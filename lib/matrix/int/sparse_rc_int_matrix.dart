/*
Copyright (C) 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
is hereby granted without fee, provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear in supporting documentation.
CERN makes no representations about the suitability of this software for any purpose.
It is provided "as is" without expressed or implied warranty.
 */
part of cern.colt.matrix;

/**
 * Sparse row-compressed 2-d matrix holding <tt>int</tt> elements. First see the
 * <a href="package-summary.html">package summary</a> and javadoc <a
 * href="package-tree.html">tree view</a> to get the broad picture.
 * <p>
 * <b>Implementation:</b>
 * <p>
 * Internally uses the standard sparse row-compressed format<br>
 * Note that this implementation is not synchronized.
 * <p>
 * <b>Memory requirements:</b>
 * <p>
 * Cells that
 * <ul>
 * <li>are never set to non-zero values do not use any memory.
 * <li>switch from zero to non-zero state do use memory.
 * <li>switch back from non-zero to zero state also do use memory. Their memory
 * is <i>not</i> automatically reclaimed (because of the lists vs. arrays).
 * Reclamation can be triggered via {@link #trimToSize()}.
 * </ul>
 * <p>
 * <tt>memory [bytes] = 4*rows + 12 * nonZeros</tt>. <br>
 * Where <tt>nonZeros = cardinality()</tt> is the number of non-zero cells.
 * Thus, a 1000 x 1000 matrix with 1000000 non-zero cells consumes 11.5 MB. The
 * same 1000 x 1000 matrix with 1000 non-zero cells consumes 15 KB.
 * <p>
 * <b>Time complexity:</b>
 * <p>
 * Getting a cell value takes time<tt> O(log nzr)</tt> where <tt>nzr</tt> is the
 * number of non-zeros of the touched row. This is usually quick, because
 * typically there are only few nonzeros per row. So, in practice, get has
 * <i>expected</i> constant time. Setting a cell value takes <i> </i>worst-case
 * time <tt>O(nz)</tt> where <tt>nzr</tt> is the total number of non-zeros in
 * the matrix. This can be extremely slow, but if you traverse coordinates
 * properly (i.e. upwards), each write is done much quicker:
 * <table>
 * <td class="PRE">
 *
 * <pre>
 * // rather quick
 * matrix.assign(0);
 * for (int row = 0; row &lt; rows; row++) {
 *     for (int column = 0; column &lt; columns; column++) {
 *         if (someCondition)
 *             matrix.setQuick(row, column, someValue);
 *     }
 * }
 *
 * // poor
 * matrix.assign(0);
 * for (int row = rows; --row &gt;= 0;) {
 *     for (int column = columns; --column &gt;= 0;) {
 *         if (someCondition)
 *             matrix.setQuick(row, column, someValue);
 *     }
 * }
 * </pre>
 *
 * </td>
 * </table>
 * If for whatever reasons you can't iterate properly, consider to create an
 * empty dense matrix, store your non-zeros in it, then call
 * <tt>sparse.assign(dense)</tt>. Under the circumstances, this is still rather
 * quick.
 * <p>
 * Fast iteration over non-zeros can be done via {@link #forEachNonZero}, which
 * supplies your function with row, column and value of each nonzero. Although
 * the internally implemented version is a bit more sophisticated, here is how a
 * quite efficient user-level matrix-vector multiplication could look like:
 * <table>
 * <td class="PRE">
 *
 * <pre>
 * // Linear algebraic y = A * x
 * A.forEachNonZero(new cern.colt.function.IntIntIntFunction() {
 *     int apply(int row, int column, int value) {
 *         y.setQuick(row, y.getQuick(row) + value * x.getQuick(column));
 *         return value;
 *     }
 * });
 * </pre>
 *
 * </td>
 * </table>
 * <p>
 * Here is how a a quite efficient user-level combined scaling operation could
 * look like:
 * <table>
 * <td class="PRE">
 *
 * <pre>
 * // Elementwise A = A + alpha*B
 * B.forEachNonZero(new cern.colt.function.IntIntIntFunction() {
 *     int apply(int row, int column, int value) {
 *         A.setQuick(row, column, A.getQuick(row, column) + alpha * value);
 *         return value;
 *     }
 * });
 * </pre>
 *
 * </td>
 * </table>
 * Method {@link #forEachWith(AbstractIntMatrix,ifunc.IntIntFunction)}
 * does just that if you supply
 * {@link ifunc.IntFunctions#plusMultSecond} as argument.
 *
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 0.9, 04/14/2000
 */
class SparseRCIntMatrix extends WrapperIntMatrix {

  static int _searchFromTo(Int32List list, int key, int from, int to) {
    while (from <= to) {
      if (list[from] == key) {
        return from;
      } else {
        from++;
        continue;
      }
    }
    return -(from + 1); // key not found.
  }

  /*
   * The elements of the matrix.
   */
  Int32List _rowPointers;

  Int32List _columnIndexes;

  Int32List _values;

  bool _columnIndexesSorted = false;

  /**
   * Constructs a matrix with a copy of the given values. <tt>values</tt> is
   * required to have the form <tt>values[row][column]</tt> and have exactly
   * the same number of columns in every row.
   * <p>
   * The values are copied. So subsequent changes in <tt>values</tt> are not
   * reflected in the matrix, and vice-versa.
   *
   * @param values
   *            The values to be filled into the new matrix.
   * @throws ArgumentError
   *             if
   *             <tt>for any 1 &lt;= row &lt; values.length: values[row].length != values[row-1].length</tt>
   *             .
   */
  factory SparseRCIntMatrix.fromList(List<Int32List> values) {
    return new SparseRCIntMatrix(values.length, values.length == 0 ? 0 : values[0].length)
      ..setAll2D(values);
  }

  /**
   * Constructs a matrix with a given number of rows and columns. All entries
   * are initially <tt>0</tt>.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @throws ArgumentError
   *             if
   *             <tt>rows<0 || columns<0 || (double)columns*rows > Integer.MAX_VALUE</tt>
   *             .
   */
  /*SparseRCIntMatrix(int rows, int columns) {
    this(rows, columns, Math.min(10 * rows, MAX_INT));
  }*/

  /**
   * Constructs a matrix with a given number of rows and columns. All entries
   * are initially <tt>0</tt>.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @param nzmax
   *            maximum number of nonzero elements
   * @throws ArgumentError
   *             if
   *             <tt>rows<0 || columns<0 || (double)columns*rows > Integer.MAX_VALUE</tt>
   *             .
   */
  factory SparseRCIntMatrix(int rows, int columns, [int nzmax=null]) {
    if (nzmax == null) {
      nzmax = Math.min(10 * rows, MAX_INT);
    }
    /*try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if (!"matrix too large".equals(exc.getMessage())) throw exc;
    }*/
    final columnIndexes = new Int32List(nzmax);
    final values = new Int32List(nzmax);
    final rowPointers = new Int32List(rows + 1);
    return new SparseRCIntMatrix._internal(rows, columns, rowPointers, columnIndexes, values);
  }

  /**
   * Constructs a matrix with indexes given in the coordinate format and
   * single value.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @param rowIndexes
   *            row indexes
   * @param columnIndexes
   *            column indexes
   * @param value
   *            numerical value, cannot be zero
   * @param removeDuplicates
   *            if true, then duplicates (if any) are removed
   * @param sortColumnIndexes
   *            if true, then column indexes are sorted
   */
  factory SparseRCIntMatrix.withValue(int rows, int columns, Int32List rowIndexes, Int32List columnIndexes, int value, bool removeDuplicates, bool sortColumnIndexes) {
    /*try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if ("matrix too large" != exc.message) {
        throw exc;
      }
    }*/
    if (rowIndexes.length != columnIndexes.length) {
      throw new ArgumentError("rowIndexes.length != columnIndexes.length");
    }
    if (value == 0) {
      throw new ArgumentError("value cannot be 0");
    }

    int nz = Math.max(rowIndexes.length, 1);
    final _columnIndexes = new Int32List(nz);
    final _values = new Int32List(nz);
    final _rowPointers = new Int32List(rows + 1);
    Int32List w = new Int32List(rows);
    int r;
    for (int k = 0; k < nz; k++) {
      w[rowIndexes[k]]++;
    }
    _cumsum(_rowPointers, w, rows);
    for (int k = 0; k < nz; k++) {
      _columnIndexes[r = w[rowIndexes[k]]++] = columnIndexes[k];
      _values[r] = value;
    }
    final m = new SparseRCIntMatrix._internal(rows, columns, _rowPointers, _columnIndexes, _values);
    if (removeDuplicates) {
      m.removeDuplicates();
    }
    if (sortColumnIndexes) {
      m.sortColumnIndexes();
    }
    return m;
  }

  /**
   * Constructs a matrix with indexes and values given in the coordinate
   * format.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @param rowIndexes
   *            row indexes
   * @param columnIndexes
   *            column indexes
   * @param values
   *            numerical values
   * @param removeDuplicates
   *            if true, then duplicates (if any) are removed
   * @param removeZeroes
   *            if true, then zeroes (if any) are removed
   * @param sortColumnIndexes
   *            if true, then column indexes are sorted
   */
  factory SparseRCIntMatrix.withValues(int rows, int columns, Int32List rowIndexes, Int32List columnIndexes, Int32List values, bool removeDuplicates, bool removeZeroes, bool sortColumnIndexes) {
    /*try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if ("matrix too large" != exc.message) {
        throw exc;
      }
    }*/
    if (rowIndexes.length != columnIndexes.length) {
      throw new ArgumentError("rowIndexes.length != columnIndexes.length");
    } else if (rowIndexes.length != values.length) {
      throw new ArgumentError("rowIndexes.length != values.length");
    }
    int nz = Math.max(rowIndexes.length, 1);
    final _columnIndexes = new Int32List(nz);
    final _values = new Int32List(nz);
    final _rowPointers = new Int32List(rows + 1);
    Int32List w = new Int32List(rows);
    int r;
    for (int k = 0; k < nz; k++) {
      w[rowIndexes[k]]++;
    }
    _cumsum(_rowPointers, w, rows);
    for (int k = 0; k < nz; k++) {
      _columnIndexes[r = w[rowIndexes[k]]++] = columnIndexes[k];
      _values[r] = values[k];
    }
    final m = new SparseRCIntMatrix._internal(rows, columns, _rowPointers, _columnIndexes, _values);
    if (removeZeroes) {
      m.removeZeroes();
    }
    if (removeDuplicates) {
      m.removeDuplicates();
    }
    if (sortColumnIndexes) {
      m.sortColumnIndexes();
    }
    return m;
  }

  /**
   * Constructs a matrix with given parameters. The arrays are not copied.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @param rowPointers
   *            row pointers
   * @param columnIndexes
   *            column indexes
   * @param values
   *            numerical values
   */
  SparseRCIntMatrix._internal(int rows, int columns, Int32List rowPointers, Int32List columnIndexes, Int32List values) : super(null) {
    try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if ("matrix too large" != exc.message) {
        throw exc;
      }
    }
    if (rowPointers.length != rows + 1) {
      throw new ArgumentError("rowPointers.length != rows + 1");
    }
    this._rowPointers = rowPointers;
    this._columnIndexes = columnIndexes;
    this._values = values;
  }

  void forEach(final ifunc.IntFunction function) {
    if (function is ifunc.IntMult) { // x[i] = mult*x[i]
      final int alpha = (function as ifunc.IntMult).multiplicator;
      if (alpha == 1) {
        return;
      }
      if (alpha == 0) {
        fill(0);
        return;
      }
      if (alpha != alpha) {
        fill(alpha); // the funny definition of isNaN(). This should better not happen.
        return;
      }

      int nz = cardinality;
      for (int j = 0; j < nz; j++) {
        _values[j] *= alpha;
      }
    } else {
      forEachNonZero((int i, int j, int value) {
        return function(value);
      });
    }
  }

  void fill(int value) {
    if (value == 0) {
      //Arrays.fill(_rowPointers, 0);
      _rowPointers.fillRange(0, _rowPointers.length, 0);
      //Arrays.fill(_columnIndexes, 0);
      _columnIndexes.fillRange(0, _columnIndexes.length, 0);
      //Arrays.fill(_values, 0);
      _values.fillRange(0, _values.length, 0);
    } else {
      int nnz = cardinality;
      for (int i = 0; i < nnz; i++) {
        _values[i] = value;
      }
    }
    return;
  }

  void copyFrom(AbstractIntMatrix source) {
    if (source == this) {
      return; // nothing to do
    }
    checkShape(source);

    if (source is SparseRCIntMatrix) {
      SparseRCIntMatrix other = source;
      //System.arraycopy(other._rowPointers, 0, _rowPointers, 0, _rows + 1);
      _rowPointers.setAll(0, other._rowPointers);
      int nzmax = other._columnIndexes.length;
      if (_columnIndexes.length < nzmax) {
        _columnIndexes = new Int32List(nzmax);
        _values = new Int32List(nzmax);
      }
      //System.arraycopy(other._columnIndexes, 0, _columnIndexes, 0, nzmax);
      _columnIndexes.setAll(0, other._columnIndexes);
      //System.arraycopy(other._values, 0, _values, 0, nzmax);
      _values.setAll(0, other._values);
      _columnIndexesSorted = other._columnIndexesSorted;
    } else if (source is SparseCCIntMatrix) {
      SparseCCIntMatrix other = source.transpose();
      _rowPointers = other.columnPointers();
      _columnIndexes = other.rowIndexes();
      _values = other.values();
      _columnIndexesSorted = true;
    } else {
      fill(0);
      source.forEachNonZero((int i, int j, int value) {
        set(i, j, value);
        return value;
      });
    }
  }

  void forEachWith(final AbstractIntMatrix y, ifunc.IntIntFunction function) {
    checkShape(y);
    if ((y is SparseRCIntMatrix) && (function == ifunc.plus)) { // x[i] = x[i] + y[i]
      SparseRCIntMatrix yy = y;

      final Int32List rowPointersY = yy._rowPointers;
      final Int32List columnIndexesY = yy._columnIndexes;
      final Int32List valuesY = yy._values;

      final Int32List rowPointersC = new Int32List(_rows + 1);
      int cnz = Math.max(_columnIndexes.length, Math.min(MAX_INT, _rowPointers[_rows] + rowPointersY[_rows]));
      final Int32List columnIndexesC = new Int32List(cnz);
      final Int32List valuesC = new Int32List(cnz);
      int nrow = _rows;
      int ncol = _columns;
      int nzmax = valuesC.length;
      if (function == ifunc.plus) { // x[i] = x[i] + y[i]
        int kc = 0;
        rowPointersC[0] = kc;
        int j1, j2;
        for (int i = 0; i < nrow; i++) {
          int ka = _rowPointers[i];
          int kb = rowPointersY[i];
          int kamax = _rowPointers[i + 1] - 1;
          int kbmax = rowPointersY[i + 1] - 1;
          while (ka <= kamax || kb <= kbmax) {
            if (ka <= kamax) {
              j1 = _columnIndexes[ka];
            } else {
              j1 = ncol + 1;
            }
            if (kb <= kbmax) {
              j2 = columnIndexesY[kb];
            } else {
              j2 = ncol + 1;
            }
            if (j1 == j2) {
              valuesC[kc] = _values[ka] + valuesY[kb];
              columnIndexesC[kc] = j1;
              ka++;
              kb++;
              kc++;
            } else if (j1 < j2) {
              columnIndexesC[kc] = j1;
              valuesC[kc] = _values[ka];
              ka++;
              kc++;
            } else if (j1 > j2) {
              columnIndexesC[kc] = j2;
              valuesC[kc] = valuesY[kb];
              kb++;
              kc++;
            }
            if (kc >= nzmax) {
              throw new ArgumentError("The number of elements in C exceeds nzmax");
            }
          }
          rowPointersC[i + 1] = kc;
        }
        this._rowPointers = rowPointersC;
        this._columnIndexes = columnIndexesC;
        this._values = valuesC;
        return;
      }
    }

    if (function is ifunc.IntPlusMultSecond) { // x[i] = x[i] + alpha*y[i]
      final int alpha = (function as ifunc.IntPlusMultSecond).multiplicator;
      if (alpha == 0) {
        return; // nothing to do
      }
      y.forEachNonZero((int i, int j, int value) {
        set(i, j, get(i, j) + alpha * value);
        return value;
      });
      return;
    }

    if (function is ifunc.IntPlusMultFirst) { // x[i] = alpha*x[i] + y[i]
      final int alpha = (function as ifunc.IntPlusMultFirst).multiplicator;
      if (alpha == 0) {
        copyFrom(y);
        return;
      }
      y.forEachNonZero((int i, int j, int value) {
        set(i, j, alpha * get(i, j) + value);
        return value;
      });
      return;
    }

    if (function == ifunc.mult) { // x[i] = x[i] * y[i]
      for (int i = _rows; --i >= 0; ) {
        int low = _rowPointers[i];
        for (int k = _rowPointers[i + 1]; --k >= low; ) {
          int j = _columnIndexes[k];
          _values[k] *= y.get(i, j);
          if (_values[k] == 0) _remove(i, j);
        }
      }
      return;
    }

    if (function == ifunc.div) { // x[i] = x[i] / y[i]

      for (int i = _rows; --i >= 0; ) {
        int low = _rowPointers[i];
        for (int k = _rowPointers[i + 1]; --k >= low; ) {
          int j = _columnIndexes[k];
          _values[k] ~/= y.get(i, j);
          if (_values[k] == 0) {
            _remove(i, j);
          }
        }
      }
      return;
    }
    super.forEachWith(y, function);
  }

  int get cardinality {
    return _rowPointers[_rows];
  }

  AbstractIntMatrix forEachNonZero(final ifunc.IntIntIntFunction function) {

    for (int i = _rows; --i >= 0; ) {
      int low = _rowPointers[i];
      for (int k = _rowPointers[i + 1]; --k >= low; ) {
        int j = _columnIndexes[k];
        int value = _values[k];
        int r = function(i, j, value);
        if (r != value) _values[k] = r;
      }
    }
    return this;
  }

  /**
   * Returns a new matrix that has the same elements as this matrix, but is in
   * a column-compressed form. This method creates a new object (not a view),
   * so changes in the returned matrix are NOT reflected in this matrix.
   *
   * @return this matrix in a column-compressed form
   */
  SparseCCIntMatrix columnCompressed() {
    SparseRCIntMatrix tr = transpose();
    SparseCCIntMatrix cc = new SparseCCIntMatrix(_rows, _columns);
    cc._rowIndexes = tr._columnIndexes;
    cc._columnPointers = tr._rowPointers;
    cc._values = tr._values;
    cc._rowIndexesSorted = true;
    return cc;
  }

  /**
   * Returns column indexes
   *
   * @return column indexes
   */
  Int32List get columnIndexes {
    return _columnIndexes;
  }

  /**
   * Returns a new matrix that has the same elements as this matrix, but is in
   * a dense form. This method creates a new object (not a view), so changes
   * in the returned matrix are NOT reflected in this matrix.
   *
   * @return this matrix in a dense form
   */
  IntMatrix dense() {
    final IntMatrix dense = new IntMatrix(_rows, _columns);
    forEachNonZero((int i, int j, int value) {
      dense.set(i, j, get(i, j));
      return value;
    });
    return dense;
  }

  int get(int row, int column) {
    //        int k = cern.colt.Sorting.binarySearchFromTo(columnIndexes, column, rowPointers[row], rowPointers[row + 1] - 1);
    int k = _searchFromTo(_columnIndexes, column, _rowPointers[row], _rowPointers[row + 1] - 1);

    int v = 0;
    if (k >= 0) v = _values[k];
    return v;
  }

  /**
   * Returns row pointers
   *
   * @return row pointers
   */
  Int32List get rowPointers {
    return _rowPointers;
  }

  /**
   * Returns a new matrix that is the transpose of this matrix. This method
   * creates a new object (not a view), so changes in the returned matrix are
   * NOT reflected in this matrix.
   *
   * @return the transpose of this matrix
   */
  SparseRCIntMatrix transpose() {
    int nnz = _rowPointers[_rows];
    Int32List w = new Int32List(_columns);
    Int32List rowPointersT = new Int32List(_columns + 1);
    Int32List columnIndexesT = new Int32List(nnz);
    Int32List valuesT = new Int32List(nnz);

    for (int p = 0; p < nnz; p++) {
      w[_columnIndexes[p]]++;
    }
    _cumsum(rowPointersT, w, _columns);
    int q;
    for (int j = 0; j < _rows; j++) {
      int high = _rowPointers[j + 1];
      for (int p = _rowPointers[j]; p < high; p++) {
        columnIndexesT[q = w[_columnIndexes[p]]++] = j;
        valuesT[q] = _values[p];
      }
    }
    SparseRCIntMatrix T = new SparseRCIntMatrix(_columns, _rows);
    T._rowPointers = rowPointersT;
    T._columnIndexes = columnIndexesT;
    T._values = valuesT;
    return T;
  }

  /**
   * Returns numerical values
   *
   * @return numerical values
   */
  Int32List get values {
    return _values;
  }

  /**
   * Returns true if column indexes are sorted, false otherwise
   *
   * @return true if column indexes are sorted, false otherwise
   */
  bool get columnIndexesSorted {
    return _columnIndexesSorted;
  }

  AbstractIntMatrix like2D(int rows, int columns) {
    return new SparseRCIntMatrix(rows, columns);
  }

  AbstractIntVector like1D(int size) {
    return new SparseIntVector(size);
  }

  /**
   * Removes (sums) duplicate entries (if any}
   */
  void removeDuplicates() {
    int nz = 0;
    int q, i;
    Int32List w = new Int32List(_columns);
    /* get workspace */
    for (i = 0; i < _columns; i++) w[i] = -1;
    /* column i not yet seen */
    for (int j = 0; j < _rows; j++) {
      q = nz;
      /* row j will start at q */
      for (int p = _rowPointers[j]; p < _rowPointers[j + 1]; p++) {
        i = _columnIndexes[p];
        /* A(i,j) is nonzero */
        if (w[i] >= q) {
          _values[w[i]] += _values[p];
          /* A(i,j) is a duplicate */
        } else {
          w[i] = nz;
          /* record where column i occurs */
          _columnIndexes[nz] = i;
          /* keep A(i,j) */
          _values[nz++] = _values[p];
        }
      }
      _rowPointers[j] = q;
      /* record start of row j */
    }
    _rowPointers[_rows] = nz;
    /* finalize A */
  }

  /**
   * Removes zero entries (if any)
   */
  void removeZeroes() {
    int nz = 0;
    for (int j = 0; j < _rows; j++) {
      int p = _rowPointers[j];
      /* get current location of row j */
      _rowPointers[j] = nz;
      /* record new location of row j */
      for ( ; p < _rowPointers[j + 1]; p++) {
        if (_values[p] != 0) {
          _values[nz] = _values[p];
          /* keep A(i,j) */
          _columnIndexes[nz++] = _columnIndexes[p];
        }
      }
    }
    _rowPointers[_rows] = nz;
    /* finalize A */
  }

  void set(int row, int column, int value) {
    //        int k = cern.colt.Sorting.binarySearchFromTo(columnIndexes, column, rowPointers[row], rowPointers[row + 1] - 1);
    int k = _searchFromTo(_columnIndexes, column, _rowPointers[row], _rowPointers[row + 1] - 1);

    if (k >= 0) { // found
      if (value == 0) _remove(row, k); else _values[k] = value;
      return;
    }

    if (value != 0) {
      k = -k - 1;
      _insert(row, column, k, value);
    }
  }

  /**
   * Sorts column indexes
   */
  void sortColumnIndexes() {
    SparseRCIntMatrix T = transpose();
    this._rows = T._rows;
    this._columns = T._columns;
    this._columnIndexes = T._columnIndexes;
    this._rowPointers = T._rowPointers;
    this._values = T._values;
    //        System.arraycopy(T.columnIndexes, 0, this.columnIndexes, 0, T.columnIndexes.length);
    //        System.arraycopy(T.rowPointers, 0, this.rowPointers, 0, T.rowPointers.length);
    //        System.arraycopy(T.values, 0, this.values, 0, T.values.length);
    T = transpose();
    this._rows = T._rows;
    this._columns = T._columns;
    this._columnIndexes = T._columnIndexes;
    this._rowPointers = T._rowPointers;
    this._values = T._values;
    _columnIndexesSorted = true;
    //        System.arraycopy(T.columnIndexes, 0, this.columnIndexes, 0, T.columnIndexes.length);
    //        System.arraycopy(T.rowPointers, 0, this.rowPointers, 0, T.rowPointers.length);
    //        System.arraycopy(T.values, 0, this.values, 0, T.values.length);
  }

  String toString() {
    StringBuffer builder = new StringBuffer();
    builder..write(_rows)..write(" x ")..write(_columns)..write(" sparse matrix, nnz = ")..write(cardinality)..write('\n');
    for (int i = 0; i < _rows; i++) {
      int high = _rowPointers[i + 1];
      for (int j = _rowPointers[i]; j < high; j++) {
        builder..write('(')..write(i)..write(',')..write(_columnIndexes[j])..write(')')..write('\t')..write(_values[j])..write('\n');
      }
    }
    return builder.toString();
  }

  void trimToSize() {
    _realloc(0);
  }

  AbstractIntVector mult(AbstractIntVector y, [AbstractIntVector z = null, final int alpha = 1, int beta = null, final bool transposeA = false]) {
    if (beta == null) {
      beta = z == null ? 1 : 0;
    }
    final int rowsA = transposeA ? _columns : _rows;
    final int columnsA = transposeA ? _rows : _columns;

    bool ignore = (z == null || !transposeA);
    if (z == null) z = new IntVector(rowsA);

    if (!(y is IntVector && z is IntVector)) {
      return super.mult(y, z, alpha, beta, transposeA);
    }

    if (columnsA != y.length || rowsA > z.length) {
      throw new ArgumentError("Incompatible args: " + ((transposeA ? dice() : this).toStringShort()) + ", " + y.toStringShort() + ", " + z.toStringShort());
    }

    IntVector zz = z as IntVector;
    final Int32List elementsZ = zz._elements;
    final int strideZ = zz.stride();
    final int zeroZ = z.index(0);

    IntVector yy = y as IntVector;
    final Int32List elementsY = yy._elements;
    final int strideY = yy.stride();
    final int zeroY = y.index(0);
    //int nthreads = ConcurrencyUtils.getNumberOfThreads();

    if (transposeA) {
      if ((!ignore) && (beta != 1.0)) {
        z.forEach(ifunc.multiply(beta));
      }

      /*if ((nthreads > 1) && (cardinality >= ConcurrencyUtils.getThreadsBeginN_2D())) {
        nthreads = 2;
        List<Future> futures = new List<Future>(nthreads);
        final Int32List result = new Int32List(rowsA);
        int k = _rows ~/ nthreads;
        for (int j = 0; j < nthreads; j++) {
          final int firstRow = j * k;
          final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
          final int threadID = j;
          futures[j] = ConcurrencyUtils.submit(() {
            if (threadID == 0) {
              for (int i = firstRow; i < lastRow; i++) {
                int high = _rowPointers[i + 1];
                int yElem = alpha * elementsY[zeroY + strideY * i];
                for (int k = _rowPointers[i]; k < high; k++) {
                  int j = _columnIndexes[k];
                  elementsZ[zeroZ + strideZ * j] += _values[k] * yElem;
                }
              }
            } else {
              for (int i = firstRow; i < lastRow; i++) {
                int high = _rowPointers[i + 1];
                int yElem = alpha * elementsY[zeroY + strideY * i];
                for (int k = _rowPointers[i]; k < high; k++) {
                  int j = _columnIndexes[k];
                  result[j] += _values[k] * yElem;
                }
              }
            }
          });
        }
        ConcurrencyUtils.waitForCompletion(futures);
        int rem = rowsA % 10;
        for (int j = rem; j < rowsA; j += 10) {
          elementsZ[zeroZ + j * strideZ] += result[j];
          elementsZ[zeroZ + (j + 1) * strideZ] += result[j + 1];
          elementsZ[zeroZ + (j + 2) * strideZ] += result[j + 2];
          elementsZ[zeroZ + (j + 3) * strideZ] += result[j + 3];
          elementsZ[zeroZ + (j + 4) * strideZ] += result[j + 4];
          elementsZ[zeroZ + (j + 5) * strideZ] += result[j + 5];
          elementsZ[zeroZ + (j + 6) * strideZ] += result[j + 6];
          elementsZ[zeroZ + (j + 7) * strideZ] += result[j + 7];
          elementsZ[zeroZ + (j + 8) * strideZ] += result[j + 8];
          elementsZ[zeroZ + (j + 9) * strideZ] += result[j + 9];
        }
        for (int j = 0; j < rem; j++) {
          elementsZ[zeroZ + j * strideZ] += result[j];
        }
      } else {*/
        for (int i = 0; i < _rows; i++) {
          int high = _rowPointers[i + 1];
          int yElem = alpha * elementsY[zeroY + strideY * i];
          for (int k = _rowPointers[i]; k < high; k++) {
            int j = _columnIndexes[k];
            elementsZ[zeroZ + strideZ * j] += _values[k] * yElem;
          }
        }
      //}

      return z;
    }

    /*if ((nthreads > 1) && (cardinality >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int zidx = zeroZ + firstRow * strideZ;
          int k = _rowPointers[firstRow];
          if (beta == 0.0) {
            for (int i = firstRow; i < lastRow; i++) {
              int sum = 0;
              int high = _rowPointers[i + 1];
              for ( ; k + 10 < high; k += 10) {
                int ind = k + 9;
                sum += _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]];
              }
              for ( ; k < high; k++) {
                sum += _values[k] * elementsY[_columnIndexes[k]];
              }
              elementsZ[zidx] = alpha * sum;
              zidx += strideZ;
            }
          } else {
            for (int i = firstRow; i < lastRow; i++) {
              int sum = 0;
              int high = _rowPointers[i + 1];
              for ( ; k + 10 < high; k += 10) {
                int ind = k + 9;
                sum += _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]];
              }
              for ( ; k < high; k++) {
                sum += _values[k] * elementsY[_columnIndexes[k]];
              }
              elementsZ[zidx] = alpha * sum + beta * elementsZ[zidx];
              zidx += strideZ;
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      int zidx = zeroZ;
      int k = _rowPointers[0];
      if (beta == 0.0) {
        for (int i = 0; i < _rows; i++) {
          int sum = 0;
          int high = _rowPointers[i + 1];
          for ( ; k + 10 < high; k += 10) {
            int ind = k + 9;
            sum += _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]];
          }
          for ( ; k < high; k++) {
            sum += _values[k] * elementsY[_columnIndexes[k]];
          }
          elementsZ[zidx] = alpha * sum;
          zidx += strideZ;
        }
      } else {
        for (int i = 0; i < _rows; i++) {
          int sum = 0;
          int high = _rowPointers[i + 1];
          for ( ; k + 10 < high; k += 10) {
            int ind = k + 9;
            sum += _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]] + _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]];
          }
          for ( ; k < high; k++) {
            sum += _values[k] * elementsY[_columnIndexes[k]];
          }
          elementsZ[zidx] = alpha * sum + beta * elementsZ[zidx];
          zidx += strideZ;
        }
      }
    //}
    return z;
  }

  AbstractIntMatrix multiply(AbstractIntMatrix B, [AbstractIntMatrix C = null, final int alpha = 1, int beta = null, final bool transposeA = false, final bool transposeB = false]) {
    if (beta == null) {
      beta = C == null ? 1 : 0;
    }
    int rowsA = _rows;
    int columnsA = _columns;
    if (transposeA) {
      rowsA = _columns;
      columnsA = _rows;
    }
    int rowsB = B.rows;
    int columnsB = B.columns;
    if (transposeB) {
      rowsB = B.columns;
      columnsB = B.rows;
    }
    int p = columnsB;
    bool ignore = (C == null);
    if (C == null) {
      if (B is SparseRCIntMatrix) {
        C = new SparseRCIntMatrix(rowsA, p, (rowsA * p));
      } else {
        C = new IntMatrix(rowsA, p);
      }
    }

    if (rowsB != columnsA) {
      throw new ArgumentError("Matrix2D inner dimensions must agree:" + toStringShort() + ", " + (transposeB ? B.dice() : B).toStringShort());
    }
    if (C.rows != rowsA || C.columns != p) {
      throw new ArgumentError("Incompatible result matrix: " + toStringShort() + ", " + (transposeB ? B.dice() : B).toStringShort() + ", " + C.toStringShort());
    }
    if (this == C || B == C) {
      throw new ArgumentError("Matrices must not be identical");
    }

    if (!ignore && beta != 1.0) {
      C.forEach(ifunc.multiply(beta));
    }

    if ((B is IntMatrix) && (C is IntMatrix)) {
      SparseRCIntMatrix AA;
      if (transposeA) {
        AA = transpose();
      } else {
        AA = this;
      }
      IntMatrix BB;
      if (transposeB) {
        BB = B.dice() as IntMatrix;
      } else {
        BB = B;
      }

      IntMatrix CC = C;
      Int32List rowPointersA = AA._rowPointers;
      Int32List columnIndexesA = AA._columnIndexes;
      Int32List valuesA = AA._values;

      for (int ii = 0; ii < rowsA; ii++) {
        int highA = rowPointersA[ii + 1];
        for (int ka = rowPointersA[ii]; ka < highA; ka++) {
          int scal = valuesA[ka] * alpha;
          int jj = columnIndexesA[ka];
          CC.row(ii).forEachWith(BB.row(jj), ifunc.plusMultSecond(scal));
        }
      }
    } else if ((B is SparseRCIntMatrix) && (C is SparseRCIntMatrix)) {
      SparseRCIntMatrix AA;
      SparseRCIntMatrix BB;
      SparseRCIntMatrix CC = C;
      if (transposeA) {
        AA = transpose();
      } else {
        AA = this;
      }
      if (transposeB) {
        BB = B.transpose();
      } else {
        BB = B;
      }

      Int32List rowPointersA = AA._rowPointers;
      Int32List columnIndexesA = AA._columnIndexes;
      Int32List valuesA = AA._values;

      Int32List rowPointersB = BB._rowPointers;
      Int32List columnIndexesB = BB._columnIndexes;
      Int32List valuesB = BB._values;

      Int32List rowPointersC = CC._rowPointers;
      Int32List columnIndexesC = CC._columnIndexes;
      Int32List valuesC = CC._values;
      int nzmax = valuesC.length;

      Int32List iw = new Int32List(columnsB + 1);
      for (int i = 0; i < iw.length; i++) {
        iw[i] = -1;
      }
      int len = -1;
      for (int ii = 0; ii < rowsA; ii++) {
        int highA = rowPointersA[ii + 1];
        for (int ka = rowPointersA[ii]; ka < highA; ka++) {
          int scal = valuesA[ka] * alpha;
          int jj = columnIndexesA[ka];
          int highB = rowPointersB[jj + 1];
          for (int kb = rowPointersB[jj]; kb < highB; kb++) {
            int jcol = columnIndexesB[kb];
            int jpos = iw[jcol];
            if (jpos == -1) {
              len++;
              if (len >= nzmax) {
                throw new ArgumentError("The max number of nonzero elements in C is too small.");
              }
              columnIndexesC[len] = jcol;
              iw[jcol] = len;
              valuesC[len] = scal * valuesB[kb];
            } else {
              valuesC[jpos] += scal * valuesB[kb];
            }
          }
        }
        for (int k = rowPointersC[ii]; k < len + 1; k++) {
          iw[columnIndexesC[k]] = -1;
        }
        rowPointersC[ii + 1] = len + 1;

        //                int length = rowPointersC[ii + 1] - rowPointersC[ii];
        //                IntMatrix1D columnIndexesCPart = columnIndexesC.viewPart(rowPointersC[ii], length);
        //                Int32List indexes = cern.colt.matrix.tint.algo.IntSorting.quickSort.sortIndex(columnIndexesCPart);
        //                Arrays.sort(columnIndexesCElements, rowPointersC[ii], rowPointersC[ii + 1]);
        //                IntMatrix1D valuesCPart = valuesC.viewPart(rowPointersC[ii], length).viewSelection(indexes);
        //                valuesC.viewPart(rowPointersC[ii], length).assign(valuesCPart);
      }
      //            CC.columnIndexes.elements((Int32List) columnIndexesC.elements());
      //            CC.columnIndexes.setSize(columnIndexesSize);
      //            CC.values.elements((Int32List) valuesC.elements());
      //            CC.values.setSize(columnIndexesSize);
    } else {
      if (transposeB) {
        B = B.dice();
      }
      // cache views
      final List<AbstractIntVector> Brows = new List<AbstractIntVector>(columnsA);
      for (int i = columnsA; --i >= 0; ) {
        Brows[i] = B.row(i);
      }
      final List<AbstractIntVector> Crows = new List<AbstractIntVector>(rowsA);
      for (int i = rowsA; --i >= 0; ) {
        Crows[i] = C.row(i);
      }

      final ifunc.IntPlusMultSecond fun = new ifunc.IntPlusMultSecond.plusMult(0);

      final Int32List columnIndexesA = _columnIndexes;
      final Int32List valuesA = _values;
      for (int i = _rows; --i >= 0; ) {
        int low = _rowPointers[i];
        for (int k = _rowPointers[i + 1]; --k >= low; ) {
          int j = columnIndexesA[k];
          fun.multiplicator = valuesA[k] * alpha;
          if (!transposeA) Crows[i].forEachWith(Brows[j], fun); else Crows[j].forEachWith(Brows[i], fun);
        }
      }
    }
    return C;
  }

  static int _cumsum(Int32List p, Int32List c, int n) {
    int nz = 0;
    int nz2 = 0;
    for (int k = 0; k < n; k++) {
      p[k] = nz;
      nz += c[k];
      nz2 += c[k];
      c[k] = p[k];
    }
    p[n] = nz;
    return (nz2);
  }

  void _realloc(int nzmax) {
    if (nzmax <= 0) nzmax = _rowPointers[_rows];
    Int32List columnIndexesNew = new Int32List(nzmax);
    int length = Math.min(nzmax, _columnIndexes.length);
    //System.arraycopy(_columnIndexes, 0, columnIndexesNew, 0, length);
    columnIndexesNew.setAll(0, _columnIndexes);
    _columnIndexes = columnIndexesNew;
    Int32List valuesNew = new Int32List(nzmax);
    length = Math.min(nzmax, _values.length);
    //System.arraycopy(_values, 0, valuesNew, 0, length);
    valuesNew.setAll(0, _values);
    _values = valuesNew;
  }

  AbstractIntMatrix _getContent() {
    return this;
  }

  void _insert(int row, int column, int index, int value) {
    List<int> columnIndexesList = new List<int>.from(_columnIndexes);
    //columnIndexesList._setSizeRaw(_rowPointers[_rows]);
    List<int> valuesList = new List<int>.from(_values);
    //valuesList._setSizeRaw(_rowPointers[_rows]);
    columnIndexesList.insert(index, column);
    valuesList.insert(index, value);
    for (int i = _rowPointers.length; --i > row; ) {
      _rowPointers[i]++;
    }
    _columnIndexes = new Int32List.fromList(columnIndexesList);
    _values = new Int32List.fromList(valuesList);
  }

  void _remove(int row, int index) {
    List<int> columnIndexesList = new List<int>.from(_columnIndexes);
    //columnIndexesList._setSizeRaw(_rowPointers[_rows]);
    List<int> valuesList = new List<int>.from(_values);
    //valuesList._setSizeRaw(_rowPointers[_rows]);
    columnIndexesList.remove(index);
    valuesList.remove(index);
    for (int i = _rowPointers.length; --i > row; ) {
      _rowPointers[i]--;
    }
    _columnIndexes = new Int32List.fromList(columnIndexesList);
    _values = new Int32List.fromList(valuesList);
  }

}
