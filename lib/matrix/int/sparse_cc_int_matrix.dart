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
 * Sparse column-compressed 2-d matrix holding <tt>int</tt> elements. First see
 * the <a href="package-summary.html">package summary</a> and javadoc <a
 * href="package-tree.html">tree view</a> to get the broad picture.
 * <p>
 * <b>Implementation:</b>
 * <p>
 * Internally uses the standard sparse column-compressed format. <br>
 * Note that this implementation is not synchronized.
 * <p>
 * Cells that
 * <ul>
 * <li>are never set to non-zero values do not use any memory.
 * <li>switch from zero to non-zero state do use memory.
 * <li>switch back from non-zero to zero state also do use memory. Their memory
 * is <i>not</i> automatically reclaimed. Reclamation can be triggered via
 * {@link #trimToSize()}.
 * </ul>
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
 * for (int column = 0; column &lt; columns; column++) {
 *     for (int row = 0; row &lt; rows; row++) {
 *         if (someCondition)
 *             matrix.setQuick(row, column, someValue);
 *     }
 * }
 *
 * // poor
 * matrix.assign(0);
 * for (int column = columns; --column &gt;= 0;) {
 *     for (int row = rows; --row &gt;= 0;) {
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
 * @author Piotr Wendykier
 *
 */
class SparseCCIntMatrix extends WrapperIntMatrix {
  /*
   * Internal storage.
   */
  Int32List _columnPointers;

  Int32List _rowIndexes;

  Int32List _values;

  bool _rowIndexesSorted = false;

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
  factory SparseCCIntMatrix.fromList(List<Int32List> values) {
    return new SparseCCIntMatrix(values.length, values[0].length)
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
   *             if <tt>rows<0 || columns<0</tt> .
   */
  /*factory SparseCCIntMatrix(int rows, int columns) {
    return new SparseCCIntMatrix(rows, columns, Math.min(10 * rows, MAX_INT));
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
   *             if <tt>rows<0 || columns<0</tt> .
   */
  factory SparseCCIntMatrix(int rows, int columns, [int nzmax = null]) {
    if (nzmax == null) {
      nzmax = Math.min(10 * rows, MAX_INT);
    }
    final rowIndexes = new Int32List(nzmax);
    final values = new Int32List(nzmax);
    final columnPointers = new Int32List(columns + 1);
    return new SparseCCIntMatrix._internal(rows, columns, rowIndexes, columnPointers, values);
  }

  SparseCCIntMatrix._internal(int rows, int columns, Int32List rowIndexes, Int32List columnPointers, Int32List values) : super(null) {
    try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if ("matrix too large" != exc.message) {
        throw exc;
      }
    }
    _rowIndexes = rowIndexes;
    _values = values;
    _columnPointers = columnPointers;
  }

  /**
   * Constructs a matrix with indexes given in the coordinate format and a
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
   *            numerical value
   * @param removeDuplicates
   *            if true, then duplicates (if any) are removed
   * @param sortRowIndexes
   *            if true, then row indexes are sorted
   */
  factory SparseCCIntMatrix.withValue(int rows, int columns, Int32List rowIndexes, Int32List columnIndexes, int value, bool removeDuplicates, bool sortRowIndexes) {//: super(null) {
    /*try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if (!"matrix too large".equals(exc.getMessage()))
        throw exc;
    }*/
    if (rowIndexes.length != columnIndexes.length) {
      throw new ArgumentError("rowIndexes.length != columnIndexes.length");
    }
    if (value == 0) {
      throw new ArgumentError("value cannot be 0");
    }

    int nz = Math.max(rowIndexes.length, 1);
    final _rowIndexes = new Int32List(nz);
    final _values = new Int32List(nz);
    final _columnPointers = new Int32List(columns + 1);
    Int32List w = new Int32List(columns);
    int r;
    for (int k = 0; k < nz; k++) {
      w[columnIndexes[k]]++;
    }
    _cumsum(_columnPointers, w, columns);
    for (int k = 0; k < nz; k++) {
      _rowIndexes[r = w[columnIndexes[k]]++] = rowIndexes[k];
      _values[r] = value;
    }
    final m = new SparseCCIntMatrix._internal(rows, columns, _rowIndexes, _columnPointers, _values);
    if (removeDuplicates) {
      m.removeDuplicates();
    }
    if (sortRowIndexes) {
      m.sortRowIndexes();
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
   * @param sortRowIndexes
   *            if true, then row indexes are sorted
   */
  factory SparseCCIntMatrix.withValues(int rows, int columns, Int32List rowIndexes, Int32List columnIndexes, Int32List values, bool removeDuplicates, bool removeZeroes, bool sortRowIndexes) {//: super(null) {
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
    final _rowIndexes = new Int32List(nz);
    final _values = new Int32List(nz);
    final _columnPointers = new Int32List(columns + 1);
    Int32List w = new Int32List(columns);
    int r;
    for (int k = 0; k < nz; k++) {
      w[columnIndexes[k]]++;
    }
    _cumsum(_columnPointers, w, columns);
    for (int k = 0; k < nz; k++) {
      _rowIndexes[r = w[columnIndexes[k]]++] = rowIndexes[k];
      _values[r] = values[k];
    }
    final m = new SparseCCIntMatrix._internal(rows, columns, _rowIndexes, _columnPointers, _values);
    if (removeDuplicates) {
      m.removeDuplicates();
    }
    if (sortRowIndexes) {
      m.sortRowIndexes();
    }
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

      final Int32List valuesE = _values;
      int nz = cardinality;
      for (int j = 0; j < nz; j++) {
        valuesE[j] *= alpha;
      }
    } else {
      forEachNonZero((int i, int j, int value) {
        return function(value);
      });
    }
  }

  void fill(int value) {
    if (value == 0) {
      //Arrays.fill(_rowIndexes, 0);
      _rowIndexes.fillRange(0, _rowIndexes.length, 0);
      //Arrays.fill(_columnPointers, 0);
      _columnPointers.fillRange(0, _columnPointers.length, 0);
      //Arrays.fill(_values, 0);
      _values.fillRange(0, _values.length, 0);
    } else {
      int nnz = cardinality;
      for (int i = 0; i < nnz; i++) {
        _values[i] = value;
      }
    }
  }

  void copyFrom(AbstractIntMatrix source) {
    if (source == this) {
      return; // nothing to do
    }
    checkShape(source);

    if (source is SparseCCIntMatrix) {
      SparseCCIntMatrix other = source;
      //System.arraycopy(other.columnPointers(), 0, _columnPointers, 0, _columns + 1);
      _columnPointers.setAll(0, other.columnPointers);
      int nzmax = other.rowIndexes.length;
      if (_rowIndexes.length < nzmax) {
        _rowIndexes = new Int32List(nzmax);
        _values = new Int32List(nzmax);
      }
      //System.arraycopy(other.rowIndexes(), 0, _rowIndexes, 0, nzmax);
      _rowIndexes.setAll(0, other.rowIndexes);
      //System.arraycopy(other.values(), 0, _values, 0, nzmax);
      _values.setAll(0, other.values);
      _rowIndexesSorted = other._rowIndexesSorted;
    } else if (source is SparseRCIntMatrix) {
      SparseRCIntMatrix other = source.transpose();
      _columnPointers = other.rowPointers;
      _rowIndexes = other.columnIndexes;
      _values = other.values;
      _rowIndexesSorted = true;
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

    if ((y is SparseCCIntMatrix) && (function == ifunc.plus)) { // x[i] = x[i] + y[i]
      SparseCCIntMatrix yy = y;
      int p,
          j,
          nz = 0,
          anz,
          m,
          n,
          bnz;
      Int32List Cp, Ci, Bp, w, x, Cx;
      m = _rows;
      anz = _columnPointers[_columns];
      n = yy._columns;
      Bp = yy._columnPointers;
      bnz = Bp[n];
      w = new Int32List(m);
      /* get workspace */
      x = new Int32List(m);
      /* get workspace */
      SparseCCIntMatrix C = new SparseCCIntMatrix(m, n, anz + bnz);
      /* allocate result*/
      Cp = C._columnPointers;
      Ci = C._rowIndexes;
      Cx = C._values;
      for (j = 0; j < n; j++) {
        Cp[j] = nz;
        /* column j of C starts here */
        nz = _scatter(this, j, 1, w, x, j + 1, C, nz);
        /* alpha*A(:,j)*/
        nz = _scatter(yy, j, 1, w, x, j + 1, C, nz);
        /* beta*B(:,j) */
        for (p = Cp[j]; p < nz; p++) Cx[p] = x[Ci[p]];
      }
      Cp[n] = nz;
      /* finalize the last column of C */
      _rowIndexes = Ci;
      _columnPointers = Cp;
      _values = Cx;
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
      final Int32List rowIndexesA = _rowIndexes;
      final Int32List columnPointersA = _columnPointers;
      final Int32List valuesA = _values;
      for (int j = _columns; --j >= 0; ) {
        int low = columnPointersA[j];
        for (int k = columnPointersA[j + 1]; --k >= low; ) {
          int i = rowIndexesA[k];
          valuesA[k] *= y.get(i, j);
          if (valuesA[k] == 0) _remove(i, j);
        }
      }
      return;
    }

    if (function == ifunc.div) { // x[i] = x[i] / y[i]
      final Int32List rowIndexesA = _rowIndexes;
      final Int32List columnPointersA = _columnPointers;
      final Int32List valuesA = _values;

      for (int j = _columns; --j >= 0; ) {
        int low = columnPointersA[j];
        for (int k = columnPointersA[j + 1]; --k >= low; ) {
          int i = rowIndexesA[k];
          valuesA[k] ~/= y.get(i, j);
          if (valuesA[k] == 0) _remove(i, j);
        }
      }
      return;
    }
    super.forEachWith(y, function);
  }

  int get cardinality {
    return _columnPointers[_columns];
  }

  void forEachNonZero(final ifunc.IntIntIntFunction function) {
    final Int32List rowIndexesA = _rowIndexes;
    final Int32List columnPointersA = _columnPointers;
    final Int32List valuesA = _values;

    for (int j = _columns; --j >= 0; ) {
      int low = columnPointersA[j];
      for (int k = columnPointersA[j + 1]; --k >= low; ) {
        int i = rowIndexesA[k];
        int value = valuesA[k];
        int r = function(i, j, value);
        valuesA[k] = r;
      }
    }
  }

  /**
   * Returns column pointers
   *
   * @return column pointers
   */
  Int32List get columnPointers => _columnPointers;

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
    //        int k = cern.colt.Sorting.binarySearchFromTo(dcs.i, row, dcs.p[column], dcs.p[column + 1] - 1);
    int k = _searchFromTo(_rowIndexes, row, _columnPointers[column], _columnPointers[column + 1] - 1);
    int v = 0;
    if (k >= 0) v = _values[k];
    return v;
  }

  /**
   * Returns a new matrix that has the same elements as this matrix, but is in
   * a row-compressed form. This method creates a new object (not a view), so
   * changes in the returned matrix are NOT reflected in this matrix.
   *
   * @return this matrix in a row-compressed form
   */
  SparseRCIntMatrix rowCompressed() {
    SparseCCIntMatrix tr = transpose();
    SparseRCIntMatrix rc = new SparseRCIntMatrix(_rows, _columns);
    rc._columnIndexes = tr._rowIndexes;
    rc._rowPointers = tr._columnPointers;
    rc._values = tr._values;
    rc._columnIndexesSorted = true;
    return rc;
  }

  /**
   * Returns row indexes;
   *
   * @return row indexes
   */
  Int32List get rowIndexes => _rowIndexes;

  /**
   * Returns a new matrix that is the transpose of this matrix. This method
   * creates a new object (not a view), so changes in the returned matrix are
   * NOT reflected in this matrix.
   *
   * @return the transpose of this matrix
   */
  SparseCCIntMatrix transpose() {
    int p, q, j, n, m;
    Int32List Cp, Ci, Ap, Ai, w, Cx, Ax;
    m = _rows;
    n = _columns;
    Ap = _columnPointers;
    Ai = _rowIndexes;
    Ax = _values;
    SparseCCIntMatrix C = new SparseCCIntMatrix(_columns, _rows, Ai.length);
    /* allocate result */
    w = new Int32List(m);
    /* get workspace */
    Cp = C._columnPointers;
    Ci = C._rowIndexes;
    Cx = C._values;
    for (p = 0; p < Ap[n]; p++) w[Ai[p]]++;
    /* row counts */
    _cumsum(Cp, w, m);
    /* row pointers */
    for (j = 0; j < n; j++) {
      for (p = Ap[j]; p < Ap[j + 1]; p++) {
        Ci[q = w[Ai[p]]++] = j;
        /* place A(i,j) as entry C(j,i) */
        Cx[q] = Ax[p];
      }
    }
    return C;
  }

  /**
   * Returns numerical values
   *
   * @return numerical values
   */
  Int32List get values => _values;

  /**
   * Returns true if row indexes are sorted, false otherwise
   *
   * @return true if row indexes are sorted, false otherwise
   */
  bool get rowIndexesSorted => _rowIndexesSorted;

  AbstractIntMatrix like2D(int rows, int columns) {
    return new SparseCCIntMatrix(rows, columns);
  }

  AbstractIntVector like1D(int size) {
    return new SparseIntVector(size);
  }

  void set(int row, int column, int value) {
    //        int k = cern.colt.Sorting.binarySearchFromTo(dcs.i, row, dcs.p[column], dcs.p[column + 1] - 1);
    int k = _searchFromTo(_rowIndexes, row, _columnPointers[column], _columnPointers[column + 1] - 1);

    if (k >= 0) { // found
      if (value == 0) _remove(column, k); else _values[k] = value;
      return;
    }

    if (value != 0) {
      k = -k - 1;
      _insert(row, column, k, value);
    }
  }

  /**
   * Sorts row indexes
   */
  void sortRowIndexes() {
    SparseCCIntMatrix tr = transpose();
    tr = tr.transpose();
    _columnPointers = tr._columnPointers;
    _rowIndexes = tr._rowIndexes;
    _values = tr._values;
    _rowIndexesSorted = true;
  }

  /**
   * Removes (sums) duplicate entries (if any}
   */
  void removeDuplicates() {
    int i,
        j,
        p,
        q,
        nz = 0,
        n,
        m;
    Int32List Ap, Ai, w, Ax;
    /* check inputs */
    m = _rows;
    n = _columns;
    Ap = _columnPointers;
    Ai = _rowIndexes;
    Ax = _values;
    w = new Int32List(m);
    /* get workspace */
    for (i = 0; i < m; i++) w[i] = -1;
    /* row i not yet seen */
    for (j = 0; j < n; j++) {
      q = nz;
      /* column j will start at q */
      for (p = Ap[j]; p < Ap[j + 1]; p++) {
        i = Ai[p];
        /* A(i,j) is nonzero */
        if (w[i] >= q) {
          Ax[w[i]] += Ax[p];
          /* A(i,j) is a duplicate */
        } else {
          w[i] = nz;
          /* record where row i occurs */
          Ai[nz] = i;
          /* keep A(i,j) */
          Ax[nz++] = Ax[p];
        }
      }
      Ap[j] = q;
      /* record start of column j */
    }
    Ap[n] = nz;
    /* finalize A */
  }

  /**
   * Removes zero entries (if any)
   */
  void removeZeroes() {
    int j,
        p,
        nz = 0,
        n;
    Int32List Ap, Ai, Ax;
    n = _columns;
    Ap = _columnPointers;
    Ai = _rowIndexes;
    Ax = _values;
    for (j = 0; j < n; j++) {
      p = Ap[j];
      /* get current location of col j */
      Ap[j] = nz;
      /* record new location of col j */
      for ( ; p < Ap[j + 1]; p++) {
        if (Ax[p] != 0) {
          Ax[nz] = Ax[p];
          /* keep A(i,j) */
          Ai[nz++] = Ai[p];
        }
      }
    }
    Ap[n] = nz;
    /* finalize A */
  }

  void trimToSize() {
    _realloc(0);
  }

  String toString() {
    StringBuffer builder = new StringBuffer();
    builder..write(_rows)..write(" x ")..write(_columns)..write(" sparse matrix, nnz = ")..write(cardinality)..write('\n');
    for (int i = 0; i < _columns; i++) {
      int high = _columnPointers[i + 1];
      for (int j = _columnPointers[i]; j < high; j++) {
        builder..write('(')..write(_rowIndexes[j])..write(',')..write(i)..write(')')..write('\t')..write(_values[j])..write('\n');
      }
    }
    return builder.toString();
  }

  AbstractIntVector mult(AbstractIntVector y, [AbstractIntVector z = null, final int alpha = 1, int beta = null, final bool transposeA = false]) {
    if (beta == null) {
      beta = z == null ? 1 : 0;
    }
    final int rowsA = transposeA ? _columns : _rows;
    final int columnsA = transposeA ? _rows : _columns;

    bool ignore = (z == null || transposeA);
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
    final int zeroZ = zz.index(0);

    IntVector yy = y as IntVector;
    final Int32List elementsY = yy._elements;
    final int strideY = yy.stride();
    final int zeroY = yy.index(0);

    final Int32List rowIndexesA = _rowIndexes;
    final Int32List columnPointersA = _columnPointers;
    final Int32List valuesA = _values;

    int zidx = zeroZ;
    //int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if (!transposeA) {
      if ((!ignore) && (beta != 1)) {
        z.forEach(ifunc.multiply(beta));
      }

      /*if ((nthreads > 1) && (cardinality() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
        nthreads = 2;
        List<Future> futures = new List<Future>(nthreads);
        final Int32List result = new Int32List(rowsA);
        int k = _columns / nthreads;
        for (int j = 0; j < nthreads; j++) {
          final int firstColumn = j * k;
          final int lastColumn = (j == nthreads - 1) ? _columns : firstColumn + k;
          final int threadID = j;
          futures[j] = ConcurrencyUtils.submit(() {
            if (threadID == 0) {
              for (int i = firstColumn; i < lastColumn; i++) {
                int high = columnPointersA[i + 1];
                int yElem = elementsY[zeroY + strideY * i];
                for (int k = columnPointersA[i]; k < high; k++) {
                  int j = rowIndexesA[k];
                  elementsZ[zeroZ + strideZ * j] += alpha * valuesA[k] * yElem;
                }
              }
            } else {
              for (int i = firstColumn; i < lastColumn; i++) {
                int high = columnPointersA[i + 1];
                int yElem = elementsY[zeroY + strideY * i];
                for (int k = columnPointersA[i]; k < high; k++) {
                  int j = rowIndexesA[k];
                  result[j] += alpha * valuesA[k] * yElem;
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
        for (int i = 0; i < _columns; i++) {
          int high = columnPointersA[i + 1];
          int yElem = elementsY[zeroY + strideY * i];
          for (int k = columnPointersA[i]; k < high; k++) {
            int j = rowIndexesA[k];
            elementsZ[zeroZ + strideZ * j] += alpha * valuesA[k] * yElem;
          }
        }
      //}
    } else {
      /*if ((nthreads > 1) && (cardinality() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
        List<Future> futures = new List<Future>(nthreads);
        int k = _columns / nthreads;
        for (int j = 0; j < nthreads; j++) {
          final int firstColumn = j * k;
          final int lastColumn = (j == nthreads - 1) ? _columns : firstColumn + k;
          futures[j] = ConcurrencyUtils.submit(() {
            int zidx = zeroZ + firstColumn * strideZ;
            int k = _columnPointers[firstColumn];
            for (int i = firstColumn; i < lastColumn; i++) {
              int sum = 0;
              int high = _columnPointers[i + 1];
              for ( ; k + 10 < high; k += 10) {
                int ind = k + 9;
                sum += valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]];
              }
              for ( ; k < high; k++) {
                sum += valuesA[k] * elementsY[_rowIndexes[k]];
              }
              elementsZ[zidx] = alpha * sum + beta * elementsZ[zidx];
              zidx += strideZ;
            }
          });
        }
        ConcurrencyUtils.waitForCompletion(futures);
      } else {*/
        int k = _columnPointers[0];
        for (int i = 0; i < _columns; i++) {
          int sum = 0;
          int high = _columnPointers[i + 1];
          for ( ; k + 10 < high; k += 10) {
            int ind = k + 9;
            sum += valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]];
          }
          for ( ; k < high; k++) {
            sum += valuesA[k] * elementsY[_rowIndexes[k]];
          }
          elementsZ[zidx] = alpha * sum + beta * elementsZ[zidx];
          zidx += strideZ;
        }
      //}
    }
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
      if (B is SparseCCIntMatrix) {
        C = new SparseCCIntMatrix(rowsA, p, (rowsA * p));
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
      SparseCCIntMatrix AA;
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
      Int32List columnPointersA = AA._columnPointers;
      Int32List rowIndexesA = AA._rowIndexes;
      Int32List valuesA = AA._values;

      int zeroB = BB.index(0, 0);
      int rowStrideB = BB.rowStride;
      int columnStrideB = BB.columnStride;
      Int32List elementsB = BB._elements;

      int zeroC = CC.index(0, 0);
      int rowStrideC = CC.rowStride;
      int columnStrideC = CC.columnStride;
      Int32List elementsC = CC._elements;

      for (int jj = 0; jj < columnsB; jj++) {
        for (int kk = 0; kk < columnsA; kk++) {
          int high = columnPointersA[kk + 1];
          int yElem = elementsB[zeroB + kk * rowStrideB + jj * columnStrideB];
          for (int ii = columnPointersA[kk]; ii < high; ii++) {
            int j = rowIndexesA[ii];
            elementsC[zeroC + j * rowStrideC + jj * columnStrideC] += valuesA[ii] * yElem;
          }
        }
      }
      if (alpha != 1.0) {
        C.forEach(ifunc.multiply(alpha));
      }

    } else if ((B is SparseCCIntMatrix) && (C is SparseCCIntMatrix)) {
      SparseCCIntMatrix AA;
      if (transposeA) {
        AA = transpose();
      } else {
        AA = this;
      }
      SparseCCIntMatrix BB = B;
      if (transposeB) {
        BB = BB.transpose();
      }
      SparseCCIntMatrix CC = C;
      int j,
          nz = 0,
          m,
          n;
      Int32List Cp, Ci, Bp, w, Bi, x, Bx, Cx;
      m = rowsA;
      n = columnsB;
      Bp = BB._columnPointers;
      Bi = BB._rowIndexes;
      Bx = BB._values;
      w = new Int32List(m);
      /* get workspace */
      x = new Int32List(m);
      /* get workspace */
      Cp = CC._columnPointers;
      Ci = CC._rowIndexes;
      Cx = CC._values;
      for (j = 0; j < n; j++) {
        int nzmaxC = CC._rowIndexes.length;
        if (nz + m > nzmaxC) {
          nzmaxC = 2 * nzmaxC + m;
          Int32List rowIndexesNew = new Int32List(nzmaxC);
          //System.arraycopy(Ci, 0, rowIndexesNew, 0, Ci.length);
          rowIndexesNew.setAll(0, Ci);
          Ci = rowIndexesNew;
          Int32List valuesNew = new Int32List(nzmaxC);
          //System.arraycopy(Cx, 0, valuesNew, 0, Cx.length);
          valuesNew.setAll(0, Cx);
          Cx = valuesNew;
        }
        Cp[j] = nz;
        /* column j of C starts here */
        for (p = Bp[j]; p < Bp[j + 1]; p++) {
          nz = _scatter(AA, Bi[p], Bx[p], w, x, j + 1, CC, nz);
        }
        for (p = Cp[j]; p < nz; p++) Cx[p] = x[Ci[p]];
      }
      Cp[n] = nz;
      /* finalize the last column of C */
      if (alpha != 1.0) {
        CC.forEach(ifunc.multiply(alpha));
      }
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

      final Int32List rowIndexesA = _rowIndexes;
      final Int32List columnPointersA = _columnPointers;
      final Int32List valuesA = _values;
      for (int i = _columns; --i >= 0; ) {
        int low = columnPointersA[i];
        for (int k = columnPointersA[i + 1]; --k >= low; ) {
          int j = rowIndexesA[k];
          fun.multiplicator = valuesA[k] * alpha;
          if (!transposeA) Crows[j].forEachWith(Brows[i], fun); else Crows[i].forEachWith(Brows[j], fun);
        }
      }
    }
    return C;
  }

  AbstractIntMatrix _getContent() {
    return this;
  }

  void _insert(int row, int column, int index, int value) {
    Int32List rowIndexesList = new Int32List.from(_rowIndexes);
    //rowIndexesList._setSizeRaw(_columnPointers[_columns]);
    Int32List valuesList = new Int32List.from(_values);
    //valuesList._setSizeRaw(_columnPointers[_columns]);
    rowIndexesList.insert(index, row);
    valuesList.insert(index, value);
    for (int i = _columnPointers.length; --i > column; ) {
      _columnPointers[i]++;
    }
    _rowIndexes = new Int32List.fromList(rowIndexesList);
    _values = new Int32List.fromList(valuesList);
  }

  void _remove(int column, int index) {
    Int32List rowIndexesList = new Int32List.from(_rowIndexes);
    Int32List valuesList = new Int32List.from(_values);
    rowIndexesList.remove(index);
    valuesList.remove(index);
    for (int i = _columnPointers.length; --i > column; ) {
      _columnPointers[i]--;
    }
    _rowIndexes = new Int32List.fromList(rowIndexesList);
    _values = new Int32List.fromList(valuesList);
  }

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
    if (nzmax <= 0) nzmax = _columnPointers[_columns];
    Int32List rowIndexesNew = new Int32List(nzmax);
    int length = Math.min(nzmax, _rowIndexes.length);
    //System.arraycopy(_rowIndexes, 0, rowIndexesNew, 0, length);
    rowIndexesNew.setAll(0, _rowIndexes);
    _rowIndexes = rowIndexesNew;
    Int32List valuesNew = new Int32List(nzmax);
    length = Math.min(nzmax, _values.length);
    //System.arraycopy(_values, 0, valuesNew, 0, length);
    valuesNew.setAll(0, _values);
    _values = valuesNew;
  }

  int _scatter(SparseCCIntMatrix A, int j, int beta, Int32List w, Int32List x, int mark, SparseCCIntMatrix C, int nz) {
    int i, p;
    Int32List Ap, Ai, Ci;
    Int32List Ax;
    Ap = A._columnPointers;
    Ai = A._rowIndexes;
    Ax = A._values;
    Ci = C._rowIndexes;
    for (p = Ap[j]; p < Ap[j + 1]; p++) {
      i = Ai[p];
      /* A(i,j) is nonzero */
      if (w[i] < mark) {
        w[i] = mark;
        /* i is new entry in column j */
        Ci[nz++] = i;
        /* add i to pattern of C(:,j) */
        if (x != null) x[i] = beta * Ax[p];
        /* x(i) = beta*A(i,j) */
      } else if (x != null) x[i] += beta * Ax[p];
      /* i exists in C(:,j) already */
    }
    return nz;
  }

  Object clone() {
    return new SparseCCIntMatrix._internal(rows, columns, _rowIndexes, _columnPointers, _values);
  }
}
