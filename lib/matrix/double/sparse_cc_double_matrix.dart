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
 * Sparse column-compressed 2-d matrix holding <tt>double</tt> elements. First
 * see the <a href="package-summary.html">package summary</a> and javadoc <a
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
 * A.forEachNonZero(new cern.colt.function.IntIntDoubleFunction() {
 *     double apply(int row, int column, double value) {
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
 * B.forEachNonZero(new cern.colt.function.IntIntDoubleFunction() {
 *     double apply(int row, int column, double value) {
 *         A.setQuick(row, column, A.getQuick(row, column) + alpha * value);
 *         return value;
 *     }
 * });
 * </pre>
 *
 * </td>
 * </table>
 * Method
 * {@link #assign(DoubleMatrix2D,func.DoubleDoubleFunction)}
 * does just that if you supply
 * {@link func#plusMultSecond} as argument.
 *
 *
 * @author Piotr Wendykier
 *
 */
class SparseCCDoubleMatrix2D extends WrapperDoubleMatrix2D {

  /*
   * Internal storage.
   */
  Dcs _dcs;

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
  factory SparseCCDoubleMatrix2D.fromList(List<Float64List> values) {
    return new SparseCCDoubleMatrix2D.sized(values.length, values[0].length)
      ..assignValues2D(values);
  }

  /**
   * Constructs a matrix with a given internal storage.
   *
   * @param dcs
   *            internal storage.
   */
  SparseCCDoubleMatrix2D(Dcs dcs) : super(null) {
    try {
      _setUp(dcs.m, dcs.n);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if ("matrix too large" != exc.message) {
        throw exc;
      }
    }
    this._dcs = dcs;
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
  /*SparseCCDoubleMatrix2D(int rows, int columns) {
        this(rows, columns, Math.min(10 * rows, Integer.MAX_VALUE));
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
  factory SparseCCDoubleMatrix2D.sized(int rows, int columns, [int nzmax = null]) {
    if (nzmax == null) {
      nzmax = 10 * rows;//Math.min(10 * rows, Integer.MAX_VALUE);
    }
    /*try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if (!"matrix too large".equals(exc.getMessage())) throw exc;
    }*/
    final dcs = cs_spalloc(rows, columns, nzmax, true, false);
    return new SparseCCDoubleMatrix2D(dcs);
  }

  /**
   * Constructs a matrix with given parameters. The arrays are not copied.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @param rowIndexes
   *            row indexes
   * @param columnPointers
   *            columns pointers
   * @param values
   *            numerical values
   */
  factory SparseCCDoubleMatrix2D.from(int rows, int columns, Int32List rowIndexes, Int32List columnPointers, Float64List values) {
    /*try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if (!"matrix too large".equals(exc.getMessage())) throw exc;
    }
    if (columnPointers.length != columns + 1) {
      throw new ArgumentError("columnsPointers.length != columns + 1");
    }*/
    Dcs dcs = new Dcs();
    dcs.m = rows;
    dcs.n = columns;
    dcs.i = rowIndexes;
    dcs.p = columnPointers;
    dcs.x = values;
    dcs.nz = -1; // column-compressed
    dcs.nzmax = values.length;
    return new SparseCCDoubleMatrix2D(dcs);
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
  factory SparseCCDoubleMatrix2D.value(int rows, int columns, Int32List rowIndexes, Int32List columnIndexes, double value, bool removeDuplicates, bool sortRowIndexes) {
    /*try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if (!"matrix too large".equals(exc.getMessage())) throw exc;
    }*/
    if (rowIndexes.length != columnIndexes.length) {
      throw new ArgumentError("rowIndexes.length != columnIndexes.length");
    }
    if (value == 0) {
      throw new ArgumentError("value cannot be 0");
    }

    int nz = Math.max(rowIndexes.length, 1);
    Dcs dcs = cs_spalloc(rows, columns, nz, true, false);
    Int32List w = new Int32List(columns);
    Int32List Cp = dcs.p;
    Int32List Ci = dcs.i;
    Float64List Cx = dcs.x;
    for (int k = 0; k < nz; k++) {
      w[columnIndexes[k]]++;
    }
    cs_cumsum(Cp, w, columns);
    int p;
    for (int k = 0; k < nz; k++) {
      Ci[p = w[columnIndexes[k]]++] = rowIndexes[k];
      if (Cx != null) Cx[p] = value;
    }
    if (removeDuplicates) {
      if (!cs_dupl(dcs)) { //remove duplicates
        throw new ArgumentError("Exception occured in cs_dupl()!");
      }
    }
    if (sortRowIndexes) {
      //sort row indexes
      dcs = cs_transpose(dcs, true);
      dcs = cs_transpose(dcs, true);
      if (dcs == null) {
        throw new ArgumentError("Exception occured in cs_transpose()!");
      }
      _rowIndexesSorted = true;
    }
    return new SparseCCDoubleMatrix2D(dcs);
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
  factory SparseCCDoubleMatrix2D.values(int rows, int columns, Int32List rowIndexes, Int32List columnIndexes, Float64List values, bool removeDuplicates, bool removeZeroes, bool sortRowIndexes) {
    /*try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if (!"matrix too large".equals(exc.getMessage())) throw exc;
    }
    if (rowIndexes.length != columnIndexes.length) {
      throw new ArgumentError("rowIndexes.length != columnIndexes.length");
    } else if (rowIndexes.length != values.length) {
      throw new ArgumentError("rowIndexes.length != values.length");
    }*/
    int nz = Math.max(rowIndexes.length, 1);
    Dcs dcs = cs_spalloc(rows, columns, nz, true, false);
    Int32List w = new Int32List(columns);
    Int32List Cp = dcs.p;
    Int32List Ci = dcs.i;
    Float64List Cx = dcs.x;
    for (int k = 0; k < nz; k++) w[columnIndexes[k]]++;
    cs_cumsum(Cp, w, columns);
    int p;
    for (int k = 0; k < nz; k++) {
      Ci[p = w[columnIndexes[k]]++] = rowIndexes[k];
      if (Cx != null) Cx[p] = values[k];
    }
    if (removeZeroes) {
      cs_dropzeros(dcs); //remove zeroes
    }
    if (removeDuplicates) {
      if (!cs_dupl(dcs)) { //remove duplicates
        throw new ArgumentError("Exception occured in cs_dupl()!");
      }
    }
    //sort row indexes
    if (sortRowIndexes) {
      dcs = cs_transpose(dcs, true);
      dcs = cs_transpose(dcs, true);
      if (dcs == null) {
        throw new ArgumentError("Exception occured in cs_transpose()!");
      }
      _rowIndexesSorted = true;
    }
    return new SparseCCDoubleMatrix2D(dcs);
  }

  DoubleMatrix2D assign(final func.DoubleFunction function) {
    if (function is DoubleMult) { // x[i] = mult*x[i]
      final double alpha = (function as DoubleMult).multiplicator;
      if (alpha == 1) return this;
      if (alpha == 0) return assignValue(0.0);
      if (alpha != alpha) {
        return assignValue(alpha); // the funny definition of isNaN(). This should better not happen.
      }

      final Float64List valuesE = _dcs.x;
      int nz = cardinality();
      for (int j = 0; j < nz; j++) {
        valuesE[j] *= alpha;
      }
    } else {
      forEachNonZero((int i, int j, double value) {
        return function(value);
      });
    }
    return this;
  }

  DoubleMatrix2D assignValue(double value) {
    if (value == 0) {
      _dcs.i.fillRange(0, _dcs.i.length, 0);
      _dcs.p.fillRange(0, _dcs.p.length, 0);
      _dcs.x.fillRange(0, _dcs.x.length, 0.0);
    } else {
      int nnz = cardinality();
      for (int i = 0; i < nnz; i++) {
        _dcs.x[i] = value;
      }
    }
    return this;
  }

  DoubleMatrix2D assignMatrix(DoubleMatrix2D source) {
    if (source == this) return this; // nothing to do
    checkShape(source);

    if (source is SparseCCDoubleMatrix2D) {
      SparseCCDoubleMatrix2D other = source;
      this._dcs.p.setAll(0, other.getColumnPointers());
      //System.arraycopy(other.getColumnPointers(), 0, this._dcs.p, 0, _columns + 1);
      int nzmax = other.getRowIndexes().length;
      if (_dcs.nzmax < nzmax) {
        _dcs.i = new Int32List(nzmax);
        _dcs.x = new Float64List(nzmax);
      }
      this._dcs.i.setAll(0, other.getRowIndexes());
      //System.arraycopy(other.getRowIndexes(), 0, this._dcs.i, 0, nzmax);
      this._dcs.x.setAll(0, other.getValues());
      //System.arraycopy(other.getValues(), 0, this._dcs.x, 0, nzmax);
      _rowIndexesSorted = other._rowIndexesSorted;
    } else if (source is SparseRCDoubleMatrix2D) {
      SparseRCDoubleMatrix2D other = source.getTranspose();
      this._dcs.p = other.getRowPointers();
      this._dcs.i = other.getColumnIndexes();
      this._dcs.x = other.getValues();
      this._dcs.nzmax = this._dcs.x.length;
      _rowIndexesSorted = true;
    } else {
      assignValue(0.0);
      source.forEachNonZero((int i, int j, double value) {
        setQuick(i, j, value);
        return value;
      });
    }
    return this;
  }

  DoubleMatrix2D assignFunc(final DoubleMatrix2D y, func.DoubleDoubleFunction function) {
    checkShape(y);

    if ((y is SparseCCDoubleMatrix2D) && (function == func.plus)) { // x[i] = x[i] + y[i]
      SparseCCDoubleMatrix2D yy = y;
      _dcs = cs_add(_dcs, yy._dcs, 1.0, 1.0);
      return this;
    }

    if (function is DoublePlusMultSecond) { // x[i] = x[i] + alpha*y[i]
      final double alpha = (function as DoublePlusMultSecond).multiplicator;
      if (alpha == 0) return this; // nothing to do
      y.forEachNonZero((int i, int j, double value) {
        setQuick(i, j, getQuick(i, j) + alpha * value);
        return value;
      });
      return this;
    }

    if (function is DoublePlusMultFirst) { // x[i] = alpha*x[i] + y[i]
      final double alpha = (function as DoublePlusMultFirst).multiplicator;
      if (alpha == 0) return assignMatrix(y);
      y.forEachNonZero((int i, int j, double value) {
        setQuick(i, j, alpha * getQuick(i, j) + value);
        return value;
      });
      return this;
    }

    if (function == func.mult) { // x[i] = x[i] * y[i]
      final Int32List rowIndexesA = _dcs.i;
      final Int32List columnPointersA = _dcs.p;
      final Float64List valuesA = _dcs.x;
      for (int j = _columns; --j >= 0; ) {
        int low = columnPointersA[j];
        for (int k = columnPointersA[j + 1]; --k >= low; ) {
          int i = rowIndexesA[k];
          valuesA[k] *= y.getQuick(i, j);
          if (valuesA[k] == 0) _remove(i, j);
        }
      }
      return this;
    }

    if (function == func.div) { // x[i] = x[i] / y[i]
      final Int32List rowIndexesA = _dcs.i;
      final Int32List columnPointersA = _dcs.p;
      final Float64List valuesA = _dcs.x;

      for (int j = _columns; --j >= 0; ) {
        int low = columnPointersA[j];
        for (int k = columnPointersA[j + 1]; --k >= low; ) {
          int i = rowIndexesA[k];
          valuesA[k] /= y.getQuick(i, j);
          if (valuesA[k] == 0) _remove(i, j);
        }
      }
      return this;
    }
    return super.assignFunc(y, function);
  }

  int cardinality() {
    return _dcs.p[_columns];
  }

  Dcs elements() {
    return _dcs;
  }

  DoubleMatrix2D forEachNonZero(final func.IntIntDoubleFunction function) {
    final Int32List rowIndexesA = _dcs.i;
    final Int32List columnPointersA = _dcs.p;
    final Float64List valuesA = _dcs.x;

    for (int j = _columns; --j >= 0; ) {
      int low = columnPointersA[j];
      for (int k = columnPointersA[j + 1]; --k >= low; ) {
        int i = rowIndexesA[k];
        double value = valuesA[k];
        double r = function(i, j, value);
        valuesA[k] = r;
      }
    }
    return this;
  }

  /**
   * Returns column pointers
   *
   * @return column pointers
   */
  Int32List getColumnPointers() {
    return _dcs.p;
  }

  /**
   * Returns a new matrix that has the same elements as this matrix, but is in
   * a dense form. This method creates a new object (not a view), so changes
   * in the returned matrix are NOT reflected in this matrix.
   *
   * @return this matrix in a dense form
   */
  DenseDoubleMatrix2D getDense() {
    final DenseDoubleMatrix2D dense = new DenseDoubleMatrix2D(_rows, _columns);
    forEachNonZero((int i, int j, double value) {
      dense.setQuick(i, j, getQuick(i, j));
      return value;
    });
    return dense;
  }

  double getQuick(int row, int column) {
    //        int k = cern.colt.Sorting.binarySearchFromTo(dcs.i, row, dcs.p[column], dcs.p[column + 1] - 1);
    int k = _searchFromTo(_dcs.i, row, _dcs.p[column], _dcs.p[column + 1] - 1);
    double v = 0.0;
    if (k >= 0) v = _dcs.x[k];
    return v;
  }

  /**
   * Returns a new matrix that has the same elements as this matrix, but is in
   * a row-compressed form. This method creates a new object (not a view), so
   * changes in the returned matrix are NOT reflected in this matrix.
   *
   * @return this matrix in a row-compressed form
   */
  SparseRCDoubleMatrix2D getRowCompressed() {
    Dcs dcst = cs_transpose(_dcs, true);
    SparseRCDoubleMatrix2D rc = new SparseRCDoubleMatrix2D.sized(_rows, _columns);
    rc._columnIndexes = dcst.i;
    rc._rowPointers = dcst.p;
    rc._values = dcst.x;
    rc._columnIndexesSorted = true;
    return rc;
  }

  /**
   * Returns row indexes;
   *
   * @return row indexes
   */
  Int32List getRowIndexes() {
    return _dcs.i;
  }

  /**
   * Returns a new matrix that is the transpose of this matrix. This method
   * creates a new object (not a view), so changes in the returned matrix are
   * NOT reflected in this matrix.
   *
   * @return the transpose of this matrix
   */
  SparseCCDoubleMatrix2D getTranspose() {
    Dcs dcst = cs_transpose(_dcs, true);
    SparseCCDoubleMatrix2D tr = new SparseCCDoubleMatrix2D.sized(_columns, _rows);
    tr._dcs = dcst;
    return tr;
  }

  /**
   * Returns numerical values
   *
   * @return numerical values
   */
  Float64List getValues() {
    return _dcs.x;
  }

  /**
   * Returns true if row indexes are sorted, false otherwise
   *
   * @return true if row indexes are sorted, false otherwise
   */
  bool hasRowIndexesSorted() {
    return _rowIndexesSorted;
  }

  DoubleMatrix2D like2D(int rows, int columns) {
    return new SparseCCDoubleMatrix2D.sized(rows, columns);
  }

  DoubleMatrix1D like1D(int size) {
    return new SparseDoubleMatrix1D(size);
  }

  void setQuick(int row, int column, double value) {
    //        int k = cern.colt.Sorting.binarySearchFromTo(dcs.i, row, dcs.p[column], dcs.p[column + 1] - 1);
    int k = _searchFromTo(_dcs.i, row, _dcs.p[column], _dcs.p[column + 1] - 1);

    if (k >= 0) { // found
      if (value == 0) _remove(column, k); else _dcs.x[k] = value;
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
    _dcs = cs_transpose(_dcs, true);
    _dcs = cs_transpose(_dcs, true);
    if (_dcs == null) {
      throw new ArgumentError("Exception occured in cs_transpose()!");
    }
    _rowIndexesSorted = true;
  }

  /**
   * Removes (sums) duplicate entries (if any}
   */
  void removeDuplicates() {
    if (!cs_dupl(_dcs)) { //remove duplicates
      throw new ArgumentError("Exception occured in cs_dupl()!");
    }
  }

  /**
   * Removes zero entries (if any)
   */
  void removeZeroes() {
    cs_dropzeros(_dcs); //remove zeroes
  }

  void trimToSize() {
    cs_sprealloc(_dcs, 0);
  }

  String toString() {
    StringBuffer builder = new StringBuffer();
    builder..write(_rows)..write(" x ")..write(_columns)..write(" sparse matrix, nnz = ")
      ..write(cardinality())..write('\n');
    for (int i = 0; i < _columns; i++) {
      int high = _dcs.p[i + 1];
      for (int j = _dcs.p[i]; j < high; j++) {
        builder..write('(')..write(_dcs.i[j])..write(',')..write(i)..write(')')
          ..write('\t')..write(_dcs.x[j])..write('\n');
      }
    }
    return builder.toString();
  }

  DoubleMatrix1D zMult(DoubleMatrix1D y, DoubleMatrix1D z, [final double alpha=1.0, final double beta=0.0, final bool transposeA=false]) {
    final int rowsA = transposeA ? _columns : _rows;
    final int columnsA = transposeA ? _rows : _columns;

    bool ignore = (z == null || transposeA);
    if (z == null) z = new DenseDoubleMatrix1D(rowsA);

    if (!(y is DenseDoubleMatrix1D && z is DenseDoubleMatrix1D)) {
      return super.zMult(y, z, alpha, beta, transposeA);
    }

    if (columnsA != y.size() || rowsA > z.size()) throw new ArgumentError("Incompatible args: " + ((transposeA ? viewDice() : this).toStringShort()) + ", " + y.toStringShort() + ", " + z.toStringShort());

    DenseDoubleMatrix1D zz = z as DenseDoubleMatrix1D;
    final Float64List elementsZ = zz._elements;
    final int strideZ = zz.stride();
    final int zeroZ = zz.index(0);

    DenseDoubleMatrix1D yy = y as DenseDoubleMatrix1D;
    final Float64List elementsY = yy._elements;
    final int strideY = yy.stride();
    final int zeroY = yy.index(0);

    final Int32List rowIndexesA = _dcs.i;
    final Int32List columnPointersA = _dcs.p;
    final Float64List valuesA = _dcs.x;

    int zidx = zeroZ;
    int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if (!transposeA) {
      if ((!ignore) && (beta / alpha != 1.0)) {
        z.assign(func.multiply(beta / alpha));
      }

      /*if ((nthreads > 1) && (cardinality() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
        nthreads = 2;
        List<Future> futures = new List<Future>(nthreads);
        final Float64List result = new Float64List(rowsA);
        int k = _columns ~/ nthreads;
        for (int j = 0; j < nthreads; j++) {
          final int firstColumn = j * k;
          final int lastColumn = (j == nthreads - 1) ? _columns : firstColumn + k;
          final int threadID = j;
          futures[j] = ConcurrencyUtils.submit(() {
            if (threadID == 0) {
              for (int i = firstColumn; i < lastColumn; i++) {
                int high = columnPointersA[i + 1];
                double yElem = elementsY[zeroY + strideY * i];
                for (int k = columnPointersA[i]; k < high; k++) {
                  int j = rowIndexesA[k];
                  elementsZ[zeroZ + strideZ * j] += valuesA[k] * yElem;
                }
              }
            } else {
              for (int i = firstColumn; i < lastColumn; i++) {
                int high = columnPointersA[i + 1];
                double yElem = elementsY[zeroY + strideY * i];
                for (int k = columnPointersA[i]; k < high; k++) {
                  int j = rowIndexesA[k];
                  result[j] += valuesA[k] * yElem;
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
          double yElem = elementsY[zeroY + strideY * i];
          for (int k = columnPointersA[i]; k < high; k++) {
            int j = rowIndexesA[k];
            elementsZ[zeroZ + strideZ * j] += valuesA[k] * yElem;
          }
        }
      //}
      if (alpha != 1.0) {
        z.assign(func.multiply(alpha));
      }
    } else {
      /*if ((nthreads > 1) && (cardinality() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
        List<Future> futures = new List<Future>(nthreads);
        int k = _columns ~/ nthreads;
        for (int j = 0; j < nthreads; j++) {
          final int firstColumn = j * k;
          final int lastColumn = (j == nthreads - 1) ? _columns : firstColumn + k;
          futures[j] = ConcurrencyUtils.submit(() {
            int zidx = zeroZ + firstColumn * strideZ;
            int k = _dcs.p[firstColumn];
            for (int i = firstColumn; i < lastColumn; i++) {
              double sum = 0;
              int high = _dcs.p[i + 1];
              for ( ; k + 10 < high; k += 10) {
                int ind = k + 9;
                sum += valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]];
              }
              for ( ; k < high; k++) {
                sum += valuesA[k] * elementsY[_dcs.i[k]];
              }
              elementsZ[zidx] = alpha * sum + beta * elementsZ[zidx];
              zidx += strideZ;
            }
          });
        }
        ConcurrencyUtils.waitForCompletion(futures);
      } else {*/
        int k = _dcs.p[0];
        for (int i = 0; i < _columns; i++) {
          double sum = 0.0;
          int high = _dcs.p[i + 1];
          for ( ; k + 10 < high; k += 10) {
            int ind = k + 9;
            sum += valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] + valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]];
          }
          for ( ; k < high; k++) {
            sum += valuesA[k] * elementsY[_dcs.i[k]];
          }
          elementsZ[zidx] = alpha * sum + beta * elementsZ[zidx];
          zidx += strideZ;
        }
      //}
    }
    return z;
  }

  DoubleMatrix2D zMult2D(DoubleMatrix2D B, DoubleMatrix2D C, [final double alpha=1.0, double beta=0.0, final bool transposeA=false, bool transposeB=false]) {
    int rowsA = _rows;
    int columnsA = _columns;
    if (transposeA) {
      rowsA = _columns;
      columnsA = _rows;
    }
    int rowsB = B.rows();
    int columnsB = B.columns();
    if (transposeB) {
      rowsB = B.columns();
      columnsB = B.rows();
    }
    int p = columnsB;
    bool ignore = (C == null);
    if (C == null) {
      if (B is SparseCCDoubleMatrix2D) {
        C = new SparseCCDoubleMatrix2D.sized(rowsA, p, (rowsA * p));
      } else {
        C = new DenseDoubleMatrix2D(rowsA, p);
      }
    }

    if (rowsB != columnsA) {
      throw new ArgumentError("Matrix2D inner dimensions must agree:" + toStringShort() + ", " + (transposeB ? B.viewDice() : B).toStringShort());
    }
    if (C.rows() != rowsA || C.columns() != p) {
      throw new ArgumentError("Incompatible result matrix: " + toStringShort() + ", " + (transposeB ? B.viewDice() : B).toStringShort() + ", " + C.toStringShort());
    }
    if (this == C || B == C) {
      throw new ArgumentError("Matrices must not be identical");
    }

    if (!ignore && beta != 1.0) {
      C.assign(func.multiply(beta));
    }

    if ((B is DenseDoubleMatrix2D) && (C is DenseDoubleMatrix2D)) {
      SparseCCDoubleMatrix2D AA;
      if (transposeA) {
        AA = getTranspose();
      } else {
        AA = this;
      }
      DenseDoubleMatrix2D BB;
      if (transposeB) {
        BB = B.viewDice() as DenseDoubleMatrix2D;
      } else {
        BB = B;
      }
      DenseDoubleMatrix2D CC = C;
      Int32List columnPointersA = AA._dcs.p;
      Int32List rowIndexesA = AA._dcs.i;
      Float64List valuesA = AA._dcs.x;

      int zeroB = BB.index(0, 0);
      int rowStrideB = BB.rowStride();
      int columnStrideB = BB.columnStride();
      Float64List elementsB = BB._elements;

      int zeroC = CC.index(0, 0);
      int rowStrideC = CC.rowStride();
      int columnStrideC = CC.columnStride();
      Float64List elementsC = CC._elements;

      for (int jj = 0; jj < columnsB; jj++) {
        for (int kk = 0; kk < columnsA; kk++) {
          int high = columnPointersA[kk + 1];
          double yElem = elementsB[zeroB + kk * rowStrideB + jj * columnStrideB];
          for (int ii = columnPointersA[kk]; ii < high; ii++) {
            int j = rowIndexesA[ii];
            elementsC[zeroC + j * rowStrideC + jj * columnStrideC] += valuesA[ii] * yElem;
          }
        }
      }
      if (alpha != 1.0) {
        C.assign(func.multiply(alpha));
      }

    } else if ((B is SparseCCDoubleMatrix2D) && (C is SparseCCDoubleMatrix2D)) {
      SparseCCDoubleMatrix2D AA;
      if (transposeA) {
        AA = getTranspose();
      } else {
        AA = this;
      }
      SparseCCDoubleMatrix2D BB = B;
      if (transposeB) {
        BB = BB.getTranspose();
      }
      SparseCCDoubleMatrix2D CC = C;
      CC._dcs = cs_multiply(AA._dcs, BB._dcs);
      if (CC._dcs == null) {
        throw new ArgumentError("Exception occured in cs_multiply()");
      }
      if (alpha != 1.0) {
        CC.assign(func.multiply(alpha));
      }
    } else {
      if (transposeB) {
        B = B.viewDice();
      }
      // cache views
      final List<DoubleMatrix1D> Brows = new List<DoubleMatrix1D>(columnsA);
      for (int i = columnsA; --i >= 0; ) Brows[i] = B.viewRow(i);
      final List<DoubleMatrix1D> Crows = new List<DoubleMatrix1D>(rowsA);
      for (int i = rowsA; --i >= 0; ) Crows[i] = C.viewRow(i);

      final DoublePlusMultSecond fun = new DoublePlusMultSecond.plusMult(0.0);

      final Int32List rowIndexesA = _dcs.i;
      final Int32List columnPointersA = _dcs.p;
      final Float64List valuesA = _dcs.x;
      for (int i = _columns; --i >= 0; ) {
        int low = columnPointersA[i];
        for (int k = columnPointersA[i + 1]; --k >= low; ) {
          int j = rowIndexesA[k];
          fun.multiplicator = valuesA[k] * alpha;
          if (!transposeA) {
            Crows[j].assignFunc(Brows[i], fun);
          } else {
            Crows[i].assignFunc(Brows[j], fun);
          }
        }
      }
    }
    return C;
  }

  DoubleMatrix2D _getContent() {
    return this;
  }

  void _insert(int row, int column, int index, double value) {
    List<int> rowIndexes = new List<int>.from(_dcs.i);
    //rowIndexes._setSizeRaw(_dcs.p[_columns]);
    List<double> values = new List<double>.from(_dcs.x);
    //values._setSizeRaw(_dcs.p[_columns]);
    rowIndexes.insert(index, row);
    values.insert(index, value);
    for (int i = _dcs.p.length; --i > column; ) {
      _dcs.p[i]++;
    }
    _dcs.i = new Int32List.fromList(rowIndexes);
    _dcs.x = new Float64List.fromList(values);
    _dcs.nzmax = rowIndexes.length;
  }

  void _remove(int column, int index) {
    List<int> rowIndexes = new List<int>.from(_dcs.i);
    List<double> values = new List<double>.from(_dcs.x);
    rowIndexes.removeAt(index);
    values.removeAt(index);
    for (int i = _dcs.p.length; --i > column; ) {
      _dcs.p[i]--;
    }
    _dcs.i = new Int32List.fromList(rowIndexes);
    _dcs.x = new Float64List.fromList(values);
    _dcs.nzmax = rowIndexes.length;
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

  Object clone() {
    return new SparseCCDoubleMatrix2D(_dcs);
  }
}
