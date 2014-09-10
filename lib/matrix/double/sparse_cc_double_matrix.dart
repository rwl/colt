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
 * {@link #assign(DoubleMatrix,func.DoubleDoubleFunction)}
 * does just that if you supply
 * {@link func#plusMultSecond} as argument.
 *
 *
 * @author Piotr Wendykier
 *
 */
class SparseCCDoubleMatrix extends WrapperDoubleMatrix {

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
  factory SparseCCDoubleMatrix.fromList(List<Float64List> values) {
    return new SparseCCDoubleMatrix(values.length, values[0].length)
      ..setAll2D(values);
  }

  /**
   * Constructs a matrix with a given internal storage.
   *
   * @param dcs
   *            internal storage.
   */
  SparseCCDoubleMatrix._internal(Dcs dcs) : super(null) {
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
  /*SparseCCDoubleMatrix(int rows, int columns) {
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
  SparseCCDoubleMatrix(int rows, int columns, [int nzmax = null]) : super(null) {
    if (nzmax == null) {
      nzmax = Math.min(10 * rows, MAX_INT);
    }
    try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if ("matrix too large" != exc.message) {
        throw exc;
      }
    }
    _dcs = cs_spalloc(rows, columns, nzmax, true, false);
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
  factory SparseCCDoubleMatrix.withPointers(int rows, int columns, Int32List rowIndexes, Int32List columnPointers, Float64List values) {
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
    return new SparseCCDoubleMatrix._internal(dcs);
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
  factory SparseCCDoubleMatrix.withValue(int rows, int columns, Int32List rowIndexes, Int32List columnIndexes, double value, bool removeDuplicates, bool sortRowIndexes) {
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
    bool rowIndexesSorted = false;
    if (sortRowIndexes) {
      //sort row indexes
      dcs = cs_transpose(dcs, true);
      dcs = cs_transpose(dcs, true);
      if (dcs == null) {
        throw new ArgumentError("Exception occured in cs_transpose()!");
      }
      rowIndexesSorted = true;
    }
    return new SparseCCDoubleMatrix._internal(dcs)
      .._rowIndexesSorted = rowIndexesSorted;
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
  factory SparseCCDoubleMatrix.withValues(int rows, int columns, Int32List rowIndexes, Int32List columnIndexes, Float64List values, [bool removeDuplicates=false, bool removeZeroes=false, bool sortRowIndexes=false]) {
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
    bool rowIndexesSorted = false;
    if (sortRowIndexes) {
      dcs = cs_transpose(dcs, true);
      dcs = cs_transpose(dcs, true);
      if (dcs == null) {
        throw new ArgumentError("Exception occured in cs_transpose()!");
      }
      rowIndexesSorted = true;
    }
    return new SparseCCDoubleMatrix._internal(dcs)
      .._rowIndexesSorted = rowIndexesSorted;
  }

  void forEach(final func.DoubleFunction function) {
    if (function is DoubleMult) { // x[i] = mult*x[i]
      final double alpha = (function as DoubleMult).multiplicator;
      if (alpha == 1) {
        return;
      }
      if (alpha == 0) {
        fill(0.0);
        return;
      }
      if (alpha != alpha) {
        fill(alpha); // the funny definition of isNaN(). This should better not happen.
        return;
      }

      final Float64List valuesE = _dcs.x;
      int nz = cardinality;
      for (int j = 0; j < nz; j++) {
        valuesE[j] *= alpha;
      }
    } else {
      forEachNonZero((int i, int j, double value) {
        return function(value);
      });
    }
  }

  void fill(double value) {
    if (value == 0) {
      _dcs.i.fillRange(0, _dcs.i.length, 0);
      _dcs.p.fillRange(0, _dcs.p.length, 0);
      _dcs.x.fillRange(0, _dcs.x.length, 0.0);
    } else {
      int nnz = cardinality;
      for (int i = 0; i < nnz; i++) {
        _dcs.x[i] = value;
      }
    }
  }

  void copyFrom(AbstractDoubleMatrix source) {
    if (source == this) {
      return; // nothing to do
    }
    checkShape(source);

    if (source is SparseCCDoubleMatrix) {
      SparseCCDoubleMatrix other = source;
      this._dcs.p.setAll(0, other.columnPointers);
      //System.arraycopy(other.getColumnPointers(), 0, this._dcs.p, 0, _columns + 1);
      int nzmax = other.rowIndexes.length;
      if (_dcs.nzmax < nzmax) {
        _dcs.i = new Int32List(nzmax);
        _dcs.x = new Float64List(nzmax);
      }
      this._dcs.i.setAll(0, other.rowIndexes);
      //System.arraycopy(other.getRowIndexes(), 0, this._dcs.i, 0, nzmax);
      this._dcs.x.setAll(0, other.values());
      //System.arraycopy(other.getValues(), 0, this._dcs.x, 0, nzmax);
      _rowIndexesSorted = other._rowIndexesSorted;
    } else if (source is SparseRCDoubleMatrix) {
      SparseRCDoubleMatrix other = source.getTranspose();
      this._dcs.p = other.getRowPointers();
      this._dcs.i = other.getColumnIndexes();
      this._dcs.x = other.getValues();
      this._dcs.nzmax = this._dcs.x.length;
      _rowIndexesSorted = true;
    } else {
      fill(0.0);
      source.forEachNonZero((int i, int j, double value) {
        set(i, j, value);
        return value;
      });
    }
  }

  void forEachMatrix(final AbstractDoubleMatrix y, func.DoubleDoubleFunction function) {
    checkShape(y);

    if ((y is SparseCCDoubleMatrix) && (function == func.plus)) { // x[i] = x[i] + y[i]
      SparseCCDoubleMatrix yy = y;
      _dcs = cs_add(_dcs, yy._dcs, 1.0, 1.0);
      return;
    }

    if (function is DoublePlusMultSecond) { // x[i] = x[i] + alpha*y[i]
      final double alpha = (function as DoublePlusMultSecond).multiplicator;
      if (alpha == 0) {
        return; // nothing to do
      }
      y.forEachNonZero((int i, int j, double value) {
        set(i, j, get(i, j) + alpha * value);
        return value;
      });
      return;
    }

    if (function is DoublePlusMultFirst) { // x[i] = alpha*x[i] + y[i]
      final double alpha = (function as DoublePlusMultFirst).multiplicator;
      if (alpha == 0) {
        copyFrom(y);
        return;
      }
      y.forEachNonZero((int i, int j, double value) {
        set(i, j, alpha * get(i, j) + value);
        return value;
      });
      return;
    }

    if (function == func.mult) { // x[i] = x[i] * y[i]
      final Int32List rowIndexesA = _dcs.i;
      final Int32List columnPointersA = _dcs.p;
      final Float64List valuesA = _dcs.x;
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

    if (function == func.div) { // x[i] = x[i] / y[i]
      final Int32List rowIndexesA = _dcs.i;
      final Int32List columnPointersA = _dcs.p;
      final Float64List valuesA = _dcs.x;

      for (int j = _columns; --j >= 0; ) {
        int low = columnPointersA[j];
        for (int k = columnPointersA[j + 1]; --k >= low; ) {
          int i = rowIndexesA[k];
          valuesA[k] /= y.get(i, j);
          if (valuesA[k] == 0) _remove(i, j);
        }
      }
      return;
    }
    super.forEachMatrix(y, function);
  }

  int get cardinality {
    return _dcs.p[_columns];
  }

  Dcs elements() {
    return _dcs;
  }

  void forEachNonZero(final func.IntIntDoubleFunction function) {
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
  }

  /**
   * Returns column pointers
   *
   * @return column pointers
   */
  Int32List get columnPointers {
    return _dcs.p;
  }

  /**
   * Returns a new matrix that has the same elements as this matrix, but is in
   * a dense form. This method creates a new object (not a view), so changes
   * in the returned matrix are NOT reflected in this matrix.
   *
   * @return this matrix in a dense form
   */
  DoubleMatrix dense() {
    final DoubleMatrix dense = new DoubleMatrix(_rows, _columns);
    forEachNonZero((int i, int j, double value) {
      dense.set(i, j, get(i, j));
      return value;
    });
    return dense;
  }

  double get(int row, int column) {
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
  SparseRCDoubleMatrix rowCompressed() {
    Dcs dcst = cs_transpose(_dcs, true);
    SparseRCDoubleMatrix rc = new SparseRCDoubleMatrix.sized(_rows, _columns);
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
  Int32List get rowIndexes {
    return _dcs.i;
  }

  /**
   * Returns a new matrix that is the transpose of this matrix. This method
   * creates a new object (not a view), so changes in the returned matrix are
   * NOT reflected in this matrix.
   *
   * @return the transpose of this matrix
   */
  SparseCCDoubleMatrix transpose() {
    Dcs dcst = cs_transpose(_dcs, true);
    SparseCCDoubleMatrix tr = new SparseCCDoubleMatrix(_columns, _rows);
    tr._dcs = dcst;
    return tr;
  }

  /**
   * Returns numerical values
   *
   * @return numerical values
   */
  Float64List values() {
    return _dcs.x;
  }

  /**
   * Returns true if row indexes are sorted, false otherwise
   *
   * @return true if row indexes are sorted, false otherwise
   */
  bool get rowIndexesSorted {
    return _rowIndexesSorted;
  }

  AbstractDoubleMatrix like2D(int rows, int columns) {
    return new SparseCCDoubleMatrix(rows, columns);
  }

  AbstractDoubleVector like1D(int size) {
    return new SparseDoubleVector(size);
  }

  void set(int row, int column, double value) {
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
      ..write(cardinality)..write('\n');
    for (int i = 0; i < _columns; i++) {
      int high = _dcs.p[i + 1];
      for (int j = _dcs.p[i]; j < high; j++) {
        builder..write('(')..write(_dcs.i[j])..write(',')..write(i)..write(')')
          ..write('\t')..write(_dcs.x[j])..write('\n');
      }
    }
    return builder.toString();
  }

  AbstractDoubleVector mult(AbstractDoubleVector y, [AbstractDoubleVector z = null, final double alpha=1.0, final double beta=0.0, final bool transposeA=false]) {
    final int rowsA = transposeA ? _columns : _rows;
    final int columnsA = transposeA ? _rows : _columns;

    bool ignore = (z == null || transposeA);
    if (z == null) z = new DoubleVector(rowsA);

    if (!(y is DoubleVector && z is DoubleVector)) {
      return super.mult(y, z, alpha, beta, transposeA);
    }

    if (columnsA != y.length || rowsA > z.length) throw new ArgumentError("Incompatible args: " + ((transposeA ? dice() : this).toStringShort()) + ", " + y.toStringShort() + ", " + z.toStringShort());

    DoubleVector zz = z as DoubleVector;
    final Float64List elementsZ = zz._elements;
    final int strideZ = zz.stride();
    final int zeroZ = zz.index(0);

    DoubleVector yy = y as DoubleVector;
    final Float64List elementsY = yy._elements;
    final int strideY = yy.stride();
    final int zeroY = yy.index(0);

    final Int32List rowIndexesA = _dcs.i;
    final Int32List columnPointersA = _dcs.p;
    final Float64List valuesA = _dcs.x;

    int zidx = zeroZ;
    //int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if (!transposeA) {
      if ((!ignore) && (beta / alpha != 1.0)) {
        z.forEach(func.multiply(beta / alpha));
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
        z.forEach(func.multiply(alpha));
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

  AbstractDoubleMatrix multiply(AbstractDoubleMatrix B, [AbstractDoubleMatrix C = null, final double alpha=1.0, double beta=0.0, final bool transposeA=false, bool transposeB=false]) {
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
      if (B is SparseCCDoubleMatrix) {
        C = new SparseCCDoubleMatrix(rowsA, p, (rowsA * p));
      } else {
        C = new DoubleMatrix(rowsA, p);
      }
    }

    if (rowsB != columnsA) {
      throw new ArgumentError("Matrix inner dimensions must agree:" + toStringShort() + ", " + (transposeB ? B.dice() : B).toStringShort());
    }
    if (C.rows != rowsA || C.columns != p) {
      throw new ArgumentError("Incompatible result matrix: " + toStringShort() + ", " + (transposeB ? B.dice() : B).toStringShort() + ", " + C.toStringShort());
    }
    if (this == C || B == C) {
      throw new ArgumentError("Matrices must not be identical");
    }

    if (!ignore && beta != 1.0) {
      C.forEach(func.multiply(beta));
    }

    if ((B is DoubleMatrix) && (C is DoubleMatrix)) {
      SparseCCDoubleMatrix AA;
      if (transposeA) {
        AA = transpose();
      } else {
        AA = this;
      }
      DoubleMatrix BB;
      if (transposeB) {
        BB = B.dice() as DoubleMatrix;
      } else {
        BB = B;
      }
      DoubleMatrix CC = C;
      Int32List columnPointersA = AA._dcs.p;
      Int32List rowIndexesA = AA._dcs.i;
      Float64List valuesA = AA._dcs.x;

      int zeroB = BB.index(0, 0);
      int rowStrideB = BB.rowStride;
      int columnStrideB = BB.columnStride;
      Float64List elementsB = BB._elements;

      int zeroC = CC.index(0, 0);
      int rowStrideC = CC.rowStride;
      int columnStrideC = CC.columnStride;
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
        C.forEach(func.multiply(alpha));
      }

    } else if ((B is SparseCCDoubleMatrix) && (C is SparseCCDoubleMatrix)) {
      SparseCCDoubleMatrix AA;
      if (transposeA) {
        AA = transpose();
      } else {
        AA = this;
      }
      SparseCCDoubleMatrix BB = B;
      if (transposeB) {
        BB = BB.transpose();
      }
      SparseCCDoubleMatrix CC = C;
      CC._dcs = cs_multiply(AA._dcs, BB._dcs);
      if (CC._dcs == null) {
        throw new ArgumentError("Exception occured in cs_multiply()");
      }
      if (alpha != 1.0) {
        CC.forEach(func.multiply(alpha));
      }
    } else {
      if (transposeB) {
        B = B.dice();
      }
      // cache views
      final List<AbstractDoubleVector> Brows = new List<AbstractDoubleVector>(columnsA);
      for (int i = columnsA; --i >= 0; ) Brows[i] = B.row(i);
      final List<AbstractDoubleVector> Crows = new List<AbstractDoubleVector>(rowsA);
      for (int i = rowsA; --i >= 0; ) Crows[i] = C.row(i);

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
            Crows[j].forEachWith(Brows[i], fun);
          } else {
            Crows[i].forEachWith(Brows[j], fun);
          }
        }
      }
    }
    return C;
  }

  AbstractDoubleMatrix _getContent() {
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
    return new SparseCCDoubleMatrix._internal(_dcs);
  }
}
