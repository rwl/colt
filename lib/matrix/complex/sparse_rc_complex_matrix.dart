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
 * Sparse row-compressed 2-d matrix holding <tt>complex</tt> elements. First see
 * the <a href="package-summary.html">package summary</a> and javadoc <a
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
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
class SparseRCComplexMatrix extends WrapperComplexMatrix {

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

  Float64List _values;

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
  factory SparseRCComplexMatrix.fromList(List<Float64List> values) {
    return new SparseRCComplexMatrix(values.length, values.length == 0 ? 0 : values[0].length)
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
  /*SparseRCComplexMatrix(int rows, int columns) {
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
   *             if
   *             <tt>rows<0 || columns<0 || (double)columns*rows > Integer.MAX_VALUE</tt>
   *             .
   */
  factory SparseRCComplexMatrix(int rows, int columns, [int nzmax = null]) {
    if (nzmax == null) {
      nzmax = 10 * rows;//Math.min(10 * rows, Integer.MAX_VALUE);
    }
    /*try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if ("matrix too large" != exc.message) {
        throw exc;
      }
    }*/
    final columnIndexes = new Int32List(nzmax);
    final values = new Float64List(2 * nzmax);
    final rowPointers = new Int32List(rows + 1);
    return new SparseRCComplexMatrix._internal(rows, columns, rowPointers, columnIndexes, values);
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
   * @param re
   *            real part of numerical value
   * @param im
   *            imaginary part of numerical value
   *
   * @param removeDuplicates
   *            if true, then duplicates (if any) are removed
   */
  factory SparseRCComplexMatrix.withValue(int rows, int columns, Int32List rowIndexes, Int32List columnIndexes, double re, double im, bool removeDuplicates) {
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
    if (re == 0 && im == 0) {
      throw new ArgumentError("value cannot be 0");
    }

    int nz = Math.max(rowIndexes.length, 1);
    final _columnIndexes = new Int32List(nz);
    final _values = new Float64List(2 * nz);
    final _rowPointers = new Int32List(rows + 1);
    Int32List w = new Int32List(rows);
    int r;
    for (int k = 0; k < nz; k++) {
      w[rowIndexes[k]]++;
    }
    _cumsum(_rowPointers, w, rows);
    for (int k = 0; k < nz; k++) {
      _columnIndexes[r = w[rowIndexes[k]]++] = columnIndexes[k];
      _values[2 * r] = re;
      _values[2 * r + 1] = im;
    }
    final m = new SparseRCComplexMatrix._internal(rows, columns, _rowPointers, _columnIndexes, _values);
    if (removeDuplicates) {
      m.removeDuplicates();
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
   */
  factory SparseRCComplexMatrix.withValues(int rows, int columns, Int32List rowIndexes, Int32List columnIndexes, Float64List values, {bool removeDuplicates: false, bool removeZeroes: false}) {
    /*try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if ("matrix too large" != exc.message) {
        throw exc;
      }
    }*/
    if (rowIndexes.length != columnIndexes.length) {
      throw new ArgumentError("rowIndexes.length != columnIndexes.length");
    } else if (2 * rowIndexes.length != values.length) {
      throw new ArgumentError("2 * rowIndexes.length != values.length");
    }
    int nz = Math.max(rowIndexes.length, 1);
    final _columnIndexes = new Int32List(nz);
    final _values = new Float64List(2 * nz);
    final _rowPointers = new Int32List(rows + 1);
    Int32List w = new Int32List(rows);
    int r;
    for (int k = 0; k < nz; k++) {
      w[rowIndexes[k]]++;
    }
    _cumsum(_rowPointers, w, rows);
    for (int k = 0; k < nz; k++) {
      _columnIndexes[r = w[rowIndexes[k]]++] = columnIndexes[k];
      _values[2 * r] = values[2 * k];
      _values[2 * r + 1] = values[2 * k + 1];
    }
    final m = new SparseRCComplexMatrix._internal(rows, columns, _rowPointers, _columnIndexes, _values);
    if (removeZeroes) {
      m.removeZeroes();
    }
    if (removeDuplicates) {
      m.removeDuplicates();
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
  SparseRCComplexMatrix._internal(int rows, int columns, Int32List rowPointers, Int32List columnIndexes, Float64List values) : super(null) {
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
    if (2 * columnIndexes.length != values.length) {
      throw new ArgumentError("2 * columnIndexes.length != values.length");
    }
    this._rowPointers = rowPointers;
    this._columnIndexes = columnIndexes;
    this._values = values;
  }

  void assign(final cfunc.ComplexComplexFunction function) {
    if (function is cfunc.ComplexMult) { // x[i] = mult*x[i]
      final Float64List alpha = (function as cfunc.ComplexMult).multiplicator;
      if (alpha[0] == 1 && alpha[1] == 0) {
        return;
      }
      if (alpha[0] == 0 && alpha[1] == 0) {
        fill(alpha[0], alpha[1]);
        return;
      }
      if (alpha[0] != alpha[0] || alpha[1] != alpha[1]) {
        fill(alpha[0], alpha[1]); // the funny definition of isNaN(). This should better not happen.
        return;
      }
      int nz = cardinality;
      Float64List elem = new Float64List(2);
      for (int j = 0; j < nz; j++) {
        elem[0] = _values[2 * j];
        elem[1] = _values[2 * j + 1];
        elem = Complex.multiply(elem, alpha);
        _values[2 * j] = elem[0];
        _values[2 * j + 1] = elem[1];
      }
    } else {
      forEachNonZero((int i, int j, Float64List value) {
        return function(value);
      });
    }
  }

  void fill(double re, double im) {
    if (re == 0 && im == 0) {
      //Arrays.fill(_rowPointers, 0);
      //Arrays.fill(_columnIndexes, 0);
      //Arrays.fill(_values, 0);
      _rowPointers.fillRange(0, _rowPointers.length, 0);
      _columnIndexes.fillRange(0, _columnIndexes.length, 0);
      _values.fillRange(0, _values.length, 0.0);
    } else {
      int nnz = cardinality;
      for (int i = 0; i < nnz; i++) {
        _values[2 * i] = re;
        _values[2 * i + 1] = im;
      }
    }
  }

  void copyFrom(AbstractComplexMatrix source) {
    if (source == this) {
      return; // nothing to do
    }
    checkShape(source);

    if (source is SparseRCComplexMatrix) {
      SparseRCComplexMatrix other = source;
      //System.arraycopy(other._rowPointers, 0, _rowPointers, 0, _rows + 1);
      _rowPointers.setAll(0, other._rowPointers);
      int nzmax = other._columnIndexes.length;
      if (_columnIndexes.length < nzmax) {
        _columnIndexes = new Int32List(nzmax);
        _values = new Float64List(2 * nzmax);
      }
      //System.arraycopy(other._columnIndexes, 0, _columnIndexes, 0, nzmax);
      _columnIndexes.setAll(0, other._columnIndexes);
      //System.arraycopy(other._values, 0, _values, 0, other._values.length);
      _values.setAll(0, other._values);
    } else if (source is SparseCCComplexMatrix) {
      SparseCCComplexMatrix other = source.conjugateTranspose();
      _rowPointers = other.columnPointers;
      _columnIndexes = other.rowIndexes;
      _values = other.values;
    } else {
      fill(0.0, 0.0);
      source.forEachNonZero((int i, int j, Float64List value) {
        set(i, j, value);
        return value;
      });
    }
  }

  void forEachMatrix(final AbstractComplexMatrix y, cfunc.ComplexComplexComplexFunction function) {
    checkShape(y);
    if ((y is SparseRCComplexMatrix) && (function == cfunc.plus)) { // x[i] = x[i] + y[i]
      SparseRCComplexMatrix yy = y;

      final Int32List rowPointersY = yy._rowPointers;
      final Int32List columnIndexesY = yy._columnIndexes;
      final Float64List valuesY = yy._values;

      final Int32List rowPointersC = new Int32List(_rows + 1);
      int cnz = Math.max(_columnIndexes.length, Math.min(MAX_INT, _rowPointers[_rows] + rowPointersY[_rows]));
      final Int32List columnIndexesC = new Int32List(cnz);
      final Float64List valuesC = new Float64List(2 * cnz);
      int nrow = _rows;
      int ncol = _columns;
      int nzmax = cnz;
      if (function == cfunc.plus) { // x[i] = x[i] + y[i]
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
              valuesC[2 * kc] = _values[2 * ka] + valuesY[2 * kb];
              valuesC[2 * kc + 1] = _values[2 * ka + 1] + valuesY[2 * kb + 1];
              columnIndexesC[kc] = j1;
              ka++;
              kb++;
              kc++;
            } else if (j1 < j2) {
              columnIndexesC[kc] = j1;
              valuesC[2 * kc] = _values[2 * ka];
              valuesC[2 * kc + 1] = _values[2 * ka + 1];
              ka++;
              kc++;
            } else if (j1 > j2) {
              columnIndexesC[kc] = j2;
              valuesC[2 * kc] = valuesY[2 * kb];
              valuesC[2 * kc + 1] = valuesY[2 * kb + 1];
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

    if (function is cfunc.ComplexPlusMultSecond) { // x[i] = x[i] + alpha*y[i]
      final Float64List alpha = (function as cfunc.ComplexPlusMultSecond).multiplicator;
      if (alpha[0] == 0 && alpha[1] == 0) {
        return; // nothing to do
      }
      y.forEachNonZero((int i, int j, Float64List value) {
        set(i, j, Complex.plus(get(i, j), Complex.multiply(alpha, value)));
        return value;
      });
      return;
    }

    if (function is cfunc.ComplexPlusMultFirst) { // x[i] = alpha*x[i] + y[i]
      final Float64List alpha = (function as cfunc.ComplexPlusMultFirst).multiplicator;
      if (alpha[0] == 0 && alpha[1] == 0) {
        copyFrom(y);
        return;
      }
      y.forEachNonZero((int i, int j, Float64List value) {
        set(i, j, Complex.plus(Complex.multiply(alpha, get(i, j)), value));
        return value;
      });
      return;
    }

    if (function == cfunc.mult) { // x[i] = x[i] * y[i]
      Float64List elem = new Float64List(2);
      for (int i = 0; i < _rows; i++) {
        int high = _rowPointers[i + 1];
        for (int k = _rowPointers[i]; k < high; k++) {
          int j = _columnIndexes[k];
          elem[0] = _values[2 * k];
          elem[1] = _values[2 * k + 1];
          elem = Complex.multiply(elem, y.get(i, j));
          _values[2 * k] = elem[0];
          _values[2 * k + 1] = elem[1];
          if (_values[2 * k] == 0 && _values[2 * k + 1] == 0) _remove(i, j);
        }
      }
      return;
    }

    if (function == cfunc.div) { // x[i] = x[i] / y[i]
      Float64List elem = new Float64List(2);
      for (int i = 0; i < _rows; i++) {
        int high = _rowPointers[i + 1];
        for (int k = _rowPointers[i]; k < high; k++) {
          int j = _columnIndexes[k];
          elem[0] = _values[2 * k];
          elem[1] = _values[2 * k + 1];
          elem = Complex.div_(elem, y.get(i, j));
          _values[2 * k] = elem[0];
          _values[2 * k + 1] = elem[1];
          if (_values[2 * k] == 0 && _values[2 * k + 1] == 0) _remove(i, j);
        }
      }
      return;
    }
    super.forEachWith(y, function);
  }

  int get cardinality {
    return _rowPointers[_rows];
  }

  void forEachNonZero(final cfunc.IntIntComplexFunction function) {
    Float64List value = new Float64List(2);
    for (int i = 0; i < _rows; i++) {
      int high = _rowPointers[i + 1];
      for (int k = _rowPointers[i]; k < high; k++) {
        int j = _columnIndexes[k];
        value[0] = _values[2 * k];
        value[1] = _values[2 * k + 1];
        Float64List r = function(i, j, value);
        if (r[0] != value[0] || r[1] != value[1]) {
          _values[2 * k] = r[0];
          _values[2 * k + 1] = r[1];
        }
      }
    }
  }

  /**
   * Returns a new matrix that has the same elements as this matrix, but is in
   * a column-compressed form. This method creates a new object (not a view),
   * so changes in the returned matrix are NOT reflected in this matrix.
   *
   * @return this matrix in a column-compressed form
   */
  SparseCCComplexMatrix columnCompressed() {
    SparseRCComplexMatrix tr = conjugateTranspose();
    SparseCCComplexMatrix cc = new SparseCCComplexMatrix(_rows, _columns);
    cc._rowIndexes = tr._columnIndexes;
    cc._columnPointers = tr._rowPointers;
    cc._values = tr._values;
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
  ComplexMatrix dense() {
    final ComplexMatrix dense = new ComplexMatrix(_rows, _columns);
    forEachNonZero((int i, int j, Float64List value) {
      dense.set(i, j, get(i, j));
      return value;
    });
    return dense;
  }

  Float64List get(int row, int column) {
    //int k = cern.colt.Sorting.binarySearchFromTo(columnIndexes, column, rowPointers[row], rowPointers[row + 1] - 1);
    int k = _searchFromTo(_columnIndexes, column, _rowPointers[row], _rowPointers[row + 1] - 1);

    Float64List v = new Float64List(2);
    if (k >= 0) {
      v[0] = _values[2 * k];
      v[1] = _values[2 * k + 1];
    }
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
  SparseRCComplexMatrix transpose() {
    int nnz = _rowPointers[_rows];
    Int32List w = new Int32List(_columns);
    Int32List rowPointersT = new Int32List(_columns + 1);
    Int32List columnIndexesT = new Int32List(nnz);
    Float64List valuesT = new Float64List(2 * nnz);

    for (int p = 0; p < nnz; p++) {
      w[_columnIndexes[p]]++;
    }
    _cumsum(rowPointersT, w, _columns);
    int q;
    for (int j = 0; j < _rows; j++) {
      int high = _rowPointers[j + 1];
      for (int p = _rowPointers[j]; p < high; p++) {
        columnIndexesT[q = w[_columnIndexes[p]]++] = j;
        valuesT[2 * q] = _values[2 * p];
        valuesT[2 * q + 1] = _values[2 * p + 1];
      }
    }
    SparseRCComplexMatrix T = new SparseRCComplexMatrix(_columns, _rows);
    T._rowPointers = rowPointersT;
    T._columnIndexes = columnIndexesT;
    T._values = valuesT;
    return T;
  }

  /**
   * Returns a new matrix that is the conjugate transpose of this matrix.
   * This method creates a new object (not a view), so changes in the
   * returned matrix are NOT reflected in this matrix.
   *
   * @return the conjugate transpose of this matrix
   */
  SparseRCComplexMatrix conjugateTranspose() {
    int nnz = _rowPointers[_rows];
    Int32List w = new Int32List(_columns);
    Int32List rowPointersT = new Int32List(_columns + 1);
    Int32List columnIndexesT = new Int32List(nnz);
    Float64List valuesT = new Float64List(2 * nnz);

    for (int p = 0; p < nnz; p++) {
      w[_columnIndexes[p]]++;
    }
    _cumsum(rowPointersT, w, _columns);
    int q;
    for (int j = 0; j < _rows; j++) {
      int high = _rowPointers[j + 1];
      for (int p = _rowPointers[j]; p < high; p++) {
        columnIndexesT[q = w[_columnIndexes[p]]++] = j;
        valuesT[2 * q] = _values[2 * p];
        valuesT[2 * q + 1] = -_values[2 * p + 1];
      }
    }
    SparseRCComplexMatrix T = new SparseRCComplexMatrix(_columns, _rows);
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
  Float64List get values {
    return _values;
  }

  AbstractComplexMatrix like2D(int rows, int columns) {
    return new SparseRCComplexMatrix(rows, columns);
  }

  AbstractComplexVector like1D(int size) {
    return new SparseComplexVector(size);
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
          _values[2 * w[i]] += _values[2 * p];
          /* A(i,j) is a duplicate */
          _values[2 * w[i] + 1] += _values[2 * p + 1];
        } else {
          w[i] = nz;
          /* record where column i occurs */
          _columnIndexes[nz] = i;
          /* keep A(i,j) */
          _values[2 * nz] = _values[2 * p];
          _values[2 * nz + 1] = _values[2 * p + 1];
          nz++;
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
    double eps = Math.pow(2, -52);
    Float64List elem = new Float64List(2);
    for (int j = 0; j < _rows; j++) {
      int p = _rowPointers[j];
      /* get current location of row j */
      _rowPointers[j] = nz;
      /* record new location of row j */
      for ( ; p < _rowPointers[j + 1]; p++) {
        elem[0] = _values[2 * p];
        elem[1] = _values[2 * p + 1];
        if (Complex.abs(elem) > eps) {
          _values[2 * nz] = _values[2 * p];
          /* keep A(i,j) */
          _values[2 * nz + 1] = _values[2 * p + 1];
          _columnIndexes[nz++] = _columnIndexes[p];
        }
      }
    }
    _rowPointers[_rows] = nz;
    /* finalize A */
  }

  void set(int row, int column, Float64List value) {
    //int k = cern.colt.Sorting.binarySearchFromTo(columnIndexes, column, rowPointers[row], rowPointers[row + 1] - 1);
    int k = _searchFromTo(_columnIndexes, column, _rowPointers[row], _rowPointers[row + 1] - 1);

    if (k >= 0) { // found
      if (value[0] == 0 && value[1] == 0) _remove(row, k); else {
        _values[2 * k] = value[0];
        _values[2 * k + 1] = value[1];
      }
      return;
    }

    if (value[0] != 0 || value[1] != 0) {
      k = -k - 1;
      _insert(row, column, k, value);
    }
  }

  void setParts(int row, int column, double re, double im) {
    //int k = cern.colt.Sorting.binarySearchFromTo(columnIndexes, column, rowPointers[row], rowPointers[row + 1] - 1);
    int k = _searchFromTo(_columnIndexes, column, _rowPointers[row], _rowPointers[row + 1] - 1);

    if (k >= 0) { // found
      if (re == 0 && im == 0) _remove(row, k); else {
        _values[2 * k] = re;
        _values[2 * k + 1] = im;
      }
      return;
    }

    if (re != 0 || im != 0) {
      k = -k - 1;
      _insertParts(row, column, k, re, im);
    }
  }

  String toString() {
    StringBuffer builder = new StringBuffer();
    builder..write(_rows)..write(" x ")..write(_columns)..write(" sparse matrix, nnz = ")..write(cardinality)..write('\n');
    for (int i = 0; i < _rows; i++) {
      int high = _rowPointers[i + 1];
      for (int j = _rowPointers[i]; j < high; j++) {
        if (_values[2 * j + 1] > 0) {
          builder..write('(')..write(i)..write(',')..write(_columnIndexes[j])..write(')')..write('\t')..write(_values[2 * j])..write('+')..write(_values[2 * j + 1])..write('i')..write('\n');
        } else if (_values[2 * j + 1] == 0) {
          builder..write('(')..write(i)..write(',')..write(_columnIndexes[j])..write(')')..write('\t')..write(_values[2 * j])..write('\n');
        } else {
          builder..write('(')..write(i)..write(',')..write(_columnIndexes[j])..write(')')..write('\t')..write(_values[2 * j])..write('-')..write(_values[2 * j + 1])..write('i')..write('\n');
        }
      }
    }
    return builder.toString();
  }

  void trimToSize() {
    _realloc(0);
  }

  AbstractComplexVector mult(AbstractComplexVector y, [AbstractComplexVector z = null, Float64List alpha = null, Float64List beta = null, bool transposeA = false]) {
    if (alpha == null) {
      alpha = new Float64List.fromList([1.0, 0.0]);
    }
    if (beta == null) {
      beta = (z == null ? new Float64List.fromList([1.0, 0.0]) : new Float64List.fromList([0.0, 0.0]));
    }
    final int rowsA = transposeA ? _columns : _rows;
    final int columnsA = transposeA ? _rows : _columns;

    bool ignore = (z == null || !transposeA);
    if (z == null) z = new ComplexVector(rowsA);

    if (!(y is ComplexVector && z is ComplexVector)) {
      return super.mult(y, z, alpha, beta, transposeA);
    }

    if (columnsA != y.length || rowsA > z.length) throw new ArgumentError("Incompatible args: " + ((transposeA ? dice() : this).toStringShort()) + ", " + y.toStringShort() + ", " + z.toStringShort());

    ComplexVector zz = z as ComplexVector;
    final Float64List elementsZ = zz._elements;
    final int strideZ = zz.stride();
    final int zeroZ = z.index(0);

    ComplexVector yy = y as ComplexVector;
    final Float64List elementsY = yy._elements;
    final int strideY = yy.stride();
    final int zeroY = y.index(0);
    //int nthreads = ConcurrencyUtils.getNumberOfThreads();

    if (transposeA) {
      if ((!ignore) && !((beta[0] == 1) && (beta[1] == 0))) {
        z.forEach(cfunc.multiply(beta));
      }

      /*if ((nthreads > 1) && (cardinality() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
        nthreads = 2;
        List<Future> futures = new List<Future>(nthreads);
        final Float64List result = new Float64List(2 * rowsA);
        int k = _rows / nthreads;
        for (int j = 0; j < nthreads; j++) {
          final int firstRow = j * k;
          final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
          final int threadID = j;
          futures[j] = ConcurrencyUtils.submit(() {
            Float64List yElem = new Float64List(2);
            Float64List val = new Float64List(2);
            if (threadID == 0) {
              for (int i = firstRow; i < lastRow; i++) {
                int high = _rowPointers[i + 1];
                yElem[0] = elementsY[zeroY + strideY * i];
                yElem[1] = elementsY[zeroY + strideY * i + 1];
                yElem = Complex.mult(alpha, yElem);
                for (int k = _rowPointers[i]; k < high; k++) {
                  int j = _columnIndexes[k];
                  val[0] = _values[2 * k];
                  val[1] = -_values[2 * k + 1];
                  val = Complex.mult(val, yElem);
                  elementsZ[zeroZ + strideZ * j] += val[0];
                  elementsZ[zeroZ + strideZ * j + 1] += val[1];
                }
              }
            } else {
              for (int i = firstRow; i < lastRow; i++) {
                int high = _rowPointers[i + 1];
                yElem[0] = elementsY[zeroY + strideY * i];
                yElem[1] = elementsY[zeroY + strideY * i + 1];
                yElem = Complex.mult(alpha, yElem);
                for (int k = _rowPointers[i]; k < high; k++) {
                  int j = _columnIndexes[k];
                  val[0] = _values[2 * k];
                  val[1] = -_values[2 * k + 1];
                  val = Complex.mult(val, yElem);
                  result[2 * j] += val[0];
                  result[2 * j + 1] += val[1];
                }
              }
            }
          });
        }
        ConcurrencyUtils.waitForCompletion(futures);
        for (int j = 0; j < rowsA; j++) {
          elementsZ[zeroZ + j * strideZ] += result[2 * j];
          elementsZ[zeroZ + j * strideZ + 1] += result[2 * j + 1];
        }
      } else {*/
        Float64List yElem = new Float64List(2);
        Float64List val = new Float64List(2);
        for (int i = 0; i < _rows; i++) {
          int high = _rowPointers[i + 1];
          yElem[0] = elementsY[zeroY + strideY * i];
          yElem[1] = elementsY[zeroY + strideY * i + 1];
          yElem = Complex.multiply(alpha, yElem);
          for (int k = _rowPointers[i]; k < high; k++) {
            int j = _columnIndexes[k];
            val[0] = _values[2 * k];
            val[1] = -_values[2 * k + 1];
            val = Complex.multiply(val, yElem);
            elementsZ[zeroZ + strideZ * j] += val[0];
            elementsZ[zeroZ + strideZ * j + 1] += val[1];
          }
        }
      //}

      return z;
    }


    /*if ((nthreads > 1) && (cardinality() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int zidx = zeroZ + firstRow * strideZ;
          Float64List yElem = new Float64List(2);
          Float64List val = new Float64List(2);
          if (beta[0] == 0.0 && beta[1] == 0) {
            for (int i = firstRow; i < lastRow; i++) {
              Float64List sum = new Float64List(2);
              int high = _rowPointers[i + 1];
              for (int k = _rowPointers[i]; k < high; k++) {
                yElem[0] = elementsY[zeroY + strideY * _columnIndexes[k]];
                yElem[1] = elementsY[zeroY + strideY * _columnIndexes[k] + 1];
                val[0] = _values[2 * k];
                val[1] = _values[2 * k + 1];
                sum = Complex.plus(sum, Complex.mult(val, yElem));
              }
              sum = Complex.mult(alpha, sum);
              elementsZ[zidx] = sum[0];
              elementsZ[zidx + 1] = sum[1];
              zidx += strideZ;
            }
          } else {
            Float64List zElem = new Float64List(2);
            for (int i = firstRow; i < lastRow; i++) {
              Float64List sum = new Float64List(2);
              int high = _rowPointers[i + 1];
              for (int k = _rowPointers[i]; k < high; k++) {
                yElem[0] = elementsY[zeroY + strideY * _columnIndexes[k]];
                yElem[1] = elementsY[zeroY + strideY * _columnIndexes[k] + 1];
                val[0] = _values[2 * k];
                val[1] = _values[2 * k + 1];
                sum = Complex.plus(sum, Complex.mult(val, yElem));
              }
              sum = Complex.mult(alpha, sum);
              zElem[0] = elementsZ[zidx];
              zElem[1] = elementsZ[zidx + 1];
              zElem = Complex.mult(beta, zElem);
              elementsZ[zidx] = sum[0] + zElem[0];
              elementsZ[zidx + 1] = sum[1] + zElem[1];
              zidx += strideZ;
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      int zidx = zeroZ;
      Float64List yElem = new Float64List(2);
      Float64List val = new Float64List(2);
      if (beta[0] == 0.0 && beta[1] == 0) {
        for (int i = 0; i < _rows; i++) {
          Float64List sum = new Float64List(2);
          int high = _rowPointers[i + 1];
          for (int k = _rowPointers[i]; k < high; k++) {
            yElem[0] = elementsY[zeroY + strideY * _columnIndexes[k]];
            yElem[1] = elementsY[zeroY + strideY * _columnIndexes[k] + 1];
            val[0] = _values[2 * k];
            val[1] = _values[2 * k + 1];
            sum = Complex.plus(sum, Complex.multiply(val, yElem));
          }
          sum = Complex.multiply(alpha, sum);
          elementsZ[zidx] = sum[0];
          elementsZ[zidx + 1] = sum[1];
          zidx += strideZ;
        }
      } else {
        Float64List zElem = new Float64List(2);
        for (int i = 0; i < _rows; i++) {
          Float64List sum = new Float64List(2);
          int high = _rowPointers[i + 1];
          for (int k = _rowPointers[i]; k < high; k++) {
            yElem[0] = elementsY[zeroY + strideY * _columnIndexes[k]];
            yElem[1] = elementsY[zeroY + strideY * _columnIndexes[k] + 1];
            val[0] = _values[2 * k];
            val[1] = _values[2 * k + 1];
            sum = Complex.plus(sum, Complex.multiply(val, yElem));
          }
          sum = Complex.multiply(alpha, sum);
          zElem[0] = elementsZ[zidx];
          zElem[1] = elementsZ[zidx + 1];
          zElem = Complex.multiply(beta, zElem);
          elementsZ[zidx] = sum[0] + zElem[0];
          elementsZ[zidx + 1] = sum[1] + zElem[1];
          zidx += strideZ;
        }
      }
    //}
    return z;
  }

  AbstractComplexMatrix multiply(AbstractComplexMatrix B, [AbstractComplexMatrix C = null, Float64List alpha = null, Float64List beta = null, bool transposeA = false, bool transposeB = false]) {
    if (alpha == null) {
      alpha = new Float64List.fromList([1.0, 0.0]);
    }
    if (beta == null) {
      beta = (C == null ? new Float64List.fromList([1.0, 0.0]) : new Float64List.fromList([0.0, 0.0]));
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
      if (B is SparseRCComplexMatrix) {
        C = new SparseRCComplexMatrix(rowsA, p, rowsA * p);
      } else {
        C = new ComplexMatrix(rowsA, p);
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

    if (!ignore && !(beta[0] == 1.0 && beta[1] == 0)) {
      C.forEach(cfunc.multiply(beta));
    }

    if ((B is ComplexMatrix) && (C is ComplexMatrix)) {
      SparseRCComplexMatrix AA;
      if (transposeA) {
        AA = conjugateTranspose();
      } else {
        AA = this;
      }
      ComplexMatrix BB;
      if (transposeB) {
        BB = B.conjugateTranspose as ComplexMatrix;
      } else {
        BB = B;
      }

      ComplexMatrix CC = C;
      Int32List rowPointersA = AA._rowPointers;
      Int32List columnIndexesA = AA._columnIndexes;
      Float64List valuesA = AA._values;
      Float64List valA = new Float64List(2);
      for (int ii = 0; ii < rowsA; ii++) {
        int highA = rowPointersA[ii + 1];
        for (int ka = rowPointersA[ii]; ka < highA; ka++) {
          valA[0] = valuesA[2 * ka];
          valA[1] = valuesA[2 * ka + 1];
          Float64List scal = Complex.multiply(alpha, valA);
          int jj = columnIndexesA[ka];
          CC.row(ii).forEachWith(BB.row(jj), cfunc.plusMultSecond(scal));
        }
      }
    } else if ((B is SparseRCComplexMatrix) && (C is SparseRCComplexMatrix)) {
      SparseRCComplexMatrix AA;
      SparseRCComplexMatrix BB;
      SparseRCComplexMatrix CC = C;
      if (transposeA) {
        AA = conjugateTranspose();
      } else {
        AA = this;
      }
      if (transposeB) {
        BB = B.conjugateTranspose();
      } else {
        BB = B;
      }

      Int32List rowPointersA = AA._rowPointers;
      Int32List columnIndexesA = AA._columnIndexes;
      Float64List valuesA = AA._values;

      Int32List rowPointersB = BB._rowPointers;
      Int32List columnIndexesB = BB._columnIndexes;
      Float64List valuesB = BB._values;

      Int32List rowPointersC = CC._rowPointers;
      Int32List columnIndexesC = CC._columnIndexes;
      Float64List valuesC = CC._values;
      int nzmax = columnIndexesC.length;

      Int32List iw = new Int32List(columnsB + 1);
      for (int i = 0; i < iw.length; i++) {
        iw[i] = -1;
      }
      int len = -1;
      Float64List valA = new Float64List(2);
      Float64List valB = new Float64List(2);
      Float64List valC = new Float64List(2);
      for (int ii = 0; ii < rowsA; ii++) {
        int highA = rowPointersA[ii + 1];
        for (int ka = rowPointersA[ii]; ka < highA; ka++) {
          valA[0] = valuesA[2 * ka];
          valA[1] = valuesA[2 * ka + 1];
          Float64List scal = Complex.multiply(alpha, valA);
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
              valB[0] = valuesB[2 * kb];
              valB[1] = valuesB[2 * kb + 1];
              valB = Complex.multiply(scal, valB);
              valuesC[2 * len] = valB[0];
              valuesC[2 * len + 1] = valB[1];
            } else {
              valB[0] = valuesB[2 * kb];
              valB[1] = valuesB[2 * kb + 1];
              valB = Complex.multiply(scal, valB);
              valuesC[2 * jpos] += valB[0];
              valuesC[2 * jpos + 1] += valB[1];
            }
          }
        }
        for (int k = rowPointersC[ii]; k < len + 1; k++) {
          iw[columnIndexesC[k]] = -1;
        }
        rowPointersC[ii + 1] = len + 1;

        //                int length = rowPointersC[ii + 1] - rowPointersC[ii];
        //                IntVector columnIndexesCPart = columnIndexesC.viewPart(rowPointersC[ii], length);
        //                Int32List indexes = cern.colt.matrix.tint.algo.IntSorting.quickSort.sortIndex(columnIndexesCPart);
        //                Arrays.sort(columnIndexesCElements, rowPointersC[ii], rowPointersC[ii + 1]);
        //                ComplexVector valuesCPart = valuesC.viewPart(rowPointersC[ii], length).viewSelection(indexes);
        //                valuesC.viewPart(rowPointersC[ii], length).assign(valuesCPart);
      }
      //            CC.columnIndexes.elements((Int32List) columnIndexesC.elements);
      //            CC.columnIndexes.setSize(columnIndexesSize);
      //            CC.values.elements((Float64List) valuesC.elements);
      //            CC.values.setSize(columnIndexesSize);
    } else {
      if (transposeB) {
        B = B.conjugateTranspose();
      }
      // cache views
      final List<AbstractComplexVector> Brows = new List<AbstractComplexVector>(columnsA);
      for (int i = columnsA; --i >= 0; ) {
        Brows[i] = B.row(i);
      }
      final List<AbstractComplexVector> Crows = new List<AbstractComplexVector>(rowsA);
      for (int i = rowsA; --i >= 0; ) {
        Crows[i] = C.row(i);
      }

      final cfunc.ComplexPlusMultSecond fun = cfunc.ComplexPlusMultSecond.plusMult(new Float64List(2));

      final Int32List columnIndexesA = _columnIndexes;
      final Float64List valuesA = _values;
      Float64List valA = new Float64List(2);
      for (int i = _rows; --i >= 0; ) {
        int low = _rowPointers[i];
        for (int k = _rowPointers[i + 1]; --k >= low; ) {
          int j = columnIndexesA[k];
          valA[0] = valuesA[2 * k];
          valA[1] = valuesA[2 * k + 1];
          fun.multiplicator = Complex.multiply(valA, alpha);
          if (!transposeA) {
            Crows[i].forEachWith(Brows[j], fun);
          } else {
            Crows[j].forEachWith(Brows[i], fun);
          }
        }
      }
    }
    return C;
  }

  static double _cumsum(Int32List p, Int32List c, int n) {
    int nz = 0;
    double nz2 = 0.0;
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
    Float64List valuesNew = new Float64List(2 * nzmax);
    length = Math.min(nzmax, _values.length);
    //System.arraycopy(_values, 0, valuesNew, 0, length);
    valuesNew.setAll(0, _values);
    _values = valuesNew;
  }

  AbstractComplexMatrix _getContent() {
    return this;
  }

  void _insert(int row, int column, int index, Float64List value) {
    var columnIndexesList = new List.from(_columnIndexes);
    //columnIndexesList._setSizeRaw(_rowPointers[_rows]);
    var valuesList = new List.from(_values);
    //valuesList._setSizeRaw(2 * _rowPointers[_rows]);
    columnIndexesList.insert(index, column);
    valuesList.insert(2 * index, value[0]);
    valuesList.insert(2 * index + 1, value[1]);
    for (int i = _rowPointers.length; --i > row; ) {
      _rowPointers[i]++;
    }
    _columnIndexes = new Int32List.fromList(columnIndexesList);
    _values = new Float64List.fromList(valuesList);
  }

  void _insertParts(int row, int column, int index, double re, double im) {
    var columnIndexesList = new List.from(_columnIndexes);
    //columnIndexesList._setSizeRaw(_rowPointers[_rows]);
    var valuesList = new List.from(_values);
    //valuesList._setSizeRaw(2 * _rowPointers[_rows]);
    columnIndexesList.insert(index, column);
    valuesList.insert(2 * index, re);
    valuesList.insert(2 * index + 1, im);
    for (int i = _rowPointers.length; --i > row; ) {
      _rowPointers[i]++;
    }
    _columnIndexes = new Int32List.fromList(columnIndexesList);
    _values = new Float64List.fromList(valuesList);
  }

  void _remove(int row, int index) {
    var columnIndexesList = new List.from(_columnIndexes);
    //columnIndexesList._setSizeRaw(_rowPointers[_rows]);
    var valuesList = new List.from(_values);
    //valuesList._setSizeRaw(_rowPointers[_rows]);
    columnIndexesList.remove(index);
    valuesList.remove(2 * index);
    valuesList.remove(2 * index + 1);
    for (int i = _rowPointers.length; --i > row; ) {
      _rowPointers[i]--;
    }
    _columnIndexes = new Int32List.fromList(columnIndexesList);
    _values = new Float64List.fromList(valuesList);
  }

  Object clone() {
    return new SparseRCComplexMatrix._internal(_rows, _columns, _rowPointers, _columnIndexes, _values);
  }

  double normInf() {
    Float64List elem = new Float64List(2);
    double norm = 0.0;
    for (int j = 0; j < _rows; j++) {
      var s = 0.0;
      for (var p = _rowPointers[j]; p < _rowPointers[j + 1]; p++) {
        elem[0] = _values[2 * p];
        elem[1] = _values[2 * p + 1];
        s += cfunc.abs(elem);
      }
      norm = Math.max(norm, s);
    }
    return norm;
  }
}
