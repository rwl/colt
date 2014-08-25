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
class SparseRCDComplexMatrix2D extends WrapperDComplexMatrix2D {

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
  factory SparseRCDComplexMatrix2D.fromValues(List<Float64List> values) {
    return new SparseRCDComplexMatrix2D.sized(values.length, values.length == 0 ? 0 : values[0].length)..assignList(values);
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
  /*SparseRCDComplexMatrix2D(int rows, int columns) {
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
  factory SparseRCDComplexMatrix2D.sized(int rows, int columns, [int nzmax = null]) {
    if (nzmax == null) {
      nzmax = 10 * rows;//Math.min(10 * rows, Integer.MAX_VALUE);
    }
    try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if ("matrix too large" != exc.message) {
        throw exc;
      }
    }
    _columnIndexes = new Int32List(nzmax);
    _values = new Float64List(2 * nzmax);
    _rowPointers = new Int32List(rows + 1);
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
  factory SparseRCDComplexMatrix2D.value(int rows, int columns, Int32List rowIndexes, Int32List columnIndexes, double re, double im, bool removeDuplicates) : super(null) {
    try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if ("matrix too large" != exc.message) {
        throw exc;
      }
    }
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
    final m = new SparseRCDComplexMatrix2D(rows, columns, _rowPointers, _columnIndexes, _values);
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
  factory SparseRCDComplexMatrix2D.values(int rows, int columns, Int32List rowIndexes, Int32List columnIndexes, Float64List values, bool removeDuplicates, bool removeZeroes) : super(null) {
    try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if ("matrix too large" != exc.message) {
        throw exc;
      }
    }
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
    final m = new SparseRCDComplexMatrix2D(rows, columns, _rowPointers, _columnIndexes, _values);
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
  SparseRCDComplexMatrix2D(int rows, int columns, Int32List rowPointers, Int32List columnIndexes, Float64List values) : super(null) {
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

  DComplexMatrix2D assign(final cfunc.DComplexDComplexFunction function) {
    if (function is cfunc.DComplexMult) { // x[i] = mult*x[i]
      final Float64List alpha = (function as cfunc.DComplexMult).multiplicator;
      if (alpha[0] == 1 && alpha[1] == 0) {
        return this;
      }
      if (alpha[0] == 0 && alpha[1] == 0) {
        return assignValues(alpha);
      }
      if (alpha[0] != alpha[0] || alpha[1] != alpha[1]) {
        return assignValues(alpha); // the funny definition of isNaN(). This should better not happen.
      }
      int nz = cardinality();
      Float64List elem = new Float64List(2);
      for (int j = 0; j < nz; j++) {
        elem[0] = _values[2 * j];
        elem[1] = _values[2 * j + 1];
        elem = DComplex.mult(elem, alpha);
        _values[2 * j] = elem[0];
        _values[2 * j + 1] = elem[1];
      }
    } else {
      forEachNonZero((int i, int j, Float64List value) {
        return function(value);
      });
    }
    return this;
  }

  DComplexMatrix2D assignValue(double re, double im) {
    if (re == 0 && im == 0) {
      Arrays.fill(_rowPointers, 0);
      Arrays.fill(_columnIndexes, 0);
      Arrays.fill(_values, 0);
    } else {
      int nnz = cardinality();
      for (int i = 0; i < nnz; i++) {
        _values[2 * i] = re;
        _values[2 * i + 1] = im;
      }
    }
    return this;
  }

  DComplexMatrix2D assignMatrix(DComplexMatrix2D source) {
    if (source == this) return this; // nothing to do
    checkShape(source);

    if (source is SparseRCDComplexMatrix2D) {
      SparseRCDComplexMatrix2D other = source as SparseRCDComplexMatrix2D;
      System.arraycopy(other._rowPointers, 0, _rowPointers, 0, _rows + 1);
      int nzmax = other._columnIndexes.length;
      if (_columnIndexes.length < nzmax) {
        _columnIndexes = new Int32List(nzmax);
        _values = new Float64List(2 * nzmax);
      }
      System.arraycopy(other._columnIndexes, 0, _columnIndexes, 0, nzmax);
      System.arraycopy(other._values, 0, _values, 0, other._values.length);
    } else if (source is SparseCCDComplexMatrix2D) {
      SparseCCDComplexMatrix2D other = (source as SparseCCDComplexMatrix2D).getConjugateTranspose();
      _rowPointers = other.getColumnPointers();
      _columnIndexes = other.getRowIndexes();
      _values = other.getValues();
    } else {
      assignValue(0, 0);
      source.forEachNonZero((int i, int j, Float64List value) {
        setQuick(i, j, value);
        return value;
      });
    }
    return this;
  }

  DComplexMatrix2D assignFunc(final DComplexMatrix2D y, cfunc.DComplexDComplexDComplexFunction function) {
    checkShape(y);
    if ((y is SparseRCDComplexMatrix2D) && (function == cfunc.DComplexFunctions.plus)) { // x[i] = x[i] + y[i]
      SparseRCDComplexMatrix2D yy = y as SparseRCDComplexMatrix2D;

      final Int32List rowPointersY = yy._rowPointers;
      final Int32List columnIndexesY = yy._columnIndexes;
      final Float64List valuesY = yy._values;

      final Int32List rowPointersC = new Int32List(_rows + 1);
      int cnz = Math.max(_columnIndexes.length, Math.min(Integer.MAX_VALUE, _rowPointers[_rows] + rowPointersY[_rows]));
      final Int32List columnIndexesC = new Int32List(cnz);
      final Float64List valuesC = new Float64List(2 * cnz);
      int nrow = _rows;
      int ncol = _columns;
      int nzmax = cnz;
      if (function == cfunc.DComplexFunctions.plus) { // x[i] = x[i] + y[i]
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
        return this;
      }
    }

    if (function is cfunc.DComplexPlusMultSecond) { // x[i] = x[i] + alpha*y[i]
      final Float64List alpha = (function as cfunc.DComplexPlusMultSecond).multiplicator;
      if (alpha[0] == 0 && alpha[1] == 0) {
        return this; // nothing to do
      }
      y.forEachNonZero((int i, int j, Float64List value) {
        setQuick(i, j, DComplex.plus(getQuick(i, j), DComplex.mult(alpha, value)));
        return value;
      });
      return this;
    }

    if (function is cfunc.DComplexPlusMultFirst) { // x[i] = alpha*x[i] + y[i]
      final Float64List alpha = (function as cfunc.DComplexPlusMultFirst).multiplicator;
      if (alpha[0] == 0 && alpha[1] == 0) {
        return assignMatrix(y);
      }
      y.forEachNonZero((int i, int j, Float64List value) {
        setQuick(i, j, DComplex.plus(DComplex.mult(alpha, getQuick(i, j)), value));
        return value;
      });
      return this;
    }

    if (function == cfunc.DComplexFunctions.mult) { // x[i] = x[i] * y[i]
      Float64List elem = new Float64List(2);
      for (int i = 0; i < _rows; i++) {
        int high = _rowPointers[i + 1];
        for (int k = _rowPointers[i]; k < high; k++) {
          int j = _columnIndexes[k];
          elem[0] = _values[2 * k];
          elem[1] = _values[2 * k + 1];
          elem = DComplex.mult(elem, y.getQuick(i, j));
          _values[2 * k] = elem[0];
          _values[2 * k + 1] = elem[1];
          if (_values[2 * k] == 0 && _values[2 * k + 1] == 0) _remove(i, j);
        }
      }
      return this;
    }

    if (function == cfunc.DComplexFunctions.div) { // x[i] = x[i] / y[i]
      Float64List elem = new Float64List(2);
      for (int i = 0; i < _rows; i++) {
        int high = _rowPointers[i + 1];
        for (int k = _rowPointers[i]; k < high; k++) {
          int j = _columnIndexes[k];
          elem[0] = _values[2 * k];
          elem[1] = _values[2 * k + 1];
          elem = DComplex.div(elem, y.getQuick(i, j));
          _values[2 * k] = elem[0];
          _values[2 * k + 1] = elem[1];
          if (_values[2 * k] == 0 && _values[2 * k + 1] == 0) _remove(i, j);
        }
      }
      return this;
    }
    return super.assignFunc(y, function);

  }

  int cardinality() {
    return _rowPointers[_rows];
  }

  DComplexMatrix2D forEachNonZero(final cfunc.IntIntDComplexFunction function) {
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
    return this;
  }

  /**
   * Returns a new matrix that has the same elements as this matrix, but is in
   * a column-compressed form. This method creates a new object (not a view),
   * so changes in the returned matrix are NOT reflected in this matrix.
   *
   * @return this matrix in a column-compressed form
   */
  SparseCCDComplexMatrix2D getColumnCompressed() {
    SparseRCDComplexMatrix2D tr = getConjugateTranspose();
    SparseCCDComplexMatrix2D cc = new SparseCCDComplexMatrix2D(_rows, _columns);
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
  Int32List getColumnIndexes() {
    return _columnIndexes;
  }

  /**
   * Returns a new matrix that has the same elements as this matrix, but is in
   * a dense form. This method creates a new object (not a view), so changes
   * in the returned matrix are NOT reflected in this matrix.
   *
   * @return this matrix in a dense form
   */
  DenseDComplexMatrix2D getDense() {
    final DenseDComplexMatrix2D dense = new DenseDComplexMatrix2D(_rows, _columns);
    forEachNonZero((int i, int j, Float64List value) {
      dense.setQuick(i, j, getQuick(i, j));
      return value;
    });
    return dense;
  }

  Float64List getQuick(int row, int column) {
    //        int k = cern.colt.Sorting.binarySearchFromTo(columnIndexes, column, rowPointers[row], rowPointers[row + 1] - 1);
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
  Int32List getRowPointers() {
    return _rowPointers;
  }

  /**
   * Returns a new matrix that is the transpose of this matrix. This method
   * creates a new object (not a view), so changes in the returned matrix are
   * NOT reflected in this matrix.
   *
   * @return the transpose of this matrix
   */
  SparseRCDComplexMatrix2D getTranspose() {
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
    SparseRCDComplexMatrix2D T = new SparseRCDComplexMatrix2D(_columns, _rows);
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
  SparseRCDComplexMatrix2D getConjugateTranspose() {
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
    SparseRCDComplexMatrix2D T = new SparseRCDComplexMatrix2D(_columns, _rows);
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
  Float64List getValues() {
    return _values;
  }

  DComplexMatrix2D like2D(int rows, int columns) {
    return new SparseRCDComplexMatrix2D(rows, columns);
  }

  DComplexMatrix1D like1D(int size) {
    return new SparseDComplexMatrix1D(size);
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
        if (DComplex.abs(elem) > eps) {
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

  void setQuick(int row, int column, Float64List value) {
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

  void setPartsQuick(int row, int column, double re, double im) {
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
    StringBuilder builder = new StringBuilder();
    builder.append(_rows).append(" x ").append(_columns).append(" sparse matrix, nnz = ").append(cardinality()).append('\n');
    for (int i = 0; i < _rows; i++) {
      int high = _rowPointers[i + 1];
      for (int j = _rowPointers[i]; j < high; j++) {
        if (_values[2 * j + 1] > 0) {
          builder.append('(').append(i).append(',').append(_columnIndexes[j]).append(')').append('\t').append(_values[2 * j]).append('+').append(_values[2 * j + 1]).append('i').append('\n');
        } else if (_values[2 * j + 1] == 0) {
          builder.append('(').append(i).append(',').append(_columnIndexes[j]).append(')').append('\t').append(_values[2 * j]).append('\n');
        } else {
          builder.append('(').append(i).append(',').append(_columnIndexes[j]).append(')').append('\t').append(_values[2 * j]).append('-').append(_values[2 * j + 1]).append('i').append('\n');
        }
      }
    }
    return builder.toString();
  }

  void trimToSize() {
    _realloc(0);
  }

  DComplexMatrix1D zMult(DComplexMatrix1D y, DComplexMatrix1D z, final Float64List alpha, final Float64List beta, final bool transposeA) {
    final int rowsA = transposeA ? _columns : _rows;
    final int columnsA = transposeA ? _rows : _columns;

    bool ignore = (z == null || !transposeA);
    if (z == null) z = new DenseDComplexMatrix1D(rowsA);

    if (!(y is DenseDComplexMatrix1D && z is DenseDComplexMatrix1D)) {
      return super.zMult(y, z, alpha, beta, transposeA);
    }

    if (columnsA != y.size() || rowsA > z.size()) throw new ArgumentError("Incompatible args: " + ((transposeA ? viewDice() : this).toStringShort()) + ", " + y.toStringShort() + ", " + z.toStringShort());

    DenseDComplexMatrix1D zz = z as DenseDComplexMatrix1D;
    final Float64List elementsZ = zz._elements;
    final int strideZ = zz.stride();
    final int zeroZ = z.index(0);

    DenseDComplexMatrix1D yy = y as DenseDComplexMatrix1D;
    final Float64List elementsY = yy._elements;
    final int strideY = yy.stride();
    final int zeroY = y.index(0);
    int nthreads = ConcurrencyUtils.getNumberOfThreads();

    if (transposeA) {
      if ((!ignore) && !((beta[0] == 1) && (beta[1] == 0))) z.assign(cfunc.DComplexFunctions.mult(beta));

      if ((nthreads > 1) && (cardinality() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
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
                yElem = DComplex.mult(alpha, yElem);
                for (int k = _rowPointers[i]; k < high; k++) {
                  int j = _columnIndexes[k];
                  val[0] = _values[2 * k];
                  val[1] = -_values[2 * k + 1];
                  val = DComplex.mult(val, yElem);
                  elementsZ[zeroZ + strideZ * j] += val[0];
                  elementsZ[zeroZ + strideZ * j + 1] += val[1];
                }
              }
            } else {
              for (int i = firstRow; i < lastRow; i++) {
                int high = _rowPointers[i + 1];
                yElem[0] = elementsY[zeroY + strideY * i];
                yElem[1] = elementsY[zeroY + strideY * i + 1];
                yElem = DComplex.mult(alpha, yElem);
                for (int k = _rowPointers[i]; k < high; k++) {
                  int j = _columnIndexes[k];
                  val[0] = _values[2 * k];
                  val[1] = -_values[2 * k + 1];
                  val = DComplex.mult(val, yElem);
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
      } else {
        Float64List yElem = new Float64List(2);
        Float64List val = new Float64List(2);
        for (int i = 0; i < _rows; i++) {
          int high = _rowPointers[i + 1];
          yElem[0] = elementsY[zeroY + strideY * i];
          yElem[1] = elementsY[zeroY + strideY * i + 1];
          yElem = DComplex.mult(alpha, yElem);
          for (int k = _rowPointers[i]; k < high; k++) {
            int j = _columnIndexes[k];
            val[0] = _values[2 * k];
            val[1] = -_values[2 * k + 1];
            val = DComplex.mult(val, yElem);
            elementsZ[zeroZ + strideZ * j] += val[0];
            elementsZ[zeroZ + strideZ * j + 1] += val[1];
          }
        }
      }

      return z;
    }


    if ((nthreads > 1) && (cardinality() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
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
                sum = DComplex.plus(sum, DComplex.mult(val, yElem));
              }
              sum = DComplex.mult(alpha, sum);
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
                sum = DComplex.plus(sum, DComplex.mult(val, yElem));
              }
              sum = DComplex.mult(alpha, sum);
              zElem[0] = elementsZ[zidx];
              zElem[1] = elementsZ[zidx + 1];
              zElem = DComplex.mult(beta, zElem);
              elementsZ[zidx] = sum[0] + zElem[0];
              elementsZ[zidx + 1] = sum[1] + zElem[1];
              zidx += strideZ;
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {
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
            sum = DComplex.plus(sum, DComplex.mult(val, yElem));
          }
          sum = DComplex.mult(alpha, sum);
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
            sum = DComplex.plus(sum, DComplex.mult(val, yElem));
          }
          sum = DComplex.mult(alpha, sum);
          zElem[0] = elementsZ[zidx];
          zElem[1] = elementsZ[zidx + 1];
          zElem = DComplex.mult(beta, zElem);
          elementsZ[zidx] = sum[0] + zElem[0];
          elementsZ[zidx + 1] = sum[1] + zElem[1];
          zidx += strideZ;
        }
      }
    }
    return z;
  }

  DComplexMatrix2D zMult2D(DComplexMatrix2D B, DComplexMatrix2D C, final Float64List alpha, Float64List beta, final bool transposeA, bool transposeB) {
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
      if (B is SparseRCDComplexMatrix2D) {
        C = new SparseRCDComplexMatrix2D(rowsA, p, (rowsA * p));
      } else {
        C = new DenseDComplexMatrix2D(rowsA, p);
      }
    }

    if (rowsB != columnsA) throw new ArgumentError("Matrix2D inner dimensions must agree:" + toStringShort() + ", " + (transposeB ? B.viewDice() : B).toStringShort());
    if (C.rows() != rowsA || C.columns() != p) throw new ArgumentError("Incompatible result matrix: " + toStringShort() + ", " + (transposeB ? B.viewDice() : B).toStringShort() + ", " + C.toStringShort());
    if (this == C || B == C) throw new ArgumentError("Matrices must not be identical");

    if (!ignore && !(beta[0] == 1.0 && beta[1] == 0)) {
      C.assign(cfunc.DComplexFunctions.mult(beta));
    }

    if ((B is DenseDComplexMatrix2D) && (C is DenseDComplexMatrix2D)) {
      SparseRCDComplexMatrix2D AA;
      if (transposeA) {
        AA = getConjugateTranspose();
      } else {
        AA = this;
      }
      DenseDComplexMatrix2D BB;
      if (transposeB) {
        BB = B.getConjugateTranspose() as DenseDComplexMatrix2D;
      } else {
        BB = B as DenseDComplexMatrix2D;
      }

      DenseDComplexMatrix2D CC = C as DenseDComplexMatrix2D;
      Int32List rowPointersA = AA._rowPointers;
      Int32List columnIndexesA = AA._columnIndexes;
      Float64List valuesA = AA._values;
      Float64List valA = new Float64List(2);
      for (int ii = 0; ii < rowsA; ii++) {
        int highA = rowPointersA[ii + 1];
        for (int ka = rowPointersA[ii]; ka < highA; ka++) {
          valA[0] = valuesA[2 * ka];
          valA[1] = valuesA[2 * ka + 1];
          Float64List scal = DComplex.mult(alpha, valA);
          int jj = columnIndexesA[ka];
          CC.viewRow(ii).assignFunc(BB.viewRow(jj), DComplexFunctions.plusMultSecond(scal));
        }
      }
    } else if ((B is SparseRCDComplexMatrix2D) && (C is SparseRCDComplexMatrix2D)) {
      SparseRCDComplexMatrix2D AA;
      SparseRCDComplexMatrix2D BB;
      SparseRCDComplexMatrix2D CC = C as SparseRCDComplexMatrix2D;
      if (transposeA) {
        AA = getConjugateTranspose();
      } else {
        AA = this;
      }
      if (transposeB) {
        BB = (B as SparseRCDComplexMatrix2D).getConjugateTranspose();
      } else {
        BB = B as SparseRCDComplexMatrix2D;
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
          Float64List scal = DComplex.mult(alpha, valA);
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
              valB = DComplex.mult(scal, valB);
              valuesC[2 * len] = valB[0];
              valuesC[2 * len + 1] = valB[1];
            } else {
              valB[0] = valuesB[2 * kb];
              valB[1] = valuesB[2 * kb + 1];
              valB = DComplex.mult(scal, valB);
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
        //                IntMatrix1D columnIndexesCPart = columnIndexesC.viewPart(rowPointersC[ii], length);
        //                Int32List indexes = cern.colt.matrix.tint.algo.IntSorting.quickSort.sortIndex(columnIndexesCPart);
        //                Arrays.sort(columnIndexesCElements, rowPointersC[ii], rowPointersC[ii + 1]);
        //                DComplexMatrix1D valuesCPart = valuesC.viewPart(rowPointersC[ii], length).viewSelection(indexes);
        //                valuesC.viewPart(rowPointersC[ii], length).assign(valuesCPart);
      }
      //            CC.columnIndexes.elements((Int32List) columnIndexesC.elements());
      //            CC.columnIndexes.setSize(columnIndexesSize);
      //            CC.values.elements((Float64List) valuesC.elements());
      //            CC.values.setSize(columnIndexesSize);
    } else {
      if (transposeB) {
        B = B.getConjugateTranspose();
      }
      // cache views
      final List<DComplexMatrix1D> Brows = new List<DComplexMatrix1D>(columnsA);
      for (int i = columnsA; --i >= 0; ) {
        Brows[i] = B.viewRow(i);
      }
      final List<DComplexMatrix1D> Crows = new List<DComplexMatrix1D>(rowsA);
      for (int i = rowsA; --i >= 0; ) {
        Crows[i] = C.viewRow(i);
      }

      final cfunc.DComplexPlusMultSecond fun = cfunc.DComplexPlusMultSecond.plusMult(new Float64List(2));

      final Int32List columnIndexesA = _columnIndexes;
      final Float64List valuesA = _values;
      Float64List valA = new Float64List(2);
      for (int i = _rows; --i >= 0; ) {
        int low = _rowPointers[i];
        for (int k = _rowPointers[i + 1]; --k >= low; ) {
          int j = columnIndexesA[k];
          valA[0] = valuesA[2 * k];
          valA[1] = valuesA[2 * k + 1];
          fun.multiplicator = DComplex.mult(valA, alpha);
          if (!transposeA) {
            Crows[i].assignFunc(Brows[j], fun);
          } else {
            Crows[j].assignFunc(Brows[i], fun);
          }
        }
      }
    }
    return C;
  }

  double _cumsum(Int32List p, Int32List c, int n) {
    int nz = 0;
    double nz2 = 0;
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
    System.arraycopy(_columnIndexes, 0, columnIndexesNew, 0, length);
    _columnIndexes = columnIndexesNew;
    Float64List valuesNew = new Float64List(2 * nzmax);
    length = Math.min(nzmax, _values.length);
    System.arraycopy(_values, 0, valuesNew, 0, length);
    _values = valuesNew;
  }

  DComplexMatrix2D _getContent() {
    return this;
  }

  void _insert(int row, int column, int index, Float64List value) {
    IntArrayList columnIndexesList = new IntArrayList(_columnIndexes);
    columnIndexesList._setSizeRaw(_rowPointers[_rows]);
    DoubleArrayList valuesList = new DoubleArrayList(_values);
    valuesList._setSizeRaw(2 * _rowPointers[_rows]);
    columnIndexesList.beforeInsert(index, column);
    valuesList.beforeInsert(2 * index, value[0]);
    valuesList.beforeInsert(2 * index + 1, value[1]);
    for (int i = _rowPointers.length; --i > row; ) _rowPointers[i]++;
    _columnIndexes = columnIndexesList.elements();
    _values = valuesList.elements();
  }

  void _insertParts(int row, int column, int index, double re, double im) {
    IntArrayList columnIndexesList = new IntArrayList(_columnIndexes);
    columnIndexesList._setSizeRaw(_rowPointers[_rows]);
    DoubleArrayList valuesList = new DoubleArrayList(_values);
    valuesList._setSizeRaw(2 * _rowPointers[_rows]);
    columnIndexesList.beforeInsert(index, column);
    valuesList.beforeInsert(2 * index, re);
    valuesList.beforeInsert(2 * index + 1, im);
    for (int i = _rowPointers.length; --i > row; ) _rowPointers[i]++;
    _columnIndexes = columnIndexesList.elements();
    _values = valuesList.elements();
  }

  void _remove(int row, int index) {
    IntArrayList columnIndexesList = new IntArrayList(_columnIndexes);
    columnIndexesList._setSizeRaw(_rowPointers[_rows]);
    DoubleArrayList valuesList = new DoubleArrayList(_values);
    valuesList._setSizeRaw(_rowPointers[_rows]);
    columnIndexesList.remove(index);
    valuesList.remove(2 * index);
    valuesList.remove(2 * index + 1);
    for (int i = _rowPointers.length; --i > row; ) _rowPointers[i]--;
    _columnIndexes = columnIndexesList.elements();
    _values = valuesList.elements();
  }

}
