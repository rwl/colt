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
 * Sparse column-compressed 2-d matrix holding <tt>complex</tt> elements. First
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
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
class SparseCCComplexMatrix extends WrapperComplexMatrix {

  /*
   * Internal storage.
   */
  Int32List _columnPointers;

  Int32List _rowIndexes;

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
  factory SparseCCComplexMatrix.fromList(List<Float64List> values) {
    return new SparseCCComplexMatrix.sized(values.length, values[0].length)
      ..setAll2D(values);
  }

  /**
   * Constructs a matrix with a given internal storage.
   *
   * @param dzcs
   *            internal storage.
   */
  /*SparseCCComplexMatrix(DZcs dzcs) {
    this(dzcs.m, dzcs.n);
    _rowIndexes = dzcs.i;
    _columnPointers = dzcs.p;
    _values = dzcs.x;
  }*/

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
  /*factory SparseCCComplexMatrix.sized(int rows, int columns) {
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
  factory SparseCCComplexMatrix.sized(int rows, int columns, [int nzmax=null]) {
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
    final rowIndexes = new Int32List(nzmax);
    final values = new Float64List(2 * nzmax);
    final columnPointers = new Int32List(columns + 1);
    return new SparseCCComplexMatrix(rows, columns, rowIndexes, columnPointers, values);
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
   *            column pointers
   * @param values
   *            numerical values
   */
  SparseCCComplexMatrix(int rows, int columns, Int32List rowIndexes, Int32List columnPointers, Float64List values) : super(null) {
    try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if ("matrix too large" != exc.message) {
        throw exc;
      }
    }
    if (columnPointers.length != columns + 1) {
      throw new ArgumentError("columnPointers.length != columns + 1");
    }
    if (2 * rowIndexes.length != values.length) {
      throw new ArgumentError("2 * rowIndexes.length != values.length");
    }
    this._columnPointers = columnPointers;
    this._rowIndexes = rowIndexes;
    this._values = values;
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
   * @param re
   *            the real part of numerical value
   * @param im
   *            the imaginary part of numerical value
   * @param removeDuplicates
   *            if true, then duplicates (if any) are removed
   */
  factory SparseCCComplexMatrix.withValue(int rows, int columns, Int32List rowIndexes, Int32List columnIndexes, double re, double im, bool removeDuplicates) {
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
    final _rowIndexes = new Int32List(nz);
    final _values = new Float64List(2 * nz);
    final _columnPointers = new Int32List(columns + 1);
    Int32List w = new Int32List(columns);
    int r;
    for (int k = 0; k < nz; k++) {
      w[columnIndexes[k]]++;
    }
    _cumsum(_columnPointers, w, columns);
    for (int k = 0; k < nz; k++) {
      _rowIndexes[r = w[columnIndexes[k]]++] = rowIndexes[k];
      _values[2 * r] = re;
      _values[2 * r + 1] = im;
    }
    final m = new SparseCCComplexMatrix(rows, columns, _rowIndexes, _columnPointers, _values);
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
  factory SparseCCComplexMatrix.withValues(int rows, int columns, Int32List rowIndexes, Int32List columnIndexes, Float64List values, bool removeDuplicates, bool removeZeroes) {
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
    final _rowIndexes = new Int32List(nz);
    final _values = new Float64List(2 * nz);
    final _columnPointers = new Int32List(columns + 1);
    Int32List w = new Int32List(columns);
    int r;
    for (int k = 0; k < nz; k++) {
      w[columnIndexes[k]]++;
    }
    _cumsum(_columnPointers, w, columns);
    for (int k = 0; k < nz; k++) {
      _rowIndexes[r = w[columnIndexes[k]]++] = rowIndexes[k];
      _values[2 * r] = values[2 * k];
      _values[2 * r + 1] = values[2 * k + 1];
    }
    final m = new SparseCCComplexMatrix(rows, columns, _rowIndexes, _columnPointers, _values);
    if (removeDuplicates) {
      m.removeDuplicates();
    }
    return m;
  }

  ComplexMatrix forEach(final cfunc.ComplexComplexFunction function) {
    if (function is cfunc.ComplexMult) { // x[i] = mult*x[i]
      final Float64List alpha = (function as cfunc.ComplexMult).multiplicator;
      if (alpha[0] == 1 && alpha[1] == 0) {
        return this;
      }
      if (alpha[0] == 0 && alpha[1] == 0) {
        return fill(0.0, 0.0);
      }
      if (alpha[0] != alpha[0] || alpha[1] != alpha[1]) {
        return fill(alpha[0], alpha[1]); // the funny definition of isNaN(). This should better not happen.
      }

      final Float64List valuesE = _values;
      int nz = cardinality;
      Float64List valE = new Float64List(2);
      for (int j = 0; j < nz; j++) {
        valE[0] = valuesE[2 * j];
        valE[1] = valuesE[2 * j + 1];
        valE = Complex.multiply(valE, alpha);
        valuesE[2 * j] = valE[0];
        valuesE[2 * j + 1] = valE[1];
      }
    } else {
      forEachNonZero((int i, int j, Float64List value) {
        return function(value);
      });
    }
    return this;
  }

  ComplexMatrix fill(double re, double im) {
    if (re == 0 && im == 0) {
      //Arrays.fill(_rowIndexes, 0);
      //Arrays.fill(_columnPointers, 0);
      //Arrays.fill(_values, 0);
      _rowIndexes.fillRange(0, _rowIndexes.length, 0);
      _columnPointers.fillRange(0, _columnPointers.length, 0);
      _values.fillRange(0, _values.length, 0.0);
    } else {
      int nnz = cardinality;
      for (int i = 0; i < nnz; i++) {
        _values[2 * i] = re;
        _values[2 * i + 1] = im;
      }
    }
    return this;
  }

  ComplexMatrix copyFrom(ComplexMatrix source) {
    if (source == this) return this; // nothing to do
    checkShape(source);

    if (source is SparseCCComplexMatrix) {
      SparseCCComplexMatrix other = source;
      //System.arraycopy(other.getColumnPointers(), 0, _columnPointers, 0, _columns + 1);
      _columnPointers.setAll(0, other.columnPointers);
      int nzmax = other.rowIndexes.length;
      if (_rowIndexes.length < nzmax) {
        _rowIndexes = new Int32List(nzmax);
        _values = new Float64List(2 * nzmax);
      }
      //System.arraycopy(other.getRowIndexes(), 0, _rowIndexes, 0, nzmax);
      _rowIndexes.setAll(0, other.rowIndexes);
      //System.arraycopy(other.getValues(), 0, _values, 0, other.getValues().length);
      _values.setAll(0, other.values);
    } else if (source is SparseRCComplexMatrix) {
      SparseRCComplexMatrix other = source.conjugateTranspose();
      _columnPointers = other.rowPointers;
      _rowIndexes = other.columnIndexes;
      _values = other.values;
    } else {
      fill(0.0, 0.0);
      source.forEachNonZero((int i, int j, Float64List value) {
        set(i, j, value);
        return value;
      });
    }
    return this;
  }

  ComplexMatrix forEachMatrix(final ComplexMatrix y, cfunc.ComplexComplexComplexFunction function) {
    checkShape(y);

    if ((y is SparseCCComplexMatrix) && (function == cfunc.plus)) { // x[i] = x[i] + y[i]
      SparseCCComplexMatrix yy = y;
      int p,
          j,
          nz = 0,
          anz;
      Int32List Cp, Ci, Bp, w;
      int m, n, bnz;
      Float64List x, Cx;
      m = _rows;
      anz = _columnPointers[_columns];
      n = yy._columns;
      Bp = yy._columnPointers;
      bnz = Bp[n];
      w = new Int32List(m);
      /* get workspace */
      x = new Float64List(2 * m);
      /* get workspace */
      SparseCCComplexMatrix C = new SparseCCComplexMatrix.sized(m, n, anz + bnz);
      /* allocate result*/
      Cp = C._columnPointers;
      Ci = C._rowIndexes;
      Cx = C._values;
      Float64List one = new Float64List.fromList([1, 0]);
      for (j = 0; j < n; j++) {
        Cp[j] = nz;
        /* column j of C starts here */
        nz = _scatter(this, j, one, w, x, j + 1, C, nz);
        /* alpha*A(:,j)*/
        nz = _scatter(yy, j, one, w, x, j + 1, C, nz);
        /* beta*B(:,j) */
        for (p = Cp[j]; p < nz; p++) {
          Cx[2 * p] = x[2 * Ci[p]];
          Cx[2 * p + 1] = x[2 * Ci[p] + 1];
        }
      }
      Cp[n] = nz;
      /* finalize the last column of C */
      _rowIndexes = Ci;
      _columnPointers = Cp;
      _values = Cx;
      return this;
    }

    if (function is cfunc.ComplexPlusMultSecond) { // x[i] = x[i] + alpha*y[i]
      final Float64List alpha = (function as cfunc.ComplexPlusMultSecond).multiplicator;
      if (alpha[0] == 0 && alpha[1] == 0) {
        return this; // nothing to do
      }
      y.forEachNonZero((int i, int j, Float64List value) {
        set(i, j, Complex.plus(get(i, j), Complex.multiply(alpha, value)));
        return value;
      });
      return this;
    }

    if (function is cfunc.ComplexPlusMultFirst) { // x[i] = alpha*x[i] + y[i]
      final Float64List alpha = (function as cfunc.ComplexPlusMultFirst).multiplicator;
      if (alpha[0] == 0 && alpha[1] == 0) {
        return copyFrom(y);
      }
      y.forEachNonZero((int i, int j, Float64List value) {
        set(i, j, Complex.plus(Complex.multiply(alpha, get(i, j)), value));
        return value;
      });
      return this;
    }

    if (function == cfunc.mult) { // x[i] = x[i] * y[i]
      final Int32List rowIndexesA = _rowIndexes;
      final Int32List columnPointersA = _columnPointers;
      final Float64List valuesA = _values;
      Float64List valA = new Float64List(2);
      for (int j = _columns; --j >= 0; ) {
        int low = columnPointersA[j];
        for (int k = columnPointersA[j + 1]; --k >= low; ) {
          int i = rowIndexesA[k];
          valA[0] = valuesA[2 * k];
          valA[1] = valuesA[2 * k + 1];
          valA = Complex.multiply(valA, y.get(i, j));
          valuesA[2 * k] = valA[0];
          valuesA[2 * k + 1] = valA[1];
          if (valuesA[2 * k] == 0 && valuesA[2 * k + 1] == 0) _remove(i, j);
        }
      }
      return this;
    }

    if (function == cfunc.div) { // x[i] = x[i] / y[i]
      final Int32List rowIndexesA = _rowIndexes;
      final Int32List columnPointersA = _columnPointers;
      final Float64List valuesA = _values;

      Float64List valA = new Float64List(2);
      for (int j = _columns; --j >= 0; ) {
        int low = columnPointersA[j];
        for (int k = columnPointersA[j + 1]; --k >= low; ) {
          int i = rowIndexesA[k];
          valA[0] = valuesA[2 * k];
          valA[1] = valuesA[2 * k + 1];
          valA = Complex.div_(valA, y.get(i, j));
          valuesA[2 * k] = valA[0];
          valuesA[2 * k + 1] = valA[1];
          if (valuesA[2 * k] == 0 && valuesA[2 * k + 1] == 0) _remove(i, j);
        }
      }
      return this;
    }
    return super.forEachMatrix(y, function);
  }

  int get cardinality {
    return _columnPointers[_columns];
  }

  ComplexMatrix forEachNonZero(final cfunc.IntIntComplexFunction function) {
    final Int32List rowIndexesA = _rowIndexes;
    final Int32List columnPointersA = _columnPointers;
    final Float64List valuesA = _values;
    Float64List valA = new Float64List(2);
    for (int j = _columns; --j >= 0; ) {
      int low = columnPointersA[j];
      for (int k = columnPointersA[j + 1]; --k >= low; ) {
        int i = rowIndexesA[k];
        valA[0] = valuesA[2 * k];
        valA[1] = valuesA[2 * k + 1];
        valA = function(i, j, valA);
        valuesA[2 * k] = valA[0];
        valuesA[2 * k + 1] = valA[1];
      }
    }
    return this;
  }

  /**
   * Returns column pointers
   *
   * @return column pointers
   */
  Int32List get columnPointers {
    return _columnPointers;
  }

  /**
   * Returns a new matrix that has the same elements as this matrix, but is in
   * a dense form. This method creates a new object (not a view), so changes
   * in the returned matrix are NOT reflected in this matrix.
   *
   * @return this matrix in a dense form
   */
  DenseComplexMatrix dense() {
    final DenseComplexMatrix dense = new DenseComplexMatrix(_rows, _columns);
    forEachNonZero((int i, int j, Float64List value) {
      dense.set(i, j, get(i, j));
      return value;
    });
    return dense;
  }

  Float64List get(int row, int column) {
    //        int k = cern.colt.Sorting.binarySearchFromTo(dcs.i, row, dcs.p[column], dcs.p[column + 1] - 1);
    int k = _searchFromTo(_rowIndexes, row, _columnPointers[column], _columnPointers[column + 1] - 1);
    Float64List v = new Float64List(2);
    if (k >= 0) {
      v[0] = _values[2 * k];
      v[1] = _values[2 * k + 1];
    }
    return v;
  }

  /**
   * Returns a new matrix that has the same elements as this matrix, but is in
   * a row-compressed form. This method creates a new object (not a view), so
   * changes in the returned matrix are NOT reflected in this matrix.
   *
   * @return this matrix in a row-compressed form
   */
  SparseRCComplexMatrix rowCompressed() {
    SparseCCComplexMatrix tr = conjugateTranspose();
    SparseRCComplexMatrix rc = new SparseRCComplexMatrix.sized(_rows, _columns);
    rc._columnIndexes = tr._rowIndexes;
    rc._rowPointers = tr._columnPointers;
    rc._values = tr._values;
    return rc;
  }

  /**
   * Returns row indexes;
   *
   * @return row indexes
   */
  Int32List get rowIndexes {
    return _rowIndexes;
  }

  DZcs elements() {
    DZcs cs = new DZcs();
    cs.m = _rows;
    cs.n = _columns;
    cs.i = _rowIndexes;
    cs.p = _columnPointers;
    cs.x = _values;
    cs.nz = -1;
    cs.nzmax = _values.length ~/ 2;
    return cs;
  }

  /**
   * Returns a new matrix that is the transpose of this matrix. This method
   * creates a new object (not a view), so changes in the returned matrix are
   * NOT reflected in this matrix.
   *
   * @return the transpose of this matrix
   */
  SparseCCComplexMatrix transpose() {
    DZcs dzcst = cxsparse.cs_transpose(elements(), true);
    SparseCCComplexMatrix tr = new SparseCCComplexMatrix(dzcst.m, dzcst.n, dzcst.i, dzcst.p, dzcst.x);
    return tr;
  }

  SparseCCComplexMatrix conjugateTranspose() {
    int p, q, j, n, m;
    Int32List Cp, Ci, Ap, Ai, w;
    Float64List Cx, Ax;
    m = _rows;
    n = _columns;
    Ap = _columnPointers;
    Ai = _rowIndexes;
    Ax = _values;
    SparseCCComplexMatrix C = new SparseCCComplexMatrix.sized(_columns, _rows, Ai.length);
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
        Cx[2 * q] = Ax[2 * p];
        Cx[2 * q + 1] = -Ax[2 * p + 1];
      }
    }
    return C;
  }

  /**
   * Returns numerical values
   *
   * @return numerical values
   */
  Float64List get values {
    return _values;
  }

  ComplexMatrix like2D(int rows, int columns) {
    return new SparseCCComplexMatrix.sized(rows, columns);
  }

  ComplexVector like1D(int size) {
    return new SparseComplexVector(size);
  }

  void set(int row, int column, Float64List value) {
    //        int k = cern.colt.Sorting.binarySearchFromTo(dcs.i, row, dcs.p[column], dcs.p[column + 1] - 1);
    int k = _searchFromTo(_rowIndexes, row, _columnPointers[column], _columnPointers[column + 1] - 1);

    if (k >= 0) { // found
      if (value[0] == 0 && value[1] == 0) {
        _remove(column, k);
      } else {
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
    //        int k = cern.colt.Sorting.binarySearchFromTo(dcs.i, row, dcs.p[column], dcs.p[column + 1] - 1);
    int k = _searchFromTo(_rowIndexes, row, _columnPointers[column], _columnPointers[column + 1] - 1);

    if (k >= 0) { // found
      if (re == 0 && im == 0) {
        _remove(column, k);
      } else {
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

  /**
   * Sorts row indexes
   */
  void sortRowIndexes() {
    SparseCCComplexMatrix tr = conjugateTranspose();
    tr = tr.conjugateTranspose();
    _columnPointers = tr._columnPointers;
    _rowIndexes = tr._rowIndexes;
    _values = tr._values;
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
    Int32List Ap, Ai, w;
    Float64List Ax;
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
          Ax[2 * nz] = Ax[2 * p];
          Ax[2 * nz + 1] = Ax[2 * p + 1];
          nz++;
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
    Int32List Ap, Ai;
    Float64List Ax;
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
          Ai[nz] = Ai[p];
          Ax[2 * nz] = Ax[2 * p];
          /* keep A(i,j) */
          Ax[2 * nz + 1] = Ax[2 * p + 1];
          nz++;
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
        if (_values[2 * j + 1] > 0) {
          builder..write('(')..write(_rowIndexes[j])..write(',')..write(i)..write(')')..write('\t')..write(_values[2 * j])..write('+')..write(_values[2 * j + 1])..write('i')..write('\n');
        } else if (_values[2 * j + 1] == 0) {
          builder..write('(')..write(_rowIndexes[j])..write(',')..write(i)..write(')')..write('\t')..write(_values[2 * j])..write('\n');
        } else {
          builder..write('(')..write(_rowIndexes[j])..write(',')..write(i)..write(')')..write('\t')..write(_values[2 * j])..write('-')..write(_values[2 * j + 1])..write('i')..write('\n');
        }
      }
    }
    return builder.toString();
  }

  ComplexVector mult(ComplexVector y, ComplexVector z, [Float64List alpha = null, Float64List beta = null, bool transposeA = false]) {
    if (alpha == null) {
      alpha = new Float64List.fromList([1.0, 0.0]);
    }
    if (beta == null) {
      beta = (z == null ? new Float64List.fromList([1.0, 0.0]) : new Float64List.fromList([0.0, 0.0]));
    }
    final int rowsA = transposeA ? _columns : _rows;
    final int columnsA = transposeA ? _rows : _columns;

    bool ignore = (z == null || transposeA);
    if (z == null) z = new DenseComplexVector(rowsA);

    if (!(y is DenseComplexVector && z is DenseComplexVector)) {
      return super.mult(y, z, alpha, beta, transposeA);
    }

    if (columnsA != y.length || rowsA > z.length) throw new ArgumentError("Incompatible args: " + ((transposeA ? dice() : this).toStringShort()) + ", " + y.toStringShort() + ", " + z.toStringShort());

    DenseComplexVector zz = z as DenseComplexVector;
    final Float64List elementsZ = zz._elements;
    final int strideZ = zz.stride();
    final int zeroZ = zz.index(0);

    DenseComplexVector yy = y as DenseComplexVector;
    final Float64List elementsY = yy._elements;
    final int strideY = yy.stride();
    final int zeroY = yy.index(0);

    final Int32List rowIndexesA = _rowIndexes;
    final Int32List columnPointersA = _columnPointers;
    final Float64List valuesA = _values;

    int zidx = zeroZ;
    //int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if (!transposeA) {
      if ((!ignore) && !(beta[0] == 1 && beta[1] == 0)) {
        z.forEach(cfunc.multiply(beta));
      }

      /*if ((nthreads > 1) && (cardinality() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
        nthreads = 2;
        List<Future> futures = new List<Future>(nthreads);
        final Float64List result = new Float64List(2 * rowsA);
        int k = _columns / nthreads;
        for (int j = 0; j < nthreads; j++) {
          final int firstColumn = j * k;
          final int lastColumn = (j == nthreads - 1) ? _columns : firstColumn + k;
          final int threadID = j;
          futures[j] = ConcurrencyUtils.submit(() {
            Float64List yElem = new Float64List(2);
            Float64List valA = new Float64List(2);
            if (threadID == 0) {
              for (int i = firstColumn; i < lastColumn; i++) {
                int high = columnPointersA[i + 1];
                yElem[0] = elementsY[zeroY + strideY * i];
                yElem[1] = elementsY[zeroY + strideY * i + 1];
                for (int k = columnPointersA[i]; k < high; k++) {
                  int j = rowIndexesA[k];
                  valA[0] = valuesA[2 * k];
                  valA[1] = valuesA[2 * k + 1];
                  valA = Complex.mult(valA, yElem);
                  valA = Complex.mult(valA, alpha);
                  elementsZ[zeroZ + strideZ * j] += valA[0];
                  elementsZ[zeroZ + strideZ * j + 1] += valA[1];
                }
              }
            } else {
              for (int i = firstColumn; i < lastColumn; i++) {
                int high = columnPointersA[i + 1];
                yElem[0] = elementsY[zeroY + strideY * i];
                yElem[1] = elementsY[zeroY + strideY * i + 1];
                for (int k = columnPointersA[i]; k < high; k++) {
                  int j = rowIndexesA[k];
                  valA[0] = valuesA[2 * k];
                  valA[1] = valuesA[2 * k + 1];
                  valA = Complex.mult(valA, yElem);
                  valA = Complex.mult(valA, alpha);
                  result[2 * j] += valA[0];
                  result[2 * j + 1] += valA[1];
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
        Float64List valA = new Float64List(2);
        for (int i = 0; i < _columns; i++) {
          int high = columnPointersA[i + 1];
          yElem[0] = elementsY[zeroY + strideY * i];
          yElem[1] = elementsY[zeroY + strideY * i + 1];
          for (int k = columnPointersA[i]; k < high; k++) {
            int j = rowIndexesA[k];
            valA[0] = valuesA[2 * k];
            valA[1] = valuesA[2 * k + 1];
            valA = Complex.multiply(valA, yElem);
            valA = Complex.multiply(valA, alpha);
            elementsZ[zeroZ + strideZ * j] += valA[0];
            elementsZ[zeroZ + strideZ * j + 1] += valA[1];
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
            Float64List valA = new Float64List(2);
            Float64List valY = new Float64List(2);
            Float64List valZ = new Float64List(2);
            for (int i = firstColumn; i < lastColumn; i++) {
              Float64List sum = new Float64List(2);
              int high = _columnPointers[i + 1];
              for (int k = _columnPointers[i]; k < high; k++) {
                valA[0] = valuesA[2 * k];
                valA[1] = -valuesA[2 * k + 1];
                valY[0] = elementsY[zeroY + strideY * _rowIndexes[k]];
                valY[1] = elementsY[zeroY + strideY * _rowIndexes[k] + 1];
                sum = Complex.plus(sum, Complex.mult(valA, valY));
              }
              sum = Complex.mult(alpha, sum);
              valZ[0] = elementsZ[zidx];
              valZ[1] = elementsZ[zidx + 1];
              valZ = Complex.mult(valZ, beta);
              elementsZ[zidx] = sum[0] + valZ[0];
              elementsZ[zidx + 1] = sum[1] + valZ[1];
              zidx += strideZ;
            }
          });
        }
        ConcurrencyUtils.waitForCompletion(futures);
      } else {*/
        Float64List valA = new Float64List(2);
        Float64List valY = new Float64List(2);
        Float64List valZ = new Float64List(2);
        for (int i = 0; i < _columns; i++) {
          Float64List sum = new Float64List(2);
          int high = _columnPointers[i + 1];
          for (int k = _columnPointers[i]; k < high; k++) {
            valA[0] = valuesA[2 * k];
            valA[1] = -valuesA[2 * k + 1];
            valY[0] = elementsY[zeroY + strideY * _rowIndexes[k]];
            valY[1] = elementsY[zeroY + strideY * _rowIndexes[k] + 1];
            sum = Complex.plus(sum, Complex.multiply(valA, valY));
          }
          sum = Complex.multiply(alpha, sum);
          valZ[0] = elementsZ[zidx];
          valZ[1] = elementsZ[zidx + 1];
          valZ = Complex.multiply(valZ, beta);
          elementsZ[zidx] = sum[0] + valZ[0];
          elementsZ[zidx + 1] = sum[1] + valZ[1];
          zidx += strideZ;
        }
      //}
    }
    return z;
  }

  ComplexMatrix multiply(ComplexMatrix B, ComplexMatrix C, [Float64List alpha = null, Float64List beta = null, bool transposeA = false, bool transposeB = false]) {
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
      if (B is SparseCCComplexMatrix) {
        C = new SparseCCComplexMatrix.sized(rowsA, p, rowsA * p);
      } else {
        C = new DenseComplexMatrix(rowsA, p);
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

    if (!ignore && !(beta[0] == 1 && beta[1] == 0)) {
      C.forEach(cfunc.multiply(beta));
    }

    if ((B is DenseComplexMatrix) && (C is DenseComplexMatrix)) {
      SparseCCComplexMatrix AA;
      if (transposeA) {
        AA = conjugateTranspose();
      } else {
        AA = this;
      }
      DenseComplexMatrix BB;
      if (transposeB) {
        BB = B.conjugateTranspose as DenseComplexMatrix;
      } else {
        BB = B;
      }
      DenseComplexMatrix CC = C;
      Int32List columnPointersA = AA._columnPointers;
      Int32List rowIndexesA = AA._rowIndexes;
      Float64List valuesA = AA._values;

      int zeroB = BB.index(0, 0);
      int rowStrideB = BB.rowStride;
      int columnStrideB = BB.columnStride;
      Float64List elementsB = BB._elements;

      int zeroC = CC.index(0, 0);
      int rowStrideC = CC.rowStride;
      int columnStrideC = CC.columnStride;
      Float64List elementsC = CC._elements;
      Float64List valA = new Float64List(2);
      Float64List valB = new Float64List(2);
      for (int jj = 0; jj < columnsB; jj++) {
        for (int kk = 0; kk < columnsA; kk++) {
          int high = columnPointersA[kk + 1];
          valB[0] = elementsB[zeroB + kk * rowStrideB + jj * columnStrideB];
          valB[1] = elementsB[zeroB + kk * rowStrideB + jj * columnStrideB + 1];
          for (int ii = columnPointersA[kk]; ii < high; ii++) {
            int j = rowIndexesA[ii];
            valA[0] = valuesA[2 * ii];
            valA[1] = valuesA[2 * ii + 1];
            valA = Complex.multiply(valA, valB);
            elementsC[zeroC + j * rowStrideC + jj * columnStrideC] += valA[0];
            elementsC[zeroC + j * rowStrideC + jj * columnStrideC + 1] += valA[1];
          }
        }
      }
      if (!(alpha[0] == 1.0 && alpha[1] == 0)) {
        C.forEach(cfunc.multiply(alpha));
      }

    } else if ((B is SparseCCComplexMatrix) && (C is SparseCCComplexMatrix)) {
      SparseCCComplexMatrix AA;
      if (transposeA) {
        AA = conjugateTranspose();
      } else {
        AA = this;
      }
      SparseCCComplexMatrix BB = B;
      if (transposeB) {
        BB = BB.conjugateTranspose();
      }
      SparseCCComplexMatrix CC = C;
      int j,
          nz = 0,
          m,
          n;
      Int32List Cp, Ci, Bp, w, Bi;
      Float64List x, Bx, Cx;
      m = rowsA;
      n = columnsB;
      Bp = BB._columnPointers;
      Bi = BB._rowIndexes;
      Bx = BB._values;
      w = new Int32List(m);
      /* get workspace */
      x = new Float64List(2 * m);
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
          Float64List valuesNew = new Float64List(2 * nzmaxC);
          //System.arraycopy(Cx, 0, valuesNew, 0, Cx.length);
          valuesNew.setAll(0, Cx);
          Cx = valuesNew;
        }
        Cp[j] = nz;
        /* column j of C starts here */
        Float64List elemB = new Float64List(2);
        for (p = Bp[j]; p < Bp[j + 1]; p++) {
          elemB[0] = Bx[2 * p];
          elemB[1] = Bx[2 * p + 1];
          nz = _scatter(AA, Bi[p], elemB, w, x, j + 1, CC, nz);
        }
        for (p = Cp[j]; p < nz; p++) {
          Cx[2 * p] = x[2 * Ci[p]];
          Cx[2 * p + 1] = x[2 * Ci[p] + 1];
        }
      }
      Cp[n] = nz;
      /* finalize the last column of C */
      if (!(alpha[0] == 1.0 && alpha[1] == 0)) {
        CC.forEach(cfunc.multiply(alpha));
      }
    } else {
      if (transposeB) {
        B = B.conjugateTranspose();
      }
      // cache views
      final List<ComplexVector> Brows = new List<ComplexVector>(columnsA);
      for (int i = columnsA; --i >= 0; ) {
        Brows[i] = B.row(i);
      }
      final List<ComplexVector> Crows = new List<ComplexVector>(rowsA);
      for (int i = rowsA; --i >= 0; ) {
        Crows[i] = C.row(i);
      }

      final cfunc.ComplexPlusMultSecond fun = cfunc.ComplexPlusMultSecond.plusMult(new Float64List(2));

      final Int32List rowIndexesA = _rowIndexes;
      final Int32List columnPointersA = _columnPointers;
      final Float64List valuesA = _values;
      Float64List valA = new Float64List(2);
      for (int i = _columns; --i >= 0; ) {
        int low = columnPointersA[i];
        for (int k = columnPointersA[i + 1]; --k >= low; ) {
          int j = rowIndexesA[k];
          valA[0] = valuesA[2 * k];
          valA[1] = valuesA[2 * k + 1];
          fun.multiplicator = Complex.multiply(valA, alpha);
          if (!transposeA) {
            Crows[j].forEachVector(Brows[i], fun);
          } else {
            Crows[i].forEachVector(Brows[j], fun);
          }
        }
      }
    }
    return C;
  }

  ComplexMatrix _getContent() {
    return this;
  }

  void _insert(int row, int column, int index, Float64List value) {
    List<int> rowIndexesList = new List<int>.from(_rowIndexes);
    //rowIndexesList._setSizeRaw(_columnPointers[_columns]);
    List<double> valuesList = new List<double>.from(_values);
    //valuesList._setSizeRaw(2 * _columnPointers[_columns]);
    rowIndexesList.insert(index, row);
    valuesList.insert(2 * index, value[0]);
    valuesList.insert(2 * index + 1, value[1]);
    for (int i = _columnPointers.length; --i > column; ) {
      _columnPointers[i]++;
    }
    _rowIndexes = new Int32List.fromList(rowIndexesList);
    _values = new Float64List.fromList(valuesList);
  }

  void _insertParts(int row, int column, int index, double re, double im) {
    List<int> rowIndexesList = new List<int>.from(_rowIndexes);
    //rowIndexesList._setSizeRaw(_columnPointers[_columns]);
    List<double> valuesList = new List<double>.from(_values);
    //valuesList._setSizeRaw(2 * _columnPointers[_columns]);
    rowIndexesList.insert(index, row);
    valuesList.insert(2 * index, re);
    valuesList.insert(2 * index + 1, im);
    for (int i = _columnPointers.length; --i > column; ) {
      _columnPointers[i]++;
    }
    _rowIndexes = new Int32List.fromList(rowIndexesList);
    _values = new Float64List.fromList(valuesList);
  }

  void _remove(int column, int index) {
    List<int> rowIndexesList = new List<int>.from(_rowIndexes);
    List<double> valuesList = new List<double>.from(_values);
    rowIndexesList.remove(index);
    valuesList.remove(2 * index);
    valuesList.remove(2 * index + 1);
    for (int i = _columnPointers.length; --i > column; ) {
      _columnPointers[i]--;
    }
    _rowIndexes = new Int32List.fromList(rowIndexesList);
    _values = new Float64List.fromList(valuesList);
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
    Float64List valuesNew = new Float64List(2 * nzmax);
    length = Math.min(nzmax, _values.length);
    //System.arraycopy(_values, 0, valuesNew, 0, length);
    valuesNew.setAll(0, _values);
    _values = valuesNew;
  }

  int _scatter(SparseCCComplexMatrix A, int j, Float64List beta, Int32List w, Float64List x, int mark, SparseCCComplexMatrix C, int nz) {
    int i, p;
    Int32List Ap, Ai, Ci;
    Float64List Ax;
    Ap = A._columnPointers;
    Ai = A._rowIndexes;
    Ax = A._values;
    Ci = C._rowIndexes;
    Float64List valX = new Float64List(2);
    Float64List valA = new Float64List(2);
    for (p = Ap[j]; p < Ap[j + 1]; p++) {
      i = Ai[p];
      /* A(i,j) is nonzero */
      if (w[i] < mark) {
        w[i] = mark;
        /* i is new entry in column j */
        Ci[nz++] = i;
        /* add i to pattern of C(:,j) */
        if (x != null) {
          valA[0] = Ax[2 * p];
          valA[1] = Ax[2 * p + 1];
          valA = Complex.multiply(beta, valA);
          x[2 * i] = valA[0];
          /* x(i) = beta*A(i,j) */
          x[2 * i + 1] = valA[1];
        }
      } else if (x != null) {
        valA[0] = Ax[2 * p];
        valA[1] = Ax[2 * p + 1];
        valA = Complex.multiply(beta, valA);
        x[2 * i] += valA[0];
        /* i exists in C(:,j) already */
        x[2 * i + 1] += valA[1];
      }
    }
    return nz;
  }

  Object clone() {
    return new SparseCCComplexMatrix(_rows, _columns, _rowIndexes, _columnPointers, _values);
  }
}
