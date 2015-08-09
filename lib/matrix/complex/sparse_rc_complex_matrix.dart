// Copyright (C) 1999 CERN - European Organization for Nuclear Research.
//
// Permission to use, copy, modify, distribute and sell this software and
// its documentation for any purpose is hereby granted without fee, provided
// that the above copyright notice appear in all copies and that both that
// copyright notice and this permission notice appear in supporting
// documentation.
//
// CERN makes no representations about the suitability of this software for
// any purpose. It is provided "as is" without expressed or implied warranty.
part of cern.colt.matrix.complex;

/// Sparse row-compressed 2-d matrix holding `complex` elements.
///
/// Internally uses the standard sparse row-compressed format.
///
/// Cells that:
/// - are never set to non-zero values do not use any memory.
/// - switch from zero to non-zero state do use memory.
/// - switch back from non-zero to zero state also do use memory. Their memory
/// is not automatically reclaimed.
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

  Int32List _rowPointers;

  Int32List _columnIndexes;

  Float64List _values;

  /// Constructs a matrix with a given number of rows and columns. All entries
  /// are initially `0`.
  factory SparseRCComplexMatrix(int rows, int columns, [int nzmax = null]) {
    if (nzmax == null) {
      nzmax = Math.min(10 * rows, MAX_INT);
    }
    final columnIndexes = new Int32List(nzmax);
    final values = new Float64List(2 * nzmax);
    final rowPointers = new Int32List(rows + 1);
    return new SparseRCComplexMatrix._internal(
        rows, columns, rowPointers, columnIndexes, values);
  }

  /// Constructs a matrix with indexes given in the coordinate format and
  /// single value.
  factory SparseRCComplexMatrix.withValue(int rows, int columns,
      Int32List rowIndexes, Int32List columnIndexes, double re, double im,
      bool removeDuplicates) {
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
    final m = new SparseRCComplexMatrix._internal(
        rows, columns, _rowPointers, _columnIndexes, _values);
    if (removeDuplicates) {
      m.removeDuplicates();
    }
    return m;
  }

  /// Constructs a matrix with indexes and values given in the coordinate
  /// format.
  factory SparseRCComplexMatrix.withValues(int rows, int columns,
      Int32List rowIndexes, Int32List columnIndexes, Float64List values,
      {bool removeDuplicates: false, bool removeZeroes: false}) {
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
    final m = new SparseRCComplexMatrix._internal(
        rows, columns, _rowPointers, _columnIndexes, _values);
    if (removeZeroes) {
      m.removeZeroes();
    }
    if (removeDuplicates) {
      m.removeDuplicates();
    }
    return m;
  }

  SparseRCComplexMatrix._internal(int rows, int columns, Int32List rowPointers,
      Int32List columnIndexes, Float64List values)
      : super._(rows, columns) {
    if (rowPointers.length != rows + 1) {
      throw new ArgumentError("rowPointers.length != rows + 1");
    }
    if (2 * columnIndexes.length != values.length) {
      throw new ArgumentError("2 * columnIndexes.length != values.length");
    }
    _rowPointers = rowPointers;
    _columnIndexes = columnIndexes;
    _values = values;
  }

  static SparseRCComplexMatrix create(int rows, int columns) {
    return new SparseRCComplexMatrix(rows, columns);
  }

  void apply(final cfunc.ComplexComplexFunction fn) {
    if (fn is cfunc.ComplexMult) {
      // x[i] = mult*x[i]
      Complex alpha = fn.multiplicator;
      if (alpha.real == 1 && alpha.imaginary == 0) {
        return;
      }
      if (alpha.real == 0 && alpha.imaginary == 0) {
        fill(alpha.real, alpha.imaginary);
        return;
      }
      if (alpha.real != alpha.real || alpha.imaginary != alpha.imaginary) {
        fill(alpha.real, alpha.imaginary); // isNaN(). This should not happen.
        return;
      }
      int nz = cardinality;
      for (int j = 0; j < nz; j++) {
        var elem = new Complex(_values[2 * j], _values[2 * j + 1]);
        elem = elem * alpha;
        _values[2 * j] = elem.real;
        _values[2 * j + 1] = elem.imaginary;
      }
    } else {
      forEachNonZero((int i, int j, Complex value) {
        return fn(value);
      });
    }
  }

  void fill(double re, double im) {
    if (re == 0 && im == 0) {
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

  void copyFrom(ComplexMatrix source) {
    if (source == this) {
      return; // nothing to do
    }
    checkShape(this, source);

    if (source is SparseRCComplexMatrix) {
      SparseRCComplexMatrix other = source;
      _rowPointers.setAll(0, other._rowPointers);
      int nzmax = other._columnIndexes.length;
      if (_columnIndexes.length < nzmax) {
        _columnIndexes = new Int32List(nzmax);
        _values = new Float64List(2 * nzmax);
      }
      _columnIndexes.setAll(0, other._columnIndexes);
      _values.setAll(0, other._values);
    } else if (source is SparseCCComplexMatrix) {
      SparseCCComplexMatrix other = source.conjugateTranspose();
      _rowPointers = other.columnPointers;
      _columnIndexes = other.rowIndexes;
      _values = other.values;
    } else {
      fill(0.0, 0.0);
      source.forEachNonZero((int i, int j, Complex value) {
        set(i, j, value);
        return value;
      });
    }
  }

  void assign(
      final ComplexMatrix y, cfunc.ComplexComplexComplexFunction fn) {
    checkShape(this, y);
    if ((y is SparseRCComplexMatrix) &&
        (fn == cfunc.plus || fn == cfunc.minus)) {
      // x[i] = x[i] + y[i]
      SparseRCComplexMatrix yy = y;

      Int32List rowPointersY = yy._rowPointers;
      Int32List columnIndexesY = yy._columnIndexes;
      Float64List valuesY = yy._values;

      Int32List rowPointersC = new Int32List(rows + 1);
      int cnz = Math.max(_columnIndexes.length,
          Math.min(MAX_INT, _rowPointers[rows] + rowPointersY[rows]));
      Int32List columnIndexesC = new Int32List(cnz);
      Float64List valuesC = new Float64List(2 * cnz);
      int nrow = rows;
      int ncol = columns;
      int nzmax = cnz;
      if (fn == cfunc.plus) {
        // x[i] = x[i] + y[i]
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
              throw new ArgumentError(
                  "The number of elements in C exceeds nzmax");
            }
          }
          rowPointersC[i + 1] = kc;
        }
        _rowPointers = rowPointersC;
        _columnIndexes = columnIndexesC;
        _values = valuesC;
        return;
      } else if (fn == cfunc.minus) {
        // x[i] = x[i] - y[i]
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
              valuesC[2 * kc] = _values[2 * ka] - valuesY[2 * kb];
              valuesC[2 * kc + 1] = _values[2 * ka + 1] - valuesY[2 * kb + 1];
              columnIndexesC[kc] = j1;
              ka++;
              kb++;
              kc++;
            } else if (j1 < j2) {
              columnIndexesC[kc] = j1;
              valuesC[2 * kc] = -_values[2 * ka];
              valuesC[2 * kc + 1] = -_values[2 * ka + 1];
              ka++;
              kc++;
            } else if (j1 > j2) {
              columnIndexesC[kc] = j2;
              valuesC[2 * kc] = -valuesY[2 * kb];
              valuesC[2 * kc + 1] = -valuesY[2 * kb + 1];
              kb++;
              kc++;
            }
            if (kc >= nzmax) {
              throw new ArgumentError(
                  "The number of elements in C exceeds nzmax");
            }
          }
          rowPointersC[i + 1] = kc;
        }
        _rowPointers = rowPointersC;
        _columnIndexes = columnIndexesC;
        _values = valuesC;
        return;
      }
    }

    if (fn is cfunc.ComplexPlusMultSecond) {
      // x[i] = x[i] + alpha*y[i]
      Complex alpha = fn.multiplicator;
      if (alpha.real == 0 && alpha.imaginary == 0) {
        return; // nothing to do
      }
      y.forEachNonZero((int i, int j, Complex value) {
        set(i, j, get(i, j) + (alpha * value));
        return value;
      });
      return;
    }

    if (fn is cfunc.ComplexPlusMultFirst) {
      // x[i] = alpha*x[i] + y[i]
      Complex alpha = fn.multiplicator;
      if (alpha.real == 0 && alpha.imaginary == 0) {
        copyFrom(y);
        return;
      }
      y.forEachNonZero((int i, int j, Complex value) {
        set(i, j, (alpha * get(i, j)) + value);
        return value;
      });
      return;
    }

    if (fn == cfunc.mult) {
      // x[i] = x[i] * y[i]
      for (int i = 0; i < rows; i++) {
        int high = _rowPointers[i + 1];
        for (int k = _rowPointers[i]; k < high; k++) {
          int j = _columnIndexes[k];
          var elem = new Complex(_values[2 * k], _values[2 * k + 1]);
          elem = elem * y.get(i, j);
          _values[2 * k] = elem.real;
          _values[2 * k + 1] = elem.imaginary;
          //if (_values[2 * k] == 0 && _values[2 * k + 1] == 0) _remove(i, j);
        }
      }
      return;
    }

    if (fn == cfunc.div) {
      // x[i] = x[i] / y[i]
      for (int i = 0; i < rows; i++) {
        int high = _rowPointers[i + 1];
        for (int k = _rowPointers[i]; k < high; k++) {
          int j = _columnIndexes[k];
          var elem = new Complex(_values[2 * k], _values[2 * k + 1]);
          elem = elem / y.get(i, j);
          _values[2 * k] = elem.real;
          _values[2 * k + 1] = elem.imaginary;
          //if (_values[2 * k] == 0 && _values[2 * k + 1] == 0) _remove(i, j);
        }
      }
      return;
    }
    super.assign(y, fn);
  }

  int get cardinality => _rowPointers[rows];

  void forEachNonZero(final cfunc.IntIntComplexFunction fn) {
    for (int i = 0; i < rows; i++) {
      int high = _rowPointers[i + 1];
      for (int k = _rowPointers[i]; k < high; k++) {
        int j = _columnIndexes[k];
        var value = new Complex(_values[2 * k], _values[2 * k + 1]);
        Complex r = fn(i, j, value);
        if (r.real != value.real || r.imaginary != value.imaginary) {
          _values[2 * k] = r.real;
          _values[2 * k + 1] = r.imaginary;
        }
      }
    }
  }

  /// Returns a new matrix that has the same elements as this matrix, but is in
  /// a column-compressed form.
  SparseCCComplexMatrix columnCompressed() {
    SparseRCComplexMatrix tr = conjugateTranspose();
    SparseCCComplexMatrix cc = new SparseCCComplexMatrix(rows, columns);
    cc._rowIndexes = tr._columnIndexes;
    cc._columnPointers = tr._rowPointers;
    cc._values = tr._values;
    return cc;
  }

  Int32List get columnIndexes => _columnIndexes;

  /// Returns a new matrix that has the same elements as this matrix, but is in
  /// a dense form.
  DenseComplexMatrix dense() {
    var dense = new DenseComplexMatrix(rows, columns);
    forEachNonZero((int i, int j, Complex value) {
      dense.set(i, j, get(i, j));
      return value;
    });
    return dense;
  }

  Complex get(int row, int column) {
    //int k = Sorting.binarySearchFromTo(columnIndexes, column, rowPointers[row], rowPointers[row + 1] - 1);
    int k = _searchFromTo(
        _columnIndexes, column, _rowPointers[row], _rowPointers[row + 1] - 1);

    if (k >= 0) {
      return new Complex(_values[2 * k], _values[2 * k + 1]);
    }
    return Complex.ZERO;
  }

  Int32List get rowPointers => _rowPointers;

  /// Returns a new matrix that is the transpose of this matrix.
  SparseRCComplexMatrix transpose() {
    int nnz = _rowPointers[rows];
    var w = new Int32List(columns);
    var rowPointersT = new Int32List(columns + 1);
    var columnIndexesT = new Int32List(nnz);
    var valuesT = new Float64List(2 * nnz);

    for (int p = 0; p < nnz; p++) {
      w[_columnIndexes[p]]++;
    }
    _cumsum(rowPointersT, w, columns);
    int q;
    for (int j = 0; j < rows; j++) {
      int high = _rowPointers[j + 1];
      for (int p = _rowPointers[j]; p < high; p++) {
        columnIndexesT[q = w[_columnIndexes[p]]++] = j;
        valuesT[2 * q] = _values[2 * p];
        valuesT[2 * q + 1] = _values[2 * p + 1];
      }
    }
    var T = new SparseRCComplexMatrix(columns, rows);
    T._rowPointers = rowPointersT;
    T._columnIndexes = columnIndexesT;
    T._values = valuesT;
    return T;
  }

  /// Returns a new matrix that is the conjugate transpose of this matrix.
  SparseRCComplexMatrix conjugateTranspose() {
    int nnz = _rowPointers[rows];
    var w = new Int32List(columns);
    var rowPointersT = new Int32List(columns + 1);
    var columnIndexesT = new Int32List(nnz);
    var valuesT = new Float64List(2 * nnz);

    for (int p = 0; p < nnz; p++) {
      w[_columnIndexes[p]]++;
    }
    _cumsum(rowPointersT, w, columns);
    int q;
    for (int j = 0; j < rows; j++) {
      int high = _rowPointers[j + 1];
      for (int p = _rowPointers[j]; p < high; p++) {
        columnIndexesT[q = w[_columnIndexes[p]]++] = j;
        valuesT[2 * q] = _values[2 * p];
        valuesT[2 * q + 1] = -_values[2 * p + 1];
      }
    }
    var T = new SparseRCComplexMatrix(columns, rows);
    T._rowPointers = rowPointersT;
    T._columnIndexes = columnIndexesT;
    T._values = valuesT;
    return T;
  }

  Float64List get values => _values;

  ComplexMatrix like2D(int rows, int columns) {
    return new SparseRCComplexMatrix(rows, columns); //, _columnIndexes.length);
  }

  ComplexVector like1D(int size) => new SparseComplexVector(size);

  /// Removes (sums) duplicate entries (if any).
  void removeDuplicates() {
    int nz = 0;
    var w = new Int32List(columns);
    for (var i = 0; i < columns; i++) {
      w[i] = -1;
    }
    for (int j = 0; j < rows; j++) {
      var q = nz;
      for (int p = _rowPointers[j]; p < _rowPointers[j + 1]; p++) {
        var i = _columnIndexes[p];
        if (w[i] >= q) {
          _values[2 * w[i]] += _values[2 * p];
          _values[2 * w[i] + 1] += _values[2 * p + 1];
        } else {
          w[i] = nz;
          _columnIndexes[nz] = i;
          _values[2 * nz] = _values[2 * p];
          _values[2 * nz + 1] = _values[2 * p + 1];
          nz++;
        }
      }
      _rowPointers[j] = q;
    }
    _rowPointers[rows] = nz;
  }

  /// Removes zero entries (if any).
  void removeZeroes() {
    int nz = 0;
    double eps = Math.pow(2, -52);
    for (int j = 0; j < rows; j++) {
      int p = _rowPointers[j];
      _rowPointers[j] = nz;
      for (; p < _rowPointers[j + 1]; p++) {
        var elem = new Complex(_values[2 * p], _values[2 * p + 1]);
        if (elem.abs() > eps) {
          _values[2 * nz] = _values[2 * p];
          _values[2 * nz + 1] = _values[2 * p + 1];
          _columnIndexes[nz++] = _columnIndexes[p];
        }
      }
    }
    _rowPointers[rows] = nz;
  }

  void set(int row, int column, Complex value) {
    //int k = Sorting.binarySearchFromTo(columnIndexes, column, rowPointers[row], rowPointers[row + 1] - 1);
    int k = _searchFromTo(
        _columnIndexes, column, _rowPointers[row], _rowPointers[row + 1] - 1);

    if (k >= 0) {
      // found
      //if (value[0] == 0 && value[1] == 0) _remove(row, k); else {
      _values[2 * k] = value.real;
      _values[2 * k + 1] = value.imaginary;
      //}
      return;
    }

    if (value.real != 0 || value.imaginary != 0) {
      k = -k - 1;
      _insert(row, column, k, value);
    }
  }

  void setParts(int row, int column, double re, double im) {
    //int k = Sorting.binarySearchFromTo(columnIndexes, column, rowPointers[row], rowPointers[row + 1] - 1);
    int k = _searchFromTo(
        _columnIndexes, column, _rowPointers[row], _rowPointers[row + 1] - 1);

    if (k >= 0) {
      // found
      //if (re == 0 && im == 0) _remove(row, k); else {
      _values[2 * k] = re;
      _values[2 * k + 1] = im;
      //}
      return;
    }

    if (re != 0 || im != 0) {
      k = -k - 1;
      _insertParts(row, column, k, re, im);
    }
  }

  String toString() {
    var buf = new StringBuffer();
    buf.write("$rows x $columns sparse matrix, nnz = $cardinality\n");
    for (int i = 0; i < rows; i++) {
      int high = _rowPointers[i + 1];
      for (int j = _rowPointers[i]; j < high; j++) {
        if (_values[2 * j + 1] > 0) {
          buf
            ..write('(')
            ..write(i)
            ..write(',')
            ..write(_columnIndexes[j])
            ..write(')')
            ..write('\t')
            ..write(_values[2 * j])
            ..write('+')
            ..write(_values[2 * j + 1])
            ..write('i')
            ..write('\n');
        } else if (_values[2 * j + 1] == 0) {
          buf
            ..write('(')
            ..write(i)
            ..write(',')
            ..write(_columnIndexes[j])
            ..write(')')
            ..write('\t')
            ..write(_values[2 * j])
            ..write('\n');
        } else {
          buf
            ..write('(')
            ..write(i)
            ..write(',')
            ..write(_columnIndexes[j])
            ..write(')')
            ..write('\t')
            ..write(_values[2 * j])
            ..write('-')
            ..write(_values[2 * j + 1])
            ..write('i')
            ..write('\n');
        }
      }
    }
    return buf.toString();
  }

  void trimToSize() => _realloc(0);

  ComplexVector mult(ComplexVector y,
      [ComplexVector z = null, Complex alpha = null,
      Complex beta = null, bool transposeA = false]) {
    if (alpha == null) {
      alpha = Complex.ONE;
    }
    if (beta == null) {
      beta = (z == null ? Complex.ONE : Complex.ZERO);
    }
    final int rowsA = transposeA ? columns : rows;
    final int columnsA = transposeA ? rows : columns;

    bool ignore = (z == null || !transposeA);
    if (z == null) {
      z = new DenseComplexVector(rowsA);
    }

    if (!(y is DenseComplexVector && z is DenseComplexVector)) {
      return super.mult(y, z, alpha, beta, transposeA);
    }

    if (columnsA != y.size || rowsA > z.size) {
      throw new ArgumentError("Incompatible args: " +
          ((transposeA ? dice() : this).toStringShort()) +
          ", " +
          y.toStringShort() +
          ", " +
          z.toStringShort());
    }

    DenseComplexVector zz = z as DenseComplexVector;
    Float64List elementsZ = zz._elements;
    int strideZ = zz.stride;
    int zeroZ = z.index(0);

    DenseComplexVector yy = y as DenseComplexVector;
    Float64List elementsY = yy._elements;
    int strideY = yy.stride;
    int zeroY = y.index(0);

    if (transposeA) {
      if ((!ignore) && !((beta.real == 1) && (beta.imaginary == 0))) {
        z.apply(cfunc.multiply(beta));
      }

      for (int i = 0; i < rows; i++) {
        int high = _rowPointers[i + 1];
        var yElem = new Complex(
            elementsY[zeroY + strideY * i], elementsY[zeroY + strideY * i + 1]);
        yElem = alpha * yElem;
        for (int k = _rowPointers[i]; k < high; k++) {
          int j = _columnIndexes[k];
          var val = new Complex(_values[2 * k], -_values[2 * k + 1]);
          val = val * yElem;
          elementsZ[zeroZ + strideZ * j] += val.real;
          elementsZ[zeroZ + strideZ * j + 1] += val.imaginary;
        }
      }

      return z;
    }

    int zidx = zeroZ;
    if (beta.real == 0.0 && beta.imaginary == 0) {
      for (int i = 0; i < rows; i++) {
        var sum = Complex.ZERO;
        int high = _rowPointers[i + 1];
        for (int k = _rowPointers[i]; k < high; k++) {
          var yElem = new Complex(
              elementsY[zeroY + strideY * _columnIndexes[k]],
              elementsY[zeroY + strideY * _columnIndexes[k] + 1]);
          var val = new Complex(_values[2 * k], _values[2 * k + 1]);
          sum += val * yElem;
        }
        sum = alpha * sum;
        elementsZ[zidx] = sum.real;
        elementsZ[zidx + 1] = sum.imaginary;
        zidx += strideZ;
      }
    } else {
      for (int i = 0; i < rows; i++) {
        var sum = Complex.ZERO;
        int high = _rowPointers[i + 1];
        for (int k = _rowPointers[i]; k < high; k++) {
          var yElem = new Complex(
              elementsY[zeroY + strideY * _columnIndexes[k]],
              elementsY[zeroY + strideY * _columnIndexes[k] + 1]);
          var val = new Complex(_values[2 * k], _values[2 * k + 1]);
          sum += val * yElem;
        }
        sum = alpha * sum;
        var zElem = new Complex(elementsZ[zidx], elementsZ[zidx + 1]);
        zElem = beta * zElem;
        elementsZ[zidx] = sum.real + zElem.real;
        elementsZ[zidx + 1] = sum.imaginary + zElem.imaginary;
        zidx += strideZ;
      }
    }
    return z;
  }

  ComplexMatrix multiply(ComplexMatrix B,
      [ComplexMatrix C = null, Complex alpha = null,
      Complex beta = null, bool transposeA = false, bool transposeB = false]) {
    if (alpha == null) {
      alpha = Complex.ONE;
    }
    if (beta == null) {
      beta = (C == null ? Complex.ONE : Complex.ZERO);
    }
    int rowsA = rows;
    int columnsA = columns;
    if (transposeA) {
      rowsA = columns;
      columnsA = rows;
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
        C = new DenseComplexMatrix(rowsA, p);
      }
    }

    if (rowsB != columnsA) {
      throw new ArgumentError("Matrix inner dimensions must agree:" +
          toStringShort() +
          ", " +
          (transposeB ? B.dice() : B).toStringShort());
    }
    if (C.rows != rowsA || C.columns != p) {
      throw new ArgumentError("Incompatible result matrix: " +
          toStringShort() +
          ", " +
          (transposeB ? B.dice() : B).toStringShort() +
          ", " +
          C.toStringShort());
    }
    if (this == C || B == C) {
      throw new ArgumentError("Matrices must not be identical");
    }

    if (!ignore && !(beta.real == 1.0 && beta.imaginary == 0)) {
      C.apply(cfunc.multiply(beta));
    }

    if ((B is DenseComplexMatrix) && (C is DenseComplexMatrix)) {
      SparseRCComplexMatrix AA;
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
      Int32List rowPointersA = AA._rowPointers;
      Int32List columnIndexesA = AA._columnIndexes;
      Float64List valuesA = AA._values;
      for (int ii = 0; ii < rowsA; ii++) {
        int highA = rowPointersA[ii + 1];
        for (int ka = rowPointersA[ii]; ka < highA; ka++) {
          var valA = new Complex(valuesA[2 * ka], valuesA[2 * ka + 1]);
          var scal = alpha * valA;
          int jj = columnIndexesA[ka];
          CC.row(ii).assign(BB.row(jj), cfunc.plusMultSecond(scal));
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

      var iw = new Int32List(columnsB + 1);
      for (int i = 0; i < iw.length; i++) {
        iw[i] = -1;
      }
      int len = -1;
      for (int ii = 0; ii < rowsA; ii++) {
        int highA = rowPointersA[ii + 1];
        for (int ka = rowPointersA[ii]; ka < highA; ka++) {
          var valA = new Complex(valuesA[2 * ka], valuesA[2 * ka + 1]);
          var scal = alpha * valA;
          int jj = columnIndexesA[ka];
          int highB = rowPointersB[jj + 1];
          for (int kb = rowPointersB[jj]; kb < highB; kb++) {
            int jcol = columnIndexesB[kb];
            int jpos = iw[jcol];
            if (jpos == -1) {
              len++;
              if (len >= nzmax) {
                throw new ArgumentError(
                    "The max number of nonzero elements in C is too small.");
              }
              columnIndexesC[len] = jcol;
              iw[jcol] = len;
              var valB = new Complex(valuesB[2 * kb], valuesB[2 * kb + 1]);
              valB = scal * valB;
              valuesC[2 * len] = valB.real;
              valuesC[2 * len + 1] = valB.imaginary;
            } else {
              var valB = new Complex(valuesB[2 * kb], valuesB[2 * kb + 1]);
              valB = scal * valB;
              valuesC[2 * jpos] += valB.real;
              valuesC[2 * jpos + 1] += valB.imaginary;
            }
          }
        }
        for (int k = rowPointersC[ii]; k < len + 1; k++) {
          iw[columnIndexesC[k]] = -1;
        }
        rowPointersC[ii + 1] = len + 1;

        //int length = rowPointersC[ii + 1] - rowPointersC[ii];
        //IntVector columnIndexesCPart = columnIndexesC.viewPart(rowPointersC[ii], length);
        //Int32List indexes = cern.colt.matrix.tint.algo.IntSorting.quickSort.sortIndex(columnIndexesCPart);
        //Arrays.sort(columnIndexesCElements, rowPointersC[ii], rowPointersC[ii + 1]);
        //ComplexVector valuesCPart = valuesC.viewPart(rowPointersC[ii], length).viewSelection(indexes);
        //valuesC.viewPart(rowPointersC[ii], length).assign(valuesCPart);
      }
      //CC.columnIndexes.elements((Int32List) columnIndexesC.elements);
      //CC.columnIndexes.setSize(columnIndexesSize);
      //CC.values.elements((Float64List) valuesC.elements);
      //CC.values.setSize(columnIndexesSize);
    } else {
      if (transposeB) {
        B = B.conjugateTranspose();
      }
      // cache views
      final List<ComplexVector> Brows =
          new List<ComplexVector>(columnsA);
      for (int i = columnsA; --i >= 0;) {
        Brows[i] = B.row(i);
      }
      final List<ComplexVector> Crows =
          new List<ComplexVector>(rowsA);
      for (int i = rowsA; --i >= 0;) {
        Crows[i] = C.row(i);
      }

      var fn = cfunc.ComplexPlusMultSecond.plusMult(Complex.ZERO);

      Int32List columnIndexesA = _columnIndexes;
      Float64List valuesA = _values;
      for (int i = rows; --i >= 0;) {
        int low = _rowPointers[i];
        for (int k = _rowPointers[i + 1]; --k >= low;) {
          int j = columnIndexesA[k];
          var valA = new Complex(valuesA[2 * k], valuesA[2 * k + 1]);
          fn.multiplicator = valA * alpha;
          if (!transposeA) {
            Crows[i].assign(Brows[j], fn);
          } else {
            Crows[j].assign(Brows[i], fn);
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
    return nz2;
  }

  void _realloc(int nzmax) {
    if (nzmax <= 0) {
      nzmax = _rowPointers[rows];
    }
    var columnIndexesNew = new Int32List(nzmax);
    columnIndexesNew.setAll(0, _columnIndexes);
    _columnIndexes = columnIndexesNew;
    var valuesNew = new Float64List(2 * nzmax);
    valuesNew.setAll(0, _values);
    _values = valuesNew;
  }

  ComplexMatrix _getContent() => this;

  void _insert(int row, int column, int index, Complex value) {
    var columnIndexesList = new List.from(_columnIndexes);
    columnIndexesList.length = _rowPointers[rows];
    var valuesList = new List.from(_values);
    valuesList.length = 2 * _rowPointers[rows];
    columnIndexesList.insert(index, column);
    valuesList.insert(2 * index, value.real);
    valuesList.insert(2 * index + 1, value.imaginary);
    for (int i = _rowPointers.length; --i > row;) {
      _rowPointers[i]++;
    }
    _columnIndexes = new Int32List.fromList(columnIndexesList);
    _values = new Float64List.fromList(valuesList);
  }

  void _insertParts(int row, int column, int index, double re, double im) {
    var columnIndexesList = new List.from(_columnIndexes);
    columnIndexesList.length = _rowPointers[rows];
    var valuesList = new List.from(_values);
    valuesList.length = 2 * _rowPointers[rows];
    columnIndexesList.insert(index, column);
    valuesList.insert(2 * index, re);
    valuesList.insert(2 * index + 1, im);
    for (int i = _rowPointers.length; --i > row;) {
      _rowPointers[i]++;
    }
    _columnIndexes = new Int32List.fromList(columnIndexesList);
    _values = new Float64List.fromList(valuesList);
  }

  /*void _remove(int row, int index) {
    var columnIndexesList = new List.from(_columnIndexes);
    columnIndexesList.length = _rowPointers[_rows];
    var valuesList = new List.from(_values);
    valuesList.length = _rowPointers[_rows];
    columnIndexesList.remove(index);
    valuesList.remove(2 * index);
    valuesList.remove(2 * index + 1);
    for (int i = _rowPointers.length; --i > row; ) {
      _rowPointers[i]--;
    }
    _columnIndexes = new Int32List.fromList(columnIndexesList);
    _values = new Float64List.fromList(valuesList);
  }*/

  Object clone() {
    return new SparseRCComplexMatrix._internal(
        rows, columns, _rowPointers, _columnIndexes, _values);
  }

  double normInf() {
    double norm = 0.0;
    for (int j = 0; j < rows; j++) {
      var s = 0.0;
      for (var p = _rowPointers[j]; p < _rowPointers[j + 1]; p++) {
        var elem = new Complex(_values[2 * p], _values[2 * p + 1]);
        s += elem.abs();
      }
      norm = Math.max(norm, s);
    }
    return norm;
  }

  DoubleMatrix real() {
    var vals = new Float64List(size);
    for (var i = 0; i < size; i++) {
      vals[i] = _values[2 * i];
    }
    var rowptr = new Int32List.fromList(_rowPointers);
    var colidx = new Int32List.fromList(_columnIndexes);
    return makeSparseRCDoubleMatrix(rows, columns, rowptr, colidx, vals);
  }

  DoubleMatrix imaginary() {
    var vals = new Float64List(size);
    for (var i = 0; i < size; i++) {
      vals[i] = _values[2 * i + 1];
    }
    var rowptr = new Int32List.fromList(_rowPointers);
    var colidx = new Int32List.fromList(_columnIndexes);
    return makeSparseRCDoubleMatrix(rows, columns, rowptr, colidx, vals);
  }
}
