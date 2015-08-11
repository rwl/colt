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

/// Sparse column-compressed 2-d matrix holding [Complex] elements.
///
/// Internally uses the standard sparse column-compressed format.
///
/// Cells that:
/// - are never set to non-zero values do not use any memory.
/// - switch from zero to non-zero state do use memory.
/// - switch back from non-zero to zero state also do use memory.
/// Reclamation can be triggered via [trimToSize].
class SparseCCComplexMatrix extends WrapperComplexMatrix {
  Int32List _columnPointers;

  Int32List _rowIndexes;

  Float64List _values;

  /// Constructs a matrix with a given number of rows and columns. All entries
  /// are initially `0`.
  factory SparseCCComplexMatrix(int rows, int columns, [int nzmax = null]) {
    if (nzmax == null) {
      nzmax = Math.min(10 * rows, MAX_INT);
    }
    var rowIndexes = new Int32List(nzmax);
    var values = new Float64List(2 * nzmax);
    var columnPointers = new Int32List(columns + 1);
    return new SparseCCComplexMatrix._internal(
        rows, columns, rowIndexes, columnPointers, values);
  }

  SparseCCComplexMatrix._internal(int rows, int columns, Int32List rowIndexes,
      Int32List columnPointers, Float64List values)
      : super._(rows, columns) {
    if (columnPointers.length != columns + 1) {
      throw new ArgumentError("columnPointers.length != columns + 1");
    }
    if (2 * rowIndexes.length != values.length) {
      throw new ArgumentError("2 * rowIndexes.length != values.length");
    }
    _columnPointers = columnPointers;
    _rowIndexes = rowIndexes;
    _values = values;
  }

  /// Constructs a matrix with indexes given in the coordinate format and a
  /// single value.
  factory SparseCCComplexMatrix.withValue(int rows, int columns,
      Int32List rowIndexes, Int32List columnIndexes, double re, double im,
      {bool removeDuplicates: false, bool sortRowIndexes: false}) {
    if (rowIndexes.length != columnIndexes.length) {
      throw new ArgumentError("rowIndexes.length != columnIndexes.length");
    }
    if (re == 0 && im == 0) {
      throw new ArgumentError("value cannot be 0");
    }

    int nz = Math.max(rowIndexes.length, 1);
    var _rowIndexes = new Int32List(nz);
    var _values = new Float64List(2 * nz);
    var _columnPointers = new Int32List(columns + 1);
    var columnCounts = new Int32List(columns);
    for (int k = 0; k < nz; k++) {
      columnCounts[columnIndexes[k]]++;
    }
    cumsum(_columnPointers, columnCounts);
    for (int k = 0; k < nz; k++) {
      var r = columnCounts[columnIndexes[k]]++;
      _rowIndexes[r] = rowIndexes[k];
      _values[2 * r] = re;
      _values[2 * r + 1] = im;
    }
    var m = new SparseCCComplexMatrix._internal(
        rows, columns, _rowIndexes, _columnPointers, _values);
    if (removeDuplicates) {
      m.removeDuplicates();
    }
    if (sortRowIndexes) {
      m.sortRowIndexes();
    }
    return m;
  }

  /// Constructs a matrix with indexes and values given in the coordinate
  /// format.
  factory SparseCCComplexMatrix.withValues(int rows, int columns,
      Int32List rowIndexes, Int32List columnIndexes, Float64List values,
      {bool removeDuplicates: false, bool removeZeroes: false,
      bool sortRowIndexes: false}) {
    if (rowIndexes.length != columnIndexes.length) {
      throw new ArgumentError("rowIndexes.length != columnIndexes.length");
    } else if (2 * rowIndexes.length != values.length) {
      throw new ArgumentError("2 * rowIndexes.length != values.length");
    }
    int nz = Math.max(rowIndexes.length, 1);
    var _rowIndexes = new Int32List(nz);
    var _values = new Float64List(2 * nz);
    var _columnPointers = new Int32List(columns + 1);
    var columnCounts = new Int32List(columns);
    for (int k = 0; k < nz; k++) {
      columnCounts[columnIndexes[k]]++;
    }
    cumsum(_columnPointers, columnCounts);
    for (int k = 0; k < nz; k++) {
      var r = columnCounts[columnIndexes[k]]++;
      _rowIndexes[r] = rowIndexes[k];
      _values[2 * r] = values[2 * k];
      _values[2 * r + 1] = values[2 * k + 1];
    }
    var m = new SparseCCComplexMatrix._internal(
        rows, columns, _rowIndexes, _columnPointers, _values);
    if (removeDuplicates) {
      m.removeDuplicates();
    }
    if (removeZeroes) {
      m.removeZeroes();
    }
    if (sortRowIndexes) {
      m.sortRowIndexes();
    }
    return m;
  }

  static SparseCCComplexMatrix create(int rows, int columns) {
    return new SparseCCComplexMatrix(rows, columns);
  }

  void apply(cfunc.ComplexComplexFunction fn) {
    if (fn is cfunc.ComplexMult) {
      // x[i] = mult*x[i]
      Complex alpha = fn.multiplicator;
      if (alpha.real == 1 && alpha.imaginary == 0) {
        return;
      }
      if (alpha.real == 0 && alpha.imaginary == 0) {
        fill(0.0, 0.0);
        return;
      }
      if (alpha.real != alpha.real || alpha.imaginary != alpha.imaginary) {
        fill(alpha.real, alpha.imaginary); // isNaN. This should not happen.
        return;
      }

      int nz = cardinality;
      for (int j = 0; j < nz; j++) {
        var value = new Complex(_values[2 * j], _values[2 * j + 1]);
        value = value * alpha;
        _values[2 * j] = value.real;
        _values[2 * j + 1] = value.imaginary;
      }
    } else {
      forEachNonZero((int i, int j, Complex value) {
        return fn(value);
      });
    }
  }

  void fill(double re, double im) {
    if (re == 0 && im == 0) {
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
  }

  void copyFrom(ComplexMatrix source) {
    if (source == this) {
      return; // nothing to do
    }
    checkShape(this, source);

    if (source is SparseCCComplexMatrix) {
      _columnPointers.setAll(0, source.columnPointers);
      int nzmax = source.rowIndexes.length;
      if (_rowIndexes.length < nzmax) {
        _rowIndexes = new Int32List(nzmax);
        _values = new Float64List(2 * nzmax);
      }
      _rowIndexes.setAll(0, source.rowIndexes);
      _values.setAll(0, source.values);
    } else if (source is SparseRCComplexMatrix) {
      SparseRCComplexMatrix tr = source.conjugateTranspose();
      _columnPointers = tr.rowPointers;
      _rowIndexes = tr.columnIndexes;
      _values = tr.values;
    } else {
      fill(0.0, 0.0);
      source.forEachNonZero((int i, int j, Complex value) {
        set(i, j, value);
        return value;
      });
    }
  }

  void assign(ComplexMatrix y, cfunc.ComplexComplexComplexFunction fn) {
    checkShape(this, y);

    if (y is SparseCCComplexMatrix && fn == cfunc.plus) {
      // x[i] = x[i] + y[i]
      SparseCCComplexMatrix yy = y;
      int nz = 0;
      var nnz = _columnPointers[columns];
      var n = yy.columns;
      var bnz = yy._columnPointers[n];
      var w = new Int32List(rows);
      var x = new Float64List(2 * rows);
      var C = new SparseCCComplexMatrix(rows, n, nnz + bnz);
      var columnPointersC = C._columnPointers;
      var rowIndexesC = C._rowIndexes;
      var valuesC = C._values;
      for (var j = 0; j < n; j++) {
        columnPointersC[j] = nz;
        nz = _toDense(this, j, Complex.ONE, w, x, j + 1, C, nz);
        nz = _toDense(yy, j, Complex.ONE, w, x, j + 1, C, nz);
        for (var p = columnPointersC[j]; p < nz; p++) {
          valuesC[2 * p] = x[2 * rowIndexesC[p]];
          valuesC[2 * p + 1] = x[2 * rowIndexesC[p] + 1];
        }
      }
      columnPointersC[n] = nz;
      _rowIndexes = rowIndexesC;
      _columnPointers = columnPointersC;
      _values = valuesC;
      return;
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
      for (int j = columns; --j >= 0;) {
        int low = _columnPointers[j];
        for (int k = _columnPointers[j + 1]; --k >= low;) {
          int i = _rowIndexes[k];
          var valA = new Complex(_values[2 * k], _values[2 * k + 1]);
          valA = valA * y.get(i, j);
          _values[2 * k] = valA.real;
          _values[2 * k + 1] = valA.imaginary;
          //if (valuesA[2 * k] == 0 && valuesA[2 * k + 1] == 0) _remove(i, j);
        }
      }
      return;
    }

    if (fn == cfunc.div) {
      // x[i] = x[i] / y[i]
      for (int j = columns; --j >= 0;) {
        int low = _columnPointers[j];
        for (int k = _columnPointers[j + 1]; --k >= low;) {
          int i = _rowIndexes[k];
          var valA = new Complex(_values[2 * k], _values[2 * k + 1]);
          valA = valA / y.get(i, j);
          _values[2 * k] = valA.real;
          _values[2 * k + 1] = valA.imaginary;
          //if (valuesA[2 * k] == 0 && valuesA[2 * k + 1] == 0) _remove(i, j);
        }
      }
      return;
    }
    super.assign(y, fn);
  }

  int get cardinality => _columnPointers[columns];

  void forEachNonZero(cfunc.IntIntComplexFunction fn) {
    for (int j = columns; --j >= 0;) {
      int low = _columnPointers[j];
      for (int k = _columnPointers[j + 1]; --k >= low;) {
        int i = _rowIndexes[k];
        var valA = new Complex(_values[2 * k], _values[2 * k + 1]);
        valA = fn(i, j, valA);
        _values[2 * k] = valA.real;
        _values[2 * k + 1] = valA.imaginary;
      }
    }
  }

  /// Column pointers (size [columns]+1).
  Int32List get columnPointers => _columnPointers;

  /// Returns a new matrix that has the same elements as this matrix, but is
  /// in a dense form.
  DenseComplexMatrix dense() {
    var dense = new DenseComplexMatrix(rows, columns);
    forEachNonZero((int i, int j, Complex value) {
      dense.set(i, j, get(i, j));
      return value;
    });
    return dense;
  }

  Complex get(int row, int column) {
    int k = searchRange(_rowIndexes, row, _columnPointers[column],
        _columnPointers[column + 1] - 1);
    if (k >= 0) {
      return new Complex(_values[2 * k], _values[2 * k + 1]);
    }
    return Complex.ZERO;
  }

  /// Returns a new matrix that has the same elements as this matrix, but is in
  /// a row-compressed form.
  SparseRCComplexMatrix rowCompressed() {
    SparseCCComplexMatrix tr = conjugateTranspose();
    var rc = new SparseRCComplexMatrix(rows, columns);
    rc._columnIndexes = tr._rowIndexes;
    rc._rowPointers = tr._columnPointers;
    rc._values = tr._values;
    return rc;
  }

  /// Row indices (size `nzmax`).
  Int32List get rowIndexes => _rowIndexes;

  dynamic get elements => _values;

  /// Returns a new matrix that is the transpose of this matrix.
  SparseCCComplexMatrix transpose() {
    var C = new SparseCCComplexMatrix(columns, rows, _rowIndexes.length);
    var rowCounts = new Int32List(rows);
    for (var p = 0; p < _columnPointers[columns]; p++) {
      rowCounts[_rowIndexes[p]]++;
    }
    cumsum(C._columnPointers, rowCounts);
    for (var j = 0; j < columns; j++) {
      for (var p = _columnPointers[j]; p < _columnPointers[j + 1]; p++) {
        var q = rowCounts[_rowIndexes[p]]++;
        C._rowIndexes[q] = j;
        C._values[q] = _values[p];
      }
    }
    return C;
  }

  SparseCCComplexMatrix conjugateTranspose() {
    var C = new SparseCCComplexMatrix(columns, rows, _rowIndexes.length);
    var rowCounts = new Int32List(rows);
    for (var p = 0; p < _columnPointers[columns]; p++) {
      rowCounts[_rowIndexes[p]]++;
    }
    cumsum(C._columnPointers, rowCounts);
    for (var j = 0; j < columns; j++) {
      for (var p = _columnPointers[j]; p < _columnPointers[j + 1]; p++) {
        var q = rowCounts[_rowIndexes[p]]++;
        C._rowIndexes[q] = j;
        C._values[2 * q] = _values[2 * p];
        C._values[2 * q + 1] = -_values[2 * p + 1];
      }
    }
    return C;
  }

  /// Numerical values (size `nzmax`).
  Float64List get values => _values;

  ComplexMatrix like2D(int rows, int columns) {
    return new SparseCCComplexMatrix(rows, columns);
  }

  ComplexVector like1D(int size) => new SparseComplexVector(size);

  void set(int row, int column, Complex value) {
    int k = searchRange(_rowIndexes, row, _columnPointers[column],
        _columnPointers[column + 1] - 1);

    if (k >= 0) {
      //if (value[0] == 0 && value[1] == 0) _remove(column, k); else
      _values[2 * k] = value.real;
      _values[2 * k + 1] = value.imaginary;
      return;
    }

    if (value.real != 0 || value.imaginary != 0) {
      k = -k - 1;
      _insert(row, column, k, value);
    }
  }

  void setParts(int row, int column, double re, double im) {
    int k = searchRange(_rowIndexes, row, _columnPointers[column],
        _columnPointers[column + 1] - 1);

    if (k >= 0) {
      //if (re == 0 && im == 0) _remove(column, k); else
      _values[2 * k] = re;
      _values[2 * k + 1] = im;
      return;
    }

    if (re != 0 || im != 0) {
      k = -k - 1;
      _insertParts(row, column, k, re, im);
    }
  }

  void sortRowIndexes() {
    SparseCCComplexMatrix tr = conjugateTranspose();
    tr = tr.conjugateTranspose();
    _columnPointers = tr._columnPointers;
    _rowIndexes = tr._rowIndexes;
    _values = tr._values;
  }

  /// Removes (sums) duplicate entries (if any).
  void removeDuplicates() {
    int nz = 0;
    var r = new Int32List(rows);
    for (var i = 0; i < rows; i++) {
      r[i] = -1;
    }
    for (var j = 0; j < columns; j++) {
      var q = nz;
      for (var p = _columnPointers[j]; p < _columnPointers[j + 1]; p++) {
        var i = _rowIndexes[p];
        if (r[i] >= q) {
          _values[r[i]] += _values[p];
        } else {
          r[i] = nz;
          _rowIndexes[nz] = i;
          _values[2 * nz] = _values[2 * p];
          _values[2 * nz + 1] = _values[2 * p + 1];
          nz++;
        }
      }
      _columnPointers[j] = q;
    }
    _columnPointers[columns] = nz;
  }

  /// Removes zero entries (if any).
  void removeZeroes() {
    int nz = 0;
    for (var j = 0; j < columns; j++) {
      var p = _columnPointers[j];
      _columnPointers[j] = nz;
      for (; p < _columnPointers[j + 1]; p++) {
        if (_values[2 * p] != 0 || _values[2 * p + 1] != 0) {
          _rowIndexes[nz] = _rowIndexes[p];
          _values[2 * nz] = _values[2 * p];
          _values[2 * nz + 1] = _values[2 * p + 1];
          nz++;
        }
      }
    }
    _columnPointers[columns] = nz;
  }

  void trimToSize() {
    var nzmax = _columnPointers[columns];
    var rowIndexesNew = new Int32List(nzmax);
    rowIndexesNew.setAll(0, _rowIndexes);
    _rowIndexes = rowIndexesNew;
    var valuesNew = new Float64List(2 * nzmax);
    valuesNew.setAll(0, _values);
    _values = valuesNew;
  }

  String toString() {
    var buf = new StringBuffer();
    buf.write("$rows x $columns sparse matrix, nnz = $cardinality\n");
    for (int i = 0; i < columns; i++) {
      int high = _columnPointers[i + 1];
      for (int j = _columnPointers[i]; j < high; j++) {
        buf.write('(${_rowIndexes[j]},$i)    ${_values[2 * j]}');
        if (_values[2 * j + 1] > 0) {
          buf.write('+${_values[2 * j + 1]}i\n');
        } else if (_values[2 * j + 1] == 0) {
          buf.write('\n');
        } else {
          buf.write('${_values[2 * j + 1]}i\n');
        }
      }
    }
    return buf.toString();
  }

  ComplexVector mult(ComplexVector y, [ComplexVector z = null,
      Complex alpha = null, Complex beta = null, bool transposeA = false]) {
    if (alpha == null) {
      alpha = Complex.ONE;
    }
    if (beta == null) {
      beta = (z == null ? Complex.ONE : Complex.ZERO);
    }
    int rowsA = transposeA ? columns : rows;
    int columnsA = transposeA ? rows : columns;

    bool ignore = (z == null || transposeA);
    if (z == null) {
      z = new DenseComplexVector(rowsA);
    }

    if (y is! DenseComplexVector || z is! DenseComplexVector) {
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
    int zeroZ = zz.index(0);

    DenseComplexVector yy = y as DenseComplexVector;
    Float64List elementsY = yy._elements;
    int strideY = yy.stride;
    int zeroY = yy.index(0);

    Int32List rowIndexesA = _rowIndexes;
    Int32List columnPointersA = _columnPointers;
    Float64List valuesA = _values;

    int zidx = zeroZ;
    if (!transposeA) {
      if ((!ignore) && !(beta.real == 1 && beta.imaginary == 0)) {
        z.apply(cfunc.multiply(beta));
      }

      for (int i = 0; i < columns; i++) {
        int high = columnPointersA[i + 1];
        var yElem = new Complex(
            elementsY[zeroY + strideY * i], elementsY[zeroY + strideY * i + 1]);
        for (int k = columnPointersA[i]; k < high; k++) {
          int j = rowIndexesA[k];
          var valA = new Complex(valuesA[2 * k], valuesA[2 * k + 1]);
          valA = valA * yElem;
          valA = valA * alpha;
          elementsZ[zeroZ + strideZ * j] += valA.real;
          elementsZ[zeroZ + strideZ * j + 1] += valA.imaginary;
        }
      }
    } else {
      for (int i = 0; i < columns; i++) {
        var sum = Complex.ZERO;
        int high = _columnPointers[i + 1];
        for (int k = _columnPointers[i]; k < high; k++) {
          var valA = new Complex(valuesA[2 * k], -valuesA[2 * k + 1]);
          var valY = new Complex(elementsY[zeroY + strideY * _rowIndexes[k]],
              elementsY[zeroY + strideY * _rowIndexes[k] + 1]);
          sum += valA * valY;
        }
        sum = alpha * sum;
        var valZ = new Complex(elementsZ[zidx], elementsZ[zidx + 1]);
        valZ = valZ * beta;
        elementsZ[zidx] = sum.real + valZ.real;
        elementsZ[zidx + 1] = sum.imaginary + valZ.imaginary;
        zidx += strideZ;
      }
    }
    return z;
  }

  ComplexMatrix multiply(ComplexMatrix B, [ComplexMatrix C = null,
      Complex alpha = null, Complex beta = null, bool transposeA = false,
      bool transposeB = false]) {
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
      if (B is SparseCCComplexMatrix) {
        C = new SparseCCComplexMatrix(rowsA, p, rowsA * p);
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

    if (!ignore && !(beta.real == 1 && beta.imaginary == 0)) {
      C.apply(cfunc.multiply(beta));
    }

    if (B is DenseComplexMatrix && C is DenseComplexMatrix) {
      SparseCCComplexMatrix AA;
      if (transposeA) {
        AA = conjugateTranspose();
      } else {
        AA = this;
      }
      DenseComplexMatrix BB;
      if (transposeB) {
        BB = B.conjugateTranspose() as DenseComplexMatrix;
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
      for (int jj = 0; jj < columnsB; jj++) {
        for (int kk = 0; kk < columnsA; kk++) {
          int high = columnPointersA[kk + 1];
          var valB = new Complex(
              elementsB[zeroB + kk * rowStrideB + jj * columnStrideB],
              elementsB[zeroB + kk * rowStrideB + jj * columnStrideB + 1]);
          for (int ii = columnPointersA[kk]; ii < high; ii++) {
            int j = rowIndexesA[ii];
            var valA = new Complex(valuesA[2 * ii], valuesA[2 * ii + 1]);
            valA = valA * valB;
            elementsC[zeroC + j * rowStrideC + jj * columnStrideC] += valA.real;
            elementsC[zeroC + j * rowStrideC + jj * columnStrideC + 1] +=
                valA.imaginary;
          }
        }
      }
      if (!(alpha.real == 1.0 && alpha.imaginary == 0)) {
        C.apply(cfunc.multiply(alpha));
      }
    } else if (B is SparseCCComplexMatrix && C is SparseCCComplexMatrix) {
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
      int nz = 0;
      var columnPointersB = BB._columnPointers;
      var rowIndexesB = BB._rowIndexes;
      var valuesB = BB._values;
      var w = new Int32List(rowsA);
      var x = new Float64List(2 * rowsA);
      var columnPointersC = CC._columnPointers;
      var rowIndexesC = CC._rowIndexes;
      var valuesC = CC._values;
      for (var j = 0; j < columnsB; j++) {
        int nzmaxC = CC._rowIndexes.length;
        if (nz + rowsA > nzmaxC) {
          nzmaxC = 2 * nzmaxC + rowsA;
          var rowIndexesNew = new Int32List(nzmaxC);
          rowIndexesNew.setAll(0, rowIndexesC);
          rowIndexesC = rowIndexesNew;
          var valuesNew = new Float64List(2 * nzmaxC);
          valuesNew.setAll(0, valuesC);
          valuesC = valuesNew;
        }
        columnPointersC[j] = nz;
        for (p = columnPointersB[j]; p < columnPointersB[j + 1]; p++) {
          var elemB = new Complex(valuesB[2 * p], valuesB[2 * p + 1]);
          nz = _toDense(AA, rowIndexesB[p], elemB, w, x, j + 1, CC, nz);
        }
        for (p = columnPointersC[j]; p < nz; p++) {
          valuesC[2 * p] = x[2 * rowIndexesC[p]];
          valuesC[2 * p + 1] = x[2 * rowIndexesC[p] + 1];
        }
      }
      columnPointersC[columnsB] = nz;
      if (!(alpha.real == 1.0 && alpha.imaginary == 0)) {
        CC.apply(cfunc.multiply(alpha));
      }
    } else {
      if (transposeB) {
        B = B.conjugateTranspose();
      }
      // cache views
      var Brows = new List<ComplexVector>(columnsA);
      for (int i = columnsA; --i >= 0;) {
        Brows[i] = B.row(i);
      }
      var Crows = new List<ComplexVector>(rowsA);
      for (int i = rowsA; --i >= 0;) {
        Crows[i] = C.row(i);
      }

      var fn = cfunc.ComplexPlusMultSecond.plusMult(Complex.ZERO);

      for (int i = columns; --i >= 0;) {
        int low = _columnPointers[i];
        for (int k = _columnPointers[i + 1]; --k >= low;) {
          int j = _rowIndexes[k];
          var valA = new Complex(_values[2 * k], _values[2 * k + 1]);
          fn.multiplicator = valA * alpha;
          if (!transposeA) {
            Crows[j].assign(Brows[i], fn);
          } else {
            Crows[i].assign(Brows[j], fn);
          }
        }
      }
    }
    return C;
  }

  ComplexMatrix _getContent() => this;

  void _insert(int row, int column, int index, Complex value) {
    var rowIndexesList = new List.from(_rowIndexes);
    rowIndexesList.length = _columnPointers[columns];
    var valuesList = new List.from(_values);
    valuesList.length = 2 * _columnPointers[columns];
    rowIndexesList.insert(index, row);
    valuesList.insert(2 * index, value.real);
    valuesList.insert(2 * index + 1, value.imaginary);
    for (int i = _columnPointers.length; --i > column;) {
      _columnPointers[i]++;
    }
    _rowIndexes = new Int32List.fromList(rowIndexesList);
    _values = new Float64List.fromList(valuesList);
  }

  void _insertParts(int row, int column, int index, double re, double im) {
    var rowIndexesList = new List.from(_rowIndexes);
    rowIndexesList.length = _columnPointers[columns];
    var valuesList = new List.from(_values);
    valuesList.length = 2 * _columnPointers[columns];
    rowIndexesList.insert(index, row);
    valuesList.insert(2 * index, re);
    valuesList.insert(2 * index + 1, im);
    for (int i = _columnPointers.length; --i > column;) {
      _columnPointers[i]++;
    }
    _rowIndexes = new Int32List.fromList(rowIndexesList);
    _values = new Float64List.fromList(valuesList);
  }

  /*void _remove(int column, int index) {
    var rowIndexesList = new List.from(_rowIndexes);
    var valuesList = new List.from(_values);
    rowIndexesList.remove(index);
    valuesList.remove(2 * index);
    valuesList.remove(2 * index + 1);
    for (int i = _columnPointers.length; --i > column; ) {
      _columnPointers[i]--;
    }
    _rowIndexes = new Int32List.fromList(rowIndexesList);
    _values = new Float64List.fromList(valuesList);
  }*/

  static int _toDense(SparseCCComplexMatrix A, int j, Complex beta, Int32List w,
      Float64List x, int mark, SparseCCComplexMatrix C, int nz) {
    for (var p = A._columnPointers[j]; p < A._columnPointers[j + 1]; p++) {
      var i = A._rowIndexes[p];
      if (w[i] < mark) {
        w[i] = mark;
        C._rowIndexes[nz++] = i;
        if (x != null) {
          var valA = new Complex(A._values[2 * p], A._values[2 * p + 1]);
          valA = beta * valA;
          x[2 * i] = valA.real;
          x[2 * i + 1] = valA.imaginary;
        }
      } else if (x != null) {
        var valA = new Complex(A._values[2 * p], A._values[2 * p + 1]);
        valA = beta * valA;
        x[2 * i] += valA.real;
        x[2 * i + 1] += valA.imaginary;
      }
    }
    return nz;
  }

  Object clone() {
    return new SparseCCComplexMatrix._internal(
        rows, columns, _rowIndexes, _columnPointers, _values);
  }
}
