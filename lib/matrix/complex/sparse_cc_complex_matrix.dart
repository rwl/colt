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

/// Sparse column-compressed 2-d matrix holding `complex` elements.
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
    final rowIndexes = new Int32List(nzmax);
    final values = new Float64List(2 * nzmax);
    final columnPointers = new Int32List(columns + 1);
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
      bool removeDuplicates) {
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
    var w = new Int32List(columns);
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
    final m = new SparseCCComplexMatrix._internal(
        rows, columns, _rowIndexes, _columnPointers, _values);
    if (removeDuplicates) {
      m.removeDuplicates();
    }
    return m;
  }

  /// Constructs a matrix with indexes and values given in the coordinate
  /// format.
  factory SparseCCComplexMatrix.withValues(int rows, int columns,
      Int32List rowIndexes, Int32List columnIndexes, Float64List values,
      {bool removeDuplicates: false, bool removeZeroes: false}) {
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
    final m = new SparseCCComplexMatrix._internal(
        rows, columns, _rowIndexes, _columnPointers, _values);
    if (removeDuplicates) {
      m.removeDuplicates();
    }
    return m;
  }

  void apply(final cfunc.ComplexComplexFunction fn) {
    if (fn is cfunc.ComplexMult) {
      // x[i] = mult*x[i]
      Float64List alpha = fn.multiplicator;
      if (alpha[0] == 1 && alpha[1] == 0) {
        return;
      }
      if (alpha[0] == 0 && alpha[1] == 0) {
        fill(0.0, 0.0);
        return;
      }
      if (alpha[0] != alpha[0] || alpha[1] != alpha[1]) {
        fill(alpha[0], alpha[1]); // isNaN(). This should not happen.
        return;
      }

      final Float64List valuesE = _values;
      int nz = cardinality;
      var valE = new Float64List(2);
      for (int j = 0; j < nz; j++) {
        valE[0] = valuesE[2 * j];
        valE[1] = valuesE[2 * j + 1];
        valE = Complex.multiply(valE, alpha);
        valuesE[2 * j] = valE[0];
        valuesE[2 * j + 1] = valE[1];
      }
    } else {
      forEachNonZero((int i, int j, Float64List value) {
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

  void copyFrom(AbstractComplexMatrix source) {
    if (source == this) {
      return; // nothing to do
    }
    checkShape(this, source);

    if (source is SparseCCComplexMatrix) {
      SparseCCComplexMatrix other = source;
      _columnPointers.setAll(0, other.columnPointers);
      int nzmax = other.rowIndexes.length;
      if (_rowIndexes.length < nzmax) {
        _rowIndexes = new Int32List(nzmax);
        _values = new Float64List(2 * nzmax);
      }
      _rowIndexes.setAll(0, other.rowIndexes);
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
  }

  void assign(
      final AbstractComplexMatrix y, cfunc.ComplexComplexComplexFunction fn) {
    checkShape(this, y);

    if ((y is SparseCCComplexMatrix) && (fn == cfunc.plus)) {
      // x[i] = x[i] + y[i]
      SparseCCComplexMatrix yy = y;
      int nz = 0;
      var m = rows;
      var anz = _columnPointers[columns];
      var n = yy.columns;
      var Bp = yy._columnPointers;
      var bnz = Bp[n];
      var w = new Int32List(m);
      var x = new Float64List(2 * m);
      var C = new SparseCCComplexMatrix(m, n, anz + bnz);
      var Cp = C._columnPointers;
      var Ci = C._rowIndexes;
      var Cx = C._values;
      var one = new Float64List.fromList([1, 0]);
      for (var j = 0; j < n; j++) {
        Cp[j] = nz;
        nz = _scatter(this, j, one, w, x, j + 1, C, nz);
        nz = _scatter(yy, j, one, w, x, j + 1, C, nz);
        for (var p = Cp[j]; p < nz; p++) {
          Cx[2 * p] = x[2 * Ci[p]];
          Cx[2 * p + 1] = x[2 * Ci[p] + 1];
        }
      }
      Cp[n] = nz;
      _rowIndexes = Ci;
      _columnPointers = Cp;
      _values = Cx;
      return;
    }

    if (fn is cfunc.ComplexPlusMultSecond) {
      // x[i] = x[i] + alpha*y[i]
      Float64List alpha =
          fn.multiplicator;
      if (alpha[0] == 0 && alpha[1] == 0) {
        return; // nothing to do
      }
      y.forEachNonZero((int i, int j, Float64List value) {
        set(i, j, Complex.plus(get(i, j), Complex.multiply(alpha, value)));
        return value;
      });
      return;
    }

    if (fn is cfunc.ComplexPlusMultFirst) {
      // x[i] = alpha*x[i] + y[i]
      Float64List alpha = (fn as cfunc.ComplexPlusMultFirst).multiplicator;
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

    if (fn == cfunc.mult) {
      // x[i] = x[i] * y[i]
      Int32List rowIndexesA = _rowIndexes;
      Int32List columnPointersA = _columnPointers;
      Float64List valuesA = _values;
      var valA = new Float64List(2);
      for (int j = columns; --j >= 0;) {
        int low = columnPointersA[j];
        for (int k = columnPointersA[j + 1]; --k >= low;) {
          int i = rowIndexesA[k];
          valA[0] = valuesA[2 * k];
          valA[1] = valuesA[2 * k + 1];
          valA = Complex.multiply(valA, y.get(i, j));
          valuesA[2 * k] = valA[0];
          valuesA[2 * k + 1] = valA[1];
          //if (valuesA[2 * k] == 0 && valuesA[2 * k + 1] == 0) _remove(i, j);
        }
      }
      return;
    }

    if (fn == cfunc.div) {
      // x[i] = x[i] / y[i]
      Int32List rowIndexesA = _rowIndexes;
      Int32List columnPointersA = _columnPointers;
      Float64List valuesA = _values;

      var valA = new Float64List(2);
      for (int j = columns; --j >= 0;) {
        int low = columnPointersA[j];
        for (int k = columnPointersA[j + 1]; --k >= low;) {
          int i = rowIndexesA[k];
          valA[0] = valuesA[2 * k];
          valA[1] = valuesA[2 * k + 1];
          valA = Complex.div_(valA, y.get(i, j));
          valuesA[2 * k] = valA[0];
          valuesA[2 * k + 1] = valA[1];
          //if (valuesA[2 * k] == 0 && valuesA[2 * k + 1] == 0) _remove(i, j);
        }
      }
      return;
    }
    super.assign(y, fn);
  }

  int get cardinality => _columnPointers[columns];

  void forEachNonZero(final cfunc.IntIntComplexFunction function) {
    Int32List rowIndexesA = _rowIndexes;
    Int32List columnPointersA = _columnPointers;
    Float64List valuesA = _values;
    var valA = new Float64List(2);
    for (int j = columns; --j >= 0;) {
      int low = columnPointersA[j];
      for (int k = columnPointersA[j + 1]; --k >= low;) {
        int i = rowIndexesA[k];
        valA[0] = valuesA[2 * k];
        valA[1] = valuesA[2 * k + 1];
        valA = function(i, j, valA);
        valuesA[2 * k] = valA[0];
        valuesA[2 * k + 1] = valA[1];
      }
    }
  }

  Int32List get columnPointers => _columnPointers;

  /// Returns a new matrix that has the same elements as this matrix, but is in
  /// a dense form.
  ComplexMatrix dense() {
    var dense = new ComplexMatrix(rows, columns);
    forEachNonZero((int i, int j, Float64List value) {
      dense.set(i, j, get(i, j));
      return value;
    });
    return dense;
  }

  Float64List get(int row, int column) {
    //int k = Sorting.binarySearchFromTo(dcs.i, row, dcs.p[column], dcs.p[column + 1] - 1);
    int k = _searchFromTo(_rowIndexes, row, _columnPointers[column],
        _columnPointers[column + 1] - 1);
    Float64List v = new Float64List(2);
    if (k >= 0) {
      v[0] = _values[2 * k];
      v[1] = _values[2 * k + 1];
    }
    return v;
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

  Int32List get rowIndexes => _rowIndexes;

  dynamic get elements {
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

  /// Returns a new matrix that is the transpose of this matrix.
  SparseCCComplexMatrix transpose() {
    DZcs dzcst = cxsparse.cs_transpose(elements(), true);
    SparseCCComplexMatrix tr = new SparseCCComplexMatrix._internal(
        dzcst.m, dzcst.n, dzcst.i, dzcst.p, dzcst.x);
    return tr;
  }

  SparseCCComplexMatrix conjugateTranspose() {
    var m = rows;
    var n = columns;
    var Ap = _columnPointers;
    var Ai = _rowIndexes;
    var Ax = _values;
    var C = new SparseCCComplexMatrix(columns, rows, Ai.length);
    var w = new Int32List(m);
    var Cp = C._columnPointers;
    var Ci = C._rowIndexes;
    var Cx = C._values;
    for (var p = 0; p < Ap[n]; p++) {
      w[Ai[p]]++;
    }
    _cumsum(Cp, w, m);
    var q;
    for (var j = 0; j < n; j++) {
      for (var p = Ap[j]; p < Ap[j + 1]; p++) {
        Ci[q = w[Ai[p]]++] = j;
        Cx[2 * q] = Ax[2 * p];
        Cx[2 * q + 1] = -Ax[2 * p + 1];
      }
    }
    return C;
  }

  Float64List get values => _values;

  AbstractComplexMatrix like2D(int rows, int columns) {
    return new SparseCCComplexMatrix(rows, columns);
  }

  AbstractComplexVector like1D(int size) => new SparseComplexVector(size);

  void set(int row, int column, Float64List value) {
    //int k = Sorting.binarySearchFromTo(dcs.i, row, dcs.p[column], dcs.p[column + 1] - 1);
    int k = _searchFromTo(_rowIndexes, row, _columnPointers[column],
        _columnPointers[column + 1] - 1);

    if (k >= 0) {
      // found
      //if (value[0] == 0 && value[1] == 0) {
      //  _remove(column, k);
      //} else {
        _values[2 * k] = value[0];
        _values[2 * k + 1] = value[1];
      //}
      return;
    }

    if (value[0] != 0 || value[1] != 0) {
      k = -k - 1;
      _insert(row, column, k, value);
    }
  }

  void setParts(int row, int column, double re, double im) {
    //int k = Sorting.binarySearchFromTo(dcs.i, row, dcs.p[column], dcs.p[column + 1] - 1);
    int k = _searchFromTo(_rowIndexes, row, _columnPointers[column],
        _columnPointers[column + 1] - 1);

    if (k >= 0) {
      // found
      //if (re == 0 && im == 0) {
      //  _remove(column, k);
      //} else {
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
    var m = rows;
    var n = columns;
    var Ap = _columnPointers;
    var Ai = _rowIndexes;
    var Ax = _values;
    var w = new Int32List(m);
    for (var i = 0; i < m; i++) {
      w[i] = -1;
    }
    for (var j = 0; j < n; j++) {
      var q = nz;
      for (var p = Ap[j]; p < Ap[j + 1]; p++) {
        var i = Ai[p];
        if (w[i] >= q) {
          Ax[w[i]] += Ax[p];
        } else {
          w[i] = nz;
          Ai[nz] = i;
          Ax[2 * nz] = Ax[2 * p];
          Ax[2 * nz + 1] = Ax[2 * p + 1];
          nz++;
        }
      }
      Ap[j] = q;
    }
    Ap[n] = nz;
  }

  /// Removes zero entries (if any).
  void removeZeroes() {
    int nz = 0;
    var n = columns;
    var Ap = _columnPointers;
    var Ai = _rowIndexes;
    var Ax = _values;
    for (var j = 0; j < n; j++) {
      var p = Ap[j];
      Ap[j] = nz;
      for (; p < Ap[j + 1]; p++) {
        if (Ax[p] != 0) {
          Ai[nz] = Ai[p];
          Ax[2 * nz] = Ax[2 * p];
          Ax[2 * nz + 1] = Ax[2 * p + 1];
          nz++;
        }
      }
    }
    Ap[n] = nz;
  }

  void trimToSize() => _realloc(0);

  String toString() {
    var buf = new StringBuffer();
    buf.write("$rows x $columns sparse matrix, nnz = $cardinality\n");
    for (int i = 0; i < columns; i++) {
      int high = _columnPointers[i + 1];
      for (int j = _columnPointers[i]; j < high; j++) {
        if (_values[2 * j + 1] > 0) {
          buf
            ..write('(')
            ..write(_rowIndexes[j])
            ..write(',')
            ..write(i)
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
            ..write(_rowIndexes[j])
            ..write(',')
            ..write(i)
            ..write(')')
            ..write('\t')
            ..write(_values[2 * j])
            ..write('\n');
        } else {
          buf
            ..write('(')
            ..write(_rowIndexes[j])
            ..write(',')
            ..write(i)
            ..write(')')
            ..write('\t')
            ..write(_values[2 * j]) /*..write('-')*/ ..write(_values[2 * j + 1])
            ..write('i')
            ..write('\n');
        }
      }
    }
    return buf.toString();
  }

  AbstractComplexVector mult(AbstractComplexVector y,
      [AbstractComplexVector z = null, Float64List alpha = null,
      Float64List beta = null, bool transposeA = false]) {
    if (alpha == null) {
      alpha = new Float64List.fromList([1.0, 0.0]);
    }
    if (beta == null) {
      beta = (z == null
          ? new Float64List.fromList([1.0, 0.0])
          : new Float64List.fromList([0.0, 0.0]));
    }
    int rowsA = transposeA ? columns : rows;
    int columnsA = transposeA ? rows : columns;

    bool ignore = (z == null || transposeA);
    if (z == null) {
      z = new ComplexVector(rowsA);
    }

    if (!(y is ComplexVector && z is ComplexVector)) {
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

    ComplexVector zz = z as ComplexVector;
    Float64List elementsZ = zz._elements;
    int strideZ = zz.stride;
    int zeroZ = zz.index(0);

    ComplexVector yy = y as ComplexVector;
    Float64List elementsY = yy._elements;
    int strideY = yy.stride;
    int zeroY = yy.index(0);

    Int32List rowIndexesA = _rowIndexes;
    Int32List columnPointersA = _columnPointers;
    Float64List valuesA = _values;

    int zidx = zeroZ;
    if (!transposeA) {
      if ((!ignore) && !(beta[0] == 1 && beta[1] == 0)) {
        z.apply(cfunc.multiply(beta));
      }

      var yElem = new Float64List(2);
      var valA = new Float64List(2);
      for (int i = 0; i < columns; i++) {
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
    } else {
      var valA = new Float64List(2);
      var valY = new Float64List(2);
      var valZ = new Float64List(2);
      for (int i = 0; i < columns; i++) {
        var sum = new Float64List(2);
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
    }
    return z;
  }

  AbstractComplexMatrix multiply(AbstractComplexMatrix B,
      [AbstractComplexMatrix C = null, Float64List alpha = null,
      Float64List beta = null, bool transposeA = false,
      bool transposeB = false]) {
    if (alpha == null) {
      alpha = new Float64List.fromList([1.0, 0.0]);
    }
    if (beta == null) {
      beta = (C == null
          ? new Float64List.fromList([1.0, 0.0])
          : new Float64List.fromList([0.0, 0.0]));
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
        C = new ComplexMatrix(rowsA, p);
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

    if (!ignore && !(beta[0] == 1 && beta[1] == 0)) {
      C.apply(cfunc.multiply(beta));
    }

    if ((B is ComplexMatrix) && (C is ComplexMatrix)) {
      SparseCCComplexMatrix AA;
      if (transposeA) {
        AA = conjugateTranspose();
      } else {
        AA = this;
      }
      ComplexMatrix BB;
      if (transposeB) {
        BB = B.conjugateTranspose() as ComplexMatrix;
      } else {
        BB = B;
      }
      ComplexMatrix CC = C;
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
      var valA = new Float64List(2);
      var valB = new Float64List(2);
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
            elementsC[zeroC + j * rowStrideC + jj * columnStrideC + 1] +=
                valA[1];
          }
        }
      }
      if (!(alpha[0] == 1.0 && alpha[1] == 0)) {
        C.apply(cfunc.multiply(alpha));
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
      int nz = 0;
      var m = rowsA;
      var n = columnsB;
      var Bp = BB._columnPointers;
      var Bi = BB._rowIndexes;
      var Bx = BB._values;
      var w = new Int32List(m);
      var x = new Float64List(2 * m);
      var Cp = CC._columnPointers;
      var Ci = CC._rowIndexes;
      var Cx = CC._values;
      for (var j = 0; j < n; j++) {
        int nzmaxC = CC._rowIndexes.length;
        if (nz + m > nzmaxC) {
          nzmaxC = 2 * nzmaxC + m;
          var rowIndexesNew = new Int32List(nzmaxC);
          rowIndexesNew.setAll(0, Ci);
          Ci = rowIndexesNew;
          var valuesNew = new Float64List(2 * nzmaxC);
          valuesNew.setAll(0, Cx);
          Cx = valuesNew;
        }
        Cp[j] = nz;
        var elemB = new Float64List(2);
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
      if (!(alpha[0] == 1.0 && alpha[1] == 0)) {
        CC.apply(cfunc.multiply(alpha));
      }
    } else {
      if (transposeB) {
        B = B.conjugateTranspose();
      }
      // cache views
      var Brows = new List<AbstractComplexVector>(columnsA);
      for (int i = columnsA; --i >= 0;) {
        Brows[i] = B.row(i);
      }
      var Crows = new List<AbstractComplexVector>(rowsA);
      for (int i = rowsA; --i >= 0;) {
        Crows[i] = C.row(i);
      }

      var fn = cfunc.ComplexPlusMultSecond.plusMult(new Float64List(2));

      Int32List rowIndexesA = _rowIndexes;
      Int32List columnPointersA = _columnPointers;
      Float64List valuesA = _values;
      var valA = new Float64List(2);
      for (int i = columns; --i >= 0;) {
        int low = columnPointersA[i];
        for (int k = columnPointersA[i + 1]; --k >= low;) {
          int j = rowIndexesA[k];
          valA[0] = valuesA[2 * k];
          valA[1] = valuesA[2 * k + 1];
          fn.multiplicator = Complex.multiply(valA, alpha);
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

  AbstractComplexMatrix _getContent() => this;

  void _insert(int row, int column, int index, Float64List value) {
    var rowIndexesList = new List.from(_rowIndexes);
    rowIndexesList.length = _columnPointers[columns];
    var valuesList = new List.from(_values);
    valuesList.length = 2 * _columnPointers[columns];
    rowIndexesList.insert(index, row);
    valuesList.insert(2 * index, value[0]);
    valuesList.insert(2 * index + 1, value[1]);
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
    return nz2;
  }

  void _realloc(int nzmax) {
    if (nzmax <= 0) {
      nzmax = _columnPointers[columns];
    }
    var rowIndexesNew = new Int32List(nzmax);
    int length = Math.min(nzmax, _rowIndexes.length);
    rowIndexesNew.setAll(0, _rowIndexes);
    _rowIndexes = rowIndexesNew;
    var valuesNew = new Float64List(2 * nzmax);
    length = Math.min(nzmax, _values.length);
    valuesNew.setAll(0, _values);
    _values = valuesNew;
  }

  int _scatter(SparseCCComplexMatrix A, int j, Float64List beta, Int32List w,
      Float64List x, int mark, SparseCCComplexMatrix C, int nz) {
    var Ap = A._columnPointers;
    var Ai = A._rowIndexes;
    var Ax = A._values;
    var Ci = C._rowIndexes;
    var valX = new Float64List(2);
    var valA = new Float64List(2);
    for (var p = Ap[j]; p < Ap[j + 1]; p++) {
      var i = Ai[p];
      if (w[i] < mark) {
        w[i] = mark;
        Ci[nz++] = i;
        if (x != null) {
          valA[0] = Ax[2 * p];
          valA[1] = Ax[2 * p + 1];
          valA = Complex.multiply(beta, valA);
          x[2 * i] = valA[0];
          x[2 * i + 1] = valA[1];
        }
      } else if (x != null) {
        valA[0] = Ax[2 * p];
        valA[1] = Ax[2 * p + 1];
        valA = Complex.multiply(beta, valA);
        x[2 * i] += valA[0];
        x[2 * i + 1] += valA[1];
      }
    }
    return nz;
  }

  Object clone() {
    return new SparseCCComplexMatrix._internal(
        rows, columns, _rowIndexes, _columnPointers, _values);
  }
}
