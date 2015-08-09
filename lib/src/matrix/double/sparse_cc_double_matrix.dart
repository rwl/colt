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
part of cern.colt.matrix.double;

/// Sparse column-compressed 2-d matrix holding [double] elements.
///
/// Internally uses the standard sparse column-compressed format.
///
/// Cells that:
/// - are never set to non-zero values do not use any memory.
/// - switch from zero to non-zero state do use memory.
/// - switch back from non-zero to zero state also do use memory.
///
/// Getting a cell value takes time `O(log nzr)` where `nzr` is the
/// number of non-zeros of the touched row. This is usually quick, because
/// typically there are only few nonzeros per row. So, in practice, get has
/// *expected* constant time. Setting a cell value takes worst-case
/// time `O(nz)` where `nzr` is the total number of non-zeros in
/// the matrix.
///
/// Fast iteration over non-zeros can be done via [forEachNonZero].
class SparseCCDoubleMatrix extends WrapperDoubleMatrix {
  Int32List _columnPointers;

  Int32List _rowIndexes;

  Float64List _values;

  bool _rowIndexesSorted = false;

  /// Constructs a matrix with a given number of rows and columns. All entries
  /// are initially `0`. [nzmax] is the maximum number of nonzero elements.
  factory SparseCCDoubleMatrix(int rows, int columns, [int nzmax = null]) {
    if (nzmax == null) {
      nzmax = Math.min(10 * rows, MAX_INT);
    }
    var rowIndexes = new Int32List(nzmax);
    var values = new Float64List(nzmax);
    var columnPointers = new Int32List(rows + 1);
    return new SparseCCDoubleMatrix._internal(
        rows, columns, rowIndexes, columnPointers, values);
  }

  SparseCCDoubleMatrix._internal(int rows, int columns, Int32List rowIndexes,
      Int32List columnPointers, Float64List values)
      : super._(rows, columns) {
    _rowIndexes = rowIndexes;
    _values = values;
    _columnPointers = columnPointers;
  }

  /// Constructs a matrix with indexes given in the coordinate format and a
  /// single value.
  factory SparseCCDoubleMatrix.withValue(int rows, int columns,
      Int32List rowIndexes, Int32List columnIndexes, double value,
      bool removeDuplicates, bool sortRowIndexes) {
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
    var _rowIndexes = new Int32List(nz);
    var _values = new Float64List(nz);
    var _columnPointers = new Int32List(columns + 1);
    var w = new Int32List(columns);
    int r;
    for (int k = 0; k < nz; k++) {
      w[columnIndexes[k]]++;
    }
    _cumsum(_columnPointers, w, columns);
    for (int k = 0; k < nz; k++) {
      _rowIndexes[r = w[columnIndexes[k]]++] = rowIndexes[k];
      _values[r] = value;
    }
    final m = new SparseCCDoubleMatrix._internal(
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
  factory SparseCCDoubleMatrix.withValues(int rows, int columns,
      Int32List rowIndexes, Int32List columnIndexes, Float64List values,
      [bool removeDuplicates = false, bool removeZeroes = false,
      bool sortRowIndexes = false]) {
    if (rowIndexes.length != columnIndexes.length) {
      throw new ArgumentError("rowIndexes.length != columnIndexes.length");
    } else if (rowIndexes.length != values.length) {
      throw new ArgumentError("rowIndexes.length != values.length");
    }
    int nz = Math.max(rowIndexes.length, 1);
    var _rowIndexes = new Int32List(nz);
    var _values = new Float64List(nz);
    var _columnPointers = new Int32List(columns + 1);
    var w = new Int32List(columns);
    int r;
    for (int k = 0; k < nz; k++) {
      w[columnIndexes[k]]++;
    }
    _cumsum(_columnPointers, w, columns);
    for (int k = 0; k < nz; k++) {
      _rowIndexes[r = w[columnIndexes[k]]++] = rowIndexes[k];
      _values[r] = values[k];
    }
    final m = new SparseCCDoubleMatrix._internal(
        rows, columns, _rowIndexes, _columnPointers, _values);
    if (removeDuplicates) {
      m.removeDuplicates();
    }
    if (sortRowIndexes) {
      m.sortRowIndexes();
    }
  }

  static SparseCCDoubleMatrix create(int rows, int columns) {
    return new SparseCCDoubleMatrix(rows, columns);
  }

  void apply(final func.DoubleFunction fn) {
    if (fn is DoubleMult) {
      // x[i] = mult*x[i]
      double alpha = fn.multiplicator;
      if (alpha == 1) {
        return;
      }
      if (alpha == 0) {
        fill(0.0);
        return;
      }
      if (alpha != alpha) {
        fill(
            alpha); // the funny definition of isNaN(). This should better not happen.
        return;
      }

      final Float64List valuesE = _values;
      int nz = cardinality;
      for (int j = 0; j < nz; j++) {
        valuesE[j] *= alpha;
      }
    } else {
      forEachNonZero((int i, int j, double value) {
        return fn(value);
      });
    }
  }

  void fill(double value) {
    if (value == 0) {
      _rowIndexes.fillRange(0, _rowIndexes.length, 0);
      _columnPointers.fillRange(0, _columnPointers.length, 0);
      _values.fillRange(0, _values.length, 0.0);
    } else {
      int nnz = cardinality;
      for (int i = 0; i < nnz; i++) {
        _values[i] = value;
      }
    }
  }

  void copyFrom(DoubleMatrix source) {
    if (source == this) {
      return; // nothing to do
    }
    checkShape(this, source);

    if (source is SparseCCDoubleMatrix) {
      SparseCCDoubleMatrix other = source;
      _columnPointers.setAll(0, other.columnPointers);
      int nzmax = other.rowIndexes.length;
      if (_rowIndexes.length < nzmax) {
        _rowIndexes = new Int32List(nzmax);
        _values = new Float64List(nzmax);
      }
      _rowIndexes.setAll(0, other.rowIndexes);
      _values.setAll(0, other.values);
      _rowIndexesSorted = other._rowIndexesSorted;
    } else if (source is SparseRCDoubleMatrix) {
      SparseRCDoubleMatrix other = source.transpose();
      _columnPointers = other.rowPointers;
      _rowIndexes = other.columnIndexes;
      _values = other.values;
      _rowIndexesSorted = true;
    } else {
      fill(0.0);
      source.forEachNonZero((int i, int j, double value) {
        set(i, j, value);
        return value;
      });
    }
  }

  void assign(final DoubleMatrix y, func.DoubleDoubleFunction fn) {
    checkShape(this, y);

    if ((y is SparseCCDoubleMatrix) && (fn == func.plus)) {
      // x[i] = x[i] + y[i]
      SparseCCDoubleMatrix yy = y;
      var nz = 0;
      var m = rows;
      var anz = _columnPointers[columns];
      var n = yy.columns;
      var Bp = yy._columnPointers;
      var bnz = Bp[n];
      var w = new Int32List(m);
      var x = new Float64List(m);
      var C = new SparseCCDoubleMatrix(m, n, anz + bnz);
      var Cp = C._columnPointers;
      var Ci = C._rowIndexes;
      var Cx = C._values;
      for (var j = 0; j < n; j++) {
        // column j of C starts here
        Cp[j] = nz;
        // alpha*A(:,j)
        nz = _scatter(this, j, 1.0, w, x, j + 1, C, nz);
        // beta*B(:,j)
        nz = _scatter(yy, j, 1.0, w, x, j + 1, C, nz);
        for (var p = Cp[j]; p < nz; p++) {
          Cx[p] = x[Ci[p]];
        }
      }
      // finalize the last column of C
      Cp[n] = nz;
      _rowIndexes = Ci;
      _columnPointers = Cp;
      _values = Cx;
      return;
    }

    if (fn is DoublePlusMultSecond) {
      // x[i] = x[i] + alpha*y[i]
      final double alpha = fn.multiplicator;
      if (alpha == 0) {
        return; // nothing to do
      }
      y.forEachNonZero((int i, int j, double value) {
        set(i, j, get(i, j) + alpha * value);
        return value;
      });
      return;
    }

    if (fn is DoublePlusMultFirst) {
      // x[i] = alpha*x[i] + y[i]
      final double alpha = fn.multiplicator;
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

    if (fn == func.mult) {
      // x[i] = x[i] * y[i]
      final Int32List rowIndexesA = _rowIndexes;
      final Int32List columnPointersA = _columnPointers;
      final Float64List valuesA = _values;
      for (int j = columns; --j >= 0;) {
        int low = columnPointersA[j];
        for (int k = columnPointersA[j + 1]; --k >= low;) {
          int i = rowIndexesA[k];
          valuesA[k] *= y.get(i, j);
          //if (valuesA[k] == 0) _remove(i, j);
        }
      }
      return;
    }

    if (fn == func.div) {
      // x[i] = x[i] / y[i]
      final Int32List rowIndexesA = _rowIndexes;
      final Int32List columnPointersA = _columnPointers;
      final Float64List valuesA = _values;

      for (int j = columns; --j >= 0;) {
        int low = columnPointersA[j];
        for (int k = columnPointersA[j + 1]; --k >= low;) {
          int i = rowIndexesA[k];
          valuesA[k] /= y.get(i, j);
          //if (valuesA[k] == 0) _remove(i, j);
        }
      }
      return;
    }
    super.assign(y, fn);
  }

  int get cardinality => _columnPointers[columns];

  Object get elements => _values;

  void forEachNonZero(final func.IntIntDoubleFunction function) {
    final Int32List rowIndexesA = _rowIndexes;
    final Int32List columnPointersA = _columnPointers;
    final Float64List valuesA = _values;

    for (int j = columns; --j >= 0;) {
      int low = columnPointersA[j];
      for (int k = columnPointersA[j + 1]; --k >= low;) {
        int i = rowIndexesA[k];
        double value = valuesA[k];
        double r = function(i, j, value);
        valuesA[k] = r;
      }
    }
  }

  Int32List get columnPointers => _columnPointers;

  /// Returns a new matrix that has the same elements as this matrix, but is
  /// in a dense form.
  DenseDoubleMatrix dense() {
    var dense = new DenseDoubleMatrix(rows, columns);
    forEachNonZero((int i, int j, double value) {
      dense.set(i, j, get(i, j));
      return value;
    });
    return dense;
  }

  double get(int row, int column) {
    //int k = Sorting.binarySearchFromTo(dcs.i, row, dcs.p[column], dcs.p[column + 1] - 1);
    int k = _searchFromTo(_rowIndexes, row, _columnPointers[column],
        _columnPointers[column + 1] - 1);
    double v = 0.0;
    if (k >= 0) {
      v = _values[k];
    }
    return v;
  }

  /// Returns a new matrix that has the same elements as this matrix, but is
  /// in a row-compressed form.
  SparseRCDoubleMatrix rowCompressed() {
    SparseCCDoubleMatrix tr = transpose();
    SparseRCDoubleMatrix rc = new SparseRCDoubleMatrix(rows, columns);
    rc._columnIndexes = tr._rowIndexes;
    rc._rowPointers = tr._columnPointers;
    rc._values = tr._values;
    rc._columnIndexesSorted = true;
    return rc;
  }

  Int32List get rowIndexes => _rowIndexes;

  /// Returns a new matrix that is the transpose of this matrix.
  SparseCCDoubleMatrix transpose() {
    var q;
    var m = rows;
    var n = columns;
    var Ap = _columnPointers;
    var Ai = _rowIndexes;
    var Ax = _values;
    var C = new SparseCCDoubleMatrix(columns, rows, Ai.length);
    var w = new Int32List(m);
    var Cp = C._columnPointers;
    var Ci = C._rowIndexes;
    var Cx = C._values;
    for (var p = 0; p < Ap[n]; p++) {
      w[Ai[p]]++;
    }
    _cumsum(Cp, w, m);
    for (var j = 0; j < n; j++) {
      for (var p = Ap[j]; p < Ap[j + 1]; p++) {
        // place A(i,j) as entry C(j,i)
        Ci[q = w[Ai[p]]++] = j;
        Cx[q] = Ax[p];
      }
    }
    return C;
  }

  Float64List get values => _values;

  bool get rowIndexesSorted => _rowIndexesSorted;

  DoubleMatrix like2D(int rows, int columns) {
    return new SparseCCDoubleMatrix(rows, columns);
  }

  DoubleVector like1D(int size) => new SparseDoubleVector(size);

  void set(int row, int column, double value) {
    //int k = Sorting.binarySearchFromTo(dcs.i, row, dcs.p[column], dcs.p[column + 1] - 1);
    int k = _searchFromTo(_rowIndexes, row, _columnPointers[column],
        _columnPointers[column + 1] - 1);

    if (k >= 0) {
      // found
      /*if (value == 0) {
        _remove(column, k);
      } else {*/
        _values[k] = value;
      //}
      return;
    }

    if (value != 0) {
      k = -k - 1;
      _insert(row, column, k, value);
    }
  }

  void sortRowIndexes() {
    SparseCCDoubleMatrix tr = transpose();
    tr = tr.transpose();
    _columnPointers = tr._columnPointers;
    _rowIndexes = tr._rowIndexes;
    _values = tr._values;
    _rowIndexesSorted = true;
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
          Ax[nz++] = Ax[p];
        }
      }
      Ap[j] = q;
    }
    Ap[n] = nz;
  }

  /// Removes zero entries (if any)
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
          Ax[nz] = Ax[p];
          Ai[nz++] = Ai[p];
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
        buf
          ..write('(')
          ..write(_rowIndexes[j])
          ..write(',')
          ..write(i)
          ..write(')')
          ..write('\t')
          ..write(_values[j])
          ..write('\n');
      }
    }
    return buf.toString();
  }

  DoubleVector mult(DoubleVector y,
      [DoubleVector z = null, final double alpha = 1.0,
      final double beta = 0.0, final bool transposeA = false]) {
    int rowsA = transposeA ? columns : rows;
    int columnsA = transposeA ? rows : columns;

    bool ignore = (z == null || transposeA);
    if (z == null) {
      z = new DenseDoubleVector(rowsA);
    }

    if (!(y is DenseDoubleVector && z is DenseDoubleVector)) {
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

    var zz = z as DenseDoubleVector;
    Float64List elementsZ = zz._elements;
    int strideZ = zz.stride;
    int zeroZ = zz.index(0);

    DenseDoubleVector yy = y as DenseDoubleVector;
    Float64List elementsY = yy._elements;
    int strideY = yy.stride;
    int zeroY = yy.index(0);

    Int32List rowIndexesA = _rowIndexes;
    Int32List columnPointersA = _columnPointers;
    Float64List valuesA = _values;

    int zidx = zeroZ;
    if (!transposeA) {
      if ((!ignore) && (beta / alpha != 1.0)) {
        z.apply(func.multiply(beta / alpha));
      }

      for (int i = 0; i < columns; i++) {
        int high = columnPointersA[i + 1];
        double yElem = elementsY[zeroY + strideY * i];
        for (int k = columnPointersA[i]; k < high; k++) {
          int j = rowIndexesA[k];
          elementsZ[zeroZ + strideZ * j] += valuesA[k] * yElem;
        }
      }

      if (alpha != 1.0) {
        z.apply(func.multiply(alpha));
      }
    } else {
      int k = _columnPointers[0];
      for (int i = 0; i < columns; i++) {
        double sum = 0.0;
        int high = _columnPointers[i + 1];
        for (; k + 10 < high; k += 10) {
          int ind = k + 9;
          sum += valuesA[ind] *
                  elementsY[zeroY + strideY * _rowIndexes[ind--]] +
              valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] +
              valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] +
              valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] +
              valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] +
              valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] +
              valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] +
              valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] +
              valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]] +
              valuesA[ind] * elementsY[zeroY + strideY * _rowIndexes[ind--]];
        }
        for (; k < high; k++) {
          sum += valuesA[k] * elementsY[_rowIndexes[k]];
        }
        elementsZ[zidx] = alpha * sum + beta * elementsZ[zidx];
        zidx += strideZ;
      }
    }
    return z;
  }

  DoubleMatrix multiply(DoubleMatrix B,
      [DoubleMatrix C = null, final double alpha = 1.0,
      double beta = 0.0, final bool transposeA = false,
      bool transposeB = false]) {
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
      if (B is SparseCCDoubleMatrix) {
        C = new SparseCCDoubleMatrix(rowsA, p, (rowsA * p));
      } else {
        C = new DenseDoubleMatrix(rowsA, p);
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

    if (!ignore && beta != 1.0) {
      C.apply(func.multiply(beta));
    }

    if ((B is DenseDoubleMatrix) && (C is DenseDoubleMatrix)) {
      SparseCCDoubleMatrix AA;
      if (transposeA) {
        AA = transpose();
      } else {
        AA = this;
      }
      DenseDoubleMatrix BB;
      if (transposeB) {
        BB = B.dice() as DenseDoubleMatrix;
      } else {
        BB = B;
      }
      DenseDoubleMatrix CC = C;
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
          double yElem =
              elementsB[zeroB + kk * rowStrideB + jj * columnStrideB];
          for (int ii = columnPointersA[kk]; ii < high; ii++) {
            int j = rowIndexesA[ii];
            elementsC[zeroC + j * rowStrideC + jj * columnStrideC] +=
                valuesA[ii] * yElem;
          }
        }
      }
      if (alpha != 1.0) {
        C.apply(func.multiply(alpha));
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
      int nz = 0;
      var m = rowsA;
      var n = columnsB;
      var Bp = BB._columnPointers;
      var Bi = BB._rowIndexes;
      var Bx = BB._values;
      var w = new Int32List(m);
      var x = new Float64List(m);
      var Cp = CC._columnPointers;
      var Ci = CC._rowIndexes;
      var Cx = CC._values;
      for (var j = 0; j < n; j++) {
        int nzmaxC = CC._rowIndexes.length;
        if (nz + m > nzmaxC) {
          nzmaxC = 2 * nzmaxC + m;
          Int32List rowIndexesNew = new Int32List(nzmaxC);
          rowIndexesNew.setAll(0, Ci);
          Ci = rowIndexesNew;
          var valuesNew = new Float64List(nzmaxC);
          valuesNew.setAll(0, Cx);
          Cx = valuesNew;
        }
        Cp[j] = nz;
        for (p = Bp[j]; p < Bp[j + 1]; p++) {
          nz = _scatter(AA, Bi[p], Bx[p], w, x, j + 1, CC, nz);
        }
        for (p = Cp[j]; p < nz; p++) Cx[p] = x[Ci[p]];
      }
      Cp[n] = nz;
      if (alpha != 1.0) {
        CC.apply(func.multiply(alpha));
      }
    } else {
      if (transposeB) {
        B = B.dice();
      }
      // cache views
      var Brows = new List<DoubleVector>(columnsA);
      for (int i = columnsA; --i >= 0;) {
        Brows[i] = B.row(i);
      }
      var Crows = new List<DoubleVector>(rowsA);
      for (int i = rowsA; --i >= 0;) {
        Crows[i] = C.row(i);
      }

      var fn = new DoublePlusMultSecond.plusMult(0.0);

      final Int32List rowIndexesA = _rowIndexes;
      final Int32List columnPointersA = _columnPointers;
      final Float64List valuesA = _values;
      for (int i = columns; --i >= 0;) {
        int low = columnPointersA[i];
        for (int k = columnPointersA[i + 1]; --k >= low;) {
          int j = rowIndexesA[k];
          fn.multiplicator = valuesA[k] * alpha;
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

  DoubleMatrix _getContent() => this;

  void _insert(int row, int column, int index, double value) {
    List rowIndexes = new List<int>.from(_rowIndexes);
    rowIndexes.length = _columnPointers[columns];
    List values = new List.from(_values);
    values.length = _columnPointers[columns];
    rowIndexes.insert(index, row);
    values.insert(index, value);
    for (int i = _columnPointers.length; --i > column;) {
      _columnPointers[i]++;
    }
    _rowIndexes = new Int32List.fromList(rowIndexes);
    _values = new Float64List.fromList(values);
  }

  /*void _remove(int column, int index) {
    List rowIndexes = new List.from(_rowIndexes);
    List values = new List.from(_values);
    rowIndexes.removeAt(index);
    values.removeAt(index);
    for (int i = _columnPointers.length; --i > column;) {
      _columnPointers[i]--;
    }
    _rowIndexes = new Int32List.fromList(rowIndexes);
    _values = new Float64List.fromList(values);
    _dcs.nzmax = rowIndexes.length;
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
    return (nz2);
  }

  void _realloc(int nzmax) {
    if (nzmax <= 0) {
      nzmax = _columnPointers[columns];
    }
    var rowIndexesNew = new Int32List(nzmax);
    rowIndexesNew.setAll(0, _rowIndexes);
    _rowIndexes = rowIndexesNew;
    var valuesNew = new Float64List(nzmax);
    valuesNew.setAll(0, _values);
    _values = valuesNew;
  }

  int _scatter(SparseCCDoubleMatrix A, int j, double beta, Int32List w, Float64List x,
               int mark, SparseCCDoubleMatrix C, int nz) {
    var Ap = A._columnPointers;
    var Ai = A._rowIndexes;
    var Ax = A._values;
    var Ci = C._rowIndexes;
    for (var p = Ap[j]; p < Ap[j + 1]; p++) {
      var i = Ai[p];
      if (w[i] < mark) {
        w[i] = mark;
        Ci[nz++] = i;
        if (x != null) {
          x[i] = beta * Ax[p];
        }
      } else if (x != null) {
        x[i] += beta * Ax[p];
      }
    }
    return nz;
  }

  Object clone() {
    return new SparseCCDoubleMatrix._internal(
        rows, columns, _rowIndexes, _columnPointers, _values);
  }
}
