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
    var columnPointers = new Int32List(columns + 1);
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
      {bool removeDuplicates: false, bool sortRowIndexes: false}) {
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
    var columnCounts = new Int32List(columns);
    for (int k = 0; k < nz; k++) {
      columnCounts[columnIndexes[k]]++;
    }
    cumsum(_columnPointers, columnCounts);
    for (int k = 0; k < nz; k++) {
      var r = columnCounts[columnIndexes[k]]++;
      _rowIndexes[r] = rowIndexes[k];
      _values[r] = value;
    }
    var m = new SparseCCDoubleMatrix._internal(
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
      {bool removeDuplicates: false, bool removeZeroes: false,
      bool sortRowIndexes: false}) {
    if (rowIndexes.length != columnIndexes.length) {
      throw new ArgumentError("rowIndexes.length != columnIndexes.length");
    } else if (rowIndexes.length != values.length) {
      throw new ArgumentError("rowIndexes.length != values.length");
    }
    int nz = Math.max(rowIndexes.length, 1);
    var _rowIndexes = new Int32List(nz);
    var _values = new Float64List(nz);
    var _columnPointers = new Int32List(columns + 1);
    var columnCounts = new Int32List(columns);
    for (int k = 0; k < nz; k++) {
      columnCounts[columnIndexes[k]]++;
    }
    cumsum(_columnPointers, columnCounts);
    for (int k = 0; k < nz; k++) {
      var r = columnCounts[columnIndexes[k]]++;
      _rowIndexes[r] = rowIndexes[k];
      _values[r] = values[k];
    }
    var m = new SparseCCDoubleMatrix._internal(
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

  static SparseCCDoubleMatrix create(int rows, int columns) {
    return new SparseCCDoubleMatrix(rows, columns);
  }

  void apply(func.DoubleFunction fn) {
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
        fill(alpha); // isNaN. This should not happen.
        return;
      }

      int nz = cardinality;
      for (int j = 0; j < nz; j++) {
        _values[j] *= alpha;
      }
    } else {
      forEachNonZero((int i, int j, double value) => fn(value));
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

  void assign(DoubleMatrix y, func.DoubleDoubleFunction fn) {
    checkShape(this, y);

    if ((y is SparseCCDoubleMatrix) && (fn == func.plus)) {
      // x[i] = x[i] + y[i]
      var nz = 0;
      var nnz = _columnPointers[columns];
      var ynz = y._columnPointers[y.columns];
      var w = new Int32List(rows);
      var x = new Float64List(rows);
      var C = new SparseCCDoubleMatrix(rows, y.columns, nnz + ynz);
      var columnPointersC = C._columnPointers;
      var rowIndexesC = C._rowIndexes;
      var valuesC = C._values;
      for (var j = 0; j < y.columns; j++) {
        columnPointersC[j] = nz;
        nz = toDense(_columnPointers, _rowIndexes, _values, j, 1.0, w, x, j + 1,
            rowIndexesC, nz);
        nz = toDense(y._columnPointers, y._rowIndexes, y._values, j, 1.0, w, x,
            j + 1, rowIndexesC, nz);
        for (var p = columnPointersC[j]; p < nz; p++) {
          valuesC[p] = x[rowIndexesC[p]];
        }
      }
      columnPointersC[y.columns] = nz;
      _rowIndexes = rowIndexesC;
      _columnPointers = columnPointersC;
      _values = valuesC;
      return;
    }

    if (fn is DoublePlusMultSecond) {
      // x[i] = x[i] + alpha*y[i]
      double alpha = fn.multiplicator;
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
      double alpha = fn.multiplicator;
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
      for (int j = columns; --j >= 0;) {
        int low = _columnPointers[j];
        for (int k = _columnPointers[j + 1]; --k >= low;) {
          int i = _rowIndexes[k];
          _values[k] *= y.get(i, j);
          //if (valuesA[k] == 0) _remove(i, j);
        }
      }
      return;
    }

    if (fn == func.div) {
      // x[i] = x[i] / y[i]
      for (int j = columns; --j >= 0;) {
        int low = _columnPointers[j];
        for (int k = _columnPointers[j + 1]; --k >= low;) {
          int i = _rowIndexes[k];
          _values[k] /= y.get(i, j);
          //if (valuesA[k] == 0) _remove(i, j);
        }
      }
      return;
    }
    super.assign(y, fn);
  }

  int get cardinality => _columnPointers[columns];

  Object get elements => _values;

  void forEachNonZero(func.IntIntDoubleFunction fn) {
    for (int j = columns; --j >= 0;) {
      int low = _columnPointers[j];
      for (int k = _columnPointers[j + 1]; --k >= low;) {
        int i = _rowIndexes[k];
        double value = _values[k];
        double r = fn(i, j, value);
        _values[k] = r;
      }
    }
  }

  /// Column pointers (size [columns]+1).
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
    int k = searchRange(_rowIndexes, row, _columnPointers[column],
        _columnPointers[column + 1] - 1);
    if (k >= 0) {
      return _values[k];
    }
    return 0.0;
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

  /// Row indices (size `nzmax`).
  Int32List get rowIndexes => _rowIndexes;

  /// Returns a new matrix that is the transpose of this matrix.
  SparseCCDoubleMatrix transpose() {
    var tr = new SparseCCDoubleMatrix(columns, rows, _rowIndexes.length);
    var rowCounts = new Int32List(rows);
    for (var p = 0; p < _columnPointers[columns]; p++) {
      rowCounts[_rowIndexes[p]]++;
    }
    cumsum(tr._columnPointers, rowCounts);
    for (var j = 0; j < columns; j++) {
      for (var p = _columnPointers[j]; p < _columnPointers[j + 1]; p++) {
        // place A(i,j) as entry C(j,i)
        var q = rowCounts[_rowIndexes[p]]++;
        tr._rowIndexes[q] = j;
        tr._values[q] = _values[p];
      }
    }
    return tr;
  }

  /// Numerical values (size `nzmax`).
  Float64List get values => _values;

  bool get rowIndexesSorted => _rowIndexesSorted;

  DoubleMatrix like2D(int rows, int columns) {
    return new SparseCCDoubleMatrix(rows, columns);
  }

  DoubleVector like1D(int size) => new SparseDoubleVector(size);

  void set(int row, int column, double value) {
    int k = searchRange(_rowIndexes, row, _columnPointers[column],
        _columnPointers[column + 1] - 1);

    if (k >= 0) {
      //if (value == 0) _remove(column, k); else
      _values[k] = value;
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
          _values[nz++] = _values[p];
        }
      }
      _columnPointers[j] = q;
    }
    _columnPointers[columns] = nz;
  }

  /// Removes zero entries (if any)
  void removeZeroes() {
    int nz = 0;
    for (var j = 0; j < columns; j++) {
      var p = _columnPointers[j];
      _columnPointers[j] = nz;
      for (; p < _columnPointers[j + 1]; p++) {
        if (_values[p] != 0) {
          _values[nz] = _values[p];
          _rowIndexes[nz++] = _rowIndexes[p];
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
    var valuesNew = new Float64List(nzmax);
    valuesNew.setAll(0, _values);
    _values = valuesNew;
  }

  String toString() {
    var buf = new StringBuffer();
    buf.write("$rows x $columns sparse matrix, nnz = $cardinality\n");
    for (int i = 0; i < columns; i++) {
      int high = _columnPointers[i + 1];
      for (int j = _columnPointers[i]; j < high; j++) {
        buf.write('(${_rowIndexes[j]},$i)\t${_values[j]}\n');
      }
    }
    return buf.toString();
  }

  DoubleVector mult(DoubleVector y, [DoubleVector z = null, double alpha = 1.0,
      double beta = 0.0, bool transposeA = false]) {
    int rowsA = transposeA ? columns : rows;
    int columnsA = transposeA ? rows : columns;

    bool ignore = (z == null || transposeA);
    if (z == null) {
      z = new DenseDoubleVector(rowsA);
    }

    if (y is! DenseDoubleVector || z is! DenseDoubleVector) {
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

  DoubleMatrix multiply(DoubleMatrix B, [DoubleMatrix C = null,
      double alpha = 1.0, double beta = 0.0, bool transposeA = false,
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
        C = new SparseCCDoubleMatrix(rowsA, p, rowsA * p);
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

    if (B is DenseDoubleMatrix && C is DenseDoubleMatrix) {
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
    } else if (B is SparseCCDoubleMatrix && C is SparseCCDoubleMatrix) {
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
      var columnPointersB = BB._columnPointers;
      var rowIndexesB = BB._rowIndexes;
      var valuesB = BB._values;
      var w = new Int32List(rowsA);
      var x = new Float64List(rowsA);
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
          var valuesNew = new Float64List(nzmaxC);
          valuesNew.setAll(0, valuesC);
          valuesC = valuesNew;
        }
        columnPointersC[j] = nz;
        for (p = columnPointersB[j]; p < columnPointersB[j + 1]; p++) {
          nz = toDense(AA._columnPointers, AA._rowIndexes, AA._values,
              rowIndexesB[p], valuesB[p], w, x, j + 1, CC._rowIndexes, nz);
        }
        for (p = columnPointersC[j]; p < nz; p++) {
          valuesC[p] = x[rowIndexesC[p]];
        }
      }
      columnPointersC[columnsB] = nz;
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

      for (int i = columns; --i >= 0;) {
        int low = _columnPointers[i];
        for (int k = _columnPointers[i + 1]; --k >= low;) {
          int j = _rowIndexes[k];
          fn.multiplicator = _values[k] * alpha;
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

  Object clone() {
    return new SparseCCDoubleMatrix._internal(
        rows, columns, _rowIndexes, _columnPointers, _values);
  }
}
