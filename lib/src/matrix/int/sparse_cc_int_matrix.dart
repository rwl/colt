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
part of cern.colt.matrix.int;

/// Sparse column-compressed 2-d matrix holding [int] elements.
///
/// Internally uses the standard sparse column-compressed format.
///
/// Cells that:
/// - are never set to non-zero values do not use any memory.
/// - switch from zero to non-zero state do use memory.
/// - switch back from non-zero to zero state also do use memory. Their
/// memory is not automatically reclaimed. Reclamation can be triggered
/// via [trimToSize].
///
/// Getting a cell value takes time` O(log nzr)` where `nzr` is the
/// number of non-zeros of the touched row. This is usually quick, because
/// typically there are only few nonzeros per row. So, in practice, get has
/// expected constant time. Setting a cell value takes worst-case time
/// `O(nz)` where `nzr` is the total number of non-zeros in the matrix.
/// This can be extremely slow.
///
/// Fast iteration over non-zeros can be done via [forEachNonZero].
class SparseCCIntMatrix extends WrapperIntMatrix {
  Int32List _columnPointers;

  Int32List _rowIndexes;

  Int32List _values;

  bool _rowIndexesSorted = false;

  /// Constructs a matrix with a given number of rows and columns. All entries
  /// are initially `0`. [nzmax] is the maximum number of nonzero elements.
  factory SparseCCIntMatrix(int rows, int columns, [int nzmax = null]) {
    if (nzmax == null) {
      nzmax = Math.min(10 * rows, MAX_INT);
    }
    var rowIndexes = new Int32List(nzmax);
    var values = new Int32List(nzmax);
    var columnPointers = new Int32List(columns + 1);
    return new SparseCCIntMatrix._internal(
        rows, columns, rowIndexes, columnPointers, values);
  }

  SparseCCIntMatrix._internal(int rows, int columns, Int32List rowIndexes,
      Int32List columnPointers, Int32List values)
      : super._(rows, columns) {
    _rowIndexes = rowIndexes;
    _values = values;
    _columnPointers = columnPointers;
  }

  /// Constructs a matrix with indexes given in the coordinate format and a
  /// single value.
  factory SparseCCIntMatrix.withValue(int rows, int columns,
      Int32List rowIndexes, Int32List columnIndexes, int value,
      {bool removeDuplicates: false, bool sortRowIndexes: false}) {
    if (rowIndexes.length != columnIndexes.length) {
      throw new ArgumentError("rowIndexes.length != columnIndexes.length");
    }
    if (value == 0) {
      throw new ArgumentError("value cannot be 0");
    }

    int nz = Math.max(rowIndexes.length, 1);
    var _rowIndexes = new Int32List(nz);
    var values = new Int32List(nz);
    var columnPointers = new Int32List(columns + 1);
    var columnCounts = new Int32List(columns);
    for (int k = 0; k < nz; k++) {
      columnCounts[columnIndexes[k]]++;
    }
    cumsum(columnPointers, columnCounts);
    for (int k = 0; k < nz; k++) {
      var r = columnCounts[columnIndexes[k]]++;
      _rowIndexes[r] = rowIndexes[k];
      values[r] = value;
    }
    var m = new SparseCCIntMatrix._internal(
        rows, columns, _rowIndexes, columnPointers, values);
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
  factory SparseCCIntMatrix.withValues(int rows, int columns,
      Int32List rowIndexes, Int32List columnIndexes, Int32List values,
      {bool removeDuplicates: false, bool removeZeroes: false,
      bool sortRowIndexes: false}) {
    if (rowIndexes.length != columnIndexes.length) {
      throw new ArgumentError("rowIndexes.length != columnIndexes.length");
    } else if (rowIndexes.length != values.length) {
      throw new ArgumentError("rowIndexes.length != values.length");
    }
    int nz = Math.max(rowIndexes.length, 1);
    var _rowIndexes = new Int32List(nz);
    var _values = new Int32List(nz);
    var columnPointers = new Int32List(columns + 1);
    var columnCounts = new Int32List(columns);
    for (int k = 0; k < nz; k++) {
      columnCounts[columnIndexes[k]]++;
    }
    cumsum(columnPointers, columnCounts);
    for (int k = 0; k < nz; k++) {
      var r = columnCounts[columnIndexes[k]]++;
      _rowIndexes[r] = rowIndexes[k];
      _values[r] = values[k];
    }
    var m = new SparseCCIntMatrix._internal(
        rows, columns, _rowIndexes, columnPointers, _values);
    if (removeDuplicates) {
      m.removeDuplicates();
    }
    if (removeZeroes) {
      m.removeZeroes();
    }
    if (sortRowIndexes) {
      m.sortRowIndexes();
    }
  }

  static SparseCCIntMatrix create(int rows, int columns) {
    return new SparseCCIntMatrix(rows, columns);
  }

  void apply(ifunc.IntFunction fn) {
    if (fn is ifunc.IntMult) {
      // x[i] = mult*x[i]
      int alpha = fn.multiplicator;
      if (alpha == 1) {
        return;
      }
      if (alpha == 0) {
        fill(0);
        return;
      }
      if (alpha != alpha) {
        fill(alpha); // isNaN. This should not happen.
        return;
      }

      for (int j = 0; j < cardinality; j++) {
        _values[j] *= alpha;
      }
    } else {
      forEachNonZero((int i, int j, int value) {
        return fn(value);
      });
    }
  }

  void fill(int value) {
    if (value == 0) {
      _rowIndexes.fillRange(0, _rowIndexes.length, 0);
      _columnPointers.fillRange(0, _columnPointers.length, 0);
      _values.fillRange(0, _values.length, 0);
    } else {
      int nnz = cardinality;
      for (int i = 0; i < nnz; i++) {
        _values[i] = value;
      }
    }
  }

  void copyFrom(IntMatrix source) {
    if (source == this) {
      return; // nothing to do
    }
    checkShape(this, source);

    if (source is SparseCCIntMatrix) {
      _columnPointers.setAll(0, source.columnPointers);
      int nzmax = source.rowIndexes.length;
      if (_rowIndexes.length < nzmax) {
        _rowIndexes = new Int32List(nzmax);
        _values = new Int32List(nzmax);
      }
      _rowIndexes.setAll(0, source.rowIndexes);
      _values.setAll(0, source.values);
      _rowIndexesSorted = source._rowIndexesSorted;
    } else if (source is SparseRCIntMatrix) {
      SparseRCIntMatrix tr = source.transpose();
      _columnPointers = tr.rowPointers;
      _rowIndexes = tr.columnIndexes;
      _values = tr.values;
      _rowIndexesSorted = true;
    } else {
      fill(0);
      source.forEachNonZero((int i, int j, int value) {
        set(i, j, value);
        return value;
      });
    }
  }

  void assign(IntMatrix y, ifunc.IntIntFunction fn) {
    checkShape(this, y);

    if (y is SparseCCIntMatrix && fn == ifunc.plus) {
      // x[i] = x[i] + y[i]
      SparseCCIntMatrix yy = y;
      var nz = 0;
      var anz = _columnPointers[columns];
      var n = yy.columns;
      var bnz = yy._columnPointers[n];
      var w = new Int32List(rows);
      var x = new Int32List(rows);
      var C = new SparseCCIntMatrix(rows, n, anz + bnz);
      var columnPointersC = C._columnPointers;
      var rowIndexesC = C._rowIndexes;
      var valuesC = C._values;
      for (var j = 0; j < n; j++) {
        columnPointersC[j] = nz;
        nz = toDense(_columnPointers, _rowIndexes, _values, j, 1, w, x, j + 1,
            rowIndexesC, nz);
        nz = toDense(yy._columnPointers, yy._rowIndexes, yy._values, j, 1, w, x,
            j + 1, rowIndexesC, nz);
        for (var p = columnPointersC[j]; p < nz; p++) {
          valuesC[p] = x[rowIndexesC[p]];
        }
      }
      columnPointersC[n] = nz;
      _rowIndexes = rowIndexesC;
      _columnPointers = columnPointersC;
      _values = valuesC;
    }

    if (fn is ifunc.IntPlusMultSecond) {
      // x[i] = x[i] + alpha*y[i]
      int alpha = fn.multiplicator;
      if (alpha == 0) {
        return; // nothing to do
      }
      y.forEachNonZero((int i, int j, int value) {
        set(i, j, get(i, j) + alpha * value);
        return value;
      });
      return;
    }

    if (fn is ifunc.IntPlusMultFirst) {
      // x[i] = alpha*x[i] + y[i]
      int alpha = fn.multiplicator;
      if (alpha == 0) {
        copyFrom(y);
        return;
      }
      y.forEachNonZero((int i, int j, int value) {
        set(i, j, alpha * get(i, j) + value);
        return value;
      });
      return;
    }

    if (fn == ifunc.mult) {
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

    if (fn == ifunc.div) {
      // x[i] = x[i] / y[i]
      for (int j = columns; --j >= 0;) {
        int low = _columnPointers[j];
        for (int k = _columnPointers[j + 1]; --k >= low;) {
          int i = _rowIndexes[k];
          _values[k] ~/= y.get(i, j);
          //if (valuesA[k] == 0) _remove(i, j);
        }
      }
      return;
    }
    super.assign(y, fn);
  }

  int get cardinality => _columnPointers[columns];

  void forEachNonZero(ifunc.IntIntIntFunction fn) {
    for (int j = columns; --j >= 0;) {
      int low = _columnPointers[j];
      for (int k = _columnPointers[j + 1]; --k >= low;) {
        int i = rowIndexes[k];
        int value = _values[k];
        int r = fn(i, j, value);
        _values[k] = r;
      }
    }
  }

  /// Column pointers (size [columns]+1).
  Int32List get columnPointers => _columnPointers;

  /// Returns a new matrix that has the same elements as this matrix, but is in
  /// a dense form.
  DenseIntMatrix dense() {
    var dense = new DenseIntMatrix(rows, columns);
    forEachNonZero((int i, int j, int value) {
      dense.set(i, j, get(i, j));
      return value;
    });
    return dense;
  }

  int get(int row, int column) {
    int k = searchRange(_rowIndexes, row, _columnPointers[column],
        _columnPointers[column + 1] - 1);
    if (k >= 0) {
      return _values[k];
    }
    return 0;
  }

  /// Returns a new matrix that has the same elements as this matrix, but is
  /// in a row-compressed form.
  SparseRCIntMatrix rowCompressed() {
    SparseCCIntMatrix tr = transpose();
    var rc = new SparseRCIntMatrix(rows, columns);
    rc._columnIndexes = tr._rowIndexes;
    rc._rowPointers = tr._columnPointers;
    rc._values = tr._values;
    rc._columnIndexesSorted = true;
    return rc;
  }

  /// Row indices (size `nzmax`).
  Int32List get rowIndexes => _rowIndexes;

  /// Returns a new matrix that is the transpose of this matrix.
  SparseCCIntMatrix transpose() {
    var n = columns;
    var tr = new SparseCCIntMatrix(columns, rows, _rowIndexes.length);
    var rowCounts = new Int32List(rows);
    for (var p = 0; p < _columnPointers[n]; p++) {
      rowCounts[_rowIndexes[p]]++;
    }
    cumsum(tr._columnPointers, rowCounts);
    for (var j = 0; j < n; j++) {
      for (var p = _columnPointers[j]; p < _columnPointers[j + 1]; p++) {
        var q = rowCounts[_rowIndexes[p]]++;
        tr._rowIndexes[q] = j;
        tr._values[q] = _values[p];
      }
    }
    return tr;
  }

  /// Numerical values (size `nzmax`).
  Int32List get values => _values;

  bool get rowIndexesSorted => _rowIndexesSorted;

  IntMatrix like2D(int rows, int columns) {
    return new SparseCCIntMatrix(rows, columns);
  }

  IntVector like1D(int size) => new SparseIntVector(size);

  void set(int row, int column, int value) {
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
    SparseCCIntMatrix tr = transpose();
    tr = tr.transpose();
    _columnPointers = tr._columnPointers;
    _rowIndexes = tr._rowIndexes;
    _values = tr._values;
    _rowIndexesSorted = true;
  }

  /// Removes (sums) duplicate entries (if any}
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

  /// Removes zero entries (if any).
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
    var valuesNew = new Int32List(nzmax);
    valuesNew.setAll(0, _values);
    _values = valuesNew;
  }

  String toString() {
    var buf = new StringBuffer();
    buf.write("$rows x $columns sparse matrix, nnz = $cardinality\n");
    for (int i = 0; i < columns; i++) {
      int high = _columnPointers[i + 1];
      for (int j = _columnPointers[i]; j < high; j++) {
        buf.write('(${_rowIndexes[j]},$i)    ${_values[j]}\n');
      }
    }
    return buf.toString();
  }

  IntVector mult(IntVector y, [IntVector z = null, int alpha = 1,
      int beta = null, bool transposeA = false]) {
    if (beta == null) {
      beta = z == null ? 1 : 0;
    }
    int rowsA = transposeA ? columns : rows;
    int columnsA = transposeA ? rows : columns;

    bool ignore = (z == null || transposeA);
    if (z == null) {
      z = new DenseIntVector(rowsA);
    }

    if (y is! DenseIntVector || z is! DenseIntVector) {
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

    DenseIntVector zz = z as DenseIntVector;
    Int32List elementsZ = zz._elements;
    int strideZ = zz.stride;
    int zeroZ = zz.index(0);

    DenseIntVector yy = y as DenseIntVector;
    Int32List elementsY = yy._elements;
    int strideY = yy.stride;
    int zeroY = yy.index(0);

    Int32List rowIndexesA = _rowIndexes;
    Int32List columnPointersA = _columnPointers;
    Int32List valuesA = _values;

    int zidx = zeroZ;
    if (!transposeA) {
      if (!ignore && beta != 1) {
        z.apply(ifunc.multiply(beta));
      }

      for (int i = 0; i < columns; i++) {
        int high = columnPointersA[i + 1];
        int yElem = elementsY[zeroY + strideY * i];
        for (int k = columnPointersA[i]; k < high; k++) {
          int j = rowIndexesA[k];
          elementsZ[zeroZ + strideZ * j] += alpha * valuesA[k] * yElem;
        }
      }
    } else {
      int k = _columnPointers[0];
      for (int i = 0; i < columns; i++) {
        int sum = 0;
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

  IntMatrix multiply(IntMatrix B, [IntMatrix C = null, int alpha = 1,
      int beta = null, bool transposeA = false,
      bool transposeB = false]) {
    if (beta == null) {
      beta = C == null ? 1 : 0;
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
      if (B is SparseCCIntMatrix) {
        C = new SparseCCIntMatrix(rowsA, p, (rowsA * p));
      } else {
        C = new DenseIntMatrix(rowsA, p);
      }
    }

    if (rowsB != columnsA) {
      throw new ArgumentError("Matrix2D inner dimensions must agree:" +
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
      C.apply(ifunc.multiply(beta));
    }

    if (B is DenseIntMatrix && C is DenseIntMatrix) {
      SparseCCIntMatrix AA;
      if (transposeA) {
        AA = transpose();
      } else {
        AA = this;
      }
      DenseIntMatrix BB;
      if (transposeB) {
        BB = B.dice() as DenseIntMatrix;
      } else {
        BB = B;
      }
      DenseIntMatrix CC = C;
      Int32List columnPointersA = AA._columnPointers;
      Int32List rowIndexesA = AA._rowIndexes;
      Int32List valuesA = AA._values;

      int zeroB = BB.index(0, 0);
      int rowStrideB = BB.rowStride;
      int columnStrideB = BB.columnStride;
      Int32List elementsB = BB._elements;

      int zeroC = CC.index(0, 0);
      int rowStrideC = CC.rowStride;
      int columnStrideC = CC.columnStride;
      Int32List elementsC = CC._elements;

      for (int jj = 0; jj < columnsB; jj++) {
        for (int kk = 0; kk < columnsA; kk++) {
          int high = columnPointersA[kk + 1];
          int yElem = elementsB[zeroB + kk * rowStrideB + jj * columnStrideB];
          for (int ii = columnPointersA[kk]; ii < high; ii++) {
            int j = rowIndexesA[ii];
            elementsC[zeroC + j * rowStrideC + jj * columnStrideC] +=
                valuesA[ii] * yElem;
          }
        }
      }
      if (alpha != 1.0) {
        C.apply(ifunc.multiply(alpha));
      }
    } else if (B is SparseCCIntMatrix && C is SparseCCIntMatrix) {
      SparseCCIntMatrix AA;
      if (transposeA) {
        AA = transpose();
      } else {
        AA = this;
      }
      SparseCCIntMatrix BB = B;
      if (transposeB) {
        BB = BB.transpose();
      }
      SparseCCIntMatrix CC = C;
      int nz = 0;
      var columnPointersB = BB._columnPointers;
      var rowIndexesB = BB._rowIndexes;
      var valuesB = BB._values;
      var w = new Int32List(rowsA);
      var x = new Int32List(rowsA);
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
          var valuesNew = new Int32List(nzmaxC);
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
        CC.apply(ifunc.multiply(alpha));
      }
    } else {
      if (transposeB) {
        B = B.dice();
      }
      // cache views
      var Brows = new List<IntVector>(columnsA);
      for (int i = columnsA; --i >= 0;) {
        Brows[i] = B.row(i);
      }
      var Crows = new List<IntVector>(rowsA);
      for (int i = rowsA; --i >= 0;) {
        Crows[i] = C.row(i);
      }

      var fn = new ifunc.IntPlusMultSecond.plusMult(0);

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

  IntMatrix _getContent() => this;

  void _insert(int row, int column, int index, int value) {
    var rowIndexesList = new List.from(_rowIndexes);
    rowIndexesList.length = _columnPointers[columns];
    var valuesList = new List.from(_values);
    valuesList.length = _columnPointers[columns];
    rowIndexesList.insert(index, row);
    valuesList.insert(index, value);
    for (int i = _columnPointers.length; --i > column;) {
      _columnPointers[i]++;
    }
    _rowIndexes = new Int32List.fromList(rowIndexesList);
    _values = new Int32List.fromList(valuesList);
  }

  /*void _remove(int column, int index) {
    var rowIndexesList = new List.from(_rowIndexes);
    var valuesList = new List.from(_values);
    rowIndexesList.remove(index);
    valuesList.remove(index);
    for (int i = _columnPointers.length; --i > column;) {
      _columnPointers[i]--;
    }
    _rowIndexes = new Int32List.fromList(rowIndexesList);
    _values = new Int32List.fromList(valuesList);
  }*/

  Object clone() {
    return new SparseCCIntMatrix._internal(
        rows, columns, _rowIndexes, _columnPointers, _values);
  }
}
