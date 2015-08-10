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

/// Sparse row-compressed 2-d matrix holding [double] elements.
///
/// Internally uses the standard sparse row-compressed format.
///
/// Cells that
/// - are never set to non-zero values do not use any memory.
/// - switch from zero to non-zero state do use memory.
/// - switch back from non-zero to zero state also do use memory. Their memory
/// is not automatically reclaimed.
///
/// Getting a cell value takes time `O(log nzr)` where `nzr` is the
/// number of non-zeros of the touched row. This is usually quick, because
/// typically there are only few nonzeros per row. So, in practice, get has
/// *expected* constant time. Setting a cell value takes worst-case
/// time `O(nz)` where `nzr` is the total number of non-zeros in
/// the matrix. This can be extremely slow, but if you traverse coordinates
/// properly (i.e. upwards), each write is done much quicker.
///
/// Fast iteration over non-zeros can be done via [forEachNonZero].
class SparseRCDoubleMatrix extends WrapperDoubleMatrix {
  Int32List _rowPointers;

  Int32List _columnIndexes;

  Float64List _values;

  bool _columnIndexesSorted = false;

  /// Constructs a matrix with a given number of rows and columns. All entries
  /// are initially `0`. [nzmax] is the maximum number of nonzero elements.
  factory SparseRCDoubleMatrix(int rows, int columns, [int nzmax = null]) {
    if (nzmax == null) {
      nzmax = Math.min(10 * rows, MAX_INT);
    }
    var columnIndexes = new Int32List(nzmax);
    var values = new Float64List(nzmax);
    var rowPointers = new Int32List(rows + 1);
    return new SparseRCDoubleMatrix._internal(
        rows, columns, rowPointers, columnIndexes, values);
  }

  /// Constructs a matrix with indexes given in the coordinate format and
  /// single value.
  factory SparseRCDoubleMatrix.withValue(int rows, int columns,
      Int32List rowIndexes, Int32List columnIndexes, double value,
      {bool removeDuplicates: false, bool sortColumnIndexes: false}) {
    if (rowIndexes.length != columnIndexes.length) {
      throw new ArgumentError("rowIndexes.length != columnIndexes.length");
    }
    if (value == 0) {
      throw new ArgumentError("value cannot be 0");
    }

    int nz = Math.max(rowIndexes.length, 1);
    var _columnIndexes = new Int32List(nz);
    var values = new Float64List(nz);
    var rowPointers = new Int32List(rows + 1);
    var rowCounts = new Int32List(rows);
    for (int k = 0; k < nz; k++) {
      rowCounts[rowIndexes[k]]++;
    }
    cumsum(rowPointers, rowCounts);
    for (int k = 0; k < nz; k++) {
      var r = rowCounts[rowIndexes[k]]++;
      _columnIndexes[r] = columnIndexes[k];
      values[r] = value;
    }
    final m = new SparseRCDoubleMatrix._internal(
        rows, columns, rowPointers, _columnIndexes, values);
    if (removeDuplicates) {
      m.removeDuplicates();
    }
    if (sortColumnIndexes) {
      m.sortColumnIndexes();
    }
    return m;
  }

  /// Constructs a matrix with indexes and values given in the coordinate
  /// format.
  factory SparseRCDoubleMatrix.withValues(int rows, int columns,
      Int32List rowIndexes, Int32List columnIndexes, Float64List values,
      {bool removeDuplicates: true, bool removeZeroes: false,
      bool sortColumnIndexes: false}) {
    if (rowIndexes.length != columnIndexes.length) {
      throw new ArgumentError("rowIndexes.length != columnIndexes.length");
    } else if (rowIndexes.length != values.length) {
      throw new ArgumentError("rowIndexes.length != values.length");
    }
    int nz = Math.max(rowIndexes.length, 1);
    var _columnIndexes = new Int32List(nz);
    var _values = new Float64List(nz);
    var rowPointers = new Int32List(rows + 1);
    var rowCounts = new Int32List(rows);
    for (int k = 0; k < nz; k++) {
      rowCounts[rowIndexes[k]]++;
    }
    cumsum(rowPointers, rowCounts);
    for (int k = 0; k < nz; k++) {
      var r = rowCounts[rowIndexes[k]]++;
      _columnIndexes[r] = columnIndexes[k];
      _values[r] = values[k];
    }
    final m = new SparseRCDoubleMatrix._internal(
        rows, columns, rowPointers, _columnIndexes, _values);
    if (removeZeroes) {
      m.removeZeroes();
    }
    if (removeDuplicates) {
      m.removeDuplicates();
    }
    if (sortColumnIndexes) {
      m.sortColumnIndexes();
    }
    return m;
  }

  SparseRCDoubleMatrix._internal(int rows, int columns, Int32List rowPointers,
      Int32List columnIndexes, Float64List values)
      : super._(rows, columns) {
    if (rowPointers.length != rows + 1) {
      throw new ArgumentError("rowPointers.length != rows + 1");
    }
    _rowPointers = rowPointers;
    _columnIndexes = columnIndexes;
    _values = values;
  }

  static SparseRCDoubleMatrix create(int rows, int columns) {
    return new SparseRCDoubleMatrix(rows, columns);
  }

  void apply(final func.DoubleFunction fn) {
    if (fn is func.DoubleMult) {
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
      forEachNonZero((i, j, double value) => fn(value));
    }
  }

  void fill(double value) {
    if (value == 0.0) {
      _rowPointers.fillRange(0, _rowPointers.length, 0);
      _columnIndexes.fillRange(0, _columnIndexes.length, 0);
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

    if (source is SparseRCDoubleMatrix) {
      _rowPointers.setAll(0, source._rowPointers);
      int nzmax = source._columnIndexes.length;
      if (_columnIndexes.length < nzmax) {
        _columnIndexes = new Int32List(nzmax);
        _values = new Float64List(nzmax);
      }
      _columnIndexes.setAll(0, source._columnIndexes);
      _values.setAll(0, source._values);
      _columnIndexesSorted = source._columnIndexesSorted;
    } else if (source is SparseCCDoubleMatrix) {
      SparseCCDoubleMatrix other = source.transpose();
      _rowPointers = other.columnPointers;
      _columnIndexes = other.rowIndexes;
      _values = other.values;
      _columnIndexesSorted = true;
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
    if (y is SparseRCDoubleMatrix && fn == func.plus) {
      // x[i] = x[i] + y[i]
      SparseRCDoubleMatrix yy = y;

      Int32List rowPointersY = yy._rowPointers;
      Int32List columnIndexesY = yy._columnIndexes;
      Float64List valuesY = yy._values;

      var rowPointersC = new Int32List(rows + 1);
      int cnz = Math.max(_columnIndexes.length,
          Math.min(MAX_INT, _rowPointers[rows] + rowPointersY[rows]));
      var columnIndexesC = new Int32List(cnz);
      var valuesC = new Float64List(cnz);
      if (fn == func.plus) {
        // x[i] = x[i] + y[i]
        int kc = 0;
        rowPointersC[0] = kc;
        int j1, j2;
        for (int i = 0; i < rows; i++) {
          int ka = _rowPointers[i];
          int kb = rowPointersY[i];
          int kamax = _rowPointers[i + 1] - 1;
          int kbmax = rowPointersY[i + 1] - 1;
          while (ka <= kamax || kb <= kbmax) {
            if (ka <= kamax) {
              j1 = _columnIndexes[ka];
            } else {
              j1 = columns + 1;
            }
            if (kb <= kbmax) {
              j2 = columnIndexesY[kb];
            } else {
              j2 = columns + 1;
            }
            if (j1 == j2) {
              valuesC[kc] = _values[ka] + valuesY[kb];
              columnIndexesC[kc] = j1;
              ka++;
              kb++;
              kc++;
            } else if (j1 < j2) {
              columnIndexesC[kc] = j1;
              valuesC[kc] = _values[ka];
              ka++;
              kc++;
            } else if (j1 > j2) {
              columnIndexesC[kc] = j2;
              valuesC[kc] = valuesY[kb];
              kb++;
              kc++;
            }
            if (kc >= cnz) {
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

    if (fn is func.DoublePlusMultSecond) {
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

    if (fn is func.DoublePlusMultFirst) {
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
      for (int i = rows; --i >= 0;) {
        int low = _rowPointers[i];
        for (int k = _rowPointers[i + 1]; --k >= low;) {
          int j = _columnIndexes[k];
          _values[k] *= y.get(i, j);
          //if (_values[k] == 0) _remove(i, j);
        }
      }
      return;
    }

    if (fn == func.div) {
      // x[i] = x[i] / y[i]
      for (int i = rows; --i >= 0;) {
        int low = _rowPointers[i];
        for (int k = _rowPointers[i + 1]; --k >= low;) {
          int j = _columnIndexes[k];
          _values[k] /= y.get(i, j);
          //if (_values[k] == 0) _remove(i, j);
        }
      }
      return;
    }
    super.assign(y, fn);
  }

  int get cardinality => _rowPointers[rows];

  void forEachNonZero(final func.IntIntDoubleFunction fn) {
    for (int i = rows; --i >= 0;) {
      int low = _rowPointers[i];
      for (int k = _rowPointers[i + 1]; --k >= low;) {
        int j = _columnIndexes[k];
        double value = _values[k];
        double r = fn(i, j, value);
        if (r != value) {
          _values[k] = r;
        }
      }
    }
  }

  /// Returns a new matrix that has the same elements as this matrix, but is in
  /// a column-compressed form.
  SparseCCDoubleMatrix columnCompressed() {
    SparseRCDoubleMatrix tr = transpose();
    var cc = new SparseCCDoubleMatrix(rows, columns);
    cc._rowIndexes = tr._columnIndexes;
    cc._columnPointers = tr._rowPointers;
    cc._values = tr._values;
    cc._rowIndexesSorted = true;
    return cc;
  }

  Int32List get columnIndexes => _columnIndexes;

  /// Returns a new matrix that has the same elements as this matrix, but is in
  /// a dense form.
  DenseDoubleMatrix dense() {
    var dense = new DenseDoubleMatrix(rows, columns);
    forEachNonZero((int i, int j, double value) {
      dense.set(i, j, get(i, j));
      return value;
    });
    return dense;
  }

  double get(int row, int column) {
    int k = searchRange(
        _columnIndexes, column, _rowPointers[row], _rowPointers[row + 1] - 1);

    if (k >= 0) {
      return _values[k];
    }
    return 0.0;
  }

  /// Row pointers (size [rows]+1).
  Int32List get rowPointers => _rowPointers;

  /// Returns a new matrix that is the transpose of this matrix.
  SparseRCDoubleMatrix transpose() {
    int nnz = _rowPointers[rows];
    var columnCounts = new Int32List(columns);
    var rowPointersT = new Int32List(columns + 1);
    var columnIndexesT = new Int32List(nnz);
    var valuesT = new Float64List(nnz);

    for (int p = 0; p < nnz; p++) {
      columnCounts[_columnIndexes[p]]++;
    }
    cumsum(rowPointersT, columnCounts);
    for (int j = 0; j < rows; j++) {
      int high = _rowPointers[j + 1];
      for (int p = _rowPointers[j]; p < high; p++) {
        var q = columnCounts[_columnIndexes[p]]++;
        columnIndexesT[q] = j;
        valuesT[q] = _values[p];
      }
    }
    var T = new SparseRCDoubleMatrix(columns, rows);
    T._rowPointers = rowPointersT;
    T._columnIndexes = columnIndexesT;
    T._values = valuesT;
    return T;
  }

  /// Numerical values (size `nzmax`).
  Float64List get values => _values;

  bool get columnIndexesSorted => _columnIndexesSorted;

  DoubleMatrix like2D(int rows, int columns) {
    return new SparseRCDoubleMatrix(rows, columns);
  }

  DoubleVector like1D(int size) => new SparseDoubleVector(size);

  /// Removes (sums) duplicate entries (if any).
  void removeDuplicates() {
    int nz = 0;
    var c = new Int32List(columns);
    for (var i = 0; i < columns; i++) {
      c[i] = -1;
    }

    for (int j = 0; j < rows; j++) {
      var q = nz;
      for (int p = _rowPointers[j]; p < _rowPointers[j + 1]; p++) {
        var i = _columnIndexes[p];
        if (c[i] >= q) {
          _values[c[i]] += _values[p];
        } else {
          c[i] = nz;
          _columnIndexes[nz] = i;
          _values[nz++] = _values[p];
        }
      }
      _rowPointers[j] = q;
    }
    _rowPointers[rows] = nz;
  }

  /// Removes zero entries (if any).
  void removeZeroes() {
    int nz = 0;
    for (int j = 0; j < rows; j++) {
      int p = _rowPointers[j];
      _rowPointers[j] = nz;
      for (; p < _rowPointers[j + 1]; p++) {
        if (_values[p].abs() > EPS) {
          _values[nz] = _values[p];
          _columnIndexes[nz++] = _columnIndexes[p];
        }
      }
    }
    _rowPointers[rows] = nz;
  }

  void set(int row, int column, double value) {
    int k = searchRange(
        _columnIndexes, column, _rowPointers[row], _rowPointers[row + 1] - 1);

    if (k >= 0) {
      //if (value == 0) _remove(row, k); else
      _values[k] = value;
      return;
    }

    if (value != 0) {
      k = -k - 1;
      _insert(row, column, k, value);
    }
  }

  void sortColumnIndexes() {
    SparseRCDoubleMatrix T = transpose();
    T = T.transpose();
    _columnIndexes = T._columnIndexes;
    _rowPointers = T._rowPointers;
    _values = T._values;
    _columnIndexesSorted = true;
  }

  String toString() {
    var buf = new StringBuffer();
    buf.write("$rows x $columns sparse matrix, nnz = $cardinality\n");
    for (int i = 0; i < rows; i++) {
      int high = _rowPointers[i + 1];
      for (int j = _rowPointers[i]; j < high; j++) {
        buf.write('($i,${_columnIndexes[j]})    ${_values[j]}\n');
      }
    }
    return buf.toString();
  }

  void trimToSize() {
    var nzmax = _rowPointers[rows];
    var columnIndexesNew = new Int32List(nzmax);
    columnIndexesNew.setAll(0, _columnIndexes);
    _columnIndexes = columnIndexesNew;
    Float64List valuesNew = new Float64List(nzmax);
    valuesNew.setAll(0, _values);
    _values = valuesNew;
  }

  DoubleVector mult(DoubleVector y, [DoubleVector z = null,
      final double alpha = 1.0, final double beta = 0.0,
      final bool transposeA = false]) {
    int rowsA = transposeA ? columns : rows;
    int columnsA = transposeA ? rows : columns;

    bool ignore = (z == null || !transposeA);
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

    DenseDoubleVector zz = z as DenseDoubleVector;
    Float64List elementsZ = zz._elements;
    int strideZ = zz.stride;
    int zeroZ = z.index(0);

    DenseDoubleVector yy = y as DenseDoubleVector;
    Float64List elementsY = yy._elements;
    int strideY = yy.stride;
    int zeroY = y.index(0);

    if (transposeA) {
      if (!ignore && beta != 1.0) {
        z.apply(func.multiply(beta));
      }
      for (int i = 0; i < rows; i++) {
        int high = _rowPointers[i + 1];
        double yElem = alpha * elementsY[zeroY + strideY * i];
        for (int k = _rowPointers[i]; k < high; k++) {
          int j = _columnIndexes[k];
          elementsZ[zeroZ + strideZ * j] += _values[k] * yElem;
        }
      }
      return z;
    }

    int zidx = zeroZ;
    int k = _rowPointers[0];
    if (beta == 0.0) {
      for (int i = 0; i < rows; i++) {
        double sum = 0.0;
        int high = _rowPointers[i + 1];
        for (; k + 10 < high; k += 10) {
          int ind = k + 9;
          sum += _values[ind] *
                  elementsY[zeroY + strideY * _columnIndexes[ind--]] +
              _values[ind] *
                  elementsY[zeroY + strideY * _columnIndexes[ind--]] +
              _values[ind] *
                  elementsY[zeroY + strideY * _columnIndexes[ind--]] +
              _values[ind] *
                  elementsY[zeroY + strideY * _columnIndexes[ind--]] +
              _values[ind] *
                  elementsY[zeroY + strideY * _columnIndexes[ind--]] +
              _values[ind] *
                  elementsY[zeroY + strideY * _columnIndexes[ind--]] +
              _values[ind] *
                  elementsY[zeroY + strideY * _columnIndexes[ind--]] +
              _values[ind] *
                  elementsY[zeroY + strideY * _columnIndexes[ind--]] +
              _values[ind] *
                  elementsY[zeroY + strideY * _columnIndexes[ind--]] +
              _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]];
        }
        for (; k < high; k++) {
          sum += _values[k] * elementsY[_columnIndexes[k]];
        }
        elementsZ[zidx] = alpha * sum;
        zidx += strideZ;
      }
    } else {
      for (int i = 0; i < rows; i++) {
        double sum = 0.0;
        int high = _rowPointers[i + 1];
        for (; k + 10 < high; k += 10) {
          int ind = k + 9;
          sum += _values[ind] *
                  elementsY[zeroY + strideY * _columnIndexes[ind--]] +
              _values[ind] *
                  elementsY[zeroY + strideY * _columnIndexes[ind--]] +
              _values[ind] *
                  elementsY[zeroY + strideY * _columnIndexes[ind--]] +
              _values[ind] *
                  elementsY[zeroY + strideY * _columnIndexes[ind--]] +
              _values[ind] *
                  elementsY[zeroY + strideY * _columnIndexes[ind--]] +
              _values[ind] *
                  elementsY[zeroY + strideY * _columnIndexes[ind--]] +
              _values[ind] *
                  elementsY[zeroY + strideY * _columnIndexes[ind--]] +
              _values[ind] *
                  elementsY[zeroY + strideY * _columnIndexes[ind--]] +
              _values[ind] *
                  elementsY[zeroY + strideY * _columnIndexes[ind--]] +
              _values[ind] * elementsY[zeroY + strideY * _columnIndexes[ind--]];
        }
        for (; k < high; k++) {
          sum += _values[k] * elementsY[_columnIndexes[k]];
        }
        elementsZ[zidx] = alpha * sum + beta * elementsZ[zidx];
        zidx += strideZ;
      }
    }

    return z;
  }

  DoubleMatrix multiply(DoubleMatrix B, [DoubleMatrix C = null,
      final double alpha = 1.0, double beta = 0.0,
      final bool transposeA = false, bool transposeB = false]) {
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
      if (B is SparseRCDoubleMatrix) {
        C = new SparseRCDoubleMatrix(rowsA, p, rowsA * p);
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
      SparseRCDoubleMatrix AA;
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
      Int32List rowPointersA = AA._rowPointers;
      Int32List columnIndexesA = AA._columnIndexes;
      Float64List valuesA = AA._values;

      for (int ii = 0; ii < rowsA; ii++) {
        int highA = rowPointersA[ii + 1];
        for (int ka = rowPointersA[ii]; ka < highA; ka++) {
          double scal = valuesA[ka] * alpha;
          int jj = columnIndexesA[ka];
          CC.row(ii).assign(BB.row(jj), func.plusMultSecond(scal));
        }
      }
    } else if ((B is SparseRCDoubleMatrix) && (C is SparseRCDoubleMatrix)) {
      SparseRCDoubleMatrix AA,
          BB,
          CC = C;
      if (transposeA) {
        AA = transpose();
      } else {
        AA = this;
      }
      if (transposeB) {
        BB = B.transpose();
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
      int nzmax = valuesC.length;

      var iw = new Int32List(columnsB + 1);
      for (int i = 0; i < iw.length; i++) {
        iw[i] = -1;
      }
      int len = -1;
      for (int ii = 0; ii < rowsA; ii++) {
        int highA = rowPointersA[ii + 1];
        for (int ka = rowPointersA[ii]; ka < highA; ka++) {
          double scal = valuesA[ka] * alpha;
          var jj = columnIndexesA[ka];
          var highB = rowPointersB[jj + 1];
          for (int kb = rowPointersB[jj]; kb < highB; kb++) {
            var jcol = columnIndexesB[kb];
            var jpos = iw[jcol];
            if (jpos == -1) {
              len++;
              if (len >= nzmax) {
                throw new ArgumentError(
                    "The max number of nonzero elements in C is too small.");
              }
              columnIndexesC[len] = jcol;
              iw[jcol] = len;
              valuesC[len] = scal * valuesB[kb];
            } else {
              valuesC[jpos] += scal * valuesB[kb];
            }
          }
        }
        for (int k = rowPointersC[ii]; k < len + 1; k++) {
          iw[columnIndexesC[k]] = -1;
        }
        rowPointersC[ii + 1] = len + 1;

        //int length = rowPointersC[ii + 1] - rowPointersC[ii];
        //var columnIndexesCPart = columnIndexesC.part(rowPointersC[ii], length);
        //Int32List indexes = QuickSort.sortIndex(columnIndexesCPart);
        //Arrays.sort(columnIndexesCElements, rowPointersC[ii], rowPointersC[ii + 1]);
        //DoubleVector valuesCPart = valuesC.part(rowPointersC[ii], length).select(indexes);
        //valuesC.part(rowPointersC[ii], length).assign(valuesCPart);
      }
      //CC.columnIndexes.elements = columnIndexesC.elements;
      //CC.columnIndexes.length = columnIndexesSize;
      //CC.values.elements = valuesC.elements;
      //CC.values.length = columnIndexesSize;
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

      var fn = new func.DoublePlusMultSecond.plusMult(0.0);

      Int32List columnIndexesA = _columnIndexes;
      Float64List valuesA = _values;
      for (int i = rows; --i >= 0;) {
        int low = _rowPointers[i];
        for (int k = _rowPointers[i + 1]; --k >= low;) {
          int j = columnIndexesA[k];
          fn.multiplicator = valuesA[k] * alpha;
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

  DoubleMatrix _getContent() => this;

  void _insert(int row, int column, int index, double value) {
    var columnIndexesList = new List.from(_columnIndexes);
    columnIndexesList.length = _rowPointers[rows];

    var valuesList = new List.from(_values);
    valuesList.length = _rowPointers[rows];

    columnIndexesList.insert(index, column);
    valuesList.insert(index, value);
    for (int i = _rowPointers.length; --i > row;) {
      _rowPointers[i]++;
    }
    _columnIndexes = new Int32List.fromList(columnIndexesList);
    _values = new Float64List.fromList(valuesList);
  }

  /*void _remove(int row, int index) {
    var columnIndexesList = new List.from(_columnIndexes);
    columnIndexesList.length = _rowPointers[rows];

    var valuesList = new List.from(_values);
    valuesList.length = _rowPointers[rows];

    columnIndexesList.removeAt(index);
    valuesList.removeAt(index);
    for (int i = _rowPointers.length; --i > row;) {
      _rowPointers[i]--;
    }
    _columnIndexes = new Int32List.fromList(columnIndexesList);
    _values = new Float64List.fromList(valuesList);
  }*/

  Object clone() {
    return new SparseRCDoubleMatrix._internal(
        rows, columns, _rowPointers, _columnIndexes, _values);
  }
}

/// For internal use.
SparseRCDoubleMatrix makeSparseRCDoubleMatrix(int rows, int columns,
    Int32List rowPointers, Int32List columnIndexes, Float64List values) {
  return new SparseRCDoubleMatrix._internal(
      rows, columns, rowPointers, columnIndexes, values);
}
