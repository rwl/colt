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

/// Sparse row-compressed 2-d matrix holding [int] elements.
///
/// Internally uses the standard sparse row-compressed format.
///
/// Cells that:
/// - are never set to non-zero values do not use any memory.
/// -switch from zero to non-zero state do use memory.
/// - switch back from non-zero to zero state also do use memory. Their memory
///  is not automatically reclaimed. Reclamation can be triggered via
///  [trimToSize].
///
///  Fast iteration over non-zeros can be done via [forEachNonZero].
class SparseRCIntMatrix extends WrapperIntMatrix {
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

  Int32List _values;

  bool _columnIndexesSorted = false;

  /// Constructs a matrix with a given number of rows and columns. All entries
  /// are initially `0`.
  factory SparseRCIntMatrix(int rows, int columns, [int nzmax = null]) {
    if (nzmax == null) {
      nzmax = Math.min(10 * rows, MAX_INT);
    }
    final columnIndexes = new Int32List(nzmax);
    final values = new Int32List(nzmax);
    final rowPointers = new Int32List(rows + 1);
    return new SparseRCIntMatrix._internal(
        rows, columns, rowPointers, columnIndexes, values);
  }

  factory SparseRCIntMatrix.withValue(int rows, int columns,
      Int32List rowIndexes, Int32List columnIndexes, int value,
      bool removeDuplicates, bool sortColumnIndexes) {
    if (rowIndexes.length != columnIndexes.length) {
      throw new ArgumentError("rowIndexes.length != columnIndexes.length");
    }
    if (value == 0) {
      throw new ArgumentError("value cannot be 0");
    }

    int nz = Math.max(rowIndexes.length, 1);
    final _columnIndexes = new Int32List(nz);
    final _values = new Int32List(nz);
    final _rowPointers = new Int32List(rows + 1);
    Int32List w = new Int32List(rows);
    int r;
    for (int k = 0; k < nz; k++) {
      w[rowIndexes[k]]++;
    }
    _cumsum(_rowPointers, w, rows);
    for (int k = 0; k < nz; k++) {
      _columnIndexes[r = w[rowIndexes[k]]++] = columnIndexes[k];
      _values[r] = value;
    }
    final m = new SparseRCIntMatrix._internal(
        rows, columns, _rowPointers, _columnIndexes, _values);
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
  factory SparseRCIntMatrix.withValues(int rows, int columns,
      Int32List rowIndexes, Int32List columnIndexes, Int32List values,
      bool removeDuplicates, bool removeZeroes, bool sortColumnIndexes) {
    if (rowIndexes.length != columnIndexes.length) {
      throw new ArgumentError("rowIndexes.length != columnIndexes.length");
    } else if (rowIndexes.length != values.length) {
      throw new ArgumentError("rowIndexes.length != values.length");
    }
    int nz = Math.max(rowIndexes.length, 1);
    final _columnIndexes = new Int32List(nz);
    final _values = new Int32List(nz);
    final _rowPointers = new Int32List(rows + 1);
    Int32List w = new Int32List(rows);
    int r;
    for (int k = 0; k < nz; k++) {
      w[rowIndexes[k]]++;
    }
    _cumsum(_rowPointers, w, rows);
    for (int k = 0; k < nz; k++) {
      _columnIndexes[r = w[rowIndexes[k]]++] = columnIndexes[k];
      _values[r] = values[k];
    }
    final m = new SparseRCIntMatrix._internal(
        rows, columns, _rowPointers, _columnIndexes, _values);
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

  SparseRCIntMatrix._internal(int rows, int columns, Int32List rowPointers,
      Int32List columnIndexes, Int32List values)
      : super._(rows, columns) {
    if (rowPointers.length != rows + 1) {
      throw new ArgumentError("rowPointers.length != rows + 1");
    }
    _rowPointers = rowPointers;
    _columnIndexes = columnIndexes;
    _values = values;
  }

  void apply(final ifunc.IntFunction fn) {
    if (fn is ifunc.IntMult) {
      // x[i] = mult*x[i]
      final int alpha = fn.multiplicator;
      if (alpha == 1) {
        return;
      }
      if (alpha == 0) {
        fill(0);
        return;
      }
      if (alpha != alpha) {
        fill(alpha); // The definition of isNaN(). This should not happen.
        return;
      }

      int nz = cardinality;
      for (int j = 0; j < nz; j++) {
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
      _rowPointers.fillRange(0, _rowPointers.length, 0);
      _columnIndexes.fillRange(0, _columnIndexes.length, 0);
      _values.fillRange(0, _values.length, 0);
    } else {
      int nnz = cardinality;
      for (int i = 0; i < nnz; i++) {
        _values[i] = value;
      }
    }
    return;
  }

  void copyFrom(AbstractIntMatrix source) {
    if (source == this) {
      return; // nothing to do
    }
    checkShape(this, source);

    if (source is SparseRCIntMatrix) {
      SparseRCIntMatrix other = source;
      _rowPointers.setAll(0, other._rowPointers);
      int nzmax = other._columnIndexes.length;
      if (_columnIndexes.length < nzmax) {
        _columnIndexes = new Int32List(nzmax);
        _values = new Int32List(nzmax);
      }
      _columnIndexes.setAll(0, other._columnIndexes);
      _values.setAll(0, other._values);
      _columnIndexesSorted = other._columnIndexesSorted;
    } else if (source is SparseCCIntMatrix) {
      SparseCCIntMatrix other = source.transpose();
      _rowPointers = other.columnPointers;
      _columnIndexes = other.rowIndexes;
      _values = other.values;
      _columnIndexesSorted = true;
    } else {
      fill(0);
      source.forEachNonZero((int i, int j, int value) {
        set(i, j, value);
        return value;
      });
    }
  }

  void assign(final AbstractIntMatrix y, ifunc.IntIntFunction fn) {
    checkShape(this, y);
    if ((y is SparseRCIntMatrix) && (fn == ifunc.plus)) {
      // x[i] = x[i] + y[i]
      SparseRCIntMatrix yy = y;

      Int32List rowPointersY = yy._rowPointers;
      Int32List columnIndexesY = yy._columnIndexes;
      Int32List valuesY = yy._values;

      var rowPointersC = new Int32List(rows + 1);
      int cnz = Math.max(_columnIndexes.length,
          Math.min(MAX_INT, _rowPointers[rows] + rowPointersY[rows]));
      var columnIndexesC = new Int32List(cnz);
      var valuesC = new Int32List(cnz);
      int nrow = rows;
      int ncol = columns;
      int nzmax = valuesC.length;
      if (fn == ifunc.plus) {
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

    if (fn is ifunc.IntPlusMultSecond) {
      // x[i] = x[i] + alpha*y[i]
      final int alpha = fn.multiplicator;
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
      final int alpha = fn.multiplicator;
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
      for (int i = rows; --i >= 0;) {
        int low = _rowPointers[i];
        for (int k = _rowPointers[i + 1]; --k >= low;) {
          int j = _columnIndexes[k];
          _values[k] *= y.get(i, j);
          if (_values[k] == 0) _remove(i, j);
        }
      }
      return;
    }

    if (fn == ifunc.div) {
      // x[i] = x[i] / y[i]

      for (int i = rows; --i >= 0;) {
        int low = _rowPointers[i];
        for (int k = _rowPointers[i + 1]; --k >= low;) {
          int j = _columnIndexes[k];
          _values[k] ~/= y.get(i, j);
          if (_values[k] == 0) {
            _remove(i, j);
          }
        }
      }
      return;
    }
    super.assign(y, fn);
  }

  int get cardinality => _rowPointers[rows];

  void forEachNonZero(final ifunc.IntIntIntFunction function) {
    for (int i = rows; --i >= 0;) {
      int low = _rowPointers[i];
      for (int k = _rowPointers[i + 1]; --k >= low;) {
        int j = _columnIndexes[k];
        int value = _values[k];
        int r = function(i, j, value);
        if (r != value) _values[k] = r;
      }
    }
  }

  /// Returns a new matrix that has the same elements as this matrix, but is
  /// in a column-compressed form.
  SparseCCIntMatrix columnCompressed() {
    SparseRCIntMatrix tr = transpose();
    var cc = new SparseCCIntMatrix(rows, columns);
    cc._rowIndexes = tr._columnIndexes;
    cc._columnPointers = tr._rowPointers;
    cc._values = tr._values;
    cc._rowIndexesSorted = true;
    return cc;
  }

  Int32List get columnIndexes => _columnIndexes;

  /// Returns a new matrix that has the same elements as this matrix, but is
  /// in a dense form.
  IntMatrix dense() {
    var dense = new IntMatrix(rows, columns);
    forEachNonZero((int i, int j, int value) {
      dense.set(i, j, get(i, j));
      return value;
    });
    return dense;
  }

  int get(int row, int column) {
    //int k = Sorting.binarySearchFromTo(columnIndexes, column, rowPointers[row], rowPointers[row + 1] - 1);
    int k = _searchFromTo(
        _columnIndexes, column, _rowPointers[row], _rowPointers[row + 1] - 1);

    int v = 0;
    if (k >= 0) {
      v = _values[k];
    }
    return v;
  }

  Int32List get rowPointers => _rowPointers;

  /// Returns a new matrix that is the transpose of this matrix.
  SparseRCIntMatrix transpose() {
    int nnz = _rowPointers[rows];
    var w = new Int32List(columns);
    var rowPointersT = new Int32List(columns + 1);
    var columnIndexesT = new Int32List(nnz);
    var valuesT = new Int32List(nnz);

    for (int p = 0; p < nnz; p++) {
      w[_columnIndexes[p]]++;
    }
    _cumsum(rowPointersT, w, columns);
    int q;
    for (int j = 0; j < rows; j++) {
      int high = _rowPointers[j + 1];
      for (int p = _rowPointers[j]; p < high; p++) {
        columnIndexesT[q = w[_columnIndexes[p]]++] = j;
        valuesT[q] = _values[p];
      }
    }
    var T = new SparseRCIntMatrix(columns, rows);
    T._rowPointers = rowPointersT;
    T._columnIndexes = columnIndexesT;
    T._values = valuesT;
    return T;
  }

  Int32List get values => _values;

  bool get columnIndexesSorted => _columnIndexesSorted;

  AbstractIntMatrix like2D(int rows, int columns) {
    return new SparseRCIntMatrix(rows, columns);
  }

  AbstractIntVector like1D(int size) {
    return new SparseIntVector(size);
  }

  /// Removes (sums) duplicate entries (if any).
  void removeDuplicates() {
    int nz = 0;
    int q, i;
    var w = new Int32List(columns);
    for (i = 0; i < columns; i++) {
      w[i] = -1; // column i not yet seen
    }
    for (int j = 0; j < rows; j++) {
      // row j will start at q
      q = nz;
      for (int p = _rowPointers[j]; p < _rowPointers[j + 1]; p++) {
        // A(i,j) is nonzero
        i = _columnIndexes[p];
        if (w[i] >= q) {
          // A(i,j) is a duplicate
          _values[w[i]] += _values[p];
        } else {
          // record where column i occurs
          w[i] = nz;
          // keep A(i,j)
          _columnIndexes[nz] = i;
          _values[nz++] = _values[p];
        }
      }
      // record start of row j
      _rowPointers[j] = q;
    }
    // finalize A
    _rowPointers[rows] = nz;
  }

  /// Removes zero entries (if any).
  void removeZeroes() {
    int nz = 0;
    for (int j = 0; j < rows; j++) {
      // get current location of row j
      int p = _rowPointers[j];
      // record new location of row j
      _rowPointers[j] = nz;
      for (; p < _rowPointers[j + 1]; p++) {
        if (_values[p] != 0) {
          // keep A(i,j)
          _values[nz] = _values[p];
          _columnIndexes[nz++] = _columnIndexes[p];
        }
      }
    }
    // finalize A
    _rowPointers[rows] = nz;
  }

  void set(int row, int column, int value) {
    //int k = Sorting.binarySearchFromTo(columnIndexes, column, rowPointers[row], rowPointers[row + 1] - 1);
    int k = _searchFromTo(
        _columnIndexes, column, _rowPointers[row], _rowPointers[row + 1] - 1);

    if (k >= 0) {
      // found
      if (value == 0) {
        _remove(row, k);
      } else {
        _values[k] = value;
      }
      return;
    }

    if (value != 0) {
      k = -k - 1;
      _insert(row, column, k, value);
    }
  }

  /// Sorts column indexes.
  void sortColumnIndexes() {
    SparseRCIntMatrix T = transpose();
    setRows(this, T.rows);
    setColumns(this, T.columns);
    _columnIndexes = T._columnIndexes;
    _rowPointers = T._rowPointers;
    _values = T._values;
    T = transpose();
    setRows(this, T.rows);
    setColumns(this, T.columns);
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
        buf
          ..write('(')
          ..write(i)
          ..write(',')
          ..write(_columnIndexes[j])
          ..write(')')
          ..write('\t')
          ..write(_values[j])
          ..write('\n');
      }
    }
    return buf.toString();
  }

  void trimToSize() => _realloc(0);

  AbstractIntVector mult(AbstractIntVector y, [AbstractIntVector z = null,
      final int alpha = 1, int beta = null, final bool transposeA = false]) {
    if (beta == null) {
      beta = z == null ? 1 : 0;
    }
    int rowsA = transposeA ? columns : rows;
    int columnsA = transposeA ? rows : columns;

    bool ignore = (z == null || !transposeA);
    if (z == null) {
      z = new IntVector(rowsA);
    }

    if (!(y is IntVector && z is IntVector)) {
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

    IntVector zz = z as IntVector;
    Int32List elementsZ = zz._elements;
    int strideZ = zz.stride;
    int zeroZ = z.index(0);

    IntVector yy = y as IntVector;
    Int32List elementsY = yy._elements;
    int strideY = yy.stride;
    int zeroY = y.index(0);

    if (transposeA) {
      if ((!ignore) && (beta != 1.0)) {
        z.apply(ifunc.multiply(beta));
      }

      for (int i = 0; i < rows; i++) {
        int high = _rowPointers[i + 1];
        int yElem = alpha * elementsY[zeroY + strideY * i];
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
        int sum = 0;
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
        int sum = 0;
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

  AbstractIntMatrix multiply(AbstractIntMatrix B, [AbstractIntMatrix C = null,
      final int alpha = 1, int beta = null, final bool transposeA = false,
      final bool transposeB = false]) {
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
      if (B is SparseRCIntMatrix) {
        C = new SparseRCIntMatrix(rowsA, p, (rowsA * p));
      } else {
        C = new IntMatrix(rowsA, p);
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

    if ((B is IntMatrix) && (C is IntMatrix)) {
      SparseRCIntMatrix AA;
      if (transposeA) {
        AA = transpose();
      } else {
        AA = this;
      }
      IntMatrix BB;
      if (transposeB) {
        BB = B.dice() as IntMatrix;
      } else {
        BB = B;
      }

      IntMatrix CC = C;
      Int32List rowPointersA = AA._rowPointers;
      Int32List columnIndexesA = AA._columnIndexes;
      Int32List valuesA = AA._values;

      for (int ii = 0; ii < rowsA; ii++) {
        int highA = rowPointersA[ii + 1];
        for (int ka = rowPointersA[ii]; ka < highA; ka++) {
          int scal = valuesA[ka] * alpha;
          int jj = columnIndexesA[ka];
          CC.row(ii).forEachWith(BB.row(jj), ifunc.plusMultSecond(scal));
        }
      }
    } else if ((B is SparseRCIntMatrix) && (C is SparseRCIntMatrix)) {
      SparseRCIntMatrix AA;
      SparseRCIntMatrix BB;
      SparseRCIntMatrix CC = C;
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
      Int32List valuesA = AA._values;

      Int32List rowPointersB = BB._rowPointers;
      Int32List columnIndexesB = BB._columnIndexes;
      Int32List valuesB = BB._values;

      Int32List rowPointersC = CC._rowPointers;
      Int32List columnIndexesC = CC._columnIndexes;
      Int32List valuesC = CC._values;
      int nzmax = valuesC.length;

      Int32List iw = new Int32List(columnsB + 1);
      for (int i = 0; i < iw.length; i++) {
        iw[i] = -1;
      }
      int len = -1;
      for (int ii = 0; ii < rowsA; ii++) {
        int highA = rowPointersA[ii + 1];
        for (int ka = rowPointersA[ii]; ka < highA; ka++) {
          int scal = valuesA[ka] * alpha;
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
        //IntMatrix1D columnIndexesCPart = columnIndexesC.viewPart(rowPointersC[ii], length);
        //Int32List indexes = cern.colt.matrix.tint.algo.IntSorting.quickSort.sortIndex(columnIndexesCPart);
        //Arrays.sort(columnIndexesCElements, rowPointersC[ii], rowPointersC[ii + 1]);
        //IntMatrix1D valuesCPart = valuesC.viewPart(rowPointersC[ii], length).viewSelection(indexes);
        //valuesC.viewPart(rowPointersC[ii], length).assign(valuesCPart);
      }
      //CC.columnIndexes.elements((Int32List) columnIndexesC.elements);
      //CC.columnIndexes.setSize(columnIndexesSize);
      //CC.values.elements((Int32List) valuesC.elements);
      //CC.values.setSize(columnIndexesSize);
    } else {
      if (transposeB) {
        B = B.dice();
      }
      // cache views
      final List<AbstractIntVector> Brows =
          new List<AbstractIntVector>(columnsA);
      for (int i = columnsA; --i >= 0;) {
        Brows[i] = B.row(i);
      }
      final List<AbstractIntVector> Crows = new List<AbstractIntVector>(rowsA);
      for (int i = rowsA; --i >= 0;) {
        Crows[i] = C.row(i);
      }

      var fn = new ifunc.IntPlusMultSecond.plusMult(0);

      final Int32List columnIndexesA = _columnIndexes;
      final Int32List valuesA = _values;
      for (int i = rows; --i >= 0;) {
        int low = _rowPointers[i];
        for (int k = _rowPointers[i + 1]; --k >= low;) {
          int j = columnIndexesA[k];
          fn.multiplicator = valuesA[k] * alpha;
          if (!transposeA) {
            Crows[i].forEachWith(Brows[j], fn);
          }
          else Crows[j].forEachWith(Brows[i], fn);
        }
      }
    }
    return C;
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
      nzmax = _rowPointers[rows];
    }
    var columnIndexesNew = new Int32List(nzmax);
    columnIndexesNew.setAll(0, _columnIndexes);
    _columnIndexes = columnIndexesNew;
    Int32List valuesNew = new Int32List(nzmax);
    valuesNew.setAll(0, _values);
    _values = valuesNew;
  }

  AbstractIntMatrix _getContent() => this;

  void _insert(int row, int column, int index, int value) {
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
    _values = new Int32List.fromList(valuesList);
  }

  void _remove(int row, int index) {
    var columnIndexesList = new List.from(_columnIndexes);
    columnIndexesList.length = _rowPointers[rows];
    var valuesList = new List.from(_values);
    valuesList.length = _rowPointers[rows];
    columnIndexesList.remove(index);
    valuesList.remove(index);
    for (int i = _rowPointers.length; --i > row;) {
      _rowPointers[i]--;
    }
    _columnIndexes = new Int32List.fromList(columnIndexesList);
    _values = new Int32List.fromList(valuesList);
  }

  Object clone() {
    return new SparseRCIntMatrix._internal(
        rows, columns, _rowPointers, _columnIndexes, _values);
  }
}
