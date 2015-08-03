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
        rows, columns, rowPointers, columnIndexes, values);
//    _dcs = cs.spalloc(rows, columns, nzmax, true, false);
  }

  factory SparseCCDoubleMatrix._internal(int rows, int columns,
      Int32List rowIndexes, Int32List columnPointers, Float64List values) {
    /*try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if (!"matrix too large".equals(exc.getMessage())) throw exc;
    }
    if (columnPointers.length != columns + 1) {
      throw new ArgumentError("columnsPointers.length != columns + 1");
    }*/
    cs.Matrix dcs = new cs.Matrix();
    dcs.m = rows;
    dcs.n = columns;
    dcs.i = rowIndexes;
    dcs.p = columnPointers;
    dcs.x = values;
    dcs.nz = -1; // column-compressed
    dcs.nzmax = values.length;
    return new SparseCCDoubleMatrix._internal(dcs);
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
    cs.Matrix dcs = cs.spalloc(rows, columns, nz, true, false);
    Int32List w = new Int32List(columns);
    Int32List Cp = dcs.p;
    Int32List Ci = dcs.i;
    Float64List Cx = dcs.x;
    for (int k = 0; k < nz; k++) {
      w[columnIndexes[k]]++;
    }
    cs.cumsum(Cp, w, columns);
    int p;
    for (int k = 0; k < nz; k++) {
      Ci[p = w[columnIndexes[k]]++] = rowIndexes[k];
      if (Cx != null) {
        Cx[p] = value;
      }
    }
    if (removeDuplicates) {
      if (!cs.dupl(dcs)) {
        //remove duplicates
        throw new ArgumentError("Exception occured in cs_dupl()!");
      }
    }
    bool rowIndexesSorted = false;
    if (sortRowIndexes) {
      //sort row indexes
      dcs = cs.transpose(dcs, true);
      dcs = cs.transpose(dcs, true);
      if (dcs == null) {
        throw new ArgumentError("Exception occured in cs_transpose()!");
      }
      rowIndexesSorted = true;
    }
    return new SparseCCDoubleMatrix._internal(dcs)
      .._rowIndexesSorted = rowIndexesSorted;
  }

  /// Constructs a matrix with indexes and values given in the coordinate
  /// format.
  factory SparseCCDoubleMatrix.withValues(int rows, int columns,
      Int32List rowIndexes, Int32List columnIndexes, Float64List values,
      [bool removeDuplicates = false, bool removeZeroes = false,
      bool sortRowIndexes = false]) {
    /*try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if (!"matrix too large".equals(exc.getMessage())) throw exc;
    }
    if (rowIndexes.length != columnIndexes.length) {
      throw new ArgumentError("rowIndexes.length != columnIndexes.length");
    } else if (rowIndexes.length != values.length) {
      throw new ArgumentError("rowIndexes.length != values.length");
    }*/
    int nz = Math.max(rowIndexes.length, 1);
    cs.Matrix dcs = cs.spalloc(rows, columns, nz, true, false);
    Int32List w = new Int32List(columns);
    Int32List Cp = dcs.p;
    Int32List Ci = dcs.i;
    Float64List Cx = dcs.x;
    for (int k = 0; k < nz; k++) w[columnIndexes[k]]++;
    cs.cumsum(Cp, w, columns);
    int p;
    for (int k = 0; k < nz; k++) {
      Ci[p = w[columnIndexes[k]]++] = rowIndexes[k];
      if (Cx != null) {
        Cx[p] = values[k];
      }
    }
    if (removeZeroes) {
      cs.dropzeros(dcs); //remove zeroes
    }
    if (removeDuplicates) {
      if (!cs.dupl(dcs)) {
        //remove duplicates
        throw new ArgumentError("Exception occured in cs_dupl()!");
      }
    }
    bool rowIndexesSorted = false;
    if (sortRowIndexes) {
      dcs = cs.transpose(dcs, true);
      dcs = cs.transpose(dcs, true);
      if (dcs == null) {
        throw new ArgumentError("Exception occured in cs_transpose()!");
      }
      rowIndexesSorted = true;
    }
    return new SparseCCDoubleMatrix._internal(dcs)
      .._rowIndexesSorted = rowIndexesSorted;
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

      final Float64List valuesE = _dcs.x;
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
      _dcs.i.fillRange(0, _dcs.i.length, 0);
      _dcs.p.fillRange(0, _dcs.p.length, 0);
      _dcs.x.fillRange(0, _dcs.x.length, 0.0);
    } else {
      int nnz = cardinality;
      for (int i = 0; i < nnz; i++) {
        _dcs.x[i] = value;
      }
    }
  }

  void copyFrom(AbstractDoubleMatrix source) {
    if (source == this) {
      return; // nothing to do
    }
    checkShape(this, source);

    if (source is SparseCCDoubleMatrix) {
      SparseCCDoubleMatrix other = source;
      this._dcs.p.setAll(0, other.columnPointers);
      //System.arraycopy(other.getColumnPointers(), 0, this._dcs.p, 0, _columns + 1);
      int nzmax = other.rowIndexes.length;
      if (_dcs.nzmax < nzmax) {
        _dcs.i = new Int32List(nzmax);
        _dcs.x = new Float64List(nzmax);
      }
      this._dcs.i.setAll(0, other.rowIndexes);
      //System.arraycopy(other.getRowIndexes(), 0, this._dcs.i, 0, nzmax);
      this._dcs.x.setAll(0, other.values);
      //System.arraycopy(other.getValues(), 0, this._dcs.x, 0, nzmax);
      _rowIndexesSorted = other._rowIndexesSorted;
    } else if (source is SparseRCDoubleMatrix) {
      SparseRCDoubleMatrix other = source.transpose();
      this._dcs.p = other.rowPointers;
      this._dcs.i = other.columnIndexes;
      this._dcs.x = other.values;
      this._dcs.nzmax = this._dcs.x.length;
      _rowIndexesSorted = true;
    } else {
      fill(0.0);
      source.forEachNonZero((int i, int j, double value) {
        set(i, j, value);
        return value;
      });
    }
  }

  void assign(final AbstractDoubleMatrix y, func.DoubleDoubleFunction fn) {
    checkShape(this, y);

    if ((y is SparseCCDoubleMatrix) && (fn == func.plus)) {
      // x[i] = x[i] + y[i]
      SparseCCDoubleMatrix yy = y;
      _dcs = cs.add(_dcs, yy._dcs, 1.0, 1.0);
      return;
    }

    if (fn is DoublePlusMultSecond) {
      // x[i] = x[i] + alpha*y[i]
      final double alpha = (fn as DoublePlusMultSecond).multiplicator;
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
      final double alpha = (fn as DoublePlusMultFirst).multiplicator;
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
      final Int32List rowIndexesA = _dcs.i;
      final Int32List columnPointersA = _dcs.p;
      final Float64List valuesA = _dcs.x;
      for (int j = _columns; --j >= 0;) {
        int low = columnPointersA[j];
        for (int k = columnPointersA[j + 1]; --k >= low;) {
          int i = rowIndexesA[k];
          valuesA[k] *= y.get(i, j);
          if (valuesA[k] == 0) _remove(i, j);
        }
      }
      return;
    }

    if (fn == func.div) {
      // x[i] = x[i] / y[i]
      final Int32List rowIndexesA = _dcs.i;
      final Int32List columnPointersA = _dcs.p;
      final Float64List valuesA = _dcs.x;

      for (int j = _columns; --j >= 0;) {
        int low = columnPointersA[j];
        for (int k = columnPointersA[j + 1]; --k >= low;) {
          int i = rowIndexesA[k];
          valuesA[k] /= y.get(i, j);
          if (valuesA[k] == 0) _remove(i, j);
        }
      }
      return;
    }
    super.forEachWith(y, fn);
  }

  int get cardinality => _dcs.p[_columns];

  Object get elements => _dcs;

  void forEachNonZero(final func.IntIntDoubleFunction function) {
    final Int32List rowIndexesA = _dcs.i;
    final Int32List columnPointersA = _dcs.p;
    final Float64List valuesA = _dcs.x;

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
  DoubleMatrix dense() {
    var dense = new DoubleMatrix(rows, columns);
    forEachNonZero((int i, int j, double value) {
      dense.set(i, j, get(i, j));
      return value;
    });
    return dense;
  }

  double get(int row, int column) {
    //int k = Sorting.binarySearchFromTo(dcs.i, row, dcs.p[column], dcs.p[column + 1] - 1);
    int k = _searchFromTo(_dcs.i, row, _dcs.p[column], _dcs.p[column + 1] - 1);
    double v = 0.0;
    if (k >= 0) {
      v = _dcs.x[k];
    }
    return v;
  }

  /// Returns a new matrix that has the same elements as this matrix, but is
  /// in a row-compressed form.
  SparseRCDoubleMatrix rowCompressed() {
    cs.Matrix dcst = cs.transpose(_dcs, true);
    SparseRCDoubleMatrix rc = new SparseRCDoubleMatrix(_rows, _columns);
    rc._columnIndexes = dcst.i;
    rc._rowPointers = dcst.p;
    rc._values = dcst.x;
    rc._columnIndexesSorted = true;
    return rc;
  }

  Int32List get rowIndexes => _rowIndexes;

  /// Returns a new matrix that is the transpose of this matrix.
  SparseCCDoubleMatrix transpose() {
    cs.Matrix dcst = cs.transpose(_dcs, true);
    SparseCCDoubleMatrix tr = new SparseCCDoubleMatrix(_columns, _rows);
    tr._dcs = dcst;
    return tr;
  }

  Float64List get values => _values;

  bool get rowIndexesSorted => _rowIndexesSorted;

  AbstractDoubleMatrix like2D(int rows, int columns) {
    return new SparseCCDoubleMatrix(rows, columns);
  }

  AbstractDoubleVector like1D(int size) => new SparseDoubleVector(size);

  void set(int row, int column, double value) {
    //int k = Sorting.binarySearchFromTo(dcs.i, row, dcs.p[column], dcs.p[column + 1] - 1);
    int k = _searchFromTo(_dcs.i, row, _dcs.p[column], _dcs.p[column + 1] - 1);

    if (k >= 0) {
      // found
      if (value == 0) {
        _remove(column, k);
      } else {
        _dcs.x[k] = value;
      }
      return;
    }

    if (value != 0) {
      k = -k - 1;
      _insert(row, column, k, value);
    }
  }

  void sortRowIndexes() {
    _dcs = cs.transpose(_dcs, true);
    _dcs = cs.transpose(_dcs, true);
    if (_dcs == null) {
      throw new ArgumentError("Exception occured in cs_transpose()!");
    }
    _rowIndexesSorted = true;
  }

  /// Removes (sums) duplicate entries (if any).
  void removeDuplicates() {
    if (!cs.dupl(_dcs)) {
      //remove duplicates
      throw new ArgumentError("Exception occured in cs_dupl()!");
    }
  }

  /// Removes zero entries (if any)
  void removeZeroes() {
    cs.dropzeros(_dcs); //remove zeroes
  }

  void trimToSize() {
    cs.sprealloc(_dcs, 0);
  }

  String toString() {
    var buf = new StringBuffer();
    buf..write("$rows x $columns sparse matrix, nnz = $cardinality\n");
    for (int i = 0; i < columns; i++) {
      int high = _dcs.p[i + 1];
      for (int j = _dcs.p[i]; j < high; j++) {
        buf
          ..write('(')
          ..write(_dcs.i[j])
          ..write(',')
          ..write(i)
          ..write(')')
          ..write('\t')
          ..write(_dcs.x[j])
          ..write('\n');
      }
    }
    return buf.toString();
  }

  AbstractDoubleVector mult(AbstractDoubleVector y,
      [AbstractDoubleVector z = null, final double alpha = 1.0,
      final double beta = 0.0, final bool transposeA = false]) {
    int rowsA = transposeA ? columns : rows;
    int columnsA = transposeA ? rows : columns;

    bool ignore = (z == null || transposeA);
    if (z == null) {
      z = new DoubleVector(rowsA);
    }

    if (!(y is DoubleVector && z is DoubleVector)) {
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

    var zz = z as DoubleVector;
    Float64List elementsZ = zz._elements;
    int strideZ = zz.stride;
    int zeroZ = zz.index(0);

    DoubleVector yy = y as DoubleVector;
    Float64List elementsY = yy._elements;
    int strideY = yy.stride;
    int zeroY = yy.index(0);

    Int32List rowIndexesA = _dcs.i;
    Int32List columnPointersA = _dcs.p;
    Float64List valuesA = _dcs.x;

    int zidx = zeroZ;
    //int nthreads = ConcurrencyUtils.getNumberOfThreads();
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
      int k = _dcs.p[0];
      for (int i = 0; i < columns; i++) {
        double sum = 0.0;
        int high = _dcs.p[i + 1];
        for (; k + 10 < high; k += 10) {
          int ind = k + 9;
          sum += valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] +
              valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] +
              valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] +
              valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] +
              valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] +
              valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] +
              valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] +
              valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] +
              valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]] +
              valuesA[ind] * elementsY[zeroY + strideY * _dcs.i[ind--]];
        }
        for (; k < high; k++) {
          sum += valuesA[k] * elementsY[_dcs.i[k]];
        }
        elementsZ[zidx] = alpha * sum + beta * elementsZ[zidx];
        zidx += strideZ;
      }
    }
    return z;
  }

  AbstractDoubleMatrix multiply(AbstractDoubleMatrix B,
      [AbstractDoubleMatrix C = null, final double alpha = 1.0,
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
        C = new DoubleMatrix(rowsA, p);
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

    if ((B is DoubleMatrix) && (C is DoubleMatrix)) {
      SparseCCDoubleMatrix AA;
      if (transposeA) {
        AA = transpose();
      } else {
        AA = this;
      }
      DoubleMatrix BB;
      if (transposeB) {
        BB = B.dice() as DoubleMatrix;
      } else {
        BB = B;
      }
      DoubleMatrix CC = C;
      Int32List columnPointersA = AA._dcs.p;
      Int32List rowIndexesA = AA._dcs.i;
      Float64List valuesA = AA._dcs.x;

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
      CC._dcs = cs.multiply(AA._dcs, BB._dcs);
      if (CC._dcs == null) {
        throw new ArgumentError("Exception occured in cs_multiply()");
      }
      if (alpha != 1.0) {
        CC.apply(func.multiply(alpha));
      }
    } else {
      if (transposeB) {
        B = B.dice();
      }
      // cache views
      final List<AbstractDoubleVector> Brows =
          new List<AbstractDoubleVector>(columnsA);
      for (int i = columnsA; --i >= 0;) Brows[i] = B.row(i);
      final List<AbstractDoubleVector> Crows =
          new List<AbstractDoubleVector>(rowsA);
      for (int i = rowsA; --i >= 0;) Crows[i] = C.row(i);

      var fn = new DoublePlusMultSecond.plusMult(0.0);

      final Int32List rowIndexesA = _dcs.i;
      final Int32List columnPointersA = _dcs.p;
      final Float64List valuesA = _dcs.x;
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

  AbstractDoubleMatrix _getContent() => this;

  void _insert(int row, int column, int index, double value) {
    List rowIndexes = new List<int>.from(_dcs.i);
    rowIndexes.lenght = _dcs.p[_columns];
    List values = new List.from(_dcs.x);
    values.length = _dcs.p[_columns];
    rowIndexes.insert(index, row);
    values.insert(index, value);
    for (int i = _dcs.p.length; --i > column;) {
      _dcs.p[i]++;
    }
    _dcs.i = new Int32List.fromList(rowIndexes);
    _dcs.x = new Float64List.fromList(values);
    _dcs.nzmax = rowIndexes.length;
  }

  void _remove(int column, int index) {
    List rowIndexes = new List.from(_dcs.i);
    List values = new List.from(_dcs.x);
    rowIndexes.removeAt(index);
    values.removeAt(index);
    for (int i = _dcs.p.length; --i > column;) {
      _dcs.p[i]--;
    }
    _dcs.i = new Int32List.fromList(rowIndexes);
    _dcs.x = new Float64List.fromList(values);
    _dcs.nzmax = rowIndexes.length;
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

  Object clone() {
    return new SparseCCDoubleMatrix._internal(_dcs);
  }
}
