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

/// Dense 2-d matrix holding [double] elements.
///
/// Internally holds one single contigous one-dimensional array, addressed in row
/// major order.
class DenseDoubleMatrix extends DoubleMatrix {
  Float64List _elements;

  /// Constructs a matrix with a given number of rows and columns. All entries
  /// are initially `0`.
  factory DenseDoubleMatrix(int rows, int columns) {
    final elements = new Float64List(rows * columns);
    return new DenseDoubleMatrix._internal(
        rows, columns, elements, 0, 0, columns, 1, false);
  }

  /// Constructs a matrix with the given parameters.
  DenseDoubleMatrix._internal(int rows, int columns, Float64List elements,
      int rowZero, int columnZero, int rowStride, int columnStride, bool isView)
      : super(rows, columns, rowZero, columnZero, rowStride, columnStride,
          !isView) {
    this._elements = elements;
  }

  static DenseDoubleMatrix create(int r, int c) => new DenseDoubleMatrix(r, c);

  factory DenseDoubleMatrix.random(int rows, int columns) {
    return dfactory.randomMatrix(rows, columns, create);
  }

  double aggregate(DoubleDoubleFunction aggr, DoubleFunction fn) {
    if (size == 0) {
      return double.NAN;
    }
    final int zero = index(0, 0);
    double a = 0.0;
    a = fn(_elements[
        zero + (rows - 1) * rowStride + (columns - 1) * columnStride]);
    int d = 1;
    for (int r = rows; --r >= 0;) {
      int ridx = zero + r * rowStride;
      for (int c = columns - d; --c >= 0;) {
        a = aggr(a, fn(_elements[ridx + c * columnStride]));
      }
      d = 0;
    }
    return a;
  }

  void apply(DoubleFunction fn) {
    final Float64List elems = this._elements;
    if (elems == null) {
      throw new Error();
    }
    final int zero = index(0, 0);
    int idx = zero + (rows - 1) * rowStride + (columns - 1) * columnStride;
    // specialization for speed
    if (fn is func.DoubleMult) {
      // x[i] =
      // mult*x[i]
      double multiplicator = fn.multiplicator;
      if (multiplicator == 1) {
        return;
      }
      if (multiplicator == 0) {
        fill(0.0);
        return;
      }
      for (int r = rows; --r >= 0;) {
        // the general case
        for (int i = idx, c = columns; --c >= 0;) {
          elems[i] *= multiplicator;
          i -= columnStride;
        }
        idx -= rowStride;
      }
    } else {
      // the general case x[i] = f(x[i])
      for (int r = rows; --r >= 0;) {
        for (int i = idx, c = columns; --c >= 0;) {
          elems[i] = fn(elems[i]);
          i -= columnStride;
        }
        idx -= rowStride;
      }
    }
  }

  void fill(final double value) {
    final Float64List elems = this._elements;
    final int zero = index(0, 0);
    int idx = zero;
    for (int r = 0; r < rows; r++) {
      for (int i = idx, c = 0; c < columns; c++) {
        elems[i] = value;
        i += columnStride;
      }
      idx += rowStride;
    }
  }

  void setAll(final Float64List values) {
    if (values.length != size) {
      throw new ArgumentError(
          "Must have same length: length=${values.length} rows()*columns()=${rows * columns}");
    }
    if (!isView) {
      _elements.setAll(0, values);
    } else {
      final int zero = index(0, 0);
      int idxOther = 0;
      int idx = zero;
      for (int r = 0; r < rows; r++) {
        for (int i = idx, c = 0; c < columns; c++) {
          _elements[i] = values[idxOther++];
          i += columnStride;
        }
        idx += rowStride;
      }
    }
  }

  void copyFrom(final DoubleMatrix source) {
    // overriden for performance only
    if (source is! DenseDoubleMatrix) {
      super.copyFrom(source);
      return;
    }
    DenseDoubleMatrix other = source as DenseDoubleMatrix;
    if (other == this) {
      return; // nothing to do
    }

    checkShape(this, other);

    if (!isView && !other.isView) {
      _elements.setAll(0, other._elements);
      return;
    }
    if (_haveSharedCells(other)) {
      DoubleMatrix c = other.copy();
      if (c is! DenseDoubleMatrix) {
        // should not happen
        super.copyFrom(other);
        return;
      }
      other = c as DenseDoubleMatrix;
    }

    Float64List elementsOther = other._elements;
    if (_elements == null || elementsOther == null) {
      throw new Error();
    }
    var zeroOther = other.index(0, 0);
    var zero = this.index(0, 0);
    var columnStrideOther = other.columnStride;
    var rowStrideOther = other.rowStride;

    int idx = zero;
    int idxOther = zeroOther;
    for (int r = 0; r < rows; r++) {
      for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
        _elements[i] = elementsOther[j];
        i += columnStride;
        j += columnStrideOther;
      }
      idx += rowStride;
      idxOther += rowStrideOther;
    }
  }

  void assign(final DoubleMatrix y, DoubleDoubleFunction fn) {
    // overriden for performance only
    if (!(y is DenseDoubleMatrix)) {
      super.assign(y, fn);
      return;
    }
    DenseDoubleMatrix other = y as DenseDoubleMatrix;
    checkShape(this, y);
    Float64List elementsOther = other._elements;
    if (_elements == null || elementsOther == null) {
      throw new Error();
    }
    var zeroOther = other.index(0, 0);
    var zero = this.index(0, 0);
    var columnStrideOther = other.columnStride;
    var rowStrideOther = other.rowStride;

    int idx;
    int idxOther;
    // specialized for speed
    if (fn == func.mult) {
      // x[i] = x[i] * y[i]
      idx = zero;
      idxOther = zeroOther;
      for (int r = 0; r < rows; r++) {
        for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
          _elements[i] *= elementsOther[j];
          i += columnStride;
          j += columnStrideOther;
        }
        idx += rowStride;
        idxOther += rowStrideOther;
      }
    } else if (fn == func.div) {
      // x[i] = x[i] / y[i]
      idx = zero;
      idxOther = zeroOther;
      for (int r = 0; r < rows; r++) {
        for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
          _elements[i] /= elementsOther[j];
          i += columnStride;
          j += columnStrideOther;
        }
        idx += rowStride;
        idxOther += rowStrideOther;
      }
    } else if (fn is func.DoublePlusMultSecond) {
      double multiplicator = fn.multiplicator;
      if (multiplicator == 0) {
        // x[i] = x[i] + 0*y[i]
        return;
      } else if (multiplicator == 1) {
        // x[i] = x[i] + y[i]
        idx = zero;
        idxOther = zeroOther;
        for (int r = 0; r < rows; r++) {
          for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
            _elements[i] += elementsOther[j];
            i += columnStride;
            j += columnStrideOther;
          }
          idx += rowStride;
          idxOther += rowStrideOther;
        }
      } else if (multiplicator == -1) {
        // x[i] = x[i] - y[i]
        idx = zero;
        idxOther = zeroOther;
        for (int r = 0; r < rows; r++) {
          for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
            _elements[i] -= elementsOther[j];
            i += columnStride;
            j += columnStrideOther;
          }
          idx += rowStride;
          idxOther += rowStrideOther;
        }
      } else {
        // the general case
        // x[i] = x[i] + mult*y[i]
        idx = zero;
        idxOther = zeroOther;
        for (int r = 0; r < rows; r++) {
          for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
            _elements[i] += multiplicator * elementsOther[j];
            i += columnStride;
            j += columnStrideOther;
          }
          idx += rowStride;
          idxOther += rowStrideOther;
        }
      }
    } else {
      // the general case x[i] = f(x[i],y[i])
      idx = zero;
      idxOther = zeroOther;
      for (int r = 0; r < rows; r++) {
        for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
          _elements[i] = fn(_elements[i], elementsOther[j]);
          i += columnStride;
          j += columnStrideOther;
        }
        idx += rowStride;
        idxOther += rowStrideOther;
      }
    }
  }

  int get cardinality {
    int cardinality = 0;
    var zero = index(0, 0);
    int idx = zero;
    for (int r = 0; r < rows; r++) {
      for (int i = idx, c = 0; c < columns; c++) {
        if (_elements[i] != 0) cardinality++;
        i += columnStride;
      }
      idx += rowStride;
    }
    return cardinality;
  }

  dynamic get elements => _elements;

  void forEachNonZero(func.IntIntDoubleFunction fn) {
    var zero = index(0, 0);
    int idx = zero;
    for (int r = 0; r < rows; r++) {
      for (int i = idx, c = 0; c < columns; c++) {
        double value = _elements[i];
        if (value != 0) {
          _elements[i] = fn(r, c, value);
        }
        i += columnStride;
      }
      idx += rowStride;
    }
  }

  DoubleMatrixLocation max() {
    int rowLocation = 0;
    int columnLocation = 0;
    var zero = index(0, 0);
    double maxValue = _elements[zero];
    int d = 1;
    for (int r = 0; r < rows; r++) {
      for (int c = d; c < columns; c++) {
        var elem = _elements[zero + r * rowStride + c * columnStride];
        if (maxValue < elem) {
          maxValue = elem;
          rowLocation = r;
          columnLocation = c;
        }
      }
      d = 0;
    }
    return new DoubleMatrixLocation._(maxValue, rowLocation, columnLocation);
  }

  DoubleMatrixLocation min() {
    int rowLocation = 0;
    int columnLocation = 0;
    var zero = index(0, 0);
    double minValue = _elements[zero];
    int d = 1;
    double elem;
    for (int r = 0; r < rows; r++) {
      for (int c = d; c < columns; c++) {
        elem = _elements[zero + r * rowStride + c * columnStride];
        if (minValue > elem) {
          minValue = elem;
          rowLocation = r;
          columnLocation = c;
        }
      }
      d = 0;
    }
    return new DoubleMatrixLocation._(minValue, rowLocation, columnLocation);
  }

  void negative(final List<int> rowList, final List<int> columnList,
      final List<double> valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    int idx = this.index(0, 0);
    for (int r = 0; r < rows; r++) {
      for (int i = idx, c = 0; c < columns; c++) {
        double value = _elements[i];
        if (value < 0) {
          rowList.add(r);
          columnList.add(c);
          valueList.add(value);
        }
        i += columnStride;
      }
      idx += rowStride;
    }
  }

  void nonzero(final List<int> rowList, final List<int> columnList,
      final List<double> valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    int idx = this.index(0, 0);
    for (int r = 0; r < rows; r++) {
      for (int i = idx, c = 0; c < columns; c++) {
        double value = _elements[i];
        if (value != 0) {
          rowList.add(r);
          columnList.add(c);
          valueList.add(value);
        }
        i += columnStride;
      }
      idx += rowStride;
    }
  }

  void positive(final List<int> rowList, final List<int> columnList,
      final List<double> valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    int idx = this.index(0, 0);
    for (int r = 0; r < rows; r++) {
      for (int i = idx, c = 0; c < columns; c++) {
        double value = _elements[i];
        if (value > 0) {
          rowList.add(r);
          columnList.add(c);
          valueList.add(value);
        }
        i += columnStride;
      }
      idx += rowStride;
    }
  }

  double get(int row, int column) {
    return _elements[
        rowZero + row * rowStride + columnZero + column * columnStride];
  }

  int index(int row, int column) {
    return rowZero + row * rowStride + columnZero + column * columnStride;
  }

  DoubleMatrix like2D(int rows, int columns) {
    return new DenseDoubleMatrix(rows, columns);
  }

  DoubleVector like1D(int size) => new DenseDoubleVector(size);

  void set(int row, int column, double value) {
    _elements[rowZero + row * rowStride + columnZero + column * columnStride] =
        value;
  }

  DoubleVector mult(final DoubleVector y,
      [DoubleVector z = null, final double alpha = 1.0,
      final double beta = 0.0, final bool transposeA = false]) {
    if (transposeA) {
      return dice().mult(y, z, alpha, beta, false);
    }
    if (z == null) {
      z = new DenseDoubleVector(rows);
    }
    if (!(y is DenseDoubleVector && z is DenseDoubleVector)) {
      return super.mult(y, z, alpha, beta, transposeA);
    }

    if (columns != y.size || rows > z.size) {
      throw new ArgumentError("Incompatible args: " +
          toStringShort() +
          ", " +
          y.toStringShort() +
          ", " +
          z.toStringShort());
    }

    var elemsY = y.elements as Float64List;
    var elemsZ = z.elements as Float64List;
    if (_elements == null || elemsY == null || elemsZ == null) {
      throw new Error();
    }
    var strideY = y.stride;
    var strideZ = z.stride;
    var zero = index(0, 0);
    var zeroY = y.index(0);
    var zeroZ = z.index(0);

    int idxZero = zero;
    int idxZeroZ = zeroZ;
    for (int r = 0; r < rows; r++) {
      double sum = 0.0;
      int idx = idxZero;
      int idxY = zeroY;
      for (int c = 0; c < columns; c++) {
        sum += _elements[idx] * elemsY[idxY];
        idx += columnStride;
        idxY += strideY;
      }
      elemsZ[idxZeroZ] = alpha * sum + beta * elemsZ[idxZeroZ];
      idxZero += rowStride;
      idxZeroZ += strideZ;
    }
    return z;
  }

  DoubleMatrix multiply(DoubleMatrix B,
      [DoubleMatrix C = null, final double alpha = 1.0,
      final double beta = 0.0, final bool transposeA = false,
      bool transposeB = false]) {
    if (transposeA) {
      return dice().multiply(B, C, alpha, beta, false, transposeB);
    }
    if (B is SparseDoubleMatrix ||
        B is SparseRCDoubleMatrix ||
        B is SparseCCDoubleMatrix) {
      // exploit quick sparse mult
      // A*B = (B' * A')'
      if (C == null) {
        return B.multiply(this, null, alpha, beta, !transposeB, true).dice();
      } else {
        B.multiply(this, C.dice(), alpha, beta, !transposeB, true);
        return C;
      }
    }
    if (transposeB) {
      return multiply(B.dice(), C, alpha, beta, transposeA, false);
    }

    int rowsA = rows;
    int columnsA = columns;
    int p = B.columns;
    if (C == null) {
      C = new DenseDoubleMatrix(rowsA, p);
    }
    if (B is! DenseDoubleMatrix || C is! DenseDoubleMatrix) {
      return super.multiply(B, C, alpha, beta, transposeA, transposeB);
    }
    if (B.rows != columnsA) {
      throw new ArgumentError("Matrix inner dimensions must agree:" +
          toStringShort() +
          ", " +
          B.toStringShort());
    }
    if (C.rows != rowsA || C.columns != p) {
      throw new ArgumentError("Incompatibel result matrix: " +
          toStringShort() +
          ", " +
          B.toStringShort() +
          ", " +
          C.toStringShort());
    }
    if (this == C || B == C) {
      throw new ArgumentError("Matrices must not be identical");
    }

    DenseDoubleMatrix BB = B as DenseDoubleMatrix;
    DenseDoubleMatrix CC = C as DenseDoubleMatrix;
    Float64List AElems = this._elements;
    Float64List BElems = BB._elements;
    Float64List CElems = CC._elements;
    if (AElems == null || BElems == null || CElems == null) {
      throw new Error();
    }

    int cA = columnStride;
    int cB = BB.columnStride;
    int cC = CC.columnStride;

    int rA = rowStride;
    int rB = BB.rowStride;
    int rC = CC.rowStride;

    // A is blocked to hide memory latency xxxxxxx B xxxxxxx xxxxxxx A xxx
    // xxxxxxx C xxx xxxxxxx --- ------- xxx xxxxxxx xxx xxxxxxx --- -------
    // xxx xxxxxxx
    const BLOCK_SIZE = 30000; // * 8 == Level 2 cache in bytes
    int m_optimal = (BLOCK_SIZE - columnsA) ~/ (columnsA + 1);
    if (m_optimal <= 0) m_optimal = 1;
    int blocks = rowsA ~/ m_optimal;
    int rr = 0;
    if (rowsA % m_optimal != 0) blocks++;
    for (; --blocks >= 0;) {
      int jB = BB.index(0, 0);
      int indexA = index(rr, 0);
      int jC = CC.index(rr, 0);
      rr += m_optimal;
      if (blocks == 0) m_optimal += rowsA - rr;

      for (int j = p; --j >= 0;) {
        int iA = indexA;
        int iC = jC;
        for (int i = m_optimal; --i >= 0;) {
          int kA = iA;
          int kB = jB;
          double s = 0.0;

          // loop unrolled
          kA -= cA;
          kB -= rB;

          for (int k = columnsA % 4; --k >= 0;) {
            s += AElems[kA += cA] * BElems[kB += rB];
          }
          for (int k = columnsA ~/ 4; --k >= 0;) {
            s += AElems[kA += cA] * BElems[kB += rB] +
                AElems[kA += cA] * BElems[kB += rB] +
                AElems[kA += cA] * BElems[kB += rB] +
                AElems[kA += cA] * BElems[kB += rB];
          }

          CElems[iC] = alpha * s + beta * CElems[iC];
          iA += rA;
          iC += rC;
        }
        jB += cB;
        jC += cC;
      }
    }
    return C;
  }

  double sum() {
    double sum = 0.0;
    if (_elements == null) {
      throw new Error();
    }
    var zero = index(0, 0);

    int idx = zero;
    for (int r = 0; r < rows; r++) {
      for (int i = idx, c = 0; c < columns; c++) {
        sum += _elements[i];
        i += columnStride;
      }
      idx += rowStride;
    }
    return sum;
  }

  bool _haveSharedCellsRaw(DoubleMatrix other) {
    if (other is SelectedDenseDoubleMatrix) {
      return _elements == other._elements;
    } else if (other is DenseDoubleMatrix) {
      return _elements == other._elements;
    }
    return false;
  }

  DoubleVector _like1D(int size, int zero, int stride) {
    return new DenseDoubleVector._internal(size, _elements, zero, stride, true);
  }

  DoubleMatrix _viewSelectionLike(
      Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedDenseDoubleMatrix._(
        _elements, rowOffsets, columnOffsets, 0);
  }

  Object clone() {
    return new DenseDoubleMatrix._internal(rows, columns, _elements, rowZero,
        columnZero, rowStride, columnStride, isView);
  }
}

/// Selection view on dense 2-d matrices holding <tt>double</tt> elements.
///
/// Objects of this class are typically constructed via `select` methods on
/// some source matrix. The interface introduced in abstract super classes
/// defines everything a user can do. From a user point of view there is
/// nothing special about this class; it presents the same functionality with the
/// same signatures and semantics as its abstract superclass(es) while
/// introducing no additional functionality.
///
/// This class uses no delegation. Its instances point directly to the data.
class SelectedDenseDoubleMatrix extends DoubleMatrix {
  Float64List _elements;

  /// The offsets of the visible cells of this matrix.
  Int32List _rowOffsets, _columnOffsets;

  int _offset;

  /// Constructs a matrix view with the given parameters.
  factory SelectedDenseDoubleMatrix._(Float64List elements,
      Int32List rowOffsets, Int32List columnOffsets, int offset) {
    return new SelectedDenseDoubleMatrix._internal(rowOffsets.length,
        columnOffsets.length, elements, 0, 0, 1, 1, rowOffsets, columnOffsets,
        offset, true);
  }

  /// Constructs a matrix view with the given parameters.
  SelectedDenseDoubleMatrix._internal(int rows, int columns,
      Float64List elements, int rowZero, int columnZero, int rowStride,
      int columnStride, Int32List rowOffsets, Int32List columnOffsets,
      int offset, bool isView)
      : super(rows, columns, rowZero, columnZero, rowStride, columnStride,
          !isView) {
    _elements = elements;
    _rowOffsets = rowOffsets;
    _columnOffsets = columnOffsets;
    _offset = offset;
  }

  Object get elements => _elements;

  double get(int row, int column) {
    return _elements[_offset +
        _rowOffsets[rowZero + row * rowStride] +
        _columnOffsets[columnZero + column * columnStride]];
  }

  int index(int row, int column) {
    return this._offset +
        _rowOffsets[rowZero + row * rowStride] +
        _columnOffsets[columnZero + column * columnStride];
  }

  DoubleMatrix like2D(int rows, int columns) {
    return new DenseDoubleMatrix(rows, columns);
  }

  DoubleVector like1D(int size) => new DenseDoubleVector(size);

  void set(int row, int column, double value) {
    _elements[_offset +
        _rowOffsets[rowZero + row * rowStride] +
        _columnOffsets[columnZero + column * columnStride]] = value;
  }

  DoubleVector column(int column) {
    checkColumn(this, column);
    int viewSize = rows;
    int viewZero = rowZero;
    int viewStride = rowStride;
    Int32List viewOffsets = _rowOffsets;
    int viewOffset = _offset + columnOffset(columnRank(column));
    return new SelectedDenseDoubleVector._internal(
        viewSize, _elements, viewZero, viewStride, viewOffsets, viewOffset);
  }

  DoubleVector row(int row) {
    checkRow(this, row);
    int viewSize = columns;
    int viewZero = columnZero;
    int viewStride = columnStride;
    Int32List viewOffsets = _columnOffsets;
    int viewOffset = _offset + rowOffset(rowRank(row));
    return new SelectedDenseDoubleVector._internal(
        viewSize, _elements, viewZero, viewStride, viewOffsets, viewOffset);
  }

  int columnOffset(int absRank) => _columnOffsets[absRank];

  int rowOffset(int absRank) => _rowOffsets[absRank];

  bool _haveSharedCellsRaw(DoubleMatrix other) {
    if (other is SelectedDenseDoubleMatrix) {
      return _elements == other._elements;
    } else if (other is DenseDoubleMatrix) {
      return _elements == other._elements;
    }
    return false;
  }

  DoubleVector _like1D(int size, int zero, int stride) {
    // This method is never called since `row` and `column` are overridden.
    throw new Error();
  }

  DoubleMatrix dice() {
    var v = _view();
    vDice(v);
    // swap
    Int32List tmp = _rowOffsets;
    _rowOffsets = _columnOffsets;
    _columnOffsets = tmp;
    return v;
  }

  DoubleMatrix _viewSelectionLike(
      Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedDenseDoubleMatrix._(
        _elements, rowOffsets, columnOffsets, _offset);
  }

  Object clone() {
    return new SelectedDenseDoubleMatrix._internal(rows, columns, _elements,
        rowZero, columnZero, rowStride, columnStride, _rowOffsets,
        _columnOffsets, _offset, isView);
  }
}
