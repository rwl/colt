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

/// Dense 2-d matrix holding [int] elements.
///
/// Internally holds one single contigous one-dimensional array, addressed in row
/// major order.
class IntMatrix extends AbstractIntMatrix {

  /// The elements of this matrix. elements are stored in row major order.
  Int32List _elements;

  /// Constructs a matrix with a given number of rows and columns. All entries
  /// are initially `0`.
  factory IntMatrix(int rows, int columns) {
    var elements = new Int32List(rows * columns);
    return new IntMatrix._internal(
        rows, columns, elements, 0, 0, columns, 1, false);
  }

  IntMatrix._internal(int rows, int columns, Int32List elements, int rowZero,
      int columnZero, int rowStride, int columnStride, bool isView)
      : super(rows, columns, rowZero, columnZero, rowStride, columnStride,
          !isView) {
    _elements = elements;
  }

  static IntMatrix create(int r, int c) => new IntMatrix(r, c);

  factory IntMatrix.random(int r, int c) {
    return ifactory.randomMatrix(r, c, create);
  }

  int aggregate(final ifunc.IntIntFunction aggr, final ifunc.IntFunction fn) {
    if (size == 0) {
      throw new ArgumentError("size == 0");
    }
    final int zero = index(0, 0);
    int a = fn(_elements[zero]);
    int d = 1; // first cell already done
    for (int r = 0; r < rows; r++) {
      for (int c = d; c < columns; c++) {
        a = aggr(a, fn(_elements[zero + r * rowStride + c * columnStride]));
      }
      d = 0;
    }
    return a;
  }

  void apply(final ifunc.IntFunction fn) {
    Int32List elems = _elements;
    if (elems == null) {
      throw new Error();
    }
    int zero = index(0, 0);
    int idx = zero;
    // specialization for speed
    if (fn is ifunc.IntMult) {
      // x[i] =
      // mult*x[i]
      int multiplicator = fn.multiplicator;
      if (multiplicator == 1) {
        return;
      }
      if (multiplicator == 0) {
        fill(0);
      }
      for (int r = 0; r < rows; r++) {
        // the general case
        for (int i = idx, c = 0; c < columns; c++) {
          elems[i] *= multiplicator;
          i += columnStride;
        }
        idx += rowStride;
      }
    } else {
      // the general case x[i] = f(x[i])
      for (int r = 0; r < rows; r++) {
        for (int i = idx, c = 0; c < columns; c++) {
          elems[i] = fn(elems[i]);
          i += columnStride;
        }
        idx += rowStride;
      }
    }
  }

  void fill(final int value) {
    Int32List elems = _elements;
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

  void setAll(final Int32List values) {
    if (values.length != size) {
      throw new ArgumentError("Must have same length: length=${values.length} "
          "rows()*columns()=${rows * columns}");
    }
    if (!isView) {
      _elements.setAll(0, values);
    } else {
      int zero = index(0, 0);
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

  void copyFrom(final AbstractIntMatrix source) {
    // overriden for performance only
    if (source is! IntMatrix) {
      super.copyFrom(source);
      return;
    }
    IntMatrix other_final = source as IntMatrix;
    if (other_final == this) {
      return; // nothing to do
    }
    checkShape(this, other_final);
    if (!isView && other_final.isView) {
      // quickest
      _elements.setAll(0, other_final._elements);
      return;
    }
    IntMatrix other = source as IntMatrix;
    if (_haveSharedCells(other)) {
      AbstractIntMatrix c = other.copy();
      if (c is! IntMatrix) {
        // should not happen
        super.copyFrom(other);
        return;
      }
      other = c as IntMatrix;
    }

    Int32List elemsOther = other._elements;
    if (_elements == null || elemsOther == null) {
      throw new Error();
    }
    int zeroOther = other.index(0, 0);
    int zero = index(0, 0);
    int columnStrideOther = other.columnStride;
    int rowStrideOther = other.rowStride;
    int idx = zero;
    int idxOther = zeroOther;
    for (int r = 0; r < rows; r++) {
      for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
        _elements[i] = elemsOther[j];
        i += columnStride;
        j += columnStrideOther;
      }
      idx += rowStride;
      idxOther += rowStrideOther;
    }
  }

  void assign(final AbstractIntMatrix y, final ifunc.IntIntFunction fn) {
    // overriden for performance only
    if (y is! IntMatrix) {
      super.assign(y, fn);
      return;
    }
    IntMatrix other = y as IntMatrix;
    checkShape(this, y);
    Int32List elemsOther = other._elements;
    if (_elements == null || elemsOther == null) {
      throw new Error();
    }
    int zeroOther = other.index(0, 0);
    int zero = index(0, 0);
    int columnStrideOther = other.columnStride;
    int rowStrideOther = other.rowStride;
    // specialized for speed
    if (fn == ifunc.mult) {
      // x[i] = x[i] * y[i]
      var idx = zero;
      var idxOther = zeroOther;
      for (int r = 0; r < rows; r++) {
        for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
          _elements[i] *= elemsOther[j];
          i += columnStride;
          j += columnStrideOther;
        }
        idx += rowStride;
        idxOther += rowStrideOther;
      }
    } else if (fn == ifunc.div) {
      // x[i] = x[i] / y[i]
      var idx = zero;
      var idxOther = zeroOther;
      for (int r = 0; r < rows; r++) {
        for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
          _elements[i] ~/= elemsOther[j];
          i += columnStride;
          j += columnStrideOther;
        }
        idx += rowStride;
        idxOther += rowStrideOther;
      }
    } else if (fn is ifunc.IntPlusMultSecond) {
      int multiplicator = fn.multiplicator;
      if (multiplicator == 0) {
        // x[i] = x[i] + 0*y[i]
        return;
      } else if (multiplicator == 1) {
        // x[i] = x[i] + y[i]
        var idx = zero;
        var idxOther = zeroOther;
        for (int r = 0; r < rows; r++) {
          for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
            _elements[i] += elemsOther[j];
            i += columnStride;
            j += columnStrideOther;
          }
          idx += rowStride;
          idxOther += rowStrideOther;
        }
      } else if (multiplicator == -1) {
        // x[i] = x[i] - y[i]
        var idx = zero;
        var idxOther = zeroOther;
        for (int r = 0; r < rows; r++) {
          for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
            _elements[i] -= elemsOther[j];
            i += columnStride;
            j += columnStrideOther;
          }
          idx += rowStride;
          idxOther += rowStrideOther;
        }
      } else {
        // the general case
        // x[i] = x[i] + mult*y[i]
        var idx = zero;
        var idxOther = zeroOther;
        for (int r = 0; r < rows; r++) {
          for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
            _elements[i] += multiplicator * elemsOther[j];
            i += columnStride;
            j += columnStrideOther;
          }
          idx += rowStride;
          idxOther += rowStrideOther;
        }
      }
    } else {
      // the general case x[i] = f(x[i],y[i])
      var idx = zero;
      var idxOther = zeroOther;
      for (int r = 0; r < rows; r++) {
        for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
          _elements[i] = fn(_elements[i], elemsOther[j]);
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
    int zero = index(0, 0);
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

  Object get elements => _elements;

  void forEachNonZero(final ifunc.IntIntIntFunction function) {
    final int zero = index(0, 0);
    int idx = zero;
    for (int r = 0; r < rows; r++) {
      for (int i = idx, c = 0; c < columns; c++) {
        int value = _elements[i];
        if (value != 0) {
          _elements[i] = function(r, c, value);
        }
        i += columnStride;
      }
      idx += rowStride;
    }
  }

  void negative(final List<int> rowList, final List<int> columnList,
      final List<int> valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    int idx = index(0, 0);
    for (int r = 0; r < rows; r++) {
      for (int i = idx, c = 0; c < columns; c++) {
        int value = _elements[i];
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
      final List<int> valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    int idx = index(0, 0);
    for (int r = 0; r < rows; r++) {
      for (int i = idx, c = 0; c < columns; c++) {
        int value = _elements[i];
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
      final List<int> valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    int idx = index(0, 0);
    for (int r = 0; r < rows; r++) {
      for (int i = idx, c = 0; c < columns; c++) {
        int value = _elements[i];
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

  int get(int row, int column) {
    return _elements[
        rowZero + row * rowStride + columnZero + column * columnStride];
  }

  int index(int row, int column) {
    return rowZero + row * rowStride + columnZero + column * columnStride;
  }

  AbstractIntMatrix like2D(int rows, int columns) {
    return new IntMatrix(rows, columns);
  }

  AbstractIntVector like1D(int size) => new IntVector(size);

  IntMatrixLocation max() {
    int rowLocation = 0;
    int columnLocation = 0;
    int zero = index(0, 0);
    int maxValue = _elements[zero];
    int d = 1;
    int elem;
    for (int r = 0; r < rows; r++) {
      for (int c = d; c < columns; c++) {
        elem = _elements[zero + r * rowStride + c * columnStride];
        if (maxValue < elem) {
          maxValue = elem;
          rowLocation = r;
          columnLocation = c;
        }
      }
      d = 0;
    }
    return new IntMatrixLocation._(maxValue, rowLocation, columnLocation);
  }

  IntMatrixLocation min() {
    int rowLocation = 0;
    int columnLocation = 0;
    final int zero = index(0, 0);
    int minValue = _elements[zero];
    int d = 1;
    int elem;
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
    return new IntMatrixLocation._(minValue, rowLocation, columnLocation);
  }

  void set(int row, int column, int value) {
    _elements[rowZero + row * rowStride + columnZero + column * columnStride] =
        value;
  }

  AbstractIntVector mult(final AbstractIntVector y, [AbstractIntVector z = null,
      final int alpha = 1, int beta = null, final bool transposeA = false]) {
    if (beta == null) {
      beta = z == null ? 1 : 0;
    }
    if (transposeA) {
      return dice().mult(y, z, alpha, beta, false);
    }
    if (z == null) {
      z = new IntVector(rows);
    }
    if (!(y is IntVector && z is IntVector)) {
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

    Int32List elemsY = y.elements as Int32List;
    Int32List elemsZ = z.elements as Int32List;
    if (_elements == null || elemsY == null || elemsZ == null) {
      throw new Error();
    }
    int strideY = y.stride;
    int strideZ = z.stride;
    int zero = index(0, 0);
    int zeroY = y.index(0);
    int zeroZ = z.index(0);
    int idxZero = zero;
    int idxZeroZ = zeroZ;
    for (int r = 0; r < rows; r++) {
      int sum = 0;
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

  AbstractIntMatrix multiply(AbstractIntMatrix B, [AbstractIntMatrix C = null,
      final int alpha = 1, int beta = null, final bool transposeA = false,
      final bool transposeB = false]) {
    if (beta == null) {
      beta = C == null ? 1 : 0;
    }
    if (transposeA) {
      return dice().multiply(B, C, alpha, beta, false, transposeB);
    }
    // TODO sparse int matrix implementations
    /*if (B is SparseIntMatrix || B is SparseRCIntMatrix || B is SparseCCIntMatrix) {
      // exploit quick sparse mult
      // A*B = (B' * A')'
      if (C == null) {
        return B.multiply(this, null, alpha, beta, !transposeB, true).dice();
      } else {
        B.multiply(this, C.dice(), alpha, beta, !transposeB, true);
        return C;
      }
    }*/
    if (transposeB) {
      return this.multiply(B.dice(), C, alpha, beta, transposeA, false);
    }

    int rowsA = rows;
    int columnsA = columns;
    int p = B.columns;
    if (C == null) {
      C = new IntMatrix(rowsA, p);
    }
    if (!(B is IntMatrix) || !(C is IntMatrix)) {
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

    IntMatrix BB = B as IntMatrix;
    IntMatrix CC = C as IntMatrix;
    final Int32List AElems = this._elements;
    final Int32List BElems = BB._elements;
    final Int32List CElems = CC._elements;
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
    if (m_optimal <= 0) {
      m_optimal = 1;
    }
    int blocks = rowsA ~/ m_optimal;
    int rr = 0;
    if (rowsA % m_optimal != 0) {
      blocks++;
    }
    for (; --blocks >= 0;) {
      int jB = BB.index(0, 0);
      int indexA = index(rr, 0);
      int jC = CC.index(rr, 0);
      rr += m_optimal;
      if (blocks == 0) {
        m_optimal += rowsA - rr;
      }

      for (int j = p; --j >= 0;) {
        int iA = indexA;
        int iC = jC;
        for (int i = m_optimal; --i >= 0;) {
          int kA = iA;
          int kB = jB;
          int s = 0;

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

  int sum() {
    int sum = 0;
    if (_elements == null) {
      throw new Error();
    }
    final int zero = index(0, 0);
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

  bool _haveSharedCellsRaw(AbstractIntMatrix other) {
    if (other is SelectedDenseIntMatrix) {
      return _elements == other._elements;
    } else if (other is IntMatrix) {
      return _elements == other._elements;
    }
    return false;
  }

  AbstractIntVector _like1D(int size, int zero, int stride) {
    return new IntVector._internal(size, _elements, zero, stride, true);
  }

  AbstractIntMatrix _viewSelectionLike(
      Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedDenseIntMatrix.offset(
        _elements, rowOffsets, columnOffsets, 0);
  }

  Object clone() {
    return new IntMatrix._internal(rows, columns, _elements, rowZero,
        columnZero, rowStride, columnStride, isView);
  }
}

/// Selection view on dense 2-d matrices holding [int] elements.
///
/// Objects of this class are typically constructed via `select` methods on
/// some source matrix. The interface introduced in abstract super classes
/// defines everything a user can do. From a user point of view there is
/// nothing special about this class; it presents the same functionality with the
/// same signatures and semantics as its abstract superclass(es) while
/// introducing no additional functionality. Thus, this class need not be visible
/// to users.
class SelectedDenseIntMatrix extends AbstractIntMatrix {
  Int32List _elements;

  /// The offsets of the visible cells of this matrix.
  Int32List _rowOffsets;

  Int32List _columnOffsets;

  int _offset;

  factory SelectedDenseIntMatrix.offset(Int32List elements,
      Int32List rowOffsets, Int32List columnOffsets, int offset) {
    return new SelectedDenseIntMatrix(rowOffsets.length, columnOffsets.length,
        elements, 0, 0, 1, 1, rowOffsets, columnOffsets, offset);
  }

  SelectedDenseIntMatrix(int rows, int columns, Int32List elements, int rowZero,
      int columnZero, int rowStride, int columnStride, Int32List rowOffsets,
      Int32List columnOffsets, int offset)
      : super(
          rows, columns, rowZero, columnZero, rowStride, columnStride, false) {
    _elements = elements;
    _rowOffsets = rowOffsets;
    _columnOffsets = columnOffsets;
    _offset = offset;
  }

  /// This method is not supported.
  Object get elements {
    throw new ArgumentError("This method is not supported.");
  }

  int get(int row, int column) {
    return _elements[_offset +
        _rowOffsets[rowZero + row * rowStride] +
        _columnOffsets[columnZero + column * columnStride]];
  }

  int index(int row, int column) {
    return _offset +
        _rowOffsets[rowZero + row * rowStride] +
        _columnOffsets[columnZero + column * columnStride];
  }

  AbstractIntMatrix like2D(int rows, int columns) {
    return new IntMatrix(rows, columns);
  }

  AbstractIntVector like1D(int size) => new IntVector(size);

  void set(int row, int column, int value) {
    _elements[_offset +
        _rowOffsets[rowZero + row * rowStride] +
        _columnOffsets[columnZero + column * columnStride]] = value;
  }

  AbstractIntVector column(int column) {
    checkColumn(this, column);
    int viewSize = rows;
    int viewZero = rowZero;
    int viewStride = rowStride;
    Int32List viewOffsets = _rowOffsets;
    int viewOffset = _offset + _columnOffset(columnRank(column));
    return new SelectedDenseIntVector._internal(
        viewSize, _elements, viewZero, viewStride, viewOffsets, viewOffset);
  }

  AbstractIntVector row(int row) {
    checkRow(this, row);
    int viewSize = columns;
    int viewZero = columnZero;
    int viewStride = columnStride;
    Int32List viewOffsets = _columnOffsets;
    int viewOffset = _offset + _rowOffset(rowRank(row));
    return new SelectedDenseIntVector._internal(
        viewSize, _elements, viewZero, viewStride, viewOffsets, viewOffset);
  }

  int _columnOffset(int absRank) => _columnOffsets[absRank];

  int _rowOffset(int absRank) => _rowOffsets[absRank];

  bool _haveSharedCellsRaw(AbstractIntMatrix other) {
    if (other is SelectedDenseIntMatrix) {
      return _elements == other._elements;
    } else if (other is IntMatrix) {
      return _elements == other._elements;
    }
    return false;
  }

  AbstractIntVector _like1D(int size, int zero, int stride) {
    // never called since [row] and [column] are overridden
    throw new Error();
  }

  AbstractIntMatrix dice() {
    var v = _view();
    vDice(v);
    // TODO: check override
    Int32List tmp = v._rowOffsets;
    v._rowOffsets = v._columnOffsets;
    v._columnOffsets = tmp;
    setIsNoView(v, false);
    return v;
  }

  AbstractIntMatrix _viewSelectionLike(
      Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedDenseIntMatrix.offset(
        this._elements, rowOffsets, columnOffsets, _offset);
  }

  Object clone() {
    return new SelectedDenseIntMatrix(rows, columns, _elements, rowZero,
        columnZero, rowStride, columnStride, _rowOffsets, _columnOffsets,
        _offset);
  }
}
