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

/// Sparse hashed 2-d matrix holding `complex` elements.
///
/// This implementation uses [Map].
class SparseComplexMatrix extends AbstractComplexMatrix {
  Map<int, Complex> _elements;

  factory SparseComplexMatrix(int rows, int columns) {
    return new SparseComplexMatrix._internal(
        rows, columns, {}, 0, 0, columns, 1, true);
  }

  SparseComplexMatrix._internal(int rows, int columns,
      Map<int, Complex> elements, int rowZero, int columnZero, int rowStride,
      int columnStride, bool isNoView)
      : super(rows, columns, rowZero, columnZero, rowStride, columnStride,
          isNoView) {
    _elements = elements;
  }

  static SparseComplexMatrix create(int rows, int columns) {
    return new SparseComplexMatrix(rows, columns);
  }

  void fill(double re, double im) {
    // overriden for performance only
    if (!isView && re == 0 && im == 0) {
      _elements.clear();
    } else {
      super.fill(re, im);
    }
  }

  void copyFrom(AbstractComplexMatrix source) {
    // overriden for performance only
    if (!(source is SparseComplexMatrix)) {
      super.copyFrom(source);
      return;
    }
    SparseComplexMatrix other = source as SparseComplexMatrix;
    if (other == this) {
      return; // nothing to do
    }
    checkShape(this, other);

    if (!isView && !other.isView) {
      // quickest
      _elements.clear();
      _elements.addAll(other._elements);
      return;
    }
    super.copyFrom(source);
  }

  void assign(
      final AbstractComplexMatrix y, cfunc.ComplexComplexComplexFunction fn) {
    if (isView) {
      super.assign(y, fn);
      return;
    }

    checkShape(this, y);

    if (fn is cfunc.ComplexPlusMultSecond) {
      // x[i] = x[i] + alpha*y[i]
      Complex alpha = fn.multiplicator;
      if (alpha.real == 0 && alpha.imaginary == 1) {
        return; // nothing to do
      }
      y.forEachNonZero((int i, int j, Complex value) {
        set(i, j, get(i, j) + (alpha * value));
        return value;
      });
      return;
    }
    super.assign(y, fn);
  }

  int get cardinality {
    if (!isView) {
      return _elements.length;
    } else {
      return super.cardinality;
    }
  }

  Complex get(int row, int column) {
    Complex elem = _elements[
        rowZero + row * rowStride + columnZero + column * columnStride];
    if (elem != null) {
      return elem;
    } else {
      return Complex.ZERO;
    }
  }

  dynamic get elements => _elements;

  bool _haveSharedCellsRaw(AbstractComplexMatrix other) {
    if (other is SelectedSparseComplexMatrix) {
      return _elements == other._elements;
    } else if (other is SparseComplexMatrix) {
      return _elements == other._elements;
    }
    return false;
  }

  int index(int row, int column) {
    return rowZero + row * rowStride + columnZero + column * columnStride;
  }

  AbstractComplexMatrix like2D(int rows, int columns) {
    return new SparseComplexMatrix(rows, columns);
  }

  AbstractComplexVector like1D(int size) => new SparseComplexVector(size);

  AbstractComplexVector _like1D(int size, int offset, int stride) {
    return new SparseComplexVector._internal(
        size, this._elements, offset, stride, true);
  }

  void set(int row, int column, Complex value) {
    int index = rowZero + row * rowStride + columnZero + column * columnStride;
    if (value.real == 0 && value.imaginary == 0) {
      _elements.remove(index);
    } else {
      _elements[index] = value;
    }
  }

  void setParts(int row, int column, double re, double im) {
    int index = rowZero + row * rowStride + columnZero + column * columnStride;
    if (re == 0 && im == 0) {
      _elements.remove(index);
    } else {
      _elements[index] = new Complex(re, im);
    }
  }

  AbstractComplexMatrix _viewSelectionLike(
      Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedSparseComplexMatrix(
        this._elements, rowOffsets, columnOffsets, 0);
  }

  AbstractDoubleMatrix imaginary() {
    AbstractDoubleMatrix Im = new SparseDoubleMatrix(rows, columns);
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        Im.set(r, c, get(r, c).imaginary);
      }
    }
    return Im;
  }

  AbstractDoubleMatrix real() {
    final AbstractDoubleMatrix Re = new SparseDoubleMatrix(rows, columns);
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        Re.set(r, c, get(r, c).real);
      }
    }
    return Re;
  }

  Object clone() {
    return new SparseComplexMatrix._internal(rows, columns, _elements, rowZero,
        columnZero, rowStride, columnStride, !isView);
  }
}

/// Selection view on sparse 2-d matrices holding `complex` elements. This
/// implementation uses [Map].
class SelectedSparseComplexMatrix extends AbstractComplexMatrix {
  Map<int, Complex> _elements;

  /// The offsets of the visible cells of this matrix.
  Int32List _rowOffsets, _columnOffsets;

  int _offset;

  SelectedSparseComplexMatrix._internal(int rows, int columns,
      Map<int, Complex> elements, int rowZero, int columnZero, int rowStride,
      int columnStride, Int32List rowOffsets, Int32List columnOffsets,
      int offset)
      : super(
          rows, columns, rowZero, columnZero, rowStride, columnStride, false) {
    _elements = elements;
    _rowOffsets = rowOffsets;
    _columnOffsets = columnOffsets;
    _offset = offset;
  }

  factory SelectedSparseComplexMatrix(Map<int, Complex> elements,
      Int32List rowOffsets, Int32List columnOffsets, int offset) {
    return new SelectedSparseComplexMatrix._internal(rowOffsets.length,
        columnOffsets.length, elements, 0, 0, 1, 1, rowOffsets, columnOffsets,
        offset);
  }

  int columnOffset(int absRank) => _columnOffsets[absRank];

  int rowOffset(int absRank) => _rowOffsets[absRank];

  Complex get(int row, int column) {
    return _elements[_offset +
        _rowOffsets[rowZero + row * rowStride] +
        _columnOffsets[columnZero + column * columnStride]];
  }

  /// This method is not supported.
  Object get elements {
    throw new UnsupportedError("This method is not supported.");
  }

  bool _haveSharedCellsRaw(AbstractComplexMatrix other) {
    if (other is SelectedSparseComplexMatrix) {
      return _elements == other._elements;
    } else if (other is SparseComplexMatrix) {
      return _elements == other._elements;
    }
    return false;
  }

  int index(int row, int column) {
    return this._offset +
        _rowOffsets[rowZero + row * rowStride] +
        _columnOffsets[columnZero + column * columnStride];
  }

  AbstractComplexMatrix like2D(int rows, int columns) {
    return new SparseComplexMatrix(rows, columns);
  }

  AbstractComplexVector like1D(int size) => new SparseComplexVector(size);

  AbstractComplexVector _like1D(int size, int zero, int stride) {
    throw new Error(); // never called since row() and column() are overridden
  }

  void set(int row, int column, Complex value) {
    int index = _offset +
        _rowOffsets[rowZero + row * rowStride] +
        _columnOffsets[columnZero + column * columnStride];

    if (value.real == 0 && value.imaginary == 0) {
      _elements.remove(index);
    } else {
      _elements[index] = value;
    }
  }

  void setParts(int row, int column, double re, double im) {
    int index = _offset +
        _rowOffsets[rowZero + row * rowStride] +
        _columnOffsets[columnZero + column * columnStride];

    if (re == 0 && im == 0) {
      _elements.remove(index);
    } else {
      _elements[index] = new Complex(re, im);
    }
  }

  AbstractComplexMatrix dice() {
    var v = _view();
    vDice(v);
    // swap
    Int32List tmp = v._rowOffsets;
    v._rowOffsets = v._columnOffsets;
    v._columnOffsets = tmp;

    // flips stay unaffected

    setIsNoView(v, false);
    return v;
  }

  AbstractComplexVector column(int column) {
    checkColumn(this, column);
    int viewSize = rows;
    int viewZero = rowZero;
    int viewStride = rowStride;
    Int32List viewOffsets = _rowOffsets;
    int viewOffset = _offset + columnOffset(columnRank(column));
    return new SelectedSparseComplexVector(
        viewSize, _elements, viewZero, viewStride, viewOffsets, viewOffset);
  }

  AbstractComplexVector row(int row) {
    checkRow(this, row);
    int viewSize = columns;
    int viewZero = columnZero;
    int viewStride = columnStride;
    Int32List viewOffsets = _columnOffsets;
    int viewOffset = _offset + rowOffset(rowRank(row));
    return new SelectedSparseComplexVector(
        viewSize, _elements, viewZero, viewStride, viewOffsets, viewOffset);
  }

  AbstractComplexMatrix _viewSelectionLike(
      Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedSparseComplexMatrix(
        _elements, rowOffsets, columnOffsets, _offset);
  }

  /// This method is not supported.
  AbstractDoubleMatrix imaginary() {
    throw new UnsupportedError("This method is not supported.");
  }

  /// This method is not supported.
  AbstractDoubleMatrix real() {
    throw new UnsupportedError("This method is not supported.");
  }

  Object clone() {
    return new SelectedSparseComplexMatrix._internal(rows, columns, _elements,
        rowZero, columnZero, rowStride, columnStride, _rowOffsets,
        _columnOffsets, _offset);
  }
}
