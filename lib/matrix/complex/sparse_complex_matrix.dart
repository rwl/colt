/*
Copyright (C) 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
is hereby granted without fee, provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear in supporting documentation.
CERN makes no representations about the suitability of this software for any purpose.
It is provided "as is" without expressed or implied warranty.
 */
part of cern.colt.matrix;

/**
 * Sparse hashed 2-d matrix holding <tt>complex</tt> elements.
 *
 * This implementation uses ConcurrentHashMap
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
class SparseComplexMatrix extends AbstractComplexMatrix {

  /*
   * The elements of the matrix.
   */
  Map<int, Float64List> _elements;

  /**
   * Constructs a matrix with a copy of the given values. <tt>values</tt> is
   * required to have the form <tt>values[row][column]</tt> and have exactly
   * the same number of columns in every row.
   * <p>
   * The values are copied. So subsequent changes in <tt>values</tt> are not
   * reflected in the matrix, and vice-versa.
   *
   * @param values
   *            The values to be filled into the new matrix.
   * @throws IllegalArgumentException
   *             if
   *             <tt>for any 1 &lt;= row &lt; values.length: values[row].length != values[row-1].length</tt>
   *             .
   */
  factory SparseComplexMatrix.fromList(List<Float64List> values) {
    return new SparseComplexMatrix(values.length, values.length == 0 ? 0 : values[0].length)
      ..setAll2D(values);
  }

  /**
   * Constructs a view with the given parameters.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @param elements
   *            the cells.
   * @param rowZero
   *            the position of the first element.
   * @param columnZero
   *            the position of the first element.
   * @param rowStride
   *            the number of elements between two rows, i.e.
   *            <tt>index(i+1,j)-index(i,j)</tt>.
   * @param columnStride
   *            the number of elements between two columns, i.e.
   *            <tt>index(i,j+1)-index(i,j)</tt>.
   * @throws IllegalArgumentException
   *             if
   *             <tt>rows<0 || columns<0 || (double)columns*rows > Integer.MAX_VALUE</tt>
   *             or flip's are illegal.
   */
  SparseComplexMatrix(int rows, int columns, [Map<int, Float64List> elements = null, int rowZero = 0, int columnZero = 0, int rowStride = null, int columnStride = 1, bool isNoView = false]) {
    if (elements == null) {
      elements = new Map<int, Float64List>();
    }
    _setUp(rows, columns, rowZero, columnZero, rowStride, columnStride);
    this._elements = elements;
    this._isNoView = isNoView;
  }

  void fill(double re, double im/*Float64List value*/) {
    // overriden for performance only
    //if (this._isNoView && value[0] == 0 && value[1] == 0) {
    if (this._isNoView && re == 0 && im == 0) {
      this._elements.clear();
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
    checkShape(other);

    if (this._isNoView && other._isNoView) { // quickest
      this._elements.clear();
      this._elements.addAll(other._elements);
      return;
    }
    super.copyFrom(source);
  }

  void forEachMatrix(final AbstractComplexMatrix y, cfunc.ComplexComplexComplexFunction function) {
    if (!this._isNoView) {
      super.forEachWith(y, function);
      return;
    }

    checkShape(y);

    if (function is cfunc.ComplexPlusMultSecond) {
      // x[i] = x[i] + alpha*y[i]
      final Float64List alpha = (function as cfunc.ComplexPlusMultSecond).multiplicator;
      if (alpha[0] == 0 && alpha[1] == 1) {
        return; // nothing to do
      }
      y.forEachNonZero((int i, int j, Float64List value) {
        set(i, j, Complex.plus(get(i, j), Complex.multiply(alpha, value)));
        return value;
      });
      return;
    }
    super.forEachWith(y, function);
  }

  int get cardinality {
    if (this._isNoView) {
      return this._elements.length;
    } else {
      return super.cardinality;
    }
  }

  Float64List get(int row, int column) {
    Float64List elem = this._elements[_rowZero + row * _rowStride + _columnZero + column * _columnStride];
    if (elem != null) {
      return new Float64List.fromList([elem[0], elem[1]]);
    } else {
      return new Float64List(2);
    }
  }

  Object get elements => _elements;

  /**
   * Returns <tt>true</tt> if both matrices share common cells. More formally,
   * returns <tt>true</tt> if at least one of the following conditions is met
   * <ul>
   * <li>the receiver is a view of the other matrix
   * <li>the other matrix is a view of the receiver
   * <li><tt>this == other</tt>
   * </ul>
   */

  bool _haveSharedCellsRaw(AbstractComplexMatrix other) {
    if (other is SelectedSparseComplexMatrix) {
      return this._elements == other._elements;
    } else if (other is SparseComplexMatrix) {
      return this._elements == other._elements;
    }
    return false;
  }

  int index(int row, int column) {
    return _rowZero + row * _rowStride + _columnZero + column * _columnStride;
  }

  AbstractComplexMatrix like2D(int rows, int columns) {
    return new SparseComplexMatrix(rows, columns);
  }

  AbstractComplexVector like1D(int size) {
    return new SparseComplexVector(size);
  }

  AbstractComplexVector _like1D(int size, int offset, int stride) {
    return new SparseComplexVector(size, this._elements, offset, stride);
  }

  void set(int row, int column, Float64List value) {
    int index = _rowZero + row * _rowStride + _columnZero + column * _columnStride;
    if (value[0] == 0 && value[1] == 0) {
      this._elements.remove(index);
    } else {
      this._elements[index] = value;
    }
  }

  AbstractComplexVector vectorize() {
    final SparseComplexVector v = new SparseComplexVector(length);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_rows * _columns >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _columns);
      List<Future> futures = new List<Future>(nthreads);
      int k = _columns ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstColumn = j * k;
        final int lastColumn = (j == nthreads - 1) ? _columns : firstColumn + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = 0;
          for (int c = firstColumn; c < lastColumn; c++) {
            idx = c * _rows;
            for (int r = 0; r < _rows; r++) {
              Float64List elem = getQuick(r, c);
              if ((elem[0] != 0) || (elem[1] != 0)) {
                v.setQuick(idx++, elem);
              }
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idx = 0;
    for (int c = 0; c < _columns; c++) {
      for (int r = 0; r < _rows; r++) {
        Float64List elem = get(r, c);
        if ((elem[0] != 0) || (elem[1] != 0)) {
          v.set(idx++, elem);
        }
      }
    }
    //}
    return v;
  }

  void setParts(int row, int column, double re, double im) {
    int index = _rowZero + row * _rowStride + _columnZero + column * _columnStride;
    if (re == 0 && im == 0) {
      this._elements.remove(index);
    } else {
      this._elements[index] = new Float64List.fromList([re, im]);
    }
  }

  AbstractComplexMatrix _viewSelectionLike(Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedSparseComplexMatrix.withOffsets(this._elements, rowOffsets, columnOffsets, 0);
  }

  AbstractDoubleMatrix imaginary() {
    final AbstractDoubleMatrix Im = new SparseDoubleMatrix(_rows, _columns);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = 0; c < _columns; c++) {
              Im.setQuick(r, c, getQuick(r, c)[1]);
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    for (int r = 0; r < _rows; r++) {
      for (int c = 0; c < _columns; c++) {
        Im.set(r, c, get(r, c)[1]);
      }
    }
    //}

    return Im;
  }

  AbstractDoubleMatrix real() {
    final AbstractDoubleMatrix Re = new SparseDoubleMatrix(_rows, _columns);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = 0; c < _columns; c++) {
              Re.setQuick(r, c, getQuick(r, c)[0]);
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    for (int r = 0; r < _rows; r++) {
      for (int c = 0; c < _columns; c++) {
        Re.set(r, c, get(r, c)[0]);
      }
    }
    //}

    return Re;
  }

  Object clone() {
    return new SparseComplexMatrix(_rows, _columns, _elements, _rowZero, _columnZero, _rowStride, _columnStride, _isNoView);
  }
}

/**
 * Selection view on sparse 2-d matrices holding <tt>complex</tt> elements. This
 * implementation uses ConcurrentHashMap
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
class SelectedSparseComplexMatrix extends AbstractComplexMatrix {

  /* The elements of the matrix. */
  Map<int, Float64List> _elements;

  /** The offsets of the visible cells of this matrix. */
  Int32List _rowOffsets;

  Int32List _columnOffsets;

  /** The offset. */
  int _offset;

  /**
   * Constructs a matrix view with the given parameters.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @param elements
   *            the cells.
   * @param rowZero
   *            the position of the first element.
   * @param columnZero
   *            the position of the first element.
   * @param rowStride
   *            the number of elements between two rows, i.e.
   *            <tt>index(i+1,j)-index(i,j)</tt>.
   * @param columnStride
   *            the number of elements between two columns, i.e.
   *            <tt>index(i,j+1)-index(i,j)</tt>.
   * @param rowOffsets
   *            The row offsets of the cells that shall be visible.
   * @param columnOffsets
   *            The column offsets of the cells that shall be visible.
   * @param offset
   */
  SelectedSparseComplexMatrix(int rows, int columns, Map<int, Float64List> elements, int rowZero, int columnZero, int rowStride, int columnStride, Int32List rowOffsets, Int32List columnOffsets, int offset) {
    // be sure parameters are valid, we do not check...
    _setUp(rows, columns, rowZero, columnZero, rowStride, columnStride);

    this._elements = elements;
    this._rowOffsets = rowOffsets;
    this._columnOffsets = columnOffsets;
    this._offset = offset;

    this._isNoView = false;
  }

  /**
   * Constructs a matrix view with the given parameters.
   *
   * @param elements
   *            the cells.
   * @param rowOffsets
   *            The row offsets of the cells that shall be visible.
   * @param columnOffsets
   *            The column offsets of the cells that shall be visible.
   * @param offset
   */
  factory SelectedSparseComplexMatrix.withOffsets(Map<int, Float64List> elements, Int32List rowOffsets, Int32List columnOffsets, int offset) {
    return new SelectedSparseComplexMatrix(rowOffsets.length, columnOffsets.length, elements, 0, 0, 1, 1, rowOffsets, columnOffsets, offset);
  }

  int _columnOffset(int absRank) {
    return _columnOffsets[absRank];
  }

  int _rowOffset(int absRank) {
    return _rowOffsets[absRank];
  }

  Float64List get(int row, int column) {
    return _elements[_offset + _rowOffsets[_rowZero + row * _rowStride] + _columnOffsets[_columnZero + column * _columnStride]];
  }

  Object get elements {
    throw new UnsupportedError("This method is not supported.");
  }

  /**
   * Returns <tt>true</tt> if both matrices share common cells. More formally,
   * returns <tt>true</tt> if <tt>other != null</tt> and at least one of the
   * following conditions is met
   * <ul>
   * <li>the receiver is a view of the other matrix
   * <li>the other matrix is a view of the receiver
   * <li><tt>this == other</tt>
   * </ul>
   */
  bool _haveSharedCellsRaw(AbstractComplexMatrix other) {
    if (other is SelectedSparseComplexMatrix) {
      return this._elements == other._elements;
    } else if (other is SparseComplexMatrix) {
      return this._elements == other._elements;
    }
    return false;
  }

  int index(int row, int column) {
    return this._offset + _rowOffsets[_rowZero + row * _rowStride] + _columnOffsets[_columnZero + column * _columnStride];
  }

  AbstractComplexMatrix like2D(int rows, int columns) {
    return new SparseComplexMatrix(rows, columns);
  }

  AbstractComplexVector like1D(int size) {
    return new SparseComplexVector(size);
  }

  AbstractComplexVector _like1D(int size, int zero, int stride) {
    throw new Error(); // this method is never called since
    // viewRow() and viewColumn are overridden
    // properly.
  }

  void set(int row, int column, Float64List value) {
    int index = _offset + _rowOffsets[_rowZero + row * _rowStride] + _columnOffsets[_columnZero + column * _columnStride];

    if (value[0] == 0 && value[1] == 0) {
      this._elements.remove(index);
    } else {
      this._elements[index] = value;
    }
  }

  AbstractComplexVector vectorize() {
    throw new UnsupportedError("This method is not supported.");
  }

  void setParts(int row, int column, double re, double im) {
    int index = _offset + _rowOffsets[_rowZero + row * _rowStride] + _columnOffsets[_columnZero + column * _columnStride];

    if (re == 0 && im == 0) this._elements.remove(index); else this._elements[index] = new Float64List.fromList([re, im]);

  }

  void _setUp(int rows, int columns, [int rowZero = 0, int columnZero = 0, int rowStride = null, int columnStride = 1]) {
    super._setUp(rows, columns);
    this._rowStride = 1;
    this._columnStride = 1;
    this._offset = 0;
  }

  void _vDice() {
    super._vDice();
    // swap
    Int32List tmp = _rowOffsets;
    _rowOffsets = _columnOffsets;
    _columnOffsets = tmp;

    // flips stay unaffected

    this._isNoView = false;
  }

  AbstractComplexVector column(int column) {
    _checkColumn(column);
    int viewSize = this._rows;
    int viewZero = this._rowZero;
    int viewStride = this._rowStride;
    Int32List viewOffsets = this._rowOffsets;
    int viewOffset = this._offset + _columnOffset(_columnRank(column));
    return new SelectedSparseComplexVector(viewSize, this._elements, viewZero, viewStride, viewOffsets, viewOffset);
  }

  AbstractComplexVector row(int row) {
    _checkRow(row);
    int viewSize = this._columns;
    int viewZero = _columnZero;
    int viewStride = this._columnStride;
    Int32List viewOffsets = this._columnOffsets;
    int viewOffset = this._offset + _rowOffset(_rowRank(row));
    return new SelectedSparseComplexVector(viewSize, this._elements, viewZero, viewStride, viewOffsets, viewOffset);
  }

  AbstractComplexMatrix _viewSelectionLike(Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedSparseComplexMatrix.withOffsets(this._elements, rowOffsets, columnOffsets, this._offset);
  }

  AbstractDoubleMatrix imaginary() {
    throw new UnsupportedError("This method is not supported.");
  }

  AbstractDoubleMatrix real() {
    throw new UnsupportedError("This method is not supported.");
  }

  Object clone() {
    return new SelectedSparseComplexMatrix(_rows, _columns, _elements, _rowZero, _columnZero, _rowStride, _columnStride, _rowOffsets, _columnOffsets, _offset);
  }

}
