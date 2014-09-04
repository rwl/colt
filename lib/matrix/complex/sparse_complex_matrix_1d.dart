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
 * Sparse hashed 1-d matrix (aka <i>vector</i>) holding <tt>complex</tt>
 * elements. This implementation uses [Map].
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
class SparseDComplexMatrix1D extends DComplexMatrix1D {

  /*
   * The elements of the matrix.
   */
  Map<int, Float64List> _elements;

  /**
   * Constructs a matrix with a copy of the given values. The values are
   * copied. So subsequent changes in <tt>values</tt> are not reflected in the
   * matrix, and vice-versa.
   *
   * @param values
   *            The values to be filled into the new matrix.
   */
  factory SparseDComplexMatrix1D.from(Float64List values) {
    return new SparseDComplexMatrix1D(values.length)..assignValues(values);
  }

  /**
   * Constructs a matrix with a given number of cells.
   *
   * @param size
   *            the number of cells the matrix shall have.
   * @throws IllegalArgumentException
   *             if <tt>size<0</tt>.
   */
  /*factory SparseDComplexMatrix1D.sized(int size) {
        _setUp(size);
        this._elements = new ConcurrentHashMap<Long, Float64List>(size / 1000);
    }*/

  /**
   * Constructs a matrix view with a given number of parameters.
   *
   * @param size
   *            the number of cells the matrix shall have.
   * @param elements
   *            the cells.
   * @param offset
   *            the index of the first element.
   * @param stride
   *            the number of indexes between any two elements, i.e.
   *            <tt>index(i+1)-index(i)</tt>.
   * @throws IllegalArgumentException
   *             if <tt>size<0</tt>.
   */
  SparseDComplexMatrix1D(int size, [Map<int, Float64List> elements = null, int offset = 0, int stride = 1, bool isNoView = false]) {
    if (elements == null) {
      elements = new Map<int, Float64List>();
    }
    _setUp(size, offset, stride);
    this._elements = elements;
    this._isNoView = isNoView;
  }

  DComplexMatrix1D assignValues(Float64List value) {
    // overriden for performance only
    if (this._isNoView && value[0] == 0 && value[1] == 0) this._elements.clear(); else super.assignValues(value);
    return this;
  }

  int cardinality() {
    if (this._isNoView) {
      return this._elements.length;
    } else {
      return super.cardinality();
    }
  }

  Float64List getQuick(int index) {
    Float64List elem = _elements[_zero + index * _stride];
    if (elem != null) {
      return new Float64List.fromList([elem[0], elem[1]]);
    } else {
      return new Float64List(2);
    }
  }

  Map<int, Float64List> elements() {
    return _elements;
  }

  /**
   * Returns <tt>true</tt> if both matrices share at least one identical cell.
   */
  bool _haveSharedCellsRaw(DComplexMatrix1D other) {
    if (other is SelectedSparseDComplexMatrix1D) {
      return this._elements == other._elements;
    } else if (other is SparseDComplexMatrix1D) {
      return this._elements == other._elements;
    }
    return false;
  }

  int index(int rank) {
    return _zero + rank * _stride;
  }

  DComplexMatrix1D like1D(int size) {
    return new SparseDComplexMatrix1D(size);
  }

  DComplexMatrix2D like2D(int rows, int columns) {
    return new SparseDComplexMatrix2D(rows, columns);
  }

  DComplexMatrix2D reshape(final int rows, final int columns) {
    if (rows * columns != _size) {
      throw new ArgumentError("rows*columns != size");
    }
    final DComplexMatrix2D M = new SparseDComplexMatrix2D(rows, columns);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = columns ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstColumn = j * k;
        final int lastColumn = (j == nthreads - 1) ? columns : firstColumn + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx;
          for (int c = firstColumn; c < lastColumn; c++) {
            idx = c * rows;
            for (int r = 0; r < rows; r++) {
              Float64List elem = getQuick(idx++);
              if ((elem[0] != 0) || (elem[1] != 0)) {
                M.setQuick(r, c, elem);
              }
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idx = 0;
    for (int c = 0; c < columns; c++) {
      for (int r = 0; r < rows; r++) {
        Float64List elem = getQuick(idx++);
        if ((elem[0] != 0) || (elem[1] != 0)) {
          M.setQuick(r, c, elem);
        }
      }
    }
    //}
    return M;
  }

  /*DComplexMatrix3D reshape3D(final int slices, final int rows, final int columns) {
    if (slices * rows * columns != _size) {
      throw new ArgumentError("slices*rows*columns != size");
    }
    final DComplexMatrix3D M = new SparseDComplexMatrix3D(slices, rows, columns);
    int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = slices / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstSlice = j * k;
        final int lastSlice = (j == nthreads - 1) ? slices : firstSlice + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx;
          for (int s = firstSlice; s < lastSlice; s++) {
            for (int c = 0; c < columns; c++) {
              idx = s * rows * columns + c * rows;
              for (int r = 0; r < rows; r++) {
                Float64List elem = getQuick(idx++);
                if ((elem[0] != 0) || (elem[1] != 0)) {
                  M.setQuick(s, r, c, elem);
                }
              }
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {
      int idx = 0;
      for (int s = 0; s < slices; s++) {
        for (int c = 0; c < columns; c++) {
          for (int r = 0; r < rows; r++) {
            Float64List elem = getQuick(idx++);
            if ((elem[0] != 0) || (elem[1] != 0)) {
              M.setQuick(s, r, c, elem);
            }
          }
        }
      }
    }
    return M;
  }*/

  void setQuick(int index, Float64List value) {
    int i = _zero + index * _stride;
    if (value[0] == 0 && value[1] == 0) {
      this._elements.remove(i);
    } else {
      this._elements[i] = value;
    }
  }

  void setPartsQuick(int index, double re, double im) {
    int i = _zero + index * _stride;
    if (re == 0 && im == 0) {
      this._elements.remove(i);
    } else {
      this._elements[i] = new Float64List.fromList([re, im]);
    }
  }

  DComplexMatrix1D _viewSelectionLike(Int32List offsets) {
    return new SelectedSparseDComplexMatrix1D.offsets(this._elements, offsets);
  }

  DoubleMatrix1D getImaginaryPart() {
    final DoubleMatrix1D Im = new SparseDoubleMatrix1D(_size);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int i = firstIdx; i < lastIdx; i++) {
            Im.setQuick(i, getQuick(i)[1]);
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    for (int i = 0; i < _size; i++) {
      Im.set(i, getQuick(i)[1]);
    }
    //}
    return Im;
  }

  DoubleMatrix1D getRealPart() {
    final DoubleMatrix1D Re = new SparseDoubleMatrix1D(_size);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int i = firstIdx; i < lastIdx; i++) {
            Re.setQuick(i, getQuick(i)[0]);
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    for (int i = 0; i < _size; i++) {
      Re.set(i, getQuick(i)[0]);
    }
    //}
    return Re;
  }
  
  Object clone() {
    return new SparseDComplexMatrix1D(_size, _elements, _zero, _stride, _isNoView);
  }
}

/**
 * Selection view on sparse 1-d matrices holding <tt>complex</tt> elements. This
 * implementation uses ConcurrentHashMap
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
class SelectedSparseDComplexMatrix1D extends DComplexMatrix1D {

  /*
   * The elements of the matrix.
   */
  Map<int, Float64List> _elements;

  /**
   * The offsets of visible indexes of this matrix.
   */
  Int32List _offsets;

  /**
   * The offset.
   */
  int __offset;

  /**
   * Constructs a matrix view with the given parameters.
   *
   * @param size
   *            the number of cells the matrix shall have.
   * @param elements
   *            the cells.
   * @param zero
   *            the index of the first element.
   * @param stride
   *            the number of indexes between any two elements, i.e.
   *            <tt>index(i+1)-index(i)</tt>.
   * @param offsets
   *            the offsets of the cells that shall be visible.
   * @param offset
   */
  SelectedSparseDComplexMatrix1D(int size, Map<int, Float64List> elements, int zero, int stride, Int32List offsets, int offset) {
    _setUp(size, zero, stride);

    this._elements = elements;
    this._offsets = offsets;
    this.__offset = offset;
    this._isNoView = false;
  }

  /**
   * Constructs a matrix view with the given parameters.
   *
   * @param elements
   *            the cells.
   * @param indexes
   *            The indexes of the cells that shall be visible.
   */
  factory SelectedSparseDComplexMatrix1D.offsets(Map<int, Float64List> elements, Int32List offsets) {
    return new SelectedSparseDComplexMatrix1D(offsets.length, elements, 0, 1, offsets, 0);
  }

  int _offset(int absRank) {
    return _offsets[absRank];
  }

  Float64List getQuick(int index) {
    return _elements[__offset + _offsets[_zero + index * _stride]];
  }

  Map<int, Float64List> elements() {
    throw new UnsupportedError("This method is not supported.");
  }

  /**
   * Returns <tt>true</tt> if both matrices share at least one identical cell.
   */

  bool _haveSharedCellsRaw(DComplexMatrix1D other) {
    if (other is SelectedSparseDComplexMatrix1D) {
      return this._elements == other._elements;
    } else if (other is SparseDComplexMatrix1D) {
      return this._elements == other._elements;
    }
    return false;
  }

  int index(int rank) {
    // return this.offset + super.index(rank);
    // manually inlined:
    return __offset + _offsets[_zero + rank * _stride];
  }

  DComplexMatrix1D like1D(int size) {
    return new SparseDComplexMatrix1D(size);
  }

  DComplexMatrix2D like2D(int rows, int columns) {
    return new SparseDComplexMatrix2D(rows, columns);
  }

  DComplexMatrix2D reshape(int rows, int columns) {
    throw new UnsupportedError("This method is not supported.");
  }

  /*DComplexMatrix3D reshape3D(int slices, int rows, int columns) {
    throw new UnsupportedError("This method is not supported.");
  }*/

  /**
   * Sets the matrix cell at coordinate <tt>index</tt> to the specified value.
   *
   * <p>
   * Provided with invalid parameters this method may access illegal indexes
   * without throwing any exception. <b>You should only use this method when
   * you are absolutely sure that the coordinate is within bounds.</b>
   * Precondition (unchecked): <tt>index&lt;0 || index&gt;=size()</tt>.
   *
   * @param index
   *            the index of the cell.
   * @param value
   *            the value to be filled into the specified cell.
   */
  void setQuick(int index, Float64List value) {
    int i = __offset + _offsets[_zero + index * _stride];
    if (value[0] == 0 && value[1] == 0) {
      this._elements.remove(i);
    } else {
      this._elements[i] = value;
    }
  }

  /**
   * Sets the matrix cell at coordinate <tt>index</tt> to the specified value.
   *
   * <p>
   * Provided with invalid parameters this method may access illegal indexes
   * without throwing any exception. <b>You should only use this method when
   * you are absolutely sure that the coordinate is within bounds.</b>
   * Precondition (unchecked): <tt>index&lt;0 || index&gt;=size()</tt>.
   *
   * @param index
   *            the index of the cell.
   * @param value
   *            the value to be filled into the specified cell.
   */
  void setPartsQuick(int index, double re, double im) {
    int i = __offset + _offsets[_zero + index * _stride];
    if (re == 0 && im == 0) this._elements.remove(i); else this._elements[i] = new Float64List.fromList([re, im]);
  }

  /**
   * Sets up a matrix with a given number of cells.
   *
   * @param size
   *            the number of cells the matrix shall have.
   */
  void _setUp(int size, [int zero = 0, int stride = 1]) {
    super._setUp(size);
    this._stride = 1;
    this.__offset = 0;
  }

  /**
   * Construct and returns a new selection view.
   *
   * @param offsets
   *            the offsets of the visible elements.
   * @return a new view.
   */
  DComplexMatrix1D _viewSelectionLike(Int32List offsets) {
    return new SelectedSparseDComplexMatrix1D.offsets(this._elements, offsets);
  }

  DoubleMatrix1D getImaginaryPart() {
    throw new UnsupportedError("This method is not supported.");
  }

  DoubleMatrix1D getRealPart() {
    throw new UnsupportedError("This method is not supported.");
  }
  
  Object clone() {
    return new SelectedSparseDComplexMatrix1D(_size, _elements, _zero, _stride, _offsets, __offset);
  }

}
