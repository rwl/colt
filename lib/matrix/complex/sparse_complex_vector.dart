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

/// Sparse hashed 1-d matrix (aka vector) holding `complex` elements in
/// a [Map].
class SparseComplexVector extends AbstractComplexVector {
  Map<int, Float64List> _elements;

  /// Constructs a matrix with a copy of the given [values]. The [values] are
  /// copied.
  factory SparseComplexVector.fromList(Float64List values) {
    return new SparseComplexVector(values.length)..setAll(values);
  }

  /// Constructs a matrix with a given number of cells.
  factory SparseComplexVector(int size) {
    return new SparseComplexVector._internal(size, {}, 0, 1, true);
  }

  SparseComplexVector._internal(int size, Map<int, Float64List> elements,
      int offset, int stride, bool isNoView)
      : super(size, offset, stride, isNoView) {
    _elements = elements;
  }

  static SparseComplexVector create(int size) {
    return new SparseComplexVector(size);
  }

  void fill(double re, double im) {
    // overriden for performance only
    if (!isView && re == 0 && im == 0) {
      _elements.clear();
    } else {
      super.fill(re, im);
    }
  }

  int get cardinality {
    if (!isView) {
      return _elements.length;
    } else {
      return super.cardinality;
    }
  }

  Float64List get(int index) {
    Float64List elem = _elements[zero + index * stride];
    if (elem != null) {
      return new Float64List.fromList([elem[0], elem[1]]);
    } else {
      return new Float64List(2);
    }
  }

  Object get elements => _elements;

  bool _haveSharedCellsRaw(AbstractComplexVector other) {
    if (other is SelectedSparseComplexVector) {
      return _elements == other._elements;
    } else if (other is SparseComplexVector) {
      return _elements == other._elements;
    }
    return false;
  }

  int index(int rank) => zero + rank * stride;

  AbstractComplexVector like1D(int size) => new SparseComplexVector(size);

  AbstractComplexMatrix like2D(int rows, int columns) {
    return new SparseComplexMatrix(rows, columns);
  }

  AbstractComplexMatrix reshape(final int rows, final int columns) {
    if (rows * columns != size) {
      throw new ArgumentError("rows*columns != size");
    }
    var M = new SparseComplexMatrix(rows, columns);
    int idx = 0;
    for (int c = 0; c < columns; c++) {
      for (int r = 0; r < rows; r++) {
        Float64List elem = get(idx++);
        if ((elem[0] != 0) || (elem[1] != 0)) {
          M.set(r, c, elem);
        }
      }
    }
    return M;
  }

  void set(int index, Float64List value) {
    int i = zero + index * stride;
    if (value[0] == 0 && value[1] == 0) {
      _elements.remove(i);
    } else {
      _elements[i] = value;
    }
  }

  void setParts(int index, double re, double im) {
    int i = zero + index * stride;
    if (re == 0 && im == 0) {
      _elements.remove(i);
    } else {
      _elements[i] = new Float64List.fromList([re, im]);
    }
  }

  AbstractComplexVector _viewSelectionLike(Int32List offsets) {
    return new SelectedSparseComplexVector.withOffsets(_elements, offsets);
  }

  AbstractDoubleVector imaginary() {
    var Im = new SparseDoubleVector(size);
    for (int i = 0; i < size; i++) {
      Im.set(i, get(i)[1]);
    }
    return Im;
  }

  AbstractDoubleVector real() {
    var Re = new SparseDoubleVector(size);
    for (int i = 0; i < size; i++) {
      Re.set(i, get(i)[0]);
    }
    return Re;
  }

  Object clone() {
    return new SparseComplexVector._internal(
        size, _elements, zero, stride, !isView);
  }
}

/// Selection view on sparse 1-d matrices holding `complex` elements. This
/// implementation uses [Map].
class SelectedSparseComplexVector extends AbstractComplexVector {
  Map<int, Float64List> _elements;

  /// The offsets of visible indexes of this matrix.
  Int32List _offsets;

  int __offset;

  SelectedSparseComplexVector(int size, Map<int, Float64List> elements,
      int zero, int stride, Int32List offsets, int offset)
      : super(size, zero, stride, false) {
    _elements = elements;
    _offsets = offsets;
    __offset = offset;
  }

  factory SelectedSparseComplexVector.withOffsets(
      Map<int, Float64List> elements, Int32List offsets) {
    return new SelectedSparseComplexVector(
        offsets.length, elements, 0, 1, offsets, 0);
  }

  int _offset(int absRank) => _offsets[absRank];

  Float64List get(int index) {
    return _elements[__offset + _offsets[zero + index * stride]];
  }

  /// This method is not supported.
  Object get elements {
    throw new UnsupportedError("This method is not supported.");
  }

  bool _haveSharedCellsRaw(AbstractComplexVector other) {
    if (other is SelectedSparseComplexVector) {
      return _elements == other._elements;
    } else if (other is SparseComplexVector) {
      return _elements == other._elements;
    }
    return false;
  }

  int index(int rank) => __offset + _offsets[zero + rank * stride];

  AbstractComplexVector like1D(int size) => new SparseComplexVector(size);

  AbstractComplexMatrix like2D(int rows, int columns) {
    return new SparseComplexMatrix(rows, columns);
  }

  /// This method is not supported.
  AbstractComplexMatrix reshape(int rows, int columns) {
    throw new UnsupportedError("This method is not supported.");
  }

  void set(int index, Float64List value) {
    int i = __offset + _offsets[zero + index * stride];
    if (value[0] == 0 && value[1] == 0) {
      _elements.remove(i);
    } else {
      _elements[i] = value;
    }
  }

  void setParts(int index, double re, double im) {
    int i = __offset + _offsets[zero + index * stride];
    if (re == 0 && im == 0) {
      _elements.remove(i);
    } else {
      _elements[i] = new Float64List.fromList([re, im]);
    }
  }

  AbstractComplexVector _viewSelectionLike(Int32List offsets) {
    return new SelectedSparseComplexVector.withOffsets(this._elements, offsets);
  }

  /// This method is not supported.
  AbstractDoubleVector imaginary() {
    throw new UnsupportedError("This method is not supported.");
  }

  /// This method is not supported.
  AbstractDoubleVector real() {
    throw new UnsupportedError("This method is not supported.");
  }

  Object clone() {
    return new SelectedSparseComplexVector(
        size, _elements, zero, stride, _offsets, __offset);
  }
}
