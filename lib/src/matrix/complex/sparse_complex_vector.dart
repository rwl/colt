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

/// Sparse hashed 1-d matrix (aka vector) holding [Complex] elements in
/// a [Map].
class SparseComplexVector extends ComplexVector {
  Map<int, Complex> _elements;

  /// Constructs a matrix with a copy of the given [values]. The [values] are
  /// copied.
  factory SparseComplexVector.fromList(Float64List values) {
    return new SparseComplexVector(values.length)..setAll(values);
  }

  /// Constructs a matrix with a given number of cells.
  factory SparseComplexVector(int size) {
    return new SparseComplexVector._internal(size, {}, 0, 1, true);
  }

  SparseComplexVector._internal(int size, Map<int, Complex> elements,
      int offset, int stride, bool isNoView)
      : super._(size, offset, stride, isNoView) {
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

  Complex get(int index) {
    var elem = _elements[zero + index * stride];
    if (elem != null) {
      return elem;
    } else {
      return Complex.ZERO;
    }
  }

  Object get elements => _elements;

  bool _haveSharedCellsRaw(ComplexVector other) {
    if (other is SelectedSparseComplexVector) {
      return _elements == other._elements;
    } else if (other is SparseComplexVector) {
      return _elements == other._elements;
    }
    return false;
  }

  int index(int rank) => zero + rank * stride;

  ComplexVector like1D(int size) => new SparseComplexVector(size);

  ComplexMatrix like2D(int rows, int columns) {
    return new SparseComplexMatrix(rows, columns);
  }

  ComplexMatrix reshape(int rows, int columns) {
    if (rows * columns != size) {
      throw new ArgumentError("rows*columns != size");
    }
    var M = new SparseComplexMatrix(rows, columns);
    int idx = 0;
    for (int c = 0; c < columns; c++) {
      for (int r = 0; r < rows; r++) {
        var elem = get(idx++);
        if ((elem.real != 0) || (elem.imaginary != 0)) {
          M.set(r, c, elem);
        }
      }
    }
    return M;
  }

  void set(int index, Complex value) {
    int i = zero + index * stride;
    if (value.real == 0 && value.imaginary == 0) {
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
      _elements[i] = new Complex(re, im);
    }
  }

  ComplexVector _viewSelectionLike(Int32List offsets) {
    return new SelectedSparseComplexVector.withOffsets(_elements, offsets);
  }

  DoubleVector imaginary() {
    var Im = new SparseDoubleVector(size);
    for (int i = 0; i < size; i++) {
      Im.set(i, get(i).imaginary);
    }
    return Im;
  }

  DoubleVector real() {
    var Re = new SparseDoubleVector(size);
    for (int i = 0; i < size; i++) {
      Re.set(i, get(i).real);
    }
    return Re;
  }

  Object clone() {
    return new SparseComplexVector._internal(
        size, _elements, zero, stride, !isView);
  }
}

/// Selection view on sparse 1-d matrices holding [Complex] elements. This
/// implementation uses [Map].
class SelectedSparseComplexVector extends ComplexVector {
  Map<int, Complex> _elements;

  /// The offsets of visible indexes of this matrix.
  Int32List _offsets;

  int __offset;

  SelectedSparseComplexVector(int size, Map<int, Complex> elements, int zero,
      int stride, Int32List offsets, int offset)
      : super._(size, zero, stride, false) {
    _elements = elements;
    _offsets = offsets;
    __offset = offset;
  }

  factory SelectedSparseComplexVector.withOffsets(
      Map<int, Complex> elements, Int32List offsets) {
    return new SelectedSparseComplexVector(
        offsets.length, elements, 0, 1, offsets, 0);
  }

  int _offset(int absRank) => _offsets[absRank];

  Complex get(int index) {
    return _elements[__offset + _offsets[zero + index * stride]];
  }

  /// This method is not supported.
  Object get elements {
    throw new UnsupportedError("This method is not supported.");
  }

  bool _haveSharedCellsRaw(ComplexVector other) {
    if (other is SelectedSparseComplexVector) {
      return _elements == other._elements;
    } else if (other is SparseComplexVector) {
      return _elements == other._elements;
    }
    return false;
  }

  int index(int rank) => __offset + _offsets[zero + rank * stride];

  ComplexVector like1D(int size) => new SparseComplexVector(size);

  ComplexMatrix like2D(int rows, int columns) {
    return new SparseComplexMatrix(rows, columns);
  }

  /// This method is not supported.
  ComplexMatrix reshape(int rows, int columns) {
    throw new UnsupportedError("This method is not supported.");
  }

  void set(int index, Complex value) {
    int i = __offset + _offsets[zero + index * stride];
    if (value.real == 0 && value.imaginary == 0) {
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
      _elements[i] = new Complex(re, im);
    }
  }

  ComplexVector _viewSelectionLike(Int32List offsets) {
    return new SelectedSparseComplexVector.withOffsets(this._elements, offsets);
  }

  /// This method is not supported.
  DoubleVector imaginary() {
    throw new UnsupportedError("This method is not supported.");
  }

  /// This method is not supported.
  DoubleVector real() {
    throw new UnsupportedError("This method is not supported.");
  }

  Object clone() {
    return new SelectedSparseComplexVector(
        size, _elements, zero, stride, _offsets, __offset);
  }
}
