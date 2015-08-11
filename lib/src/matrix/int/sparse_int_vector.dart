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

/// Sparse hashed 1-d matrix (aka vector) holding [int] elements in a [Map].
///
/// Cells that:
/// - are never set to non-zero values do not use any memory.
/// - switch from zero to non-zero state do use memory.
/// - switch back from non-zero to zero state also do use memory.
///
/// This class offers expected time complexity `O(1)` (i.e. constant time)
/// for the basic operations `get`, `set`, and `size`. As such this sparse
/// class is expected to have no worse time complexity than its dense
/// counterpart [DenseIntVector]. However, constant factors are considerably
/// larger.
class SparseIntVector extends IntVector {
  Map<int, int> _elements;

  /// Constructs a matrix with a copy of the given values.
  factory SparseIntVector.fromList(Int32List values) {
    return new SparseIntVector(values.length)..setAll(values);
  }

  /// Constructs a matrix with a given number of cells. All entries are
  /// initially `0`.
  factory SparseIntVector(int size) {
    return new SparseIntVector._internal(size, {}, 0, 1);
  }

  /// Constructs a matrix view with a given number of parameters.
  SparseIntVector._internal(
      int size, Map<int, int> elements, int offset, int stride)
      : super._(size, offset, stride, false) {
    _elements = elements;
  }

  static SparseIntVector create(int size) {
    return new SparseIntVector(size);
  }

  void fill(int value) {
    if (!isView && value == 0) {
      _elements.clear();
    } else {
      super.fill(value);
    }
  }

  int get cardinality {
    if (!isView) {
      return _elements.length;
    } else {
      return super.cardinality;
    }
  }

  Object get elements => _elements;

  int get(int index) {
    var i = zero + index * stride;
    if (_elements.containsKey(i)) {
      return _elements[i];
    } else {
      return 0;
    }
  }

  int index(int rank) => zero + rank * stride;

  IntVector like1D(int size) => new SparseIntVector(size);

  IntMatrix like2D(int rows, int columns) {
    return new SparseIntMatrix(rows, columns);
  }

  IntMatrix reshape(int rows, int columns) {
    if (rows * columns != size) {
      throw new ArgumentError("rows*columns != size");
    }
    var M = new SparseIntMatrix(rows, columns);
    int idx = 0;
    for (int c = 0; c < columns; c++) {
      for (int r = 0; r < rows; r++) {
        int elem = get(idx++);
        if (elem != 0) {
          M.set(r, c, elem);
        }
      }
    }
    return M;
  }

  void set(int index, int value) {
    int i = zero + index * stride;
    if (value == 0) {
      _elements.remove(i);
    } else {
      _elements[i] = value;
    }
  }

  String toString() {
    var buf = new StringBuffer();
    buf.write("1 x $size sparse matrix, nnz = $cardinality\n");
    for (int i = 0; i < size; i++) {
      int elem = get(i);
      if (elem != 0) {
        buf.write('($i)    $elem\n');
      }
    }
    return buf.toString();
  }

  bool _haveSharedCellsRaw(IntVector other) {
    if (other is SelectedSparseIntVector) {
      return _elements == other._elements;
    } else if (other is SparseIntVector) {
      return _elements == other._elements;
    }
    return false;
  }

  IntVector _viewSelectionLike(Int32List offsets) {
    return new SelectedSparseIntVector(_elements, offsets);
  }

  Object clone() {
    return new SparseIntVector._internal(size, _elements, zero, stride);
  }
}

/// Selection view on sparse 1-d matrices holding [int] elements.
///
/// Objects of this class are typically constructed via `select` methods on
/// some source matrix. The interface introduced in abstract super classes
/// defines everything a user can do. From a user point of view there is
/// nothing special about this class; it presents the same functionality with the
/// same signatures and semantics as its abstract superclass(es) while
/// introducing no additional functionality. Thus, this class need not be visible
/// to users.
class SelectedSparseIntVector extends IntVector {
  Map<int, int> _elements;

  /// The offsets of visible indexes of this matrix.
  Int32List _offsets;

  int __offset;

  factory SelectedSparseIntVector(Map<int, int> elements, Int32List offsets) {
    return new SelectedSparseIntVector._internal(
        offsets.length, elements, 0, 1, offsets, 0);
  }

  SelectedSparseIntVector._internal(int size, Map<int, int> elements, int zero,
      int stride, Int32List offsets, int offset)
      : super._(size, zero, stride, false) {
    _elements = elements;
    _offsets = offsets;
    __offset = offset;
  }

  /// This method is not supported.
  Object get elements {
    throw new ArgumentError("This method is not supported.");
  }

  int get(int index) => _elements[__offset + _offsets[zero + index * stride]];

  int index(int rank) => __offset + _offsets[zero + rank * stride];

  IntVector like1D(int size) => new SparseIntVector(size);

  IntMatrix like2D(int rows, int columns) {
    return new SparseIntMatrix(rows, columns);
  }

  // This method is not supported.
  IntMatrix reshape(int rows, int columns) {
    throw new ArgumentError("This method is not supported.");
  }

  void set(int index, int value) {
    int i = __offset + _offsets[zero + index * stride];
    if (value == 0) {
      _elements.remove(i);
    } else {
      _elements[i] = value;
    }
  }

  int _offset(int absRank) => _offsets[absRank];

  bool _haveSharedCellsRaw(IntVector other) {
    if (other is SelectedSparseIntVector) {
      return _elements == other._elements;
    } else if (other is SparseIntVector) {
      return _elements == other._elements;
    }
    return false;
  }

  IntVector _viewSelectionLike(Int32List offsets) {
    return new SelectedSparseIntVector(_elements, offsets);
  }

  Object clone() {
    return new SelectedSparseIntVector._internal(
        size, _elements, zero, stride, _offsets, __offset);
  }
}
