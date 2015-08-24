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
part of cern.colt.matrix.internal;

/// Abstract base class for 1-d matrices (aka vectors) holding objects or
/// primitive data types such as [int], [double], etc.
abstract class AbstractVector<E> extends ListBase<E> {
  bool _isNoView = true;

  int _size;

  int _zero;

  int _stride;

  /// Returns whether the receiver is a view or not.
  bool get isView => !_isNoView;

  /// Returns the position of the element with the given relative rank within
  /// the (virtual or non-virtual) internal 1-dimensional array. Override this
  /// method for performance.
  int index(int relRank) => offset(rank(relRank));

  /// Returns the position of the given absolute [rank] within the (virtual or
  /// non-virtual) internal 1-dimensional array. Default implementation.
  /// Override, if necessary.
  int offset(int rank) => rank;

  /// Returns the absolute rank of the given relative [rank].
  int rank(int rank) => _zero + rank * _stride;

  /// Sets up a matrix with the given parameters.
//  void _setUp(int size, [int zero = 0, int stride = 1]) {
  AbstractVector(int size,
      [int zero = 0, int stride = 1, bool isNoView = true]) {
    if (size < 0) {
      throw new ArgumentError("negative size");
    }

    _size = size;
    _zero = zero;
    _stride = stride;
    _isNoView = isNoView;
  }

  /// Returns the number of cells.
  int get size => _size;

  /// Throws an [UnsupportedError].

  /// The number of indexes between any two elements.
  int get stride => _stride;

  /// The index of the first element.
  int get zero => _zero;

  /// Returns a string representation of the receiver's shape.
  String toStringShort() {
    return AbstractFormatter.shapeVector(this);
  }

  /// Throws an [UnsupportedError].
  void set length(int _) {
    throw new UnsupportedError('fixed size');
  }

  /// Synonym for [size] ([List] API).
  int get length => size;
}

/// Sanity check for operations requiring an index to be within bounds.
void checkIndex(AbstractVector vector, int index) {
  if (index < 0 || index >= vector._size) {
    throw new RangeError(
        "Attempted to access " + vector.toStringShort() + " at index=$index");
  }
}

/// Checks whether indexes are legal.
void checkIndexes(AbstractVector vector, List<int> indexes) {
  for (int index in indexes) {
    if (index < 0 || index >= vector._size) {
      checkIndex(vector, index);
    }
  }
}

/// Checks whether the receiver contains the given range.
void checkRange(AbstractVector vector, int index, int width) {
  if (index < 0 || index + width > vector._size) {
    throw new RangeError("index: $index, width: $width, size=${vector._size}");
  }
}

/// Sanity check for operations requiring two matrices with the same size.
void checkSize(AbstractVector A, AbstractVector B) {
  if (A._size != B._size) {
    throw new ArgumentError("Incompatible sizes: " +
        A.toStringShort() +
        " and " +
        B.toStringShort());
  }
}

/// Self modifying version of viewFlip(). What used to be index `0` is
/// now index `size()-1`, ..., what used to be index `size()-1` is now
/// index `0`.
void vFlip(AbstractVector vector) {
  if (vector._size > 0) {
    vector._zero += (vector.size - 1) * vector._stride;
    vector._stride = -vector._stride;
    vector._isNoView = false;
  }
}

/// Self modifying version of [viewPart].
void vPart(AbstractVector vector, int index, int width) {
  checkRange(vector, index, width);
  vector._zero += vector._stride * index;
  vector._size = width;
  vector._isNoView = false;
}

/// Self modifying version of [viewStrides].
void vStride(AbstractVector vector, int stride) {
  if (stride <= 0) {
    throw new RangeError("illegal stride: $stride");
  }
  vector._stride *= stride;
  if (vector._size != 0) {
    vector._size = (vector._size - 1) ~/ stride + 1;
  }
  vector._isNoView = false;
}
