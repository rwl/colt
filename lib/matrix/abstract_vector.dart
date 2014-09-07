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
 * Abstract base class for 1-d matrices (aka <i>vectors</i>) holding objects or
 * primitive data types such as <code>int</code>, <code>double</code>, etc.
 * First see the <a href="package-summary.html">package summary</a> and javadoc
 * <a href="package-tree.html">tree view</a> to get the broad picture.
 * <p>
 * <b>Note that this implementation is not synchronized.</b>
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 09/24/99
 */
abstract class AbstractVector {

  bool _isNoView = true;

  /** the number of cells this matrix (view) has */
  int _size;

  /** the index of the first element */
  int _zero;

  /**
   * the number of indexes between any two elements, i.e.
   * <tt>index(i+1) - index(i)</tt>.
   */
  int _stride;

  /**
   * Makes this class non instantiable, but still let's others inherit from
   * it.
   */
  AbstractVector();

  /**
   * Returns whether the receiver is a view or not.
   */
  bool get isView {
    return !this._isNoView;
  }

  /**
   * Returns the number of cells.
   */
  //int get length;

  /**
   * Returns the position of the given absolute rank within the (virtual or
   * non-virtual) internal 1-dimensional array. Default implementation.
   * Override, if necessary.
   *
   * @param rank
   *            the absolute rank of the element.
   * @return the position.
   */
  int _offset(int absRank) {
    return absRank;
  }

  /**
   * Returns the absolute rank of the given relative rank.
   *
   * @param rank
   *            the relative rank of the element.
   * @return the absolute rank of the element.
   */
  int _rank(int rank) {
    return _zero + rank * _stride;
  }

  /**
   * Sanity check for operations requiring an index to be within bounds.
   *
   * @throws IndexOutOfBoundsException
   *             if <tt>index < 0 || index >= size()</tt>.
   */
  void _checkIndex(int index) {
    if (index < 0 || index >= _size) {
      throw new RangeError("Attempted to access " + toStringShort() + " at index=$index");
    }
  }

  /**
   * Checks whether indexes are legal and throws an exception, if necessary.
   *
   * @throws IndexOutOfBoundsException
   *             if <tt>! (0 <= indexes[i] < size())</tt> for any
   *             i=0..indexes.length()-1.
   */
  void _checkIndexes(List<int> indexes) {
    for (int i = indexes.length; --i >= 0; ) {
      int index = indexes[i];
      if (index < 0 || index >= _size) _checkIndex(index);
    }
  }

  /**
   * Checks whether the receiver contains the given range and throws an
   * exception, if necessary.
   *
   * @throws IndexOutOfBoundsException
   *             if <tt>index<0 || index+width>size()</tt>.
   */
  void _checkRange(int index, int width) {
    if (index < 0 || index + width > _size) {
      throw new RangeError("index: $index, width: $width, size=$_size");
    }
  }

  /**
   * Sanity check for operations requiring two matrices with the same size.
   *
   * @throws IllegalArgumentException
   *             if <tt>size() != B.size()</tt>.
   */
  void checkSize(AbstractVector B) {
    if (_size != B._size) {
      throw new ArgumentError("Incompatible sizes: " + toStringShort() + " and " + B.toStringShort());
    }
  }

  /**
   * Returns the position of the element with the given relative rank within
   * the (virtual or non-virtual) internal 1-dimensional array. You may want
   * to override this method for performance.
   *
   * @param rank
   *            the rank of the element.
   */
  int index(int rank) {
    return _offset(_rank(rank));
  }

  /**
   * Sets up a matrix with the given parameters.
   *
   * @param size
   *            the number of elements the matrix shall have.
   * @param zero
   *            the index of the first element.
   * @param stride
   *            the number of indexes between any two elements, i.e.
   *            <tt>index(i+1)-index(i)</tt>.
   * @throws IllegalArgumentException
   *             if <tt>size<0</tt>.
   */
  void _setUp(int size, [int zero = 0, int stride = 1]) {
    if (size < 0) throw new ArgumentError("negative size");

    this._size = size;
    this._zero = zero;
    this._stride = stride;
    this._isNoView = true;
  }

  /**
   * Returns the number of cells.
   */
  int get length {
    return _size;
  }

  /**
   * Returns the stride of the given dimension (axis, rank).
   *
   * @dimension the index of the dimension.
   * @return the stride in the given dimension.
   * @throws IllegalArgumentException if <tt>dimension != 0</tt>.
   */
  int stride([int dimension = 0]) {
    if (dimension != 0) {
      throw new ArgumentError("invalid dimension: $dimension used to access" + toStringShort());
    }
    return this._stride;
  }

  /**
   * Returns a string representation of the receiver's shape.
   */
  String toStringShort() {
    return AbstractFormatter.shape(this);
  }

  /**
   * Self modifying version of viewFlip(). What used to be index <tt>0</tt> is
   * now index <tt>size()-1</tt>, ..., what used to be index <tt>size()-1</tt>
   * is now index <tt>0</tt>.
   */
  AbstractVector _vFlip() {
    if (_size > 0) {
      this._zero += (this.length - 1) * this._stride;
      this._stride = -this._stride;
      this._isNoView = false;
    }
    return this;
  }

  /**
   * Self modifying version of viewPart().
   *
   * @throws IndexOutOfBoundsException
   *             if <tt>index<0 || index+width>size()</tt>.
   */
  AbstractVector _vPart(int index, int width) {
    _checkRange(index, width);
    this._zero += this._stride * index;
    this._size = width;
    this._isNoView = false;
    return this;
  }

  /**
   * Self modifying version of viewStrides().
   *
   * @throws IndexOutOfBoundsException
   *             if <tt>stride <= 0</tt>.
   */
  AbstractVector _vStrides(int stride) {
    if (stride <= 0) {
      throw new RangeError("illegal stride: $stride");
    }
    this._stride *= stride;
    if (this._size != 0) {
      this._size = (this._size - 1) ~/ stride + 1;
    }
    this._isNoView = false;
    return this;
  }
}
