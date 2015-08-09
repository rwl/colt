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

/**
 * 1-d matrix holding <tt>int</tt> elements; either a view wrapping another 2-d
 * matrix and therefore delegating calls to it.
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 09/24/99
 */
//class DelegateIntVector extends WrapperIntVector {
//
//  /*
//   * The elements of the matrix.
//   */
//  IntMatrix _content;
//
//  /*
//   * The row this view is bound to.
//   */
//  int _row;
//
//  DelegateIntVector(IntMatrix newContent, int row) : super(null) {
//    if (row < 0 || row >= newContent.rows()) {
//      throw new ArgumentError();
//    }
//    _setUp(newContent.columns);
//    this._row = row;
//    this._content = newContent;
//  }
//
//  /**
//   * Returns the matrix cell value at coordinate <tt>index</tt>.
//   *
//   * <p>
//   * Provided with invalid parameters this method may return invalid objects
//   * without throwing any exception. <b>You should only use this method when
//   * you are absolutely sure that the coordinate is within bounds.</b>
//   * Precondition (unchecked): <tt>index&lt;0 || index&gt;=size()</tt>.
//   *
//   * @param index
//   *            the index of the cell.
//   * @return the value of the specified cell.
//   */
//  int get(int index) {
//    return _content.get(_row, index);
//  }
//
//  /**
//   * Construct and returns a new empty matrix <i>of the same dynamic type</i>
//   * as the receiver, having the specified size. For example, if the receiver
//   * is an instance of type <tt>DenseIntMatrix1D</tt> the new matrix must also
//   * be of type <tt>DenseIntMatrix1D</tt>, if the receiver is an instance of
//   * type <tt>SparseIntMatrix1D</tt> the new matrix must also be of type
//   * <tt>SparseIntMatrix1D</tt>, etc. In general, the new matrix should have
//   * internal parametrization as similar as possible.
//   *
//   * @param size
//   *            the number of cell the matrix shall have.
//   * @return a new empty matrix of the same dynamic type.
//   */
//  IntVector like(int size) {
//    return _content.like1D(size);
//  }
//
//  /**
//   * Construct and returns a new 2-d matrix <i>of the corresponding dynamic
//   * type</i>, entirelly independent of the receiver. For example, if the
//   * receiver is an instance of type <tt>DenseIntMatrix1D</tt> the new matrix
//   * must be of type <tt>DenseIntMatrix2D</tt>, if the receiver is an instance
//   * of type <tt>SparseIntMatrix1D</tt> the new matrix must be of type
//   * <tt>SparseIntMatrix2D</tt>, etc.
//   *
//   * @param rows
//   *            the number of rows the matrix shall have.
//   * @param columns
//   *            the number of columns the matrix shall have.
//   * @return a new matrix of the corresponding dynamic type.
//   */
//  IntMatrix like2D(int rows, int columns) {
//    return _content.like2D(rows, columns);
//  }
//
//  /**
//   * Sets the matrix cell at coordinate <tt>index</tt> to the specified value.
//   *
//   * <p>
//   * Provided with invalid parameters this method may access illegal indexes
//   * without throwing any exception. <b>You should only use this method when
//   * you are absolutely sure that the coordinate is within bounds.</b>
//   * Precondition (unchecked): <tt>index&lt;0 || index&gt;=size()</tt>.
//   *
//   * @param index
//   *            the index of the cell.
//   * @param value
//   *            the value to be filled into the specified cell.
//   */
//  void set(int index, int value) {
//    _content.set(_row, index, value);
//  }
//}


class DelegateIntVector extends IntVector {
  IntMatrix _content;

  /// The row this view is bound to.
  int _row;

  /// Constructs a matrix view with a given content and row.
  DelegateIntVector(IntMatrix newContent, int row)
      : super(newContent.columns) {
    if (row < 0 || row >= newContent.rows) {
      throw new ArgumentError();
    }
    _row = row;
    _content = newContent;
  }

  int get(int index) => _content.get(_row, index);

  IntVector like1D(int size) => _content.like1D(size);

  IntMatrix like2D(int rows, int columns) {
    return _content.like2D(rows, columns);
  }

  void set(int index, int value) => _content.set(_row, index, value);

  Object get elements => _content.elements;

  IntMatrix reshape(int rows, int columns) {
    throw new ArgumentError("This method is not supported.");
  }

  IntVector _viewSelectionLike(Int32List offsets) {
    throw new Error(); // should never get called
  }

  Object clone() => new DelegateIntVector(_content, _row);
}
