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
 * 1-d matrix holding <tt>double</tt> elements; a view wrapping another 2-d
 * matrix and therefore delegating calls to it.
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 09/24/99
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
class DelegateDoubleMatrix1D extends DoubleMatrix1D {

  /*
   * The elements of the matrix.
   */
  DoubleMatrix2D _content;

  /*
   * The row this view is bound to.
   */
  int _row;

  /**
   * Constructs a matrix view with a given content and row
   *
   * @param newContent
   *            the content of this view
   * @param row
   *            the row this view is bound to
   */
  DelegateDoubleMatrix1D(DoubleMatrix2D newContent, int row) {
    if (row < 0 || row >= newContent.rows()) {
      throw new ArgumentError();
    }
    _setUp(newContent.columns());
    this._row = row;
    this._content = newContent;
  }

  double getQuick(int index) {
    return _content.getQuick(_row, index);
  }

  DoubleMatrix1D like1D(int size) {
    return _content.like1D(size);
  }

  DoubleMatrix2D like2D(int rows, int columns) {
    return _content.like2D(rows, columns);
  }

  void setQuick(int index, double value) {
    _content.setQuick(_row, index, value);
  }

  Object elements() {
    return _content.elements();
  }

  DoubleMatrix2D reshape(int rows, int columns) {
    throw new ArgumentError("This method is not supported.");
  }

  /*DoubleMatrix3D reshape(int slices, int rows, int columns) {
    throw new IllegalArgumentException("This method is not supported.");
  }*/

  DoubleMatrix1D _viewSelectionLike(Int32List offsets) {
    throw new Error(); // should never get called
  }

}
