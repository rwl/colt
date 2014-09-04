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
 * 1-d matrix holding <tt>complex</tt> elements; either a view wrapping another
 * 2-d matrix and therefore delegating calls to it.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
class DelegateDComplexMatrix1D extends DComplexMatrix1D {

  /*
   * The elements of the matrix.
   */
  DComplexMatrix2D _content;

  /*
   * The row this view is bound to.
   */
  int _row;

  /**
   * Creates new instance of DelegateDComplexMatrix1D
   *
   * @param newContent
   *            the content
   * @param row
   *            the row this view is bound to
   */
  DelegateDComplexMatrix1D(DComplexMatrix2D newContent, int row) {
    if (row < 0 || row >= newContent.rows) {
      throw new ArgumentError();
    }
    _setUp(newContent.columns);
    this._row = row;
    this._content = newContent;
  }

  Float64List getQuick(int index) {
    return _content.getQuick(_row, index);
  }

  DComplexMatrix1D like1D(int size) {
    return _content.like1D(size);
  }

  DComplexMatrix2D like2D(int rows, int columns) {
    return _content.like2D(rows, columns);
  }

  void setQuick(int index, Float64List value) {
    _content.setQuick(_row, index, value);
  }

  void setPartsQuick(int index, double re, double im) {
    _content.setPartsQuick(_row, index, re, im);
  }

  Object elements() {
    return _content.elements();
  }

  DComplexMatrix2D reshape(int rows, int columns) {
    throw new ArgumentError("This method is not supported.");
  }

  /*DComplexMatrix3D reshape3D(int slices, int rows, int columns) {
    throw new ArgumentError("This method is not supported.");
  }*/

  DComplexMatrix1D _viewSelectionLike(Int32List offsets) {
    throw new ArgumentError("This method is not supported.");
  }

  DoubleMatrix1D getImaginaryPart() {
    return _content.viewRow(_row).getImaginaryPart();
  }

  DoubleMatrix1D getRealPart() {
    return _content.viewRow(_row).getRealPart();
  }
  
  Object clone() {
    return new DelegateDComplexMatrix1D(_content, _row);
  }

}

/**
 * 2-d matrix holding <tt>complex</tt> elements; either a view wrapping another
 * 3-d matrix and therefore delegating calls to it.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
/*class DelegateDComplexMatrix2D extends DComplexMatrix2D {

  /*
     * The elements of the matrix.
     */
  DComplexMatrix3D _content;

  /*
   * The index this view is bound to.
   */
  int _index;

  int _axis; //0-2

  /**
   * Creates new instance of DelegateDComplexMatrix2D
   *
   * @param newContent
   *            the content
   * @param axis
   *            the axis this view is bound to
   * @param index
   *            the index this view is bound to
   */
  DelegateDComplexMatrix2D(DComplexMatrix3D newContent, int axis, int index) {
    switch (axis) {
      case 0:
        if (index < 0 || index >= newContent.slices()) throw new ArgumentError();
        _setUp(newContent.rows(), newContent.columns());
        break;
      case 1:
        if (index < 0 || index >= newContent.rows()) throw new ArgumentError();
        _setUp(newContent.slices(), newContent.columns());
        break;
      case 2:
        if (index < 0 || index >= newContent.columns()) throw new ArgumentError();
        _setUp(newContent.slices(), newContent.rows());
        break;
      default:
        throw new ArgumentError();
    }
    this._axis = axis;
    this._index = index;
    this._content = newContent;
  }

  Float64List getQuick(int row, int column) {
    switch (_axis) {
      case 0:
        return _content.getQuick(_index, row, column);
      case 1:
        return _content.getQuick(row, _index, column);
      case 2:
        return _content.getQuick(row, column, _index);
      default:
        throw new ArgumentError();
    }
  }

  DComplexMatrix2D like2D(int rows, int columns) {
    return _content.like2D(rows, columns);
  }

  void setQuick(int row, int column, Float64List value) {
    switch (_axis) {
      case 0:
        _content.setQuick(_index, row, column, value);
        break;
      case 1:
        _content.setQuick(row, _index, column, value);
        break;
      case 2:
        _content.setQuick(row, column, _index, value);
        break;
      default:
        throw new ArgumentError();
    }
  }

  void setPartsQuick(int row, int column, double re, double im) {
    switch (_axis) {
      case 0:
        _content.setQuick(_index, row, column, re, im);
        break;
      case 1:
        _content.setQuick(row, _index, column, re, im);
        break;
      case 2:
        _content.setQuick(row, column, _index, re, im);
        break;
      default:
        throw new ArgumentError();
    }
  }

  DComplexMatrix1D viewColumn(int column) {
    _checkColumn(column);
    return new WrapperDComplexMatrix2D(this).viewColumn(column);
  }

  Object elements() {
    return _content.elements();
  }

  DComplexMatrix2D _viewSelectionLike(Int32List rowOffsets, Int32List columnOffsets) {
    throw new ArgumentError("This method is not supported.");
  }

  DComplexMatrix1D like1D(int size) {
    throw new ArgumentError("This method is not supported.");
  }

  DComplexMatrix1D _like1D(int size, int zero, int stride) {
    throw new ArgumentError("This method is not supported.");
  }

  DComplexMatrix1D vectorize() {
    DComplexMatrix1D v = new DenseDComplexMatrix1D(_rows * _columns);
    int idx = 0;
    for (int c = 0; c < _columns; c++) {
      for (int r = 0; r < _rows; r++) {
        v.setQuick(idx++, getQuick(r, c));
      }
    }
    return v;
  }

  DoubleMatrix2D getImaginaryPart() {
    switch (_axis) {
      case 0:
        return _content.viewSlice(_index).getImaginaryPart();
      case 1:
        return _content.viewRow(_index).getImaginaryPart();
      case 2:
        return _content.viewColumn(_index).getImaginaryPart();
      default:
        throw new ArgumentError();
    }
  }

  DoubleMatrix2D getRealPart() {
    switch (_axis) {
      case 0:
        return _content.viewSlice(_index).getRealPart();
      case 1:
        return _content.viewRow(_index).getRealPart();
      case 2:
        return _content.viewColumn(_index).getRealPart();
      default:
        throw new ArgumentError();
    }
  }

}*/
