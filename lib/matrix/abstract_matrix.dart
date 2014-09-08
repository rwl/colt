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
 * Abstract base class for 2-d matrices holding objects or primitive data types
 * such as <code>int</code>, <code>double</code>, etc. First see the <a
 * href="package-summary.html">package summary</a> and javadoc <a
 * href="package-tree.html">tree view</a> to get the broad picture.
 * <p>
 * <b>Note that this implementation is not synchronized.</b>
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 09/24/99
 */
abstract class AbstractMatrix {

  bool _isNoView = true;

  /** the number of colums and rows this matrix (view) has */
  int _columns, _rows;

  /**
   * the number of elements between two rows, i.e.
   * <tt>index(i+1,j,k) - index(i,j,k)</tt>.
   */
  int _rowStride;

  /**
   * the number of elements between two columns, i.e.
   * <tt>index(i,j+1,k) - index(i,j,k)</tt>.
   */
  int _columnStride;

  /** the index of the first element */
  int _rowZero, _columnZero;

  AbstractMatrix();

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
  int _columnOffset(int absRank) {
    return absRank;
  }

  /**
   * Returns the absolute rank of the given relative rank.
   *
   * @param rank
   *            the relative rank of the element.
   * @return the absolute rank of the element.
   */
  int _columnRank(int rank) {
    return _columnZero + rank * _columnStride;
  }

  /**
   * Returns the position of the given absolute rank within the (virtual or
   * non-virtual) internal 1-dimensional array. Default implementation.
   * Override, if necessary.
   *
   * @param rank
   *            the absolute rank of the element.
   * @return the position.
   */
  int _rowOffset(int absRank) {
    return absRank;
  }

  /**
   * Returns the absolute rank of the given relative rank.
   *
   * @param rank
   *            the relative rank of the element.
   * @return the absolute rank of the element.
   */
  int _rowRank(int rank) {
    return _rowZero + rank * _rowStride;
  }

  /**
   * Checks whether the receiver contains the given box and throws an
   * exception, if necessary.
   *
   * @throws IndexOutOfBoundsException
   *             if
   *             <tt>column<0 || width<0 || column+width>columns() || row<0 || height<0 || row+height>rows()</tt>
   */
  void _checkBox(int row, int column, int height, int width) {
    if (column < 0 || width < 0 || column + width > _columns || row < 0 || height < 0 || row + height > _rows) throw new RangeError(toStringShort() + ", column:$column, row:$row ,width:$width, height:$height");
  }

  /**
   * Sanity check for operations requiring a column index to be within bounds.
   *
   * @throws IndexOutOfBoundsException
   *             if <tt>column < 0 || column >= columns()</tt>.
   */
  void _checkColumn(int column) {
    if (column < 0 || column >= _columns) throw new RangeError("Attempted to access " + toStringShort() + " at column=$column");
  }

  /**
   * Checks whether indexes are legal and throws an exception, if necessary.
   *
   * @throws IndexOutOfBoundsException
   *             if <tt>! (0 <= indexes[i] < columns())</tt> for any
   *             i=0..indexes.length()-1.
   */
  void _checkColumnIndexes(List<int> indexes) {
    for (int i = indexes.length; --i >= 0; ) {
      int index = indexes[i];
      if (index < 0 || index >= _columns) _checkColumn(index);
    }
  }

  /**
   * Sanity check for operations requiring a row index to be within bounds.
   *
   * @throws IndexOutOfBoundsException
   *             if <tt>row < 0 || row >= rows()</tt>.
   */
  void _checkRow(int row) {
    if (row < 0 || row >= _rows) {
      throw new RangeError("Attempted to access " + toStringShort() + " at row=$row");
    }
  }

  /**
   * Checks whether indexes are legal and throws an exception, if necessary.
   *
   * @throws IndexOutOfBoundsException
   *             if <tt>! (0 <= indexes[i] < rows())</tt> for any
   *             i=0..indexes.length()-1.
   */
  void _checkRowIndexes(List<int> indexes) {
    for (int i = indexes.length; --i >= 0; ) {
      int index = indexes[i];
      if (index < 0 || index >= _rows) _checkRow(index);
    }
  }

  /**
   * Sanity check for operations requiring matrices with the same number of
   * columns and rows.
   *
   * @throws IllegalArgumentException
   *             if
   *             <tt>columns() != B.columns() || rows() != B.rows() || columns() != C.columns() || rows() != C.rows()</tt>
   *             .
   */
  void checkShape(AbstractMatrix B, [AbstractMatrix C = null]) {
    if (C == null) {
      if (_columns != B._columns || _rows != B._rows) {
        throw new ArgumentError("Incompatible dimensions: " + toStringShort() + " and " + B.toStringShort());
      }
    } else {
      if (_columns != B._columns || _rows != B._rows || _columns != C._columns || _rows != C._rows) {
        throw new ArgumentError("Incompatible dimensions: " + toStringShort() + ", " + B.toStringShort() + ", " + C.toStringShort());
      }
    }
  }

  /**
   * Returns the number of columns.
   */
  int get columns {
    return _columns;
  }

  /**
   * Returns the column stride.
   */
  int get columnStride {
    return _columnStride;
  }

  /**
   * Returns the position of the given coordinate within the (virtual or
   * non-virtual) internal 1-dimensional array.
   *
   * @param row
   *            the index of the row-coordinate.
   * @param column
   *            the index of the column-coordinate.
   */
  int index(int row, int column) {
    return _rowOffset(_rowRank(row)) + _columnOffset(_columnRank(column));
  }

  /**
   * Returns the number of rows.
   */
  int get rows {
    return _rows;
  }

  /**
   * Returns the row stride.
   */
  int get rowStride {
    return _rowStride;
  }

  /**
   * Sets up a matrix with a given number of rows and columns and the given
   * strides.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
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
  void _setUp(int rows, int columns, [int rowZero = 0, int columnZero = 0, int rowStride = null, int columnStride = 1]) {
  //void _setUp(int rows, int columns, int rowZero, int columnZero, int rowStride, int columnStride) {
    if (rowStride == null) {
      rowStride = columns;
    }
    if (rows < 0 || columns < 0) {
      throw new ArgumentError("negative size");
    }
    this._rows = rows;
    this._columns = columns;

    this._rowZero = rowZero;
    this._columnZero = columnZero;

    this._rowStride = rowStride;
    this._columnStride = columnStride;

    this._isNoView = true;
    if (columns * rows > MAX_INT) {
      throw new ArgumentError("matrix too large");
    }
  }

  /**
   * Returns the number of cells which is <tt>rows()*columns()</tt>.
   */
  int get length {
    return _rows * _columns;
  }

  /**
   * Returns a string representation of the receiver's shape.
   */
  String toStringShort() {
    return AbstractFormatter.shape2D(this);
  }

  /**
   * Self modifying version of viewColumnFlip().
   */
  void _vColumnFlip() {
    if (_columns > 0) {
      _columnZero += (_columns - 1) * _columnStride;
      _columnStride = -_columnStride;
      this._isNoView = false;
    }
  }

  /**
   * Self modifying version of viewDice().
   */
  void _vDice() {
    int tmp;
    // swap;
    tmp = _rows;
    _rows = _columns;
    _columns = tmp;
    tmp = _rowStride;
    _rowStride = _columnStride;
    _columnStride = tmp;
    tmp = _rowZero;
    _rowZero = _columnZero;
    _columnZero = tmp;

    // flips stay unaffected

    this._isNoView = false;
  }

  /**
   * Self modifying version of viewPart().
   *
   * @throws IndexOutOfBoundsException
   *             if
   *             <tt>column<0 || width<0 || column+width>columns() || row<0 || height<0 || row+height>rows()</tt>
   */
  void _vPart(int row, int column, int height, int width) {
    _checkBox(row, column, height, width);
    this._rowZero += this._rowStride * row;
    this._columnZero += this._columnStride * column;
    this._rows = height;
    this._columns = width;
    this._isNoView = false;
  }

  /**
   * Self modifying version of viewRowFlip().
   */
  void _vRowFlip() {
    if (_rows > 0) {
      _rowZero += (_rows - 1) * _rowStride;
      _rowStride = -_rowStride;
      this._isNoView = false;
    }
  }

  /**
   * Self modifying version of viewStrides().
   *
   * @throws IndexOutOfBoundsException
   *             if <tt>rowStride<=0 || columnStride<=0</tt>.
   */
  void _vStrides(int rowStride, int columnStride) {
    if (rowStride <= 0 || columnStride <= 0) throw new RangeError("illegal strides: $rowStride, $columnStride");
    this._rowStride *= rowStride;
    this._columnStride *= columnStride;
    if (this._rows != 0) this._rows = (this._rows - 1) ~/ rowStride + 1;
    if (this._columns != 0) this._columns = (this._columns - 1) ~/ columnStride + 1;
    this._isNoView = false;
  }
}
