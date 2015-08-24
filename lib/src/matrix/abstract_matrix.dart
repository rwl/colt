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

/// Abstract base class for 2-d matrices holding objects or primitive data
/// types such as [int], [double], etc.
abstract class AbstractMatrix {
  bool _isNoView = true;

  /// The number of colums and rows this matrix (view) has.
  int _columns, _rows;

  /// The number of elements between two rows, i.e.
  /// `index(i+1,j,k) - index(i,j,k)`.
  int _rowStride;

  /// The number of elements between two columns, i.e.
  /// `index(i,j+1,k) - index(i,j,k)`.
  int _columnStride;

  /// The index of the first element.
  int _rowZero, _columnZero;

  AbstractMatrix(int rows, int columns, [int rowZero = 0, int columnZero = 0,
      int rowStride = null, int columnStride = 1, bool isNoView = true]) {
    //void _setUp(int rows, int columns, int rowZero, int columnZero, int rowStride, int columnStride) {
    if (rowStride == null) {
      rowStride = columns;
    }
    if (rows < 0 || columns < 0) {
      throw new ArgumentError("negative size");
    }
    _rows = rows;
    _columns = columns;

    _rowZero = rowZero;
    _columnZero = columnZero;

    _rowStride = rowStride;
    _columnStride = columnStride;

    _isNoView = isNoView;
    if (columns * rows > MAX_INT) {
      throw new ArgumentError("matrix too large");
    }
  }

  /// Returns whether the receiver is a view or not.
  bool get isView => !_isNoView;

  /// Returns the position of the given absolute rank within the (virtual or
  /// non-virtual) internal 1-dimensional array. Override, if necessary.
  int columnOffset(int absRank) => absRank;

  /// Returns the absolute rank of the given relative rank.
  int columnRank(int rank) => _columnZero + rank * _columnStride;

  /// Returns the position of the given absolute rank within the (virtual or
  /// non-virtual) internal 1-dimensional array. Override, if necessary.
  int rowOffset(int absRank) => absRank;

  /// Returns the absolute rank of the given relative rank.
  int rowRank(int rank) => _rowZero + rank * _rowStride;

  /// The number of columns.
  int get columns => _columns;

  /// The column stride.
  int get columnStride => _columnStride;

  /// Returns the position of the given coordinate within the (virtual or
  /// non-virtual) internal 1-dimensional array.
  int index(int row, int column) {
    return rowOffset(rowRank(row)) + columnOffset(columnRank(column));
  }

  /// The number of rows.
  int get rows => _rows;

  /// The row stride.
  int get rowStride => _rowStride;

  /// The number of cells which is `rows*columns`.
  int get size => _rows * _columns;

  int get rowZero => _rowZero;

  int get columnZero => _columnZero;

  /// A string representation of the receiver's shape.
  String toStringShort() {
    return AbstractFormatter.shapeMatrix(this);
  }
}

/// Checks whether the receiver contains the given box and throws an
/// exception, if necessary.
void checkBox(AbstractMatrix m, int row, int column, int height, int width) {
  if (column < 0 ||
      width < 0 ||
      column + width > m._columns ||
      row < 0 ||
      height < 0 ||
      row + height > m._rows) {
    throw new RangeError(m.toStringShort() +
        ", column:$column, row:$row ,width:$width, height:$height");
  }
}

/// Sanity check for operations requiring a column index to be within bounds.
void checkColumn(AbstractMatrix m, int column) {
  if (column < 0 || column >= m._columns) {
    throw new RangeError(
        "Attempted to access " + m.toStringShort() + " at column=$column");
  }
}

/// Checks whether indexes are legal and throws an exception, if necessary.
void checkColumnIndexes(AbstractMatrix m, List<int> indexes) {
  for (int index in indexes) {
    checkColumn(m, index);
  }
}

/// Sanity check for operations requiring a row index to be within bounds.
void checkRow(AbstractMatrix m, int row) {
  if (row < 0 || row >= m._rows) {
    throw new RangeError(
        "Attempted to access " + m.toStringShort() + " at row=$row");
  }
}

/// Checks whether indexes are legal and throws an exception, if necessary.
void checkRowIndexes(AbstractMatrix m, List<int> indexes) {
  for (int index in indexes) {
    checkRow(m, index);
  }
}

/// Sanity check for operations requiring matrices with the same number of
/// columns and rows.
void checkShape(AbstractMatrix A, AbstractMatrix B, [AbstractMatrix C = null]) {
  if (C == null) {
    if (A._columns != B._columns || A._rows != B._rows) {
      throw new ArgumentError("Incompatible dimensions: " +
          A.toStringShort() +
          " and " +
          B.toStringShort());
    }
  } else {
    if (A._columns != B._columns ||
        A._rows != B._rows ||
        A._columns != C._columns ||
        A._rows != C._rows) {
      throw new ArgumentError("Incompatible dimensions: " +
          A.toStringShort() +
          ", " +
          B.toStringShort() +
          ", " +
          C.toStringShort());
    }
  }
}

/// Self modifying version of viewColumnFlip().
void vColumnFlip(AbstractMatrix m) {
  if (m._columns > 0) {
    m._columnZero += (m._columns - 1) * m._columnStride;
    m._columnStride = -m._columnStride;
    m._isNoView = false;
  }
}

/// Self modifying version of viewDice().
void vDice(AbstractMatrix m) {
  // swap;
  int tmp = m._rows;
  m._rows = m._columns;
  m._columns = tmp;
  tmp = m._rowStride;
  m._rowStride = m._columnStride;
  m._columnStride = tmp;
  tmp = m._rowZero;
  m._rowZero = m._columnZero;
  m._columnZero = tmp;

  // flips stay unaffected

  m._isNoView = false;
}

/// Self modifying version of viewPart().
void vBox(AbstractMatrix m, int row, int column, int height, int width) {
  checkBox(m, row, column, height, width);
  m._rowZero += m._rowStride * row;
  m._columnZero += m._columnStride * column;
  m._rows = height;
  m._columns = width;
  m._isNoView = false;
}

/// Self modifying version of viewRowFlip().
void vRowFlip(AbstractMatrix m) {
  if (m._rows > 0) {
    m._rowZero += (m._rows - 1) * m._rowStride;
    m._rowStride = -m._rowStride;
    m._isNoView = false;
  }
}

/// Self modifying version of viewStrides().
void vStrides(AbstractMatrix m, int rowStride, int columnStride) {
  if (rowStride <= 0 || columnStride <= 0) {
    throw new RangeError("illegal strides: $rowStride, $columnStride");
  }
  m._rowStride *= rowStride;
  m._columnStride *= columnStride;
  if (m._rows != 0) {
    m._rows = (m._rows - 1) ~/ rowStride + 1;
  }
  if (m._columns != 0) {
    m._columns = (m._columns - 1) ~/ columnStride + 1;
  }
  m._isNoView = false;
}

setIsNoView(AbstractMatrix m, bool isNoView) => m._isNoView = isNoView;

setRows(AbstractMatrix m, int rows) => m._rows = rows;

setColumns(AbstractMatrix m, int columns) => m._columns = columns;
