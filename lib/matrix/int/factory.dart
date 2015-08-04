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
library cern.colt.matrix.int.factory;

import 'dart:math' as Math;
import 'dart:typed_data';
import 'package:colt/function/int.dart' as ifunc;
import 'matrix.dart';

typedef AbstractIntVector VectorFn(int);
typedef AbstractIntMatrix MatrixFn(int r, int c);

/// `C = A||B`; Constructs a new matrix which is the concatenation of two
/// other matrices.
AbstractIntVector append(
    AbstractIntVector A, AbstractIntVector B, VectorFn make) {
  // concatenate
  AbstractIntVector matrix = make(A.size + B.size);
  matrix.part(0, A.size).copyFrom(A);
  matrix.part(A.size, B.size).copyFrom(B);
  return matrix;
}

/// Constructs a matrix with cells having ascending values.
AbstractIntVector ascending(int size, VectorFn make) {
  return descending(size, make)
    ..apply(ifunc.chain2(ifunc.neg, ifunc.subtract(size)));
}

/// Constructs a matrix with cells having descending values.
AbstractIntVector descending(int size, VectorFn make) {
  AbstractIntVector matrix = make(size);
  int v = 0;
  for (int i = size; --i >= 0;) {
    matrix.set(i, v++);
  }
  return matrix;
}

/// Constructs a matrix which is the concatenation of all given parts.
AbstractIntVector concat(List<AbstractIntVector> parts, VectorFn make) {
  if (parts.length == 0) {
    return make(0);
  }

  int size = 0;
  for (int i = 0; i < parts.length; i++) {
    size += parts[i].size;
  }

  AbstractIntVector vector = make(size);
  size = 0;
  for (int i = 0; i < parts.length; i++) {
    vector.part(size, parts[i].size).copyFrom(parts[i]);
    size += parts[i].size;
  }

  return vector;
}

/// Constructs a matrix with the given shape, each cell initialized with the
/// given value.
AbstractIntVector fill(int size, int initialValue, VectorFn make) {
  return make(size)..fill(initialValue);
}

/// Constructs a matrix with uniformly distributed values in `(0,1)`
/// (exclusive).
AbstractIntVector random(int size, VectorFn make) {
  return make(size)..apply(ifunc.random());
}

/// `C = A||A||..||A`; Constructs a new matrix which is concatenated
/// `repeat` times. Example:
///
///     0 1
///     repeat(3) -->
///     0 1 0 1 0 1
AbstractIntVector repeat(AbstractIntVector A, int repeat, VectorFn make) {
  int size = A.size;
  AbstractIntVector matrix = make(repeat * size);
  for (int i = repeat; --i >= 0;) {
    matrix.part(size * i, size).copyFrom(A);
  }
  return matrix;
}

/// `C = A||B`; Constructs a new matrix which is the column-wise concatenation
/// of two other matrices.
///
///     0 1 2
///     3 4 5
///     appendColumns
///     6 7
///     8 9
///     -->
///     0 1 2 6 7
///     3 4 5 8 9
AbstractIntMatrix appendColumns(
    AbstractIntMatrix A, AbstractIntMatrix B, MatrixFn make) {
  // force both to have maximal shared number of rows.
  if (B.rows > A.rows) {
    B = B.part(0, 0, A.rows, B.columns);
  } else if (B.rows < A.rows) {
    A = A.part(0, 0, B.rows, A.columns);
  }

  // concatenate
  int ac = A.columns;
  int bc = B.columns;
  int r = A.rows;
  AbstractIntMatrix matrix = make(r, ac + bc);
  matrix.part(0, 0, r, ac).copyFrom(A);
  matrix.part(0, ac, r, bc).copyFrom(B);
  return matrix;
}

AbstractIntMatrix appendColumn(
    AbstractIntMatrix A, AbstractIntVector b, MatrixFn make) {
  // force both to have maximal shared number of rows.
  if (b.size > A.rows) {
    b = b.part(0, A.rows);
  } else if (b.size < A.rows) {
    A = A.part(0, 0, b.size, A.columns);
  }

  // concatenate
  int ac = A.columns;
  int bc = 1;
  int r = A.rows;
  AbstractIntMatrix matrix = make(r, ac + bc);
  matrix.part(0, 0, r, ac).copyFrom(A);
  matrix.column(ac).copyFrom(b);
  return matrix;
}

/// `C = A||B`; Constructs a new matrix which is the row-wise concatenation of
/// two other matrices.
///
///     0 1
///     2 3
///     4 5
///     appendRows
///     6 7
///     8 9
///     -->
///     0 1
///     2 3
///     4 5
///     6 7
///     8 9
AbstractIntMatrix appendRows(
    AbstractIntMatrix A, AbstractIntMatrix B, MatrixFn make) {
  // force both to have maximal shared number of columns.
  if (B.columns > A.columns) B = B.part(0, 0, B.rows, A.columns);
  else if (B.columns < A.columns) A = A.part(0, 0, A.rows, B.columns);

  // concatenate
  int ar = A.rows;
  int br = B.rows;
  int c = A.columns;
  AbstractIntMatrix matrix = make(ar + br, c);
  matrix.part(0, 0, ar, c).copyFrom(A);
  matrix.part(ar, 0, br, c).copyFrom(B);
  return matrix;
}

AbstractIntMatrix appendRow(
    AbstractIntMatrix A, AbstractIntVector b, MatrixFn make) {
  // force both to have maximal shared number of columns.
  if (b.size > A.columns) {
    b = b.part(0, A.columns);
  } else if (b.size < A.columns) {
    A = A.part(0, 0, A.rows, b.size);
  }

  // concatenate
  int ar = A.rows;
  int br = 1;
  int c = A.columns;
  AbstractIntMatrix matrix = make(ar + br, c);
  matrix.part(0, 0, ar, c).copyFrom(A);
  matrix.row(ar).copyFrom(b);
  return matrix;
}

/// Constructs a matrix with cells having ascending values. Example:
///
///     0 1 2
///     3 4 5
AbstractIntMatrix ascendingMatrix(int rows, int columns, MatrixFn make) {
  return descendingMatrix(rows, columns, make)
    ..apply(ifunc.chain2(ifunc.neg, ifunc.subtract(columns * rows)));
}

/// Checks whether the given array is rectangular, that is, whether all rows
/// have the same number of columns.
void _checkRectangularShape(List<Int32List> array) {
  int columns = -1;
  for (int row = array.length; --row >= 0;) {
    if (array[row] != null) {
      if (columns == -1) {
        columns = array[row].length;
      }
      if (array[row].length != columns) {
        throw new ArgumentError(
            "All rows of array must have same number of columns.");
      }
    }
  }
}

/// Checks whether the given array is rectangular, that is, whether all rows
/// have the same number of columns.
void _checkRectShape(List<List<AbstractIntMatrix>> array) {
  int columns = -1;
  for (int row = array.length; --row >= 0;) {
    if (array[row] != null) {
      if (columns == -1) {
        columns = array[row].length;
      }
      if (array[row].length != columns) {
        throw new ArgumentError(
            "All rows of array must have same number of columns.");
      }
    }
  }
}

AbstractIntMatrix reshape(
    AbstractIntVector a, int rows, int columns, MatrixFn make) {
  if (a.size != rows * columns) {
    throw new ArgumentError("a.length != rows*columns");
  }
  AbstractIntMatrix A = make(rows, columns);
  for (int c = 0; c < columns; c++) {
    A.column(c).copyFrom(a.part(c * rows, rows));
  }
  return A;
}

/// Constructs a block matrix made from the given parts.
///
/// All matrices of a given column within [parts] must have the same
/// number of columns. All matrices of a given row within [parts] must
/// have the same number of rows. Otherwise an [ArgumentError] is thrown.
/// Note that `null`s within `parts[row,col]` are an exception to this
/// rule: they are ignored.
AbstractIntMatrix compose(List<List<AbstractIntMatrix>> parts, MatrixFn make) {
  _checkRectShape(parts);
  int rows = parts.length;
  int columns = 0;
  if (parts.length > 0) {
    columns = parts[0].length;
  }
  AbstractIntMatrix empty = make(0, 0);

  if (rows == 0 || columns == 0) {
    return empty;
  }

  // determine maximum column width of each column
  var maxWidths = new Int32List(columns);
  for (int column = columns; --column >= 0;) {
    int maxWidth = 0;
    for (int row = rows; --row >= 0;) {
      AbstractIntMatrix part = parts[row][column];
      if (part != null) {
        int width = part.columns;
        if (maxWidth > 0 && width > 0 && width != maxWidth) {
          throw new ArgumentError("Different number of columns.");
        }
        maxWidth = Math.max(maxWidth, width);
      }
    }
    maxWidths[column] = maxWidth;
  }

  // determine row height of each row
  var maxHeights = new Int32List(rows);
  for (int row = rows; --row >= 0;) {
    int maxHeight = 0;
    for (int column = columns; --column >= 0;) {
      AbstractIntMatrix part = parts[row][column];
      if (part != null) {
        int height = part.rows;
        if (maxHeight > 0 && height > 0 && height != maxHeight) {
          throw new ArgumentError("Different number of rows.");
        }
        maxHeight = Math.max(maxHeight, height);
      }
    }
    maxHeights[row] = maxHeight;
  }

  // shape of result
  int resultRows = 0;
  for (int row = rows; --row >= 0;) {
    resultRows += maxHeights[row];
  }
  int resultCols = 0;
  for (int column = columns; --column >= 0;) {
    resultCols += maxWidths[column];
  }

  AbstractIntMatrix matrix = make(resultRows, resultCols);

  // copy
  int r = 0;
  for (int row = 0; row < rows; row++) {
    int c = 0;
    for (int column = 0; column < columns; column++) {
      AbstractIntMatrix part = parts[row][column];
      if (part != null) {
        matrix.part(r, c, part.rows, part.columns).copyFrom(part);
      }
      c += maxWidths[column];
    }
    r += maxHeights[row];
  }

  return matrix;
}

/// Constructs a diagonal block matrix from the given parts (the direct
/// sum of two matrices). That is the concatenation
///
///     A 0
///     0 B
AbstractIntMatrix composeDiagonal(
    AbstractIntMatrix A, AbstractIntMatrix B, MatrixFn make) {
  int ar = A.rows;
  int ac = A.columns;
  int br = B.rows;
  int bc = B.columns;
  AbstractIntMatrix sum = make(ar + br, ac + bc);
  sum.part(0, 0, ar, ac).copyFrom(A);
  sum.part(ar, ac, br, bc).copyFrom(B);
  return sum;
}

/// Constructs a diagonal block matrix from the given parts. The
/// concatenation has the form
///
///     A 0 0
///     0 B 0
///     0 0 C
AbstractIntMatrix composeDiag(AbstractIntMatrix A, AbstractIntMatrix B,
    AbstractIntMatrix C, MatrixFn make) {
  AbstractIntMatrix diag =
      make(A.rows + B.rows + C.rows, A.columns + B.columns + C.columns);
  diag.part(0, 0, A.rows, A.columns).copyFrom(A);
  diag.part(A.rows, A.columns, B.rows, B.columns).copyFrom(B);
  diag
      .part(A.rows + B.rows, A.columns + B.columns, C.rows, C.columns)
      .copyFrom(C);
  return diag;
}

AbstractIntMatrix composeBidiagonal(
    AbstractIntMatrix A, AbstractIntMatrix B, MatrixFn make) {
  int ar = A.rows;
  int ac = A.columns;
  int br = B.rows;
  int bc = B.columns;
  AbstractIntMatrix sum = make(ar + br - 1, ac + bc);
  sum.part(0, 0, ar, ac).copyFrom(A);
  sum.part(ar - 1, ac, br, bc).copyFrom(B);
  return sum;
}

/// Splits a block matrix into its constituent blocks; Copies blocks of a
/// matrix into the given parts. The inverse to method [compose].
///
/// All matrices of a given column within [parts] must have the same
/// number of columns. All matrices of a given row within [parts] must
/// have the same number of rows. Otherwise an [ArgumentError] is thrown.
/// Note that `null`s within `parts[row,col]` are an exception to this
/// rule: they are ignored.
void decompose(List<List<AbstractIntMatrix>> parts, AbstractIntMatrix matrix) {
  _checkRectShape(parts);
  int rows = parts.length;
  int columns = 0;
  if (parts.length > 0) {
    columns = parts[0].length;
  }
  if (rows == 0 || columns == 0) {
    return;
  }

  // determine maximum column width of each column
  var maxWidths = new Int32List(columns);
  for (int column = columns; --column >= 0;) {
    int maxWidth = 0;
    for (int row = rows; --row >= 0;) {
      AbstractIntMatrix part = parts[row][column];
      if (part != null) {
        int width = part.columns;
        if (maxWidth > 0 && width > 0 && width != maxWidth) {
          throw new ArgumentError("Different number of columns.");
        }
        maxWidth = Math.max(maxWidth, width);
      }
    }
    maxWidths[column] = maxWidth;
  }

  // determine row height of each row
  var maxHeights = new Int32List(rows);
  for (int row = rows; --row >= 0;) {
    int maxHeight = 0;
    for (int column = columns; --column >= 0;) {
      AbstractIntMatrix part = parts[row][column];
      if (part != null) {
        int height = part.rows;
        if (maxHeight > 0 && height > 0 && height != maxHeight) {
          throw new ArgumentError("Different number of rows.");
        }
        maxHeight = Math.max(maxHeight, height);
      }
    }
    maxHeights[row] = maxHeight;
  }

  // shape of result parts
  int resultRows = 0;
  for (int row = rows; --row >= 0;) {
    resultRows += maxHeights[row];
  }
  int resultCols = 0;
  for (int column = columns; --column >= 0;) {
    resultCols += maxWidths[column];
  }

  if (matrix.rows < resultRows || matrix.columns < resultCols) {
    throw new ArgumentError("Parts larger than matrix.");
  }

  // copy
  int r = 0;
  for (int row = 0; row < rows; row++) {
    int c = 0;
    for (int column = 0; column < columns; column++) {
      AbstractIntMatrix part = parts[row][column];
      if (part != null) {
        part.copyFrom(matrix.part(r, c, part.rows, part.columns));
      }
      c += maxWidths[column];
    }
    r += maxHeights[row];
  }
}

/// Constructs a matrix with cells having descending values.
AbstractIntMatrix descendingMatrix(int rows, int columns, MatrixFn make) {
  AbstractIntMatrix matrix = make(rows, columns);
  int v = 0;
  for (int row = rows; --row >= 0;) {
    for (int column = columns; --column >= 0;) {
      matrix.set(row, column, v++);
    }
  }
  return matrix;
}

/// Constructs a new diagonal matrix whose diagonal elements are the elements
/// of [vector]. Example:
///
///     5 4 3 -->
///     5 0 0
///     0 4 0
///     0 0 3
AbstractIntMatrix diagonal(AbstractIntVector vector, MatrixFn make) {
  int size = vector.size;
  AbstractIntMatrix diag = make(size, size);
  for (int i = size; --i >= 0;) {
    diag.set(i, i, vector.get(i));
  }
  return diag;
}

/// Constructs a new vector consisting of the diagonal elements of `A`.
/// Example:
///
///     5 0 0 9
///     0 4 0 9
///     0 0 3 9
///     --> 5 4 3
AbstractIntVector diagonal2D(AbstractIntMatrix A, MatrixFn make) {
  int min = Math.min(A.rows, A.columns);
  AbstractIntVector diag = make(0, 0).like1D(min);
  for (int i = min; --i >= 0;) {
    diag.set(i, A.get(i, i));
  }
  return diag;
}

/// Constructs an identity matrix (having ones on the diagonal and zeros
/// elsewhere).
AbstractIntMatrix identity(int rowsAndColumns, MatrixFn make) {
  AbstractIntMatrix matrix = make(rowsAndColumns, rowsAndColumns);
  for (int i = rowsAndColumns; --i >= 0;) {
    matrix.set(i, i, 1);
  }
  return matrix;
}

/// Construct a matrix from a one-dimensional column-major packed array, ala
/// Fortran. Has the form
/// `matrix.get(row,column) == values[row + column*rows]`.
AbstractIntMatrix columnMajor(Int32List values, int rows, MatrixFn make) {
  int columns = (rows != 0 ? values.length / rows : 0);
  if (rows * columns != values.length) {
    throw new ArgumentError("Array length must be a multiple of m.");
  }

  AbstractIntMatrix matrix = make(rows, columns);
  for (int row = 0; row < rows; row++) {
    for (int column = 0; column < columns; column++) {
      matrix.set(row, column, values[row + column * rows]);
    }
  }
  return matrix;
}

/// Constructs a matrix with the given shape, each cell initialized with the
/// given value.
AbstractIntMatrix fillMatrix(
    int rows, int columns, int initialValue, MatrixFn make) {
  if (initialValue == 0) return make(rows, columns);
  return make(rows, columns)..fill(initialValue);
}

/// Constructs a matrix with uniformly distributed values in `(0,1)`
/// (exclusive).
AbstractIntMatrix randomMatrix(int rows, int columns, MatrixFn make) {
  return make(rows, columns)..apply(ifunc.random());
}

/// `C = A||A||..||A`; Constructs a new matrix which is duplicated both along
/// the row and column dimension. Example:
///
///     0 1
///     2 3
///     repeat(2,3) -->
///     0 1 0 1 0 1
///     2 3 2 3 2 3
///     0 1 0 1 0 1
///     2 3 2 3 2 3
AbstractIntMatrix repeatMatrix(
    AbstractIntMatrix A, int rowRepeat, int columnRepeat, MatrixFn make) {
  int r = A.rows;
  int c = A.columns;
  AbstractIntMatrix matrix = make(r * rowRepeat, c * columnRepeat);
  for (int i = rowRepeat; --i >= 0;) {
    for (int j = columnRepeat; --j >= 0;) {
      matrix.part(r * i, c * j, r, c).copyFrom(A);
    }
  }
  return matrix;
}
