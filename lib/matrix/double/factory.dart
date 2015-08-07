library cern.colt.matrix.double.factory;

import 'dart:math' as Math;
import 'dart:typed_data';
import 'package:colt/function/double.dart' as func;
import 'matrix.dart';

typedef AbstractDoubleVector VectorFn(int);
typedef AbstractDoubleMatrix MatrixFn(int rows, int columns);

/// `C = A||B`; Constructs a new matrix which is the concatenation of two
/// other matrices.
AbstractDoubleVector append(AbstractDoubleVector A, AbstractDoubleVector B, VectorFn make) {
  // concatenate
  AbstractDoubleVector matrix = make(A.size + B.size);
  matrix.part(0, A.size).copyFrom(A);
  matrix.part(A.size, B.size).copyFrom(B);
  return matrix;
}

/// Constructs a matrix with cells having ascending values.
AbstractDoubleVector ascending(int size, VectorFn make) {
  return descending(size, make)..apply(func.chainGH(func.neg, func.subtract(size)));
}

/// Constructs a matrix with cells having descending values.
AbstractDoubleVector descending(int size, VectorFn make) {
  AbstractDoubleVector matrix = make(size);
  double v = 0.0;
  for (int i = size; --i >= 0; ) {
    matrix.set(i, v);
    v += 1;
  }
  return matrix;
}

/// Constructs a matrix which is the concatenation of all given parts. Cells
/// are copied.
AbstractDoubleVector concat(List<AbstractDoubleVector> parts, VectorFn make) {
  if (parts.length == 0) {
    return make(0);
  }

  int size = 0;
  for (int i = 0; i < parts.length; i++) {
    size += parts[i].size;
  }

  AbstractDoubleVector vector = make(size);
  size = 0;
  for (int i = 0; i < parts.length; i++) {
    vector.part(size, parts[i].size).copyFrom(parts[i]);
    size += parts[i].size;
  }

  return vector;
}

/// Constructs a matrix with the given shape, each cell initialized with the
/// given value.
AbstractDoubleVector fill(int size, double initialValue, VectorFn make) {
  return make(size)..fill(initialValue);
}

/// Constructs a matrix with uniformly distributed values in `(0,1)`
/// (exclusive).
AbstractDoubleVector random(int size, VectorFn make) {
  return make(size)..apply(func.random());
}

/// `C = A||A||..||A`; Constructs a new matrix which is concatenated
/// [repeat} times. Example:
///
///     0 1
///     repeat(3) -->
///     0 1 0 1 0 1
AbstractDoubleVector repeat(AbstractDoubleVector A, int repeat, VectorFn make) {
  int size = A.size;
  AbstractDoubleVector matrix = make(repeat * size);
  for (int i = repeat; --i >= 0; ) {
    matrix.part(size * i, size).copyFrom(A);
  }
  return matrix;
}

/// Checks whether the given array is rectangular, that is, whether all rows
/// have the same number of columns.
/*void _checkRectangularShape(List<Float64List> array) {
  int columns = -1;
  for (int row = array.length; --row >= 0; ) {
    if (array[row] != null) {
      if (columns == -1) {
        columns = array[row].length;
      }
      if (array[row].length != columns) {
        throw new ArgumentError("All rows of array must have same number of columns.");
      }
    }
  }
}*/

/// Checks whether the given array is rectangular, that is, whether all rows
/// have the same number of columns.
void _checkRectangularShapeParts(List<List<AbstractDoubleMatrix>> array) {
  int columns = -1;
  for (int row = array.length; --row >= 0; ) {
    if (array[row] != null) {
      if (columns == -1) {
        columns = array[row].length;
      }
      if (array[row].length != columns) {
        throw new ArgumentError("All rows of array must have same number of columns.");
      }
    }
  }
}

/// `C = A||b`; Constructs a new matrix which is the column-wise concatenation
/// of two other matrices.
///
///     0 1 2
///     3 4 5
///     appendColumn
///     6
///     8
///     -->
///     0 1 2 6
///     3 4 5 8
AbstractDoubleMatrix appendColumn(AbstractDoubleMatrix A, AbstractDoubleVector b, MatrixFn make) {
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
  AbstractDoubleMatrix matrix = make(r, ac + bc);
  matrix.part(0, 0, r, ac).copyFrom(A);
  matrix.column(ac).copyFrom(b);
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
AbstractDoubleMatrix appendColumns(AbstractDoubleMatrix A, AbstractDoubleMatrix B, MatrixFn make) {
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
  AbstractDoubleMatrix matrix = make(r, ac + bc);
  matrix.part(0, 0, r, ac).copyFrom(A);
  matrix.part(0, ac, r, bc).copyFrom(B);
  return matrix;
}

/// `C = A||b`; Constructs a new matrix which is the row-wise concatenation of
/// two other matrices.
///
///     0 1
///     2 3
///     4 5
///     appendRow
///     6 7
///     -->
///     0 1
///     2 3
///     4 5
///     6 7
AbstractDoubleMatrix appendRow(AbstractDoubleMatrix A, AbstractDoubleVector b, MatrixFn make) {
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
  AbstractDoubleMatrix matrix = make(ar + br, c);
  matrix.part(0, 0, ar, c).copyFrom(A);
  matrix.row(ar).copyFrom(b);
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
AbstractDoubleMatrix appendRows(AbstractDoubleMatrix A, AbstractDoubleMatrix B, MatrixFn make) {
  // force both to have maximal shared number of columns.
  if (B.columns > A.columns) {
    B = B.part(0, 0, B.rows, A.columns);
  } else if (B.columns < A.columns) {
    A = A.part(0, 0, A.rows, B.columns);
  }

  // concatenate
  int ar = A.rows;
  int br = B.rows;
  int c = A.columns;
  AbstractDoubleMatrix matrix = make(ar + br, c);
  matrix.part(0, 0, ar, c).copyFrom(A);
  matrix.part(ar, 0, br, c).copyFrom(B);
  return matrix;
}

/// Constructs a matrix with cells having ascending values.
///
/// Example:
///
///     0 1 2
///     3 4 5
AbstractDoubleMatrix ascendingMatrix(int rows, int columns, MatrixFn make) {
  return descendingMatrix(rows, columns, make)..apply(func.chainGH(func.neg, func.subtract(columns * rows)));
}

/// Constructs a block matrix made from the given parts. The inverse to
/// method [decompose].
///
/// All matrices of a given column within [parts] must have the same
/// number of columns. All matrices of a given row within [parts] must
/// have the same number of rows. Otherwise an [ArgumentError] is thrown.
/// Note that `null` within `parts[row,col]` are an exception to this rule:
/// they are ignored. Cells are copied.
AbstractDoubleMatrix compose(List<List<AbstractDoubleMatrix>> parts, MatrixFn make) {
  _checkRectangularShapeParts(parts);
  int rows = parts.length;
  int columns = 0;
  if (parts.length > 0) {
    columns = parts[0].length;
  }
  AbstractDoubleMatrix empty = make(0, 0);

  if (rows == 0 || columns == 0) {
    return empty;
  }

  // determine maximum column width of each column
  var maxWidths = new Int32List(columns);
  for (int column = columns; --column >= 0; ) {
    int maxWidth = 0;
    for (int row = rows; --row >= 0; ) {
      AbstractDoubleMatrix part = parts[row][column];
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
  for (int row = rows; --row >= 0; ) {
    int maxHeight = 0;
    for (int column = columns; --column >= 0; ) {
      AbstractDoubleMatrix part = parts[row][column];
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
  for (int row = rows; --row >= 0; ) {
    resultRows += maxHeights[row];
  }
  int resultCols = 0;
  for (int column = columns; --column >= 0; ) {
    resultCols += maxWidths[column];
  }

  AbstractDoubleMatrix matrix = make(resultRows, resultCols);

  // copy
  int r = 0;
  for (int row = 0; row < rows; row++) {
    int c = 0;
    for (int column = 0; column < columns; column++) {
      AbstractDoubleMatrix part = parts[row][column];
      if (part != null) {
        matrix.part(r, c, part.rows, part.columns).copyFrom(part);
      }
      c += maxWidths[column];
    }
    r += maxHeights[row];
  }

  return matrix;
}

/// Constructs a matrix with cells having descending values. Example:
///     5 4 3
///     2 1 0
AbstractDoubleMatrix descendingMatrix(int rows, int columns, MatrixFn make) {
  AbstractDoubleMatrix matrix = make(rows, columns);
  double v = 0.0;
  for (int row = rows; --row >= 0; ) {
    for (int column = columns; --column >= 0; ) {
      matrix.set(row, column, v);
      v += 1;
    }
  }
  return matrix;
}

/// Constructs a new diagonal matrix whose diagonal elements are the elements
/// of [vector]. Cells values are copied. The new matrix is not a view.
///
/// Example:
///
///     5 4 3 -->
///     5 0 0
///     0 4 0
///     0 0 3
AbstractDoubleMatrix diagonalMatrix(AbstractDoubleVector vector, MatrixFn make) {
  int size = vector.size;
  AbstractDoubleMatrix diag = make(size, size);
  for (int i = 0; i < size; i++) {
    diag.set(i, i, vector[i]);
  }
  return diag;
}

/// Constructs a new vector consisting of the diagonal elements of [A].
/// Cells values are copied. The new vector is not a view. Example:
///
///     5 0 0 9
///     0 4 0 9
///     0 0 3 9
///     --> 5 4 3
AbstractDoubleVector diagonal(AbstractDoubleMatrix A, MatrixFn make) {
  int min = Math.min(A.rows, A.columns);
  AbstractDoubleVector diag = make(0, 0).like1D(min);
  for (int i = min; --i >= 0; ) {
    diag.set(i, A.get(i, i));
  }
  return diag;
}

/// Constructs an identity matrix (having ones on the diagonal and zeros
/// elsewhere).
AbstractDoubleMatrix identity(int rowsAndColumns, MatrixFn make) {
  AbstractDoubleMatrix matrix = make(rowsAndColumns, rowsAndColumns);
  for (int i = rowsAndColumns; --i >= 0; ) {
    matrix.set(i, i, 1.0);
  }
  return matrix;
}

/// Construct a matrix from a one-dimensional column-major packed array, ala
/// Fortran. Has the form
/// `matrix.get(row,column) == values[row + column*rows]`. The values
/// are copied.
AbstractDoubleMatrix makeColumn(Float64List values, int rows, MatrixFn make) {
  int columns = (rows != 0 ? values.length / rows : 0);
  if (rows * columns != values.length) {
    throw new ArgumentError("Array length must be a multiple of m.");
  }

  AbstractDoubleMatrix matrix = make(rows, columns);
  for (int row = 0; row < rows; row++) {
    for (int column = 0; column < columns; column++) {
      matrix.set(row, column, values[row + column * rows]);
    }
  }
  return matrix;
}

/// Constructs a matrix with the given shape, each cell initialized with the
/// given value.
AbstractDoubleMatrix fillMatrix(int rows, int columns, double initialValue, MatrixFn make) {
  if (initialValue == 0) {
    return make(rows, columns);
  }
  return make(rows, columns)..fill(initialValue);
}

/// Constructs a matrix with uniformly distributed values in `(0,1)`
/// (exclusive).
AbstractDoubleMatrix randomMatrix(int rows, int columns, MatrixFn make) {
  return make(rows, columns)..apply(func.random());
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
AbstractDoubleMatrix repeatMatrix(AbstractDoubleMatrix A, int rowRepeat, int columnRepeat, MatrixFn make) {
  int r = A.rows;
  int c = A.columns;
  AbstractDoubleMatrix matrix = make(r * rowRepeat, c * columnRepeat);
  for (int i = rowRepeat; --i >= 0; ) {
    for (int j = columnRepeat; --j >= 0; ) {
      matrix.part(r * i, c * j, r, c).copyFrom(A);
    }
  }
  return matrix;
}
