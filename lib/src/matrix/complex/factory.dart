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
library cern.colt.matrix.complex.factory;

import 'dart:math' as Math;
import 'dart:typed_data';
import 'package:colt/function/complex.dart' as cfunc;
import 'matrix.dart';
import 'package:complex/complex.dart';

typedef ComplexVector VectorFn(int);
typedef ComplexMatrix MatrixFn(int r, int c);

/// `C = A||B`; Constructs a new matrix which is the concatenation of two other
/// matrices.
ComplexVector append(ComplexVector A, ComplexVector B, VectorFn make) {
  ComplexVector matrix = make((A.size + B.size));
  matrix.part(0, A.size).copyFrom(A);
  matrix.part(A.size, B.size).copyFrom(B);
  return matrix;
}

/// Constructs a matrix which is the concatenation of all given parts. Cells
/// are copied.
ComplexVector concat(List<ComplexVector> parts, VectorFn make) {
  if (parts.length == 0) {
    return make(0);
  }

  int size = 0;
  for (int i = 0; i < parts.length; i++) {
    size += parts[i].size;
  }

  ComplexVector vector = make(size);
  size = 0;
  for (int i = 0; i < parts.length; i++) {
    vector.part(size, parts[i].size).copyFrom(parts[i]);
    size += parts[i].size;
  }

  return vector;
}

/// Constructs a matrix with the given shape, each cell initialized with the
/// given value.
ComplexVector fill(int size, Complex initialValue, VectorFn make) {
  return make(size)..fill(initialValue.real, initialValue.imaginary);
}

/// Constructs a matrix with uniformly distributed values in `(0,1)`
/// (exclusive).
ComplexVector random(int size, VectorFn make) {
  return make(size)..apply(cfunc.random);
}

/// `C = A||A||..||A`; Constructs a new matrix which is concatenated
/// `repeat` times.
ComplexVector repeat(ComplexVector A, int repeat, VectorFn make) {
  int size = A.size;
  ComplexVector matrix = make(repeat * size);
  for (int i = 0; i < repeat; i++) {
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
ComplexMatrix appendColumns(ComplexMatrix A, ComplexMatrix B, MatrixFn make) {
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
  ComplexMatrix matrix = make(r, ac + bc);
  matrix.part(0, 0, r, ac).copyFrom(A);
  matrix.part(0, ac, r, bc).copyFrom(B);
  return matrix;
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
ComplexMatrix appendColumn(ComplexMatrix A, ComplexVector b, MatrixFn make) {
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
  ComplexMatrix matrix = make(r, ac + bc);
  matrix.part(0, 0, r, ac).copyFrom(A);
  matrix.column(ac).copyFrom(b);
  return matrix;
}

/// `C = A||B`; Constructs a new matrix which is the row-wise concatenation of
/// two other matrices.
ComplexMatrix appendRows(ComplexMatrix A, ComplexMatrix B, MatrixFn make) {
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
  ComplexMatrix matrix = make(ar + br, c);
  matrix.part(0, 0, ar, c).copyFrom(A);
  matrix.part(ar, 0, br, c).copyFrom(B);
  return matrix;
}

/// `C = A||b`; Constructs a new matrix which is the row-wise concatenation of
/// two other matrices.
ComplexMatrix appendRow(ComplexMatrix A, ComplexVector b, MatrixFn make) {
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
  ComplexMatrix matrix = make(ar + br, c);
  matrix.part(0, 0, ar, c).copyFrom(A);
  matrix.row(ar).copyFrom(b);
  return matrix;
}

/// Checks whether the given array is rectangular, that is, whether all rows
/// have the same number of columns.
/*void _checkRectangularShape(List<Float64List> array) {
  int columns = -1;
  for (int r = 0; r < array.length; r++) {
    if (array[r] != null) {
      if (columns == -1) columns = array[r].length;
      if (array[r].length != columns) {
        throw new ArgumentError("All rows of array must have same number of columns.");
      }
    }
  }
}*/

/// Checks whether the given array is rectangular, that is, whether all rows
/// have the same number of columns.
void _checkRectangularShape2D(List<List<ComplexMatrix>> array) {
  int columns = -1;
  for (int r = 0; r < array.length; r++) {
    if (array[r] != null) {
      if (columns == -1) columns = array[r].length;
      if (array[r].length != columns) {
        throw new ArgumentError("All rows of array must have same number of columns.");
      }
    }
  }
}

/// Constructs a block matrix made from the given parts. The inverse to
/// method [decompose].
///
/// All matrices of a given column within [parts] must have the same
/// number of columns. All matrices of a given row within [parts] must
/// have the same number of rows. Otherwise an [ArgumentError] is thrown.
/// Note that `null`s within `parts[row,col]` are an exception to this rule:
/// they are ignored.
ComplexMatrix compose(List<List<ComplexMatrix>> parts, MatrixFn make) {
  _checkRectangularShape2D(parts);
  int rows = parts.length;
  int columns = 0;
  if (parts.length > 0) {
    columns = parts[0].length;
  }
  ComplexMatrix empty = make(0, 0);

  if (rows == 0 || columns == 0) {
    return empty;
  }

  // determine maximum column width of each column
  var maxWidths = new Int32List(columns);
  for (int c = 0; c < columns; c++) {
    int maxWidth = 0;
    for (int r = 0; r < rows; r++) {
      ComplexMatrix part = parts[r][c];
      if (part != null) {
        int width = part.columns;
        if (maxWidth > 0 && width > 0 && width != maxWidth) {
          throw new ArgumentError("Different number of columns.");
        }
        maxWidth = Math.max(maxWidth, width);
      }
    }
    maxWidths[c] = maxWidth;
  }

  // determine row height of each row
  var maxHeights = new Int32List(rows);
  for (int r = 0; r < rows; r++) {
    int maxHeight = 0;
    for (int c = 0; c < columns; c++) {
      ComplexMatrix part = parts[r][c];
      if (part != null) {
        int height = part.rows;
        if (maxHeight > 0 && height > 0 && height != maxHeight) {
          throw new ArgumentError("Different number of rows.");
        }
        maxHeight = Math.max(maxHeight, height);
      }
    }
    maxHeights[r] = maxHeight;
  }

  // shape of result
  int resultRows = 0;
  for (int r = 0; r < rows; r++) {
    resultRows += maxHeights[r];
  }
  int resultCols = 0;
  for (int c = 0; c < columns; c++) {
    resultCols += maxWidths[c];
  }

  ComplexMatrix matrix = make(resultRows, resultCols);

  // copy
  int idxr = 0;
  for (int r = 0; r < rows; r++) {
    int idxc = 0;
    for (int c = 0; c < columns; c++) {
      ComplexMatrix part = parts[r][c];
      if (part != null) {
        matrix.part(idxr, idxc, part.rows, part.columns).copyFrom(part);
      }
      idxc += maxWidths[c];
    }
    idxr += maxHeights[r];
  }

  return matrix;
}

/// Constructs a diagonal block matrix from the given parts (the direct
/// sum of two matrices). That is the concatenation:
///     A 0
///     0 B
ComplexMatrix composeDiagonal(ComplexMatrix A, ComplexMatrix B, MatrixFn make) {
  int ar = A.rows;
  int ac = A.columns;
  int br = B.rows;
  int bc = B.columns;
  ComplexMatrix sum = make(ar + br, ac + bc);
  sum.part(0, 0, ar, ac).copyFrom(A);
  sum.part(ar, ac, br, bc).copyFrom(B);
  return sum;
}

/// Constructs a diagonal block matrix from the given parts. The
/// concatenation has the form:
///     A 0 0
///     0 B 0
///     0 0 C
ComplexMatrix composeDiagonal3(ComplexMatrix A, ComplexMatrix B, ComplexMatrix C, MatrixFn make) {
  ComplexMatrix diag = make(A.rows + B.rows + C.rows, A.columns + B.columns + C.columns);
  diag.part(0, 0, A.rows, A.columns).copyFrom(A);
  diag.part(A.rows, A.columns, B.rows, B.columns).copyFrom(B);
  diag.part(A.rows + B.rows, A.columns + B.columns, C.rows, C.columns).copyFrom(C);
  return diag;
}

/// Constructs a bidiagonal block matrix from the given parts.
ComplexMatrix composeBidiagonal(ComplexMatrix A, ComplexMatrix B, MatrixFn make) {
  int ar = A.rows;
  int ac = A.columns;
  int br = B.rows;
  int bc = B.columns;
  ComplexMatrix sum = make(ar + br - 1, ac + bc);
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
void decompose(List<List<ComplexMatrix>> parts, ComplexMatrix matrix) {
  _checkRectangularShape2D(parts);
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
  for (int c = 0; c < columns; c++) {
    int maxWidth = 0;
    for (int r = 0; r < rows; r++) {
      ComplexMatrix part = parts[r][c];
      if (part != null) {
        int width = part.columns;
        if (maxWidth > 0 && width > 0 && width != maxWidth) {
          throw new ArgumentError("Different number of columns.");
        }
        maxWidth = Math.max(maxWidth, width);
      }
    }
    maxWidths[c] = maxWidth;
  }

  // determine row height of each row
  var maxHeights = new Int32List(rows);
  for (int r = 0; r < rows; r++) {
    int maxHeight = 0;
    for (int c = 0; c < columns; c++) {
      ComplexMatrix part = parts[r][c];
      if (part != null) {
        int height = part.rows;
        if (maxHeight > 0 && height > 0 && height != maxHeight) {
          throw new ArgumentError("Different number of rows.");
        }
        maxHeight = Math.max(maxHeight, height);
      }
    }
    maxHeights[r] = maxHeight;
  }

  // shape of result parts
  int resultRows = 0;
  for (int r = 0; r < rows; r++) {
    resultRows += maxHeights[r];
  }
  int resultCols = 0;
  for (int c = 0; c < columns; c++) {
    resultCols += maxWidths[c];
  }

  if (matrix.rows < resultRows || matrix.columns < resultCols) {
    throw new ArgumentError("Parts larger than matrix.");
  }

  // copy
  int idxr = 0;
  for (int r = 0; r < rows; r++) {
    int idxc = 0;
    for (int c = 0; c < columns; c++) {
      ComplexMatrix part = parts[r][c];
      if (part != null) {
        part.copyFrom(matrix.part(idxr, idxc, part.rows, part.columns));
      }
      idxc += maxWidths[c];
    }
    idxr += maxHeights[r];
  }
}

/// Constructs a new diagonal matrix whose diagonal elements are the elements
/// of [vector].
ComplexMatrix diag(ComplexVector vector, MatrixFn make) {
  int size = vector.size;
  ComplexMatrix diag = make(size, size);
  for (int i = 0; i < size; i++) {
    diag.set(i, i, vector.get(i));
  }
  return diag;
}

/// Constructs a new vector consisting of the diagonal elements of [A].
ComplexVector diagonal(ComplexMatrix A, MatrixFn make) {
  int min = Math.min(A.rows, A.columns);
  ComplexVector diag = make(0, 0).like1D(min);
  for (int i = 0; i < min; i++) {
    diag.set(i, A.get(i, i));
  }
  return diag;
}

/// Constructs an identity matrix (having ones on the diagonal and zeros
/// elsewhere).
ComplexMatrix identity(int rowsAndColumns, MatrixFn make) {
  ComplexMatrix matrix = make(rowsAndColumns, rowsAndColumns);
  for (int i = rowsAndColumns; --i >= 0; ) {
    matrix.set(i, i, Complex.ONE);
  }
  return matrix;
}

/// Constructs a matrix with the given shape, each cell initialized with the
/// given value.
ComplexMatrix fillMatrix(int rows, int columns, Complex initialValue, MatrixFn make) {
  if (initialValue.real == 0 && initialValue.imaginary == 0) {
    return make(rows, columns);
  }
  return make(rows, columns)..fill(initialValue.real, initialValue.imaginary);
}

/// Constructs a matrix with uniformly distributed values in `(0,1)`
/// (exclusive).
ComplexMatrix randomMatrix(int rows, int columns, MatrixFn make) {
  return make(rows, columns)..apply(cfunc.random);
}

/// `C = A||A||..||A`; Constructs a new matrix which is duplicated both along
/// the row and column dimension.
ComplexMatrix repeatMatrix(ComplexMatrix A, int rowRepeat, int columnRepeat, MatrixFn make) {
  int r = A.rows;
  int c = A.columns;
  ComplexMatrix matrix = make(r * rowRepeat, c * columnRepeat);
  for (int i = 0; i < rowRepeat; i++) {
    for (int j = 0; j < columnRepeat; j++) {
      matrix.part(r * i, c * j, r, c).copyFrom(A);
    }
  }
  return matrix;
}
