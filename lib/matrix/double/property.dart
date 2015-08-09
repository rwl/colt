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
library cern.colt.matrix.tdouble.property;

import 'matrix.dart';
import '../../math.dart';

/// Returns a String with [length] blanks.
/*String _blanks(int length) {
  if (length < 0) {
    length = 0;
  }
  var buf = new StringBuffer(length);
  for (int k = 0; k < length; k++) {
    buf.write(' ');
  }
  return buf.toString();
}*/

/// Checks whether the given matrix [A] is rectangular.
/*void _checkRectangular(AbstractDoubleMatrix A) {
  if (A.rows < A.columns) {
    throw new ArgumentError(
        "Matrix must be rectangular: " + AbstractFormatter.shapeMatrix(A));
  }
}*/

/// Checks whether the given matrix `A` is square.
/*void _checkSquare(AbstractDoubleMatrix A) {
  if (A.rows != A.columns) {
    throw new ArgumentError(
        "Matrix must be square: " + AbstractFormatter.shapeMatrix(A));
  }
}*/

/*void _checkDense(AbstractDoubleVector A) {
  if (A is! DoubleVector) {
    throw new ArgumentError("Matrix must be dense");
  }
}

void _checkSparse(AbstractDoubleVector A) {
  if (A is! SparseDoubleVector) {
    throw new ArgumentError("Matrix must be sparse");
  }
}

void _checkSparseMatrix(AbstractDoubleMatrix A) {
  if (A is! SparseCCDoubleMatrix &&
      A is! SparseRCDoubleMatrix &&
      A is! SparseDoubleMatrix) {
    throw new ArgumentError("Matrix must be sparse");
  }
}*/

/// Returns the matrix's fraction of non-zero cells;
/// `A.cardinality / A.size`.
double density(AbstractDoubleMatrix A) => A.cardinality / A.size.toDouble();

/// Returns whether all cells of the given matrix [A] are equal to the
/// given value. The result is `true` if and only if
/// `A != null` and `! (Math.abs(value - A[i]) > tolerance())`
/// holds for all coordinates.
bool allVector(final AbstractDoubleVector A, final double value,
    [double epsilon = EPSILON]) {
  if (A == null) {
    return false;
  }
  int size = A.size;
  for (int i = 0; i < size; i++) {
    double x = A.get(i);
    double diff = (value - x).abs();
    if ((diff != diff) && ((value != value && x != x) || value == x)) {
      diff = 0.0;
    }
    if (!(diff <= epsilon)) {
      return false;
    }
  }
  return true;
}

/// Returns whether both given matrices [A] and [B] are equal.
/// The result is `true` if `A==B`. Otherwise, the result is
/// `true` if and only if both arguments are `!= null`, have
/// the same size and `! (Math.abs(A[i] - B[i]) > tolerance())`
/// holds for all indexes.
bool equalsVector(final AbstractDoubleVector A, final AbstractDoubleVector B,
    [double epsilon = EPSILON]) {
  if (identical(A, B)) {
    return true;
  }
  if (!(A != null && B != null)) {
    return false;
  }
  int size = A.size;
  if (size != B.size) {
    return false;
  }

  for (int i = 0; i < size; i++) {
    double x = A.get(i);
    double value = B.get(i);
    double diff = (value - x).abs();
    if ((diff != diff) && ((value != value && x != x) || value == x)) {
      diff = 0.0;
    }
    if (!(diff <= epsilon)) {
      return false;
    }
  }
  return true;
}

/// Returns whether all cells of the given matrix [A] are equal to the
/// given value. The result is `true` if and only if `A != null` and
/// `! (Math.abs(value - A[row,col]) > tolerance())` holds for all
/// coordinates.
bool allMatrix(final AbstractDoubleMatrix A, final double value,
    [double epsilon = EPSILON]) {
  if (A == null) {
    return false;
  }
  final int rows = A.rows;
  final int columns = A.columns;
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < columns; c++) {
      double x = A.get(r, c);
      double diff = (value - x).abs();
      if ((diff != diff) && ((value != value && x != x) || value == x)) {
        diff = 0.0;
      }
      if (!(diff <= epsilon)) {
        return false;
      }
    }
  }
  return true;
}

/// Returns whether both given matrices [A] and [B] are equal.
/// The result is `true` if `A==B`. Otherwise, the result is
/// `true` if and only if both arguments are `!= null`, have
/// the same number of columns and rows and
/// `! (Math.abs(A[row,col] - B[row,col]) > tolerance())` holds
/// for all coordinates.
bool equalsMatrix(final AbstractDoubleMatrix A, final AbstractDoubleMatrix B,
    [double epsilon = EPSILON]) {
  if (identical(A, B)) {
    return true;
  }
  if (A == null || B == null) {
    return false;
  }
  final int rows = A.rows;
  final int columns = A.columns;
  if (columns != B.columns || rows != B.rows) {
    return false;
  }
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < columns; c++) {
      double x = A.get(r, c);
      double value = B.get(r, c);
      double diff = (value - x).abs();
      if ((diff != diff) && ((value != value && x != x) || value == x)) diff =
          0.0;
      if (!(diff <= epsilon)) {
        return false;
      }
    }
  }
  return true;
}
