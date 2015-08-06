library cern.colt.matrix.int.property;

import 'matrix.dart';
import '../former.dart';

/// Checks whether the given matrix [A] is rectangular.
void _checkRectangular(AbstractIntMatrix A) {
  if (A.rows < A.columns) {
    throw new ArgumentError(
        "Matrix must be rectangular: " + AbstractFormatter.shape2D(A));
  }
}

/// Checks whether the given matrix [A] is square.
void _checkSquare(AbstractIntMatrix A) {
  if (A.rows != A.columns) {
    throw new ArgumentError(
        "Matrix must be square: " + AbstractFormatter.shape2D(A));
  }
}

/// Returns the matrix's fraction of non-zero cells;
/// `A.cardinality / A.size`.
int density(AbstractIntMatrix A) {
  return A.cardinality ~/ A.size;
}

/// Returns whether all cells of the given matrix [A] are equal to the
/// given value. The result is `true` if and only if
/// `A != null` and `! (Math.abs(value - A[i]) > tolerance())`
/// holds for all coordinates.
bool allVector(final AbstractIntVector A, final int value) {
  if (A == null) {
    return false;
  }
  int size = A.size;
  bool result = false;
  for (int i = 0; i < size; i++) {
    if (!(A.get(i) == value)) {
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
bool equalsVector(final AbstractIntVector A, final AbstractIntVector B) {
  if (identical(A, B)) {
    return true;
  }
  if (A == null || B == null) {
    return false;
  }
  int size = A.size;
  if (size != B.size) {
    return false;
  }

  bool result = false;
  for (int i = 0; i < size; i++) {
    if (!(A.get(i) == B.get(i))) return false;
  }
  return true;
}

/// Returns whether all cells of the given matrix [A] are equal to the
/// given value. The result is `true` if and only if `A != null` and
/// `! (Math.abs(value - A[row,col]) > tolerance())` holds for all
/// coordinates.
bool allMatrix(final AbstractIntMatrix A, final int value) {
  if (A == null) {
    return false;
  }
  final int rows = A.rows;
  final int columns = A.columns;
  bool result = false;
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < columns; c++) {
      if (!(A.get(r, c) == value)) {
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
bool equalsMatrix(final AbstractIntMatrix A, final AbstractIntMatrix B) {
  if (identical(A, B)) {
    return true;
  }
  if (!(A != null && B != null)) {
    return false;
  }
  final int rows = A.rows;
  final int columns = A.columns;
  if (columns != B.columns || rows != B.rows) {
    return false;
  }
  bool result = false;
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < columns; c++) {
      if (!(A.get(r, c) == B.get(r, c))) {
        return false;
      }
    }
  }
  return true;
}
