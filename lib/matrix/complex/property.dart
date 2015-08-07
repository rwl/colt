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
library cern.colt.matrix.complex.property;

import 'dart:typed_data';

import 'matrix.dart';
//import '../formatter.dart';

import '../../math.dart' show EPSILON;
import '../../math.dart' as cmath;

/*void _checkDense(AbstractComplexVector A) {
  if (A is! ComplexVector) {
    throw new ArgumentError("Matrix must be dense");
  }
}*/

/// Checks whether the given matrix `A` is square.
/*void _checkSquare(AbstractComplexMatrix A) {
  if (A.rows != A.columns) {
    throw new ArgumentError(
        "Matrix must be square: " + AbstractFormatter.shapeMatrix(A));
  }
}*/

/*void _checkSparse(AbstractComplexMatrix A) {
  if (A is! SparseCCComplexMatrix && A is! SparseRCComplexMatrix) {
    throw new ArgumentError("Matrix must be sparse");
  }
}*/

/// Returns whether all cells of the given matrix [A] are equal to the
/// given value.
bool allVector(final AbstractComplexVector A, final Float64List value,
    [double epsilon = EPSILON]) {
  if (A == null) {
    return false;
  }
  var diff = new Float64List(2);
  for (int i = 0; i < A.size; i++) {
    Float64List x = A.get(i);
    diff[0] = (value[0] - x[0]).abs();
    diff[1] = (value[1] - x[1]).abs();
    if (((diff[0] != diff[0]) || (diff[1] != diff[1])) &&
            ((((value[0] != value[0]) || (value[1] != value[1])) &&
                ((x[0] != x[0]) || (x[1] != x[1])))) ||
        (cmath.isEqual(value, x, epsilon))) {
      diff[0] = 0.0;
      diff[1] = 0.0;
    }
    if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
      return false;
    }
  }
  return true;
}

/// Returns whether both given matrices [A] and [B] are equal.
bool equalsVector(final AbstractComplexVector A, final AbstractComplexVector B,
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

  var diff = new Float64List(2);
  for (int i = 0; i < size; i++) {
    Float64List x = A.get(i);
    Float64List value = B.get(i);
    diff[0] = (value[0] - x[0]).abs();
    diff[1] = (value[1] - x[1]).abs();
    if (((diff[0] != diff[0]) || (diff[1] != diff[1])) &&
            ((((value[0] != value[0]) || (value[1] != value[1])) &&
                ((x[0] != x[0]) || (x[1] != x[1])))) ||
        (cmath.isEqual(value, x, epsilon))) {
      diff[0] = 0.0;
      diff[1] = 0.0;
    }
    if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
      return false;
    }
  }
  return true;
}

/// Returns whether all cells of the given matrix [A] are equal to the
/// given value.
bool allMatrix(final AbstractComplexMatrix A, final Float64List value,
    [double epsilon = EPSILON]) {
  if (A == null) {
    return false;
  }
  int rows = A.rows;
  int columns = A.columns;
  var diff = new Float64List(2);
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < columns; c++) {
      Float64List x = A.get(r, c);
      diff[0] = (value[0] - x[0]).abs();
      diff[1] = (value[1] - x[1]).abs();
      if (((diff[0] != diff[0]) || (diff[1] != diff[1])) &&
              ((((value[0] != value[0]) || (value[1] != value[1])) &&
                  ((x[0] != x[0]) || (x[1] != x[1])))) ||
          (cmath.isEqual(value, x, epsilon))) {
        diff[0] = 0.0;
        diff[1] = 0.0;
      }
      if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
        return false;
      }
    }
  }
  return true;
}

/// Returns whether both given matrices [A] and [B] are equal.
bool equalsMatrix(final AbstractComplexMatrix A, final AbstractComplexMatrix B,
    [double epsilon = EPSILON]) {
  if (identical(A, B)) {
    return true;
  }
  if (A == null || B == null) {
    return false;
  }
  int rows = A.rows;
  int columns = A.columns;
  if (columns != B.columns || rows != B.rows) {
    return false;
  }
  var diff = new Float64List(2);
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < columns; c++) {
      Float64List x = A.get(r, c);
      Float64List value = B.get(r, c);
      diff[0] = (value[0] - x[0]).abs();
      diff[1] = (value[1] - x[1]).abs();
      if (((diff[0] != diff[0]) || (diff[1] != diff[1])) &&
              ((((value[0] != value[0]) || (value[1] != value[1])) &&
                  ((x[0] != x[0]) || (x[1] != x[1])))) ||
          (cmath.isEqual(value, x, epsilon))) {
        diff[0] = 0.0;
        diff[1] = 0.0;
      }
      if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
        return false;
      }
    }
  }
  return true;
}
