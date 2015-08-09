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

import 'package:complex/complex.dart';

import 'matrix.dart';
//import '../formatter.dart';

import '../../math.dart' show EPSILON;
import '../../math.dart' as cmath;

/*void _checkDense(ComplexVector A) {
  if (A is! ComplexVector) {
    throw new ArgumentError("Matrix must be dense");
  }
}*/

/// Checks whether the given matrix `A` is square.
/*void _checkSquare(ComplexMatrix A) {
  if (A.rows != A.columns) {
    throw new ArgumentError(
        "Matrix must be square: " + AbstractFormatter.shapeMatrix(A));
  }
}*/

/*void _checkSparse(ComplexMatrix A) {
  if (A is! SparseCCComplexMatrix && A is! SparseRCComplexMatrix) {
    throw new ArgumentError("Matrix must be sparse");
  }
}*/

/// Returns whether all cells of the given matrix [A] are equal to the
/// given value.
bool allVector(final ComplexVector A, final Complex value,
    [double epsilon = EPSILON]) {
  if (A == null) {
    return false;
  }
  for (int i = 0; i < A.size; i++) {
    var x = A.get(i);
    var diff = new Complex((value.real - x.real).abs(), (value.imaginary - x.imaginary).abs());
    if (((diff.real != diff.real) || (diff.imaginary != diff.imaginary)) &&
            ((((value.real != value.real) || (value.imaginary != value.imaginary)) &&
                ((x.real != x.real) || (x.imaginary != x.imaginary)))) ||
        (cmath.isEqual(value, x, epsilon))) {
      diff = Complex.ZERO;
    }
    if ((diff.real > epsilon) || (diff.imaginary > epsilon)) {
      return false;
    }
  }
  return true;
}

/// Returns whether both given matrices [A] and [B] are equal.
bool equalsVector(final ComplexVector A, final ComplexVector B,
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
    var x = A.get(i);
    var value = B.get(i);
    var diff = new Complex((value.real - x.real).abs(), (value.imaginary - x.imaginary).abs());
    if (((diff.real != diff.real) || (diff.imaginary != diff.imaginary)) &&
            ((((value.real != value.real) || (value.imaginary != value.imaginary)) &&
                ((x.real != x.real) || (x.imaginary != x.imaginary)))) ||
        (cmath.isEqual(value, x, epsilon))) {
      diff = Complex.ZERO;
    }
    if ((diff.real > epsilon) || (diff.imaginary > epsilon)) {
      return false;
    }
  }
  return true;
}

/// Returns whether all cells of the given matrix [A] are equal to the
/// given value.
bool allMatrix(final ComplexMatrix A, final Complex value,
    [double epsilon = EPSILON]) {
  if (A == null) {
    return false;
  }
  int rows = A.rows;
  int columns = A.columns;
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < columns; c++) {
      var x = A.get(r, c);
      var diff = new Complex((value.real - x.real).abs(), (value.imaginary - x.imaginary).abs());
      if (((diff.real != diff.real) || (diff.imaginary != diff.imaginary)) &&
              ((((value.real != value.real) || (value.imaginary != value.imaginary)) &&
                  ((x.real != x.real) || (x.imaginary != x.imaginary)))) ||
          (cmath.isEqual(value, x, epsilon))) {
        diff = Complex.ZERO;
      }
      if ((diff.real > epsilon) || (diff.imaginary > epsilon)) {
        return false;
      }
    }
  }
  return true;
}

/// Returns whether both given matrices [A] and [B] are equal.
bool equalsMatrix(final ComplexMatrix A, final ComplexMatrix B,
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
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < columns; c++) {
      var x = A.get(r, c);
      var value = B.get(r, c);
      var diff = new Complex((value.real - x.real).abs(), (value.imaginary - x.imaginary).abs());
      if (((diff.real != diff.real) || (diff.imaginary != diff.imaginary)) &&
              ((((value.real != value.real) || (value.imaginary != value.imaginary)) &&
                  ((x.real != x.real) || (x.imaginary != x.imaginary)))) ||
          (cmath.isEqual(value, x, epsilon))) {
        diff = Complex.ZERO;
      }
      if ((diff.real > epsilon) || (diff.imaginary > epsilon)) {
        return false;
      }
    }
  }
  return true;
}
