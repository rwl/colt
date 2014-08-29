library cern.colt.matrix.double.test;

import 'dart:math' as math;
import 'dart:typed_data';

import 'package:unittest/unittest.dart';
import 'package:colt/colt.dart';
import 'package:colt/function/double.dart' hide equals;

part 'double_matrix_1d_test.dart';
part 'double_matrix_2d_test.dart';
part 'dense_double_matrix_1d_test.dart';
part 'dense_double_matrix_2d_test.dart';
part 'diagonal_double_matrix_2d_test.dart';

final math.Random random = new math.Random(0);

doubleMatrix1DTests() {
  testDoubleMatrix1D(new DenseDoubleMatrix1DTest());
  testDoubleMatrix1D(new DenseDoubleMatrix1DViewTest());
}

doubleMatrix2DTests() {
  testDenseDoubleMatrix2D(new DenseDoubleMatrix2DTest());
  testDenseDoubleMatrix2D(new DenseDoubleMatrix2DViewTest());
  testDoubleMatrix2D(new DiagonalDoubleMatrix2DTest());
  testDoubleMatrix2D(new DiagonalDoubleMatrix2DViewTest());
}
