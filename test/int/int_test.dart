library cern.colt.matrix.int.test;

import 'dart:math' as math;
import 'dart:typed_data';

import 'package:unittest/unittest.dart';
import 'package:colt/math.dart' show MAX_INT;
import 'package:colt/colt.dart';
import 'package:colt/function/int.dart' as ifunc;

part 'int_matrix_1d_test.dart';

final math.Random random = new math.Random(0);

intMatrix1DTests() {
  testIntMatrix1D('DenseIntMatrix1D', new DenseIntMatrix1DTest());
  testIntMatrix1D('DenseIntMatrix1D viewFlip', new DenseIntMatrix1DViewTest());
}

intMatrix2DTests() {
  testIntDoubleMatrix2D('DenseIntMatrix2D', new DenseIntMatrix2DTest());
  testIntDoubleMatrix2D('DenseIntMatrix2D viewDice', new DenseIntMatrix2DViewTest());
}
