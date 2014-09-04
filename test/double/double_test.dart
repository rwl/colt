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
part 'sparse_double_matrix_test.dart';

final math.Random random = new math.Random(0);

doubleMatrix1DTests() {
  testDoubleMatrix1D('DenseDoubleMatrix1D', new DenseDoubleMatrix1DTest());
  testDoubleMatrix1D('DenseDoubleMatrix1D viewFlip', new DenseDoubleMatrix1DViewTest());

  testDoubleMatrix1D('SparseDoubleMatrix1D', new SparseDoubleMatrix1DTest());
  testDoubleMatrix1D('SparseDoubleMatrix1D viewFlip', new SparseDoubleMatrix1DViewTest());
}

doubleMatrix2DTests() {
  testDenseDoubleMatrix2D('DenseDoubleMatrix2D', new DenseDoubleMatrix2DTest());
  testDenseDoubleMatrix2D('DenseDoubleMatrix2D viewDice', new DenseDoubleMatrix2DViewTest());
  testDoubleMatrix2D('DiagonalDoubleMatrix2D', new DiagonalDoubleMatrix2DTest());
  testDoubleMatrix2D('DiagonalDoubleMatrix2D viewDice', new DiagonalDoubleMatrix2DViewTest());

  testDoubleMatrix2D('SparseCCDoubleMatrix2D', new SparseCCDoubleMatrix2DTest());
  testDoubleMatrix2D('SparseCCDoubleMatrix2D viewDice', new SparseCCDoubleMatrix2DViewTest());
  testDoubleMatrix2D('SparseDoubleMatrix2D', new SparseDoubleMatrix2DTest());
  testDoubleMatrix2D('SparseDoubleMatrix2D viewDice', new SparseDoubleMatrix2DViewTest());
  testDoubleMatrix2D('SparseRCDoubleMatrix2D', new SparseRCDoubleMatrix2DTest());
  testDoubleMatrix2D('SparseRCDoubleMatrix2D viewDice', new SparseRCDoubleMatrix2DViewTest());
}
