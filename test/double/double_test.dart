library cern.colt.matrix.double.test;

import 'dart:math' as math;
import 'dart:typed_data';

import 'package:unittest/unittest.dart';
import 'package:colt/colt.dart';
import 'package:colt/function/double.dart' hide equals;

part 'double_vector_test.dart';
part 'double_matrix_test.dart';
part 'dense_double_vector_test.dart';
part 'dense_double_matrix_test.dart';
part 'diagonal_double_matrix_test.dart';
part 'sparse_double_matrix_test.dart';

final math.Random random = new math.Random(0);

doubleVectorTests() {
  testDoubleVector('DenseDoubleVector', new DenseDoubleVectorTest());
  testDoubleVector('DenseDoubleVector viewFlip', new DenseDoubleVectorViewTest());

  testDoubleVector('SparseDoubleVector', new SparseDoubleVectorTest());
  testDoubleVector('SparseDoubleVector viewFlip', new SparseDoubleVectorViewTest());
}

doubleMatrixTests() {
  testDenseDoubleMatrix('DenseDoubleMatrix', new DenseDoubleMatrixTest());
  testDenseDoubleMatrix('DenseDoubleMatrix viewDice', new DenseDoubleMatrixViewTest());
  testDoubleMatrix('DiagonalDoubleMatrix', new DiagonalDoubleMatrixTest());
  testDoubleMatrix('DiagonalDoubleMatrix viewDice', new DiagonalDoubleMatrixViewTest());

  testDoubleMatrix('SparseCCDoubleMatrix', new SparseCCDoubleMatrixTest());
  testDoubleMatrix('SparseCCDoubleMatrix viewDice', new SparseCCDoubleMatrixViewTest());
  testDoubleMatrix('SparseDoubleMatrix', new SparseDoubleMatrixTest());
  testDoubleMatrix('SparseDoubleMatrix viewDice', new SparseDoubleMatrixViewTest());
  testDoubleMatrix('SparseRCDoubleMatrix', new SparseRCDoubleMatrixTest());
  testDoubleMatrix('SparseRCDoubleMatrix viewDice', new SparseRCDoubleMatrixViewTest());
}
