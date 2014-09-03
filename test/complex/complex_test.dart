library cern.colt.matrix.complex.test;

import 'dart:math' as math;
import 'dart:typed_data';

import 'package:unittest/unittest.dart';
import 'package:colt/colt.dart';
import 'package:colt/math.dart';
import 'package:colt/function/complex.dart' hide equals;

part 'complex_matrix_1d_test.dart';
part 'complex_matrix_2d_test.dart';
part 'dense_complex_matrix_test.dart';
part 'sparse_complex_matrix_test.dart';
part 'diagonal_complex_matrix_test.dart';

final math.Random random = new math.Random(0);

complexMatrix1DTests() {
  testDComplexMatrix1D('DenseDComplexMatrix1D', new DenseDComplexMatrix1DTest());
  testDComplexMatrix1D('DenseDComplexMatrix1D viewFlip', new DenseDComplexMatrix1DViewTest());

  testDComplexMatrix1D('SparseDComplexMatrix1D', new SparseDComplexMatrix1DTest());
  testDComplexMatrix1D('SparseDComplexMatrix1D viewFlip', new SparseDComplexMatrix1DViewTest());
}

complexMatrix2DTests() {
  testDComplexMatrix2D('DenseDComplexMatrix2D', new DenseDComplexMatrix2DTest());
  testDComplexMatrix2D('DenseDComplexMatrix2D viewDice', new DenseDComplexMatrix2DViewTest());
  testDComplexMatrix2D('DiagonalDComplexMatrix2D', new DiagonalDComplexMatrix2DTest());
  testDComplexMatrix2D('DiagonalDComplexMatrix2D viewDice', new DiagonalDComplexMatrix2DViewTest());

  testDComplexMatrix2D('SparseCCDComplexMatrix2D', new SparseCCDComplexMatrix2DTest());
  testDComplexMatrix2D('SparseCCDComplexMatrix2D viewDice', new SparseCCDComplexMatrix2DViewTest());
  testDComplexMatrix2D('SparseDComplexMatrix2D', new SparseDComplexMatrix2DTest());
  testDComplexMatrix2D('SparseDComplexMatrix2D viewDice', new SparseDComplexMatrix2DViewTest());
  testDComplexMatrix2D('SparseRCDComplexMatrix2D', new SparseRCDComplexMatrix2DTest());
  testDComplexMatrix2D('SparseRCDComplexMatrix2D viewDice', new SparseRCDComplexMatrix2DViewTest());
}
