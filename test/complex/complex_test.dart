library cern.colt.matrix.complex.test;

import 'dart:math' as math;
import 'dart:typed_data';

import 'package:test/test.dart';
import 'package:colt/colt.dart';
import 'package:colt/math.dart';
import 'package:colt/function/complex.dart' hide equals;

part 'complex_vector_test.dart';
part 'complex_matrix_test.dart';
part 'dense_complex_matrix_test.dart';
part 'sparse_complex_matrix_test.dart';
part 'diagonal_complex_matrix_test.dart';

final math.Random random = new math.Random(0);

complexVectorTests() {
  testComplexVector('DenseComplexVector', new DenseComplexVectorTest());
  testComplexVector('DenseComplexVector viewFlip', new DenseComplexVectorViewTest());

  testComplexVector('SparseComplexVector', new SparseComplexVectorTest());
  testComplexVector('SparseComplexVector viewFlip', new SparseComplexVectorViewTest());
}

complexMatrixTests() {
  testComplexMatrix('DenseComplexMatrix', new DenseComplexMatrixTest());
  testComplexMatrix('DenseComplexMatrix viewDice', new DenseComplexMatrixViewTest());
  testComplexMatrix('DiagonalComplexMatrix', new DiagonalComplexMatrixTest());
  testComplexMatrix('DiagonalComplexMatrix viewDice', new DiagonalComplexMatrixViewTest());

  testComplexMatrix('SparseCCComplexMatrix', new SparseCCComplexMatrixTest());
  testComplexMatrix('SparseCCComplexMatrix viewDice', new SparseCCComplexMatrixViewTest());
  testComplexMatrix('SparseComplexMatrix', new SparseComplexMatrixTest());
  testComplexMatrix('SparseComplexMatrix viewDice', new SparseComplexMatrixViewTest());
  testComplexMatrix('SparseRCComplexMatrix', new SparseRCComplexMatrixTest());
  testComplexMatrix('SparseRCComplexMatrix viewDice', new SparseRCComplexMatrixViewTest());
}
