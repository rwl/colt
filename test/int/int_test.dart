library cern.colt.matrix.int.test;

import 'dart:math' as math;
import 'dart:typed_data';

import 'package:unittest/unittest.dart';
import 'package:colt/math.dart' show MAX_INT;
import 'package:colt/colt.dart';
import 'package:colt/function/int.dart' as ifunc;

part 'int_vector_test.dart';
part 'int_matrix_test.dart';
part 'dense_int_vector_test.dart';
part 'dense_int_matrix_test.dart';

final math.Random random = new math.Random(0);

intVectorTests() {
  testIntVector('DenseIntVector', new DenseIntVectorTest());
  testIntVector('DenseIntVector flip', new DenseIntVectorViewTest());
}

intMatrixTests() {
  testDenseIntMatrix('DenseIntMatrix', new DenseIntMatrixTest());
  testDenseIntMatrix('DenseIntMatrix dice', new DenseIntMatrixViewTest());
}
