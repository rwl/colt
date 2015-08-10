library cern.colt.matrix.double.test;

import 'dart:math' as math;
import 'dart:typed_data';

import 'package:test/test.dart';
import 'package:colt/colt.dart';
import 'package:colt/function/double.dart' hide equals;
import 'package:colt/src/math.dart' show MAX_INT;

part 'double_vector_test.dart';
part 'double_matrix_test.dart';
part 'dense_double_matrix_test.dart';
part 'diagonal_double_matrix_test.dart';
part 'sparse_double_matrix_test.dart';

final math.Random _r = new math.Random(0);

doubleVectorTests() {
  testDoubleVector('dense', DenseDoubleVector.create);
  testDoubleVector('dense view', _flip(DenseDoubleVector.create));

  testDoubleVector('sparse', SparseDoubleVector.create);
  testDoubleVector('sparse view', _flip(SparseDoubleVector.create));
}

doubleMatrixTests() {
  testDoubleMatrix('dense', DenseDoubleMatrix.create);
  testDoubleMatrix('dense view', _dice(DenseDoubleMatrix.create));
  testDenseDoubleMatrix('raw', DenseDoubleMatrix.create);
  testDenseDoubleMatrix('view', _dice(DenseDoubleMatrix.create));

//  testDiagonalDoubleMatrix(true);
//  testDiagonalDoubleMatrix(false);

  testDoubleMatrix('sparse', SparseDoubleMatrix.create);
  testDoubleMatrix('sparse view', _dice(SparseDoubleMatrix.create));
//  testSparseDoubleMatrix();
//  testSparseDoubleMatrix();

  testDoubleMatrix('sparse rc', SparseRCDoubleMatrix.create);
  testDoubleMatrix('sparse rc view', _dice(SparseRCDoubleMatrix.create));

//  testDoubleMatrix('sparse cc', SparseCCDoubleMatrix.create);
//  testDoubleMatrix('sparse cc view', _view(SparseCCDoubleMatrix.create));
}

_flip(make) => (sz) => make(sz).flip();

_dice(make) => (r, c) => make(r, c).dice();

List<Float64List> toList(DoubleMatrix m) {
  var values = new List.generate(m.rows, (_) => new Float64List(m.columns));
  for (int r = 0; r < m.rows; r++) {
    Float64List currentRow = values[r];
    for (int c = 0; c < m.columns; c++) {
      currentRow[c] = m.get(r, c);
    }
  }
  return values;
}
