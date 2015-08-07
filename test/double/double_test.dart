library cern.colt.matrix.double.test;

import 'dart:math' as math;
import 'dart:typed_data';

import 'package:test/test.dart';
import 'package:colt/colt.dart';
import 'package:colt/function/double.dart' hide equals;
import 'package:colt/math.dart' show MAX_INT;

part 'double_vector_test.dart';
part 'double_matrix_test.dart';
part 'dense_double_matrix_test.dart';
part 'diagonal_double_matrix_test.dart';
part 'sparse_double_matrix_test.dart';

final math.Random _r = new math.Random(0);

doubleVectorTests() {
  testAbstractDoubleVector('dense', DoubleVector.create);
  testAbstractDoubleVector('dense view', _flip(DoubleVector.create));

  testAbstractDoubleVector('sparse', SparseDoubleVector.create);
  testAbstractDoubleVector('sparse view', _flip(SparseDoubleVector.create));
}

doubleMatrixTests() {
  testAbstractDoubleMatrix('dense', DoubleMatrix.create);
  testAbstractDoubleMatrix('dense view', _dice(DoubleMatrix.create));
  testDenseDoubleMatrix('raw', DoubleMatrix.create);
  testDenseDoubleMatrix('view', _dice(DoubleMatrix.create));

//  testDiagonalDoubleMatrix(true);
//  testDiagonalDoubleMatrix(false);

  testAbstractDoubleMatrix('sparse', SparseDoubleMatrix.create);
  testAbstractDoubleMatrix('sparse view', _dice(SparseDoubleMatrix.create));
//  testSparseDoubleMatrix();
//  testSparseDoubleMatrix();

  testAbstractDoubleMatrix('sparse rc', SparseRCDoubleMatrix.create);
  testAbstractDoubleMatrix('sparse rc view', _dice(SparseRCDoubleMatrix.create));

//  testAbstractDoubleMatrix('sparse cc', SparseCCDoubleMatrix.create);
//  testAbstractDoubleMatrix('sparse cc view', _view(SparseCCDoubleMatrix.create));
}

_flip(make) => (sz) => make(sz).flip();

_dice(make) => (r, c) => make(r, c).dice();

List<Float64List> toList(AbstractDoubleMatrix m) {
  var values = new List.generate(m.rows, (_) => new Float64List(m.columns));
  for (int r = 0; r < m.rows; r++) {
    Float64List currentRow = values[r];
    for (int c = 0; c < m.columns; c++) {
      currentRow[c] = m.get(r, c);
    }
  }
  return values;
}
