library cern.colt.matrix.complex.test;

import 'dart:math' as math;
import 'dart:typed_data';

import 'package:test/test.dart';
import 'package:complex/complex.dart';
import 'package:colt/colt.dart';
import 'package:colt/src/math.dart' as cmath;
import 'package:colt/function/complex.dart' hide equals;

part 'complex_vector_test.dart';
part 'complex_matrix_test.dart';
part 'diagonal_complex_matrix_test.dart';

final math.Random random = new math.Random(0);

complexVectorTests() {
  testComplexVector('dense', DenseComplexVector.create);
  testComplexVector('dense view', _flip(DenseComplexVector.create));

  testComplexVector('sparse', SparseComplexVector.create);
  testComplexVector('sparse view', _flip(SparseComplexVector.create));
}

complexMatrixTests() {
  testComplexMatrix('dense', DenseComplexMatrix.create);
  testComplexMatrix('dense view', _dice(DenseComplexMatrix.create));

  testDiagonalComplexMatrix(false);
  testDiagonalComplexMatrix(true);

  testComplexMatrix('sparse', SparseComplexMatrix.create);
  testComplexMatrix('sparse view', _dice(SparseComplexMatrix.create));

  testComplexMatrix('sparse cc', SparseCCComplexMatrix.create);
  testComplexMatrix('sparse cc view', _dice(SparseCCComplexMatrix.create));

  testComplexMatrix('sparse rc ', SparseRCComplexMatrix.create);
  testComplexMatrix('sparse rc view', _dice(SparseRCComplexMatrix.create));
}

_flip(make) => (sz) => make(sz).flip();

_dice(make) => (r, c) => make(r, c).dice();

List<Float64List> toList(ComplexMatrix m) {
  var values = new List.generate(m.rows, (_) => new Float64List(2 * m.columns));
  for (int r = 0; r < m.rows; r++) {
    for (int c = 0; c < m.columns; c++) {
      var tmp = m.get(r, c);
      values[r][2 * c] = tmp.real;
      values[r][2 * c + 1] = tmp.imaginary;
    }
  }
  return values;
}
