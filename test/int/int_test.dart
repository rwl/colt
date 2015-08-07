library cern.colt.matrix.int.test;

import 'dart:math' as math;
import 'dart:typed_data';

import 'package:test/test.dart';
import 'package:colt/math.dart' show MAX_INT;
import 'package:colt/colt.dart';
import 'package:colt/function/int.dart' as ifunc;

part 'int_vector_test.dart';
part 'int_matrix_test.dart';
part 'dense_int_matrix_test.dart';

final math.Random random = new math.Random(0);

intVectorTests() {
  testAbstractIntVector('dense', IntVector.create);
  testAbstractIntVector('dense view', _flip(IntVector.create));

  testAbstractIntVector('sparse', SparseIntVector.create);
  testAbstractIntVector('sparse view', _flip(SparseIntVector.create));
}

intMatrixTests() {
  testAbstractIntMatrix('dense', IntMatrix.create);
  testAbstractIntMatrix('dense view', _dice(IntMatrix.create));
  testIntMatrix('raw', IntMatrix.create);
  testIntMatrix('view', _dice(IntMatrix.create));

  testAbstractIntMatrix('sparse', SparseIntMatrix.create);
  testAbstractIntMatrix('sparse view', _dice(SparseIntMatrix.create));

  testAbstractIntMatrix('sparse rc', SparseRCIntMatrix.create);
  testAbstractIntMatrix('sparse rc view', _dice(SparseRCIntMatrix.create));

  testAbstractIntMatrix('sparse cc', SparseCCIntMatrix.create);
  testAbstractIntMatrix('sparse cc view', _dice(SparseCCIntMatrix.create));
}

_flip(make) => (sz) => make(sz).flip();

_dice(make) => (r, c) => make(r, c).dice();

List<Int32List> toList(AbstractIntMatrix m) {
  var values = new List.generate(m.rows, (_) => new Int32List(m.columns));
  for (int r = 0; r < m.rows; r++) {
    Int32List currentRow = values[r];
    for (int c = 0; c < m.columns; c++) {
      currentRow[c] = m.get(r, c);
    }
  }
  return values;
}
