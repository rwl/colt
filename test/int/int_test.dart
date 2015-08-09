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
  testIntVector('dense', DenseIntVector.create);
  testIntVector('dense view', _flip(DenseIntVector.create));

  testIntVector('sparse', SparseIntVector.create);
  testIntVector('sparse view', _flip(SparseIntVector.create));
}

intMatrixTests() {
  testIntMatrix('dense', DenseIntMatrix.create);
  testIntMatrix('dense view', _dice(DenseIntMatrix.create));
  testIntMatrix('raw', DenseIntMatrix.create);
  testIntMatrix('view', _dice(DenseIntMatrix.create));

  testIntMatrix('sparse', SparseIntMatrix.create);
  testIntMatrix('sparse view', _dice(SparseIntMatrix.create));

  testIntMatrix('sparse rc', SparseRCIntMatrix.create);
  testIntMatrix('sparse rc view', _dice(SparseRCIntMatrix.create));

  testIntMatrix('sparse cc', SparseCCIntMatrix.create);
  testIntMatrix('sparse cc view', _dice(SparseCCIntMatrix.create));
}

_flip(make) => (sz) => make(sz).flip();

_dice(make) => (r, c) => make(r, c).dice();

List<Int32List> toList(IntMatrix m) {
  var values = new List.generate(m.rows, (_) => new Int32List(m.columns));
  for (int r = 0; r < m.rows; r++) {
    Int32List currentRow = values[r];
    for (int c = 0; c < m.columns; c++) {
      currentRow[c] = m.get(r, c);
    }
  }
  return values;
}
