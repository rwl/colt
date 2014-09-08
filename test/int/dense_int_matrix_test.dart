part of cern.colt.matrix.int.test;

testDenseIntMatrix(String name, DenseIntMatrixTest t) {
  testIntMatrix(name, t);
  group('DenseIntMatrix', () {
    setUp(t.setUp);
    tearDown(t.tearDown);
    test('setAll', t.testSetAll);
  });
}

class DenseIntMatrixTest extends IntMatrixTest {

  void createMatrices() {
    A = new IntMatrix(NROWS, NCOLUMNS);
    B = new IntMatrix(NROWS, NCOLUMNS);
    Bt = new IntMatrix(NCOLUMNS, NROWS);
  }

  testSetAll() {
    final expected = new Int32List(A.length);
    for (int i = 0; i < A.length; i++) {
      expected[i] = random.nextInt(MAX_INT);
    }
    A.setAll(expected);
    int idx = 0;
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(0, equals((expected[idx++] - A.get(r, c)).abs()));
      }
    }
  }
}

class DenseIntMatrixViewTest extends DenseIntMatrixTest {

  void createMatrices() {
    A = new IntMatrix(NCOLUMNS, NROWS).dice();
    B = new IntMatrix(NCOLUMNS, NROWS).dice();
    Bt = new IntMatrix(NROWS, NCOLUMNS).dice();
  }
}
