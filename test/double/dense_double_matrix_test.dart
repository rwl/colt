part of cern.colt.matrix.double.test;

testDenseDoubleMatrix(String name, DenseDoubleMatrixTest t) {
  testDoubleMatrix(name, t);
  group('DenseDoubleMatrix', () {
    setUp(t.setUp);
    tearDown(t.tearDown);
    test('setAll', t.testSetAll);
  });
}

class DenseDoubleMatrixTest extends DoubleMatrixTest {
  void createMatrices() {
    A = new DoubleMatrix(NROWS, NCOLUMNS);
    B = new DoubleMatrix(NROWS, NCOLUMNS);
    Bt = new DoubleMatrix(NCOLUMNS, NROWS);
  }

  void testSetAll() {
    var expected = new Float64List(A.size);
    for (int i = 0; i < A.size; i++) {
      expected[i] = random.nextDouble();
    }
    A.setAll(expected);
    int idx = 0;
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(0, closeTo((expected[idx++] - A.get(r, c)).abs(), TOL));
      }
    }
  }
}

class DenseDoubleMatrixViewTest extends DenseDoubleMatrixTest {
  void createMatrices() {
    A = new DoubleMatrix(NCOLUMNS, NROWS).dice();
    B = new DoubleMatrix(NCOLUMNS, NROWS).dice();
    Bt = new DoubleMatrix(NROWS, NCOLUMNS).dice();
  }
}
