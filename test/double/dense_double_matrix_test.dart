part of cern.colt.matrix.double.test;

testDenseDoubleMatrix(String name, DenseDoubleMatrixTest t) {
  testDoubleMatrix(name, t);
  group('DenseDoubleMatrix2D', () {
    setUp(t.setUp);
    tearDown(t.tearDown);
    test('assignValues', t.testAssignValues);
  });
}

class DenseDoubleMatrixTest extends DoubleMatrixTest {

  void createMatrices() {
    A = new DenseDoubleMatrix(NROWS, NCOLUMNS);
    B = new DenseDoubleMatrix(NROWS, NCOLUMNS);
    Bt = new DenseDoubleMatrix(NCOLUMNS, NROWS);
  }

  void testAssignValues() {
    Float64List expected = new Float64List(A.length);
    for (int i = 0; i < A.length; i++) {
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
    A = new DenseDoubleMatrix(NCOLUMNS, NROWS).dice();
    B = new DenseDoubleMatrix(NCOLUMNS, NROWS).dice();
    Bt = new DenseDoubleMatrix(NROWS, NCOLUMNS).dice();
  }

}
