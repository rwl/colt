part of cern.colt.matrix.double.test;

testDenseDoubleMatrix2D(String name, DenseDoubleMatrix2DTest t) {
  testDoubleMatrix2D(name, t);
  group('DenseDoubleMatrix2D', () {
    setUp(t.setUp);
    tearDown(t.tearDown);
    test('assignValues', t.testAssignValues);
  });
}

class DenseDoubleMatrix2DTest extends DoubleMatrix2DTest {

  void createMatrices() {
    A = new DenseDoubleMatrix2D(NROWS, NCOLUMNS);
    B = new DenseDoubleMatrix2D(NROWS, NCOLUMNS);
    Bt = new DenseDoubleMatrix2D(NCOLUMNS, NROWS);
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

class DenseDoubleMatrix2DViewTest extends DenseDoubleMatrix2DTest {

  void createMatrices() {
    A = new DenseDoubleMatrix2D(NCOLUMNS, NROWS).dice();
    B = new DenseDoubleMatrix2D(NCOLUMNS, NROWS).dice();
    Bt = new DenseDoubleMatrix2D(NROWS, NCOLUMNS).dice();
  }

}
