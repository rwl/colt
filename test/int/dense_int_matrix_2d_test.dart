part of cern.colt.matrix.int.test;

class DenseIntMatrix2DTest extends IntMatrix2DTest {

  void createMatrices() {
    A = new DenseIntMatrix2D(NROWS, NCOLUMNS);
    B = new DenseIntMatrix2D(NROWS, NCOLUMNS);
    Bt = new DenseIntMatrix2D(NCOLUMNS, NROWS);
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

class DenseIntMatrix2DViewTest extends DenseIntMatrix2DTest {

  void createMatrices() {
    A = new DenseIntMatrix2D(NCOLUMNS, NROWS).dice();
    B = new DenseIntMatrix2D(NCOLUMNS, NROWS).dice();
    Bt = new DenseIntMatrix2D(NROWS, NCOLUMNS).dice();
  }
}
