part of cern.colt.matrix.double.test;

const MAX_INT = (1 << 32) - 1;

testSparseDoubleMatrix(String name, SparseDoubleMatrixTest t) {
  testDoubleMatrix(name, t);
  group('SparseDoubleMatrix', () {
    setUp(t.setUp);
    tearDown(t.tearDown);
    test('rowCompressed', t.testRowCompressed);
    test('columnCompressed', t.testColumnCompressed);
  });
}

class SparseDoubleVectorTest extends DoubleVectorTest {
  void createMatrices() {
    A = new SparseDoubleVector(SIZE);
    B = new SparseDoubleVector(SIZE);
  }
}

class SparseDoubleVectorViewTest extends SparseDoubleVectorTest {
  void createMatrices() {
    A = new SparseDoubleVector(SIZE).flip();
    B = new SparseDoubleVector(SIZE).flip();
  }
}

class SparseCCDoubleMatrixTest extends DoubleMatrixTest {
  void createMatrices() {
    A = new SparseCCDoubleMatrix(NROWS, NCOLUMNS);
    B = new SparseCCDoubleMatrix(NROWS, NCOLUMNS);
    Bt = new SparseCCDoubleMatrix(NCOLUMNS, NROWS);
  }
}

class SparseCCDoubleMatrixViewTest extends SparseCCDoubleMatrixTest {
  void createMatrices() {
    A = new SparseCCDoubleMatrix(NCOLUMNS, NROWS).dice();
    B = new SparseCCDoubleMatrix(NCOLUMNS, NROWS).dice();
    Bt = new SparseCCDoubleMatrix(NROWS, NCOLUMNS).dice();
  }
}

class SparseDoubleMatrixTest extends DoubleMatrixTest {

  void createMatrices() {
    A = new SparseDoubleMatrix(NROWS, NCOLUMNS);
    B = new SparseDoubleMatrix(NROWS, NCOLUMNS);
    Bt = new SparseDoubleMatrix(NCOLUMNS, NROWS);
  }

  void testRowCompressed() {
    int SIZE = NROWS * NCOLUMNS;
    Int32List rowindexes = new Int32List(SIZE);
    Int32List columnindexes = new Int32List(SIZE);
    Float64List values = new Float64List(SIZE);
    for (int i = 0; i < SIZE; i++) {
      rowindexes[i] = (random.nextInt(MAX_INT) % NROWS).abs();
      columnindexes[i] = (random.nextInt(MAX_INT) % NCOLUMNS).abs();
      values[i] = random.nextDouble();
    }
    SparseDoubleMatrix A = new SparseDoubleMatrix.withValues(NROWS, NCOLUMNS, rowindexes, columnindexes, values);
    SparseRCDoubleMatrix B = A.rowCompressed(false);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(A.get(r, c), equals(B.get(r, c)));
      }
    }
    B = A.rowCompressed(true);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(A.get(r, c), equals(B.get(r, c)));
      }
    }
  }

  void testColumnCompressed() {
    int SIZE = A.rows * A.columns;
    Int32List rowindexes = new Int32List(SIZE);
    Int32List columnindexes = new Int32List(SIZE);
    Float64List values = new Float64List(SIZE);
    for (int i = 0; i < SIZE; i++) {
      rowindexes[i] = (random.nextInt(MAX_INT) % NROWS).abs();
      columnindexes[i] = (random.nextInt(MAX_INT) % NCOLUMNS).abs();
      values[i] = random.nextDouble();
    }
    SparseDoubleMatrix S = new SparseDoubleMatrix.withValues(A.rows, A.columns, rowindexes, columnindexes, values);
    SparseCCDoubleMatrix B = S.columnCompressed(false);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(S.get(r, c), equals(B.get(r, c)));
      }
    }
    B = S.columnCompressed(true);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(S.get(r, c), equals(B.get(r, c)));
      }
    }
  }
}

class SparseDoubleMatrixViewTest extends SparseDoubleMatrixTest {
  void createMatrices() {
    A = new SparseDoubleMatrix(NCOLUMNS, NROWS).dice();
    B = new SparseDoubleMatrix(NCOLUMNS, NROWS).dice();
    Bt = new SparseDoubleMatrix(NROWS, NCOLUMNS).dice();
  }
}

class SparseRCDoubleMatrixTest extends SparseDoubleMatrixTest {
  void createMatrices() {
    A = new SparseRCDoubleMatrix(NROWS, NCOLUMNS);
    B = new SparseRCDoubleMatrix(NROWS, NCOLUMNS);
    Bt = new SparseRCDoubleMatrix(NCOLUMNS, NROWS);
  }
}

class SparseRCDoubleMatrixViewTest extends SparseDoubleMatrixTest {
  void createMatrices() {
    A = new SparseRCDoubleMatrix(NCOLUMNS, NROWS).dice();
    B = new SparseRCDoubleMatrix(NCOLUMNS, NROWS).dice();
    Bt = new SparseRCDoubleMatrix(NROWS, NCOLUMNS).dice();
  }
}
