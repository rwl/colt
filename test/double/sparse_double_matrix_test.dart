part of cern.colt.matrix.double.test;

const MAX_INT = (1 << 32) - 1;

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

  void testGetRowCompressed() {
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
    SparseRCDoubleMatrix B = A.getRowCompressed(false);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(A.get(r, c), equals(B.get(r, c)));
      }
    }
    B = A.getRowCompressed(true);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(A.get(r, c), equals(B.get(r, c)));
      }
    }
  }

  /*void testGetRowCompressedModified() {
    int SIZE = A.rows() * A.columns();
    Int32List rowindexes = new Int32List(SIZE);
    Int32List columnindexes = new Int32List(SIZE);
    Float64List values = new Float64List(SIZE);
    for (int i = 0; i < SIZE; i++) {
      rowindexes[i] = (random.nextInt(MAX_INT) % NROWS).abs();
      columnindexes[i] = (random.nextInt(MAX_INT) % NCOLUMNS).abs();
      values[i] = random.nextDouble();
    }
    SparseDoubleMatrix2D S = new SparseDoubleMatrix2D.values(A.rows(), A.columns(), rowindexes, columnindexes, values);
    SparseRCMDoubleMatrix2D B = S.getRowCompressedModified();
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        expect(S.getQuick(r, c), equals(B.getQuick(r, c)));
      }
    }
  }*/

  void testGetColumnCompressed() {
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
    SparseCCDoubleMatrix B = S.getColumnCompressed(false);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(S.get(r, c), equals(B.get(r, c)));
      }
    }
    B = S.getColumnCompressed(true);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(S.get(r, c), equals(B.get(r, c)));
      }
    }

  }

  /*void testGetColumnCompressedModified() {
    int SIZE = A.rows() * A.columns();
    Int32List rowindexes = new Int32List(SIZE);
    Int32List columnindexes = new Int32List(SIZE);
    Float64List values = new Float64List(SIZE);
    for (int i = 0; i < SIZE; i++) {
      rowindexes[i] = (random.nextInt(MAX_INT) % NROWS).abs();
      columnindexes[i] = (random.nextInt(MAX_INT) % NCOLUMNS).abs();
      values[i] = random.nextDouble();
    }
    SparseDoubleMatrix2D S = new SparseDoubleMatrix2D.values(A.rows(), A.columns(), rowindexes, columnindexes, values);
    SparseCCMDoubleMatrix2D B = S.getColumnCompressedModified();
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        expect(S.getQuick(r, c), equals(B.getQuick(r, c)));
      }
    }
  }*/

  void testAssignIntArrayIntArrayDoubleArrayDoubleDoubleFunction() {
    int SIZE = A.rows * A.columns;
    Int32List rowindexes = new Int32List(SIZE);
    Int32List columnindexes = new Int32List(SIZE);
    Float64List values = new Float64List(SIZE);
    AbstractDoubleMatrix Adense = new DoubleMatrix(A.rows, A.columns);
    for (int i = 0; i < SIZE; i++) {
      rowindexes[i] = i % A.rows;
      columnindexes[i] = i % A.columns;
      values[i] = random.nextDouble();
      Adense.set(rowindexes[i], columnindexes[i], values[i]);
    }
    SparseDoubleMatrix S = new SparseDoubleMatrix(A.rows, A.columns);
    S.assignIndexValuesFunc(rowindexes, columnindexes, values, multSecond(2.0));
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(2 * Adense.get(r, c), equals(S.get(r, c)));
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

class SparseRCDoubleMatrixTest extends DoubleMatrixTest {
  void createMatrices() {
    A = new SparseRCDoubleMatrix.sized(NROWS, NCOLUMNS);
    B = new SparseRCDoubleMatrix.sized(NROWS, NCOLUMNS);
    Bt = new SparseRCDoubleMatrix.sized(NCOLUMNS, NROWS);
  }
}

class SparseRCDoubleMatrixViewTest extends SparseRCDoubleMatrixTest {
  void createMatrices() {
    A = new SparseRCDoubleMatrix.sized(NCOLUMNS, NROWS).dice();
    B = new SparseRCDoubleMatrix.sized(NCOLUMNS, NROWS).dice();
    Bt = new SparseRCDoubleMatrix.sized(NROWS, NCOLUMNS).dice();
  }
}
