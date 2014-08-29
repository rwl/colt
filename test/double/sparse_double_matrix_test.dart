part of cern.colt.matrix.double.test;

const MAX_INT = (1 << 32) - 1;

class SparseDoubleMatrix1DTest extends DoubleMatrix1DTest {
  void createMatrices() {
    A = new SparseDoubleMatrix1D(SIZE);
    B = new SparseDoubleMatrix1D(SIZE);
  }
}

class SparseDoubleMatrix1DViewTest extends SparseDoubleMatrix1DTest {
  void createMatrices() {
    A = new SparseDoubleMatrix1D(SIZE).viewFlip();
    B = new SparseDoubleMatrix1D(SIZE).viewFlip();
  }
}

class SparseCCDoubleMatrix2DTest extends DoubleMatrix2DTest {
  void createMatrices() {
    A = new SparseCCDoubleMatrix2D.sized(NROWS, NCOLUMNS);
    B = new SparseCCDoubleMatrix2D.sized(NROWS, NCOLUMNS);
    Bt = new SparseCCDoubleMatrix2D.sized(NCOLUMNS, NROWS);
  }
}

class SparseCCDoubleMatrix2DViewTest extends SparseCCDoubleMatrix2DTest {
  void createMatrices() {
    A = new SparseCCDoubleMatrix2D.sized(NCOLUMNS, NROWS).viewDice();
    B = new SparseCCDoubleMatrix2D.sized(NCOLUMNS, NROWS).viewDice();
    Bt = new SparseCCDoubleMatrix2D.sized(NROWS, NCOLUMNS).viewDice();
  }
}

class SparseDoubleMatrix2DTest extends DoubleMatrix2DTest {

  void createMatrices() {
    A = new SparseDoubleMatrix2D(NROWS, NCOLUMNS);
    B = new SparseDoubleMatrix2D(NROWS, NCOLUMNS);
    Bt = new SparseDoubleMatrix2D(NCOLUMNS, NROWS);
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
    SparseDoubleMatrix2D A = new SparseDoubleMatrix2D.values(NROWS, NCOLUMNS, rowindexes, columnindexes, values);
    SparseRCDoubleMatrix2D B = A.getRowCompressed(false);
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        expect(A.getQuick(r, c), equals(B.getQuick(r, c)));
      }
    }
    B = A.getRowCompressed(true);
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        expect(A.getQuick(r, c), equals(B.getQuick(r, c)));
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
    SparseCCDoubleMatrix2D B = S.getColumnCompressed(false);
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        expect(S.getQuick(r, c), equals(B.getQuick(r, c)));
      }
    }
    B = S.getColumnCompressed(true);
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        expect(S.getQuick(r, c), equals(B.getQuick(r, c)));
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
    int SIZE = A.rows() * A.columns();
    Int32List rowindexes = new Int32List(SIZE);
    Int32List columnindexes = new Int32List(SIZE);
    Float64List values = new Float64List(SIZE);
    DoubleMatrix2D Adense = new DenseDoubleMatrix2D(A.rows(), A.columns());
    for (int i = 0; i < SIZE; i++) {
      rowindexes[i] = i % A.rows();
      columnindexes[i] = i % A.columns();
      values[i] = random.nextDouble();
      Adense.setQuick(rowindexes[i], columnindexes[i], values[i]);
    }
    SparseDoubleMatrix2D S = new SparseDoubleMatrix2D(A.rows(), A.columns());
    S.assignIndexValuesFunc(rowindexes, columnindexes, values, multSecond(2.0));
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        expect(2 * Adense.getQuick(r, c), equals(S.getQuick(r, c)));
      }
    }
  }
}

class SparseDoubleMatrix2DViewTest extends SparseDoubleMatrix2DTest {
  void createMatrices() {
    A = new SparseDoubleMatrix2D(NCOLUMNS, NROWS).viewDice();
    B = new SparseDoubleMatrix2D(NCOLUMNS, NROWS).viewDice();
    Bt = new SparseDoubleMatrix2D(NROWS, NCOLUMNS).viewDice();
  }
}

class SparseRCDoubleMatrix2DTest extends DoubleMatrix2DTest {
  void createMatrices() {
    A = new SparseRCDoubleMatrix2D.sized(NROWS, NCOLUMNS);
    B = new SparseRCDoubleMatrix2D.sized(NROWS, NCOLUMNS);
    Bt = new SparseRCDoubleMatrix2D.sized(NCOLUMNS, NROWS);
  }
}

class SparseRCDoubleMatrix2DViewTest extends SparseRCDoubleMatrix2DTest {
  void createMatrices() {
    A = new SparseRCDoubleMatrix2D.sized(NCOLUMNS, NROWS).viewDice();
    B = new SparseRCDoubleMatrix2D.sized(NCOLUMNS, NROWS).viewDice();
    Bt = new SparseRCDoubleMatrix2D.sized(NROWS, NCOLUMNS).viewDice();
  }
}
