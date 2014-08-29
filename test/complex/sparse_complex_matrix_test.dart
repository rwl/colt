part of cern.colt.matrix.complex.test;

class SparseCCDComplexMatrix2DTest extends DComplexMatrix2DTest {
  void createMatrices() {
    A = new SparseCCDComplexMatrix2D.sized(NROWS, NCOLUMNS);
    B = new SparseCCDComplexMatrix2D.sized(NROWS, NCOLUMNS);
    Bt = new SparseCCDComplexMatrix2D.sized(NCOLUMNS, NROWS);
  }
}

class SparseCCDComplexMatrix2DViewTest extends SparseCCDComplexMatrix2DTest {
  void createMatrices() {
    A = new SparseCCDComplexMatrix2D.sized(NCOLUMNS, NROWS).viewDice();
    B = new SparseCCDComplexMatrix2D.sized(NCOLUMNS, NROWS).viewDice();
    Bt = new SparseCCDComplexMatrix2D.sized(NROWS, NCOLUMNS).viewDice();
  }
}

class SparseDComplexMatrix1DTest extends DComplexMatrix1DTest {
  void createMatrices() {
    A = new SparseDComplexMatrix1D(SIZE);
    B = new SparseDComplexMatrix1D(SIZE);
  }
}

class SparseDComplexMatrix1DViewTest extends SparseDComplexMatrix1DTest {
  void createMatrices() {
    A = new SparseDComplexMatrix1D(SIZE).viewFlip();
    B = new SparseDComplexMatrix1D(SIZE).viewFlip();
  }
}

class SparseDComplexMatrix2DTest extends DComplexMatrix2DTest {
  void createMatrices() {
    A = new SparseDComplexMatrix2D(NROWS, NCOLUMNS);
    B = new SparseDComplexMatrix2D(NROWS, NCOLUMNS);
    Bt = new SparseDComplexMatrix2D(NCOLUMNS, NROWS);
  }
}

class SparseDComplexMatrix2DViewTest extends SparseDComplexMatrix2DTest {
  void createMatrices() {
    A = new SparseDComplexMatrix2D(NCOLUMNS, NROWS).viewDice();
    B = new SparseDComplexMatrix2D(NCOLUMNS, NROWS).viewDice();
    Bt = new SparseDComplexMatrix2D(NROWS, NCOLUMNS).viewDice();
  }
}

class SparseRCDComplexMatrix2DTest extends DComplexMatrix2DTest {
  void createMatrices() {
    A = new SparseRCDComplexMatrix2D.sized(NROWS, NCOLUMNS);
    B = new SparseRCDComplexMatrix2D.sized(NROWS, NCOLUMNS);
    Bt = new SparseRCDComplexMatrix2D.sized(NCOLUMNS, NROWS);
  }
}

class SparseRCDComplexMatrix2DViewTest extends SparseRCDComplexMatrix2DTest {
  void createMatrices() {
    A = new SparseRCDComplexMatrix2D.sized(NCOLUMNS, NROWS).viewDice();
    B = new SparseRCDComplexMatrix2D.sized(NCOLUMNS, NROWS).viewDice();
    Bt = new SparseRCDComplexMatrix2D.sized(NROWS, NCOLUMNS).viewDice();
  }
}
