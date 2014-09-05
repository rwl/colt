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
    A = new SparseCCDComplexMatrix2D.sized(NCOLUMNS, NROWS).dice();
    B = new SparseCCDComplexMatrix2D.sized(NCOLUMNS, NROWS).dice();
    Bt = new SparseCCDComplexMatrix2D.sized(NROWS, NCOLUMNS).dice();
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
    A = new SparseDComplexMatrix1D(SIZE).flip();
    B = new SparseDComplexMatrix1D(SIZE).flip();
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
    A = new SparseDComplexMatrix2D(NCOLUMNS, NROWS).dice();
    B = new SparseDComplexMatrix2D(NCOLUMNS, NROWS).dice();
    Bt = new SparseDComplexMatrix2D(NROWS, NCOLUMNS).dice();
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
    A = new SparseRCDComplexMatrix2D.sized(NCOLUMNS, NROWS).dice();
    B = new SparseRCDComplexMatrix2D.sized(NCOLUMNS, NROWS).dice();
    Bt = new SparseRCDComplexMatrix2D.sized(NROWS, NCOLUMNS).dice();
  }
}
