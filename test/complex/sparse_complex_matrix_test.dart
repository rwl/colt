part of cern.colt.matrix.complex.test;

class SparseCCComplexMatrixTest extends ComplexMatrixTest {
  void createMatrices() {
    A = new SparseCCComplexMatrix(NROWS, NCOLUMNS);
    B = new SparseCCComplexMatrix(NROWS, NCOLUMNS);
    Bt = new SparseCCComplexMatrix(NCOLUMNS, NROWS);
  }
}

class SparseCCComplexMatrixViewTest extends SparseCCComplexMatrixTest {
  void createMatrices() {
    A = new SparseCCComplexMatrix(NCOLUMNS, NROWS).dice();
    B = new SparseCCComplexMatrix(NCOLUMNS, NROWS).dice();
    Bt = new SparseCCComplexMatrix(NROWS, NCOLUMNS).dice();
  }
}

class SparseComplexVectorTest extends ComplexVectorTest {
  void createMatrices() {
    A = new SparseComplexVector(SIZE);
    B = new SparseComplexVector(SIZE);
  }
}

class SparseComplexVectorViewTest extends SparseComplexVectorTest {
  void createMatrices() {
    A = new SparseComplexVector(SIZE).flip();
    B = new SparseComplexVector(SIZE).flip();
  }
}

class SparseComplexMatrixTest extends ComplexMatrixTest {
  void createMatrices() {
    A = new SparseComplexMatrix(NROWS, NCOLUMNS);
    B = new SparseComplexMatrix(NROWS, NCOLUMNS);
    Bt = new SparseComplexMatrix(NCOLUMNS, NROWS);
  }
}

class SparseComplexMatrixViewTest extends SparseComplexMatrixTest {
  void createMatrices() {
    A = new SparseComplexMatrix(NCOLUMNS, NROWS).dice();
    B = new SparseComplexMatrix(NCOLUMNS, NROWS).dice();
    Bt = new SparseComplexMatrix(NROWS, NCOLUMNS).dice();
  }
}

class SparseRCComplexMatrixTest extends ComplexMatrixTest {
  void createMatrices() {
    A = new SparseRCComplexMatrix(NROWS, NCOLUMNS);
    B = new SparseRCComplexMatrix(NROWS, NCOLUMNS);
    Bt = new SparseRCComplexMatrix(NCOLUMNS, NROWS);
  }
}

class SparseRCComplexMatrixViewTest extends SparseRCComplexMatrixTest {
  void createMatrices() {
    A = new SparseRCComplexMatrix(NCOLUMNS, NROWS).dice();
    B = new SparseRCComplexMatrix(NCOLUMNS, NROWS).dice();
    Bt = new SparseRCComplexMatrix(NROWS, NCOLUMNS).dice();
  }
}
