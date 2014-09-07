part of cern.colt.matrix.complex.test;

class SparseCCComplexMatrixTest extends ComplexMatrixTest {
  void createMatrices() {
    A = new SparseCCComplexMatrix.sized(NROWS, NCOLUMNS);
    B = new SparseCCComplexMatrix.sized(NROWS, NCOLUMNS);
    Bt = new SparseCCComplexMatrix.sized(NCOLUMNS, NROWS);
  }
}

class SparseCCComplexMatrixViewTest extends SparseCCComplexMatrixTest {
  void createMatrices() {
    A = new SparseCCComplexMatrix.sized(NCOLUMNS, NROWS).dice();
    B = new SparseCCComplexMatrix.sized(NCOLUMNS, NROWS).dice();
    Bt = new SparseCCComplexMatrix.sized(NROWS, NCOLUMNS).dice();
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
    A = new SparseRCComplexMatrix.sized(NROWS, NCOLUMNS);
    B = new SparseRCComplexMatrix.sized(NROWS, NCOLUMNS);
    Bt = new SparseRCComplexMatrix.sized(NCOLUMNS, NROWS);
  }
}

class SparseRCComplexMatrixViewTest extends SparseRCComplexMatrixTest {
  void createMatrices() {
    A = new SparseRCComplexMatrix.sized(NCOLUMNS, NROWS).dice();
    B = new SparseRCComplexMatrix.sized(NCOLUMNS, NROWS).dice();
    Bt = new SparseRCComplexMatrix.sized(NROWS, NCOLUMNS).dice();
  }
}
