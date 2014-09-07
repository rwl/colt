part of cern.colt.matrix.complex.test;

class DenseComplexVectorTest extends ComplexVectorTest {
    void createMatrices() {
        A = new DenseComplexVector(SIZE);
        B = new DenseComplexVector(SIZE);
    }
}

class DenseComplexVectorViewTest extends DenseComplexVectorTest {
    void createMatrices() {
        A = new DenseComplexVector(SIZE).flip();
        B = new DenseComplexVector(SIZE).flip();
    }
}

class DenseComplexMatrixTest extends ComplexMatrixTest {
    void createMatrices() {
        A = new DenseComplexMatrix(NROWS, NCOLUMNS);
        B = new DenseComplexMatrix(NROWS, NCOLUMNS);
        Bt = new DenseComplexMatrix(NCOLUMNS, NROWS);
    }
}

class DenseComplexMatrixViewTest extends DenseComplexMatrixTest {
    void createMatrices() {
        A = new DenseComplexMatrix(NCOLUMNS, NROWS).dice();
        B = new DenseComplexMatrix(NCOLUMNS, NROWS).dice();
        Bt = new DenseComplexMatrix(NROWS, NCOLUMNS).dice();
    }
}
