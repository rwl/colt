part of cern.colt.matrix.complex.test;

class DenseComplexVectorTest extends ComplexVectorTest {
    void createMatrices() {
        A = new ComplexVector(SIZE);
        B = new ComplexVector(SIZE);
    }
}

class DenseComplexVectorViewTest extends DenseComplexVectorTest {
    void createMatrices() {
        A = new ComplexVector(SIZE).flip();
        B = new ComplexVector(SIZE).flip();
    }
}

class DenseComplexMatrixTest extends ComplexMatrixTest {
    void createMatrices() {
        A = new ComplexMatrix(NROWS, NCOLUMNS);
        B = new ComplexMatrix(NROWS, NCOLUMNS);
        Bt = new ComplexMatrix(NCOLUMNS, NROWS);
    }
}

class DenseComplexMatrixViewTest extends DenseComplexMatrixTest {
    void createMatrices() {
        A = new ComplexMatrix(NCOLUMNS, NROWS).dice();
        B = new ComplexMatrix(NCOLUMNS, NROWS).dice();
        Bt = new ComplexMatrix(NROWS, NCOLUMNS).dice();
    }
}
