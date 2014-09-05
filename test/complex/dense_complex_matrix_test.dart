part of cern.colt.matrix.complex.test;

class DenseDComplexMatrix1DTest extends DComplexMatrix1DTest {
    void createMatrices() {
        A = new DenseDComplexMatrix1D(SIZE);
        B = new DenseDComplexMatrix1D(SIZE);
    }
}

class DenseDComplexMatrix1DViewTest extends DenseDComplexMatrix1DTest {
    void createMatrices() {
        A = new DenseDComplexMatrix1D(SIZE).flip();
        B = new DenseDComplexMatrix1D(SIZE).flip();
    }
}

class DenseDComplexMatrix2DTest extends DComplexMatrix2DTest {
    void createMatrices() {
        A = new DenseDComplexMatrix2D(NROWS, NCOLUMNS);
        B = new DenseDComplexMatrix2D(NROWS, NCOLUMNS);
        Bt = new DenseDComplexMatrix2D(NCOLUMNS, NROWS);
    }
}

class DenseDComplexMatrix2DViewTest extends DenseDComplexMatrix2DTest {
    void createMatrices() {
        A = new DenseDComplexMatrix2D(NCOLUMNS, NROWS).dice();
        B = new DenseDComplexMatrix2D(NCOLUMNS, NROWS).dice();
        Bt = new DenseDComplexMatrix2D(NROWS, NCOLUMNS).dice();
    }
}
