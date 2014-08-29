part of cern.colt.matrix.double.test;

class DenseDoubleMatrix1DTest extends DoubleMatrix1DTest {

  void createMatrices() {
    A = new DenseDoubleMatrix1D(SIZE);
    B = new DenseDoubleMatrix1D(SIZE);
  }

}

class DenseDoubleMatrix1DViewTest extends DenseDoubleMatrix1DTest {

  void createMatrices() {
    A = new DenseDoubleMatrix1D(SIZE).viewFlip();
    B = new DenseDoubleMatrix1D(SIZE).viewFlip();
  }

}
