part of cern.colt.matrix.int.test;

class DenseIntMatrix1DTest extends IntMatrix1DTest {

  void createMatrices() {
    A = new DenseIntMatrix1D(SIZE);
    B = new DenseIntMatrix1D(SIZE);
  }
}

class DenseIntMatrix1DViewTest extends DenseIntMatrix1DTest {

  void createMatrices() {
    A = new DenseIntMatrix1D(SIZE).flip();
    B = new DenseIntMatrix1D(SIZE).flip();
  }
}
