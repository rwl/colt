part of cern.colt.matrix.double.test;

class DenseDoubleVectorTest extends DoubleVectorTest {

  void createMatrices() {
    A = new DenseDoubleVector(SIZE);
    B = new DenseDoubleVector(SIZE);
  }

}

class DenseDoubleVectorViewTest extends DenseDoubleVectorTest {

  void createMatrices() {
    A = new DenseDoubleVector(SIZE).flip();
    B = new DenseDoubleVector(SIZE).flip();
  }

}
