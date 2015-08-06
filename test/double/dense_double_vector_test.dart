part of cern.colt.matrix.double.test;

class DenseDoubleVectorTest extends DoubleVectorTest {
  void createMatrices() {
    A = new DoubleVector(SIZE);
    B = new DoubleVector(SIZE);
  }
}

class DenseDoubleVectorViewTest extends DenseDoubleVectorTest {
  void createMatrices() {
    A = new DoubleVector(SIZE).flip();
    B = new DoubleVector(SIZE).flip();
  }
}
