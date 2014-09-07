part of cern.colt.matrix.int.test;

class DenseIntVectorTest extends IntVectorTest {

  void createMatrices() {
    A = new DenseIntVector(SIZE);
    B = new DenseIntVector(SIZE);
  }
}

class DenseIntVectorViewTest extends DenseIntVectorTest {

  void createMatrices() {
    A = new DenseIntVector(SIZE).flip();
    B = new DenseIntVector(SIZE).flip();
  }
}
