part of cern.colt.matrix.int.test;

class DenseIntVectorTest extends IntVectorTest {

  void createMatrices() {
    A = new IntVector(SIZE);
    B = new IntVector(SIZE);
  }
}

class DenseIntVectorViewTest extends DenseIntVectorTest {

  void createMatrices() {
    A = new IntVector(SIZE).flip();
    B = new IntVector(SIZE).flip();
  }
}
