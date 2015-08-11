part of cern.colt.matrix.double.test;

testSparseDoubleMatrix() {
  group('SparseDoubleMatrix', () {
    test('rowCompressed', testRowCompressed);
    test('columnCompressed', testColumnCompressed);
  });
}

void testRowCompressed() {
  int SIZE = NROWS * NCOLUMNS;
  var rowindexes = new Int32List(SIZE);
  var columnindexes = new Int32List(SIZE);
  Float64List values = new Float64List(SIZE);
  for (int i = 0; i < SIZE; i++) {
    rowindexes[i] = (_r.nextInt(MAX_INT) % NROWS).abs();
    columnindexes[i] = (_r.nextInt(MAX_INT) % NCOLUMNS).abs();
    values[i] = _r.nextDouble();
  }
  var A = new SparseDoubleMatrix.withValues(
      NROWS, NCOLUMNS, rowindexes, columnindexes, values);
  SparseRCDoubleMatrix B = A.rowCompressed(false);
  for (int r = 0; r < A.rows; r++) {
    for (int c = 0; c < A.columns; c++) {
      expect(A.get(r, c), equals(B.get(r, c)));
    }
  }
  B = A.rowCompressed(true);
  for (int r = 0; r < A.rows; r++) {
    for (int c = 0; c < A.columns; c++) {
      expect(A.get(r, c), equals(B.get(r, c)));
    }
  }
}

void testColumnCompressed() {
  int SIZE = NROWS * NCOLUMNS;
  var rowindexes = new Int32List(SIZE);
  var columnindexes = new Int32List(SIZE);
  Float64List values = new Float64List(SIZE);
  for (int i = 0; i < SIZE; i++) {
    rowindexes[i] = (_r.nextInt(MAX_INT) % NROWS).abs();
    columnindexes[i] = (_r.nextInt(MAX_INT) % NCOLUMNS).abs();
    values[i] = _r.nextDouble();
  }
  var A = new SparseDoubleMatrix.withValues(
      NROWS, NCOLUMNS, rowindexes, columnindexes, values);
  SparseCCDoubleMatrix B = A.columnCompressed(false);
  for (int r = 0; r < A.rows; r++) {
    for (int c = 0; c < A.columns; c++) {
      expect(A.get(r, c), equals(B.get(r, c)));
    }
  }
  B = A.columnCompressed(true);
  for (int r = 0; r < A.rows; r++) {
    for (int c = 0; c < A.columns; c++) {
      expect(A.get(r, c), equals(B.get(r, c)));
    }
  }
}
