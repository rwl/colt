part of cern.colt.matrix.int.test;

testDenseIntMatrix(String kind, IntMatrix make(int rows, int columns)) {
  group('IntMatrix ($kind)', () {
    IntMatrix A;
    setUp(() {
      A = make(NROWS, NCOLUMNS);
    });
    test('setAll', () {
      final expected = new Int32List(A.size);
      for (int i = 0; i < A.size; i++) {
        expected[i] = random.nextInt(MAX_INT);
      }
      A.setAll(expected);
      int idx = 0;
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(0, equals((expected[idx++] - A.get(r, c)).abs()));
        }
      }
    });
  });
}
