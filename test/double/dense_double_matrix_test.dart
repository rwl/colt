part of cern.colt.matrix.double.test;

testDenseDoubleMatrix(String kind, DoubleMatrix make(int rows, int columns)) {
  group('DoubleMatrix ($kind)', () {
    DoubleMatrix A;
    setUp(() {
      A = make(NROWS, NCOLUMNS);
    });

    test('setValues', () {
      var expected = new Float64List(A.size);
      for (int i = 0; i < A.size; i++) {
        expected[i] = _r.nextDouble();
      }
      A.setValues(expected);
      int idx = 0;
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(expected[idx++], closeTo(A.get(r, c), TOL));
        }
      }
    });
  });
}
