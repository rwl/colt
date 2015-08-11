part of cern.colt.matrix.double.test;

const int DINDEX = 3;

testDiagonalDoubleMatrix(bool view) {
  group("DiagonalDoubleMatrix (${view ? 'view' : 'raw'})", () {
    DoubleMatrix A, B, Bt;
    int dlength;

    setUp(() {
      if (!view) {
        A = new DiagonalDoubleMatrix(NROWS, NCOLUMNS, DINDEX);
        B = new DiagonalDoubleMatrix(NROWS, NCOLUMNS, DINDEX);
        Bt = new DiagonalDoubleMatrix(NCOLUMNS, NROWS, -DINDEX);
        dlength = (A as DiagonalDoubleMatrix).diagonalLength;
      } else {
        A = new DiagonalDoubleMatrix(NCOLUMNS, NROWS, -DINDEX);
        dlength = (A as DiagonalDoubleMatrix).diagonalLength;
        A = A.dice();
        B = new DiagonalDoubleMatrix(NCOLUMNS, NROWS, -DINDEX).dice();
        Bt = new DiagonalDoubleMatrix(NROWS, NCOLUMNS, DINDEX).dice();
      }

      if (DINDEX >= 0) {
        for (int r = 0; r < dlength; r++) {
          A.set(r, r + DINDEX, _r.nextDouble());
          B.set(r, r + DINDEX, _r.nextDouble());
          Bt.set(r - DINDEX, r, _r.nextDouble());
        }
      } else {
        for (int r = 0; r < dlength; r++) {
          A.set(r - DINDEX, r, _r.nextDouble());
          B.set(r - DINDEX, r, _r.nextDouble());
          Bt.set(r, r + DINDEX, _r.nextDouble());
        }
      }
    });

    test('fill', () {
      double value = _r.nextDouble();
      A.fill(value);
      if (DINDEX >= 0) {
        for (int r = 0; r < dlength; r++) {
          expect(value, closeTo(A.get(r, r + DINDEX), TOL));
        }
      } else {
        for (int r = 0; r < dlength; r++) {
          expect(value, closeTo(A.get(r - DINDEX, r), TOL));
        }
      }
    });

    test('setAll', () {
      Float64List expected = new Float64List(dlength);
      for (int i = 0; i < dlength; i++) {
        expected[i] = _r.nextDouble();
      }
      A.setAll(expected);
      if (DINDEX >= 0) {
        for (int r = 0; r < dlength; r++) {
          expect(expected[r], closeTo(A.get(r, r + DINDEX), TOL));
        }
      } else {
        for (int r = 0; r < dlength; r++) {
          expect(expected[r], closeTo(A.get(r - DINDEX, r), TOL));
        }
      }
    });

    test('apply', () {
      DoubleMatrix Acopy = A.copy();
      A.apply(acos);
      if (DINDEX >= 0) {
        for (int r = 0; r < dlength; r++) {
          double expected = math.acos(Acopy.get(r, r + DINDEX));
          expect(expected, closeTo(A.get(r, r + DINDEX), TOL));
        }
      } else {
        for (int r = 0; r < dlength; r++) {
          double expected = math.acos(Acopy.get(r - DINDEX, r));
          expect(expected, closeTo(A.get(r - DINDEX, r), TOL));
        }
      }
    });

    test('assign', () {
      DoubleMatrix Acopy = A.copy();
      A.assign(B, div);
      if (DINDEX >= 0) {
        for (int r = 0; r < dlength; r++) {
          expect(Acopy.get(r, r + DINDEX) / B.get(r, r + DINDEX),
              closeTo(A.get(r, r + DINDEX), TOL));
        }
      } else {
        for (int r = 0; r < dlength; r++) {
          expect(Acopy.get(r - DINDEX, r) / B.get(r - DINDEX, r),
              closeTo(A.get(r - DINDEX, r), TOL));
        }
      }
    });

    test('cardinality', () {
      int card = A.cardinality;
      expect(dlength, equals(card));
    });

    test('max', () {
      A.fill(0.0);
      if (DINDEX >= 0) {
        A.set(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, 0.7);
        A.set(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, 0.1);
        var maxAndLoc = A.max();
        expect(0.7, closeTo(maxAndLoc.value, TOL));
        expect(NROWS ~/ 3, equals(maxAndLoc.row));
        expect(NROWS ~/ 3 + DINDEX, equals(maxAndLoc.column));
      } else {
        A.set(NROWS ~/ 3 - DINDEX, NROWS ~/ 3, 0.7);
        A.set(NROWS ~/ 2 - DINDEX, NROWS ~/ 2, 0.1);
        var maxAndLoc = A.max();
        expect(0.7, closeTo(maxAndLoc.value, TOL));
        expect(NROWS ~/ 3 - DINDEX, equals(maxAndLoc.row));
        expect(NROWS ~/ 3, equals(maxAndLoc.column));
      }
    });

    test('min', () {
      A.fill(0.0);
      if (DINDEX >= 0) {
        A.set(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, -0.7);
        A.set(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, -0.1);
        var minAndLoc = A.min();
        expect(-0.7, closeTo(minAndLoc.value, TOL));
        expect(NROWS ~/ 3, equals(minAndLoc.row));
        expect(NROWS ~/ 3 + DINDEX, equals(minAndLoc.column));
      } else {
        A.set(NROWS ~/ 3 - DINDEX, NROWS ~/ 3, -0.7);
        A.set(NROWS ~/ 2 - DINDEX, NROWS ~/ 2, -0.1);
        var minAndLoc = A.min();
        expect(-0.7, closeTo(minAndLoc.value, TOL));
        expect(NROWS ~/ 3 - DINDEX, equals(minAndLoc.row));
        expect(NROWS ~/ 3, equals(minAndLoc.column));
      }
    });

    test('negative', () {
      A.fill(0.0);
      if (DINDEX >= 0) {
        A.set(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, -0.7);
        A.set(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, -0.1);
        List<int> rowList = new List<int>();
        List<int> columnList = new List<int>();
        List<double> valueList = new List<double>();
        A.negative(rowList, columnList, valueList);
        expect(2, equals(rowList.length));
        expect(2, equals(columnList.length));
        expect(2, equals(valueList.length));
        expect(rowList.contains(NROWS ~/ 3), isTrue);
        expect(rowList.contains(NROWS ~/ 2), isTrue);
        expect(columnList.contains(NROWS ~/ 3 + DINDEX), isTrue);
        expect(columnList.contains(NROWS ~/ 2 + DINDEX), isTrue);
        expect(valueList.contains(-0.7), isTrue);
        expect(valueList.contains(-0.1), isTrue);
      } else {
        A.set(NROWS ~/ 3 - DINDEX, NROWS ~/ 3, -0.7);
        A.set(NROWS ~/ 2 - DINDEX, NROWS ~/ 2, -0.1);
        List<int> rowList = new List<int>();
        List<int> columnList = new List<int>();
        List<double> valueList = new List<double>();
        A.negative(rowList, columnList, valueList);
        expect(2, equals(rowList.length));
        expect(2, equals(columnList.length));
        expect(2, equals(valueList.length));
        expect(rowList.contains(NROWS ~/ 3 - DINDEX), isTrue);
        expect(rowList.contains(NROWS ~/ 2 - DINDEX), isTrue);
        expect(columnList.contains(NROWS ~/ 3), isTrue);
        expect(columnList.contains(NROWS ~/ 2), isTrue);
        expect(valueList.contains(-0.7), isTrue);
        expect(valueList.contains(-0.1), isTrue);
      }
    });

    test('nonzero', () {
      A.fill(0.0);
      if (DINDEX >= 0) {
        A.set(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, 0.7);
        A.set(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, 0.1);
        List<int> rowList = new List<int>();
        List<int> columnList = new List<int>();
        List<double> valueList = new List<double>();
        A.nonzero(rowList, columnList, valueList);
        expect(2, equals(rowList.length));
        expect(2, equals(columnList.length));
        expect(2, equals(valueList.length));
        expect(rowList.contains(NROWS ~/ 3), isTrue);
        expect(rowList.contains(NROWS ~/ 2), isTrue);
        expect(columnList.contains(NROWS ~/ 3 + DINDEX), isTrue);
        expect(columnList.contains(NROWS ~/ 2 + DINDEX), isTrue);
        expect(valueList.contains(0.7), isTrue);
        expect(valueList.contains(0.1), isTrue);
      } else {
        A.set(NROWS ~/ 3 - DINDEX, NROWS ~/ 3, 0.7);
        A.set(NROWS ~/ 2 - DINDEX, NROWS ~/ 2, 0.1);
        List<int> rowList = new List<int>();
        List<int> columnList = new List<int>();
        List<double> valueList = new List<double>();
        A.nonzero(rowList, columnList, valueList);
        expect(2, equals(rowList.length));
        expect(2, equals(columnList.length));
        expect(2, equals(valueList.length));
        expect(rowList.contains(NROWS ~/ 3 - DINDEX), isTrue);
        expect(rowList.contains(NROWS ~/ 2 - DINDEX), isTrue);
        expect(columnList.contains(NROWS ~/ 3), isTrue);
        expect(columnList.contains(NROWS ~/ 2), isTrue);
        expect(valueList.contains(0.7), isTrue);
        expect(valueList.contains(0.1), isTrue);
      }
    });

    test('positive', () {
      A.fill(0.0);
      if (DINDEX >= 0) {
        A.set(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, 0.7);
        A.set(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, 0.1);
        List<int> rowList = new List<int>();
        List<int> columnList = new List<int>();
        List<double> valueList = new List<double>();
        A.positive(rowList, columnList, valueList);
        expect(2, equals(rowList.length));
        expect(2, equals(columnList.length));
        expect(2, equals(valueList.length));
        expect(rowList.contains(NROWS ~/ 3), isTrue);
        expect(rowList.contains(NROWS ~/ 2), isTrue);
        expect(columnList.contains(NROWS ~/ 3 + DINDEX), isTrue);
        expect(columnList.contains(NROWS ~/ 2 + DINDEX), isTrue);
        expect(valueList.contains(0.7), isTrue);
        expect(valueList.contains(0.1), isTrue);
      } else {
        A.set(NROWS ~/ 3 - DINDEX, NROWS ~/ 3, 0.7);
        A.set(NROWS ~/ 2 - DINDEX, NROWS ~/ 2, 0.1);
        List<int> rowList = new List<int>();
        List<int> columnList = new List<int>();
        List<double> valueList = new List<double>();
        A.positive(rowList, columnList, valueList);
        expect(2, equals(rowList.length));
        expect(2, equals(columnList.length));
        expect(2, equals(valueList.length));
        expect(rowList.contains(NROWS ~/ 3 - DINDEX), isTrue);
        expect(rowList.contains(NROWS ~/ 2 - DINDEX), isTrue);
        expect(columnList.contains(NROWS ~/ 3), isTrue);
        expect(columnList.contains(NROWS ~/ 2), isTrue);
        expect(valueList.contains(0.7), isTrue);
        expect(valueList.contains(0.1), isTrue);
      }
    });

    test('column', () {
      DoubleVector col = A.column(NCOLUMNS ~/ 2);
      expect(NROWS, equals(col.size));
      for (int r = 0; r < NROWS; r++) {
        expect(A.get(r, NCOLUMNS ~/ 2), closeTo(col.get(r), TOL));
      }
    });

    test('columnFlip', () {
      DoubleMatrix B = A.columnFlip();
      expect(A.size, equals(B.size));
      for (int r = 0; r < NROWS; r++) {
        for (int c = 0; c < NCOLUMNS; c++) {
          expect(A.get(r, NCOLUMNS - 1 - c), closeTo(B.get(r, c), TOL));
        }
      }
    });

    test('dice', () {
      DoubleMatrix B = A.dice();
      expect(NROWS, equals(B.columns));
      expect(NCOLUMNS, equals(B.rows));
      for (int r = 0; r < NROWS; r++) {
        for (int c = 0; c < NCOLUMNS; c++) {
          expect(A.get(r, c), closeTo(B.get(c, r), TOL));
        }
      }
    });

    test('part', () {
      DoubleMatrix B =
          A.part(NROWS ~/ 2, NCOLUMNS ~/ 2, NROWS ~/ 3, NCOLUMNS ~/ 3);
      expect(NROWS ~/ 3, equals(B.rows));
      expect(NCOLUMNS ~/ 3, equals(B.columns));
      for (int r = 0; r < NROWS ~/ 3; r++) {
        for (int c = 0; c < NCOLUMNS ~/ 3; c++) {
          expect(A.get(NROWS ~/ 2 + r, NCOLUMNS ~/ 2 + c),
              closeTo(B.get(r, c), TOL));
        }
      }
    });

    test('row', () {
      DoubleVector B = A.row(NROWS ~/ 2);
      expect(NCOLUMNS, equals(B.size));
      for (int r = 0; r < NCOLUMNS; r++) {
        expect(A.get(NROWS ~/ 2, r), closeTo(B.get(r), TOL));
      }
    });

    test('rowFlip', () {
      DoubleMatrix B = A.rowFlip();
      expect(A.size, equals(B.size));
      for (int r = 0; r < NROWS; r++) {
        for (int c = 0; c < NCOLUMNS; c++) {
          expect(A.get(NROWS - 1 - r, c), closeTo(B.get(r, c), TOL));
        }
      }
    });

    test('select', () {
      var rowIndexes = new Int32List.fromList(
          [NROWS ~/ 6, NROWS ~/ 5, NROWS ~/ 4, NROWS ~/ 3, NROWS ~/ 2]);
      var colIndexes = new Int32List.fromList([
        NROWS ~/ 6,
        NROWS ~/ 5,
        NROWS ~/ 4,
        NROWS ~/ 3,
        NROWS ~/ 2,
        NROWS - 1
      ]);
      DoubleMatrix B = A.select(rowIndexes, colIndexes);
      expect(rowIndexes.length, equals(B.rows));
      expect(colIndexes.length, equals(B.columns));
      for (int r = 0; r < rowIndexes.length; r++) {
        for (int c = 0; c < colIndexes.length; c++) {
          expect(
              A.get(rowIndexes[r], colIndexes[c]), closeTo(B.get(r, c), TOL));
        }
      }
    });

    test('strides', () {
      int rowStride = 3;
      int colStride = 5;
      DoubleMatrix B = A.strides(rowStride, colStride);
      for (int r = 0; r < B.rows; r++) {
        for (int c = 0; c < B.columns; c++) {
          expect(
              A.get(r * rowStride, c * colStride), closeTo(B.get(r, c), TOL));
        }
      }
    });

    test('multiply', () {
      double alpha = 3.0;
      double beta = 5.0;
      DoubleMatrix C = new DiagonalDoubleMatrix(NROWS, NROWS, 0);
      for (int i = 0; i < dlength; i++) {
        C.set(i, i, _r.nextDouble());
      }
      List<Float64List> expected = toList(C);
      C = A.multiply(Bt, C, alpha, beta, false, false);
      for (int j = 0; j < NROWS; j++) {
        for (int i = 0; i < NROWS; i++) {
          double s = 0.0;
          for (int k = 0; k < NCOLUMNS; k++) {
            s += A.get(i, k) * Bt.get(k, j);
          }
          expected[i][j] = s * alpha + expected[i][j] * beta;
        }
      }
      for (int r = 0; r < NROWS; r++) {
        for (int c = 0; c < NROWS; c++) {
          expect(expected[r][c], closeTo(C.get(r, c), TOL));
        }
      }

      //---
      C = null;
      C = A.multiply(Bt, C, alpha, beta, false, false);
      expected =
          new List<Float64List>.generate(NROWS, (_) => new Float64List(NROWS));
      for (int j = 0; j < NROWS; j++) {
        for (int i = 0; i < NROWS; i++) {
          double s = 0.0;
          for (int k = 0; k < NCOLUMNS; k++) {
            s += A.get(i, k) * Bt.get(k, j);
          }
          expected[i][j] = s * alpha;
        }
      }
      for (int r = 0; r < NROWS; r++) {
        for (int c = 0; c < NROWS; c++) {
          expect(expected[r][c], closeTo(C.get(r, c), TOL));
        }
      }

      //transposeA
      C = new DiagonalDoubleMatrix(NCOLUMNS, NCOLUMNS, 0);
      for (int i = 0; i < dlength; i++) {
        C.set(i, i, _r.nextDouble());
      }
      expected = toList(C);
      C = A.multiply(B, C, alpha, beta, true, false);
      for (int j = 0; j < NCOLUMNS; j++) {
        for (int i = 0; i < NCOLUMNS; i++) {
          double s = 0.0;
          for (int k = 0; k < NROWS; k++) {
            s += A.get(k, i) * B.get(k, j);
          }
          expected[i][j] = s * alpha + expected[i][j] * beta;
        }
      }
      for (int r = 0; r < NCOLUMNS; r++) {
        for (int c = 0; c < NCOLUMNS; c++) {
          expect(expected[r][c], closeTo(C.get(r, c), TOL));
        }
      }
      //---
      C = null;
      C = A.multiply(B, C, alpha, beta, true, false);
      expected = new List<Float64List>.generate(
          NCOLUMNS, (_) => new Float64List(NCOLUMNS));
      for (int j = 0; j < NCOLUMNS; j++) {
        for (int i = 0; i < NCOLUMNS; i++) {
          double s = 0.0;
          for (int k = 0; k < NROWS; k++) {
            s += A.get(k, i) * B.get(k, j);
          }
          expected[i][j] = s * alpha;
        }
      }
      for (int r = 0; r < NCOLUMNS; r++) {
        for (int c = 0; c < NCOLUMNS; c++) {
          expect(expected[r][c], closeTo(C.get(r, c), TOL));
        }
      }

      //transposeB
      C = new DiagonalDoubleMatrix(NROWS, NROWS, 0);
      for (int i = 0; i < dlength; i++) {
        C.set(i, i, _r.nextDouble());
      }
      expected = toList(C);
      C = A.multiply(B, C, alpha, beta, false, true);
      for (int j = 0; j < NROWS; j++) {
        for (int i = 0; i < NROWS; i++) {
          double s = 0.0;
          for (int k = 0; k < NCOLUMNS; k++) {
            s += A.get(i, k) * B.get(j, k);
          }
          expected[i][j] = s * alpha + expected[i][j] * beta;
        }
      }
      for (int r = 0; r < NROWS; r++) {
        for (int c = 0; c < NROWS; c++) {
          expect(expected[r][c], closeTo(C.get(r, c), TOL));
        }
      }
      //---
      C = null;
      C = A.multiply(B, C, alpha, beta, false, true);
      expected =
          new List<Float64List>.generate(NROWS, (_) => new Float64List(NROWS));
      for (int j = 0; j < NROWS; j++) {
        for (int i = 0; i < NROWS; i++) {
          double s = 0.0;
          for (int k = 0; k < NCOLUMNS; k++) {
            s += A.get(i, k) * B.get(j, k);
          }
          expected[i][j] = s * alpha;
        }
      }
      for (int r = 0; r < NROWS; r++) {
        for (int c = 0; c < NROWS; c++) {
          expect(expected[r][c], closeTo(C.get(r, c), TOL));
        }
      }
      //transposeA and transposeB
      C = new DiagonalDoubleMatrix(NCOLUMNS, NCOLUMNS, 0);
      for (int i = 0; i < dlength; i++) {
        C.set(i, i, _r.nextDouble());
      }
      expected = toList(C);
      C = A.multiply(Bt, C, alpha, beta, true, true);
      for (int j = 0; j < NCOLUMNS; j++) {
        for (int i = 0; i < NCOLUMNS; i++) {
          double s = 0.0;
          for (int k = 0; k < NROWS; k++) {
            s += A.get(k, i) * Bt.get(j, k);
          }
          expected[i][j] = s * alpha + expected[i][j] * beta;
        }
      }
      for (int r = 0; r < NCOLUMNS; r++) {
        for (int c = 0; c < NCOLUMNS; c++) {
          expect(expected[r][c], closeTo(C.get(r, c), TOL));
        }
      }
      //---
      C = null;
      C = A.multiply(Bt, C, alpha, beta, true, true);
      expected = new List<Float64List>.generate(
          NCOLUMNS, (_) => new Float64List(NCOLUMNS));
      for (int j = 0; j < NCOLUMNS; j++) {
        for (int i = 0; i < NCOLUMNS; i++) {
          double s = 0.0;
          for (int k = 0; k < NROWS; k++) {
            s += A.get(k, i) * Bt.get(j, k);
          }
          expected[i][j] = s * alpha;
        }
      }
      for (int r = 0; r < NCOLUMNS; r++) {
        for (int c = 0; c < NCOLUMNS; c++) {
          expect(expected[r][c], closeTo(C.get(r, c), TOL));
        }
      }
    });
  });
}
