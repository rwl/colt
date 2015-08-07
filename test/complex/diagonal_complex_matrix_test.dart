part of cern.colt.matrix.complex.test;

const DINDEX = 3;

testDiagonalComplexMatrix(bool view) {
  group('DiagonalComplexMatrix', () {
    AbstractComplexMatrix A, B, Bt;
    int DLENGTH;
    setUp(() {
      if (!view) {
        A = new DiagonalComplexMatrix(NROWS, NCOLUMNS, DINDEX);
        B = new DiagonalComplexMatrix(NROWS, NCOLUMNS, DINDEX);
        Bt = new DiagonalComplexMatrix(NCOLUMNS, NROWS, -DINDEX);
        DLENGTH = (A as DiagonalComplexMatrix).diagonalLength;
      } else {
        A = new DiagonalComplexMatrix(NCOLUMNS, NROWS, -DINDEX);
        DLENGTH = (A as DiagonalComplexMatrix).diagonalLength;
        A = A.dice();
        B = new DiagonalComplexMatrix(NCOLUMNS, NROWS, -DINDEX).dice();
        Bt = new DiagonalComplexMatrix(NROWS, NCOLUMNS, DINDEX).dice();
      }

      if (DINDEX >= 0) {
        for (int r = 0; r < DLENGTH; r++) {
          A.setParts(r, r + DINDEX, random.nextDouble(), random.nextDouble());
          B.setParts(r, r + DINDEX, random.nextDouble(), random.nextDouble());
          Bt.setParts(r - DINDEX, r, random.nextDouble(), random.nextDouble());
        }
      } else {
        for (int r = 0; r < DLENGTH; r++) {
          A.setParts(r - DINDEX, r, random.nextDouble(), random.nextDouble());
          B.setParts(r - DINDEX, r, random.nextDouble(), random.nextDouble());
          Bt.setParts(r, r + DINDEX, random.nextDouble(), random.nextDouble());
        }
      }
    });

    test('fill', () {
      Float64List value =
          new Float64List.fromList([random.nextDouble(), random.nextDouble()]);
      A.fill(value[0], value[1]);
      if (DINDEX >= 0) {
        for (int r = 0; r < DLENGTH; r++) {
          assertEquals(value, A.get(r, r + DINDEX), TOL);
        }
      } else {
        for (int r = 0; r < DLENGTH; r++) {
          assertEquals(value, A.get(r - DINDEX, r), TOL);
        }
      }
    });

    test('setAll', () {
      Float64List expected = new Float64List(2 * DLENGTH);
      for (int i = 0; i < 2 * DLENGTH; i++) {
        expected[i] = random.nextDouble();
      }
      A.setAll(expected);
      if (DINDEX >= 0) {
        for (int r = 0; r < DLENGTH; r++) {
          expect(expected[2 * r], closeTo(A.get(r, r + DINDEX)[0], TOL));
          expect(expected[2 * r + 1], closeTo(A.get(r, r + DINDEX)[1], TOL));
        }
      } else {
        for (int r = 0; r < DLENGTH; r++) {
          expect(expected[2 * r], closeTo(A.get(r - DINDEX, r)[0], TOL));
          expect(expected[2 * r + 1], closeTo(A.get(r - DINDEX, r)[1], TOL));
        }
      }
    });

    test('setImaginary', () {
      AbstractDoubleMatrix Im = new DoubleMatrix.random(A.rows, A.columns);
      AbstractComplexMatrix Acopy = A.copy();
      A.setImaginary(Im);
      if (DINDEX >= 0) {
        for (int r = 0; r < DLENGTH; r++) {
          expect(Acopy.get(r, r + DINDEX)[0],
              closeTo(A.get(r, r + DINDEX)[0], TOL));
          expect(Im.get(r, r + DINDEX), closeTo(A.get(r, r + DINDEX)[1], TOL));
        }
      } else {
        for (int r = 0; r < DLENGTH; r++) {
          expect(Acopy.get(r - DINDEX, r)[0],
              closeTo(A.get(r - DINDEX, r)[0], TOL));
          expect(Im.get(r - DINDEX, r), closeTo(A.get(r - DINDEX, r)[1], TOL));
        }
      }
    });

    test('setReal', () {
      AbstractDoubleMatrix Re = new DoubleMatrix.random(A.rows, A.columns);
      AbstractComplexMatrix Acopy = A.copy();
      A.setReal(Re);
      if (DINDEX >= 0) {
        for (int r = 0; r < DLENGTH; r++) {
          expect(Acopy.get(r, r + DINDEX)[1],
              closeTo(A.get(r, r + DINDEX)[1], TOL));
          expect(Re.get(r, r + DINDEX), closeTo(A.get(r, r + DINDEX)[0], TOL));
        }
      } else {
        for (int r = 0; r < DLENGTH; r++) {
          expect(Acopy.get(r - DINDEX, r)[1],
              closeTo(A.get(r - DINDEX, r)[1], TOL));
          expect(Re.get(r - DINDEX, r), closeTo(A.get(r - DINDEX, r)[0], TOL));
        }
      }
    });

    test('apply', () {
      AbstractComplexMatrix Acopy = A.copy();
      A.apply(acos);
      if (DINDEX >= 0) {
        for (int r = 0; r < DLENGTH; r++) {
          Float64List expected = cmath.acos(Acopy.get(r, r + DINDEX));
          expect(expected[0], closeTo(A.get(r, r + DINDEX)[0], TOL));
          expect(expected[1], closeTo(A.get(r, r + DINDEX)[1], TOL));
        }
      } else {
        for (int r = 0; r < DLENGTH; r++) {
          Float64List expected = cmath.acos(Acopy.get(r - DINDEX, r));
          expect(expected[0], closeTo(A.get(r - DINDEX, r)[0], TOL));
          expect(expected[1], closeTo(A.get(r - DINDEX, r)[1], TOL));
        }
      }
    });

    test('assign', () {
      AbstractComplexMatrix Acopy = A.copy();
      A.assign(B, div);
      if (DINDEX >= 0) {
        for (int r = 0; r < DLENGTH; r++) {
          expect(cmath.div_(Acopy.get(r, r + DINDEX), B.get(r, r + DINDEX))[0],
              closeTo(A.get(r, r + DINDEX)[0], TOL));
          expect(cmath.div_(Acopy.get(r, r + DINDEX), B.get(r, r + DINDEX))[1],
              closeTo(A.get(r, r + DINDEX)[1], TOL));
        }
      } else {
        for (int r = 0; r < DLENGTH; r++) {
          expect(cmath.div_(Acopy.get(r - DINDEX, r), B.get(r - DINDEX, r))[0],
              closeTo(A.get(r - DINDEX, r)[0], TOL));
          expect(cmath.div_(Acopy.get(r - DINDEX, r), B.get(r - DINDEX, r))[1],
              closeTo(A.get(r - DINDEX, r)[1], TOL));
        }
      }
    });

    test('cardinality', () {
      int card = A.cardinality;
      expect(DLENGTH, equals(card));
    });

    test('nonzero', () {
      A.fill(0.0, 0.0);
      Float64List elem1 = new Float64List.fromList([0.7, 0.8]);
      Float64List elem2 = new Float64List.fromList([0.1, 0.2]);
      if (DINDEX >= 0) {
        A.set(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, elem1);
        A.set(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, elem2);
        List<int> rowList = new List<int>();
        List<int> columnList = new List<int>();
        List<Float64List> valueList = new List<Float64List>();
        A.nonzero(rowList, columnList, valueList);
        expect(2, equals(rowList.length));
        expect(2, equals(columnList.length));
        expect(2, equals(valueList.length));
        expect(rowList.contains(NROWS ~/ 3), isTrue);
        expect(rowList.contains(NROWS ~/ 2), isTrue);
        expect(columnList.contains(NROWS ~/ 3 + DINDEX), isTrue);
        expect(columnList.contains(NROWS ~/ 2 + DINDEX), isTrue);
        assertEquals(A.get(rowList[0], columnList[0]), valueList[0], TOL);
        assertEquals(A.get(rowList[1], columnList[1]), valueList[1], TOL);
      } else {
        A.set(NROWS ~/ 3 - DINDEX, NROWS ~/ 3, elem1);
        A.set(NROWS ~/ 2 - DINDEX, NROWS ~/ 2, elem2);
        List<int> rowList = new List<int>();
        List<int> columnList = new List<int>();
        List<Float64List> valueList = new List<Float64List>();
        A.nonzero(rowList, columnList, valueList);
        expect(2, equals(rowList.length));
        expect(2, equals(columnList.length));
        expect(2, equals(valueList.length));
        expect(rowList.contains(NROWS ~/ 3 - DINDEX), isTrue);
        expect(rowList.contains(NROWS ~/ 2 - DINDEX), isTrue);
        expect(columnList.contains(NROWS ~/ 3), isTrue);
        expect(columnList.contains(NROWS ~/ 2), isTrue);
        assertEquals(A.get(rowList[0], columnList[0]), valueList[0], TOL);
        assertEquals(A.get(rowList[1], columnList[1]), valueList[1], TOL);
      }
    });

    test('column', () {
      AbstractComplexVector col = A.column(NCOLUMNS ~/ 2);
      expect(NROWS, equals(col.size));
      for (int r = 0; r < NROWS; r++) {
        assertEquals(A.get(r, NCOLUMNS ~/ 2), col.get(r), TOL);
      }
    });

    test('columnFlip', () {
      AbstractComplexMatrix B = A.columnFlip();
      expect(A.size, equals(B.size));
      for (int r = 0; r < NROWS; r++) {
        for (int c = 0; c < NCOLUMNS; c++) {
          assertEquals(A.get(r, NCOLUMNS - 1 - c), B.get(r, c), TOL);
        }
      }
    });

    test('dice', () {
      AbstractComplexMatrix B = A.dice();
      expect(NROWS, equals(B.columns));
      expect(NCOLUMNS, equals(B.rows));
      for (int r = 0; r < NROWS; r++) {
        for (int c = 0; c < NCOLUMNS; c++) {
          assertEquals(A.get(r, c), B.get(c, r), TOL);
        }
      }
    });

    test('part', () {
      AbstractComplexMatrix B =
          A.part(NROWS ~/ 2, NCOLUMNS ~/ 2, NROWS ~/ 3, NCOLUMNS ~/ 3);
      expect(NROWS ~/ 3, equals(B.rows));
      expect(NCOLUMNS ~/ 3, equals(B.columns));
      for (int r = 0; r < NROWS ~/ 3; r++) {
        for (int c = 0; c < NCOLUMNS ~/ 3; c++) {
          assertEquals(
              A.get(NROWS ~/ 2 + r, NCOLUMNS ~/ 2 + c), B.get(r, c), TOL);
        }
      }
    });

    test('row', () {
      AbstractComplexVector B = A.row(NROWS ~/ 2);
      expect(NCOLUMNS, equals(B.size));
      for (int r = 0; r < NCOLUMNS; r++) {
        assertEquals(A.get(NROWS ~/ 2, r), B.get(r), TOL);
      }
    });

    test('rowFlip', () {
      AbstractComplexMatrix B = A.rowFlip();
      expect(A.size, equals(B.size));
      for (int r = 0; r < NROWS; r++) {
        for (int c = 0; c < NCOLUMNS; c++) {
          assertEquals(A.get(NROWS - 1 - r, c), B.get(r, c), TOL);
        }
      }
    });

    test('select', () {
      Int32List rowIndexes = new Int32List.fromList(
          [NROWS ~/ 6, NROWS ~/ 5, NROWS ~/ 4, NROWS ~/ 3, NROWS ~/ 2]);
      Int32List colIndexes = new Int32List.fromList([
        NROWS ~/ 6,
        NROWS ~/ 5,
        NROWS ~/ 4,
        NROWS ~/ 3,
        NROWS ~/ 2,
        NROWS - 1
      ]);
      AbstractComplexMatrix B = A.select(rowIndexes, colIndexes);
      expect(rowIndexes.length, equals(B.rows));
      expect(colIndexes.length, equals(B.columns));
      for (int r = 0; r < rowIndexes.length; r++) {
        for (int c = 0; c < colIndexes.length; c++) {
          assertEquals(A.get(rowIndexes[r], colIndexes[c]), B.get(r, c), TOL);
        }
      }
    });

    test('strides', () {
      int rowStride = 3;
      int colStride = 5;
      AbstractComplexMatrix B = A.strides(rowStride, colStride);
      for (int r = 0; r < B.rows; r++) {
        for (int c = 0; c < B.columns; c++) {
          assertEquals(A.get(r * rowStride, c * colStride), B.get(r, c), TOL);
        }
      }
    });

    test('multiply', () {
      Float64List alpha = new Float64List.fromList([3.0, 4.0]);
      Float64List beta = new Float64List.fromList([5.0, 6.0]);
      AbstractComplexMatrix C = new DiagonalComplexMatrix(NROWS, NROWS, 0);
      for (int i = 0; i < DLENGTH; i++) {
        C.setParts(i, i, random.nextDouble(), random.nextDouble());
      }
      List<Float64List> expected = toList(C);
      C = A.multiply(Bt, C, alpha, beta, false, false);
      Float64List elem = new Float64List(2);
      for (int j = 0; j < NROWS; j++) {
        for (int i = 0; i < NROWS; i++) {
          Float64List s = new Float64List(2);
          for (int k = 0; k < NCOLUMNS; k++) {
            s = cmath.plus(s, cmath.multiply(A.get(i, k), Bt.get(k, j)));
          }
          elem[0] = expected[i][2 * j];
          elem[1] = expected[i][2 * j + 1];
          elem = cmath.multiply(beta, elem);
          s = cmath.multiply(alpha, s);
          expected[i][2 * j] = s[0] + elem[0];
          expected[i][2 * j + 1] = s[1] + elem[1];
        }
      }
      for (int r = 0; r < NROWS; r++) {
        for (int c = 0; c < NROWS; c++) {
          expect(expected[r][2 * c], closeTo(C.at(r, c)[0], TOL));
          expect(expected[r][2 * c + 1], closeTo(C.at(r, c)[1], TOL));
        }
      }

      //---
      C = null;
      C = A.multiply(Bt, C, alpha, beta, false, false);
      expected = new List<Float64List>.generate(
          NROWS, (_) => new Float64List(2 * NROWS));
      for (int j = 0; j < NROWS; j++) {
        for (int i = 0; i < NROWS; i++) {
          Float64List s = new Float64List(2);
          for (int k = 0; k < NCOLUMNS; k++) {
            s = cmath.plus(s, cmath.multiply(A.get(i, k), Bt.get(k, j)));
          }
          s = cmath.multiply(alpha, s);
          expected[i][2 * j] = s[0];
          expected[i][2 * j + 1] = s[1];
        }
      }
      for (int r = 0; r < NROWS; r++) {
        for (int c = 0; c < NROWS; c++) {
          expect(expected[r][2 * c], closeTo(C.at(r, c)[0], TOL));
          expect(expected[r][2 * c + 1], closeTo(C.at(r, c)[1], TOL));
        }
      }

      //transposeA
      C = new DiagonalComplexMatrix(NCOLUMNS, NCOLUMNS, 0);
      for (int i = 0; i < DLENGTH; i++) {
        C.setParts(i, i, random.nextDouble(), random.nextDouble());
      }
      expected = toList(C);
      C = A.multiply(B, C, alpha, beta, true, false);
      for (int j = 0; j < NCOLUMNS; j++) {
        for (int i = 0; i < NCOLUMNS; i++) {
          Float64List s = new Float64List(2);
          for (int k = 0; k < NROWS; k++) {
            s = cmath.plus(
                s, cmath.multiply(cmath.conj(A.get(k, i)), B.get(k, j)));
          }
          elem[0] = expected[i][2 * j];
          elem[1] = expected[i][2 * j + 1];
          elem = cmath.multiply(beta, elem);
          s = cmath.multiply(alpha, s);
          expected[i][2 * j] = s[0] + elem[0];
          expected[i][2 * j + 1] = s[1] + elem[1];
        }
      }
      for (int r = 0; r < NCOLUMNS; r++) {
        for (int c = 0; c < NCOLUMNS; c++) {
          expect(expected[r][2 * c], closeTo(C.at(r, c)[0], TOL));
          expect(expected[r][2 * c + 1], closeTo(C.at(r, c)[1], TOL));
        }
      }
      //---
      C = null;
      C = A.multiply(B, C, alpha, beta, true, false);
      expected = new List<Float64List>.generate(
          NCOLUMNS, (_) => new Float64List(2 * NCOLUMNS));
      for (int j = 0; j < NCOLUMNS; j++) {
        for (int i = 0; i < NCOLUMNS; i++) {
          Float64List s = new Float64List(2);
          for (int k = 0; k < NROWS; k++) {
            s = cmath.plus(
                s, cmath.multiply(cmath.conj(A.get(k, i)), B.get(k, j)));
          }
          s = cmath.multiply(alpha, s);
          expected[i][2 * j] = s[0];
          expected[i][2 * j + 1] = s[1];
        }
      }
      for (int r = 0; r < NCOLUMNS; r++) {
        for (int c = 0; c < NCOLUMNS; c++) {
          expect(expected[r][2 * c], closeTo(C.at(r, c)[0], TOL));
          expect(expected[r][2 * c + 1], closeTo(C.at(r, c)[1], TOL));
        }
      }

      //transposeB
      C = new DiagonalComplexMatrix(NROWS, NROWS, 0);
      for (int i = 0; i < DLENGTH; i++) {
        C.setParts(i, i, random.nextDouble(), random.nextDouble());
      }
      expected = toList(C);
      C = A.multiply(B, C, alpha, beta, false, true);
      for (int j = 0; j < NROWS; j++) {
        for (int i = 0; i < NROWS; i++) {
          Float64List s = new Float64List(2);
          for (int k = 0; k < NCOLUMNS; k++) {
            s = cmath.plus(
                s, cmath.multiply(A.get(i, k), cmath.conj(B.get(j, k))));
          }
          elem[0] = expected[i][2 * j];
          elem[1] = expected[i][2 * j + 1];
          elem = cmath.multiply(beta, elem);
          s = cmath.multiply(alpha, s);
          expected[i][2 * j] = s[0] + elem[0];
          expected[i][2 * j + 1] = s[1] + elem[1];
        }
      }
      for (int r = 0; r < NROWS; r++) {
        for (int c = 0; c < NROWS; c++) {
          expect(expected[r][2 * c], closeTo(C.at(r, c)[0], TOL));
          expect(expected[r][2 * c + 1], closeTo(C.at(r, c)[1], TOL));
        }
      }
      //---
      C = null;
      C = A.multiply(B, C, alpha, beta, false, true);
      expected = new List<Float64List>.generate(
          NROWS, (_) => new Float64List(2 * NROWS));
      for (int j = 0; j < NROWS; j++) {
        for (int i = 0; i < NROWS; i++) {
          Float64List s = new Float64List(2);
          for (int k = 0; k < NCOLUMNS; k++) {
            s = cmath.plus(
                s, cmath.multiply(A.get(i, k), cmath.conj(B.get(j, k))));
          }
          s = cmath.multiply(alpha, s);
          expected[i][2 * j] = s[0];
          expected[i][2 * j + 1] = s[1];
        }
      }
      for (int r = 0; r < NROWS; r++) {
        for (int c = 0; c < NROWS; c++) {
          expect(expected[r][2 * c], closeTo(C.at(r, c)[0], TOL));
          expect(expected[r][2 * c + 1], closeTo(C.at(r, c)[1], TOL));
        }
      }
      //transposeA and transposeB
      C = new DiagonalComplexMatrix(NCOLUMNS, NCOLUMNS, 0);
      for (int i = 0; i < DLENGTH; i++) {
        C.setParts(i, i, random.nextDouble(), random.nextDouble());
      }
      expected = toList(C);
      C = A.multiply(Bt, C, alpha, beta, true, true);
      for (int j = 0; j < NCOLUMNS; j++) {
        for (int i = 0; i < NCOLUMNS; i++) {
          Float64List s = new Float64List(2);
          for (int k = 0; k < NROWS; k++) {
            s = cmath.plus(s, cmath.multiply(
                cmath.conj(A.get(k, i)), cmath.conj(Bt.get(j, k))));
          }
          elem[0] = expected[i][2 * j];
          elem[1] = expected[i][2 * j + 1];
          elem = cmath.multiply(beta, elem);
          s = cmath.multiply(alpha, s);
          expected[i][2 * j] = s[0] + elem[0];
          expected[i][2 * j + 1] = s[1] + elem[1];
        }
      }
      for (int r = 0; r < NCOLUMNS; r++) {
        for (int c = 0; c < NCOLUMNS; c++) {
          expect(expected[r][2 * c], closeTo(C.at(r, c)[0], TOL));
          expect(expected[r][2 * c + 1], closeTo(C.at(r, c)[1], TOL));
        }
      }
      //---
      C = null;
      C = A.multiply(Bt, C, alpha, beta, true, true);
      expected = new List<Float64List>.generate(
          NCOLUMNS, (_) => new Float64List(2 * NCOLUMNS));
      for (int j = 0; j < NCOLUMNS; j++) {
        for (int i = 0; i < NCOLUMNS; i++) {
          Float64List s = new Float64List(2);
          for (int k = 0; k < NROWS; k++) {
            s = cmath.plus(s, cmath.multiply(A.get(k, i), Bt.get(j, k)));
          }
          s = cmath.multiply(alpha, s);
          expected[i][2 * j] = s[0];
          expected[i][2 * j + 1] = s[1];
        }
      }
      for (int r = 0; r < NCOLUMNS; r++) {
        for (int c = 0; c < NCOLUMNS; c++) {
          expect(expected[r][2 * c], closeTo(C.at(r, c)[0], TOL));
          expect(expected[r][2 * c + 1], closeTo(C.at(r, c)[1], TOL));
        }
      }
    });
  });
}
