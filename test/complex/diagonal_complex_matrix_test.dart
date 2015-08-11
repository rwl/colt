part of cern.colt.matrix.complex.test;

const DINDEX = 3;

testDiagonalComplexMatrix(bool view) {
  group('DiagonalComplexMatrix', () {
    ComplexMatrix A, B, Bt;
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
      var value = new Complex(random.nextDouble(), random.nextDouble());
      A.fill(value.real, value.imaginary);
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
          expect(expected[2 * r], closeTo(A.get(r, r + DINDEX).real, TOL));
          expect(expected[2 * r + 1],
              closeTo(A.get(r, r + DINDEX).imaginary, TOL));
        }
      } else {
        for (int r = 0; r < DLENGTH; r++) {
          expect(expected[2 * r], closeTo(A.get(r - DINDEX, r).real, TOL));
          expect(expected[2 * r + 1],
              closeTo(A.get(r - DINDEX, r).imaginary, TOL));
        }
      }
    });

    test('setImaginary', () {
      DoubleMatrix Im = new DenseDoubleMatrix.random(A.rows, A.columns);
      ComplexMatrix Acopy = A.copy();
      A.setImaginary(Im);
      if (DINDEX >= 0) {
        for (int r = 0; r < DLENGTH; r++) {
          expect(Acopy.get(r, r + DINDEX).real,
              closeTo(A.get(r, r + DINDEX).real, TOL));
          expect(Im.get(r, r + DINDEX),
              closeTo(A.get(r, r + DINDEX).imaginary, TOL));
        }
      } else {
        for (int r = 0; r < DLENGTH; r++) {
          expect(Acopy.get(r - DINDEX, r).real,
              closeTo(A.get(r - DINDEX, r).real, TOL));
          expect(Im.get(r - DINDEX, r),
              closeTo(A.get(r - DINDEX, r).imaginary, TOL));
        }
      }
    });

    test('setReal', () {
      DoubleMatrix Re = new DenseDoubleMatrix.random(A.rows, A.columns);
      ComplexMatrix Acopy = A.copy();
      A.setReal(Re);
      if (DINDEX >= 0) {
        for (int r = 0; r < DLENGTH; r++) {
          expect(Acopy.get(r, r + DINDEX).imaginary,
              closeTo(A.get(r, r + DINDEX).imaginary, TOL));
          expect(
              Re.get(r, r + DINDEX), closeTo(A.get(r, r + DINDEX).real, TOL));
        }
      } else {
        for (int r = 0; r < DLENGTH; r++) {
          expect(Acopy.get(r - DINDEX, r).imaginary,
              closeTo(A.get(r - DINDEX, r).imaginary, TOL));
          expect(
              Re.get(r - DINDEX, r), closeTo(A.get(r - DINDEX, r).real, TOL));
        }
      }
    });

    test('apply', () {
      ComplexMatrix Acopy = A.copy();
      A.apply(acos);
      if (DINDEX >= 0) {
        for (int r = 0; r < DLENGTH; r++) {
          var expected = Acopy.get(r, r + DINDEX).acos();
          expect(expected.real, closeTo(A.get(r, r + DINDEX).real, TOL));
          expect(
              expected.imaginary, closeTo(A.get(r, r + DINDEX).imaginary, TOL));
        }
      } else {
        for (int r = 0; r < DLENGTH; r++) {
          var expected = Acopy.get(r - DINDEX, r).acos();
          expect(expected.real, closeTo(A.get(r - DINDEX, r).real, TOL));
          expect(
              expected.imaginary, closeTo(A.get(r - DINDEX, r).imaginary, TOL));
        }
      }
    });

    test('assign', () {
      ComplexMatrix Acopy = A.copy();
      A.assign(B, div);
      if (DINDEX >= 0) {
        for (int r = 0; r < DLENGTH; r++) {
          expect((Acopy.get(r, r + DINDEX) / B.get(r, r + DINDEX)).real,
              closeTo(A.get(r, r + DINDEX).real, TOL));
          expect((Acopy.get(r, r + DINDEX) / B.get(r, r + DINDEX)).imaginary,
              closeTo(A.get(r, r + DINDEX).imaginary, TOL));
        }
      } else {
        for (int r = 0; r < DLENGTH; r++) {
          expect((Acopy.get(r - DINDEX, r) / B.get(r - DINDEX, r)).real,
              closeTo(A.get(r - DINDEX, r).real, TOL));
          expect((Acopy.get(r - DINDEX, r) / B.get(r - DINDEX, r)).imaginary,
              closeTo(A.get(r - DINDEX, r).imaginary, TOL));
        }
      }
    });

    test('cardinality', () {
      int card = A.cardinality;
      expect(DLENGTH, equals(card));
    });

    test('nonzero', () {
      A.fill(0.0, 0.0);
      var elem1 = new Complex(0.7, 0.8);
      var elem2 = new Complex(0.1, 0.2);
      if (DINDEX >= 0) {
        A.set(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, elem1);
        A.set(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, elem2);
        var rowList = new List<int>();
        var columnList = new List<int>();
        var valueList = new List<Complex>();
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
        var rowList = new List<int>();
        var columnList = new List<int>();
        var valueList = new List<Complex>();
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
      ComplexVector col = A.column(NCOLUMNS ~/ 2);
      expect(NROWS, equals(col.size));
      for (int r = 0; r < NROWS; r++) {
        assertEquals(A.get(r, NCOLUMNS ~/ 2), col.get(r), TOL);
      }
    });

    test('columnFlip', () {
      ComplexMatrix B = A.columnFlip();
      expect(A.size, equals(B.size));
      for (int r = 0; r < NROWS; r++) {
        for (int c = 0; c < NCOLUMNS; c++) {
          assertEquals(A.get(r, NCOLUMNS - 1 - c), B.get(r, c), TOL);
        }
      }
    });

    test('dice', () {
      ComplexMatrix B = A.dice();
      expect(NROWS, equals(B.columns));
      expect(NCOLUMNS, equals(B.rows));
      for (int r = 0; r < NROWS; r++) {
        for (int c = 0; c < NCOLUMNS; c++) {
          assertEquals(A.get(r, c), B.get(c, r), TOL);
        }
      }
    });

    test('part', () {
      ComplexMatrix B =
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
      ComplexVector B = A.row(NROWS ~/ 2);
      expect(NCOLUMNS, equals(B.size));
      for (int r = 0; r < NCOLUMNS; r++) {
        assertEquals(A.get(NROWS ~/ 2, r), B.get(r), TOL);
      }
    });

    test('rowFlip', () {
      ComplexMatrix B = A.rowFlip();
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
      ComplexMatrix B = A.select(rowIndexes, colIndexes);
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
      ComplexMatrix B = A.strides(rowStride, colStride);
      for (int r = 0; r < B.rows; r++) {
        for (int c = 0; c < B.columns; c++) {
          assertEquals(A.get(r * rowStride, c * colStride), B.get(r, c), TOL);
        }
      }
    });

    test('multiply', () {
      var alpha = new Complex(3.0, 4.0);
      var beta = new Complex(5.0, 6.0);
      ComplexMatrix C = new DiagonalComplexMatrix(NROWS, NROWS, 0);
      for (int i = 0; i < DLENGTH; i++) {
        C.setParts(i, i, random.nextDouble(), random.nextDouble());
      }
      List<Float64List> expected = toList(C);
      C = A.multiply(Bt, C, alpha, beta, false, false);
      for (int j = 0; j < NROWS; j++) {
        for (int i = 0; i < NROWS; i++) {
          var s = Complex.ZERO;
          for (int k = 0; k < NCOLUMNS; k++) {
            s += A.get(i, k) * Bt.get(k, j);
          }
          var elem = new Complex(expected[i][2 * j], expected[i][2 * j + 1]);
          elem = beta * elem;
          s = alpha * s;
          expected[i][2 * j] = s.real + elem.real;
          expected[i][2 * j + 1] = s.imaginary + elem.imaginary;
        }
      }
      for (int r = 0; r < NROWS; r++) {
        for (int c = 0; c < NROWS; c++) {
          expect(expected[r][2 * c], closeTo(C.at(r, c).real, TOL));
          expect(expected[r][2 * c + 1], closeTo(C.at(r, c).imaginary, TOL));
        }
      }

      //---
      C = null;
      C = A.multiply(Bt, C, alpha, beta, false, false);
      expected = new List<Float64List>.generate(
          NROWS, (_) => new Float64List(2 * NROWS));
      for (int j = 0; j < NROWS; j++) {
        for (int i = 0; i < NROWS; i++) {
          var s = Complex.ZERO;
          for (int k = 0; k < NCOLUMNS; k++) {
            s += A.get(i, k) * Bt.get(k, j);
          }
          s = alpha * s;
          expected[i][2 * j] = s.real;
          expected[i][2 * j + 1] = s.imaginary;
        }
      }
      for (int r = 0; r < NROWS; r++) {
        for (int c = 0; c < NROWS; c++) {
          expect(expected[r][2 * c], closeTo(C.at(r, c).real, TOL));
          expect(expected[r][2 * c + 1], closeTo(C.at(r, c).imaginary, TOL));
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
          var s = Complex.ZERO;
          for (int k = 0; k < NROWS; k++) {
            s += A.get(k, i).conjugate() * B.get(k, j);
          }
          var elem = new Complex(expected[i][2 * j], expected[i][2 * j + 1]);
          elem = beta * elem;
          s = alpha * s;
          expected[i][2 * j] = s.real + elem.real;
          expected[i][2 * j + 1] = s.imaginary + elem.imaginary;
        }
      }
      for (int r = 0; r < NCOLUMNS; r++) {
        for (int c = 0; c < NCOLUMNS; c++) {
          expect(expected[r][2 * c], closeTo(C.at(r, c).real, TOL));
          expect(expected[r][2 * c + 1], closeTo(C.at(r, c).imaginary, TOL));
        }
      }
      //---
      C = null;
      C = A.multiply(B, C, alpha, beta, true, false);
      expected = new List<Float64List>.generate(
          NCOLUMNS, (_) => new Float64List(2 * NCOLUMNS));
      for (int j = 0; j < NCOLUMNS; j++) {
        for (int i = 0; i < NCOLUMNS; i++) {
          var s = Complex.ZERO;
          for (int k = 0; k < NROWS; k++) {
            s += A.get(k, i).conjugate() * B.get(k, j);
          }
          s = alpha * s;
          expected[i][2 * j] = s.real;
          expected[i][2 * j + 1] = s.imaginary;
        }
      }
      for (int r = 0; r < NCOLUMNS; r++) {
        for (int c = 0; c < NCOLUMNS; c++) {
          expect(expected[r][2 * c], closeTo(C.at(r, c).real, TOL));
          expect(expected[r][2 * c + 1], closeTo(C.at(r, c).imaginary, TOL));
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
          var s = Complex.ZERO;
          for (int k = 0; k < NCOLUMNS; k++) {
            s += A.get(i, k) * B.get(j, k).conjugate();
          }
          var elem = new Complex(expected[i][2 * j], expected[i][2 * j + 1]);
          elem = beta * elem;
          s = alpha * s;
          expected[i][2 * j] = s.real + elem.real;
          expected[i][2 * j + 1] = s.imaginary + elem.imaginary;
        }
      }
      for (int r = 0; r < NROWS; r++) {
        for (int c = 0; c < NROWS; c++) {
          expect(expected[r][2 * c], closeTo(C.at(r, c).real, TOL));
          expect(expected[r][2 * c + 1], closeTo(C.at(r, c).imaginary, TOL));
        }
      }
      //---
      C = null;
      C = A.multiply(B, C, alpha, beta, false, true);
      expected = new List<Float64List>.generate(
          NROWS, (_) => new Float64List(2 * NROWS));
      for (int j = 0; j < NROWS; j++) {
        for (int i = 0; i < NROWS; i++) {
          var s = Complex.ZERO;
          for (int k = 0; k < NCOLUMNS; k++) {
            s += A.get(i, k) * B.get(j, k).conjugate();
          }
          s = alpha * s;
          expected[i][2 * j] = s.real;
          expected[i][2 * j + 1] = s.imaginary;
        }
      }
      for (int r = 0; r < NROWS; r++) {
        for (int c = 0; c < NROWS; c++) {
          expect(expected[r][2 * c], closeTo(C.at(r, c).real, TOL));
          expect(expected[r][2 * c + 1], closeTo(C.at(r, c).imaginary, TOL));
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
          var s = Complex.ZERO;
          for (int k = 0; k < NROWS; k++) {
            s += A.get(k, i).conjugate() * Bt.get(j, k).conjugate();
          }
          var elem = new Complex(expected[i][2 * j], expected[i][2 * j + 1]);
          elem = beta * elem;
          s = alpha * s;
          expected[i][2 * j] = s.real + elem.real;
          expected[i][2 * j + 1] = s.imaginary + elem.imaginary;
        }
      }
      for (int r = 0; r < NCOLUMNS; r++) {
        for (int c = 0; c < NCOLUMNS; c++) {
          expect(expected[r][2 * c], closeTo(C.at(r, c).real, TOL));
          expect(expected[r][2 * c + 1], closeTo(C.at(r, c).imaginary, TOL));
        }
      }
      //---
      C = null;
      C = A.multiply(Bt, C, alpha, beta, true, true);
      expected = new List<Float64List>.generate(
          NCOLUMNS, (_) => new Float64List(2 * NCOLUMNS));
      for (int j = 0; j < NCOLUMNS; j++) {
        for (int i = 0; i < NCOLUMNS; i++) {
          var s = Complex.ZERO;
          for (int k = 0; k < NROWS; k++) {
            s += A.get(k, i) * Bt.get(j, k);
          }
          s = alpha * s;
          expected[i][2 * j] = s.real;
          expected[i][2 * j + 1] = s.imaginary;
        }
      }
      for (int r = 0; r < NCOLUMNS; r++) {
        for (int c = 0; c < NCOLUMNS; c++) {
          expect(expected[r][2 * c], closeTo(C.at(r, c).real, TOL));
          expect(expected[r][2 * c + 1], closeTo(C.at(r, c).imaginary, TOL));
        }
      }
    });
  });
}
