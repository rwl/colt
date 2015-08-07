part of cern.colt.matrix.complex.test;

const NROWS = 13;
const NCOLUMNS = 17;
//const TOL = 1e-10;

testAbstractComplexMatrix(
    String kind, AbstractComplexMatrix make(int rows, int columns)) {
  group("AbstractComplexMatrix ($kind)", () {
    AbstractComplexMatrix A, B, Bt;
    setUp(() {
      A = make(NROWS, NCOLUMNS);
      B = make(NROWS, NCOLUMNS);
      Bt = make(NCOLUMNS, NROWS);

      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          A.setParts(r, c, random.nextDouble(), random.nextDouble());
          B.setParts(r, c, random.nextDouble(), random.nextDouble());
          Bt.setParts(c, r, random.nextDouble(), random.nextDouble());
        }
      }
    });

    test('aggregate', () {
      Float64List actual = A.aggregate(plus, square);
      Float64List expected = new Float64List(2);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expected = cmath.plus(expected, cmath.square(A.get(r, c)));
        }
      }
      assertEquals(expected, actual, TOL);
    });

    test('apply', () {
      AbstractComplexMatrix Acopy = A.copy();
      A.apply(acos);
      Float64List tmp;
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          tmp = cmath.acos(Acopy.get(r, c));
          assertEquals(tmp, A.get(r, c), TOL);
        }
      }
    });

    test('copyFrom', () {
      A.copyFrom(B);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          assertEquals(B.get(r, c), A.get(r, c), TOL);
        }
      }
    });

    test('assign', () {
      AbstractComplexMatrix Acopy = A.copy();
      A.assign(B, div);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          assertEquals(
              cmath.div_(Acopy.get(r, c), B.get(r, c)), A.get(r, c), TOL);
        }
      }
    });

    test('applyReal', () {
      AbstractComplexMatrix Acopy = A.copy();
      A.applyReal(abs);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          var tmp = A.get(r, c);
          expect(cmath.abs(Acopy.get(r, c)), closeTo(tmp[0], TOL));
          expect(0, closeTo(tmp[1], TOL));
        }
      }
    });

    test('setAll', () {
      var expected = new Float64List(2 * A.size);
      for (int i = 0; i < 2 * A.size; i++) {
        expected[i] = random.nextDouble();
      }
      A.setAll(expected);
      int idx = 0;
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          Float64List elem = A.get(r, c);
          expect(expected[idx], closeTo(elem[0], TOL));
          expect(expected[idx + 1], closeTo(elem[1], TOL));
          idx += 2;
        }
      }
    });

    test('fill', () {
      final value = [random.nextDouble(), random.nextDouble()];
      A.fill(value[0], value[1]);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          Float64List elem = A.get(r, c);
          assertEquals(value, elem, TOL);
        }
      }
    });

    test('setImaginary', () {
      AbstractDoubleMatrix Im = new DoubleMatrix.random(A.rows, A.columns);
      AbstractComplexMatrix Acopy = A.copy();
      A.setImaginary(Im);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(Acopy.get(r, c)[0], closeTo(A.get(r, c)[0], TOL));
          expect(Im.get(r, c), closeTo(A.get(r, c)[1], TOL));
        }
      }
    });

    test('setReal', () {
      AbstractDoubleMatrix Re = new DoubleMatrix.random(A.rows, A.columns);
      AbstractComplexMatrix Acopy = A.copy();
      A.setReal(Re);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(Acopy.get(r, c)[1], closeTo(A.get(r, c)[1], TOL));
          expect(Re.get(r, c), closeTo(A.get(r, c)[0], TOL));
        }
      }
    });

    test('cardinality', () {
      int card = A.cardinality;
      expect(A.size, equals(card));
    });

    test('equals', () {
      expect(A.equals(A), isTrue);
      expect(A.equals(B), isFalse);
    });

    test('forEachNonZero', () {
      AbstractComplexMatrix Acopy = A.copy();
      Float64List function(int first, int second, Float64List third) {
        return cmath.sqrt(third);
      }
      A.forEachNonZero(function);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          assertEquals(cmath.sqrt(Acopy.get(r, c)), A.get(r, c), TOL);
        }
      }
    });

    test('conjugateTranspose', () {
      AbstractComplexMatrix Aconj = A.conjugateTranspose();
      expect(A.rows, equals(Aconj.columns));
      expect(A.columns, equals(Aconj.rows));
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(A.get(r, c)[0], closeTo(Aconj.get(c, r)[0], TOL));
          expect(-A.get(r, c)[1], closeTo(Aconj.get(c, r)[1], TOL));
        }
      }
    });

    test('imaginary', () {
      AbstractDoubleMatrix Im = A.imaginary();
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(A.get(r, c)[1], closeTo(Im.get(r, c), TOL));
        }
      }
    });

    test('nonzero', () {
      var rowList = <int>[];
      var colList = <int>[];
      var valueList = <Float64List>[];
      A.nonzero(rowList, colList, valueList);
      expect(A.size, equals(rowList.length));
      expect(A.size, equals(colList.length));
      expect(A.size, equals(valueList.length));
      int idx = 0;
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          assertEquals(A.get(rowList[idx], colList[idx]), valueList[idx], TOL);
          idx++;
        }
      }
    });

    test('real', () {
      AbstractDoubleMatrix Re = A.real();
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(A.get(r, c)[0], closeTo(Re.get(r, c), TOL));
        }
      }
    });

    test('toList', () {
      List<Float64List> array = toList(A);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(A.get(r, c)[0], closeTo(array[r][2 * c], TOL));
          expect(A.get(r, c)[1], closeTo(array[r][2 * c + 1], TOL));
        }
      }
    });

    test('column', () {
      AbstractComplexVector B = A.column(A.columns ~/ 2);
      expect(A.rows, equals(B.size));
      for (int r = 0; r < A.rows; r++) {
        assertEquals(A.get(r, A.columns ~/ 2), B.get(r), TOL);
      }
    });

    test('columnFlip', () {
      AbstractComplexMatrix B = A.columnFlip();
      expect(A.size, equals(B.size));
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          assertEquals(A.get(r, A.columns - 1 - c), B.get(r, c), TOL);
        }
      }
    });

    test('dice', () {
      AbstractComplexMatrix B = A.dice();
      expect(A.rows, equals(B.columns));
      expect(A.columns, equals(B.rows));
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          assertEquals(A.get(r, c), B.get(c, r), TOL);
        }
      }
    });

    test('part', () {
      AbstractComplexMatrix B =
          A.part(A.rows ~/ 2, A.columns ~/ 2, A.rows ~/ 3, A.columns ~/ 3);
      for (int r = 0; r < A.rows / 3; r++) {
        for (int c = 0; c < A.columns / 3; c++) {
          assertEquals(
              A.get(A.rows ~/ 2 + r, A.columns ~/ 2 + c), B.get(r, c), TOL);
        }
      }
    });

    test('row', () {
      AbstractComplexVector B = A.row(A.rows ~/ 2);
      expect(A.columns, equals(B.size));
      for (int c = 0; c < A.columns; c++) {
        assertEquals(A.get(A.rows ~/ 2, c), B.get(c), TOL);
      }
    });

    test('rowFlip', () {
      AbstractComplexMatrix B = A.rowFlip();
      expect(A.size, equals(B.size));
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          assertEquals(A.get(A.rows - 1 - r, c), B.get(r, c), TOL);
        }
      }
    });

    test('select', () {
      final rowIndexes = new Int32List.fromList(
          [A.rows ~/ 6, A.rows ~/ 5, A.rows ~/ 4, A.rows ~/ 3, A.rows ~/ 2]);
      final colIndexes = new Int32List.fromList([
        A.columns ~/ 6,
        A.columns ~/ 5,
        A.columns ~/ 4,
        A.columns ~/ 3,
        A.columns ~/ 2,
        A.columns - 1
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

    test('mult', () {
      AbstractComplexVector y = new ComplexVector(A.columns);
      for (int i = 0; i < y.size; i++) {
        y.set(i, new Float64List.fromList(
            [random.nextDouble(), random.nextDouble()]));
      }
      final alpha = new Float64List.fromList([3.0, 2.0]);
      final beta = new Float64List.fromList([5.0, 4.0]);
      AbstractComplexVector z = null;
      z = A.mult(y, z, alpha, beta, false);
      Float64List expected = new Float64List(2 * A.rows);
      Float64List tmp = new Float64List(2);
      for (int r = 0; r < A.rows; r++) {
        Float64List s = new Float64List(2);
        for (int c = 0; c < A.columns; c++) {
          s = cmath.plus(s, cmath.multiply(A.get(r, c), y.get(c)));
        }
        tmp[0] = expected[2 * r];
        tmp[1] = expected[2 * r + 1];
        tmp = cmath.multiply(beta, tmp);
        tmp = cmath.plus(tmp, cmath.multiply(alpha, s));
        expected[2 * r] = tmp[0];
        expected[2 * r + 1] = tmp[1];
      }

      for (int r = 0; r < A.rows; r++) {
        expect(expected[2 * r], closeTo(z.get(r)[0], TOL));
        expect(expected[2 * r + 1], closeTo(z.get(r)[1], TOL));
      }
      //transpose
      y = new ComplexVector(A.rows);
      for (int i = 0; i < y.size; i++) {
        y.set(i, new Float64List.fromList(
            [random.nextDouble(), random.nextDouble()]));
      }
      z = null;
      z = A.mult(y, z, alpha, beta, true);
      expected = new Float64List(2 * A.columns);
      for (int r = 0; r < A.columns; r++) {
        Float64List s = new Float64List(2);
        for (int c = 0; c < A.rows; c++) {
          s = cmath.plus(s, cmath.multiply(cmath.conj(A.get(c, r)), y.get(c)));
        }
        tmp[0] = expected[2 * r];
        tmp[1] = expected[2 * r + 1];
        tmp = cmath.multiply(beta, tmp);
        tmp = cmath.plus(tmp, cmath.multiply(alpha, s));
        expected[2 * r] = tmp[0];
        expected[2 * r + 1] = tmp[1];
      }
      for (int r = 0; r < A.columns; r++) {
        expect(expected[2 * r], closeTo(z.get(r)[0], TOL));
        expect(expected[2 * r + 1], closeTo(z.get(r)[1], TOL));
      }
    });

    test('multiply', () {
      final alpha = new Float64List.fromList([3.0, 2.0]);
      final beta = new Float64List.fromList([5.0, 4.0]);
      Float64List tmp = new Float64List(2);
      AbstractComplexMatrix C = null;
      C = A.multiply(Bt, C, alpha, beta, false, false);
      List<Float64List> expected = new List<Float64List>.generate(
          A.rows, (_) => new Float64List(2 * A.rows));
      for (int j = 0; j < A.rows; j++) {
        for (int i = 0; i < A.rows; i++) {
          Float64List s = new Float64List(2);
          for (int k = 0; k < A.columns; k++) {
            s = cmath.plus(s, cmath.multiply(A.get(i, k), Bt.get(k, j)));
          }
          tmp[0] = expected[i][2 * j];
          tmp[1] = expected[i][2 * j + 1];
          tmp = cmath.multiply(tmp, beta);
          tmp = cmath.plus(tmp, cmath.multiply(s, alpha));
          expected[i][2 * j] = tmp[0];
          expected[i][2 * j + 1] = tmp[1];
        }
      }
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.rows; c++) {
          expect(expected[r][2 * c], closeTo(C.get(r, c)[0], TOL));
          expect(expected[r][2 * c + 1], closeTo(C.get(r, c)[1], TOL));
        }
      }

      //transposeA
      C = null;
      C = A.multiply(B, C, alpha, beta, true, false);
      expected = new List<Float64List>.generate(
          A.columns, (_) => new Float64List(2 * A.columns));
      for (int j = 0; j < A.columns; j++) {
        for (int i = 0; i < A.columns; i++) {
          Float64List s = new Float64List(2);
          for (int k = 0; k < A.rows; k++) {
            s = cmath.plus(
                s, cmath.multiply(cmath.conj(A.get(k, i)), B.get(k, j)));
          }
          tmp[0] = expected[i][2 * j];
          tmp[1] = expected[i][2 * j + 1];
          tmp = cmath.multiply(tmp, beta);
          tmp = cmath.plus(tmp, cmath.multiply(s, alpha));
          expected[i][2 * j] = tmp[0];
          expected[i][2 * j + 1] = tmp[1];
        }
      }
      for (int r = 0; r < A.columns; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(expected[r][2 * c], closeTo(C.get(r, c)[0], TOL));
          expect(expected[r][2 * c + 1], closeTo(C.get(r, c)[1], TOL));
        }
      }
      //transposeB
      C = null;
      C = A.multiply(B, C, alpha, beta, false, true);
      expected = new List<Float64List>.generate(
          A.rows, (_) => new Float64List(2 * A.rows));
      for (int j = 0; j < A.rows; j++) {
        for (int i = 0; i < A.rows; i++) {
          Float64List s = new Float64List(2);
          for (int k = 0; k < A.columns; k++) {
            s = cmath.plus(
                s, cmath.multiply(A.get(i, k), cmath.conj(B.get(j, k))));
          }
          tmp[0] = expected[i][2 * j];
          tmp[1] = expected[i][2 * j + 1];
          tmp = cmath.multiply(tmp, beta);
          tmp = cmath.plus(tmp, cmath.multiply(s, alpha));
          expected[i][2 * j] = tmp[0];
          expected[i][2 * j + 1] = tmp[1];
        }
      }
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.rows; c++) {
          expect(expected[r][2 * c], closeTo(C.get(r, c)[0], TOL));
          expect(expected[r][2 * c + 1], closeTo(C.get(r, c)[1], TOL));
        }
      }
      //transposeA and transposeB
      C = null;
      C = A.multiply(Bt, C, alpha, beta, true, true);
      expected = new List<Float64List>.generate(
          A.columns, (_) => new Float64List(2 * A.columns));
      for (int j = 0; j < A.columns; j++) {
        for (int i = 0; i < A.columns; i++) {
          Float64List s = new Float64List(2);
          for (int k = 0; k < A.rows; k++) {
            s = cmath.plus(s, cmath.multiply(
                cmath.conj(A.get(k, i)), cmath.conj(Bt.get(j, k))));
          }
          tmp[0] = expected[i][2 * j];
          tmp[1] = expected[i][2 * j + 1];
          tmp = cmath.multiply(tmp, beta);
          tmp = cmath.plus(tmp, cmath.multiply(s, alpha));
          expected[i][2 * j] = tmp[0];
          expected[i][2 * j + 1] = tmp[1];
        }
      }
      for (int r = 0; r < A.columns; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(expected[r][2 * c], closeTo(C.get(r, c)[0], TOL));
          expect(expected[r][2 * c + 1], closeTo(C.get(r, c)[1], TOL));
        }
      }
    });

    test('sum', () {
      Float64List actual = A.sum();
      var expected = new Float64List(2);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expected = cmath.plus(expected, A.get(r, c));
        }
      }
      assertEquals(expected, actual, TOL);
    });
  });
}

//void assertEquals(List<double> expected, List<double> actual, double tol) {
//  for (int i = 0; i < actual.length; i++) {
//    expect(expected[i], closeTo(actual[i], tol));
//  }
//}
