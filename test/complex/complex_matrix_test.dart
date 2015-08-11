part of cern.colt.matrix.complex.test;

const NROWS = 13;
const NCOLUMNS = 17;
//const TOL = 1e-10;

testComplexMatrix(String kind, ComplexMatrix make(int rows, int columns)) {
  group("ComplexMatrix ($kind)", () {
    ComplexMatrix A, B, Bt;
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
      var actual = A.aggregate(plus, square);
      var expected = Complex.ZERO;
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expected += A.get(r, c).pow(2);
        }
      }
      assertEquals(expected, actual, TOL);
    });

    test('apply', () {
      ComplexMatrix Acopy = A.copy();
      A.apply(acos);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          var tmp = Acopy.get(r, c).acos();
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
      ComplexMatrix Acopy = A.copy();
      A.assign(B, div);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          assertEquals(Acopy.get(r, c) / B.get(r, c), A.get(r, c), TOL);
        }
      }
    });

    test('applyReal', () {
      ComplexMatrix Acopy = A.copy();
      A.applyReal(abs);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          var tmp = A.get(r, c);
          expect(Acopy.get(r, c).abs(), closeTo(tmp.real, TOL));
          expect(0, closeTo(tmp.imaginary, TOL));
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
          var elem = A.get(r, c);
          expect(expected[idx], closeTo(elem.real, TOL));
          expect(expected[idx + 1], closeTo(elem.imaginary, TOL));
          idx += 2;
        }
      }
    });

    test('fill', () {
      var value = new Complex(random.nextDouble(), random.nextDouble());
      A.fill(value.real, value.imaginary);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          var elem = A.get(r, c);
          assertEquals(value, elem, TOL);
        }
      }
    });

    test('setImaginary', () {
      DoubleMatrix Im = new DenseDoubleMatrix.random(A.rows, A.columns);
      ComplexMatrix Acopy = A.copy();
      A.setImaginary(Im);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(Acopy.get(r, c).real, closeTo(A.get(r, c).real, TOL));
          expect(Im.get(r, c), closeTo(A.get(r, c).imaginary, TOL));
        }
      }
    });

    test('setReal', () {
      DoubleMatrix Re = new DenseDoubleMatrix.random(A.rows, A.columns);
      ComplexMatrix Acopy = A.copy();
      A.setReal(Re);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(
              Acopy.get(r, c).imaginary, closeTo(A.get(r, c).imaginary, TOL));
          expect(Re.get(r, c), closeTo(A.get(r, c).real, TOL));
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
      ComplexMatrix Acopy = A.copy();
      Complex fn(int first, int second, Complex third) {
        return third.sqrt();
      }
      A.forEachNonZero(fn);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          assertEquals(Acopy.get(r, c).sqrt(), A.get(r, c), TOL);
        }
      }
    });

    test('conjugateTranspose', () {
      ComplexMatrix Aconj = A.conjugateTranspose();
      expect(A.rows, equals(Aconj.columns));
      expect(A.columns, equals(Aconj.rows));
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(A.get(r, c).real, closeTo(Aconj.get(c, r).real, TOL));
          expect(
              -A.get(r, c).imaginary, closeTo(Aconj.get(c, r).imaginary, TOL));
        }
      }
    });

    test('imaginary', () {
      DoubleMatrix Im = A.imaginary();
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(A.get(r, c).imaginary, closeTo(Im.get(r, c), TOL));
        }
      }
    });

    test('nonzero', () {
      var rowList = <int>[];
      var colList = <int>[];
      var valueList = <Complex>[];
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
      DoubleMatrix Re = A.real();
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(A.get(r, c).real, closeTo(Re.get(r, c), TOL));
        }
      }
    });

    test('column', () {
      ComplexVector B = A.column(A.columns ~/ 2);
      expect(A.rows, equals(B.size));
      for (int r = 0; r < A.rows; r++) {
        assertEquals(A.get(r, A.columns ~/ 2), B.get(r), TOL);
      }
    });

    test('columnFlip', () {
      ComplexMatrix B = A.columnFlip();
      expect(A.size, equals(B.size));
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          assertEquals(A.get(r, A.columns - 1 - c), B.get(r, c), TOL);
        }
      }
    });

    test('dice', () {
      ComplexMatrix B = A.dice();
      expect(A.rows, equals(B.columns));
      expect(A.columns, equals(B.rows));
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          assertEquals(A.get(r, c), B.get(c, r), TOL);
        }
      }
    });

    test('part', () {
      ComplexMatrix B =
          A.part(A.rows ~/ 2, A.columns ~/ 2, A.rows ~/ 3, A.columns ~/ 3);
      for (int r = 0; r < A.rows / 3; r++) {
        for (int c = 0; c < A.columns / 3; c++) {
          assertEquals(
              A.get(A.rows ~/ 2 + r, A.columns ~/ 2 + c), B.get(r, c), TOL);
        }
      }
    });

    test('row', () {
      ComplexVector B = A.row(A.rows ~/ 2);
      expect(A.columns, equals(B.size));
      for (int c = 0; c < A.columns; c++) {
        assertEquals(A.get(A.rows ~/ 2, c), B.get(c), TOL);
      }
    });

    test('rowFlip', () {
      ComplexMatrix B = A.rowFlip();
      expect(A.size, equals(B.size));
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          assertEquals(A.get(A.rows - 1 - r, c), B.get(r, c), TOL);
        }
      }
    });

    test('select', () {
      var rowIndexes = new Int32List.fromList(
          [A.rows ~/ 6, A.rows ~/ 5, A.rows ~/ 4, A.rows ~/ 3, A.rows ~/ 2]);
      var colIndexes = new Int32List.fromList([
        A.columns ~/ 6,
        A.columns ~/ 5,
        A.columns ~/ 4,
        A.columns ~/ 3,
        A.columns ~/ 2,
        A.columns - 1
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

    test('mult', () {
      ComplexVector y = new DenseComplexVector(A.columns);
      for (int i = 0; i < y.size; i++) {
        y.set(i, new Complex(random.nextDouble(), random.nextDouble()));
      }
      var alpha = new Complex(3.0, 2.0);
      var beta = new Complex(5.0, 4.0);
      ComplexVector z = null;
      z = A.mult(y, z, alpha, beta, false);
      var expected = new Float64List(2 * A.rows);
      for (int r = 0; r < A.rows; r++) {
        var s = Complex.ZERO;
        for (int c = 0; c < A.columns; c++) {
          s = s + (A.get(r, c) * y.get(c));
        }
        var tmp = new Complex(expected[2 * r], expected[2 * r + 1]);
        tmp = beta * tmp;
        tmp = tmp + (alpha * s);
        expected[2 * r] = tmp.real;
        expected[2 * r + 1] = tmp.imaginary;
      }

      for (int r = 0; r < A.rows; r++) {
        expect(expected[2 * r], closeTo(z.get(r).real, TOL));
        expect(expected[2 * r + 1], closeTo(z.get(r).imaginary, TOL));
      }
      //transpose
      y = new DenseComplexVector(A.rows);
      for (int i = 0; i < y.size; i++) {
        y.set(i, new Complex(random.nextDouble(), random.nextDouble()));
      }
      z = null;
      z = A.mult(y, z, alpha, beta, true);
      expected = new Float64List(2 * A.columns);
      for (int r = 0; r < A.columns; r++) {
        var s = Complex.ZERO;
        for (int c = 0; c < A.rows; c++) {
          s += A.get(c, r).conjugate() * y.get(c);
        }
        var tmp = new Complex(expected[2 * r], expected[2 * r + 1]);
        tmp = beta * tmp;
        tmp = tmp + (alpha * s);
        expected[2 * r] = tmp.real;
        expected[2 * r + 1] = tmp.imaginary;
      }
      for (int r = 0; r < A.columns; r++) {
        expect(expected[2 * r], closeTo(z.get(r).real, TOL));
        expect(expected[2 * r + 1], closeTo(z.get(r).imaginary, TOL));
      }
    });

    test('multiply', () {
      var alpha = new Complex(3.0, 2.0);
      var beta = new Complex(5.0, 4.0);
      ComplexMatrix C = null;
      C = A.multiply(Bt, C, alpha, beta, false, false);
      var expected =
          new List.generate(A.rows, (_) => new Float64List(2 * A.rows));
      for (int j = 0; j < A.rows; j++) {
        for (int i = 0; i < A.rows; i++) {
          var s = Complex.ZERO;
          for (int k = 0; k < A.columns; k++) {
            s += A.get(i, k) * Bt.get(k, j);
          }
          var tmp = new Complex(expected[i][2 * j], expected[i][2 * j + 1]);
          tmp = tmp * beta;
          tmp += s * alpha;
          expected[i][2 * j] = tmp.real;
          expected[i][2 * j + 1] = tmp.imaginary;
        }
      }
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.rows; c++) {
          expect(expected[r][2 * c], closeTo(C.get(r, c).real, TOL));
          expect(expected[r][2 * c + 1], closeTo(C.get(r, c).imaginary, TOL));
        }
      }

      //transposeA
      C = null;
      C = A.multiply(B, C, alpha, beta, true, false);
      expected =
          new List.generate(A.columns, (_) => new Float64List(2 * A.columns));
      for (int j = 0; j < A.columns; j++) {
        for (int i = 0; i < A.columns; i++) {
          var s = Complex.ZERO;
          for (int k = 0; k < A.rows; k++) {
            s += A.get(k, i).conjugate() * B.get(k, j);
          }
          var tmp = new Complex(expected[i][2 * j], expected[i][2 * j + 1]);
          tmp = tmp * beta;
          tmp += s * alpha;
          expected[i][2 * j] = tmp.real;
          expected[i][2 * j + 1] = tmp.imaginary;
        }
      }
      for (int r = 0; r < A.columns; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(expected[r][2 * c], closeTo(C.get(r, c).real, TOL));
          expect(expected[r][2 * c + 1], closeTo(C.get(r, c).imaginary, TOL));
        }
      }
      //transposeB
      C = null;
      C = A.multiply(B, C, alpha, beta, false, true);
      expected = new List<Float64List>.generate(
          A.rows, (_) => new Float64List(2 * A.rows));
      for (int j = 0; j < A.rows; j++) {
        for (int i = 0; i < A.rows; i++) {
          var s = Complex.ZERO;
          for (int k = 0; k < A.columns; k++) {
            s += A.get(i, k) * B.get(j, k).conjugate();
          }
          var tmp = new Complex(expected[i][2 * j], expected[i][2 * j + 1]);
          tmp = tmp * beta;
          tmp += s * alpha;
          expected[i][2 * j] = tmp.real;
          expected[i][2 * j + 1] = tmp.imaginary;
        }
      }
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.rows; c++) {
          expect(expected[r][2 * c], closeTo(C.get(r, c).real, TOL));
          expect(expected[r][2 * c + 1], closeTo(C.get(r, c).imaginary, TOL));
        }
      }
      //transposeA and transposeB
      C = null;
      C = A.multiply(Bt, C, alpha, beta, true, true);
      expected =
          new List.generate(A.columns, (_) => new Float64List(2 * A.columns));
      for (int j = 0; j < A.columns; j++) {
        for (int i = 0; i < A.columns; i++) {
          var s = Complex.ZERO;
          for (int k = 0; k < A.rows; k++) {
            s += A.get(k, i).conjugate() * Bt.get(j, k).conjugate();
          }
          var tmp = new Complex(expected[i][2 * j], expected[i][2 * j + 1]);
          tmp = tmp * beta;
          tmp += s * alpha;
          expected[i][2 * j] = tmp.real;
          expected[i][2 * j + 1] = tmp.imaginary;
        }
      }
      for (int r = 0; r < A.columns; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(expected[r][2 * c], closeTo(C.get(r, c).real, TOL));
          expect(expected[r][2 * c + 1], closeTo(C.get(r, c).imaginary, TOL));
        }
      }
    });

    test('sum', () {
      var actual = A.sum();
      var expected = Complex.ZERO;
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expected += A.get(r, c);
        }
      }
      assertEquals(expected, actual, TOL);
    });
  });
}
