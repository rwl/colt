part of cern.colt.matrix.double.test;

const int NROWS = 13;
const int NCOLUMNS = 17;

void testAbstractDoubleMatrix(
    String kind, AbstractDoubleMatrix make(int rows, int columns)) {
  group("AbstractDoubleMatrix ($kind)", () {
    AbstractDoubleMatrix A, B, Bt;

    setUp(() {
      A = make(NROWS, NCOLUMNS);
      B = make(NROWS, NCOLUMNS);
      Bt = make(NCOLUMNS, NROWS);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          A.set(r, c, _r.nextDouble());
          B.set(r, c, _r.nextDouble());
          Bt.set(c, r, _r.nextDouble());
        }
      }
    });

    test('aggregate', () {
      double expected = 0.0;
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          double elem = A.get(r, c);
          expected += elem * elem;
        }
      }
      double result = A.aggregate(plus, square);
      expect(expected, closeTo(result, TOL));
    });

    test('fill', () {
      double value = _r.nextDouble();
      A.fill(value);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(value, closeTo(A.get(r, c), TOL));
        }
      }
    });

    test('apply', () {
      AbstractDoubleMatrix Acopy = A.copy();
      A.apply(acos);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          double expected = math.acos(Acopy.get(r, c));
          expect(expected, closeTo(A.get(r, c), TOL));
        }
      }
    });

    test('copyFrom', () {
      A.copyFrom(B);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(B.get(r, c), closeTo(A.get(r, c), TOL));
        }
      }
    });

    test('assign', () {
      AbstractDoubleMatrix Acopy = A.copy();
      A.assign(B, plus);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(Acopy.get(r, c) + B.get(r, c), closeTo(A.get(r, c), TOL));
        }
      }
    });

    test('cardinality', () {
      int card = A.cardinality;
      expect(A.rows * A.columns, equals(card));
    });

    test('equals', () {
      expect(A.equals(A), isTrue);
      expect(A.equals(B), isFalse);
    });

    test('forEachNonZero', () {
      AbstractDoubleMatrix Acopy = A.copy();
      double fn(int first, int second, double third) {
        return math.sqrt(third);
      }
      A.forEachNonZero(fn);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(math.sqrt(Acopy.get(r, c)), closeTo(A.get(r, c), TOL));
        }
      }
    });

    test('max', () {
      A.fill(0.0);
      A.set(A.rows ~/ 3, A.columns ~/ 3, 0.7);
      A.set(A.rows ~/ 2, A.columns ~/ 2, 0.1);
      final maxAndLoc = A.max();
      expect(0.7, closeTo(maxAndLoc.value, TOL));
      expect(A.rows ~/ 3, equals(maxAndLoc.row));
      expect(A.columns ~/ 3, equals(maxAndLoc.column));
    });

    test('min', () {
      A.fill(0.0);
      A.set(A.rows ~/ 3, A.columns ~/ 3, -0.7);
      A.set(A.rows ~/ 2, A.columns ~/ 2, -0.1);
      final minAndLoc = A.min();
      expect(-0.7, closeTo(minAndLoc.value, TOL));
      expect(A.rows ~/ 3, equals(minAndLoc.row));
      expect(A.columns ~/ 3, equals(minAndLoc.column));
    });

    test('negative', () {
      A.fill(0.0);
      A.set(A.rows ~/ 3, A.columns ~/ 3, -0.7);
      A.set(A.rows ~/ 2, A.columns ~/ 2, -0.1);
      var rowList = new List<int>();
      var columnList = new List<int>();
      var valueList = new List<double>();
      A.negative(rowList, columnList, valueList);
      expect(2, equals(rowList.length));
      expect(2, equals(columnList.length));
      expect(2, equals(valueList.length));
      expect(rowList.contains(A.rows ~/ 3), isTrue);
      expect(rowList.contains(A.rows ~/ 2), isTrue);
      expect(columnList.contains(A.columns ~/ 3), isTrue);
      expect(columnList.contains(A.columns ~/ 2), isTrue);
      expect(valueList.contains(-0.7), isTrue);
      expect(valueList.contains(-0.1), isTrue);
    });

    test('nonzero', () {
      A.fill(0.0);
      A.set(A.rows ~/ 3, A.columns ~/ 3, 0.7);
      A.set(A.rows ~/ 2, A.columns ~/ 2, 0.1);
      var rowList = new List<int>();
      var columnList = new List<int>();
      var valueList = new List<double>();
      A.nonzero(rowList, columnList, valueList);
      expect(2, equals(rowList.length));
      expect(2, equals(columnList.length));
      expect(2, equals(valueList.length));
      expect(rowList.contains(A.rows ~/ 3), isTrue);
      expect(rowList.contains(A.rows ~/ 2), isTrue);
      expect(columnList.contains(A.columns ~/ 3), isTrue);
      expect(columnList.contains(A.columns ~/ 2), isTrue);
      expect(valueList.contains(0.7), isTrue);
      expect(valueList.contains(0.1), isTrue);
    });

    test('positive', () {
      A.fill(0.0);
      A.set(A.rows ~/ 3, A.columns ~/ 3, 0.7);
      A.set(A.rows ~/ 2, A.columns ~/ 2, 0.1);
      var rowList = new List<int>();
      var columnList = new List<int>();
      var valueList = new List<double>();
      A.positive(rowList, columnList, valueList);
      expect(2, equals(rowList.length));
      expect(2, equals(columnList.length));
      expect(2, equals(valueList.length));
      expect(rowList.contains(A.rows ~/ 3), isTrue);
      expect(rowList.contains(A.rows ~/ 2), isTrue);
      expect(columnList.contains(A.columns ~/ 3), isTrue);
      expect(columnList.contains(A.columns ~/ 2), isTrue);
      expect(valueList.contains(0.7), isTrue);
      expect(valueList.contains(0.1), isTrue);
    });

    test('column', () {
      AbstractDoubleVector col = A.column(A.columns ~/ 2);
      expect(A.rows, equals(col.size));
      for (int r = 0; r < A.rows; r++) {
        expect(A.get(r, A.columns ~/ 2), closeTo(col.get(r), TOL));
      }
    });

    test('columnFlip', () {
      AbstractDoubleMatrix B = A.columnFlip();
      expect(A.size, equals(B.size));
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(A.get(r, A.columns - 1 - c), closeTo(B.get(r, c), TOL));
        }
      }
    });

    test('dice', () {
      AbstractDoubleMatrix B = A.dice();
      expect(A.rows, equals(B.columns));
      expect(A.columns, equals(B.rows));
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(A.get(r, c), closeTo(B.get(c, r), TOL));
        }
      }
    });

    test('part', () {
      AbstractDoubleMatrix B =
          A.part(A.rows ~/ 2, A.columns ~/ 2, A.rows ~/ 3, A.columns ~/ 3);
      expect(A.rows ~/ 3, equals(B.rows));
      expect(A.columns ~/ 3, equals(B.columns));
      for (int r = 0; r < A.rows / 3; r++) {
        for (int c = 0; c < A.columns / 3; c++) {
          expect(A.get(A.rows ~/ 2 + r, A.columns ~/ 2 + c),
              closeTo(B.get(r, c), TOL));
        }
      }
    });

    test('row', () {
      AbstractDoubleVector B = A.row(A.rows ~/ 2);
      expect(A.columns, equals(B.size));
      for (int r = 0; r < A.columns; r++) {
        expect(A.get(A.rows ~/ 2, r), closeTo(B.get(r), TOL));
      }
    });

    test('rowFlip', () {
      AbstractDoubleMatrix B = A.rowFlip();
      expect(A.size, equals(B.size));
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(A.get(A.rows - 1 - r, c), closeTo(B.get(r, c), TOL));
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
      AbstractDoubleMatrix B = A.select(rowIndexes, colIndexes);
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
      AbstractDoubleMatrix B = A.strides(rowStride, colStride);
      for (int r = 0; r < B.rows; r++) {
        for (int c = 0; c < B.columns; c++) {
          expect(
              A.get(r * rowStride, c * colStride), closeTo(B.get(r, c), TOL));
        }
      }
    });

    test('mult', () {
      AbstractDoubleVector y = new DoubleVector(A.columns);
      for (int i = 0; i < y.size; i++) {
        y.set(i, _r.nextDouble());
      }
      double alpha = 3.0;
      double beta = 5.0;
      AbstractDoubleVector z = new DoubleVector.random(A.rows);
      Float64List expected = z.toList();
      z = A.mult(y, z, alpha, beta, false);
      for (int r = 0; r < A.rows; r++) {
        double s = 0.0;
        for (int c = 0; c < A.columns; c++) {
          s += A.get(r, c) * y.get(c);
        }
        expected[r] = s * alpha + expected[r] * beta;
      }

      for (int r = 0; r < A.rows; r++) {
        expect(expected[r], closeTo(z.get(r), TOL));
      }
      //---
      z = null;
      z = A.mult(y, z, alpha, beta, false);
      expected = new Float64List(A.rows);
      for (int r = 0; r < A.rows; r++) {
        double s = 0.0;
        for (int c = 0; c < A.columns; c++) {
          s += A.get(r, c) * y.get(c);
        }
        expected[r] = s * alpha;
      }
      for (int r = 0; r < A.rows; r++) {
        expect(expected[r], closeTo(z.get(r), TOL));
      }

      //transpose
      y = new DoubleVector(A.rows);
      for (int i = 0; i < y.size; i++) {
        y.set(i, _r.nextDouble());
      }
      z = new DoubleVector.random(A.columns);
      expected = z.toList();
      z = A.mult(y, z, alpha, beta, true);
      for (int r = 0; r < A.columns; r++) {
        double s = 0.0;
        for (int c = 0; c < A.rows; c++) {
          s += A.get(c, r) * y.get(c);
        }
        expected[r] = s * alpha + expected[r] * beta;
      }
      for (int r = 0; r < A.columns; r++) {
        expect(expected[r], closeTo(z.get(r), TOL));
      }
      //---
      z = null;
      z = A.mult(y, z, alpha, beta, true);
      expected = new Float64List(A.columns);
      for (int r = 0; r < A.columns; r++) {
        double s = 0.0;
        for (int c = 0; c < A.rows; c++) {
          s += A.get(c, r) * y.get(c);
        }
        expected[r] = s * alpha;
      }
      for (int r = 0; r < A.columns; r++) {
        expect(expected[r], closeTo(z.get(r), TOL));
      }
    });

    test('multiply', () {
      double alpha = 3.0;
      double beta = 5.0;
      AbstractDoubleMatrix C = new DoubleMatrix.random(A.rows, A.rows);
      List<Float64List> expected = toList(C);
      C = A.multiply(Bt, C, alpha, beta, false, false);
      for (int j = 0; j < A.rows; j++) {
        for (int i = 0; i < A.rows; i++) {
          double s = 0.0;
          for (int k = 0; k < A.columns; k++) {
            s += A.get(i, k) * Bt.get(k, j);
          }
          expected[i][j] = s * alpha + expected[i][j] * beta;
        }
      }
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.rows; c++) {
          expect(expected[r][c], closeTo(C.get(r, c), TOL));
        }
      }

      //---
      C = null;
      C = A.multiply(Bt, C, alpha, beta, false, false);
      expected = new List<Float64List>.generate(
          A.rows, (_) => new Float64List(A.rows));
      for (int j = 0; j < A.rows; j++) {
        for (int i = 0; i < A.rows; i++) {
          double s = 0.0;
          for (int k = 0; k < A.columns; k++) {
            s += A.get(i, k) * Bt.get(k, j);
          }
          expected[i][j] = s * alpha;
        }
      }
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.rows; c++) {
          expect(expected[r][c], closeTo(C.get(r, c), TOL));
        }
      }

      //transposeA
      C = new DoubleMatrix.random(A.columns, A.columns);
      expected = toList(C);
      C = A.multiply(B, C, alpha, beta, true, false);
      for (int j = 0; j < A.columns; j++) {
        for (int i = 0; i < A.columns; i++) {
          double s = 0.0;
          for (int k = 0; k < A.rows; k++) {
            s += A.get(k, i) * B.get(k, j);
          }
          expected[i][j] = s * alpha + expected[i][j] * beta;
        }
      }
      for (int r = 0; r < A.columns; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(expected[r][c], closeTo(C.get(r, c), TOL));
        }
      }
      //---
      C = null;
      C = A.multiply(B, C, alpha, beta, true, false);
      expected = new List<Float64List>.generate(
          A.columns, (_) => new Float64List(A.columns));
      for (int j = 0; j < A.columns; j++) {
        for (int i = 0; i < A.columns; i++) {
          double s = 0.0;
          for (int k = 0; k < A.rows; k++) {
            s += A.get(k, i) * B.get(k, j);
          }
          expected[i][j] = s * alpha;
        }
      }
      for (int r = 0; r < A.columns; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(expected[r][c], closeTo(C.get(r, c), TOL));
        }
      }

      //transposeB
      C = new DoubleMatrix.random(A.rows, A.rows);
      expected = toList(C);
      C = A.multiply(B, C, alpha, beta, false, true);
      for (int j = 0; j < A.rows; j++) {
        for (int i = 0; i < A.rows; i++) {
          double s = 0.0;
          for (int k = 0; k < A.columns; k++) {
            s += A.get(i, k) * B.get(j, k);
          }
          expected[i][j] = s * alpha + expected[i][j] * beta;
        }
      }
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.rows; c++) {
          expect(expected[r][c], closeTo(C.get(r, c), TOL));
        }
      }
      //---
      C = null;
      C = A.multiply(B, C, alpha, beta, false, true);
      expected = new List<Float64List>.generate(
          A.rows, (_) => new Float64List(A.rows));
      for (int j = 0; j < A.rows; j++) {
        for (int i = 0; i < A.rows; i++) {
          double s = 0.0;
          for (int k = 0; k < A.columns; k++) {
            s += A.get(i, k) * B.get(j, k);
          }
          expected[i][j] = s * alpha;
        }
      }
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.rows; c++) {
          expect(expected[r][c], closeTo(C.get(r, c), TOL));
        }
      }
      //transposeA and transposeB
      C = new DoubleMatrix.random(A.columns, A.columns);
      expected = toList(C);
      C = A.multiply(Bt, C, alpha, beta, true, true);
      for (int j = 0; j < A.columns; j++) {
        for (int i = 0; i < A.columns; i++) {
          double s = 0.0;
          for (int k = 0; k < A.rows; k++) {
            s += A.get(k, i) * Bt.get(j, k);
          }
          expected[i][j] = s * alpha + expected[i][j] * beta;
        }
      }
      for (int r = 0; r < A.columns; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(expected[r][c], closeTo(C.get(r, c), TOL));
        }
      }
      //---
      C = null;
      C = A.multiply(Bt, C, alpha, beta, true, true);
      expected = new List<Float64List>.generate(
          A.columns, (_) => new Float64List(A.columns));
      for (int j = 0; j < A.columns; j++) {
        for (int i = 0; i < A.columns; i++) {
          double s = 0.0;
          for (int k = 0; k < A.rows; k++) {
            s += A.get(k, i) * Bt.get(j, k);
          }
          expected[i][j] = s * alpha;
        }
      }
      for (int r = 0; r < A.columns; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(expected[r][c], closeTo(C.get(r, c), TOL));
        }
      }
    });

    test('sum', () {
      double sum = A.sum();
      double expected = 0.0;
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expected += A.get(r, c);
        }
      }
      expect(expected, closeTo(sum, TOL));
    });
  });
}
