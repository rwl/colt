part of cern.colt.matrix.int.test;

const int NROWS = 13;

const int NCOLUMNS = 17;

testIntMatrix(String kind, IntMatrix make(int rows, int columns)) {
  group("IntMatrix ($kind)", () {
    IntMatrix A, B, Bt;
    setUp(() {
      A = make(NROWS, NCOLUMNS);
      B = make(NROWS, NCOLUMNS);
      Bt = make(NCOLUMNS, NROWS);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          A.set(r, c, math.max(1, random.nextInt(MAX_INT) % A.rows));
          B.set(r, c, math.max(1, random.nextInt(MAX_INT) % A.rows));
          Bt.set(c, r, math.max(1, random.nextInt(MAX_INT) % A.rows));
        }
      }
    });

    test('aggregate', () {
      int expected = 0;
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          int elem = A.get(r, c);
          expected += elem * elem;
        }
      }
      int result = A.aggregate(ifunc.plus, ifunc.square);
      expect(expected, equals(result));
    });

    test('fill', () {
      int value = random.nextInt(MAX_INT);
      A.fill(value);
      List<int> rowList = new List<int>();
      List<int> columnList = new List<int>();
      List<int> valueList = new List<int>();
      A.nonzero(rowList, columnList, valueList);
      for (int i = 0; i < valueList.length; i++) {
        expect(value, equals(valueList[i]));
      }
    });

    test('apply', () {
      IntMatrix Acopy = A.copy();
      A.apply(ifunc.neg);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          int expected = -Acopy.get(r, c);
          expect(expected, equals(A.get(r, c)));
        }
      }
    });

    test('copyFrom', () {
      A.copyFrom(B);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(B.get(r, c), equals(A.get(r, c)));
        }
      }
    });

    test('assign', () {
      IntMatrix Acopy = A.copy();
      A.assign(B, ifunc.plus);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(Acopy.get(r, c) + B.get(r, c), equals(A.get(r, c)));
        }
      }
    });

    test('cardinality', () {
      int card = A.cardinality;
      int expected = 0;
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          if (A.get(r, c) != 0) expected++;
        }
      }
      expect(expected, equals(card));
    });

    test('equals', () {
      expect(A.equals(A), isTrue);
      expect(A.equals(B), isFalse);
    });

    test('forEachNonZero', () {
      IntMatrix Acopy = A.copy();
      int fn(int first, int second, int third) {
        return -third;
      }
      A.forEachNonZero(fn);
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(-Acopy.get(r, c), equals(A.get(r, c)));
        }
      }
    });

    test('max', () {
      A.fill(0);
      A.set(A.rows ~/ 3, A.columns ~/ 3, 7);
      A.set(A.rows ~/ 2, A.columns ~/ 2, 1);
      var maxAndLoc = A.max();
      expect(7, equals(maxAndLoc.value));
      expect(A.rows ~/ 3, equals(maxAndLoc.row));
      expect(A.columns ~/ 3, equals(maxAndLoc.column));
    });

    test('min', () {
      A.fill(0);
      A.set(A.rows ~/ 3, A.columns ~/ 3, -7);
      A.set(A.rows ~/ 2, A.columns ~/ 2, -1);
      var minAndLoc = A.min();
      expect(-7, equals(minAndLoc.value));
      expect(A.rows ~/ 3, equals(minAndLoc.row));
      expect(A.columns ~/ 3, equals(minAndLoc.column));
    });

    test('negative', () {
      A.fill(0);
      A.set(A.rows ~/ 3, A.columns ~/ 3, -7);
      A.set(A.rows ~/ 2, A.columns ~/ 2, -1);
      var rowList = <int>[];
      var columnList = <int>[];
      var valueList = <int>[];
      A.negative(rowList, columnList, valueList);
      expect(rowList.length, equals(2));
      expect(columnList.length, equals(2));
      expect(valueList.length, equals(2));
      expect(rowList.contains(A.rows ~/ 3), isTrue);
      expect(rowList.contains(A.rows ~/ 2), isTrue);
      expect(columnList.contains(A.columns ~/ 3), isTrue);
      expect(columnList.contains(A.columns ~/ 2), isTrue);
      expect(valueList.contains(-7), isTrue);
      expect(valueList.contains(-1), isTrue);
    });

    test('nonzero', () {
      A.fill(0);
      A.set(A.rows ~/ 3, A.columns ~/ 3, 7);
      A.set(A.rows ~/ 2, A.columns ~/ 2, 1);
      var rowList = <int>[];
      var columnList = <int>[];
      var valueList = <int>[];
      A.nonzero(rowList, columnList, valueList);
      expect(rowList.length, equals(2));
      expect(columnList.length, equals(2));
      expect(valueList.length, equals(2));
      expect(rowList.contains(A.rows ~/ 3), isTrue);
      expect(rowList.contains(A.rows ~/ 2), isTrue);
      expect(columnList.contains(A.columns ~/ 3), isTrue);
      expect(columnList.contains(A.columns ~/ 2), isTrue);
      expect(valueList.contains(7), isTrue);
      expect(valueList.contains(1), isTrue);
    });

    test('positive', () {
      A.fill(0);
      A.set(A.rows ~/ 3, A.columns ~/ 3, 7);
      A.set(A.rows ~/ 2, A.columns ~/ 2, 1);
      var rowList = <int>[];
      var columnList = <int>[];
      var valueList = <int>[];
      A.positive(rowList, columnList, valueList);
      expect(rowList.length, equals(2));
      expect(columnList.length, equals(2));
      expect(valueList.length, equals(2));
      expect(rowList.contains(A.rows ~/ 3), isTrue);
      expect(rowList.contains(A.rows ~/ 2), isTrue);
      expect(columnList.contains(A.columns ~/ 3), isTrue);
      expect(columnList.contains(A.columns ~/ 2), isTrue);
      expect(valueList.contains(7), isTrue);
      expect(valueList.contains(1), isTrue);
    });

    test('column', () {
      IntVector col = A.column(A.columns ~/ 2);
      expect(A.rows, col.size);
      for (int r = 0; r < A.rows; r++) {
        expect(A.get(r, A.columns ~/ 2), equals(col.get(r)));
      }
    });

    test('columnFlip', () {
      IntMatrix B = A.columnFlip();
      expect(A.size, equals(B.size));
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(A.get(r, A.columns - 1 - c), equals(B.get(r, c)));
        }
      }
    });

    test('dice', () {
      IntMatrix B = A.dice();
      expect(A.rows, equals(B.columns));
      expect(A.columns, equals(B.rows));
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(A.get(r, c), equals(B.get(c, r)));
        }
      }
    });

    test('part', () {
      IntMatrix B =
      A.part(A.rows ~/ 2, A.columns ~/ 2, A.rows ~/ 3, A.columns ~/ 3);
      expect(A.rows ~/ 3, equals(B.rows));
      expect(A.columns ~/ 3, equals(B.columns));
      for (int r = 0; r < A.rows / 3; r++) {
        for (int c = 0; c < A.columns / 3; c++) {
          expect(A.get(A.rows ~/ 2 + r, A.columns ~/ 2 + c), equals(B.get(r, c)));
        }
      }
    });

    test('row', () {
      IntVector B = A.row(A.rows ~/ 2);
      expect(A.columns, equals(B.size));
      for (int r = 0; r < A.columns; r++) {
        expect(A.get(A.rows ~/ 2, r), equals(B.get(r)));
      }
    });

    test('rowFlip', () {
      IntMatrix B = A.rowFlip();
      expect(A.size, equals(B.size));
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(A.get(A.rows - 1 - r, c), equals(B.get(r, c)));
        }
      }
    });

    test('select', () {
      var rowIndexes = [
        A.rows ~/ 6,
        A.rows ~/ 5,
        A.rows ~/ 4,
        A.rows ~/ 3,
        A.rows ~/ 2
      ];
      var colIndexes = [
        A.columns ~/ 6,
        A.columns ~/ 5,
        A.columns ~/ 4,
        A.columns ~/ 3,
        A.columns ~/ 2,
        A.columns - 1
      ];
      IntMatrix B = A.select(
          new Int32List.fromList(rowIndexes), new Int32List.fromList(colIndexes));
      expect(rowIndexes.length, equals(B.rows));
      expect(colIndexes.length, equals(B.columns));
      for (int r = 0; r < rowIndexes.length; r++) {
        for (int c = 0; c < colIndexes.length; c++) {
          expect(A.get(rowIndexes[r], colIndexes[c]), equals(B.get(r, c)));
        }
      }
    });

    test('strides', () {
      int rowStride = 3;
      int colStride = 5;
      IntMatrix B = A.strides(rowStride, colStride);
      for (int r = 0; r < B.rows; r++) {
        for (int c = 0; c < B.columns; c++) {
          expect(A.get(r * rowStride, c * colStride), equals(B.get(r, c)));
        }
      }
    });

    test('mult', () {
      IntVector y = new DenseIntVector(A.columns);
      for (int i = 0; i < y.size; i++) {
        y.set(i, random.nextInt(MAX_INT) % A.rows);
      }
      int alpha = 3;
      int beta = 5;
      IntVector z = new DenseIntVector.random(A.rows);
      z.apply(ifunc.modulus(A.rows));
      Int32List expected = z.toList();
      z = A.mult(y, z, alpha, beta, false);
      for (int r = 0; r < A.rows; r++) {
        int s = 0;
        for (int c = 0; c < A.columns; c++) {
          s += A.get(r, c) * y.get(c);
        }
        expected[r] = s * alpha + expected[r] * beta;
      }

      for (int r = 0; r < A.rows; r++) {
        expect(expected[r], equals(z.get(r)));
      }
      //---
      z = null;
      z = A.mult(y, z, alpha, beta, false);
      expected = new Int32List(A.rows);
      for (int r = 0; r < A.rows; r++) {
        int s = 0;
        for (int c = 0; c < A.columns; c++) {
          s += A.get(r, c) * y.get(c);
        }
        expected[r] = s * alpha;
      }
      for (int r = 0; r < A.rows; r++) {
        expect(expected[r], equals(z.get(r)));
      }

      //transpose
      y = new DenseIntVector(A.rows);
      for (int i = 0; i < y.size; i++) {
        y.set(i, random.nextInt(MAX_INT) % A.rows);
      }
      z = new DenseIntVector.random(A.columns);
      z.apply(ifunc.modulus(A.rows));
      expected = z.toList();
      z = A.mult(y, z, alpha, beta, true);
      for (int r = 0; r < A.columns; r++) {
        int s = 0;
        for (int c = 0; c < A.rows; c++) {
          s += A.get(c, r) * y.get(c);
        }
        expected[r] = s * alpha + expected[r] * beta;
      }
      for (int r = 0; r < A.columns; r++) {
        expect(expected[r], equals(z.get(r)));
      }
      //---
      z = null;
      z = A.mult(y, z, alpha, beta, true);
      expected = new Int32List(A.columns);
      for (int r = 0; r < A.columns; r++) {
        int s = 0;
        for (int c = 0; c < A.rows; c++) {
          s += A.get(c, r) * y.get(c);
        }
        expected[r] = s * alpha;
      }
      for (int r = 0; r < A.columns; r++) {
        expect(expected[r], equals(z.get(r)));
      }
    });

    test('multiply', () {
      int alpha = 3;
      int beta = 5;
      IntMatrix C = new DenseIntMatrix.random(A.rows, A.rows);
      C.apply(ifunc.modulus(A.rows));
      List<Int32List> expected = toList(C);
      C = A.multiply(Bt, C, alpha, beta, false, false);
      for (int j = 0; j < A.rows; j++) {
        for (int i = 0; i < A.rows; i++) {
          int s = 0;
          for (int k = 0; k < A.columns; k++) {
            s += A.get(i, k) * Bt.get(k, j);
          }
          expected[i][j] = s * alpha + expected[i][j] * beta;
        }
      }
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.rows; c++) {
          expect(expected[r][c], equals(C.get(r, c)));
        }
      }

      //---
      C = null;
      C = A.multiply(Bt, C, alpha, beta, false, false);
      expected =
      new List<Int32List>.generate(A.rows, (_) => new Int32List(A.rows));
      for (int j = 0; j < A.rows; j++) {
        for (int i = 0; i < A.rows; i++) {
          int s = 0;
          for (int k = 0; k < A.columns; k++) {
            s += A.get(i, k) * Bt.get(k, j);
          }
          expected[i][j] = s * alpha;
        }
      }
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.rows; c++) {
          expect(expected[r][c], equals(C.get(r, c)));
        }
      }

      //transposeA
      C = new DenseIntMatrix.random(A.columns, A.columns);
      C.apply(ifunc.modulus(A.rows));
      expected = toList(C);
      C = A.multiply(B, C, alpha, beta, true, false);
      for (int j = 0; j < A.columns; j++) {
        for (int i = 0; i < A.columns; i++) {
          int s = 0;
          for (int k = 0; k < A.rows; k++) {
            s += A.get(k, i) * B.get(k, j);
          }
          expected[i][j] = s * alpha + expected[i][j] * beta;
        }
      }
      for (int r = 0; r < A.columns; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(expected[r][c], equals(C.get(r, c)));
        }
      }
      //---
      C = null;
      C = A.multiply(B, C, alpha, beta, true, false);
      expected = new List<Int32List>.generate(
          A.columns, (_) => new Int32List(A.columns));
      for (int j = 0; j < A.columns; j++) {
        for (int i = 0; i < A.columns; i++) {
          int s = 0;
          for (int k = 0; k < A.rows; k++) {
            s += A.get(k, i) * B.get(k, j);
          }
          expected[i][j] = s * alpha;
        }
      }
      for (int r = 0; r < A.columns; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(expected[r][c], equals(C.get(r, c)));
        }
      }

      //transposeB
      C = new DenseIntMatrix.random(A.rows, A.rows);
      C.apply(ifunc.modulus(A.rows));
      expected = toList(C);
      C = A.multiply(B, C, alpha, beta, false, true);
      for (int j = 0; j < A.rows; j++) {
        for (int i = 0; i < A.rows; i++) {
          int s = 0;
          for (int k = 0; k < A.columns; k++) {
            s += A.get(i, k) * B.get(j, k);
          }
          expected[i][j] = s * alpha + expected[i][j] * beta;
        }
      }
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.rows; c++) {
          expect(expected[r][c], equals(C.get(r, c)));
        }
      }
      //---
      C = null;
      C = A.multiply(B, C, alpha, beta, false, true);
      expected =
      new List<Int32List>.generate(A.rows, (_) => new Int32List(A.rows));
      for (int j = 0; j < A.rows; j++) {
        for (int i = 0; i < A.rows; i++) {
          int s = 0;
          for (int k = 0; k < A.columns; k++) {
            s += A.get(i, k) * B.get(j, k);
          }
          expected[i][j] = s * alpha;
        }
      }
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.rows; c++) {
          expect(expected[r][c], equals(C.get(r, c)));
        }
      }
      //transposeA and transposeB
      C = new DenseIntMatrix.random(A.columns, A.columns);
      C.apply(ifunc.modulus(A.rows));
      expected = toList(C);
      C = A.multiply(Bt, C, alpha, beta, true, true);
      for (int j = 0; j < A.columns; j++) {
        for (int i = 0; i < A.columns; i++) {
          int s = 0;
          for (int k = 0; k < A.rows; k++) {
            s += A.get(k, i) * Bt.get(j, k);
          }
          expected[i][j] = s * alpha + expected[i][j] * beta;
        }
      }
      for (int r = 0; r < A.columns; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(expected[r][c], equals(C.get(r, c)));
        }
      }
      //---
      C = null;
      C = A.multiply(Bt, C, alpha, beta, true, true);
      expected = new List<Int32List>.generate(
          A.columns, (_) => new Int32List(A.columns));
      for (int j = 0; j < A.columns; j++) {
        for (int i = 0; i < A.columns; i++) {
          int s = 0;
          for (int k = 0; k < A.rows; k++) {
            s += A.get(k, i) * Bt.get(j, k);
          }
          expected[i][j] = s * alpha;
        }
      }
      for (int r = 0; r < A.columns; r++) {
        for (int c = 0; c < A.columns; c++) {
          expect(expected[r][c], equals(C.get(r, c)));
        }
      }
    });

    test('sum', () {
      int sum = A.sum();
      int expected = 0;
      for (int r = 0; r < A.rows; r++) {
        for (int c = 0; c < A.columns; c++) {
          expected += A.get(r, c);
        }
      }
      expect(expected, equals(sum));
    });
  });
}
