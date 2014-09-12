part of cern.colt.matrix.int.test;

testIntMatrix(String name, IntMatrixTest t) {
  group(name, () {
    setUp(t.setUp);
    tearDown(t.tearDown);
    test('reduce', t.testReduce);
    test('reduceWhere', t.testReduceWhere);
    test('reduceRange', t.testReduceRange);
    test('reduceWith', t.testReduceWith);
    test('fill', t.testFill);
    test('setAll2D', t.testSetAll2D);
    test('forEach', t.testForEach);
    test('copyFrom', t.testCopyFrom);
    test('forEachWith', t.testForEachWith);
    test('forEachWithRange', t.testForEachWithRange);
    test('fillWhere', t.testFillWhere);
    test('forEachWhere', t.testForEachWhere);
    test('cardinality', t.testCardinality);
    test('all', t.testAll);
    test('equals', t.testEquals);
    test('forEachNonZero', t.testForEachNonZero);
    test('max', t.testMax);
    test('min', t.testMin);
    test('negativeValues', t.testNegativeValues);
    test('nonZeros', t.testNonZeros);
    test('positiveValues', t.testPositiveValues);
    test('toList', t.testToList);
    test('vectorize', t.testVectorize);
    test('column', t.testColumn);
    test('columnFlip', t.testColumnFlip);
    test('dice', t.testDice);
    test('part', t.testPart);
    test('row', t.testRow);
    test('rowFlip', t.testRowFlip);
    test('where', t.testWhere);
    test('select', t.testSelect);
    test('strides', t.testStrides);
    test('mult', t.testMult);
    test('multiply', t.testMultiply);
    test('sum', t.testSum);
  });
}

abstract class IntMatrixTest {
  /** Matrix to test. */
  AbstractIntMatrix A;

  /** Matrix of the same size as [A]. */
  AbstractIntMatrix B;

  /** Matrix of the size `A.columns` x `A.rows`. */
  AbstractIntMatrix Bt;

  int NROWS = 13;

  int NCOLUMNS = 17;

  void setUp() {
    createMatrices();
    populateMatrices();
  }

  void createMatrices();

  void populateMatrices() {
    //ConcurrencyUtils.setThreadsBeginN_2D(1);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        A.set(r, c, math.max(1, random.nextInt(MAX_INT) % A.rows));
      }
    }

    for (int r = 0; r < B.rows; r++) {
      for (int c = 0; c < B.columns; c++) {
        B.set(r, c, math.max(1, random.nextInt(MAX_INT) % B.rows));
      }
    }

    for (int r = 0; r < Bt.rows; r++) {
      for (int c = 0; c < Bt.columns; c++) {
        Bt.set(r, c, math.max(1, random.nextInt(MAX_INT) % Bt.rows));
      }
    }
  }

  void tearDown() {
    A = B = Bt = null;
  }

  void testReduce() {
    int expected = 0;
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        int elem = A.get(r, c);
        expected += elem * elem;
      }
    }
    int result = A.reduce(ifunc.plus, ifunc.square);
    expect(expected, equals(result));
  }

  void testReduceWhere() {
    bool procedure(int element) {
      if (element.abs() > 0.2) {
        return true;
      } else {
        return false;
      }
    }
    int expected = 0;
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        int elem = A.get(r, c);
        if (elem.abs() > 0.2) {
          expected += elem * elem;
        }
      }
    }

    int result = A.reduceWhere(ifunc.plus, ifunc.square, procedure);
    expect(expected, equals(result));
  }

  void testReduceRange() {
    List<int> rowList = new List<int>();
    List<int> columnList = new List<int>();
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        rowList.add(r);
        columnList.add(c);
      }
    }
    int expected = 0;
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        int elem = A.get(r, c);
        expected += elem * elem;
      }
    }
    int result = A.reduceRange(ifunc.plus, ifunc.square,
        new Int32List.fromList(rowList),
        new Int32List.fromList(columnList));
    expect(expected, equals(result));
  }

  void testReduceWith() {
    int expected = 0;
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        int elemA = A.get(r, c);
        int elemB = B.get(r, c);
        expected += elemA * elemB;
      }
    }
    int result = A.reduceWith(B, ifunc.plus, ifunc.mult);
    expect(expected, equals(result));
  }

  void testFill() {
    int value = random.nextInt(MAX_INT);
    A.fill(value);
    List<int> rowList = new List<int>();
    List<int> columnList = new List<int>();
    List<int> valueList = new List<int>();
    A.nonZeros(rowList, columnList, valueList);
    for (int i = 0; i < valueList.length; i++) {
      expect(value, equals(valueList[i]));
    }
  }

  void testSetAll2D() {
    List<Int32List> expected = new List<Int32List>.generate(A.rows, (_) => new Int32List(A.columns));
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expected[r][c] = random.nextInt(MAX_INT);
      }
    }
    A.setAll2D(expected);
    for (int r = 0; r < A.rows; r++) {
      expect(A.columns == expected[r].length, isTrue);
      for (int c = 0; c < A.columns; c++) expect(expected[r][c], equals(A.get(r, c)));
    }
  }

  void testForEach() {
    AbstractIntMatrix Acopy = A.copy();
    A.forEach(ifunc.neg);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        int expected = -Acopy.get(r, c);
        expect(expected, equals(A.get(r, c)));
      }
    }
  }

  void testCopyFrom() {
    A.copyFrom(B);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(B.get(r, c), equals(A.get(r, c)));
      }
    }
  }

  void testForEachWith() {
    AbstractIntMatrix Acopy = A.copy();
    A.forEachWith(B, ifunc.plus);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(Acopy.get(r, c) + B.get(r, c), equals(A.get(r, c)));
      }
    }
  }

  void testForEachWithRange() {
    List<int> rowList = new List<int>();
    List<int> columnList = new List<int>();
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        rowList.add(r);
        columnList.add(c);
      }
    }
    AbstractIntMatrix Acopy = A.copy();
    A.forEachWithRange(B, ifunc.plus,
        new Int32List.fromList(rowList),
        new Int32List.fromList(columnList));
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(Acopy.get(r, c) + B.get(r, c), equals(A.get(r, c)));
      }
    }
  }

  void testFillWhere() {
    bool procedure(int element) {
      if (element.abs() > 1) {
        return true;
      } else {
        return false;
      }
    }
    AbstractIntMatrix Acopy = A.copy();
    A.fillWhere(procedure, -1);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        if (Acopy.get(r, c).abs() > 1) {
          expect(-1, equals(A.get(r, c)));
        } else {
          expect(Acopy.get(r, c), equals(A.get(r, c)));
        }
      }
    }
  }

  void testForEachWhere() {
    bool procedure(int element) {
      if (element.abs() > 1) {
        return true;
      } else {
        return false;
      }
    }
    AbstractIntMatrix Acopy = A.copy();
    A.forEachWhere(procedure, ifunc.neg);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        if (Acopy.get(r, c).abs() > 1) {
          expect(-Acopy.get(r, c), equals(A.get(r, c)));
        } else {
          expect(Acopy.get(r, c), equals(A.get(r, c)));
        }
      }
    }
  }

  void testCardinality() {
    int card = A.cardinality;
    int expected = 0;
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        if (A.get(r, c) != 0) expected++;
      }
    }
    expect(expected, equals(card));
  }

  void testAll() {
    int value = 1;
    A.fill(value);
    expect(A.all(value), isTrue);
    expect(A.all(2), isFalse);
  }

  void testEquals() {
    expect(A.equals(A), isTrue);
    expect(A.equals(B), isFalse);
  }

  void testForEachNonZero() {
    AbstractIntMatrix Acopy = A.copy();
    int function(int first, int second, int third) {
      return -third;
    }
    A.forEachNonZero(function);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(-Acopy.get(r, c), equals(A.get(r, c)));
      }
    }
  }

  void testMax() {
    A.fill(0);
    A.set(A.rows ~/ 3, A.columns ~/ 3, 7);
    A.set(A.rows ~/ 2, A.columns ~/ 2, 1);
    final maxAndLoc = A.max();
    expect(7, equals(maxAndLoc.value));
    expect(A.rows ~/ 3, equals(maxAndLoc.row));
    expect(A.columns ~/ 3, equals(maxAndLoc.column));
  }

  void testMin() {
    A.fill(0);
    A.set(A.rows ~/ 3, A.columns ~/ 3, -7);
    A.set(A.rows ~/ 2, A.columns ~/ 2, -1);
    final minAndLoc = A.min();
    expect(-7, equals(minAndLoc.value));
    expect(A.rows ~/ 3, equals(minAndLoc.row));
    expect(A.columns ~/ 3, equals(minAndLoc.column));
  }

  void testNegativeValues() {
    A.fill(0);
    A.set(A.rows ~/ 3, A.columns ~/ 3, -7);
    A.set(A.rows ~/ 2, A.columns ~/ 2, -1);
    List<int> rowList = new List<int>();
    List<int> columnList = new List<int>();
    List<int> valueList = new List<int>();
    A.negativeValues(rowList, columnList, valueList);
    expect(2, equals(rowList.length));
    expect(2, equals(columnList.length));
    expect(2, equals(valueList.length));
    expect(rowList.contains(A.rows ~/ 3), isTrue);
    expect(rowList.contains(A.rows ~/ 2), isTrue);
    expect(columnList.contains(A.columns ~/ 3), isTrue);
    expect(columnList.contains(A.columns ~/ 2), isTrue);
    expect(valueList.contains(-7), isTrue);
    expect(valueList.contains(-1), isTrue);
  }

  void testNonZeros() {
    A.fill(0);
    A.set(A.rows ~/ 3, A.columns ~/ 3, 7);
    A.set(A.rows ~/ 2, A.columns ~/ 2, 1);
    List<int> rowList = new List<int>();
    List<int> columnList = new List<int>();
    List<int> valueList = new List<int>();
    A.nonZeros(rowList, columnList, valueList);
    expect(2, equals(rowList.length));
    expect(2, equals(columnList.length));
    expect(2, equals(valueList.length));
    expect(rowList.contains(A.rows ~/ 3), isTrue);
    expect(rowList.contains(A.rows ~/ 2), isTrue);
    expect(columnList.contains(A.columns ~/ 3), isTrue);
    expect(columnList.contains(A.columns ~/ 2), isTrue);
    expect(valueList.contains(7), isTrue);
    expect(valueList.contains(1), isTrue);
  }

  void testPositiveValues() {
    A.fill(0);
    A.set(A.rows ~/ 3, A.columns ~/ 3, 7);
    A.set(A.rows ~/ 2, A.columns ~/ 2, 1);
    List<int> rowList = new List<int>();
    List<int> columnList = new List<int>();
    List<int> valueList = new List<int>();
    A.positiveValues(rowList, columnList, valueList);
    expect(2, equals(rowList.length));
    expect(2, equals(columnList.length));
    expect(2, equals(valueList.length));
    expect(rowList.contains(A.rows ~/ 3), isTrue);
    expect(rowList.contains(A.rows ~/ 2), isTrue);
    expect(columnList.contains(A.columns ~/ 3), isTrue);
    expect(columnList.contains(A.columns ~/ 2), isTrue);
    expect(valueList.contains(7), isTrue);
    expect(valueList.contains(1), isTrue);
  }

  void testToList() {
    List<Int32List> array = A.toList();
    expect(A.rows == array.length, isTrue);
    for (int r = 0; r < A.rows; r++) {
      expect(A.columns == array[r].length, isTrue);
      for (int c = 0; c < A.columns; c++) {
        expect(0, equals((array[r][c] - A.get(r, c)).abs()));
      }
    }
  }

  void testVectorize() {
    AbstractIntVector Avec = A.vectorize();
    int idx = 0;
    for (int c = 0; c < A.columns; c++) {
      for (int r = 0; r < A.rows; r++) {
        expect(A.get(r, c), equals(Avec.get(idx++)));
      }
    }
  }

  void testColumn() {
    AbstractIntVector col = A.column(A.columns ~/ 2);
    expect(A.rows, col.length);
    for (int r = 0; r < A.rows; r++) {
      expect(A.get(r, A.columns ~/ 2), equals(col.get(r)));
    }
  }

  void testColumnFlip() {
    AbstractIntMatrix B = A.columnFlip();
    expect(A.length, equals(B.length));
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(A.get(r, A.columns - 1 - c), equals(B.get(r, c)));
      }
    }
  }

  void testDice() {
    AbstractIntMatrix B = A.dice();
    expect(A.rows, equals(B.columns));
    expect(A.columns, equals(B.rows));
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(A.get(r, c), equals(B.get(c, r)));
      }
    }
  }

  void testPart() {
    AbstractIntMatrix B = A.part(A.rows ~/ 2, A.columns ~/ 2, A.rows ~/ 3, A.columns ~/ 3);
    expect(A.rows ~/ 3, equals(B.rows));
    expect(A.columns ~/ 3, equals(B.columns));
    for (int r = 0; r < A.rows / 3; r++) {
      for (int c = 0; c < A.columns / 3; c++) {
        expect(A.get(A.rows ~/ 2 + r, A.columns ~/ 2 + c), equals(B.get(r, c)));
      }
    }
  }

  void testRow() {
    AbstractIntVector B = A.row(A.rows ~/ 2);
    expect(A.columns, equals(B.length));
    for (int r = 0; r < A.columns; r++) {
      expect(A.get(A.rows ~/ 2, r), equals(B.get(r)));
    }
  }

  void testRowFlip() {
    AbstractIntMatrix B = A.rowFlip();
    expect(A.length, equals(B.length));
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(A.get(A.rows - 1 - r, c), equals(B.get(r, c)));
      }
    }
  }

  void testWhere() {
    final int value = 2;
    A.fill(0);
    A.set(A.rows ~/ 4, 0, value);
    A.set(A.rows ~/ 2, 0, value);
    AbstractIntMatrix B = A.where((IntVector element) {
      if ((element.get(0) - value).abs() == 0) {
        return true;
      } else {
        return false;
      }
    });
    expect(2, equals(B.rows));
    expect(A.columns, equals(B.columns));
    expect(A.get(A.rows ~/ 4, 0), equals(B.get(0, 0)));
    expect(A.get(A.rows ~/ 2, 0), equals(B.get(1, 0)));
  }

  void testSelect() {
    final rowIndexes = [A.rows ~/ 6, A.rows ~/ 5, A.rows ~/ 4, A.rows ~/ 3, A.rows ~/ 2];
    final colIndexes = [A.columns ~/ 6, A.columns ~/ 5, A.columns ~/ 4, A.columns ~/ 3, A.columns ~/ 2, A.columns - 1];
    AbstractIntMatrix B = A.select(new Int32List.fromList(rowIndexes),
        new Int32List.fromList(colIndexes));
    expect(rowIndexes.length, equals(B.rows));
    expect(colIndexes.length, equals(B.columns));
    for (int r = 0; r < rowIndexes.length; r++) {
      for (int c = 0; c < colIndexes.length; c++) {
        expect(A.get(rowIndexes[r], colIndexes[c]), equals(B.get(r, c)));
      }
    }
  }

  /*void testSorted() {
    IntMatrix B = A.sorted(1);
    for (int r = 0; r < A.rows - 1; r++) {
      expect(B.get(r + 1, 1) >= B.get(r, 1), isTrue);
    }
  }*/

  void testStrides() {
    int rowStride = 3;
    int colStride = 5;
    AbstractIntMatrix B = A.strides(rowStride, colStride);
    for (int r = 0; r < B.rows; r++) {
      for (int c = 0; c < B.columns; c++) {
        expect(A.get(r * rowStride, c * colStride), equals(B.get(r, c)));
      }
    }
  }

  void testMult() {
    AbstractIntVector y = new IntVector(A.columns);
    for (int i = 0; i < y.length; i++) {
      y.set(i, random.nextInt(MAX_INT) % A.rows);
    }
    int alpha = 3;
    int beta = 5;
    AbstractIntVector z = IntFactory1D.dense.random(A.rows);
    z.forEach(ifunc.modulus(A.rows));
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
    y = new IntVector(A.rows);
    for (int i = 0; i < y.length; i++) {
      y.set(i, random.nextInt(MAX_INT) % A.rows);
    }
    z = IntFactory1D.dense.random(A.columns);
    z.forEach(ifunc.modulus(A.rows));
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
  }

  void testMultiply() {
    int alpha = 3;
    int beta = 5;
    AbstractIntMatrix C = IntFactory2D.dense.random(A.rows, A.rows);
    C.forEach(ifunc.modulus(A.rows));
    List<Int32List> expected = C.toList();
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
    expected = new List<Int32List>.generate(A.rows, (_) => new Int32List(A.rows));
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
    C = IntFactory2D.dense.random(A.columns, A.columns);
    C.forEach(ifunc.modulus(A.rows));
    expected = C.toList();
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
    expected = new List<Int32List>.generate(A.columns, (_) => new Int32List(A.columns));
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
    C = IntFactory2D.dense.random(A.rows, A.rows);
    C.forEach(ifunc.modulus(A.rows));
    expected = C.toList();
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
    expected = new List<Int32List>.generate(A.rows, (_) => new Int32List(A.rows));
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
    C = IntFactory2D.dense.random(A.columns, A.columns);
    C.forEach(ifunc.modulus(A.rows));
    expected = C.toList();
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
    expected = new List<Int32List>.generate(A.columns, (_) => new Int32List(A.columns));
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

  }

  void testSum() {
    int sum = A.sum();
    int expected = 0;
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expected += A.get(r, c);
      }
    }
    expect(expected, equals(sum));
  }

}
