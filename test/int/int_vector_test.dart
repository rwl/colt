part of cern.colt.matrix.int.test;

testIntVector(String name, IntVectorTest t) {
  group(name, () {
    setUp(t.setUp);
    tearDown(t.tearDown);
    test('reduce', t.testReduce);
    test('reduceRange', t.testReduceRange);
    test('reduceWith', t.testReduceWith);
    test('fill', t.testFill);
    test('setAll', t.testSetAll);
    test('forEach', t.testForEach);
    test('copyFrom', t.testCopyFrom);
    test('forEachWith', t.testForEachWith);
    test('fillWhere', t.testFillWhere);
    test('forEachWhere', t.testForEachWhere);
    test('cardinality', t.testCardinality);
    test('all', t.testAll);
    test('equals', t.testEquals);
    test('max', t.testMax);
    test('min', t.testMin);
    test('negativeValues', t.testNegativeValues);
    test('nonZeros', t.testNonZeros);
    test('positiveValues', t.testPositiveValues);
    test('toList', t.testToList);
    test('fillList', t.testFillList);
    test('reshape', t.testReshape);
    test('swap', t.testSwap);
    test('flip', t.testFlip);
    test('part', t.testPart);
    test('where', t.testWhere);
    test('select', t.testSelect);
    test('strides', t.testStrides);
    test('dot', t.testDot);
    test('dotRange', t.testDotRange);
    test('dotNonZero', t.testDotNonZero);
    test('sum', t.testSum);
  });
}

abstract class IntVectorTest {
  /** Matrix to test. */
  AbstractIntVector A;

  /** Matrix of the same size as [A]. */
  AbstractIntVector B;

  int SIZE = 2 * 17 * 5;

  void setUp() {
    createMatrices();
    populateMatrices();
  }

  void createMatrices();

  void populateMatrices() {
    //ConcurrencyUtils.setThreadsBeginN_1D(1);

    for (int i = 0; i < A.length; i++) {
      A.set(i, math.max(1, random.nextInt(MAX_INT) % A.length));
    }

    for (int i = 0; i < B.length; i++) {
      B.set(i, math.max(1, random.nextInt(MAX_INT) % A.length));
    }
  }

  void tearDown() {
    A = B = null;
  }

  void testReduce() {
    int expected = 0;
    for (int i = 0; i < A.length; i++) {
      int elem = A.get(i);
      expected += elem * elem;
    }
    int result = A.aggregate(ifunc.plus, ifunc.square);
    expect(expected, equals(result));
  }

  void testReduceRange() {
    List<int> indexList = new List<int>();
    for (int i = 0; i < A.length; i++) {
      indexList.add(i);
    }
    int expected = 0;
    for (int i = 0; i < A.length; i++) {
      int elem = A.get(i);
      expected += elem * elem;
    }
    int result = A.reduceRange(ifunc.plus, ifunc.square,
        new Int32List.fromList(indexList));
    expect(expected, equals(result));
  }

  void testReduceWith() {
    int expected = 0;
    for (int i = 0; i < A.length; i++) {
      int elemA = A.get(i);
      int elemB = B.get(i);
      expected += elemA * elemB;
    }
    int result = A.reduceWith(B, ifunc.plus, ifunc.mult);
    expect(expected, equals(result));
  }

  void testFill() {
    int value = random.nextInt(MAX_INT);
    A.fill(value);
    for (int i = 0; i < A.length; i++) {
      expect(value, equals(A.get(i)));
    }
  }

  void testSetAll() {
    Int32List expected = new Int32List(A.length);
    for (int i = 0; i < A.length; i++) {
      expected[i] = random.nextInt(MAX_INT);
    }
    A.setAll(0, expected);
    for (int i = 0; i < A.length; i++) {
      expect(expected[i], A.get(i));
    }
  }

  void testForEach() {
    AbstractIntVector Acopy = A.copy();
    A.forEach(ifunc.neg);
    for (int i = 0; i < A.length; i++) {
      int expected = -Acopy.get(i);
      expect(expected, equals(A.get(i)));
    }
  }

  void testCopyFrom() {
    A.copyFrom(B);
    expect(A.length == B.length, isTrue);
    for (int i = 0; i < A.length; i++) {
      expect(B.get(i), equals(A.get(i)));
    }
  }

  void testForEachWith() {
    AbstractIntVector Acopy = A.copy();
    A.forEachWith(B, ifunc.plus);
    for (int i = 0; i < A.length; i++) {
      expect(Acopy.get(i) + B.get(i), equals(A.get(i)));
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
    AbstractIntVector Acopy = A.copy();
    A.fillWhere(procedure, -1);
    for (int i = 0; i < A.length; i++) {
      if (Acopy.get(i).abs() > 1) {
        expect(-1, equals(A.get(i)));
      } else {
        expect(Acopy.get(i), equals(A.get(i)));
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
    AbstractIntVector Acopy = A.copy();
    A.forEachWhere(procedure, ifunc.neg);
    for (int i = 0; i < A.length; i++) {
      if (Acopy.get(i).abs() > 1) {
        expect(-Acopy.get(i), equals(A.get(i)));
      } else {
        expect(Acopy.get(i), equals(A.get(i)));
      }
    }
  }

  void testCardinality() {
    int card = A.cardinality;
    int expected = 0;
    for (int i = 0; i < A.length; i++) {
      if (A.get(i) != 0) expected++;
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

  void testMax() {
    A.fill(0);
    A.set(A.length ~/ 3, 7);
    A.set(A.length ~/ 2, 1);
    final maxAndLoc = A.max();
    expect(7, equals(maxAndLoc.value));
    expect(A.length ~/ 3, equals(maxAndLoc.location));
  }

  void testMin() {
    A.fill(0);
    A.set(A.length ~/ 3, -7);
    A.set(A.length ~/ 2, -1);
    final minAndLoc = A.min();
    expect(-7, equals(minAndLoc.value));
    expect(A.length ~/ 3, equals(minAndLoc.location));
  }

  void testNegativeValues() {
    A.fill(0);
    A.set(A.length ~/ 3, -7);
    A.set(A.length ~/ 2, -1);
    List<int> indexList = new List<int>();
    List<int> valueList = new List<int>();
    A.negativeValues(indexList: indexList, valueList: valueList);
    expect(2, equals(indexList.length));
    expect(2, equals(valueList.length));
    expect(indexList.contains(A.length ~/ 3), isTrue);
    expect(indexList.contains(A.length ~/ 2), isTrue);
    expect(valueList.contains(-7), isTrue);
    expect(valueList.contains(-1), isTrue);
  }

  void testNonZeros() {
    A.fill(0);
    A.set(A.length ~/ 3, 7);
    A.set(A.length ~/ 2, 1);
    List<int> indexList = new List<int>();
    List<int> valueList = new List<int>();
    A.nonZeros(indexList: indexList, valueList: valueList);
    expect(2, equals(indexList.length));
    expect(2, equals(valueList.length));
    expect(indexList.contains(A.length ~/ 3), isTrue);
    expect(indexList.contains(A.length ~/ 2), isTrue);
    expect(valueList.contains(7), isTrue);
    expect(valueList.contains(1), isTrue);
  }

  void testPositiveValues() {
    A.fill(0);
    A.set(A.length ~/ 3, 7);
    A.set(A.length ~/ 2, 1);
    List<int> indexList = new List<int>();
    List<int> valueList = new List<int>();
    A.positiveValues(indexList: indexList, valueList: valueList);
    expect(2, equals(indexList.length));
    expect(2, equals(valueList.length));
    expect(indexList.contains(A.length ~/ 3), isTrue);
    expect(indexList.contains(A.length ~/ 2), isTrue);
    expect(valueList.contains(7), isTrue);
    expect(valueList.contains(1), isTrue);
  }

  void testToList() {
    final array = A.toList();
    expect(A.length == array.length, isTrue);
    for (int i = 0; i < A.length; i++) {
      expect(array[i], equals(A.get(i)));
    }
  }

  void testFillList() {
    final array = new Int32List(A.length);
    A.fillList(array);
    for (int i = 0; i < A.length; i++) {
      expect(A.get(i), equals(array[i]));
    }
  }

  void testReshape() {
    int rows = 10;
    int columns = 17;
    AbstractIntMatrix B = A.reshape(rows, columns);
    int idx = 0;
    for (int c = 0; c < columns; c++) {
      for (int r = 0; r < rows; r++) {
        expect(A.get(idx++), equals(B.get(r, c)));
      }
    }
  }

  /*void testReshape3D() {
    int slices = 2;
    int rows = 5;
    int columns = 17;
    final B = A.reshape3D(slices, rows, columns);
    int idx = 0;
    for (int s = 0; s < slices; s++) {
      for (int c = 0; c < columns; c++) {
        for (int r = 0; r < rows; r++) {
          expect(A.get(idx++), equals(B.getQuick(s, r, c)));
        }
      }
    }
  }*/

  void testSwap() {
    AbstractIntVector Acopy = A.copy();
    AbstractIntVector Bcopy = B.copy();
    A.swap(B);
    for (int i = 0; i < A.length; i++) {
      expect(Bcopy.get(i), equals(A.get(i)));
      expect(Acopy.get(i), equals(B.get(i)));
    }
  }

  void testFlip() {
    AbstractIntVector b = A.flip();
    expect(A.length, b.length);
    for (int i = 0; i < A.length; i++) {
      expect(A.get(i), equals(b.get(A.length - 1 - i)));
    }
  }

  void testPart() {
    AbstractIntVector b = A.part(15, 11);
    for (int i = 0; i < 11; i++) {
      expect(A.get(15 + i), equals(b.get(i)));
    }
  }

  void testWhere() {
    AbstractIntVector b = A.where((int element) {
      return element % 2 == 0;
    });
    for (int i = 0; i < b.length; i++) {
      int el = b.get(i);
      if (el % 2 != 0) {
        fail('');
      }
    }
  }

  void testSelect() {
    final indexes = new Int32List.fromList([5, 11, 22, 37, 101]);
    AbstractIntVector b = A.select(indexes);
    for (int i = 0; i < indexes.length; i++) {
      expect(A.get(indexes[i]), equals(b.get(i)));
    }
  }

  /*void testSorted() {
    IntVector b = A.sorted();
    for (int i = 0; i < A.length - 1; i++) {
      expect(b.get(i + 1) >= b.get(i), isTrue);
    }
  }*/

  void testStrides() {
    int stride = 3;
    AbstractIntVector b = A.strides(stride);
    for (int i = 0; i < b.length; i++) {
      expect(A.get(i * stride), equals(b.get(i)));
    }
  }

  void testDot() {
    int product = A.dot(B);
    int expected = 0;
    for (int i = 0; i < A.length; i++) {
      expected += A.get(i) * B.get(i);
    }
    expect(expected, equals(product));
  }

  void testDotRange() {
    int product = A.dot(B, 5, B.length - 10);
    int expected = 0;
    for (int i = 5; i < A.length - 5; i++) {
      expected += A.get(i) * B.get(i);
    }
    expect(expected, equals(product));
  }

  void testDotNonZero() {
    List<int> indexList = new List<int>();
    List<int> valueList = new List<int>();
    B.nonZeros(indexList: indexList, valueList: valueList);
    int product = A.dotNonZero(B, indexList, 5, B.length - 10);
    int expected = 0;
    for (int i = 5; i < A.length - 5; i++) {
      expected += A.get(i) * B.get(i);
    }
    expect(expected, equals(product));
  }

  void testSum() {
    int sum = A.sum();
    int expected = 0;
    for (int i = 0; i < A.length; i++) {
      expected += A.get(i);
    }
    expect(expected, equals(sum));
  }

}
