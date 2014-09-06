part of cern.colt.matrix.int.test;

abstract class IntMatrix1DTest {
  /** Matrix to test. */
  IntMatrix1D A;

  /** Matrix of the same size as [A]. */
  IntMatrix1D B;

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
    int result = A.reduce(ifunc.plus, ifunc.square);
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
    int result = A.reduceRange(ifunc.plus, ifunc.square, indexList);
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
    A.setAll(expected);
    for (int i = 0; i < A.length; i++) {
      expect(expected[i], A.get(i));
    }
  }

  void testForEach() {
    IntMatrix1D Acopy = A.copy();
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
    IntMatrix1D Acopy = A.copy();
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
    IntMatrix1D Acopy = A.copy();
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
    IntMatrix1D Acopy = A.copy();
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
    int card = A.cardinality();
    int expected = 0;
    for (int i = 0; i < A.length; i++) {
      if (A.get(i) != 0) expected++;
    }
    expect(expected, equals(card));
  }

  void testEqualsInt() {
    int value = 1;
    A.fill(value);
    expect(A == value, isTrue);
    expect(A == 2, isFalse);
  }

  void testEqualsObject() {
    expect(A == A, isTrue);
    expect(A == B, isFalse);
  }

  void testMax() {
    A.fill(0);
    A.set(A.length ~/ 3, 7);
    A.set(A.length ~/ 2, 1);
    final maxAndLoc = A.max();
    expect(7, equals(maxAndLoc[0]));
    expect(A.length / 3, equals(maxAndLoc[1]));
  }

  void testMin() {
    A.fill(0);
    A.set(A.length ~/ 3, -7);
    A.set(A.length ~/ 2, -1);
    final minAndLoc = A.min();
    expect(-7, equals(minAndLoc[0]));
    expect(A.length / 3, equals(minAndLoc[1]));
  }

  void testNegativeValues() {
    A.fill(0);
    A.set(A.length ~/ 3, -7);
    A.set(A.length ~/ 2, -1);
    List<int> indexList = new List<int>();
    List<int> valueList = new List<int>();
    A.negativeValues(indexList, valueList);
    expect(2, equals(indexList.length));
    expect(2, equals(valueList.length));
    expect(indexList.contains(A.length / 3), isTrue);
    expect(indexList.contains(A.length / 2), isTrue);
    expect(valueList.contains(-7), isTrue);
    expect(valueList.contains(-1), isTrue);
  }

  void testNonZeros() {
    A.fill(0);
    A.set(A.length ~/ 3, 7);
    A.set(A.length ~/ 2, 1);
    List<int> indexList = new List<int>();
    List<int> valueList = new List<int>();
    A.nonZeros(indexList, valueList);
    expect(2, equals(indexList.length));
    expect(2, equals(valueList.length));
    expect(indexList.contains(A.length / 3), isTrue);
    expect(indexList.contains(A.length / 2), isTrue);
    expect(valueList.contains(7), isTrue);
    expect(valueList.contains(1), isTrue);
  }

  void testPositiveValues() {
    A.fill(0);
    A.set(A.length ~/ 3, 7);
    A.set(A.length ~/ 2, 1);
    List<int> indexList = new List<int>();
    List<int> valueList = new List<int>();
    A.positiveValues(indexList, valueList);
    expect(2, equals(indexList.length));
    expect(2, equals(valueList.length));
    expect(indexList.contains(A.length / 3), isTrue);
    expect(indexList.contains(A.length / 2), isTrue);
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
    IntMatrix2D B = A.reshape(rows, columns);
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
    IntMatrix1D Acopy = A.copy();
    IntMatrix1D Bcopy = B.copy();
    A.swap(B);
    for (int i = 0; i < A.length; i++) {
      expect(Bcopy.get(i), equals(A.get(i)));
      expect(Acopy.get(i), equals(B.get(i)));
    }
  }

  void testFlip() {
    IntMatrix1D b = A.flip();
    expect(A.length, b.length);
    for (int i = 0; i < A.length; i++) {
      expect(A.get(i), equals(b.get(A.length - 1 - i)));
    }
  }

  void testPart() {
    IntMatrix1D b = A.part(15, 11);
    for (int i = 0; i < 11; i++) {
      expect(A.get(15 + i), equals(b.get(i)));
    }
  }

  void testWhere() {
    IntMatrix1D b = A.where((int element) {
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
    IntMatrix1D b = A.select(indexes);
    for (int i = 0; i < indexes.length; i++) {
      expect(A.get(indexes[i]), equals(b.get(i)));
    }
  }

  /*void testSorted() {
    IntMatrix1D b = A.sorted();
    for (int i = 0; i < A.length - 1; i++) {
      expect(b.get(i + 1) >= b.get(i), isTrue);
    }
  }*/

  void testStrides() {
    int stride = 3;
    IntMatrix1D b = A.strides(stride);
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
    B.nonZeros(indexList, valueList);
    int product = A.dotNonZero(B, new Int32List.fromList(indexList), 5, B.length - 10);
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
