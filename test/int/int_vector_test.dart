part of cern.colt.matrix.int.test;

const SIZE = 2 * 17 * 5;

testAbstractIntVector(String kind, AbstractIntVector make(int size)) {
  group("AbstractIntVector ($kind)", () {
    AbstractIntVector A, B;
    setUp(() {
      A = make(SIZE);
      B = make(SIZE);
      for (int i = 0; i < A.size; i++) {
        A.set(i, math.max(1, random.nextInt(MAX_INT) % A.size));
        B.set(i, math.max(1, random.nextInt(MAX_INT) % A.size));
      }
    });
    test('aggregate', () {
      int expected = 0;
      for (int i = 0; i < A.size; i++) {
        int elem = A.get(i);
        expected += elem * elem;
      }
      int result = A.aggregate(ifunc.plus, ifunc.square);
      expect(expected, equals(result));
    });

    test('fill', () {
      int value = random.nextInt(MAX_INT);
      A.fill(value);
      for (int i = 0; i < A.size; i++) {
        expect(value, equals(A.get(i)));
      }
    });

    test('setAll', () {
      var expected = new Int32List(A.size);
      for (int i = 0; i < A.size; i++) {
        expected[i] = random.nextInt(MAX_INT);
      }
      A.setAll(expected);
      for (int i = 0; i < A.size; i++) {
        expect(expected[i], A.get(i));
      }
    });

    test('apply', () {
      AbstractIntVector Acopy = A.copy();
      A.apply(ifunc.neg);
      for (int i = 0; i < A.size; i++) {
        int expected = -Acopy.get(i);
        expect(expected, equals(A.get(i)));
      }
    });

    test('copyFrom', () {
      A.copyFrom(B);
      expect(A.size == B.size, isTrue);
      for (int i = 0; i < A.size; i++) {
        expect(B.get(i), equals(A.get(i)));
      }
    });

    test('assign', () {
      AbstractIntVector Acopy = A.copy();
      A.assign(B, ifunc.plus);
      for (int i = 0; i < A.size; i++) {
        expect(Acopy.get(i) + B.get(i), equals(A.get(i)));
      }
    });

    test('cardinality', () {
      int card = A.cardinality;
      int expected = 0;
      for (int i = 0; i < A.size; i++) {
        if (A.get(i) != 0) expected++;
      }
      expect(expected, equals(card));
    });

    test('equals', () {
      expect(A.equals(A), isTrue);
      expect(A.equals(B), isFalse);
    });

    test('max', () {
      A.fill(0);
      A.set(A.size ~/ 3, 7);
      A.set(A.size ~/ 2, 1);
      final maxAndLoc = A.max();
      expect(7, equals(maxAndLoc.value));
      expect(A.size ~/ 3, equals(maxAndLoc.location));
    });

    test('min', () {
      A.fill(0);
      A.set(A.size ~/ 3, -7);
      A.set(A.size ~/ 2, -1);
      final minAndLoc = A.min();
      expect(-7, equals(minAndLoc.value));
      expect(A.size ~/ 3, equals(minAndLoc.location));
    });

    test('negative', () {
      A.fill(0);
      A.set(A.size ~/ 3, -7);
      A.set(A.size ~/ 2, -1);
      List<int> indexList = new List<int>();
      List<int> valueList = new List<int>();
      A.negative(indexList: indexList, valueList: valueList);
      expect(2, equals(indexList.length));
      expect(2, equals(valueList.length));
      expect(indexList.contains(A.size ~/ 3), isTrue);
      expect(indexList.contains(A.size ~/ 2), isTrue);
      expect(valueList.contains(-7), isTrue);
      expect(valueList.contains(-1), isTrue);
    });

    test('nonzero', () {
      A.fill(0);
      A.set(A.size ~/ 3, 7);
      A.set(A.size ~/ 2, 1);
      List<int> indexList = new List<int>();
      List<int> valueList = new List<int>();
      A.nonzero(indexList: indexList, valueList: valueList);
      expect(2, equals(indexList.length));
      expect(2, equals(valueList.length));
      expect(indexList.contains(A.size ~/ 3), isTrue);
      expect(indexList.contains(A.size ~/ 2), isTrue);
      expect(valueList.contains(7), isTrue);
      expect(valueList.contains(1), isTrue);
    });

    test('positive', () {
      A.fill(0);
      A.set(A.size ~/ 3, 7);
      A.set(A.size ~/ 2, 1);
      List<int> indexList = new List<int>();
      List<int> valueList = new List<int>();
      A.positive(indexList: indexList, valueList: valueList);
      expect(2, equals(indexList.length));
      expect(2, equals(valueList.length));
      expect(indexList.contains(A.size ~/ 3), isTrue);
      expect(indexList.contains(A.size ~/ 2), isTrue);
      expect(valueList.contains(7), isTrue);
      expect(valueList.contains(1), isTrue);
    });

    test('toList', () {
      final array = A.toList();
      expect(A.size == array.length, isTrue);
      for (int i = 0; i < A.size; i++) {
        expect(array[i], equals(A.get(i)));
      }
    });

    test('fillList', () {
      final array = new Int32List(A.size);
      A.fillList(array);
      for (int i = 0; i < A.size; i++) {
        expect(A.get(i), equals(array[i]));
      }
    });

    test('reshape', () {
      int rows = 10;
      int columns = 17;
      AbstractIntMatrix B = A.reshape(rows, columns);
      int idx = 0;
      for (int c = 0; c < columns; c++) {
        for (int r = 0; r < rows; r++) {
          expect(A.get(idx++), equals(B.get(r, c)));
        }
      }
    });

    test('flip', () {
      AbstractIntVector b = A.flip();
      expect(A.size, b.size);
      for (int i = 0; i < A.size; i++) {
        expect(A.get(i), equals(b.get(A.size - 1 - i)));
      }
    });

    test('part', () {
      AbstractIntVector b = A.part(15, 11);
      for (int i = 0; i < 11; i++) {
        expect(A.get(15 + i), equals(b.get(i)));
      }
    });

    test('select', () {
      final indexes = new Int32List.fromList([5, 11, 22, 37, 101]);
      AbstractIntVector b = A.select(indexes);
      for (int i = 0; i < indexes.length; i++) {
        expect(A.get(indexes[i]), equals(b.get(i)));
      }
    });

    test('strides', () {
      int stride = 3;
      AbstractIntVector b = A.strides(stride);
      for (int i = 0; i < b.size; i++) {
        expect(A.get(i * stride), equals(b.get(i)));
      }
    });

    test('dot', () {
      int product = A.dot(B);
      int expected = 0;
      for (int i = 0; i < A.size; i++) {
        expected += A.get(i) * B.get(i);
      }
      expect(expected, equals(product));
    });

    test('dot range', () {
      int product = A.dot(B, 5, B.size - 10);
      int expected = 0;
      for (int i = 5; i < A.size - 5; i++) {
        expected += A.get(i) * B.get(i);
      }
      expect(expected, equals(product));
    });

    test('sum', () {
      int sum = A.sum();
      int expected = 0;
      for (int i = 0; i < A.size; i++) {
        expected += A.get(i);
      }
      expect(expected, equals(sum));
    });
  });
}
