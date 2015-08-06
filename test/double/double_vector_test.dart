part of cern.colt.matrix.double.test;

testDoubleVector(String name, DoubleVectorTest t) {
  group(name, () {
    setUp(t.setUp);
    tearDown(t.tearDown);
    test('aggregate', t.testAggregate);
    test('fill', t.testFill);
    test('setAll', t.testSetAll);
    test('apply', t.testApply);
    test('copyFrom', t.testCopyFrom);
    test('assign', t.testAssign);
    test('cardinality', t.testCardinality);
    test('equals', t.testEquals);
    test('maxLocation', t.testMax);
    test('minLocation', t.testMin);
    test('negative', t.testNegative);
    test('nonzero', t.testNonzero);
    test('positive', t.testPositive);
    test('toList', t.testToList);
    test('fillList', t.testFillList);
    test('reshape', t.testReshape);
    test('flip', t.testFlip);
    test('part', t.testPart);
    test('select', t.testSelect);
    test('strides', t.testStrides);
    test('dot', t.testDot);
    test('dot range', t.testDotRange);
  });
}

abstract class DoubleVectorTest {
  /// Matrix to test.
  AbstractDoubleVector A;

  /// Matrix of the same size as [A].
  AbstractDoubleVector B;

  double TOL = 1e-10;

  final int SIZE = 2 * 17 * 5;

  void populateMatrices() {
    for (int i = 0; i < A.size; i++) {
      A.set(i, random.nextDouble());
    }

    for (int i = 0; i < B.size; i++) {
      B.set(i, random.nextDouble());
    }
  }

  void createMatrices();

  void setUp() {
    createMatrices();
    populateMatrices();
  }

  void tearDown() {
    A = null;
    B = null;
  }

  testAggregate() {
    double expected = 0.0;
    for (int i = 0; i < A.size; i++) {
      double elem = A.get(i);
      expected += elem * elem;
    }
    double result = A.aggregate(plus, square);
    expect(result, closeTo(expected, TOL));
  }

  testFill() {
    double value = random.nextDouble();
    A.fill(value);
    for (int i = 0; i < A.size; i++) {
      expect(value, closeTo(A.get(i), TOL));
    }
  }

  testSetAll() {
    Float64List expected = new Float64List(A.size);
    for (int i = 0; i < A.size; i++) {
      expected[i] = random.nextDouble();
    }
    A.setAll(expected);
    for (int i = 0; i < A.size; i++) {
      expect(expected[i], closeTo(A.get(i), TOL));
    }
  }

  testApply() {
    AbstractDoubleVector Acopy = A.copy();
    A.apply(acos);
    for (int i = 0; i < A.size; i++) {
      double expected = math.acos(Acopy.get(i));
      expect(expected, closeTo(A.get(i), TOL));
    }
  }

  testCopyFrom() {
    A.copyFrom(B);
    expect(A.size, equals(B.size));
    for (int i = 0; i < A.size; i++) {
      expect(B.get(i), closeTo(A.get(i), TOL));
    }
  }

  testAssign() {
    AbstractDoubleVector Acopy = A.copy();
    A.assign(B, div);
    for (int i = 0; i < A.size; i++) {
      expect(Acopy.get(i) / B.get(i), closeTo(A.get(i), TOL));
    }
  }

  testCardinality() {
    int card = A.cardinality;
    expect(A.size, equals(card));
  }

  testEquals() {
    expect(A.equals(A), isTrue);
    expect(A.equals(B), isFalse);
  }

  testMax() {
    A.fill(0.0);
    A.set(A.size ~/ 3, 0.7);
    A.set(A.size ~/ 2, 0.1);
    final maxAndLoc = A.max();
    expect(0.7, closeTo(maxAndLoc.value, TOL));
    expect(A.size ~/ 3, equals(maxAndLoc.location));
  }

  testMin() {
    A.fill(0.0);
    A.set(A.size ~/ 3, -0.7);
    A.set(A.size ~/ 2, -0.1);
    final minAndLoc = A.min();
    expect(-0.7, closeTo(minAndLoc.value, TOL));
    expect(A.size ~/ 3, equals(minAndLoc.location));
  }

  testNegative() {
    A.fill(0.0);
    A.set(A.size ~/ 3, -0.7);
    A.set(A.size ~/ 2, -0.1);
    List<int> indexList = new List<int>();
    List<double> valueList = new List<double>();
    A.negative(indexList: indexList, valueList: valueList);
    expect(2, indexList.length);
    expect(2, valueList.length);
    expect(indexList.contains(A.size ~/ 3), isTrue);
    expect(indexList.contains(A.size ~/ 2), isTrue);
    expect(valueList.contains(-0.7), isTrue);
    expect(valueList.contains(-0.1), isTrue);
  }

  testNonzero() {
    A.fill(0.0);
    A.set(A.size ~/ 3, 0.7);
    A.set(A.size ~/ 2, 0.1);
    List<int> indexList = new List<int>();
    List<double> valueList = new List<double>();
    A.nonzero(indexList: indexList, valueList: valueList);
    expect(2, indexList.length);
    expect(2, valueList.length);
    expect(indexList.contains(A.size ~/ 3), isTrue);
    expect(indexList.contains(A.size ~/ 2), isTrue);
    expect(valueList.contains(0.7), isTrue);
    expect(valueList.contains(0.1), isTrue);
  }

  testPositive() {
    A.fill(0.0);
    A.set(A.size ~/ 3, 0.7);
    A.set(A.size ~/ 2, 0.1);
    List<int> indexList = new List<int>();
    List<double> valueList = new List<double>();
    A.positive(indexList: indexList, valueList: valueList);
    expect(2, indexList.length);
    expect(2, valueList.length);
    expect(indexList.contains(A.size ~/ 3), isTrue);
    expect(indexList.contains(A.size ~/ 2), isTrue);
    expect(valueList.contains(0.7), isTrue);
    expect(valueList.contains(0.1), isTrue);
  }

  testToList() {
    Float64List array = A.toList();
    expect(A.size == array.length, isTrue);
    for (int i = 0; i < A.size; i++) {
      expect(array[i], closeTo(A.get(i), TOL));
    }
  }

  testFillList() {
    var array = new Float64List(A.size);
    A.fillList(array);
    for (int i = 0; i < A.size; i++) {
      expect(A.get(i), closeTo(array[i], TOL));
    }
  }

  testReshape() {
    int rows = 10;
    int columns = 17;
    AbstractDoubleMatrix B = A.reshape(rows, columns);
    int idx = 0;
    for (int c = 0; c < columns; c++) {
      for (int r = 0; r < rows; r++) {
        expect(A.get(idx++), closeTo(B.get(r, c), TOL));
      }
    }
  }

  testFlip() {
    AbstractDoubleVector b = A.flip();
    expect(A.size, equals(b.size));
    for (int i = 0; i < A.size; i++) {
      expect(A.get(i), closeTo(b.get(A.size - 1 - i), TOL));
    }
  }

  testPart() {
    AbstractDoubleVector b = A.part(15, 11);
    for (int i = 0; i < 11; i++) {
      expect(A.get(15 + i), closeTo(b.get(i), TOL));
    }
  }

  testSelect() {
    var indexes = new Int32List.fromList([5, 11, 22, 37, 101]);
    AbstractDoubleVector b = A.select(indexes);
    for (int i = 0; i < indexes.length; i++) {
      expect(A.get(indexes[i]), closeTo(b.get(i), TOL));
    }
  }

  testStrides() {
    int stride = 3;
    AbstractDoubleVector b = A.strides(stride);
    for (int i = 0; i < b.size; i++) {
      expect(A.get(i * stride), closeTo(b.get(i), TOL));
    }
  }

  testDot() {
    double product = A.dot(B);
    double expected = 0.0;
    for (int i = 0; i < A.size; i++) {
      expected += A.get(i) * B.get(i);
    }
    expect(expected, closeTo(product, TOL));
  }

  testDotRange() {
    double product = A.dot(B, 5, B.size - 10);
    double expected = 0.0;
    for (int i = 5; i < A.size - 5; i++) {
      expected += A.get(i) * B.get(i);
    }
    expect(expected, closeTo(product, TOL));
  }

  testSum() {
    double sum = A.sum();
    double expected = 0.0;
    for (int i = 0; i < A.size; i++) {
      expected += A.get(i);
    }
    expect(expected, closeTo(sum, TOL));
  }
}
