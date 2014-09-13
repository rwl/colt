part of cern.colt.matrix.double.test;

testDoubleVector(String name, DoubleVectorTest t) {
  group(name, () {
    setUp(t.setUp);
    tearDown(t.tearDown);
    test('aggregate', t.testAggregate);
    test('aggregateIndex', t.testAggregateIndex);
    test('aggregateMatrix', t.testAggregateMatrix);
    test('assignValue', t.testAssignValue);
    test('assignValues', t.testAssignValues);
    test('assign', t.testAssign);
    test('assignMatrix', t.testAssignMatrix);
    test('assignFunc', t.testAssignFunc);
    test('assignProcValue', t.testAssignProceValue);
    test('assignProc', t.testAssignProc);
    test('cardinality', t.testCardinality);
    test('all', t.testAll);
    test('equals', t.testEquals);
    test('maxLocation', t.testMaxLocation);
    test('minLocation', t.testMinLocation);
    test('getNegativeValues', t.testGetNegativeValues);
    test('getNonZeros', t.testGetNonZeros);
    test('getPositiveValues', t.testGetPositiveValues);
    test('toArray', t.testToArray);
    test('toArrayFill', t.testToArrayFill);
    test('reshape', t.testReshape);
    test('swap', t.testSwap);
    test('viewFlip', t.testViewFlip);
    test('viewPart', t.testViewPart);
    test('viewSelection', t.testViewSelection);
    test('viewSelectionIndex', t.testViewSelectionIndex);
    test('viewStrides', t.testViewStrides);
    test('zDotProduct', t.testDot);
    test('zDotProductRange', t.testDotRange);
    test('zDotProductIndex', t.testDotNonZero);
  });
}

abstract class DoubleVectorTest {
  /** Matrix to test. */
  AbstractDoubleVector A;

  /** Matrix of the same size as [A]. */
  AbstractDoubleVector B;

  double TOL = 1e-10;

  final int SIZE = 2 * 17 * 5;

  void populateMatrices() {
    //ConcurrencyUtils.setThreadsBeginN_1D(1);

    final r = new math.Random();

    for (int i = 0; i < A.length; i++) {
      A.set(i, r.nextDouble());
    }

    for (int i = 0; i < B.length; i++) {
      B.set(i, r.nextDouble());
    }
  }

  void createMatrices();

  void setUp() {//throws Exception {
    createMatrices();
    populateMatrices();
  }

  void tearDown() {//throws Exception {
    A = null;
    B = null;
  }

  testAggregate() {
    double expected = 0.0;
    for (int i = 0; i < A.length; i++) {
      double elem = A.get(i);
      expected += elem * elem;
    }
    double result = A.aggregate(plus, square);
    expect(result, closeTo(expected, TOL));
  }

  testAggregateIndex() {
    List<int> indexList = new List<int>();
    for (int i = 0; i < A.length; i++) {
      indexList.add(i);
    }
    double expected = 0.0;
    for (int i = 0; i < A.length; i++) {
      double elem = A.get(i);
      expected += elem * elem;
    }
    double result = A.reduceRange(plus, square,
        new Int32List.fromList(indexList));
    expect(result, closeTo(expected, TOL));
  }

  testAggregateMatrix() {
    double expected = 0.0;
    for (int i = 0; i < A.length; i++) {
      double elemA = A.get(i);
      double elemB = B.get(i);
      expected += elemA * elemB;
    }
    double result = A.reduceWith(B, plus, mult);
    expect(result, closeTo(expected, TOL));
  }

  testAssignValue() {
    double value = random.nextDouble();
    A.fill(value);
    for (int i = 0; i < A.length; i++) {
      expect(value, closeTo(A.get(i), TOL));
    }
  }

  testAssignValues() {
    Float64List expected = new Float64List(A.length);
    for (int i = 0; i < A.length; i++) {
      expected[i] = random.nextDouble();
    }
    A.setValues(expected);
    for (int i = 0; i < A.length; i++) {
      expect(expected[i], closeTo(A.get(i), TOL));
    }
  }

  testAssign() {
    AbstractDoubleVector Acopy = A.copy();
    A.forEach(acos);
    for (int i = 0; i < A.length; i++) {
      double expected = math.acos(Acopy.get(i));
      expect(expected, closeTo(A.get(i), TOL));
    }
  }

  testAssignMatrix() {
    A.copyFrom(B);
    expect(A.length == B.length, isTrue);
    for (int i = 0; i < A.length; i++) {
      expect(B.get(i), closeTo(A.get(i), TOL));
    }
  }

  testAssignFunc() {
    AbstractDoubleVector Acopy = A.copy();
    A.forEachWith(B, div);
    for (int i = 0; i < A.length; i++) {
      expect(Acopy.get(i) / B.get(i), closeTo(A.get(i), TOL));
    }
  }

  testAssignProceValue() {
    bool procedure(double element) {
      if (element.abs() > 0.1) {
        return true;
      } else {
        return false;
      }
    }
    AbstractDoubleVector Acopy = A.copy();
    A.fillWhere(procedure, -1.0);
    for (int i = 0; i < A.length; i++) {
      if (Acopy.get(i).abs() > 0.1) {
        expect(-1.0, closeTo(A.get(i), TOL));
      } else {
        expect(Acopy.get(i), closeTo(A.get(i), TOL));
      }
    }
  }

  testAssignProc() {
    bool procedure(double element) {
      if (element.abs() > 0.1) {
        return true;
      } else {
        return false;
      }
    }
    AbstractDoubleVector Acopy = A.copy();
    A.forEachWhere(procedure, tan);
    for (int i = 0; i < A.length; i++) {
      if (Acopy.get(i).abs() > 0.1) {
        expect(math.tan(Acopy.get(i)), closeTo(A.get(i), TOL));
      } else {
        expect(Acopy.get(i), closeTo(A.get(i), TOL));
      }
    }
  }

  testCardinality() {
    int card = A.cardinality;
    expect(A.length, equals(card));
  }

  testAll() {
    double value = 1.0;
    A.fill(value);
    expect(A.all(value), isTrue);
    expect(A.all(2.0), isFalse);
  }

  testEquals() {
    expect(A.equals(A), isTrue);
    expect(A.equals(B), isFalse);
  }

  testMaxLocation() {
    A.fill(0.0);
    A.set(A.length ~/ 3, 0.7);
    A.set(A.length ~/ 2, 0.1);
    final maxAndLoc = A.max();
    expect(0.7, closeTo(maxAndLoc.value, TOL));
    expect(A.length ~/ 3, equals(maxAndLoc.location));
  }

  testMinLocation() {
    A.fill(0.0);
    A.set(A.length ~/ 3, -0.7);
    A.set(A.length ~/ 2, -0.1);
    final minAndLoc = A.min();
    expect(-0.7, closeTo(minAndLoc.value, TOL));
    expect(A.length ~/ 3, equals(minAndLoc.location));
  }

  testGetNegativeValues() {
    A.fill(0.0);
    A.set(A.length ~/ 3, -0.7);
    A.set(A.length ~/ 2, -0.1);
    List<int> indexList = new List<int>();
    List<double> valueList = new List<double>();
    A.negativeValues(indexList: indexList, valueList: valueList);
    expect(2, indexList.length);
    expect(2, valueList.length);
    expect(indexList.contains(A.length ~/ 3), isTrue);
    expect(indexList.contains(A.length ~/ 2), isTrue);
    expect(valueList.contains(-0.7), isTrue);
    expect(valueList.contains(-0.1), isTrue);
  }

  testGetNonZeros() {
    A.fill(0.0);
    A.set(A.length ~/ 3, 0.7);
    A.set(A.length ~/ 2, 0.1);
    List<int> indexList = new List<int>();
    List<double> valueList = new List<double>();
    A.nonZeros(indexList: indexList, valueList: valueList);
    expect(2, indexList.length);
    expect(2, valueList.length);
    expect(indexList.contains(A.length ~/ 3), isTrue);
    expect(indexList.contains(A.length ~/ 2), isTrue);
    expect(valueList.contains(0.7), isTrue);
    expect(valueList.contains(0.1), isTrue);
  }

  testGetPositiveValues() {
    A.fill(0.0);
    A.set(A.length ~/ 3, 0.7);
    A.set(A.length ~/ 2, 0.1);
    List<int> indexList = new List<int>();
    List<double> valueList = new List<double>();
    A.positiveValues(indexList: indexList, valueList: valueList);
    expect(2, indexList.length);
    expect(2, valueList.length);
    expect(indexList.contains(A.length ~/ 3), isTrue);
    expect(indexList.contains(A.length ~/ 2), isTrue);
    expect(valueList.contains(0.7), isTrue);
    expect(valueList.contains(0.1), isTrue);
  }

  testToArray() {
    Float64List array = A.toList();
    expect(A.length == array.length, isTrue);
    for (int i = 0; i < A.length; i++) {
      expect(array[i], closeTo(A.get(i), TOL));
    }
  }

  testToArrayFill() {
    Float64List array = new Float64List(A.length);
    A.fillList(array);
    for (int i = 0; i < A.length; i++) {
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

  /*void testReshapeIntIntInt() {
    int slices = 2;
    int rows = 5;
    int columns = 17;
    DoubleMatrix3D B = A.reshape(slices, rows, columns);
    int idx = 0;
    for (int s = 0; s < slices; s++) {
      for (int c = 0; c < columns; c++) {
        for (int r = 0; r < rows; r++) {
          expect(A.getQuick(idx++), closeTo(B.getQuick(s, r, c), TOL));
        }
      }
    }
  }*/

  testSwap() {
    AbstractDoubleVector Acopy = A.copy();
    AbstractDoubleVector Bcopy = B.copy();
    A.swap(B);
    for (int i = 0; i < A.length; i++) {
      expect(Bcopy.get(i), closeTo(A.get(i), TOL));
      expect(Acopy.get(i), closeTo(B.get(i), TOL));
    }
  }

  testViewFlip() {
    AbstractDoubleVector b = A.flip();
    expect(A.length, equals(b.length));
    for (int i = 0; i < A.length; i++) {
      expect(A.get(i), closeTo(b.get(A.length - 1 - i), TOL));
    }
  }

  testViewPart() {
    AbstractDoubleVector b = A.part(15, 11);
    for (int i = 0; i < 11; i++) {
      expect(A.get(15 + i), closeTo(b.get(i), TOL));
    }
  }

  testViewSelection() {
    AbstractDoubleVector b = A.where((double element) {
      return element % 2 == 0;
    });
    for (int i = 0; i < b.length; i++) {
      double el = b.get(i);
      if (el % 2 != 0) {
        fail("viewSelection");
      }
    }
  }

  testViewSelectionIndex() {
    final indexes = [5, 11, 22, 37, 101];
    AbstractDoubleVector b = A.select(indexes);
    for (int i = 0; i < indexes.length; i++) {
      expect(A.get(indexes[i]), closeTo(b.get(i), TOL));
    }
  }

  /*testViewSorted() {
    DoubleVector b = A.viewSorted();
    for (int i = 0; i < A.size() - 1; i++) {
      expect(b.getQuick(i + 1) >= b.getQuick(i), isTrue);
    }
  }*/

  testViewStrides() {
    int stride = 3;
    AbstractDoubleVector b = A.strides(stride);
    for (int i = 0; i < b.length; i++) {
      expect(A.get(i * stride), closeTo(b.get(i), TOL));
    }
  }

  testDot() {
    double product = A.dot(B);
    double expected = 0.0;
    for (int i = 0; i < A.length; i++) {
      expected += A.get(i) * B.get(i);
    }
    expect(expected, closeTo(product, TOL));
  }

  testDotRange() {
    double product = A.dot(B, 5, B.length - 10);
    double expected = 0.0;
    for (int i = 5; i < A.length - 5; i++) {
      expected += A.get(i) * B.get(i);
    }
    expect(expected, closeTo(product, TOL));
  }

  testDotNonZero() {
    List<int> indexList = new List<int>();
    List<double> valueList = new List<double>();
    B.nonZeros(indexList: indexList, valueList: valueList);
    double product = A.dotNonZero(B, indexList, 5, B.length - 10);
    double expected = 0.0;
    for (int i = 5; i < A.length - 5; i++) {
      expected += A.get(i) * B.get(i);
    }
    expect(expected, closeTo(product, TOL));
  }

  testSum() {
    double sum = A.sum();
    double expected = 0.0;
    for (int i = 0; i < A.length; i++) {
      expected += A.get(i);
    }
    expect(expected, closeTo(sum, TOL));
  }

}
