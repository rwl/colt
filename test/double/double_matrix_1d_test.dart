part of cern.colt.matrix.double.test;

testDoubleMatrix1D(String name, DoubleMatrix1DTest t) {
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
    test('equalsValue', t.testEqualsValue);
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
    test('zDotProduct', t.testZDotProduct);
    test('zDotProductRange', t.testZDotProductRange);
    test('zDotProductIndex', t.testZDotProductIndex);
  });
}

abstract class DoubleMatrix1DTest {
  /** Matrix to test. */
  DoubleMatrix1D A;

  /** Matrix of the same size as [A]. */
  DoubleMatrix1D B;

  double TOL = 1e-10;

  final int SIZE = 2 * 17 * 5;

  void populateMatrices() {
    //ConcurrencyUtils.setThreadsBeginN_1D(1);

    final r = new math.Random();

    for (int i = 0; i < A.size(); i++) {
      A.setQuick(i, r.nextDouble());
    }

    for (int i = 0; i < B.size(); i++) {
      B.setQuick(i, r.nextDouble());
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
    for (int i = 0; i < A.size(); i++) {
      double elem = A.getQuick(i);
      expected += elem * elem;
    }
    double result = A.aggregate(plus, square);
    expect(result, closeTo(expected, TOL));
  }

  testAggregateIndex() {
    List<int> indexList = new List<int>();
    for (int i = 0; i < A.size(); i++) {
      indexList.add(i);
    }
    double expected = 0.0;
    for (int i = 0; i < A.size(); i++) {
      double elem = A.getQuick(i);
      expected += elem * elem;
    }
    double result = A.aggregateIndex(plus, square,
        new Int32List.fromList(indexList));
    expect(result, closeTo(expected, TOL));
  }

  testAggregateMatrix() {
    double expected = 0.0;
    for (int i = 0; i < A.size(); i++) {
      double elemA = A.getQuick(i);
      double elemB = B.getQuick(i);
      expected += elemA * elemB;
    }
    double result = A.aggregateMatrix(B, plus, mult);
    expect(result, closeTo(expected, TOL));
  }

  testAssignValue() {
    double value = random.nextDouble();
    A.assignValue(value);
    for (int i = 0; i < A.size(); i++) {
      expect(value, closeTo(A.getQuick(i), TOL));
    }
  }

  testAssignValues() {
    Float64List expected = new Float64List(A.size());
    for (int i = 0; i < A.size(); i++) {
      expected[i] = random.nextDouble();
    }
    A.assignValues(expected);
    for (int i = 0; i < A.size(); i++) {
      expect(expected[i], closeTo(A.getQuick(i), TOL));
    }
  }

  testAssign() {
    DoubleMatrix1D Acopy = A.copy();
    A.assign(acos);
    for (int i = 0; i < A.size(); i++) {
      double expected = math.acos(Acopy.getQuick(i));
      expect(expected, closeTo(A.getQuick(i), TOL));
    }
  }

  testAssignMatrix() {
    A.assignMatrix(B);
    expect(A.size() == B.size(), isTrue);
    for (int i = 0; i < A.size(); i++) {
      expect(B.getQuick(i), closeTo(A.getQuick(i), TOL));
    }
  }

  testAssignFunc() {
    DoubleMatrix1D Acopy = A.copy();
    A.assignFunc(B, div);
    for (int i = 0; i < A.size(); i++) {
      expect(Acopy.getQuick(i) / B.getQuick(i), closeTo(A.getQuick(i), TOL));
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
    DoubleMatrix1D Acopy = A.copy();
    A.assignProcValue(procedure, -1.0);
    for (int i = 0; i < A.size(); i++) {
      if (Acopy.getQuick(i).abs() > 0.1) {
        expect(-1.0, closeTo(A.getQuick(i), TOL));
      } else {
        expect(Acopy.getQuick(i), closeTo(A.getQuick(i), TOL));
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
    DoubleMatrix1D Acopy = A.copy();
    A.assignProc(procedure, tan);
    for (int i = 0; i < A.size(); i++) {
      if (Acopy.getQuick(i).abs() > 0.1) {
        expect(math.tan(Acopy.getQuick(i)), closeTo(A.getQuick(i), TOL));
      } else {
        expect(Acopy.getQuick(i), closeTo(A.getQuick(i), TOL));
      }
    }
  }

  testCardinality() {
    int card = A.cardinality();
    expect(A.size(), equals(card));
  }

  testEqualsValue() {
    double value = 1.0;
    A.assignValue(value);
    bool eq = A.equalsValue(value);
    expect(eq, isTrue);
    eq = A.equals(2);
    expect(eq, isFalse);
  }

  testEquals() {
    bool eq = A.equals(A);
    expect(eq, isTrue);
    eq = A.equals(B);
    expect(eq, isFalse);
  }

  testMaxLocation() {
    A.assignValue(0.0);
    A.setQuick(A.size() ~/ 3, 0.7);
    A.setQuick(A.size() ~/ 2, 0.1);
    List<num> maxAndLoc = A.getMaxLocation();
    expect(0.7, closeTo(maxAndLoc[0], TOL));
    expect(A.size() ~/ 3, equals(maxAndLoc[1]));
  }

  testMinLocation() {
    A.assignValue(0.0);
    A.setQuick(A.size() ~/ 3, -0.7);
    A.setQuick(A.size() ~/ 2, -0.1);
    List<num> minAndLoc = A.getMinLocation();
    expect(-0.7, closeTo(minAndLoc[0], TOL));
    expect(A.size() ~/ 3, equals(minAndLoc[1]));
  }

  testGetNegativeValues() {
    A.assignValue(0.0);
    A.setQuick(A.size() ~/ 3, -0.7);
    A.setQuick(A.size() ~/ 2, -0.1);
    List<int> indexList = new List<int>();
    List<double> valueList = new List<double>();
    A.getNegativeValues(indexList, valueList);
    expect(2, indexList.length);
    expect(2, valueList.length);
    expect(indexList.contains(A.size() ~/ 3), isTrue);
    expect(indexList.contains(A.size() ~/ 2), isTrue);
    expect(valueList.contains(-0.7), isTrue);
    expect(valueList.contains(-0.1), isTrue);
  }

  testGetNonZeros() {
    A.assignValue(0.0);
    A.setQuick(A.size() ~/ 3, 0.7);
    A.setQuick(A.size() ~/ 2, 0.1);
    List<int> indexList = new List<int>();
    List<double> valueList = new List<double>();
    A.getNonZeros(indexList, valueList);
    expect(2, indexList.length);
    expect(2, valueList.length);
    expect(indexList.contains(A.size() ~/ 3), isTrue);
    expect(indexList.contains(A.size() ~/ 2), isTrue);
    expect(valueList.contains(0.7), isTrue);
    expect(valueList.contains(0.1), isTrue);
  }

  testGetPositiveValues() {
    A.assignValue(0.0);
    A.setQuick(A.size() ~/ 3, 0.7);
    A.setQuick(A.size() ~/ 2, 0.1);
    List<int> indexList = new List<int>();
    List<double> valueList = new List<double>();
    A.getPositiveValues(indexList, valueList);
    expect(2, indexList.length);
    expect(2, valueList.length);
    expect(indexList.contains(A.size() ~/ 3), isTrue);
    expect(indexList.contains(A.size() ~/ 2), isTrue);
    expect(valueList.contains(0.7), isTrue);
    expect(valueList.contains(0.1), isTrue);
  }

  testToArray() {
    Float64List array = A.toArray();
    expect(A.size() == array.length, isTrue);
    for (int i = 0; i < A.size(); i++) {
      expect(array[i], closeTo(A.getQuick(i), TOL));
    }
  }

  testToArrayFill() {
    Float64List array = new Float64List(A.size());
    A.toArrayValues(array);
    for (int i = 0; i < A.size(); i++) {
      expect(A.getQuick(i), closeTo(array[i], TOL));
    }
  }

  testReshape() {
    int rows = 10;
    int columns = 17;
    DoubleMatrix2D B = A.reshape(rows, columns);
    int idx = 0;
    for (int c = 0; c < columns; c++) {
      for (int r = 0; r < rows; r++) {
        expect(A.getQuick(idx++), closeTo(B.getQuick(r, c), TOL));
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
    DoubleMatrix1D Acopy = A.copy();
    DoubleMatrix1D Bcopy = B.copy();
    A.swap(B);
    for (int i = 0; i < A.size(); i++) {
      expect(Bcopy.getQuick(i), closeTo(A.getQuick(i), TOL));
      expect(Acopy.getQuick(i), closeTo(B.getQuick(i), TOL));
    }
  }

  testViewFlip() {
    DoubleMatrix1D b = A.viewFlip();
    expect(A.size(), equals(b.size()));
    for (int i = 0; i < A.size(); i++) {
      expect(A.getQuick(i), closeTo(b.getQuick(A.size() - 1 - i), TOL));
    }
  }

  testViewPart() {
    DoubleMatrix1D b = A.viewPart(15, 11);
    for (int i = 0; i < 11; i++) {
      expect(A.getQuick(15 + i), closeTo(b.getQuick(i), TOL));
    }
  }

  testViewSelection() {
    DoubleMatrix1D b = A.viewSelection((double element) {
      return element % 2 == 0;
    });
    for (int i = 0; i < b.size(); i++) {
      double el = b.getQuick(i);
      if (el % 2 != 0) {
        fail("viewSelection");
      }
    }
  }

  testViewSelectionIndex() {
    Int32List indexes = new Int32List.fromList([5, 11, 22, 37, 101]);
    DoubleMatrix1D b = A.viewSelectionIndex(indexes);
    for (int i = 0; i < indexes.length; i++) {
      expect(A.getQuick(indexes[i]), closeTo(b.getQuick(i), TOL));
    }
  }

  /*testViewSorted() {
    DoubleMatrix1D b = A.viewSorted();
    for (int i = 0; i < A.size() - 1; i++) {
      expect(b.getQuick(i + 1) >= b.getQuick(i), isTrue);
    }
  }*/

  testViewStrides() {
    int stride = 3;
    DoubleMatrix1D b = A.viewStrides(stride);
    for (int i = 0; i < b.size(); i++) {
      expect(A.getQuick(i * stride), closeTo(b.getQuick(i), TOL));
    }
  }

  testZDotProduct() {
    double product = A.zDotProduct(B);
    double expected = 0.0;
    for (int i = 0; i < A.size(); i++) {
      expected += A.getQuick(i) * B.getQuick(i);
    }
    expect(expected, closeTo(product, TOL));
  }

  testZDotProductRange() {
    double product = A.zDotProductRange(B, 5, B.size() - 10);
    double expected = 0.0;
    for (int i = 5; i < A.size() - 5; i++) {
      expected += A.getQuick(i) * B.getQuick(i);
    }
    expect(expected, closeTo(product, TOL));
  }

  //@Test
  testZDotProductIndex() {
    List<int> indexList = new List<int>();
    List<double> valueList = new List<double>();
    B.getNonZeros(indexList, valueList);
    double product = A.zDotProductIndex(B, 5, B.size() - 10, new Int32List.fromList(indexList));
    double expected = 0.0;
    for (int i = 5; i < A.size() - 5; i++) {
      expected += A.getQuick(i) * B.getQuick(i);
    }
    expect(expected, closeTo(product, TOL));
  }

  testZSum() {
    double sum = A.zSum();
    double expected = 0.0;
    for (int i = 0; i < A.size(); i++) {
      expected += A.getQuick(i);
    }
    expect(expected, closeTo(sum, TOL));
  }

}
