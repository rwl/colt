part of cern.colt.matrix.complex.test;

testComplexVector(String name, ComplexVectorTest t) {
  group(name, () {
    setUp(t.setUp);
    tearDown(t.tearDown);
    test('aggregate', t.testAggregate);
    test('aggregateMatrix', t.testAggregateMatrix);
    test('assign', t.testAssign);
    test('assignMatrix', t.testAssignMatrix);
    test('assignMatrixFunc', t.testAssignMatrixFunc);
    test('assignProc', t.testAssignProc);
    test('assignProcValue', t.testAssignProcValue);
    test('assignRealFunc', t.testAssignRealFunc);
    test('assignValues', t.testAssignValues);
    test('assignValue', t.testAssignValue);
    test('assignImaginary', t.testAssignImaginary);
    test('assignReal', t.testAssignReal);
    test('cardinality', t.testCardinality);
    test('equalsValue', t.testEqualsValue);
    test('equals', t.testEquals);
    test('getImaginaryPart', t.testGetImaginaryPart);
    test('getRealPart', t.testGetRealPart);
    test('getNonZeros', t.testGetNonZeros);
    test('reshape', t.testReshape);
    test('swap', t.testSwap);
    test('toArray', t.testToArray);
    test('toArrayFill', t.testToArrayFill);
    test('viewFlip', t.testViewFlip);
    test('viewPart', t.testViewPart);
    test('viewSelectionProc', t.testViewSelectionProc);
    test('viewSelection', t.testViewSelection);
    test('viewStrides', t.testViewStrides);
    test('zDotProduct', t.testZDotProduct);
    test('zDotProductRange', t.testZDotProductRange);
    test('zDotProductIndex', t.testZDotProductIndex);
    test('zSum', t.testZSum);
  });
}

abstract class ComplexVectorTest {

  /** Matrix to test. */
  AbstractComplexVector A;

  /** Matrix of the same size as [A]. */
  AbstractComplexVector B;

  int SIZE = 2 * 17 * 5;

  double TOL = 1e-10;

  void setUp() {
    createMatrices();
    populateMatrices();
  }

  void createMatrices();

  void populateMatrices() {
    //ConcurrencyUtils.setThreadsBeginN_1D(1);

    for (int i = 0; i < A.length; i++) {
      A.set(i, new Float64List.fromList([random.nextDouble(), random.nextDouble()]));
    }

    for (int i = 0; i < A.length; i++) {
      B.set(i, new Float64List.fromList([random.nextDouble(), random.nextDouble()]));
    }
  }

  void tearDown() {
    A = B = null;
  }

  void testAggregate() {
    Float64List expected = new Float64List(2);
    for (int i = 0; i < A.length; i++) {
      expected = Complex.plus(expected, Complex.square(A.get(i)));
    }
    Float64List result = A.reduce(plus, square);
    assertEquals(expected, result, TOL);
  }

  void testAggregateMatrix() {
    Float64List actual = A.reduceVector(B, plus, mult);
    Float64List expected = new Float64List(2);
    for (int i = 0; i < A.length; i++) {
      expected = Complex.plus(expected, Complex.multiply(A.get(i), B.get(i)));
    }
    assertEquals(expected, actual, TOL);
  }

  void testAssign() {
    AbstractComplexVector Acopy = A.copy();
    A.forEach(acos);
    for (int i = 0; i < A.length; i++) {
      Float64List expected = Complex.acos(Acopy.get(i));
      assertEquals(expected, A.get(i), TOL);
    }
  }

  void testAssignMatrix() {
    A.copyFrom(B);
    expect(A.length == B.length, isTrue);
    for (int i = 0; i < A.length; i++) {
      assertEquals(B.get(i), A.get(i), TOL);
    }
  }

  void testAssignMatrixFunc() {
    AbstractComplexVector Acopy = A.copy();
    A.forEachWith(B, div);
    for (int i = 0; i < A.length; i++) {
      assertEquals(Complex.div_(Acopy.get(i), B.get(i)), A.get(i), TOL);
    }
  }

  void testAssignProc() {
    bool procedure(Float64List element) {
      if (Complex.abs(element) > 0.1) {
        return true;
      } else {
        return false;
      }
    }
    AbstractComplexVector Acopy = A.copy();
    A.forEachWhere(procedure, tan);
    for (int i = 0; i < A.length; i++) {
      if (Complex.abs(Acopy.get(i)) > 0.1) {
        assertEquals(Complex.tan(Acopy.get(i)), A.get(i), TOL);
      } else {
        assertEquals(Acopy.get(i), A.get(i), TOL);
      }
    }
  }

  void testAssignProcValue() {
    bool procedure(Float64List element) {
      if (Complex.abs(element) > 0.1) {
        return true;
      } else {
        return false;
      }
    }
    AbstractComplexVector Acopy = A.copy();
    Float64List value = new Float64List.fromList([-1.0, -1.0]);
    A.fillWhere(procedure, value);
    for (int i = 0; i < A.length; i++) {
      if (Complex.abs(Acopy.get(i)) > 0.1) {
        assertEquals(value, A.get(i), TOL);
      } else {
        assertEquals(Acopy.get(i), A.get(i), TOL);
      }
    }
  }

  void testAssignRealFunc() {
    AbstractComplexVector Acopy = A.copy();
    A.forEachReal(abs);
    for (int i = 0; i < A.length; i++) {
      Float64List elem = A.get(i);
      expect(Complex.abs(Acopy.get(i)), closeTo(elem[0], TOL));
      expect(0, closeTo(elem[1], TOL));
    }
  }

  void testAssignValues() {
    Float64List expected = new Float64List(2 * A.length);
    for (int i = 0; i < 2 * A.length; i++) {
      expected[i] = random.nextDouble();
    }
    A.setAll(expected);
    for (int i = 0; i < A.length; i++) {
      Float64List elem = A.get(i);
      expect(expected[2 * i], closeTo(elem[0], TOL));
      expect(expected[2 * i + 1], closeTo(elem[1], TOL));
    }

  }

  void testAssignValue() {
    double re = random.nextDouble();
    double im = random.nextDouble();
    A.fill(re, im);
    for (int i = 0; i < A.length; i++) {
      Float64List elem = A.get(i);
      expect(re, closeTo(elem[0], TOL));
      expect(im, closeTo(elem[1], TOL));
    }
  }

  void testAssignImaginary() {
    AbstractComplexVector Acopy = A.copy();
    AbstractDoubleVector Im = DoubleFactory1D.dense.random(A.length);
    A.setImaginary(Im);
    for (int i = 0; i < A.length; i++) {
      Float64List elem = A.get(i);
      expect(Acopy.get(i)[0], closeTo(elem[0], TOL));
      expect(Im.get(i), closeTo(elem[1], TOL));
    }
  }

  void testAssignReal() {
    AbstractComplexVector Acopy = A.copy();
    AbstractDoubleVector Re = DoubleFactory1D.dense.random(A.length);
    A.setReal(Re);
    for (int i = 0; i < A.length; i++) {
      Float64List elem = A.get(i);
      expect(Acopy.get(i)[1], closeTo(elem[1], TOL));
      expect(Re.get(i), closeTo(elem[0], TOL));
    }
  }

  void testCardinality() {
    int card = A.cardinality;
    expect(A.length, equals(card));
  }

  void testEqualsValue() {
    Float64List value = new Float64List.fromList([1.0, 2.0]);
    A.fill(value[0], value[1]);
    expect(A == value, isTrue);
    final eq = A == new Float64List.fromList([2.0, 2.0]);
    expect(eq, isFalse);
  }

  void testEquals() {
    expect(A == A, isTrue);
    expect(A == B, isFalse);
  }

  void testGetImaginaryPart() {
    AbstractDoubleVector Im = A.imaginary();
    for (int i = 0; i < A.length; i++) {
      expect(A.get(i)[1], closeTo(Im.get(i), TOL));
    }
  }

  void testGetRealPart() {
    AbstractDoubleVector Re = A.real();
    for (int i = 0; i < A.length; i++) {
      expect(A.get(i)[0], closeTo(Re.get(i), TOL));
    }
  }

  void testGetNonZeros() {
    List<int> indexList = new List<int>();
    List<Float64List> valueList = new List<Float64List>();
    A.nonZeros(indexList, valueList);
    expect(A.length, equals(indexList.length));
    expect(A.length, equals(valueList.length));
    for (int i = 0; i < A.length; i++) {
      assertEquals(A.get(indexList[i]), valueList[i], TOL);
      expect(valueList[i][0] != 0 || valueList[i][1] != 0, isTrue);
    }
  }

  void testReshape() {
    int rows = 10;
    int columns = 17;
    AbstractComplexMatrix B = A.reshape(rows, columns);
    int idx = 0;
    for (int c = 0; c < columns; c++) {
      for (int r = 0; r < rows; r++) {
        assertEquals(A.get(idx++), B.get(r, c), TOL);
      }
    }
  }

  /*void testReshape3D() {
    int slices = 2;
    int rows = 5;
    int columns = 17;
    ComplexMatrix3D B = A.reshape3D(slices, rows, columns);
    int idx = 0;
    for (int s = 0; s < slices; s++) {
      for (int c = 0; c < columns; c++) {
        for (int r = 0; r < rows; r++) {
          assertEquals(A.getQuick(idx++), B.getQuick(s, r, c), TOL);
        }
      }
    }
  }*/

  void testSwap() {
    AbstractComplexVector Acopy = A.copy();
    AbstractComplexVector Bcopy = B.copy();
    A.swap(B);
    for (int i = 0; i < A.length; i++) {
      assertEquals(Bcopy.get(i), A.get(i), TOL);
      assertEquals(Acopy.get(i), B.get(i), TOL);
    }
  }

  void testToArray() {
    Float64List array = A.toList();
    for (int i = 0; i < A.length; i++) {
      Float64List elem = A.get(i);
      expect(elem[0], closeTo(array[2 * i], TOL));
      expect(elem[1], closeTo(array[2 * i + 1], TOL));
    }
  }

  void testToArrayFill() {
    Float64List array = new Float64List(2 * A.length);
    A.fillList(array);
    for (int i = 0; i < A.length; i++) {
      Float64List elem = A.get(i);
      expect(elem[0], closeTo(array[2 * i], TOL));
      expect(elem[1], closeTo(array[2 * i + 1], TOL));
    }
  }

  void testViewFlip() {
    AbstractComplexVector B = A.flip();
    for (int i = 0; i < A.length; i++) {
      assertEquals(A.get(A.length - 1 - i), B.get(i), TOL);
    }
  }

  void testViewPart() {
    AbstractComplexVector B = A.part(A.length ~/ 2, A.length ~/ 3);
    for (int i = 0; i < A.length / 3; i++) {
      assertEquals(A.get(A.length ~/ 2 + i), B.get(i), TOL);
    }
  }

  void testViewSelectionProc() {
    AbstractComplexVector B = A.where((Float64List element) {
      if (element[0] < element[1]) {
        return true;
      } else {
        return false;
      }
    });
    for (int i = 0; i < B.length; i++) {
      Float64List el = B.get(i);
      if (el[0] >= el[1]) {
        fail('viewSelectionProc');
      }
    }
  }

  void testViewSelection() {
    Int32List indexes = new Int32List.fromList([A.length ~/ 6, A.length ~/ 5, A.length ~/ 4, A.length ~/ 3, A.length ~/ 2]);
    AbstractComplexVector B = A.select(indexes);
    for (int i = 0; i < indexes.length; i++) {
      assertEquals(A.get(indexes[i]), B.get(i), TOL);
    }
  }

  void testViewStrides() {
    int stride = 3;
    AbstractComplexVector B = A.strides(stride);
    for (int i = 0; i < B.length; i++) {
      assertEquals(A.get(i * stride), B.get(i), TOL);
    }
  }

  void testZDotProduct() {
    Float64List actual = A.dot(B);
    Float64List expected = new Float64List(2);
    for (int i = 0; i < A.length; i++) {
      expected = Complex.plus(expected, Complex.multiply(Complex.conj(B.get(i)), A.get(i)));
    }
    assertEquals(expected, actual, TOL);
  }

  void testZDotProductRange() {
    Float64List actual = A.dot(B, 5, B.length - 10);
    Float64List expected = new Float64List(2);
    for (int i = 5; i < A.length - 5; i++) {
      expected = Complex.plus(expected, Complex.multiply(Complex.conj(B.get(i)), A.get(i)));
    }
    assertEquals(expected, actual, TOL);
  }

  void testZDotProductIndex() {
    List<int> indexList = new List<int>();
    List<Float64List> valueList = new List<Float64List>();
    B.nonZeros(indexList, valueList);
    Float64List actual = A.dotNonZero(B,
        new Int32List.fromList(indexList), 5, B.length - 10);
    Float64List expected = new Float64List(2);
    for (int i = 5; i < A.length - 5; i++) {
      expected = Complex.plus(expected, Complex.multiply(A.get(i), Complex.conj(B.get(i))));
    }
    assertEquals(expected, actual, TOL);
  }

  void testZSum() {
    Float64List actual = A.sum();
    Float64List expected = new Float64List(2);
    for (int i = 0; i < A.length; i++) {
      expected = Complex.plus(expected, A.get(i));
    }
    assertEquals(expected, actual, TOL);
  }

  void assertEquals(Float64List expected, Float64List actual, double tol) {
    for (int i = 0; i < actual.length; i++) {
      expect(expected[i], closeTo(actual[i], tol));
    }
  }

}
