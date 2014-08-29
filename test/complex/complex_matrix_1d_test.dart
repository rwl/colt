part of cern.colt.matrix.complex.test;

abstract class DComplexMatrix1DTest {

  /** Matrix to test. */
  DComplexMatrix1D A;

  /** Matrix of the same size as [A]. */
  DComplexMatrix1D B;

  int SIZE = 2 * 17 * 5;

  double TOL = 1e-10;

  void setUp() {
    createMatrices();
    populateMatrices();
  }

  void createMatrices();

  void populateMatrices() {
    //ConcurrencyUtils.setThreadsBeginN_1D(1);

    for (int i = 0; i < A.size(); i++) {
      A.setQuick(i, new Float64List.fromList([random.nextDouble(), random.nextDouble()]));
    }

    for (int i = 0; i < A.size(); i++) {
      B.setQuick(i, new Float64List.fromList([random.nextDouble(), random.nextDouble()]));
    }
  }

  void tearDown() {
    A = B = null;
  }

  void testAggregateDoubleDoubleFunctionDoubleFunction() {
    Float64List expected = new Float64List(2);
    for (int i = 0; i < A.size(); i++) {
      expected = DComplex.plus(expected, DComplex.square(A.getQuick(i)));
    }
    Float64List result = A.aggregate(plus, square);
    assertEquals(expected, result, TOL);
  }

  void testAggregateComplexMatrix1DComplexComplexFunctionComplexComplexFunction() {
    Float64List actual = A.aggregateMatrix(B, plus, mult);
    Float64List expected = new Float64List(2);
    for (int i = 0; i < A.size(); i++) {
      expected = DComplex.plus(expected, DComplex.multiply(A.getQuick(i), B.getQuick(i)));
    }
    assertEquals(expected, actual, TOL);
  }

  void testAssignComplexComplexFunction() {
    DComplexMatrix1D Acopy = A.copy();
    A.assign(acos);
    for (int i = 0; i < A.size(); i++) {
      Float64List expected = DComplex.acos(Acopy.getQuick(i));
      assertEquals(expected, A.getQuick(i), TOL);
    }
  }

  void testAssignComplexMatrix1D() {
    A.assignMatrix(B);
    expect(A.size() == B.size(), isTrue);
    for (int i = 0; i < A.size(); i++) {
      assertEquals(B.getQuick(i), A.getQuick(i), TOL);
    }
  }

  void testAssignComplexMatrix1DComplexComplexFunction() {
    DComplexMatrix1D Acopy = A.copy();
    A.assignMatrixFunc(B, div);
    for (int i = 0; i < A.size(); i++) {
      assertEquals(DComplex.div_(Acopy.getQuick(i), B.getQuick(i)), A.getQuick(i), TOL);
    }
  }

  void testAssignComplexProcedureComplexComplexFunction() {
    bool procedure(Float64List element) {
      if (DComplex.abs(element) > 0.1) {
        return true;
      } else {
        return false;
      }
    }
    ;
    DComplexMatrix1D Acopy = A.copy();
    A.assignProc(procedure, tan);
    for (int i = 0; i < A.size(); i++) {
      if (DComplex.abs(Acopy.getQuick(i)) > 0.1) {
        assertEquals(DComplex.tan(Acopy.getQuick(i)), A.getQuick(i), TOL);
      } else {
        assertEquals(Acopy.getQuick(i), A.getQuick(i), TOL);
      }
    }
  }

  void testAssignComplexProcedureDoubleArray() {
    bool procedure(Float64List element) {
      if (DComplex.abs(element) > 0.1) {
        return true;
      } else {
        return false;
      }
    }
    DComplexMatrix1D Acopy = A.copy();
    Float64List value = new Float64List.fromList([-1, -1]);
    A.assignProcValue(procedure, value);
    for (int i = 0; i < A.size(); i++) {
      if (DComplex.abs(Acopy.getQuick(i)) > 0.1) {
        assertEquals(value, A.getQuick(i), TOL);
      } else {
        assertEquals(Acopy.getQuick(i), A.getQuick(i), TOL);
      }
    }
  }

  void testAssignComplexRealFunction() {
    DComplexMatrix1D Acopy = A.copy();
    A.assignRealFunc(abs);
    for (int i = 0; i < A.size(); i++) {
      Float64List elem = A.getQuick(i);
      expect(DComplex.abs(Acopy.getQuick(i)), closeTo(elem[0], TOL));
      expect(0, closeTo(elem[1], TOL));
    }
  }

  void testAssignDoubleArray() {
    Float64List expected = new Float64List(2 * A.size());
    for (int i = 0; i < 2 * A.size(); i++) {
      expected[i] = random.nextDouble();
    }
    A.assignValues(expected);
    for (int i = 0; i < A.size(); i++) {
      Float64List elem = A.getQuick(i);
      expect(expected[2 * i], closeTo(elem[0], TOL));
      expect(expected[2 * i + 1], closeTo(elem[1], TOL));
    }

  }

  void testAssignDoubleDouble() {
    double re = random.nextDouble();
    double im = random.nextDouble();
    A.assignValue(re, im);
    for (int i = 0; i < A.size(); i++) {
      Float64List elem = A.getQuick(i);
      expect(re, closeTo(elem[0], TOL));
      expect(im, closeTo(elem[1], TOL));
    }
  }

  void testAssignImaginary() {
    DComplexMatrix1D Acopy = A.copy();
    DoubleMatrix1D Im = DoubleFactory1D.dense.random(A.size());
    A.assignImaginary(Im);
    for (int i = 0; i < A.size(); i++) {
      Float64List elem = A.getQuick(i);
      expect(Acopy.getQuick(i)[0], closeTo(elem[0], TOL));
      expect(Im.getQuick(i), closeTo(elem[1], TOL));
    }
  }

  void testAssignReal() {
    DComplexMatrix1D Acopy = A.copy();
    DoubleMatrix1D Re = DoubleFactory1D.dense.random(A.size());
    A.assignReal(Re);
    for (int i = 0; i < A.size(); i++) {
      Float64List elem = A.getQuick(i);
      expect(Acopy.getQuick(i)[1], closeTo(elem[1], TOL));
      expect(Re.getQuick(i), closeTo(elem[0], TOL));
    }
  }

  void testCardinality() {
    int card = A.cardinality();
    expect(A.size(), equals(card));
  }

  void testEqualsDouble() {
    Float64List value = new Float64List.fromList([1, 2]);
    A.assignValue(value[0], value[1]);
    bool eq = A.equals(value);
    expect(eq, isTrue);
    eq = A.equals(new Float64List.fromList([2, 2]));
    expect(eq, isFalse);
  }

  void testEqualsObject() {
    bool eq = A.equals(A);
    expect(eq, isTrue);
    eq = A.equals(B);
    expect(eq, isFalse);
  }

  void testGetImaginaryPart() {
    DoubleMatrix1D Im = A.getImaginaryPart();
    for (int i = 0; i < A.size(); i++) {
      expect(A.getQuick(i)[1], closeTo(Im.getQuick(i), TOL));
    }
  }

  void testGetRealPart() {
    DoubleMatrix1D Re = A.getRealPart();
    for (int i = 0; i < A.size(); i++) {
      expect(A.getQuick(i)[0], closeTo(Re.getQuick(i), TOL));
    }
  }

  void testGetNonZerosIntArrayListArrayListOfdouble() {
    List<int> indexList = new List<int>();
    List<Float64List> valueList = new List<Float64List>();
    A.getNonZeros(indexList, valueList);
    expect(A.size(), equals(indexList.length));
    expect(A.size(), equals(valueList.length));
    for (int i = 0; i < A.size(); i++) {
      assertEquals(A.getQuick(indexList[i]), valueList[i], TOL);
      expect(valueList[i][0] != 0 || valueList[i][1] != 0, isTrue);
    }
  }

  void testReshapeIntInt() {
    int rows = 10;
    int columns = 17;
    DComplexMatrix2D B = A.reshape(rows, columns);
    int idx = 0;
    for (int c = 0; c < columns; c++) {
      for (int r = 0; r < rows; r++) {
        assertEquals(A.getQuick(idx++), B.getQuick(r, c), TOL);
      }
    }
  }

  /*void testReshape3D() {
    int slices = 2;
    int rows = 5;
    int columns = 17;
    DComplexMatrix3D B = A.reshape3D(slices, rows, columns);
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
    DComplexMatrix1D Acopy = A.copy();
    DComplexMatrix1D Bcopy = B.copy();
    A.swap(B);
    for (int i = 0; i < A.size(); i++) {
      assertEquals(Bcopy.getQuick(i), A.getQuick(i), TOL);
      assertEquals(Acopy.getQuick(i), B.getQuick(i), TOL);
    }
  }

  void testToArray() {
    Float64List array = A.toArray();
    for (int i = 0; i < A.size(); i++) {
      Float64List elem = A.getQuick(i);
      expect(elem[0], closeTo(array[2 * i], TOL));
      expect(elem[1], closeTo(array[2 * i + 1], TOL));
    }
  }

  void testToArrayDoubleArray() {
    Float64List array = new Float64List(2 * A.size());
    A.toArrayFill(array);
    for (int i = 0; i < A.size(); i++) {
      Float64List elem = A.getQuick(i);
      expect(elem[0], closeTo(array[2 * i], TOL));
      expect(elem[1], closeTo(array[2 * i + 1], TOL));
    }
  }

  void testViewFlip() {
    DComplexMatrix1D B = A.viewFlip();
    for (int i = 0; i < A.size(); i++) {
      assertEquals(A.getQuick(A.size() - 1 - i), B.getQuick(i), TOL);
    }
  }

  void testViewPart() {
    DComplexMatrix1D B = A.viewPart(A.size() ~/ 2, A.size() ~/ 3);
    for (int i = 0; i < A.size() / 3; i++) {
      assertEquals(A.getQuick(A.size() ~/ 2 + i), B.getQuick(i), TOL);
    }
  }

  void testViewSelectionComplexProcedure() {
    DComplexMatrix1D B = A.viewSelectionProc((Float64List element) {
      if (element[0] < element[1]) {
        return true;
      } else {
        return false;
      }
    });
    for (int i = 0; i < B.size(); i++) {
      Float64List el = B.getQuick(i);
      if (el[0] >= el[1]) {
        fail('viewSelectionProc');
      }
    }
  }

  void testViewSelectionIntArray() {
    Int32List indexes = new Int32List.fromList([A.size() / 6, A.size() / 5, A.size() / 4, A.size() / 3, A.size() / 2]);
    DComplexMatrix1D B = A.viewSelection(indexes);
    for (int i = 0; i < indexes.length; i++) {
      assertEquals(A.getQuick(indexes[i]), B.getQuick(i), TOL);
    }
  }

  void testViewStrides() {
    int stride = 3;
    DComplexMatrix1D B = A.viewStrides(stride);
    for (int i = 0; i < B.size(); i++) {
      assertEquals(A.getQuick(i * stride), B.getQuick(i), TOL);
    }
  }

  void testZDotProductComplexMatrix1D() {
    Float64List actual = A.zDotProduct(B);
    Float64List expected = new Float64List(2);
    for (int i = 0; i < A.size(); i++) {
      expected = DComplex.plus(expected, DComplex.multiply(DComplex.conj(B.getQuick(i)), A.getQuick(i)));
    }
    assertEquals(expected, actual, TOL);
  }

  void testZDotProductComplexMatrix1DIntInt() {
    Float64List actual = A.zDotProduct(B, 5, B.size() - 10);
    Float64List expected = new Float64List(2);
    for (int i = 5; i < A.size() - 5; i++) {
      expected = DComplex.plus(expected, DComplex.multiply(DComplex.conj(B.getQuick(i)), A.getQuick(i)));
    }
    assertEquals(expected, actual, TOL);
  }

  void testZDotProductComplexMatrix1DIntIntIntArrayList() {
    List<int> indexList = new List<int>();
    List<Float64List> valueList = new List<Float64List>();
    B.getNonZeros(indexList, valueList);
    Float64List actual = A.zDotProductIndex(B, indexList, 5, B.size() - 10);
    Float64List expected = new Float64List(2);
    for (int i = 5; i < A.size() - 5; i++) {
      expected = DComplex.plus(expected, DComplex.multiply(A.getQuick(i), DComplex.conj(B.getQuick(i))));
    }
    assertEquals(expected, actual, TOL);
  }

  void testZSum() {
    Float64List actual = A.zSum();
    Float64List expected = new Float64List(2);
    for (int i = 0; i < A.size(); i++) {
      expected = DComplex.plus(expected, A.getQuick(i));
    }
    assertEquals(expected, actual, TOL);
  }

  void assertEquals(Float64List expected, Float64List actual, double tol) {
    for (int i = 0; i < actual.length; i++) {
      expect(expected[i], closeTo(actual[i], tol));
    }
  }

}
