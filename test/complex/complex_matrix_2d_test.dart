part of cern.colt.matrix.complex.test;

testDComplexMatrix2D(String name, DComplexMatrix2DTest t) {
  group(name, () {
    setUp(t.setUp);
    tearDown(t.tearDown);
    test('aggregate', t.testAggregate);
    test('aggregateMatrix', t.testAggregateMatrix);
    test('assign', t.testAssign);
    test('assignMatrix', t.testAssignMatrix);
    test('assignFunc', t.testAssignFunc);
    test('assignProc', t.testAssignProc);
    test('assignProcValue', t.testProcValue);
    test('assignRealFunc', t.testRealFunc);
    test('assignValues', t.testAssignValues);
    test('assignList', t.testAssignList);
    test('assignValue', t.testAssignValue);
    test('assignImaginary', t.testAssignImaginary);
    test('assignReal', t.testAssignReal);
    test('cardinality', t.testCardinality);
    test('equalsValue', t.testEqualsValue);
    test('equals', t.testEquals);
    test('forEachNonZero', t.testForEachNonZero);
    test('getConjugateTranspose', t.testGetConjugateTranspose);
    test('getImaginaryPart', t.testGetImaginaryPart);
    test('getNonZeros', t.testGetNonZeros);
    test('getRealPart', t.testGetRealPart);
    test('toArray', t.testToArray);
    test('vectorize', t.testVectorize);
    test('viewColumn', t.testViewColumn);
    test('viewColumnFlip', t.testViewColumnFlip);
    test('viewDice', t.testViewDice);
    test('viewPart', t.testViewPart);
    test('viewRow', t.testViewRow);
    test('viewRowFlip', t.testViewRowFlip);
    test('viewSelectionProc', t.testViewSelectionProc);
    test('viewSelection', t.testViewSelection);
    test('viewStrides', t.testViewStrides);
    test('zMult', t.testZMult);
    test('zMult2D', t.testZMult2D);
    test('zSum', t.testZSum);
  });
}

abstract class DComplexMatrix2DTest {

  /** Matrix to test. */
  DComplexMatrix2D A;

  /** Matrix of the same size as [A]. */
  DComplexMatrix2D B;

  /** Matrix of the size `A.columns() x A.rows()`. */
  DComplexMatrix2D Bt;

  int NROWS = 13;

  int NCOLUMNS = 17;

  double TOL = 1e-10;

  void setUp() {
    createMatrices();
    populateMatrices();
  }

  void createMatrices();

  void populateMatrices() {
    //ConcurrencyUtils.setThreadsBeginN_2D(1);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        A.set(r, c, new Float64List.fromList([random.nextDouble(), random.nextDouble()]));
      }
    }

    for (int r = 0; r < B.rows; r++) {
      for (int c = 0; c < B.columns; c++) {
        B.set(r, c, new Float64List.fromList([random.nextDouble(), random.nextDouble()]));
      }
    }

    for (int r = 0; r < Bt.rows; r++) {
      for (int c = 0; c < Bt.columns; c++) {
        Bt.set(r, c, new Float64List.fromList([random.nextDouble(), random.nextDouble()]));
      }
    }
  }

  void tearDown() {
    A = B = Bt = null;
  }

  void testAggregate() {
    Float64List actual = A.reduce(plus, square);
    Float64List expected = new Float64List(2);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expected = DComplex.plus(expected, DComplex.square(A.get(r, c)));
      }
    }
    assertEquals(expected, actual, TOL);
  }

  void testAggregateMatrix() {
    Float64List actual = A.reduceMatrix(B, plus, mult);
    Float64List expected = new Float64List(2);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expected = DComplex.plus(expected, DComplex.multiply(A.get(r, c), B.get(r, c)));
      }
    }
    assertEquals(expected, actual, TOL);
  }

  void testAssign() {
    DComplexMatrix2D Acopy = A.copy();
    A.forEach(acos);
    Float64List tmp;
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        tmp = DComplex.acos(Acopy.get(r, c));
        assertEquals(tmp, A.get(r, c), TOL);
      }
    }
  }

  void testAssignMatrix() {
    A.copyFrom(B);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        assertEquals(B.get(r, c), A.get(r, c), TOL);
      }
    }
  }

  void testAssignFunc() {
    DComplexMatrix2D Acopy = A.copy();
    A.forEachMatrix(B, div);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        assertEquals(DComplex.div_(Acopy.get(r, c), B.get(r, c)), A.get(r, c), TOL);
      }
    }
  }

  void testAssignProc() {
    DComplexMatrix2D Acopy = A.copy();
    A.forEachWhere((Float64List element) {
      if (DComplex.abs(element) > 3) {
        return true;
      } else {
        return false;
      }
    }, tan);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        if (DComplex.abs(Acopy.get(r, c)) > 3) {
          assertEquals(DComplex.tan(Acopy.get(r, c)), A.get(r, c), TOL);
        } else {
          assertEquals(Acopy.get(r, c), A.get(r, c), TOL);
        }
      }
    }
  }

  void testProcValue() {
    DComplexMatrix2D Acopy = A.copy();
    Float64List value = new Float64List.fromList([random.nextDouble(), random.nextDouble()]);
    A.fillWhere((Float64List element) {
      if (DComplex.abs(element) > 3) {
        return true;
      } else {
        return false;
      }
    }, value);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        if (DComplex.abs(A.get(r, c)) > 3) {
          assertEquals(value, A.get(r, c), TOL);
        } else {
          assertEquals(Acopy.get(r, c), A.get(r, c), TOL);
        }
      }
    }
  }

  void testRealFunc() {
    DComplexMatrix2D Acopy = A.copy();
    A.forEachReal(abs);
    Float64List tmp;
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        tmp = A.get(r, c);
        expect(DComplex.abs(Acopy.get(r, c)), closeTo(tmp[0], TOL));
        expect(0, closeTo(tmp[1], TOL));
      }
    }
  }

  void testAssignValues() {
    Float64List expected = new Float64List(2 * A.length);
    for (int i = 0; i < 2 * A.length; i++) {
      expected[i] = random.nextDouble();
    }
    A.setAll(expected);
    int idx = 0;
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        Float64List elem = A.get(r, c);
        expect(expected[idx], closeTo(elem[0], TOL));
        expect(expected[idx + 1], closeTo(elem[1], TOL));
        idx += 2;
      }
    }
  }

  void testAssignList() {
    List<Float64List> expected = new List<Float64List>.generate(A.rows,
        (_) => new Float64List(2 * A.columns));
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < 2 * A.columns; c++) {
        expected[r][c] = random.nextDouble();
      }
    }
    A.setAll2D(expected);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        Float64List elem = A.get(r, c);
        expect(expected[r][2 * c], closeTo(elem[0], TOL));
        expect(expected[r][2 * c + 1], closeTo(elem[1], TOL));
      }
    }
  }

  void testAssignValue() {
    Float64List value = new Float64List.fromList([random.nextDouble(), random.nextDouble()]);
    A.fill(value[0], value[1]);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        Float64List elem = A.get(r, c);
        assertEquals(value, elem, TOL);
      }
    }
  }

  void testAssignImaginary() {
    DoubleMatrix2D Im = DoubleFactory2D.dense.random(A.rows, A.columns);
    DComplexMatrix2D Acopy = A.copy();
    A.setImaginary(Im);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(Acopy.get(r, c)[0], closeTo(A.get(r, c)[0], TOL));
        expect(Im.get(r, c), closeTo(A.get(r, c)[1], TOL));
      }
    }
  }

  void testAssignReal() {
    DoubleMatrix2D Re = DoubleFactory2D.dense.random(A.rows, A.columns);
    DComplexMatrix2D Acopy = A.copy();
    A.setReal(Re);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(Acopy.get(r, c)[1], closeTo(A.get(r, c)[1], TOL));
        expect(Re.get(r, c), closeTo(A.get(r, c)[0], TOL));
      }
    }
  }

  void testCardinality() {
    int card = A.cardinality;
    expect(A.length, equals(card));
  }

  void testEqualsValue() {
    Float64List value = new Float64List.fromList([random.nextDouble(), random.nextDouble()]);
    A.fill(value[0], value[1]);
    expect(A == value, isTrue);
    final eq = A == new Float64List.fromList([value[0] + 1, value[1] + 1]);
    expect(eq, isFalse);
  }

  void testEquals() {
    expect(A == A, isTrue);
    expect(A == B, isFalse);
  }

  void testForEachNonZero() {
    DComplexMatrix2D Acopy = A.copy();
    Float64List function(int first, int second, Float64List third) {
      return DComplex.sqrt(third);
    }
    A.forEachNonZero(function);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        assertEquals(DComplex.sqrt(Acopy.get(r, c)), A.get(r, c), TOL);
      }
    }
  }

  void testGetConjugateTranspose() {
    DComplexMatrix2D Aconj = A.conjugateTranspose();
    expect(A.rows, equals(Aconj.columns));
    expect(A.columns, equals(Aconj.rows));
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(A.get(r, c)[0], closeTo(Aconj.get(c, r)[0], TOL));
        expect(-A.get(r, c)[1], closeTo(Aconj.get(c, r)[1], TOL));
      }
    }
  }

  void testGetImaginaryPart() {
    DoubleMatrix2D Im = A.imaginary();
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(A.get(r, c)[1], closeTo(Im.get(r, c), TOL));
      }
    }
  }

  void testGetNonZeros() {
    List<int> rowList = new List<int>();
    List<int> colList = new List<int>();
    List<Float64List> valueList = new List<Float64List>();
    A.nonZeros(rowList, colList, valueList);
    expect(A.length, equals(rowList.length));
    expect(A.length, equals(colList.length));
    expect(A.length, equals(valueList.length));
    int idx = 0;
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        assertEquals(A.get(rowList[idx], colList[idx]), valueList[idx], TOL);
        idx++;
      }
    }
  }

  void testGetRealPart() {
    DoubleMatrix2D Re = A.real();
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(A.get(r, c)[0], closeTo(Re.get(r, c), TOL));
      }
    }
  }

  void testToArray() {
    List<Float64List> array = A.toList();
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(A.get(r, c)[0], closeTo(array[r][2 * c], TOL));
        expect(A.get(r, c)[1], closeTo(array[r][2 * c + 1], TOL));
      }
    }
  }

  void testVectorize() {
    DComplexMatrix1D B = A.vectorize();
    int idx = 0;
    for (int c = 0; c < A.columns; c++) {
      for (int r = 0; r < A.rows; r++) {
        assertEquals(A.get(r, c), B.get(idx++), TOL);
      }
    }
  }

  void testViewColumn() {
    DComplexMatrix1D B = A.column(A.columns ~/ 2);
    expect(A.rows, equals(B.length));
    for (int r = 0; r < A.rows; r++) {
      assertEquals(A.get(r, A.columns ~/ 2), B.get(r), TOL);
    }
  }

  void testViewColumnFlip() {
    DComplexMatrix2D B = A.columnFlip();
    expect(A.length, equals(B.length));
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        assertEquals(A.get(r, A.columns - 1 - c), B.get(r, c), TOL);
      }
    }
  }

  void testViewDice() {
    DComplexMatrix2D B = A.dice();
    expect(A.rows, equals(B.columns));
    expect(A.columns, equals(B.rows));
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        assertEquals(A.get(r, c), B.get(c, r), TOL);
      }
    }
  }

  void testViewPart() {
    DComplexMatrix2D B = A.part(A.rows ~/ 2, A.columns ~/ 2, A.rows ~/ 3, A.columns ~/ 3);
    for (int r = 0; r < A.rows / 3; r++) {
      for (int c = 0; c < A.columns / 3; c++) {
        assertEquals(A.get(A.rows ~/ 2 + r, A.columns ~/ 2 + c), B.get(r, c), TOL);
      }
    }
  }

  void testViewRow() {
    DComplexMatrix1D B = A.row(A.rows ~/ 2);
    expect(A.columns, equals(B.length));
    for (int c = 0; c < A.columns; c++) {
      assertEquals(A.get(A.rows ~/ 2, c), B.get(c), TOL);
    }
  }

  void testViewRowFlip() {
    DComplexMatrix2D B = A.rowFlip();
    expect(A.length, equals(B.length));
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        assertEquals(A.get(A.rows - 1 - r, c), B.get(r, c), TOL);
      }
    }
  }

  void testViewSelectionProc() {
    final Float64List value = new Float64List.fromList([random.nextDouble(), random.nextDouble()]);
    A.set(A.rows ~/ 3, 0, value);
    A.set(A.rows ~/ 2, 0, value);
    DComplexMatrix2D B = A.where((DComplexMatrix1D element) {
      return DComplex.isEqual(element.get(0), value, TOL);
    });
    expect(2, equals(B.rows));
    expect(A.columns, equals(B.columns));
    assertEquals(A.get(A.rows ~/ 3, 0), B.get(0, 0), TOL);
    assertEquals(A.get(A.rows ~/ 2, 0), B.get(1, 0), TOL);
  }

  void testViewSelection() {
    Int32List rowIndexes = new Int32List.fromList([A.rows ~/ 6, A.rows ~/ 5, A.rows ~/ 4, A.rows ~/ 3, A.rows ~/ 2]);
    Int32List colIndexes = new Int32List.fromList([A.columns ~/ 6, A.columns ~/ 5, A.columns ~/ 4, A.columns ~/ 3, A.columns ~/ 2, A.columns - 1]);
    DComplexMatrix2D B = A.select(rowIndexes, colIndexes);
    expect(rowIndexes.length, equals(B.rows));
    expect(colIndexes.length, equals(B.columns));
    for (int r = 0; r < rowIndexes.length; r++) {
      for (int c = 0; c < colIndexes.length; c++) {
        assertEquals(A.get(rowIndexes[r], colIndexes[c]), B.get(r, c), TOL);
      }
    }
  }

  void testViewStrides() {
    int rowStride = 3;
    int colStride = 5;
    DComplexMatrix2D B = A.strides(rowStride, colStride);
    for (int r = 0; r < B.rows; r++) {
      for (int c = 0; c < B.columns; c++) {
        assertEquals(A.get(r * rowStride, c * colStride), B.get(r, c), TOL);
      }
    }
  }

  void testZMult() {
    DComplexMatrix1D y = new DenseDComplexMatrix1D(A.columns);
    for (int i = 0; i < y.length; i++) {
      y.set(i, new Float64List.fromList([random.nextDouble(), random.nextDouble()]));
    }
    Float64List alpha = new Float64List.fromList([3.0, 2.0]);
    Float64List beta = new Float64List.fromList([5.0, 4.0]);
    DComplexMatrix1D z = null;
    z = A.mult(y, z, alpha, beta, false);
    Float64List expected = new Float64List(2 * A.rows);
    Float64List tmp = new Float64List(2);
    for (int r = 0; r < A.rows; r++) {
      Float64List s = new Float64List(2);
      for (int c = 0; c < A.columns; c++) {
        s = DComplex.plus(s, DComplex.multiply(A.get(r, c), y.get(c)));
      }
      tmp[0] = expected[2 * r];
      tmp[1] = expected[2 * r + 1];
      tmp = DComplex.multiply(beta, tmp);
      tmp = DComplex.plus(tmp, DComplex.multiply(alpha, s));
      expected[2 * r] = tmp[0];
      expected[2 * r + 1] = tmp[1];
    }

    for (int r = 0; r < A.rows; r++) {
      expect(expected[2 * r], closeTo(z.get(r)[0], TOL));
      expect(expected[2 * r + 1], closeTo(z.get(r)[1], TOL));
    }
    //transpose
    y = new DenseDComplexMatrix1D(A.rows);
    for (int i = 0; i < y.length; i++) {
      y.set(i, new Float64List.fromList([random.nextDouble(), random.nextDouble()]));
    }
    z = null;
    z = A.mult(y, z, alpha, beta, true);
    expected = new Float64List(2 * A.columns);
    for (int r = 0; r < A.columns; r++) {
      Float64List s = new Float64List(2);
      for (int c = 0; c < A.rows; c++) {
        s = DComplex.plus(s, DComplex.multiply(DComplex.conj(A.get(c, r)), y.get(c)));
      }
      tmp[0] = expected[2 * r];
      tmp[1] = expected[2 * r + 1];
      tmp = DComplex.multiply(beta, tmp);
      tmp = DComplex.plus(tmp, DComplex.multiply(alpha, s));
      expected[2 * r] = tmp[0];
      expected[2 * r + 1] = tmp[1];
    }
    for (int r = 0; r < A.columns; r++) {
      expect(expected[2 * r], closeTo(z.get(r)[0], TOL));
      expect(expected[2 * r + 1], closeTo(z.get(r)[1], TOL));
    }
  }

  void testZMult2D() {
    Float64List alpha = new Float64List.fromList([3.0, 2.0]);
    Float64List beta = new Float64List.fromList([5.0, 4.0]);
    Float64List tmp = new Float64List(2);
    DComplexMatrix2D C = null;
    C = A.multiply(Bt, C, alpha, beta, false, false);
    List<Float64List> expected = new List<Float64List>.generate(A.rows,
        (_) => new Float64List(2 * A.rows));
    for (int j = 0; j < A.rows; j++) {
      for (int i = 0; i < A.rows; i++) {
        Float64List s = new Float64List(2);
        for (int k = 0; k < A.columns; k++) {
          s = DComplex.plus(s, DComplex.multiply(A.get(i, k), Bt.get(k, j)));
        }
        tmp[0] = expected[i][2 * j];
        tmp[1] = expected[i][2 * j + 1];
        tmp = DComplex.multiply(tmp, beta);
        tmp = DComplex.plus(tmp, DComplex.multiply(s, alpha));
        expected[i][2 * j] = tmp[0];
        expected[i][2 * j + 1] = tmp[1];
      }
    }
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.rows; c++) {
        expect(expected[r][2 * c], closeTo(C.get(r, c)[0], TOL));
        expect(expected[r][2 * c + 1], closeTo(C.get(r, c)[1], TOL));
      }
    }

    //transposeA
    C = null;
    C = A.multiply(B, C, alpha, beta, true, false);
    expected = new List<Float64List>.generate(A.columns,
        (_) => new Float64List(2 * A.columns));
    for (int j = 0; j < A.columns; j++) {
      for (int i = 0; i < A.columns; i++) {
        Float64List s = new Float64List(2);
        for (int k = 0; k < A.rows; k++) {
          s = DComplex.plus(s, DComplex.multiply(DComplex.conj(A.get(k, i)), B.get(k, j)));
        }
        tmp[0] = expected[i][2 * j];
        tmp[1] = expected[i][2 * j + 1];
        tmp = DComplex.multiply(tmp, beta);
        tmp = DComplex.plus(tmp, DComplex.multiply(s, alpha));
        expected[i][2 * j] = tmp[0];
        expected[i][2 * j + 1] = tmp[1];
      }
    }
    for (int r = 0; r < A.columns; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(expected[r][2 * c], closeTo(C.get(r, c)[0], TOL));
        expect(expected[r][2 * c + 1], closeTo(C.get(r, c)[1], TOL));
      }
    }
    //transposeB
    C = null;
    C = A.multiply(B, C, alpha, beta, false, true);
    expected = new List<Float64List>.generate(A.rows,
        (_) => new Float64List(2 * A.rows));
    for (int j = 0; j < A.rows; j++) {
      for (int i = 0; i < A.rows; i++) {
        Float64List s = new Float64List(2);
        for (int k = 0; k < A.columns; k++) {
          s = DComplex.plus(s, DComplex.multiply(A.get(i, k), DComplex.conj(B.get(j, k))));
        }
        tmp[0] = expected[i][2 * j];
        tmp[1] = expected[i][2 * j + 1];
        tmp = DComplex.multiply(tmp, beta);
        tmp = DComplex.plus(tmp, DComplex.multiply(s, alpha));
        expected[i][2 * j] = tmp[0];
        expected[i][2 * j + 1] = tmp[1];
      }
    }
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.rows; c++) {
        expect(expected[r][2 * c], closeTo(C.get(r, c)[0], TOL));
        expect(expected[r][2 * c + 1], closeTo(C.get(r, c)[1], TOL));
      }
    }
    //transposeA and transposeB
    C = null;
    C = A.multiply(Bt, C, alpha, beta, true, true);
    expected = new List<Float64List>.generate(A.columns,
        (_) => new Float64List(2 * A.columns));
    for (int j = 0; j < A.columns; j++) {
      for (int i = 0; i < A.columns; i++) {
        Float64List s = new Float64List(2);
        for (int k = 0; k < A.rows; k++) {
          s = DComplex.plus(s, DComplex.multiply(DComplex.conj(A.get(k, i)), DComplex.conj(Bt.get(j, k))));
        }
        tmp[0] = expected[i][2 * j];
        tmp[1] = expected[i][2 * j + 1];
        tmp = DComplex.multiply(tmp, beta);
        tmp = DComplex.plus(tmp, DComplex.multiply(s, alpha));
        expected[i][2 * j] = tmp[0];
        expected[i][2 * j + 1] = tmp[1];
      }
    }
    for (int r = 0; r < A.columns; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(expected[r][2 * c], closeTo(C.get(r, c)[0], TOL));
        expect(expected[r][2 * c + 1], closeTo(C.get(r, c)[1], TOL));
      }
    }

  }

  void testZSum() {
    Float64List actual = A.sum();
    Float64List expected = new Float64List(2);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expected = DComplex.plus(expected, A.get(r, c));
      }
    }
    assertEquals(expected, actual, TOL);
  }

  void assertEquals(Float64List expected, Float64List actual, double tol) {
    for (int i = 0; i < actual.length; i++) {
      expect(expected[i], closeTo(actual[i], tol));
    }
  }

}
