part of cern.colt.matrix.double.test;

void testDoubleMatrix2D(String name, DoubleMatrix2DTest t) {
  group('DoubleMatrix2D ($name)', () {
    setUp() => t.setUp;
    tearDown() => t.tearDown();
    test('aggregate', t.testAggregate);
    test('aggregateProc', t.testAggregateProc);
    test('aggregateIndex', t.testAggregateIndex);
    test('aggregateFunc', t.testAggregateFunc);
    test('assignValue', t.testAssignValue);
    test('assignValues2D', t.testAssignValues2D);
    test('assign', t.testAssign);
    test('assignMatrix', t.testAssignMatrix);
    test('assignFunc', t.testAssignFunc);
    test('assignFuncIndex', t.testAssignFuncIndex);
    test('assignProc', t.testAssignProc);
    test('assignProcFunc', t.testAssignProcFunc);
    test('cardinality', t.testCardinality);
    test('equalsValue', t.testEqualsValue);
    test('equals', t.testEquals);
    test('forEachNonZero', t.testForEachNonZero);
    test('maxLocation', t.testMaxLocation);
    test('minLocation', t.testMinLocation);
    test('getNegativeValues', t.testGetNegativeValues);
    test('getNonZeros', t.testGetNonZeros);
    test('getPositiveValues', t.testGetPositiveValues);
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
    test('testZSum', t.testZSum);
  });
}

abstract class DoubleMatrix2DTest {

  /** Matrix to test. */
  DoubleMatrix2D A;

  /** Matrix of the same size as [A]. */
  DoubleMatrix2D B;

  /** Matrix of the size A.columns() x A.rows(). */
  DoubleMatrix2D Bt;

  final int NROWS = 13;
  final int NCOLUMNS = 17;

  final double TOL = 1e-10;

  void populateMatrices() {
    //ConcurrencyUtils.setThreadsBeginN_2D(1);
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        A.setQuick(r, c, random.nextDouble());
      }
    }

    for (int r = 0; r < B.rows(); r++) {
      for (int c = 0; c < B.columns(); c++) {
        B.setQuick(r, c, random.nextDouble());
      }
    }

    for (int r = 0; r < Bt.rows(); r++) {
      for (int c = 0; c < Bt.columns(); c++) {
        Bt.setQuick(r, c, random.nextDouble());
      }
    }
  }

  void createMatrices();

  void runTests() {
    group('DoubleMatrix2D', () {
      setUp() => this.setUp;
      tearDown() => this.tearDown();
      test('aggregate', testAggregate);
      test('aggregateProc', testAggregateProc);
      test('aggregateIndex', testAggregateIndex);
      test('aggregateFunc', testAggregateFunc);
      test('assignValue', testAssignValue);
      test('assignValues2D', testAssignValues2D);
      test('assign', testAssign);
      test('assignMatrix', testAssignMatrix);
      test('assignFunc', testAssignFunc);
      test('assignFuncIndex', testAssignFuncIndex);
      test('assignProc', testAssignProc);
      test('assignProcFunc', testAssignProcFunc);
      test('cardinality', testCardinality);
      test('equalsValue', testEqualsValue);
      test('equals', testEquals);
      test('forEachNonZero', testForEachNonZero);
      test('maxLocation', testMaxLocation);
      test('minLocation', testMinLocation);
      test('getNegativeValues', testGetNegativeValues);
      test('getNonZeros', testGetNonZeros);
      test('getPositiveValues', testGetPositiveValues);
      test('toArray', testToArray);
      test('vectorize', testVectorize);
      test('viewColumn', testViewColumn);
      test('viewColumnFlip', testViewColumnFlip);
      test('viewDice', testViewDice);
      test('viewPart', testViewPart);
      test('viewRow', testViewRow);
      test('viewRowFlip', testViewRowFlip);
      test('viewSelectionProc', testViewSelectionProc);
      test('viewSelection', testViewSelection);
      test('viewStrides', testViewStrides);
      test('zMult', testZMult);
      test('zMult2D', testZMult2D);
      test('testZSum', testZSum);
    });
  }

  void setUp() {
    createMatrices();
    populateMatrices();
  }

  void tearDown() {
    A = B = Bt = null;
  }

  testAggregate() {
    double expected = 0.0;
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        double elem = A.getQuick(r, c);
        expected += elem * elem;
      }
    }
    double result = A.aggregate(plus, square);
    expect(expected, closeTo(result, TOL));
  }

  testAggregateProc() {
    bool procedure(double element) {
      if (element.abs() > 0.2) {
        return true;
      } else {
        return false;
      }
    }
    double expected = 0.0;
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        double elem = A.getQuick(r, c);
        if (elem.abs() > 0.2) {
          expected += elem * elem;
        }
      }
    }

    double result = A.aggregateProc(plus, square, procedure);
    expect(expected, closeTo(result, TOL));
  }

  testAggregateIndex() {
    List<int> rowList = new List<int>();
    List<int> columnList = new List<int>();
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        rowList.add(r);
        columnList.add(c);
      }
    }
    double expected = 0.0;
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        double elem = A.getQuick(r, c);
        expected += elem * elem;
      }
    }
    double result = A.aggregateIndex(plus, square, rowList, columnList);
    expect(expected, closeTo(result, TOL));
  }

  testAggregateFunc() {
    double expected = 0.0;
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        double elemA = A.getQuick(r, c);
        double elemB = B.getQuick(r, c);
        expected += elemA * elemB;
      }
    }
    double result = A.aggregateFunc(B, plus, mult);
    expect(expected, closeTo(result, TOL));
  }

  testAssignValue() {
    double value = random.nextDouble();
    A.assignValue(value);
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        expect(value, closeTo(A.getQuick(r, c), TOL));
      }
    }
  }

  testAssignValues2D() {
    List<Float64List> expected = new List<Float64List>(A.rows());//[A.columns()];
    for (int r = 0; r < A.rows(); r++) {
      expected[r] = new Float64List(A.columns());
      for (int c = 0; c < A.columns(); c++) {
        expected[r][c] = random.nextDouble();
      }
    }
    A.assignValues2D(expected);
    for (int r = 0; r < A.rows(); r++) {
      expect(A.columns() == expected[r].length, isTrue);
      for (int c = 0; c < A.columns(); c++) {
        expect(expected[r][c], closeTo(A.getQuick(r, c), TOL));
      }
    }
  }

  testAssign() {
    DoubleMatrix2D Acopy = A.copy();
    A.assign(acos);
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        double expected = math.acos(Acopy.getQuick(r, c));
        expect(expected, closeTo(A.getQuick(r, c), TOL));
      }
    }
  }

  testAssignMatrix() {
    A.assignMatrix(B);
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        expect(B.getQuick(r, c), closeTo(A.getQuick(r, c), TOL));
      }
    }
  }

  testAssignFunc() {
    DoubleMatrix2D Acopy = A.copy();
    A.assignFunc(B, plus);
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        expect(Acopy.getQuick(r, c) + B.getQuick(r, c), closeTo(A.getQuick(r, c), TOL));
      }
    }
  }

  testAssignFuncIndex() {
    List<int> rowList = new List<int>();
    List<int> columnList = new List<int>();
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        rowList.add(r);
        columnList.add(c);
      }
    }
    DoubleMatrix2D Acopy = A.copy();
    A.assignFuncIndex(B, div, rowList, columnList);
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        expect(Acopy.getQuick(r, c) / B.getQuick(r, c), closeTo(A.getQuick(r, c), TOL));
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
    DoubleMatrix2D Acopy = A.copy();
    A.assignProc(procedure, -1.0);
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        if (Acopy.getQuick(r, c).abs() > 0.1) {
          expect(-1.0, closeTo(A.getQuick(r, c), TOL));
        } else {
          expect(Acopy.getQuick(r, c), closeTo(A.getQuick(r, c), TOL));
        }
      }
    }
  }

  testAssignProcFunc() {
    bool procedure(double element) {
      if (element.abs() > 0.1) {
        return true;
      } else {
        return false;
      }
    }
    DoubleMatrix2D Acopy = A.copy();
    A.assignProcFunc(procedure, tan);
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        if (Acopy.getQuick(r, c).abs() > 0.1) {
          expect(math.tan(Acopy.getQuick(r, c)), closeTo(A.getQuick(r, c), TOL));
        } else {
          expect(Acopy.getQuick(r, c), closeTo(A.getQuick(r, c), TOL));
        }
      }
    }
  }

  testCardinality() {
    int card = A.cardinality();
    expect(A.rows() * A.columns(), equals(card));
  }

  testEqualsValue() {
    double value = 1.0;
    A.assignValue(value);
    bool eq = A.equals(value);
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

  testForEachNonZero() {
    DoubleMatrix2D Acopy = A.copy();
    double function(int first, int second, double third) {
      return math.sqrt(third);
    }
    ;
    A.forEachNonZero(function);
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        expect(math.sqrt(Acopy.getQuick(r, c)), closeTo(A.getQuick(r, c), TOL));
      }
    }
  }

  testMaxLocation() {
    A.assignValue(0.0);
    A.setQuick(A.rows() ~/ 3, A.columns() ~/ 3, 0.7);
    A.setQuick(A.rows() ~/ 2, A.columns() ~/ 2, 0.1);
    Float64List maxAndLoc = A.getMaxLocation();
    expect(0.7, closeTo(maxAndLoc[0], TOL));
    expect(A.rows() / 3, equals(maxAndLoc[1]));
    expect(A.columns() / 3, equals(maxAndLoc[2]));
  }

  testMinLocation() {
    A.assignValue(0.0);
    A.setQuick(A.rows() ~/ 3, A.columns() ~/ 3, -0.7);
    A.setQuick(A.rows() ~/ 2, A.columns() ~/ 2, -0.1);
    Float64List minAndLoc = A.getMinLocation();
    expect(-0.7, closeTo(minAndLoc[0], TOL));
    expect(A.rows() ~/ 3, equals(minAndLoc[1]));
    expect(A.columns() ~/ 3, equals(minAndLoc[2]));
  }

  testGetNegativeValues() {
    A.assignValue(0.0);
    A.setQuick(A.rows() ~/ 3, A.columns() ~/ 3, -0.7);
    A.setQuick(A.rows() ~/ 2, A.columns() ~/ 2, -0.1);
    List<int> rowList = new List<int>();
    List<int> columnList = new List<int>();
    List<double> valueList = new List<double>();
    A.getNegativeValues(rowList, columnList, valueList);
    expect(2, equals(rowList.length));
    expect(2, equals(columnList.length));
    expect(2, equals(valueList.length));
    expect(rowList.contains(A.rows() / 3), isTrue);
    expect(rowList.contains(A.rows() / 2), isTrue);
    expect(columnList.contains(A.columns() / 3), isTrue);
    expect(columnList.contains(A.columns() / 2), isTrue);
    expect(valueList.contains(-0.7), isTrue);
    expect(valueList.contains(-0.1), isTrue);
  }

  testGetNonZeros() {
    A.assignValue(0.0);
    A.setQuick(A.rows() ~/ 3, A.columns() ~/ 3, 0.7);
    A.setQuick(A.rows() ~/ 2, A.columns() ~/ 2, 0.1);
    List<int> rowList = new List<int>();
    List<int> columnList = new List<int>();
    List<double> valueList = new List<double>();
    A.getNonZeros(rowList, columnList, valueList);
    expect(2, equals(rowList.length));
    expect(2, equals(columnList.length));
    expect(2, equals(valueList.length));
    expect(rowList.contains(A.rows() / 3), isTrue);
    expect(rowList.contains(A.rows() / 2), isTrue);
    expect(columnList.contains(A.columns() / 3), isTrue);
    expect(columnList.contains(A.columns() / 2), isTrue);
    expect(valueList.contains(0.7), isTrue);
    expect(valueList.contains(0.1), isTrue);
  }

  testGetPositiveValues() {
    A.assignValue(0.0);
    A.setQuick(A.rows() ~/ 3, A.columns() ~/ 3, 0.7);
    A.setQuick(A.rows() ~/ 2, A.columns() ~/ 2, 0.1);
    List<int> rowList = new List<int>();
    List<int> columnList = new List<int>();
    List<double> valueList = new List<double>();
    A.getPositiveValues(rowList, columnList, valueList);
    expect(2, equals(rowList.length));
    expect(2, equals(columnList.length));
    expect(2, equals(valueList.length));
    expect(rowList.contains(A.rows() / 3), isTrue);
    expect(rowList.contains(A.rows() / 2), isTrue);
    expect(columnList.contains(A.columns() / 3), isTrue);
    expect(columnList.contains(A.columns() / 2), isTrue);
    expect(valueList.contains(0.7), isTrue);
    expect(valueList.contains(0.1), isTrue);
  }

  testToArray() {
    List<Float64List> array = A.toArray();
    expect(A.rows() == array.length, isTrue);
    for (int r = 0; r < A.rows(); r++) {
      expect(A.columns() == array[r].length, isTrue);
      for (int c = 0; c < A.columns(); c++) {
        expect(0, closeTo((array[r][c] - A.getQuick(r, c)).abs(), TOL));
      }
    }
  }

  testVectorize() {
    DoubleMatrix1D Avec = A.vectorize();
    int idx = 0;
    for (int c = 0; c < A.columns(); c++) {
      for (int r = 0; r < A.rows(); r++) {
        expect(A.getQuick(r, c), closeTo(Avec.getQuick(idx++), TOL));
      }
    }
  }

  testViewColumn() {
    DoubleMatrix1D col = A.viewColumn(A.columns() ~/ 2);
    expect(A.rows(), equals(col.size()));
    for (int r = 0; r < A.rows(); r++) {
      expect(A.getQuick(r, A.columns() ~/ 2), closeTo(col.getQuick(r), TOL));
    }
  }

  testViewColumnFlip() {
    DoubleMatrix2D B = A.viewColumnFlip();
    expect(A.size(), equals(B.size()));
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        expect(A.getQuick(r, A.columns() - 1 - c), closeTo(B.getQuick(r, c), TOL));
      }
    }
  }

  testViewDice() {
    DoubleMatrix2D B = A.viewDice();
    expect(A.rows(), equals(B.columns()));
    expect(A.columns(), equals(B.rows()));
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        expect(A.getQuick(r, c), closeTo(B.getQuick(c, r), TOL));
      }
    }
  }

  testViewPart() {
    DoubleMatrix2D B = A.viewPart(A.rows() ~/ 2, A.columns() ~/ 2, A.rows() ~/ 3, A.columns() ~/ 3);
    expect(A.rows() / 3, equals(B.rows()));
    expect(A.columns() / 3, equals(B.columns()));
    for (int r = 0; r < A.rows() / 3; r++) {
      for (int c = 0; c < A.columns() / 3; c++) {
        expect(A.getQuick(A.rows() ~/ 2 + r, A.columns() ~/ 2 + c), closeTo(B.getQuick(r, c), TOL));
      }
    }
  }

  testViewRow() {
    DoubleMatrix1D B = A.viewRow(A.rows() ~/ 2);
    expect(A.columns(), equals(B.size()));
    for (int r = 0; r < A.columns(); r++) {
      expect(A.getQuick(A.rows() ~/ 2, r), closeTo(B.getQuick(r), TOL));
    }
  }

  testViewRowFlip() {
    DoubleMatrix2D B = A.viewRowFlip();
    expect(A.size(), equals(B.size()));
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        expect(A.getQuick(A.rows() - 1 - r, c), closeTo(B.getQuick(r, c), TOL));
      }
    }
  }

  testViewSelectionProc() {
    final double value = 2.0;
    A.assignValue(0.0);
    A.setQuick(A.rows() ~/ 4, 0, value);
    A.setQuick(A.rows() ~/ 2, 0, value);
    DoubleMatrix2D B = A.viewSelectionProc((DoubleMatrix1D element) {
      if ((element.getQuick(0) - value).abs() < TOL) {
        return true;
      } else {
        return false;
      }
    });
    expect(2, equals(B.rows()));
    expect(A.columns(), B.columns());
    expect(A.getQuick(A.rows() ~/ 4, 0), closeTo(B.getQuick(0, 0), TOL));
    expect(A.getQuick(A.rows() ~/ 2, 0), closeTo(B.getQuick(1, 0), TOL));
  }

  testViewSelection() {
    Int32List rowIndexes = new Int32List.fromList([A.rows() / 6, A.rows() / 5, A.rows() / 4, A.rows() / 3, A.rows() / 2]);
    Int32List colIndexes = new Int32List.fromList([A.columns() / 6, A.columns() / 5, A.columns() / 4, A.columns() / 3, A.columns() / 2, A.columns() - 1]);
    DoubleMatrix2D B = A.viewSelection(rowIndexes, colIndexes);
    expect(rowIndexes.length, equals(B.rows()));
    expect(colIndexes.length, equals(B.columns()));
    for (int r = 0; r < rowIndexes.length; r++) {
      for (int c = 0; c < colIndexes.length; c++) {
        expect(A.getQuick(rowIndexes[r], colIndexes[c]), closeTo(B.getQuick(r, c), TOL));
      }
    }
  }

  /*test('viewSorted', () {
    DoubleMatrix2D B = A.viewSorted(1);
    for (int r = 0; r < A.rows() - 1; r++) {
      expect(B.getQuick(r + 1, 1) >= B.getQuick(r, 1), isTrue);
    }
  }*/

  testViewStrides() {
    int rowStride = 3;
    int colStride = 5;
    DoubleMatrix2D B = A.viewStrides(rowStride, colStride);
    for (int r = 0; r < B.rows(); r++) {
      for (int c = 0; c < B.columns(); c++) {
        expect(A.getQuick(r * rowStride, c * colStride), closeTo(B.getQuick(r, c), TOL));
      }
    }
  }

  testZMult() {
    DoubleMatrix1D y = new DenseDoubleMatrix1D(A.columns());
    for (int i = 0; i < y.size(); i++) {
      y.setQuick(i, random.nextDouble());
    }
    double alpha = 3.0;
    double beta = 5.0;
    DoubleMatrix1D z = DoubleFactory1D.dense.random(A.rows());
    Float64List expected = z.toArray();
    z = A.zMult(y, z, alpha, beta, false);
    for (int r = 0; r < A.rows(); r++) {
      double s = 0.0;
      for (int c = 0; c < A.columns(); c++) {
        s += A.getQuick(r, c) * y.getQuick(c);
      }
      expected[r] = s * alpha + expected[r] * beta;
    }

    for (int r = 0; r < A.rows(); r++) {
      expect(expected[r], closeTo(z.getQuick(r), TOL));
    }
    //---
    z = null;
    z = A.zMult(y, z, alpha, beta, false);
    expected = new Float64List(A.rows());
    for (int r = 0; r < A.rows(); r++) {
      double s = 0.0;
      for (int c = 0; c < A.columns(); c++) {
        s += A.getQuick(r, c) * y.getQuick(c);
      }
      expected[r] = s * alpha;
    }
    for (int r = 0; r < A.rows(); r++) {
      expect(expected[r], closeTo(z.getQuick(r), TOL));
    }

    //transpose
    y = new DenseDoubleMatrix1D(A.rows());
    for (int i = 0; i < y.size(); i++) {
      y.setQuick(i, random.nextDouble());
    }
    z = DoubleFactory1D.dense.random(A.columns());
    expected = z.toArray();
    z = A.zMult(y, z, alpha, beta, true);
    for (int r = 0; r < A.columns(); r++) {
      double s = 0.0;
      for (int c = 0; c < A.rows(); c++) {
        s += A.getQuick(c, r) * y.getQuick(c);
      }
      expected[r] = s * alpha + expected[r] * beta;
    }
    for (int r = 0; r < A.columns(); r++) {
      expect(expected[r], closeTo(z.getQuick(r), TOL));
    }
    //---
    z = null;
    z = A.zMult(y, z, alpha, beta, true);
    expected = new Float64List(A.columns());
    for (int r = 0; r < A.columns(); r++) {
      double s = 0.0;
      for (int c = 0; c < A.rows(); c++) {
        s += A.getQuick(c, r) * y.getQuick(c);
      }
      expected[r] = s * alpha;
    }
    for (int r = 0; r < A.columns(); r++) {
      expect(expected[r], closeTo(z.getQuick(r), TOL));
    }
  }

  testZMult2D() {
    double alpha = 3.0;
    double beta = 5.0;
    DoubleMatrix2D C = DoubleFactory2D.dense.random(A.rows(), A.rows());
    List<Float64List> expected = C.toArray();
    C = A.zMult2D(Bt, C, alpha, beta, false, false);
    for (int j = 0; j < A.rows(); j++) {
      for (int i = 0; i < A.rows(); i++) {
        double s = 0.0;
        for (int k = 0; k < A.columns(); k++) {
          s += A.getQuick(i, k) * Bt.getQuick(k, j);
        }
        expected[i][j] = s * alpha + expected[i][j] * beta;
      }
    }
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.rows(); c++) {
        expect(expected[r][c], closeTo(C.getQuick(r, c), TOL));
      }
    }

    //---
    C = null;
    C = A.zMult2D(Bt, C, alpha, beta, false, false);
    expected = new List<Float64List>(A.rows());//[A.rows()];
    for (int j = 0; j < A.rows(); j++) {
      expected[j] = new Float64List(A.rows());
      for (int i = 0; i < A.rows(); i++) {
        double s = 0.0;
        for (int k = 0; k < A.columns(); k++) {
          s += A.getQuick(i, k) * Bt.getQuick(k, j);
        }
        expected[i][j] = s * alpha;
      }
    }
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.rows(); c++) {
        expect(expected[r][c], closeTo(C.getQuick(r, c), TOL));
      }
    }

    //transposeA
    C = DoubleFactory2D.dense.random(A.columns(), A.columns());
    expected = C.toArray();
    C = A.zMult2D(B, C, alpha, beta, true, false);
    for (int j = 0; j < A.columns(); j++) {
      for (int i = 0; i < A.columns(); i++) {
        double s = 0.0;
        for (int k = 0; k < A.rows(); k++) {
          s += A.getQuick(k, i) * B.getQuick(k, j);
        }
        expected[i][j] = s * alpha + expected[i][j] * beta;
      }
    }
    for (int r = 0; r < A.columns(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        expect(expected[r][c], closeTo(C.getQuick(r, c), TOL));
      }
    }
    //---
    C = null;
    C = A.zMult2D(B, C, alpha, beta, true, false);
    expected = new List<Float64List>(A.columns());//[A.columns()];
    for (int j = 0; j < A.columns(); j++) {
      expected[j] = new Float64List(A.columns());
      for (int i = 0; i < A.columns(); i++) {
        double s = 0.0;
        for (int k = 0; k < A.rows(); k++) {
          s += A.getQuick(k, i) * B.getQuick(k, j);
        }
        expected[i][j] = s * alpha;
      }
    }
    for (int r = 0; r < A.columns(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        expect(expected[r][c], closeTo(C.getQuick(r, c), TOL));
      }
    }

    //transposeB
    C = DoubleFactory2D.dense.random(A.rows(), A.rows());
    expected = C.toArray();
    C = A.zMult2D(B, C, alpha, beta, false, true);
    for (int j = 0; j < A.rows(); j++) {
      for (int i = 0; i < A.rows(); i++) {
        double s = 0.0;
        for (int k = 0; k < A.columns(); k++) {
          s += A.getQuick(i, k) * B.getQuick(j, k);
        }
        expected[i][j] = s * alpha + expected[i][j] * beta;
      }
    }
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.rows(); c++) {
        expect(expected[r][c], closeTo(C.getQuick(r, c), TOL));
      }
    }
    //---
    C = null;
    C = A.zMult2D(B, C, alpha, beta, false, true);
    expected = new List<Float64List>(A.rows());//[A.rows()];
    for (int j = 0; j < A.rows(); j++) {
      expected[j] = new Float64List(A.rows());
      for (int i = 0; i < A.rows(); i++) {
        double s = 0.0;
        for (int k = 0; k < A.columns(); k++) {
          s += A.getQuick(i, k) * B.getQuick(j, k);
        }
        expected[i][j] = s * alpha;
      }
    }
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.rows(); c++) {
        expect(expected[r][c], closeTo(C.getQuick(r, c), TOL));
      }
    }
    //transposeA and transposeB
    C = DoubleFactory2D.dense.random(A.columns(), A.columns());
    expected = C.toArray();
    C = A.zMult2D(Bt, C, alpha, beta, true, true);
    for (int j = 0; j < A.columns(); j++) {
      for (int i = 0; i < A.columns(); i++) {
        double s = 0.0;
        for (int k = 0; k < A.rows(); k++) {
          s += A.getQuick(k, i) * Bt.getQuick(j, k);
        }
        expected[i][j] = s * alpha + expected[i][j] * beta;
      }
    }
    for (int r = 0; r < A.columns(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        expect(expected[r][c], closeTo(C.getQuick(r, c), TOL));
      }
    }
    //---
    C = null;
    C = A.zMult2D(Bt, C, alpha, beta, true, true);
    expected = new List<Float64List>(A.columns());//[A.columns()];
    for (int j = 0; j < A.columns(); j++) {
      expected[j] = new Float64List(A.columns());
      for (int i = 0; i < A.columns(); i++) {
        double s = 0.0;
        for (int k = 0; k < A.rows(); k++) {
          s += A.getQuick(k, i) * Bt.getQuick(j, k);
        }
        expected[i][j] = s * alpha;
      }
    }
    for (int r = 0; r < A.columns(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        expect(expected[r][c], closeTo(C.getQuick(r, c), TOL));
      }
    }

  }

  testZSum() {
    double sum = A.zSum();
    double expected = 0.0;
    for (int r = 0; r < A.rows(); r++) {
      for (int c = 0; c < A.columns(); c++) {
        expected += A.getQuick(r, c);
      }
    }
    expect(expected, closeTo(sum, TOL));
  }
}
