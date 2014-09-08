part of cern.colt.matrix.double.test;

void testDoubleMatrix(String name, DoubleMatrixTest t) {
  group(name, () {
    setUp(t.setUp);
    tearDown(t.tearDown);
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

abstract class DoubleMatrixTest {

  /** Matrix to test. */
  AbstractDoubleMatrix A;

  /** Matrix of the same size as [A]. */
  AbstractDoubleMatrix B;

  /** Matrix of the size A.columns() x A.rows(). */
  AbstractDoubleMatrix Bt;

  final int NROWS = 13;
  final int NCOLUMNS = 17;

  final double TOL = 1e-10;

  void populateMatrices() {
    //ConcurrencyUtils.setThreadsBeginN_2D(1);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        A.set(r, c, random.nextDouble());
      }
    }

    for (int r = 0; r < B.rows; r++) {
      for (int c = 0; c < B.columns; c++) {
        B.set(r, c, random.nextDouble());
      }
    }

    for (int r = 0; r < Bt.rows; r++) {
      for (int c = 0; c < Bt.columns; c++) {
        Bt.set(r, c, random.nextDouble());
      }
    }
  }

  void createMatrices();

  void setUp() {
    createMatrices();
    populateMatrices();
  }

  void tearDown() {
    A = B = Bt = null;
  }

  testAggregate() {
    double expected = 0.0;
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        double elem = A.get(r, c);
        expected += elem * elem;
      }
    }
    double result = A.reduce(plus, square);
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
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        double elem = A.get(r, c);
        if (elem.abs() > 0.2) {
          expected += elem * elem;
        }
      }
    }

    double result = A.reduceWhere(plus, square, procedure);
    expect(expected, closeTo(result, TOL));
  }

  testAggregateIndex() {
    List<int> rowList = new List<int>();
    List<int> columnList = new List<int>();
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        rowList.add(r);
        columnList.add(c);
      }
    }
    double expected = 0.0;
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        double elem = A.get(r, c);
        expected += elem * elem;
      }
    }
    double result = A.reduceRange(plus, square,
        new Int32List.fromList(rowList), 
        new Int32List.fromList(columnList));
    expect(expected, closeTo(result, TOL));
  }

  testAggregateFunc() {
    double expected = 0.0;
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        double elemA = A.get(r, c);
        double elemB = B.get(r, c);
        expected += elemA * elemB;
      }
    }
    double result = A.reduceMatrix(B, plus, mult);
    expect(expected, closeTo(result, TOL));
  }

  testAssignValue() {
    double value = random.nextDouble();
    A.fill(value);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(value, closeTo(A.get(r, c), TOL));
      }
    }
  }

  testAssignValues2D() {
    List<Float64List> expected = new List<Float64List>(A.rows);//[A.columns()];
    for (int r = 0; r < A.rows; r++) {
      expected[r] = new Float64List(A.columns);
      for (int c = 0; c < A.columns; c++) {
        expected[r][c] = random.nextDouble();
      }
    }
    A.setAll2D(expected);
    for (int r = 0; r < A.rows; r++) {
      expect(A.columns == expected[r].length, isTrue);
      for (int c = 0; c < A.columns; c++) {
        expect(expected[r][c], closeTo(A.get(r, c), TOL));
      }
    }
  }

  testAssign() {
    AbstractDoubleMatrix Acopy = A.copy();
    A.forEach(acos);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        double expected = math.acos(Acopy.get(r, c));
        expect(expected, closeTo(A.get(r, c), TOL));
      }
    }
  }

  testAssignMatrix() {
    A.copyFrom(B);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(B.get(r, c), closeTo(A.get(r, c), TOL));
      }
    }
  }

  testAssignFunc() {
    AbstractDoubleMatrix Acopy = A.copy();
    A.forEachMatrix(B, plus);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(Acopy.get(r, c) + B.get(r, c), closeTo(A.get(r, c), TOL));
      }
    }
  }

  testAssignFuncIndex() {
    List<int> rowList = new List<int>();
    List<int> columnList = new List<int>();
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        rowList.add(r);
        columnList.add(c);
      }
    }
    AbstractDoubleMatrix Acopy = A.copy();
    A.forEachMatrixRange(B, div, new Int32List.fromList(rowList), new Int32List.fromList(columnList));
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(Acopy.get(r, c) / B.get(r, c), closeTo(A.get(r, c), TOL));
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
    AbstractDoubleMatrix Acopy = A.copy();
    A.fillWhere(procedure, -1.0);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        if (Acopy.get(r, c).abs() > 0.1) {
          expect(-1.0, closeTo(A.get(r, c), TOL));
        } else {
          expect(Acopy.get(r, c), closeTo(A.get(r, c), TOL));
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
    AbstractDoubleMatrix Acopy = A.copy();
    A.forEachWhere(procedure, tan);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        if (Acopy.get(r, c).abs() > 0.1) {
          expect(math.tan(Acopy.get(r, c)), closeTo(A.get(r, c), TOL));
        } else {
          expect(Acopy.get(r, c), closeTo(A.get(r, c), TOL));
        }
      }
    }
  }

  testCardinality() {
    int card = A.cardinality;
    expect(A.rows * A.columns, equals(card));
  }

  testEqualsValue() {
    double value = 1.0;
    A.fill(value);
    expect(A == value, isTrue);
    expect(A == 2, isFalse);
  }

  testEquals() {
    expect(A == A, isTrue);
    expect(A == B, isFalse);
  }

  testForEachNonZero() {
    AbstractDoubleMatrix Acopy = A.copy();
    double function(int first, int second, double third) {
      return math.sqrt(third);
    }
    A.forEachNonZero(function);
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(math.sqrt(Acopy.get(r, c)), closeTo(A.get(r, c), TOL));
      }
    }
  }

  testMaxLocation() {
    A.fill(0.0);
    A.set(A.rows ~/ 3, A.columns ~/ 3, 0.7);
    A.set(A.rows ~/ 2, A.columns ~/ 2, 0.1);
    final maxAndLoc = A.max();
    expect(0.7, closeTo(maxAndLoc.value, TOL));
    expect(A.rows ~/ 3, equals(maxAndLoc.row));
    expect(A.columns ~/ 3, equals(maxAndLoc.column));
  }

  testMinLocation() {
    A.fill(0.0);
    A.set(A.rows ~/ 3, A.columns ~/ 3, -0.7);
    A.set(A.rows ~/ 2, A.columns ~/ 2, -0.1);
    final minAndLoc = A.min();
    expect(-0.7, closeTo(minAndLoc.value, TOL));
    expect(A.rows ~/ 3, equals(minAndLoc.row));
    expect(A.columns ~/ 3, equals(minAndLoc.column));
  }

  testGetNegativeValues() {
    A.fill(0.0);
    A.set(A.rows ~/ 3, A.columns ~/ 3, -0.7);
    A.set(A.rows ~/ 2, A.columns ~/ 2, -0.1);
    List<int> rowList = new List<int>();
    List<int> columnList = new List<int>();
    List<double> valueList = new List<double>();
    A.negativeValues(rowList, columnList, valueList);
    expect(2, equals(rowList.length));
    expect(2, equals(columnList.length));
    expect(2, equals(valueList.length));
    expect(rowList.contains(A.rows ~/ 3), isTrue);
    expect(rowList.contains(A.rows ~/ 2), isTrue);
    expect(columnList.contains(A.columns ~/ 3), isTrue);
    expect(columnList.contains(A.columns ~/ 2), isTrue);
    expect(valueList.contains(-0.7), isTrue);
    expect(valueList.contains(-0.1), isTrue);
  }

  testGetNonZeros() {
    A.fill(0.0);
    A.set(A.rows ~/ 3, A.columns ~/ 3, 0.7);
    A.set(A.rows ~/ 2, A.columns ~/ 2, 0.1);
    List<int> rowList = new List<int>();
    List<int> columnList = new List<int>();
    List<double> valueList = new List<double>();
    A.nonZeros(rowList, columnList, valueList);
    expect(2, equals(rowList.length));
    expect(2, equals(columnList.length));
    expect(2, equals(valueList.length));
    expect(rowList.contains(A.rows ~/ 3), isTrue);
    expect(rowList.contains(A.rows ~/ 2), isTrue);
    expect(columnList.contains(A.columns ~/ 3), isTrue);
    expect(columnList.contains(A.columns ~/ 2), isTrue);
    expect(valueList.contains(0.7), isTrue);
    expect(valueList.contains(0.1), isTrue);
  }

  testGetPositiveValues() {
    A.fill(0.0);
    A.set(A.rows ~/ 3, A.columns ~/ 3, 0.7);
    A.set(A.rows ~/ 2, A.columns ~/ 2, 0.1);
    List<int> rowList = new List<int>();
    List<int> columnList = new List<int>();
    List<double> valueList = new List<double>();
    A.positiveValues(rowList, columnList, valueList);
    expect(2, equals(rowList.length));
    expect(2, equals(columnList.length));
    expect(2, equals(valueList.length));
    expect(rowList.contains(A.rows ~/ 3), isTrue);
    expect(rowList.contains(A.rows ~/ 2), isTrue);
    expect(columnList.contains(A.columns ~/ 3), isTrue);
    expect(columnList.contains(A.columns ~/ 2), isTrue);
    expect(valueList.contains(0.7), isTrue);
    expect(valueList.contains(0.1), isTrue);
  }

  testToArray() {
    List<Float64List> array = A.toList();
    expect(A.rows == array.length, isTrue);
    for (int r = 0; r < A.rows; r++) {
      expect(A.columns == array[r].length, isTrue);
      for (int c = 0; c < A.columns; c++) {
        expect(0, closeTo((array[r][c] - A.get(r, c)).abs(), TOL));
      }
    }
  }

  testVectorize() {
    AbstractDoubleVector Avec = A.vectorize();
    int idx = 0;
    for (int c = 0; c < A.columns; c++) {
      for (int r = 0; r < A.rows; r++) {
        expect(A.get(r, c), closeTo(Avec.get(idx++), TOL));
      }
    }
  }

  testViewColumn() {
    AbstractDoubleVector col = A.column(A.columns ~/ 2);
    expect(A.rows, equals(col.length));
    for (int r = 0; r < A.rows; r++) {
      expect(A.get(r, A.columns ~/ 2), closeTo(col.get(r), TOL));
    }
  }

  testViewColumnFlip() {
    AbstractDoubleMatrix B = A.columnFlip();
    expect(A.length, equals(B.length));
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(A.get(r, A.columns - 1 - c), closeTo(B.get(r, c), TOL));
      }
    }
  }

  testViewDice() {
    AbstractDoubleMatrix B = A.dice();
    expect(A.rows, equals(B.columns));
    expect(A.columns, equals(B.rows));
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(A.get(r, c), closeTo(B.get(c, r), TOL));
      }
    }
  }

  testViewPart() {
    AbstractDoubleMatrix B = A.part(A.rows ~/ 2, A.columns ~/ 2, A.rows ~/ 3, A.columns ~/ 3);
    expect(A.rows ~/ 3, equals(B.rows));
    expect(A.columns ~/ 3, equals(B.columns));
    for (int r = 0; r < A.rows / 3; r++) {
      for (int c = 0; c < A.columns / 3; c++) {
        expect(A.get(A.rows ~/ 2 + r, A.columns ~/ 2 + c), closeTo(B.get(r, c), TOL));
      }
    }
  }

  testViewRow() {
    AbstractDoubleVector B = A.row(A.rows ~/ 2);
    expect(A.columns, equals(B.length));
    for (int r = 0; r < A.columns; r++) {
      expect(A.get(A.rows ~/ 2, r), closeTo(B.get(r), TOL));
    }
  }

  testViewRowFlip() {
    AbstractDoubleMatrix B = A.rowFlip();
    expect(A.length, equals(B.length));
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(A.get(A.rows - 1 - r, c), closeTo(B.get(r, c), TOL));
      }
    }
  }

  testViewSelectionProc() {
    final double value = 2.0;
    A.fill(0.0);
    A.set(A.rows ~/ 4, 0, value);
    A.set(A.rows ~/ 2, 0, value);
    AbstractDoubleMatrix B = A.where((DoubleVector element) {
      if ((element.get(0) - value).abs() < TOL) {
        return true;
      } else {
        return false;
      }
    });
    expect(2, equals(B.rows));
    expect(A.columns, B.columns);
    expect(A.get(A.rows ~/ 4, 0), closeTo(B.get(0, 0), TOL));
    expect(A.get(A.rows ~/ 2, 0), closeTo(B.get(1, 0), TOL));
  }

  testViewSelection() {
    Int32List rowIndexes = new Int32List.fromList([A.rows ~/ 6, A.rows ~/ 5, A.rows ~/ 4, A.rows ~/ 3, A.rows ~/ 2]);
    Int32List colIndexes = new Int32List.fromList([A.columns ~/ 6, A.columns ~/ 5, A.columns ~/ 4, A.columns ~/ 3, A.columns ~/ 2, A.columns - 1]);
    AbstractDoubleMatrix B = A.select(rowIndexes, colIndexes);
    expect(rowIndexes.length, equals(B.rows));
    expect(colIndexes.length, equals(B.columns));
    for (int r = 0; r < rowIndexes.length; r++) {
      for (int c = 0; c < colIndexes.length; c++) {
        expect(A.get(rowIndexes[r], colIndexes[c]), closeTo(B.get(r, c), TOL));
      }
    }
  }

  /*test('viewSorted', () {
    DoubleMatrix B = A.viewSorted(1);
    for (int r = 0; r < A.rows() - 1; r++) {
      expect(B.getQuick(r + 1, 1) >= B.getQuick(r, 1), isTrue);
    }
  }*/

  testViewStrides() {
    int rowStride = 3;
    int colStride = 5;
    AbstractDoubleMatrix B = A.strides(rowStride, colStride);
    for (int r = 0; r < B.rows; r++) {
      for (int c = 0; c < B.columns; c++) {
        expect(A.get(r * rowStride, c * colStride), closeTo(B.get(r, c), TOL));
      }
    }
  }

  testZMult() {
    AbstractDoubleVector y = new DoubleVector(A.columns);
    for (int i = 0; i < y.length; i++) {
      y.set(i, random.nextDouble());
    }
    double alpha = 3.0;
    double beta = 5.0;
    AbstractDoubleVector z = DoubleFactory1D.dense.random(A.rows);
    Float64List expected = z.toList();
    z = A.mult(y, z, alpha, beta, false);
    for (int r = 0; r < A.rows; r++) {
      double s = 0.0;
      for (int c = 0; c < A.columns; c++) {
        s += A.get(r, c) * y.get(c);
      }
      expected[r] = s * alpha + expected[r] * beta;
    }

    for (int r = 0; r < A.rows; r++) {
      expect(expected[r], closeTo(z.get(r), TOL));
    }
    //---
    z = null;
    z = A.mult(y, z, alpha, beta, false);
    expected = new Float64List(A.rows);
    for (int r = 0; r < A.rows; r++) {
      double s = 0.0;
      for (int c = 0; c < A.columns; c++) {
        s += A.get(r, c) * y.get(c);
      }
      expected[r] = s * alpha;
    }
    for (int r = 0; r < A.rows; r++) {
      expect(expected[r], closeTo(z.get(r), TOL));
    }

    //transpose
    y = new DoubleVector(A.rows);
    for (int i = 0; i < y.length; i++) {
      y.set(i, random.nextDouble());
    }
    z = DoubleFactory1D.dense.random(A.columns);
    expected = z.toList();
    z = A.mult(y, z, alpha, beta, true);
    for (int r = 0; r < A.columns; r++) {
      double s = 0.0;
      for (int c = 0; c < A.rows; c++) {
        s += A.get(c, r) * y.get(c);
      }
      expected[r] = s * alpha + expected[r] * beta;
    }
    for (int r = 0; r < A.columns; r++) {
      expect(expected[r], closeTo(z.get(r), TOL));
    }
    //---
    z = null;
    z = A.mult(y, z, alpha, beta, true);
    expected = new Float64List(A.columns);
    for (int r = 0; r < A.columns; r++) {
      double s = 0.0;
      for (int c = 0; c < A.rows; c++) {
        s += A.get(c, r) * y.get(c);
      }
      expected[r] = s * alpha;
    }
    for (int r = 0; r < A.columns; r++) {
      expect(expected[r], closeTo(z.get(r), TOL));
    }
  }

  testZMult2D() {
    double alpha = 3.0;
    double beta = 5.0;
    AbstractDoubleMatrix C = DoubleFactory2D.dense.random(A.rows, A.rows);
    List<Float64List> expected = C.toList();
    C = A.multiply(Bt, C, alpha, beta, false, false);
    for (int j = 0; j < A.rows; j++) {
      for (int i = 0; i < A.rows; i++) {
        double s = 0.0;
        for (int k = 0; k < A.columns; k++) {
          s += A.get(i, k) * Bt.get(k, j);
        }
        expected[i][j] = s * alpha + expected[i][j] * beta;
      }
    }
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.rows; c++) {
        expect(expected[r][c], closeTo(C.get(r, c), TOL));
      }
    }

    //---
    C = null;
    C = A.multiply(Bt, C, alpha, beta, false, false);
    expected = new List<Float64List>.generate(A.rows, (_) => new Float64List(A.rows));
    for (int j = 0; j < A.rows; j++) {
      for (int i = 0; i < A.rows; i++) {
        double s = 0.0;
        for (int k = 0; k < A.columns; k++) {
          s += A.get(i, k) * Bt.get(k, j);
        }
        expected[i][j] = s * alpha;
      }
    }
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.rows; c++) {
        expect(expected[r][c], closeTo(C.get(r, c), TOL));
      }
    }

    //transposeA
    C = DoubleFactory2D.dense.random(A.columns, A.columns);
    expected = C.toList();
    C = A.multiply(B, C, alpha, beta, true, false);
    for (int j = 0; j < A.columns; j++) {
      for (int i = 0; i < A.columns; i++) {
        double s = 0.0;
        for (int k = 0; k < A.rows; k++) {
          s += A.get(k, i) * B.get(k, j);
        }
        expected[i][j] = s * alpha + expected[i][j] * beta;
      }
    }
    for (int r = 0; r < A.columns; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(expected[r][c], closeTo(C.get(r, c), TOL));
      }
    }
    //---
    C = null;
    C = A.multiply(B, C, alpha, beta, true, false);
    expected = new List<Float64List>.generate(A.columns,
        (_) => new Float64List(A.columns));
    for (int j = 0; j < A.columns; j++) {
      for (int i = 0; i < A.columns; i++) {
        double s = 0.0;
        for (int k = 0; k < A.rows; k++) {
          s += A.get(k, i) * B.get(k, j);
        }
        expected[i][j] = s * alpha;
      }
    }
    for (int r = 0; r < A.columns; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(expected[r][c], closeTo(C.get(r, c), TOL));
      }
    }

    //transposeB
    C = DoubleFactory2D.dense.random(A.rows, A.rows);
    expected = C.toList();
    C = A.multiply(B, C, alpha, beta, false, true);
    for (int j = 0; j < A.rows; j++) {
      for (int i = 0; i < A.rows; i++) {
        double s = 0.0;
        for (int k = 0; k < A.columns; k++) {
          s += A.get(i, k) * B.get(j, k);
        }
        expected[i][j] = s * alpha + expected[i][j] * beta;
      }
    }
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.rows; c++) {
        expect(expected[r][c], closeTo(C.get(r, c), TOL));
      }
    }
    //---
    C = null;
    C = A.multiply(B, C, alpha, beta, false, true);
    expected = new List<Float64List>.generate(A.rows,
        (_) => new Float64List(A.rows));
    for (int j = 0; j < A.rows; j++) {
      for (int i = 0; i < A.rows; i++) {
        double s = 0.0;
        for (int k = 0; k < A.columns; k++) {
          s += A.get(i, k) * B.get(j, k);
        }
        expected[i][j] = s * alpha;
      }
    }
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.rows; c++) {
        expect(expected[r][c], closeTo(C.get(r, c), TOL));
      }
    }
    //transposeA and transposeB
    C = DoubleFactory2D.dense.random(A.columns, A.columns);
    expected = C.toList();
    C = A.multiply(Bt, C, alpha, beta, true, true);
    for (int j = 0; j < A.columns; j++) {
      for (int i = 0; i < A.columns; i++) {
        double s = 0.0;
        for (int k = 0; k < A.rows; k++) {
          s += A.get(k, i) * Bt.get(j, k);
        }
        expected[i][j] = s * alpha + expected[i][j] * beta;
      }
    }
    for (int r = 0; r < A.columns; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(expected[r][c], closeTo(C.get(r, c), TOL));
      }
    }
    //---
    C = null;
    C = A.multiply(Bt, C, alpha, beta, true, true);
    expected = new List<Float64List>.generate(A.columns,
        (_) => new Float64List(A.columns));
    for (int j = 0; j < A.columns; j++) {
      for (int i = 0; i < A.columns; i++) {
        double s = 0.0;
        for (int k = 0; k < A.rows; k++) {
          s += A.get(k, i) * Bt.get(j, k);
        }
        expected[i][j] = s * alpha;
      }
    }
    for (int r = 0; r < A.columns; r++) {
      for (int c = 0; c < A.columns; c++) {
        expect(expected[r][c], closeTo(C.get(r, c), TOL));
      }
    }

  }

  testZSum() {
    double sum = A.sum();
    double expected = 0.0;
    for (int r = 0; r < A.rows; r++) {
      for (int c = 0; c < A.columns; c++) {
        expected += A.get(r, c);
      }
    }
    expect(expected, closeTo(sum, TOL));
  }
}
