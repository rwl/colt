part of cern.colt.matrix.double.test;

class DiagonalDoubleMatrixTest extends DoubleMatrixTest {
  int DLENGTH;

  int DINDEX;

  void createMatrices() {
    DINDEX = 3;
    A = new DiagonalDoubleMatrix(NROWS, NCOLUMNS, DINDEX);
    B = new DiagonalDoubleMatrix(NROWS, NCOLUMNS, DINDEX);
    Bt = new DiagonalDoubleMatrix(NCOLUMNS, NROWS, -DINDEX);
    DLENGTH = (A as DiagonalDoubleMatrix).diagonalLength;
  }

  void populateMatrices() {
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        A.set(r, r + DINDEX, random.nextDouble());
      }

      for (int r = 0; r < DLENGTH; r++) {
        B.set(r, r + DINDEX, random.nextDouble());
      }

      for (int r = 0; r < DLENGTH; r++) {
        Bt.set(r - DINDEX, r, random.nextDouble());
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        A.set(r - DINDEX, r, random.nextDouble());
      }

      for (int r = 0; r < DLENGTH; r++) {
        B.set(r - DINDEX, r, random.nextDouble());
      }
      for (int r = 0; r < DLENGTH; r++) {
        Bt.set(r, r + DINDEX, random.nextDouble());
      }
    }
  }

  void testFill() {
    double value = random.nextDouble();
    A.fill(value);
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        expect(value, closeTo(A.get(r, r + DINDEX), TOL));
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        expect(value, closeTo(A.get(r - DINDEX, r), TOL));
      }
    }
  }

  void testSetAll() {
    Float64List expected = new Float64List(DLENGTH);
    for (int i = 0; i < DLENGTH; i++) {
      expected[i] = random.nextDouble();
    }
    A.setAll(expected);
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        expect(expected[r], closeTo(A.get(r, r + DINDEX), TOL));
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        expect(expected[r], closeTo(A.get(r - DINDEX, r), TOL));
      }
    }
  }

  void testApply() {
    AbstractDoubleMatrix Acopy = A.copy();
    A.apply(acos);
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        double expected = math.acos(Acopy.get(r, r + DINDEX));
        expect(expected, closeTo(A.get(r, r + DINDEX), TOL));
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        double expected = math.acos(Acopy.get(r - DINDEX, r));
        expect(expected, closeTo(A.get(r - DINDEX, r), TOL));
      }
    }
  }

  void testAssign() {
    AbstractDoubleMatrix Acopy = A.copy();
    A.assign(B, div);
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        expect(Acopy.get(r, r + DINDEX) / B.get(r, r + DINDEX),
            closeTo(A.get(r, r + DINDEX), TOL));
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        expect(Acopy.get(r - DINDEX, r) / B.get(r - DINDEX, r),
            closeTo(A.get(r - DINDEX, r), TOL));
      }
    }
  }

  void testCardinality() {
    int card = A.cardinality;
    expect(DLENGTH, equals(card));
  }

  void testMax() {
    A.fill(0.0);
    if (DINDEX >= 0) {
      A.set(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, 0.7);
      A.set(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, 0.1);
      final maxAndLoc = A.max();
      expect(0.7, closeTo(maxAndLoc.value, TOL));
      expect(NROWS ~/ 3, equals(maxAndLoc.row));
      expect(NROWS ~/ 3 + DINDEX, equals(maxAndLoc.column));
    } else {
      A.set(NROWS ~/ 3 - DINDEX, NROWS ~/ 3, 0.7);
      A.set(NROWS ~/ 2 - DINDEX, NROWS ~/ 2, 0.1);
      final maxAndLoc = A.max();
      expect(0.7, closeTo(maxAndLoc.value, TOL));
      expect(NROWS ~/ 3 - DINDEX, equals(maxAndLoc.row));
      expect(NROWS ~/ 3, equals(maxAndLoc.column));
    }
  }

  void testMin() {
    A.fill(0.0);
    if (DINDEX >= 0) {
      A.set(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, -0.7);
      A.set(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, -0.1);
      final minAndLoc = A.min();
      expect(-0.7, closeTo(minAndLoc.value, TOL));
      expect(NROWS ~/ 3, equals(minAndLoc.row));
      expect(NROWS ~/ 3 + DINDEX, equals(minAndLoc.column));
    } else {
      A.set(NROWS ~/ 3 - DINDEX, NROWS ~/ 3, -0.7);
      A.set(NROWS ~/ 2 - DINDEX, NROWS ~/ 2, -0.1);
      final minAndLoc = A.min();
      expect(-0.7, closeTo(minAndLoc.value, TOL));
      expect(NROWS ~/ 3 - DINDEX, equals(minAndLoc.row));
      expect(NROWS ~/ 3, equals(minAndLoc.column));
    }
  }

  void testNegative() {
    A.fill(0.0);
    if (DINDEX >= 0) {
      A.set(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, -0.7);
      A.set(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, -0.1);
      List<int> rowList = new List<int>();
      List<int> columnList = new List<int>();
      List<double> valueList = new List<double>();
      A.negative(rowList, columnList, valueList);
      expect(2, equals(rowList.length));
      expect(2, equals(columnList.length));
      expect(2, equals(valueList.length));
      expect(rowList.contains(NROWS ~/ 3), isTrue);
      expect(rowList.contains(NROWS ~/ 2), isTrue);
      expect(columnList.contains(NROWS ~/ 3 + DINDEX), isTrue);
      expect(columnList.contains(NROWS ~/ 2 + DINDEX), isTrue);
      expect(valueList.contains(-0.7), isTrue);
      expect(valueList.contains(-0.1), isTrue);
    } else {
      A.set(NROWS ~/ 3 - DINDEX, NROWS ~/ 3, -0.7);
      A.set(NROWS ~/ 2 - DINDEX, NROWS ~/ 2, -0.1);
      List<int> rowList = new List<int>();
      List<int> columnList = new List<int>();
      List<double> valueList = new List<double>();
      A.negative(rowList, columnList, valueList);
      expect(2, equals(rowList.length));
      expect(2, equals(columnList.length));
      expect(2, equals(valueList.length));
      expect(rowList.contains(NROWS ~/ 3 - DINDEX), isTrue);
      expect(rowList.contains(NROWS ~/ 2 - DINDEX), isTrue);
      expect(columnList.contains(NROWS ~/ 3), isTrue);
      expect(columnList.contains(NROWS ~/ 2), isTrue);
      expect(valueList.contains(-0.7), isTrue);
      expect(valueList.contains(-0.1), isTrue);
    }
  }

  void testNonzero() {
    A.fill(0.0);
    if (DINDEX >= 0) {
      A.set(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, 0.7);
      A.set(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, 0.1);
      List<int> rowList = new List<int>();
      List<int> columnList = new List<int>();
      List<double> valueList = new List<double>();
      A.nonzero(rowList, columnList, valueList);
      expect(2, equals(rowList.length));
      expect(2, equals(columnList.length));
      expect(2, equals(valueList.length));
      expect(rowList.contains(NROWS ~/ 3), isTrue);
      expect(rowList.contains(NROWS ~/ 2), isTrue);
      expect(columnList.contains(NROWS ~/ 3 + DINDEX), isTrue);
      expect(columnList.contains(NROWS ~/ 2 + DINDEX), isTrue);
      expect(valueList.contains(0.7), isTrue);
      expect(valueList.contains(0.1), isTrue);
    } else {
      A.set(NROWS ~/ 3 - DINDEX, NROWS ~/ 3, 0.7);
      A.set(NROWS ~/ 2 - DINDEX, NROWS ~/ 2, 0.1);
      List<int> rowList = new List<int>();
      List<int> columnList = new List<int>();
      List<double> valueList = new List<double>();
      A.nonzero(rowList, columnList, valueList);
      expect(2, equals(rowList.length));
      expect(2, equals(columnList.length));
      expect(2, equals(valueList.length));
      expect(rowList.contains(NROWS ~/ 3 - DINDEX), isTrue);
      expect(rowList.contains(NROWS ~/ 2 - DINDEX), isTrue);
      expect(columnList.contains(NROWS ~/ 3), isTrue);
      expect(columnList.contains(NROWS ~/ 2), isTrue);
      expect(valueList.contains(0.7), isTrue);
      expect(valueList.contains(0.1), isTrue);
    }
  }

  void testPositive() {
    A.fill(0.0);
    if (DINDEX >= 0) {
      A.set(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, 0.7);
      A.set(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, 0.1);
      List<int> rowList = new List<int>();
      List<int> columnList = new List<int>();
      List<double> valueList = new List<double>();
      A.positive(rowList, columnList, valueList);
      expect(2, equals(rowList.length));
      expect(2, equals(columnList.length));
      expect(2, equals(valueList.length));
      expect(rowList.contains(NROWS ~/ 3), isTrue);
      expect(rowList.contains(NROWS ~/ 2), isTrue);
      expect(columnList.contains(NROWS ~/ 3 + DINDEX), isTrue);
      expect(columnList.contains(NROWS ~/ 2 + DINDEX), isTrue);
      expect(valueList.contains(0.7), isTrue);
      expect(valueList.contains(0.1), isTrue);
    } else {
      A.set(NROWS ~/ 3 - DINDEX, NROWS ~/ 3, 0.7);
      A.set(NROWS ~/ 2 - DINDEX, NROWS ~/ 2, 0.1);
      List<int> rowList = new List<int>();
      List<int> columnList = new List<int>();
      List<double> valueList = new List<double>();
      A.positive(rowList, columnList, valueList);
      expect(2, equals(rowList.length));
      expect(2, equals(columnList.length));
      expect(2, equals(valueList.length));
      expect(rowList.contains(NROWS ~/ 3 - DINDEX), isTrue);
      expect(rowList.contains(NROWS ~/ 2 - DINDEX), isTrue);
      expect(columnList.contains(NROWS ~/ 3), isTrue);
      expect(columnList.contains(NROWS ~/ 2), isTrue);
      expect(valueList.contains(0.7), isTrue);
      expect(valueList.contains(0.1), isTrue);
    }
  }

  void testColumn() {
    AbstractDoubleVector col = A.column(NCOLUMNS ~/ 2);
    expect(NROWS, equals(col.size));
    for (int r = 0; r < NROWS; r++) {
      expect(A.get(r, NCOLUMNS ~/ 2), closeTo(col.get(r), TOL));
    }
  }

  void testColumnFlip() {
    AbstractDoubleMatrix B = A.columnFlip();
    expect(A.size, equals(B.size));
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(A.get(r, NCOLUMNS - 1 - c), closeTo(B.get(r, c), TOL));
      }
    }
  }

  void testDice() {
    AbstractDoubleMatrix B = A.dice();
    expect(NROWS, equals(B.columns));
    expect(NCOLUMNS, equals(B.rows));
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(A.get(r, c), closeTo(B.get(c, r), TOL));
      }
    }
  }

  void testPart() {
    AbstractDoubleMatrix B =
        A.part(NROWS ~/ 2, NCOLUMNS ~/ 2, NROWS ~/ 3, NCOLUMNS ~/ 3);
    expect(NROWS ~/ 3, equals(B.rows));
    expect(NCOLUMNS ~/ 3, equals(B.columns));
    for (int r = 0; r < NROWS ~/ 3; r++) {
      for (int c = 0; c < NCOLUMNS ~/ 3; c++) {
        expect(A.get(NROWS ~/ 2 + r, NCOLUMNS ~/ 2 + c),
            closeTo(B.get(r, c), TOL));
      }
    }
  }

  void testRow() {
    AbstractDoubleVector B = A.row(NROWS ~/ 2);
    expect(NCOLUMNS, equals(B.size));
    for (int r = 0; r < NCOLUMNS; r++) {
      expect(A.get(NROWS ~/ 2, r), closeTo(B.get(r), TOL));
    }
  }

  void testRowFlip() {
    AbstractDoubleMatrix B = A.rowFlip();
    expect(A.size, equals(B.size));
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(A.get(NROWS - 1 - r, c), closeTo(B.get(r, c), TOL));
      }
    }
  }

  void testSelect() {
    final rowIndexes = new Int32List.fromList(
        [NROWS ~/ 6, NROWS ~/ 5, NROWS ~/ 4, NROWS ~/ 3, NROWS ~/ 2]);
    final colIndexes = new Int32List.fromList([
      NROWS ~/ 6,
      NROWS ~/ 5,
      NROWS ~/ 4,
      NROWS ~/ 3,
      NROWS ~/ 2,
      NROWS - 1
    ]);
    AbstractDoubleMatrix B = A.select(rowIndexes, colIndexes);
    expect(rowIndexes.size, equals(B.rows));
    expect(colIndexes.size, equals(B.columns));
    for (int r = 0; r < rowIndexes.size; r++) {
      for (int c = 0; c < colIndexes.size; c++) {
        expect(A.get(rowIndexes[r], colIndexes[c]), closeTo(B.get(r, c), TOL));
      }
    }
  }

  void testStrides() {
    int rowStride = 3;
    int colStride = 5;
    AbstractDoubleMatrix B = A.strides(rowStride, colStride);
    for (int r = 0; r < B.rows; r++) {
      for (int c = 0; c < B.columns; c++) {
        expect(A.get(r * rowStride, c * colStride), closeTo(B.get(r, c), TOL));
      }
    }
  }

  void testMultiply() {
    double alpha = 3.0;
    double beta = 5.0;
    AbstractDoubleMatrix C = new DiagonalDoubleMatrix(NROWS, NROWS, 0);
    for (int i = 0; i < DLENGTH; i++) {
      C.set(i, i, random.nextDouble());
    }
    List<Float64List> expected = C.toList();
    C = A.multiply(Bt, C, alpha, beta, false, false);
    for (int j = 0; j < NROWS; j++) {
      for (int i = 0; i < NROWS; i++) {
        double s = 0.0;
        for (int k = 0; k < NCOLUMNS; k++) {
          s += A.get(i, k) * Bt.get(k, j);
        }
        expected[i][j] = s * alpha + expected[i][j] * beta;
      }
    }
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NROWS; c++) {
        expect(expected[r][c], closeTo(C.get(r, c), TOL));
      }
    }

    //---
    C = null;
    C = A.multiply(Bt, C, alpha, beta, false, false);
    expected =
        new List<Float64List>.generate(NROWS, (_) => new Float64List(NROWS));
    for (int j = 0; j < NROWS; j++) {
      for (int i = 0; i < NROWS; i++) {
        double s = 0.0;
        for (int k = 0; k < NCOLUMNS; k++) {
          s += A.get(i, k) * Bt.get(k, j);
        }
        expected[i][j] = s * alpha;
      }
    }
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NROWS; c++) {
        expect(expected[r][c], closeTo(C.get(r, c), TOL));
      }
    }

    //transposeA
    C = new DiagonalDoubleMatrix(NCOLUMNS, NCOLUMNS, 0);
    for (int i = 0; i < DLENGTH; i++) {
      C.set(i, i, random.nextDouble());
    }
    expected = C.toList();
    C = A.multiply(B, C, alpha, beta, true, false);
    for (int j = 0; j < NCOLUMNS; j++) {
      for (int i = 0; i < NCOLUMNS; i++) {
        double s = 0.0;
        for (int k = 0; k < NROWS; k++) {
          s += A.get(k, i) * B.get(k, j);
        }
        expected[i][j] = s * alpha + expected[i][j] * beta;
      }
    }
    for (int r = 0; r < NCOLUMNS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(expected[r][c], closeTo(C.get(r, c), TOL));
      }
    }
    //---
    C = null;
    C = A.multiply(B, C, alpha, beta, true, false);
    expected = new List<Float64List>.generate(
        NCOLUMNS, (_) => new Float64List(NCOLUMNS));
    for (int j = 0; j < NCOLUMNS; j++) {
      for (int i = 0; i < NCOLUMNS; i++) {
        double s = 0.0;
        for (int k = 0; k < NROWS; k++) {
          s += A.get(k, i) * B.get(k, j);
        }
        expected[i][j] = s * alpha;
      }
    }
    for (int r = 0; r < NCOLUMNS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(expected[r][c], closeTo(C.get(r, c), TOL));
      }
    }

    //transposeB
    C = new DiagonalDoubleMatrix(NROWS, NROWS, 0);
    for (int i = 0; i < DLENGTH; i++) {
      C.set(i, i, random.nextDouble());
    }
    expected = C.toList();
    C = A.multiply(B, C, alpha, beta, false, true);
    for (int j = 0; j < NROWS; j++) {
      for (int i = 0; i < NROWS; i++) {
        double s = 0.0;
        for (int k = 0; k < NCOLUMNS; k++) {
          s += A.get(i, k) * B.get(j, k);
        }
        expected[i][j] = s * alpha + expected[i][j] * beta;
      }
    }
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NROWS; c++) {
        expect(expected[r][c], closeTo(C.get(r, c), TOL));
      }
    }
    //---
    C = null;
    C = A.multiply(B, C, alpha, beta, false, true);
    expected =
        new List<Float64List>.generate(NROWS, (_) => new Float64List(NROWS));
    for (int j = 0; j < NROWS; j++) {
      for (int i = 0; i < NROWS; i++) {
        double s = 0.0;
        for (int k = 0; k < NCOLUMNS; k++) {
          s += A.get(i, k) * B.get(j, k);
        }
        expected[i][j] = s * alpha;
      }
    }
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NROWS; c++) {
        expect(expected[r][c], closeTo(C.get(r, c), TOL));
      }
    }
    //transposeA and transposeB
    C = new DiagonalDoubleMatrix(NCOLUMNS, NCOLUMNS, 0);
    for (int i = 0; i < DLENGTH; i++) {
      C.set(i, i, random.nextDouble());
    }
    expected = C.toList();
    C = A.multiply(Bt, C, alpha, beta, true, true);
    for (int j = 0; j < NCOLUMNS; j++) {
      for (int i = 0; i < NCOLUMNS; i++) {
        double s = 0.0;
        for (int k = 0; k < NROWS; k++) {
          s += A.get(k, i) * Bt.get(j, k);
        }
        expected[i][j] = s * alpha + expected[i][j] * beta;
      }
    }
    for (int r = 0; r < NCOLUMNS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(expected[r][c], closeTo(C.get(r, c), TOL));
      }
    }
    //---
    C = null;
    C = A.multiply(Bt, C, alpha, beta, true, true);
    expected = new List<Float64List>.generate(
        NCOLUMNS, (_) => new Float64List(NCOLUMNS));
    for (int j = 0; j < NCOLUMNS; j++) {
      for (int i = 0; i < NCOLUMNS; i++) {
        double s = 0.0;
        for (int k = 0; k < NROWS; k++) {
          s += A.get(k, i) * Bt.get(j, k);
        }
        expected[i][j] = s * alpha;
      }
    }
    for (int r = 0; r < NCOLUMNS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(expected[r][c], closeTo(C.get(r, c), TOL));
      }
    }
  }
}

class DiagonalDoubleMatrixViewTest extends DiagonalDoubleMatrixTest {
  void createMatrices() {
    DINDEX = 3;
    A = new DiagonalDoubleMatrix(NCOLUMNS, NROWS, -DINDEX);
    DLENGTH = (A as DiagonalDoubleMatrix).diagonalLength;
    A = A.dice();
    B = new DiagonalDoubleMatrix(NCOLUMNS, NROWS, -DINDEX).dice();
    Bt = new DiagonalDoubleMatrix(NROWS, NCOLUMNS, DINDEX).dice();
  }
}
