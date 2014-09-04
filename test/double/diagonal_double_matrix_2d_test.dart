part of cern.colt.matrix.double.test;

class DiagonalDoubleMatrix2DTest extends DoubleMatrix2DTest {

  int DLENGTH;

  int DINDEX;

  void createMatrices() {
    DINDEX = 3;
    A = new DiagonalDoubleMatrix2D(NROWS, NCOLUMNS, DINDEX);
    B = new DiagonalDoubleMatrix2D(NROWS, NCOLUMNS, DINDEX);
    Bt = new DiagonalDoubleMatrix2D(NCOLUMNS, NROWS, -DINDEX);
    DLENGTH = (A as DiagonalDoubleMatrix2D).diagonalLength();
  }

  void populateMatrices() {
    //ConcurrencyUtils.setThreadsBeginN_2D(1);
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

  void testAssignValue() {
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

  void testAssignValues() {
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

  void testAssignValues2D() {
    List<Float64List> expected = new List<Float64List>(NROWS);//[NCOLUMNS];
    for (int r = 0; r < NROWS; r++) {
      expected[r] = new Float64List(NCOLUMNS);
      for (int c = 0; c < NCOLUMNS; c++) {
        expected[r][c] = random.nextDouble();
      }
    }
    A.setAll2D(expected);
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        expect(expected[r][r + DINDEX], closeTo(A.get(r, r + DINDEX), TOL));
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        expect(expected[r - DINDEX][r], closeTo(A.get(r - DINDEX, r), TOL));
      }
    }
  }

  void testAssign() {
    DoubleMatrix2D Acopy = A.copy();
    A.forEach(acos);
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

  void testAssignFunc() {
    DoubleMatrix2D Acopy = A.copy();
    A.forEachMatrix(B, div);
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        expect(Acopy.get(r, r + DINDEX) / B.get(r, r + DINDEX), closeTo(A.get(r, r + DINDEX), TOL));
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        expect(Acopy.get(r - DINDEX, r) / B.get(r - DINDEX, r), closeTo(A.get(r - DINDEX, r), TOL));
      }
    }
  }

  void testAssignFuncIndex() {
    List<int> rowList = new List<int>();
    List<int> columnList = new List<int>();
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        rowList.add(r);
        columnList.add(r + DINDEX);
      }
      DoubleMatrix2D Acopy = A.copy();
      A.forEachMatrixRange(B, div,
          new Int32List.fromList(rowList),
          new Int32List.fromList(columnList));
      for (int r = 0; r < DLENGTH; r++) {
        expect(Acopy.get(r, r + DINDEX) / B.get(r, r + DINDEX), closeTo(A.get(r, r + DINDEX), TOL));
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        rowList.add(r - DINDEX);
        columnList.add(r);
      }
      DoubleMatrix2D Acopy = A.copy();
      A.forEachMatrixRange(B, div, new Int32List.fromList(rowList),
          new Int32List.fromList(columnList));
      for (int r = 0; r < DLENGTH; r++) {
        expect(Acopy.get(r - DINDEX, r) / B.get(r - DINDEX, r), closeTo(A.get(r - DINDEX, r), TOL));
      }
    }
  }

  void testCardinality() {
    int card = A.cardinality;
    expect(DLENGTH, equals(card));
  }

  void testMaxLocation() {
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

  void testMinLocation() {
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

  void testGetNegativeValues() {
    A.fill(0.0);
    if (DINDEX >= 0) {
      A.set(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, -0.7);
      A.set(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, -0.1);
      List<int> rowList = new List<int>();
      List<int> columnList = new List<int>();
      List<double> valueList = new List<double>();
      A.negativeValues(rowList, columnList, valueList);
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
      A.negativeValues(rowList, columnList, valueList);
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

  void testGetNonZeros() {
    A.fill(0.0);
    if (DINDEX >= 0) {
      A.set(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, 0.7);
      A.set(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, 0.1);
      List<int> rowList = new List<int>();
      List<int> columnList = new List<int>();
      List<double> valueList = new List<double>();
      A.nonZeros(rowList, columnList, valueList);
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
      A.nonZeros(rowList, columnList, valueList);
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

  void testGetPositiveValues() {
    A.fill(0.0);
    if (DINDEX >= 0) {
      A.set(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, 0.7);
      A.set(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, 0.1);
      List<int> rowList = new List<int>();
      List<int> columnList = new List<int>();
      List<double> valueList = new List<double>();
      A.positiveValues(rowList, columnList, valueList);
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
      A.positiveValues(rowList, columnList, valueList);
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

  void testToArray() {
    List<Float64List> array = A.toList();
    expect(NROWS == array.length, isTrue);
    for (int r = 0; r < NROWS; r++) {
      expect(NCOLUMNS == array[r].length, isTrue);
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(array[r][c], closeTo(A.get(r, c), TOL));
      }
    }
  }

  void testVectorize() {
    DoubleMatrix1D Avec = A.vectorize();
    int idx = 0;
    for (int c = 0; c < NCOLUMNS; c++) {
      for (int r = 0; r < NROWS; r++) {
        expect(A.get(r, c), closeTo(Avec.get(idx++), TOL));
      }
    }
  }

  void testViewColumn() {
    DoubleMatrix1D col = A.column(NCOLUMNS ~/ 2);
    expect(NROWS, equals(col.length));
    for (int r = 0; r < NROWS; r++) {
      expect(A.get(r, NCOLUMNS ~/ 2), closeTo(col.get(r), TOL));
    }
  }

  void testViewColumnFlip() {
    DoubleMatrix2D B = A.columnFlip();
    expect(A.length, equals(B.length));
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(A.get(r, NCOLUMNS - 1 - c), closeTo(B.get(r, c), TOL));
      }
    }
  }

  void testViewDice() {
    DoubleMatrix2D B = A.dice();
    expect(NROWS, equals(B.columns));
    expect(NCOLUMNS, equals(B.rows));
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(A.get(r, c), closeTo(B.get(c, r), TOL));
      }
    }
  }

  void testViewPart() {
    DoubleMatrix2D B = A.part(NROWS ~/ 2, NCOLUMNS ~/ 2, NROWS ~/ 3, NCOLUMNS ~/ 3);
    expect(NROWS ~/ 3, equals(B.rows));
    expect(NCOLUMNS ~/ 3, equals(B.columns));
    for (int r = 0; r < NROWS ~/ 3; r++) {
      for (int c = 0; c < NCOLUMNS ~/ 3; c++) {
        expect(A.get(NROWS ~/ 2 + r, NCOLUMNS ~/ 2 + c), closeTo(B.get(r, c), TOL));
      }
    }
  }

  void testViewRow() {
    DoubleMatrix1D B = A.row(NROWS ~/ 2);
    expect(NCOLUMNS, equals(B.length));
    for (int r = 0; r < NCOLUMNS; r++) {
      expect(A.get(NROWS ~/ 2, r), closeTo(B.get(r), TOL));
    }
  }

  void testViewRowFlip() {
    DoubleMatrix2D B = A.rowFlip();
    expect(A.length, equals(B.length));
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(A.get(NROWS - 1 - r, c), closeTo(B.get(r, c), TOL));
      }
    }
  }

  void testViewSelectionProc() {
    final double value = 2.0;
    A.fill(0.0);
    if (DINDEX >= 0) {
      A.set(NROWS ~/ 4, NROWS ~/ 4 + DINDEX, value);
      A.set(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, value);
      DoubleMatrix2D B = A.where((DoubleMatrix1D element) {
        if ((element.get(NROWS ~/ 4 + DINDEX) - value).abs() < TOL) {
          return true;
        } else {
          return false;
        }
      });
      expect(1, equals(B.rows));
      expect(NCOLUMNS, equals(B.columns));
      expect(A.get(NROWS ~/ 4, NROWS ~/ 4 + DINDEX), closeTo(B.get(0, NROWS ~/ 4 + DINDEX), TOL));
    } else {
      A.set(NROWS ~/ 4 - DINDEX, NROWS ~/ 4, value);
      A.set(NROWS ~/ 2 - DINDEX, NROWS ~/ 2, value);
      DoubleMatrix2D B = A.where((DoubleMatrix1D element) {
        if ((element.get(NROWS ~/ 4) - value).abs() < TOL) {
          return true;
        } else {
          return false;
        }
      });
      expect(1, equals(B.rows));
      expect(NCOLUMNS, equals(B.columns));
      expect(A.get(NROWS ~/ 4 - DINDEX, NROWS ~/ 4), closeTo(B.get(0, NROWS ~/ 4), TOL));
    }
  }

  void testViewSelection() {
    final rowIndexes = new Int32List.fromList([NROWS ~/ 6, NROWS ~/ 5, NROWS ~/ 4, NROWS ~/ 3, NROWS ~/ 2]);
    final colIndexes = new Int32List.fromList([NROWS ~/ 6, NROWS ~/ 5, NROWS ~/ 4, NROWS ~/ 3, NROWS ~/ 2, NROWS - 1]);
    DoubleMatrix2D B = A.select(rowIndexes, colIndexes);
    expect(rowIndexes.length, equals(B.rows));
    expect(colIndexes.length, equals(B.columns));
    for (int r = 0; r < rowIndexes.length; r++) {
      for (int c = 0; c < colIndexes.length; c++) {
        expect(A.get(rowIndexes[r], colIndexes[c]), closeTo(B.get(r, c), TOL));
      }
    }
  }

  /*void testViewSorted() {
    DoubleMatrix2D B = A.viewSorted(1);
    for (int r = 0; r < NROWS - 1; r++) {
      expect(B.getQuick(r + 1, 1) >= B.getQuick(r, 1), isTrue);
    }
  }*/

  void testViewStrides() {
    int rowStride = 3;
    int colStride = 5;
    DoubleMatrix2D B = A.strides(rowStride, colStride);
    for (int r = 0; r < B.rows; r++) {
      for (int c = 0; c < B.columns; c++) {
        expect(A.get(r * rowStride, c * colStride), closeTo(B.get(r, c), TOL));
      }
    }
  }

  void testZMult2D() {
    double alpha = 3.0;
    double beta = 5.0;
    DoubleMatrix2D C = new DiagonalDoubleMatrix2D(NROWS, NROWS, 0);
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
    expected = new List<Float64List>.generate(NROWS,
        (_) => new Float64List(NROWS));
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
    C = new DiagonalDoubleMatrix2D(NCOLUMNS, NCOLUMNS, 0);
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
    expected = new List<Float64List>.generate(NCOLUMNS, 
        (_) => new Float64List(NCOLUMNS));
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
    C = new DiagonalDoubleMatrix2D(NROWS, NROWS, 0);
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
    expected = new List<Float64List>.generate(NROWS, 
        (_) => new Float64List(NROWS));
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
    C = new DiagonalDoubleMatrix2D(NCOLUMNS, NCOLUMNS, 0);
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
    expected = new List<Float64List>.generate(NCOLUMNS, 
        (_) => new Float64List(NCOLUMNS));
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

class DiagonalDoubleMatrix2DViewTest extends DiagonalDoubleMatrix2DTest {

  void createMatrices() {
    DINDEX = 3;
    A = new DiagonalDoubleMatrix2D(NCOLUMNS, NROWS, -DINDEX);
    DLENGTH = (A as DiagonalDoubleMatrix2D).diagonalLength();
    A = A.dice();
    B = new DiagonalDoubleMatrix2D(NCOLUMNS, NROWS, -DINDEX).dice();
    Bt = new DiagonalDoubleMatrix2D(NROWS, NCOLUMNS, DINDEX).dice();
  }

}
