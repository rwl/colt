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
        A.setQuick(r, r + DINDEX, random.nextDouble());
      }

      for (int r = 0; r < DLENGTH; r++) {
        B.setQuick(r, r + DINDEX, random.nextDouble());
      }

      for (int r = 0; r < DLENGTH; r++) {
        Bt.setQuick(r - DINDEX, r, random.nextDouble());
      }

    } else {
      for (int r = 0; r < DLENGTH; r++) {
        A.setQuick(r - DINDEX, r, random.nextDouble());
      }

      for (int r = 0; r < DLENGTH; r++) {
        B.setQuick(r - DINDEX, r, random.nextDouble());
      }
      for (int r = 0; r < DLENGTH; r++) {
        Bt.setQuick(r, r + DINDEX, random.nextDouble());
      }

    }
  }

  void testAssignValue() {
    double value = random.nextDouble();
    A.assignValue(value);
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        expect(value, closeTo(A.getQuick(r, r + DINDEX), TOL));
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        expect(value, closeTo(A.getQuick(r - DINDEX, r), TOL));
      }
    }
  }

  void testAssignValues() {
    Float64List expected = new Float64List(DLENGTH);
    for (int i = 0; i < DLENGTH; i++) {
      expected[i] = random.nextDouble();
    }
    A.assignValues(expected);
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        expect(expected[r], closeTo(A.getQuick(r, r + DINDEX), TOL));
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        expect(expected[r], closeTo(A.getQuick(r - DINDEX, r), TOL));
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
    A.assignValues2D(expected);
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        expect(expected[r][r + DINDEX], closeTo(A.getQuick(r, r + DINDEX), TOL));
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        expect(expected[r - DINDEX][r], closeTo(A.getQuick(r - DINDEX, r), TOL));
      }
    }
  }

  void testAssign() {
    DoubleMatrix2D Acopy = A.copy();
    A.assign(acos);
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        double expected = math.acos(Acopy.getQuick(r, r + DINDEX));
        expect(expected, closeTo(A.getQuick(r, r + DINDEX), TOL));
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        double expected = math.acos(Acopy.getQuick(r - DINDEX, r));
        expect(expected, closeTo(A.getQuick(r - DINDEX, r), TOL));
      }
    }
  }

  void testAssignFunc() {
    DoubleMatrix2D Acopy = A.copy();
    A.assignFunc(B, div);
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        expect(Acopy.getQuick(r, r + DINDEX) / B.getQuick(r, r + DINDEX), closeTo(A.getQuick(r, r + DINDEX), TOL));
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        expect(Acopy.getQuick(r - DINDEX, r) / B.getQuick(r - DINDEX, r), closeTo(A.getQuick(r - DINDEX, r), TOL));
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
      A.assignFuncIndex(B, div, rowList, columnList);
      for (int r = 0; r < DLENGTH; r++) {
        expect(Acopy.getQuick(r, r + DINDEX) / B.getQuick(r, r + DINDEX), closeTo(A.getQuick(r, r + DINDEX), TOL));
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        rowList.add(r - DINDEX);
        columnList.add(r);
      }
      DoubleMatrix2D Acopy = A.copy();
      A.assignFuncIndex(B, div, rowList, columnList);
      for (int r = 0; r < DLENGTH; r++) {
        expect(Acopy.getQuick(r - DINDEX, r) / B.getQuick(r - DINDEX, r), closeTo(A.getQuick(r - DINDEX, r), TOL));
      }
    }
  }

  void testCardinality() {
    int card = A.cardinality();
    expect(DLENGTH, equals(card));
  }

  void testMaxLocation() {
    A.assignValue(0.0);
    if (DINDEX >= 0) {
      A.setQuick(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, 0.7);
      A.setQuick(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, 0.1);
      Float64List maxAndLoc = A.getMaxLocation();
      expect(0.7, closeTo(maxAndLoc[0], TOL));
      expect(NROWS / 3, equals(maxAndLoc[1]));
      expect(NROWS / 3 + DINDEX, equals(maxAndLoc[2]));
    } else {
      A.setQuick(NROWS ~/ 3 - DINDEX, NROWS ~/ 3, 0.7);
      A.setQuick(NROWS ~/ 2 - DINDEX, NROWS ~/ 2, 0.1);
      Float64List maxAndLoc = A.getMaxLocation();
      expect(0.7, closeTo(maxAndLoc[0], TOL));
      expect(NROWS ~/ 3 - DINDEX, equals(maxAndLoc[1]));
      expect(NROWS ~/ 3, equals(maxAndLoc[2]));
    }
  }

  void testMinLocation() {
    A.assignValue(0.0);
    if (DINDEX >= 0) {
      A.setQuick(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, -0.7);
      A.setQuick(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, -0.1);
      Float64List minAndLoc = A.getMinLocation();
      expect(-0.7, closeTo(minAndLoc[0], TOL));
      expect(NROWS / 3, equals(minAndLoc[1]));
      expect(NROWS / 3 + DINDEX, equals(minAndLoc[2]));
    } else {
      A.setQuick(NROWS ~/ 3 - DINDEX, NROWS ~/ 3, -0.7);
      A.setQuick(NROWS ~/ 2 - DINDEX, NROWS ~/ 2, -0.1);
      Float64List minAndLoc = A.getMinLocation();
      expect(-0.7, closeTo(minAndLoc[0], TOL));
      expect(NROWS ~/ 3 - DINDEX, equals(minAndLoc[1]));
      expect(NROWS ~/ 3, equals(minAndLoc[2]));
    }
  }

  void testGetNegativeValues() {
    A.assignValue(0.0);
    if (DINDEX >= 0) {
      A.setQuick(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, -0.7);
      A.setQuick(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, -0.1);
      List<int> rowList = new List<int>();
      List<int> columnList = new List<int>();
      List<double> valueList = new List<double>();
      A.getNegativeValues(rowList, columnList, valueList);
      expect(2, equals(rowList.length));
      expect(2, equals(columnList.length));
      expect(2, equals(valueList.length));
      expect(rowList.contains(NROWS / 3), isTrue);
      expect(rowList.contains(NROWS / 2), isTrue);
      expect(columnList.contains(NROWS / 3 + DINDEX), isTrue);
      expect(columnList.contains(NROWS / 2 + DINDEX), isTrue);
      expect(valueList.contains(-0.7), isTrue);
      expect(valueList.contains(-0.1), isTrue);
    } else {
      A.setQuick(NROWS ~/ 3 - DINDEX, NROWS ~/ 3, -0.7);
      A.setQuick(NROWS ~/ 2 - DINDEX, NROWS ~/ 2, -0.1);
      List<int> rowList = new List<int>();
      List<int> columnList = new List<int>();
      List<double> valueList = new List<double>();
      A.getNegativeValues(rowList, columnList, valueList);
      expect(2, equals(rowList.length));
      expect(2, equals(columnList.length));
      expect(2, equals(valueList.length));
      expect(rowList.contains(NROWS / 3 - DINDEX), isTrue);
      expect(rowList.contains(NROWS / 2 - DINDEX), isTrue);
      expect(columnList.contains(NROWS / 3), isTrue);
      expect(columnList.contains(NROWS / 2), isTrue);
      expect(valueList.contains(-0.7), isTrue);
      expect(valueList.contains(-0.1), isTrue);
    }
  }

  void testGetNonZeros() {
    A.assignValue(0.0);
    if (DINDEX >= 0) {
      A.setQuick(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, 0.7);
      A.setQuick(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, 0.1);
      List<int> rowList = new List<int>();
      List<int> columnList = new List<int>();
      List<double> valueList = new List<double>();
      A.getNonZeros(rowList, columnList, valueList);
      expect(2, equals(rowList.length));
      expect(2, equals(columnList.length));
      expect(2, equals(valueList.length));
      expect(rowList.contains(NROWS / 3), isTrue);
      expect(rowList.contains(NROWS / 2), isTrue);
      expect(columnList.contains(NROWS / 3 + DINDEX), isTrue);
      expect(columnList.contains(NROWS / 2 + DINDEX), isTrue);
      expect(valueList.contains(0.7), isTrue);
      expect(valueList.contains(0.1), isTrue);
    } else {
      A.setQuick(NROWS ~/ 3 - DINDEX, NROWS ~/ 3, 0.7);
      A.setQuick(NROWS ~/ 2 - DINDEX, NROWS ~/ 2, 0.1);
      List<int> rowList = new List<int>();
      List<int> columnList = new List<int>();
      List<double> valueList = new List<double>();
      A.getNonZeros(rowList, columnList, valueList);
      expect(2, equals(rowList.length));
      expect(2, equals(columnList.length));
      expect(2, equals(valueList.length));
      expect(rowList.contains(NROWS / 3 - DINDEX), isTrue);
      expect(rowList.contains(NROWS / 2 - DINDEX), isTrue);
      expect(columnList.contains(NROWS / 3), isTrue);
      expect(columnList.contains(NROWS / 2), isTrue);
      expect(valueList.contains(0.7), isTrue);
      expect(valueList.contains(0.1), isTrue);
    }
  }

  void testGetPositiveValues() {
    A.assignValue(0.0);
    if (DINDEX >= 0) {
      A.setQuick(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, 0.7);
      A.setQuick(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, 0.1);
      List<int> rowList = new List<int>();
      List<int> columnList = new List<int>();
      List<double> valueList = new List<double>();
      A.getPositiveValues(rowList, columnList, valueList);
      expect(2, equals(rowList.length));
      expect(2, equals(columnList.length));
      expect(2, equals(valueList.length));
      expect(rowList.contains(NROWS / 3), isTrue);
      expect(rowList.contains(NROWS / 2), isTrue);
      expect(columnList.contains(NROWS / 3 + DINDEX), isTrue);
      expect(columnList.contains(NROWS / 2 + DINDEX), isTrue);
      expect(valueList.contains(0.7), isTrue);
      expect(valueList.contains(0.1), isTrue);
    } else {
      A.setQuick(NROWS ~/ 3 - DINDEX, NROWS ~/ 3, 0.7);
      A.setQuick(NROWS ~/ 2 - DINDEX, NROWS ~/ 2, 0.1);
      List<int> rowList = new List<int>();
      List<int> columnList = new List<int>();
      List<double> valueList = new List<double>();
      A.getPositiveValues(rowList, columnList, valueList);
      expect(2, equals(rowList.length));
      expect(2, equals(columnList.length));
      expect(2, equals(valueList.length));
      expect(rowList.contains(NROWS / 3 - DINDEX), isTrue);
      expect(rowList.contains(NROWS / 2 - DINDEX), isTrue);
      expect(columnList.contains(NROWS / 3), isTrue);
      expect(columnList.contains(NROWS / 2), isTrue);
      expect(valueList.contains(0.7), isTrue);
      expect(valueList.contains(0.1), isTrue);
    }
  }

  void testToArray() {
    List<Float64List> array = A.toArray();
    expect(NROWS == array.length, isTrue);
    for (int r = 0; r < NROWS; r++) {
      expect(NCOLUMNS == array[r].length, isTrue);
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(array[r][c], closeTo(A.getQuick(r, c), TOL));
      }
    }
  }

  void testVectorize() {
    DoubleMatrix1D Avec = A.vectorize();
    int idx = 0;
    for (int c = 0; c < NCOLUMNS; c++) {
      for (int r = 0; r < NROWS; r++) {
        expect(A.getQuick(r, c), closeTo(Avec.getQuick(idx++), TOL));
      }
    }
  }

  void testViewColumn() {
    DoubleMatrix1D col = A.viewColumn(NCOLUMNS ~/ 2);
    expect(NROWS, equals(col.size()));
    for (int r = 0; r < NROWS; r++) {
      expect(A.getQuick(r, NCOLUMNS ~/ 2), closeTo(col.getQuick(r), TOL));
    }
  }

  void testViewColumnFlip() {
    DoubleMatrix2D B = A.viewColumnFlip();
    expect(A.size(), equals(B.size()));
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(A.getQuick(r, NCOLUMNS - 1 - c), closeTo(B.getQuick(r, c), TOL));
      }
    }
  }

  void testViewDice() {
    DoubleMatrix2D B = A.viewDice();
    expect(NROWS, equals(B.columns()));
    expect(NCOLUMNS, equals(B.rows()));
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(A.getQuick(r, c), closeTo(B.getQuick(c, r), TOL));
      }
    }
  }

  void testViewPart() {
    DoubleMatrix2D B = A.viewPart(NROWS ~/ 2, NCOLUMNS ~/ 2, NROWS ~/ 3, NCOLUMNS ~/ 3);
    expect(NROWS ~/ 3, equals(B.rows()));
    expect(NCOLUMNS ~/ 3, equals(B.columns()));
    for (int r = 0; r < NROWS ~/ 3; r++) {
      for (int c = 0; c < NCOLUMNS ~/ 3; c++) {
        expect(A.getQuick(NROWS ~/ 2 + r, NCOLUMNS ~/ 2 + c), closeTo(B.getQuick(r, c), TOL));
      }
    }
  }

  void testViewRow() {
    DoubleMatrix1D B = A.viewRow(NROWS ~/ 2);
    expect(NCOLUMNS, equals(B.size()));
    for (int r = 0; r < NCOLUMNS; r++) {
      expect(A.getQuick(NROWS ~/ 2, r), closeTo(B.getQuick(r), TOL));
    }
  }

  void testViewRowFlip() {
    DoubleMatrix2D B = A.viewRowFlip();
    expect(A.size(), equals(B.size()));
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(A.getQuick(NROWS - 1 - r, c), closeTo(B.getQuick(r, c), TOL));
      }
    }
  }

  void testViewSelectionProc() {
    final double value = 2.0;
    A.assignValue(0.0);
    if (DINDEX >= 0) {
      A.setQuick(NROWS ~/ 4, NROWS ~/ 4 + DINDEX, value);
      A.setQuick(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, value);
      DoubleMatrix2D B = A.viewSelectionProc((DoubleMatrix1D element) {
        if ((element.getQuick(NROWS ~/ 4 + DINDEX) - value).abs() < TOL) {
          return true;
        } else {
          return false;
        }
      });
      expect(1, equals(B.rows()));
      expect(NCOLUMNS, equals(B.columns()));
      expect(A.getQuick(NROWS ~/ 4, NROWS ~/ 4 + DINDEX), closeTo(B.getQuick(0, NROWS ~/ 4 + DINDEX), TOL));
    } else {
      A.setQuick(NROWS ~/ 4 - DINDEX, NROWS ~/ 4, value);
      A.setQuick(NROWS ~/ 2 - DINDEX, NROWS ~/ 2, value);
      DoubleMatrix2D B = A.viewSelectionProc((DoubleMatrix1D element) {
        if ((element.getQuick(NROWS ~/ 4) - value).abs() < TOL) {
          return true;
        } else {
          return false;
        }
      });
      expect(1, equals(B.rows()));
      expect(NCOLUMNS, equals(B.columns()));
      expect(A.getQuick(NROWS ~/ 4 - DINDEX, NROWS ~/ 4), closeTo(B.getQuick(0, NROWS ~/ 4), TOL));
    }
  }

  void testViewSelection() {
    final rowIndexes = new Int32List.fromList([NROWS / 6, NROWS / 5, NROWS / 4, NROWS / 3, NROWS / 2]);
    final colIndexes = new Int32List.fromList([NROWS / 6, NROWS / 5, NROWS / 4, NROWS / 3, NROWS / 2, NROWS - 1]);
    DoubleMatrix2D B = A.viewSelection(rowIndexes, colIndexes);
    expect(rowIndexes.length, equals(B.rows()));
    expect(colIndexes.length, equals(B.columns()));
    for (int r = 0; r < rowIndexes.length; r++) {
      for (int c = 0; c < colIndexes.length; c++) {
        expect(A.getQuick(rowIndexes[r], colIndexes[c]), closeTo(B.getQuick(r, c), TOL));
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
    DoubleMatrix2D B = A.viewStrides(rowStride, colStride);
    for (int r = 0; r < B.rows(); r++) {
      for (int c = 0; c < B.columns(); c++) {
        expect(A.getQuick(r * rowStride, c * colStride), closeTo(B.getQuick(r, c), TOL));
      }
    }
  }

  void testZMult2D() {
    double alpha = 3.0;
    double beta = 5.0;
    DoubleMatrix2D C = new DiagonalDoubleMatrix2D(NROWS, NROWS, 0);
    for (int i = 0; i < DLENGTH; i++) {
      C.setQuick(i, i, random.nextDouble());
    }
    List<Float64List> expected = C.toArray();
    C = A.zMult2D(Bt, C, alpha, beta, false, false);
    for (int j = 0; j < NROWS; j++) {
      for (int i = 0; i < NROWS; i++) {
        double s = 0.0;
        for (int k = 0; k < NCOLUMNS; k++) {
          s += A.getQuick(i, k) * Bt.getQuick(k, j);
        }
        expected[i][j] = s * alpha + expected[i][j] * beta;
      }
    }
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NROWS; c++) {
        expect(expected[r][c], closeTo(C.getQuick(r, c), TOL));
      }
    }

    //---
    C = null;
    C = A.zMult2D(Bt, C, alpha, beta, false, false);
    expected = new List<Float64List>(NROWS);//[NROWS];
    for (int j = 0; j < NROWS; j++) {
      expected[j] = new Float64List(NROWS);
      for (int i = 0; i < NROWS; i++) {
        double s = 0.0;
        for (int k = 0; k < NCOLUMNS; k++) {
          s += A.getQuick(i, k) * Bt.getQuick(k, j);
        }
        expected[i][j] = s * alpha;
      }
    }
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NROWS; c++) {
        expect(expected[r][c], closeTo(C.getQuick(r, c), TOL));
      }
    }

    //transposeA
    C = new DiagonalDoubleMatrix2D(NCOLUMNS, NCOLUMNS, 0);
    for (int i = 0; i < DLENGTH; i++) {
      C.setQuick(i, i, random.nextDouble());
    }
    expected = C.toArray();
    C = A.zMult2D(B, C, alpha, beta, true, false);
    for (int j = 0; j < NCOLUMNS; j++) {
      for (int i = 0; i < NCOLUMNS; i++) {
        double s = 0.0;
        for (int k = 0; k < NROWS; k++) {
          s += A.getQuick(k, i) * B.getQuick(k, j);
        }
        expected[i][j] = s * alpha + expected[i][j] * beta;
      }
    }
    for (int r = 0; r < NCOLUMNS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(expected[r][c], closeTo(C.getQuick(r, c), TOL));
      }
    }
    //---
    C = null;
    C = A.zMult2D(B, C, alpha, beta, true, false);
    expected = new List<Float64List>(NCOLUMNS);//[NCOLUMNS];
    for (int j = 0; j < NCOLUMNS; j++) {
      expected[j] = new Float64List(NCOLUMNS);
      for (int i = 0; i < NCOLUMNS; i++) {
        double s = 0.0;
        for (int k = 0; k < NROWS; k++) {
          s += A.getQuick(k, i) * B.getQuick(k, j);
        }
        expected[i][j] = s * alpha;
      }
    }
    for (int r = 0; r < NCOLUMNS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(expected[r][c], closeTo(C.getQuick(r, c), TOL));
      }
    }

    //transposeB
    C = new DiagonalDoubleMatrix2D(NROWS, NROWS, 0);
    for (int i = 0; i < DLENGTH; i++) {
      C.setQuick(i, i, random.nextDouble());
    }
    expected = C.toArray();
    C = A.zMult2D(B, C, alpha, beta, false, true);
    for (int j = 0; j < NROWS; j++) {
      for (int i = 0; i < NROWS; i++) {
        double s = 0.0;
        for (int k = 0; k < NCOLUMNS; k++) {
          s += A.getQuick(i, k) * B.getQuick(j, k);
        }
        expected[i][j] = s * alpha + expected[i][j] * beta;
      }
    }
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NROWS; c++) {
        expect(expected[r][c], closeTo(C.getQuick(r, c), TOL));
      }
    }
    //---
    C = null;
    C = A.zMult2D(B, C, alpha, beta, false, true);
    expected = new List<Float64List>(NROWS);//[NROWS];
    for (int j = 0; j < NROWS; j++) {
      expected[j] = new Float64List(NROWS);
      for (int i = 0; i < NROWS; i++) {
        double s = 0.0;
        for (int k = 0; k < NCOLUMNS; k++) {
          s += A.getQuick(i, k) * B.getQuick(j, k);
        }
        expected[i][j] = s * alpha;
      }
    }
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NROWS; c++) {
        expect(expected[r][c], closeTo(C.getQuick(r, c), TOL));
      }
    }
    //transposeA and transposeB
    C = new DiagonalDoubleMatrix2D(NCOLUMNS, NCOLUMNS, 0);
    for (int i = 0; i < DLENGTH; i++) {
      C.setQuick(i, i, random.nextDouble());
    }
    expected = C.toArray();
    C = A.zMult2D(Bt, C, alpha, beta, true, true);
    for (int j = 0; j < NCOLUMNS; j++) {
      for (int i = 0; i < NCOLUMNS; i++) {
        double s = 0.0;
        for (int k = 0; k < NROWS; k++) {
          s += A.getQuick(k, i) * Bt.getQuick(j, k);
        }
        expected[i][j] = s * alpha + expected[i][j] * beta;
      }
    }
    for (int r = 0; r < NCOLUMNS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(expected[r][c], closeTo(C.getQuick(r, c), TOL));
      }
    }
    //---
    C = null;
    C = A.zMult2D(Bt, C, alpha, beta, true, true);
    expected = new List<Float64List>(NCOLUMNS);//[NCOLUMNS];
    for (int j = 0; j < NCOLUMNS; j++) {
      expected[j] = new Float64List(NCOLUMNS);
      for (int i = 0; i < NCOLUMNS; i++) {
        double s = 0.0;
        for (int k = 0; k < NROWS; k++) {
          s += A.getQuick(k, i) * Bt.getQuick(j, k);
        }
        expected[i][j] = s * alpha;
      }
    }
    for (int r = 0; r < NCOLUMNS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(expected[r][c], closeTo(C.getQuick(r, c), TOL));
      }
    }

  }
}

class DiagonalDoubleMatrix2DViewTest extends DiagonalDoubleMatrix2DTest {

    void createMatrices() {
        DINDEX = 3;
        A = new DiagonalDoubleMatrix2D(NCOLUMNS, NROWS, -DINDEX);
        DLENGTH = (A as DiagonalDoubleMatrix2D).diagonalLength();
        A = A.viewDice();
        B = new DiagonalDoubleMatrix2D(NCOLUMNS, NROWS, -DINDEX).viewDice();
        Bt = new DiagonalDoubleMatrix2D(NROWS, NCOLUMNS, DINDEX).viewDice();
    }

}
