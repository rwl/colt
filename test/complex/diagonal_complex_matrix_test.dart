part of cern.colt.matrix.complex.test;

class DiagonalDComplexMatrix2DTest extends DComplexMatrix2DTest {

  int DLENGTH;

  int DINDEX;

  void createMatrices() {
    DINDEX = 3;
    A = new DiagonalDComplexMatrix2D(NROWS, NCOLUMNS, DINDEX);
    B = new DiagonalDComplexMatrix2D(NROWS, NCOLUMNS, DINDEX);
    Bt = new DiagonalDComplexMatrix2D(NCOLUMNS, NROWS, -DINDEX);
    DLENGTH = (A as DiagonalDComplexMatrix2D).diagonalLength();
  }

  void populateMatrices() {
    //ConcurrencyUtils.setThreadsBeginN_2D(1);
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        A.setPartsQuick(r, r + DINDEX, random.nextDouble(), random.nextDouble());
      }

      for (int r = 0; r < DLENGTH; r++) {
        B.setPartsQuick(r, r + DINDEX, random.nextDouble(), random.nextDouble());
      }

      for (int r = 0; r < DLENGTH; r++) {
        Bt.setPartsQuick(r - DINDEX, r, random.nextDouble(), random.nextDouble());
      }

    } else {
      for (int r = 0; r < DLENGTH; r++) {
        A.setPartsQuick(r - DINDEX, r, random.nextDouble(), random.nextDouble());
      }

      for (int r = 0; r < DLENGTH; r++) {
        B.setPartsQuick(r - DINDEX, r, random.nextDouble(), random.nextDouble());
      }
      for (int r = 0; r < DLENGTH; r++) {
        Bt.setPartsQuick(r, r + DINDEX, random.nextDouble(), random.nextDouble());
      }

    }
  }

  void testAssignValue() {
    Float64List value = new Float64List.fromList([random.nextDouble(), random.nextDouble()]);
    A.assignValue(value[0], value[1]);
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        assertEquals(value, A.getQuick(r, r + DINDEX), TOL);
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        assertEquals(value, A.getQuick(r - DINDEX, r), TOL);
      }
    }
  }

  void testAssignValues() {
    Float64List expected = new Float64List(2 * DLENGTH);
    for (int i = 0; i < 2 * DLENGTH; i++) {
      expected[i] = random.nextDouble();
    }
    A.assignValues(expected);
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        expect(expected[2 * r], closeTo(A.getQuick(r, r + DINDEX)[0], TOL));
        expect(expected[2 * r + 1], closeTo(A.getQuick(r, r + DINDEX)[1], TOL));
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        expect(expected[2 * r], closeTo(A.getQuick(r - DINDEX, r)[0], TOL));
        expect(expected[2 * r + 1], closeTo(A.getQuick(r - DINDEX, r)[1], TOL));
      }
    }
  }

  void testAssignImaginary() {
    DoubleMatrix2D Im = DoubleFactory2D.dense.random(A.rows(), A.columns());
    DComplexMatrix2D Acopy = A.copy();
    A.assignImaginary(Im);
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        expect(Acopy.getQuick(r, r + DINDEX)[0], closeTo(A.getQuick(r, r + DINDEX)[0], TOL));
        expect(Im.getQuick(r, r + DINDEX), closeTo(A.getQuick(r, r + DINDEX)[1], TOL));
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        expect(Acopy.getQuick(r - DINDEX, r)[0], closeTo(A.getQuick(r - DINDEX, r)[0], TOL));
        expect(Im.getQuick(r - DINDEX, r), closeTo(A.getQuick(r - DINDEX, r)[1], TOL));
      }
    }
  }

  void testAssignReal() {
    DoubleMatrix2D Re = DoubleFactory2D.dense.random(A.rows(), A.columns());
    DComplexMatrix2D Acopy = A.copy();
    A.assignReal(Re);
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        expect(Acopy.getQuick(r, r + DINDEX)[1], closeTo(A.getQuick(r, r + DINDEX)[1], TOL));
        expect(Re.getQuick(r, r + DINDEX), closeTo(A.getQuick(r, r + DINDEX)[0], TOL));
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        expect(Acopy.getQuick(r - DINDEX, r)[1], closeTo(A.getQuick(r - DINDEX, r)[1], TOL));
        expect(Re.getQuick(r - DINDEX, r), closeTo(A.getQuick(r - DINDEX, r)[0], TOL));
      }
    }
  }

  void testAssignList() {
    List<Float64List> expected = new List<Float64List>.generate(NROWS,
        (_) => new Float64List(2 * NCOLUMNS));
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expected[r][2 * c] = random.nextDouble();
        expected[r][2 * c + 1] = random.nextDouble();
      }
    }
    A.assignList(expected);
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        expect(expected[r][2 * (r + DINDEX)], closeTo(A.getQuick(r, r + DINDEX)[0], TOL));
        expect(expected[r][2 * (r + DINDEX) + 1], closeTo(A.getQuick(r, r + DINDEX)[1], TOL));
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        expect(expected[r - DINDEX][2 * r], closeTo(A.getQuick(r - DINDEX, r)[0], TOL));
        expect(expected[r - DINDEX][2 * r + 1], closeTo(A.getQuick(r - DINDEX, r)[1], TOL));
      }
    }
  }

  void testAssign() {
    DComplexMatrix2D Acopy = A.copy();
    A.assign(acos);
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        Float64List expected = DComplex.acos(Acopy.getQuick(r, r + DINDEX));
        expect(expected[0], closeTo(A.getQuick(r, r + DINDEX)[0], TOL));
        expect(expected[1], closeTo(A.getQuick(r, r + DINDEX)[1], TOL));
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        Float64List expected = DComplex.acos(Acopy.getQuick(r - DINDEX, r));
        expect(expected[0], closeTo(A.getQuick(r - DINDEX, r)[0], TOL));
        expect(expected[1], closeTo(A.getQuick(r - DINDEX, r)[1], TOL));
      }
    }
  }

  void testAssignFunc() {
    DComplexMatrix2D Acopy = A.copy();
    A.assignFunc(B, div);
    if (DINDEX >= 0) {
      for (int r = 0; r < DLENGTH; r++) {
        expect(DComplex.div_(Acopy.getQuick(r, r + DINDEX), B.getQuick(r, r + DINDEX))[0], closeTo(A.getQuick(r, r + DINDEX)[0], TOL));
        expect(DComplex.div_(Acopy.getQuick(r, r + DINDEX), B.getQuick(r, r + DINDEX))[1], closeTo(A.getQuick(r, r + DINDEX)[1], TOL));
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        expect(DComplex.div_(Acopy.getQuick(r - DINDEX, r), B.getQuick(r - DINDEX, r))[0], closeTo(A.getQuick(r - DINDEX, r)[0], TOL));
        expect(DComplex.div_(Acopy.getQuick(r - DINDEX, r), B.getQuick(r - DINDEX, r))[1], closeTo(A.getQuick(r - DINDEX, r)[1], TOL));
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
      DComplexMatrix2D Acopy = A.copy();
      A.assignFuncIndex(B, div, rowList, columnList);
      for (int r = 0; r < DLENGTH; r++) {
        expect(DComplex.div_(Acopy.getQuick(r, r + DINDEX), B.getQuick(r, r + DINDEX))[0], closeTo(A.getQuick(r, r + DINDEX)[0], TOL));
        expect(DComplex.div_(Acopy.getQuick(r, r + DINDEX), B.getQuick(r, r + DINDEX))[1], closeTo(A.getQuick(r, r + DINDEX)[1], TOL));
      }
    } else {
      for (int r = 0; r < DLENGTH; r++) {
        rowList.add(r - DINDEX);
        columnList.add(r);
      }
      DComplexMatrix2D Acopy = A.copy();
      A.assignFuncIndex(B, div, rowList, columnList);
      for (int r = 0; r < DLENGTH; r++) {
        expect(DComplex.div_(Acopy.getQuick(r - DINDEX, r), B.getQuick(r - DINDEX, r))[0], closeTo(A.getQuick(r - DINDEX, r)[0], TOL));
        expect(DComplex.div_(Acopy.getQuick(r - DINDEX, r), B.getQuick(r - DINDEX, r))[1], closeTo(A.getQuick(r - DINDEX, r)[1], TOL));
      }
    }
  }

  void testCardinality() {
    int card = A.cardinality();
    expect(DLENGTH, equals(card));
  }

  void testGetNonZeros() {
    A.assignValue(0.0, 0.0);
    Float64List elem1 = new Float64List.fromList([0.7, 0.8]);
    Float64List elem2 = new Float64List.fromList([0.1, 0.2]);
    if (DINDEX >= 0) {
      A.setQuick(NROWS ~/ 3, NROWS ~/ 3 + DINDEX, elem1);
      A.setQuick(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, elem2);
      List<int> rowList = new List<int>();
      List<int> columnList = new List<int>();
      List<Float64List> valueList = new List<Float64List>();
      A.getNonZeros(rowList, columnList, valueList);
      expect(2, equals(rowList.length));
      expect(2, equals(columnList.length));
      expect(2, equals(valueList.length));
      expect(rowList.contains(NROWS ~/ 3), isTrue);
      expect(rowList.contains(NROWS ~/ 2), isTrue);
      expect(columnList.contains(NROWS ~/ 3 + DINDEX), isTrue);
      expect(columnList.contains(NROWS ~/ 2 + DINDEX), isTrue);
      assertEquals(A.getQuick(rowList[0], columnList[0]), valueList[0], TOL);
      assertEquals(A.getQuick(rowList[1], columnList[1]), valueList[1], TOL);
    } else {
      A.setQuick(NROWS ~/ 3 - DINDEX, NROWS ~/ 3, elem1);
      A.setQuick(NROWS ~/ 2 - DINDEX, NROWS ~/ 2, elem2);
      List<int> rowList = new List<int>();
      List<int> columnList = new List<int>();
      List<Float64List> valueList = new List<Float64List>();
      A.getNonZeros(rowList, columnList, valueList);
      expect(2, equals(rowList.length));
      expect(2, equals(columnList.length));
      expect(2, equals(valueList.length));
      expect(rowList.contains(NROWS ~/ 3 - DINDEX), isTrue);
      expect(rowList.contains(NROWS ~/ 2 - DINDEX), isTrue);
      expect(columnList.contains(NROWS ~/ 3), isTrue);
      expect(columnList.contains(NROWS ~/ 2), isTrue);
      assertEquals(A.getQuick(rowList[0], columnList[0]), valueList[0], TOL);
      assertEquals(A.getQuick(rowList[1], columnList[1]), valueList[1], TOL);
    }
  }

  void testToArray() {
    List<Float64List> array = A.toArray();
    expect(NROWS == array.length, isTrue);
    for (int r = 0; r < NROWS; r++) {
      expect(2 * NCOLUMNS == array[r].length, isTrue);
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(array[r][2 * c], closeTo(A.getQuick(r, c)[0], TOL));
        expect(array[r][2 * c + 1], closeTo(A.getQuick(r, c)[1], TOL));
      }
    }
  }

  void testVectorize() {
    DComplexMatrix1D Avec = A.vectorize();
    int idx = 0;
    for (int c = 0; c < NCOLUMNS; c++) {
      for (int r = 0; r < NROWS; r++) {
        assertEquals(A.getQuick(r, c), Avec.getQuick(idx++), TOL);
      }
    }
  }

  void testViewColumn() {
    DComplexMatrix1D col = A.viewColumn(NCOLUMNS ~/ 2);
    expect(NROWS, equals(col.size()));
    for (int r = 0; r < NROWS; r++) {
      assertEquals(A.getQuick(r, NCOLUMNS ~/ 2), col.getQuick(r), TOL);
    }
  }

  void testViewColumnFlip() {
    DComplexMatrix2D B = A.viewColumnFlip();
    expect(A.size(), equals(B.size()));
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        assertEquals(A.getQuick(r, NCOLUMNS - 1 - c), B.getQuick(r, c), TOL);
      }
    }
  }

  void testViewDice() {
    DComplexMatrix2D B = A.viewDice();
    expect(NROWS, equals(B.columns()));
    expect(NCOLUMNS, equals(B.rows()));
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        assertEquals(A.getQuick(r, c), B.getQuick(c, r), TOL);
      }
    }
  }

  void testViewPart() {
    DComplexMatrix2D B = A.viewPart(NROWS ~/ 2, NCOLUMNS ~/ 2, NROWS ~/ 3, NCOLUMNS ~/ 3);
    expect(NROWS ~/ 3, equals(B.rows()));
    expect(NCOLUMNS ~/ 3, equals(B.columns()));
    for (int r = 0; r < NROWS ~/ 3; r++) {
      for (int c = 0; c < NCOLUMNS ~/ 3; c++) {
        assertEquals(A.getQuick(NROWS ~/ 2 + r, NCOLUMNS ~/ 2 + c), B.getQuick(r, c), TOL);
      }
    }
  }

  void testViewRow() {
    DComplexMatrix1D B = A.viewRow(NROWS ~/ 2);
    expect(NCOLUMNS, equals(B.size()));
    for (int r = 0; r < NCOLUMNS; r++) {
      assertEquals(A.getQuick(NROWS ~/ 2, r), B.getQuick(r), TOL);
    }
  }

  void testViewRowFlip() {
    DComplexMatrix2D B = A.viewRowFlip();
    expect(A.size(), equals(B.size()));
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        assertEquals(A.getQuick(NROWS - 1 - r, c), B.getQuick(r, c), TOL);
      }
    }
  }

  void testViewSelectionProc() {
    final Float64List value = new Float64List.fromList([2.0, 3.0]);
    A.assignValue(0.0, 0.0);
    if (DINDEX >= 0) {
      A.setQuick(NROWS ~/ 4, NROWS ~/ 4 + DINDEX, value);
      A.setQuick(NROWS ~/ 2, NROWS ~/ 2 + DINDEX, value);
      DComplexMatrix2D B = A.viewSelectionProc((DComplexMatrix1D element) {
        if (DComplex.abs(DComplex.minus(element.getQuick(NROWS ~/ 4 + DINDEX), value)) < TOL) {
          return true;
        } else {
          return false;
        }
      });
      expect(1, equals(B.rows()));
      expect(NCOLUMNS, equals(B.columns()));
      assertEquals(A.getQuick(NROWS ~/ 4, NROWS ~/ 4 + DINDEX), B.getQuick(0, NROWS ~/ 4 + DINDEX), TOL);
    } else {
      A.setQuick(NROWS ~/ 4 - DINDEX, NROWS ~/ 4, value);
      A.setQuick(NROWS ~/ 2 - DINDEX, NROWS ~/ 2, value);
      DComplexMatrix2D B = A.viewSelectionProc((DComplexMatrix1D element) {
        if (DComplex.abs(DComplex.minus(element.getQuick(NROWS ~/ 4), value)) < TOL) {
          return true;
        } else {
          return false;
        }
      });
      expect(1, equals(B.rows()));
      expect(NCOLUMNS, equals(B.columns()));
      assertEquals(A.getQuick(NROWS ~/ 4 - DINDEX, NROWS ~/ 4), B.getQuick(0, NROWS ~/ 4), TOL);
    }
  }

  void testViewSelection() {
    Int32List rowIndexes = new Int32List.fromList([NROWS ~/ 6, NROWS ~/ 5, NROWS ~/ 4, NROWS ~/ 3, NROWS ~/ 2]);
    Int32List colIndexes = new Int32List.fromList([NROWS ~/ 6, NROWS ~/ 5, NROWS ~/ 4, NROWS ~/ 3, NROWS ~/ 2, NROWS - 1]);
    DComplexMatrix2D B = A.viewSelection(rowIndexes, colIndexes);
    expect(rowIndexes.length, equals(B.rows()));
    expect(colIndexes.length, equals(B.columns()));
    for (int r = 0; r < rowIndexes.length; r++) {
      for (int c = 0; c < colIndexes.length; c++) {
        assertEquals(A.getQuick(rowIndexes[r], colIndexes[c]), B.getQuick(r, c), TOL);
      }
    }
  }

  void testViewStrides() {
    int rowStride = 3;
    int colStride = 5;
    DComplexMatrix2D B = A.viewStrides(rowStride, colStride);
    for (int r = 0; r < B.rows(); r++) {
      for (int c = 0; c < B.columns(); c++) {
        assertEquals(A.getQuick(r * rowStride, c * colStride), B.getQuick(r, c), TOL);
      }
    }
  }

  void testZMult2D() {
    Float64List alpha = new Float64List.fromList([3.0, 4.0]);
    Float64List beta = new Float64List.fromList([5.0, 6.0]);
    DComplexMatrix2D C = new DiagonalDComplexMatrix2D(NROWS, NROWS, 0);
    for (int i = 0; i < DLENGTH; i++) {
      C.setPartsQuick(i, i, random.nextDouble(), random.nextDouble());
    }
    List<Float64List> expected = C.toArray();
    C = A.zMult2D(Bt, C, alpha, beta, false, false);
    Float64List elem = new Float64List(2);
    for (int j = 0; j < NROWS; j++) {
      for (int i = 0; i < NROWS; i++) {
        Float64List s = new Float64List(2);
        for (int k = 0; k < NCOLUMNS; k++) {
          s = DComplex.plus(s, DComplex.multiply(A.getQuick(i, k), Bt.getQuick(k, j)));
        }
        elem[0] = expected[i][2 * j];
        elem[1] = expected[i][2 * j + 1];
        elem = DComplex.multiply(beta, elem);
        s = DComplex.multiply(alpha, s);
        expected[i][2 * j] = s[0] + elem[0];
        expected[i][2 * j + 1] = s[1] + elem[1];
      }
    }
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NROWS; c++) {
        expect(expected[r][2 * c], closeTo(C.getQuick(r, c)[0], TOL));
        expect(expected[r][2 * c + 1], closeTo(C.getQuick(r, c)[1], TOL));
      }
    }

    //---
    C = null;
    C = A.zMult2D(Bt, C, alpha, beta, false, false);
    expected = new List<Float64List>.generate(NROWS,
        (_) => new Float64List(2 * NROWS));
    for (int j = 0; j < NROWS; j++) {
      for (int i = 0; i < NROWS; i++) {
        Float64List s = new Float64List(2);
        for (int k = 0; k < NCOLUMNS; k++) {
          s = DComplex.plus(s, DComplex.multiply(A.getQuick(i, k), Bt.getQuick(k, j)));
        }
        s = DComplex.multiply(alpha, s);
        expected[i][2 * j] = s[0];
        expected[i][2 * j + 1] = s[1];
      }
    }
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NROWS; c++) {
        expect(expected[r][2 * c], closeTo(C.getQuick(r, c)[0], TOL));
        expect(expected[r][2 * c + 1], closeTo(C.getQuick(r, c)[1], TOL));
      }
    }

    //transposeA
    C = new DiagonalDComplexMatrix2D(NCOLUMNS, NCOLUMNS, 0);
    for (int i = 0; i < DLENGTH; i++) {
      C.setPartsQuick(i, i, random.nextDouble(), random.nextDouble());
    }
    expected = C.toArray();
    C = A.zMult2D(B, C, alpha, beta, true, false);
    for (int j = 0; j < NCOLUMNS; j++) {
      for (int i = 0; i < NCOLUMNS; i++) {
        Float64List s = new Float64List(2);
        for (int k = 0; k < NROWS; k++) {
          s = DComplex.plus(s, DComplex.multiply(DComplex.conj(A.getQuick(k, i)), B.getQuick(k, j)));
        }
        elem[0] = expected[i][2 * j];
        elem[1] = expected[i][2 * j + 1];
        elem = DComplex.multiply(beta, elem);
        s = DComplex.multiply(alpha, s);
        expected[i][2 * j] = s[0] + elem[0];
        expected[i][2 * j + 1] = s[1] + elem[1];
      }
    }
    for (int r = 0; r < NCOLUMNS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(expected[r][2 * c], closeTo(C.getQuick(r, c)[0], TOL));
        expect(expected[r][2 * c + 1], closeTo(C.getQuick(r, c)[1], TOL));
      }
    }
    //---
    C = null;
    C = A.zMult2D(B, C, alpha, beta, true, false);
    expected = new List<Float64List>.generate(NCOLUMNS,
        (_) => new Float64List(2 * NCOLUMNS));
    for (int j = 0; j < NCOLUMNS; j++) {
      for (int i = 0; i < NCOLUMNS; i++) {
        Float64List s = new Float64List(2);
        for (int k = 0; k < NROWS; k++) {
          s = DComplex.plus(s, DComplex.multiply(DComplex.conj(A.getQuick(k, i)), B.getQuick(k, j)));
        }
        s = DComplex.multiply(alpha, s);
        expected[i][2 * j] = s[0];
        expected[i][2 * j + 1] = s[1];
      }
    }
    for (int r = 0; r < NCOLUMNS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(expected[r][2 * c], closeTo(C.getQuick(r, c)[0], TOL));
        expect(expected[r][2 * c + 1], closeTo(C.getQuick(r, c)[1], TOL));
      }
    }

    //transposeB
    C = new DiagonalDComplexMatrix2D(NROWS, NROWS, 0);
    for (int i = 0; i < DLENGTH; i++) {
      C.setPartsQuick(i, i, random.nextDouble(), random.nextDouble());
    }
    expected = C.toArray();
    C = A.zMult2D(B, C, alpha, beta, false, true);
    for (int j = 0; j < NROWS; j++) {
      for (int i = 0; i < NROWS; i++) {
        Float64List s = new Float64List(2);
        for (int k = 0; k < NCOLUMNS; k++) {
          s = DComplex.plus(s, DComplex.multiply(A.getQuick(i, k), DComplex.conj(B.getQuick(j, k))));
        }
        elem[0] = expected[i][2 * j];
        elem[1] = expected[i][2 * j + 1];
        elem = DComplex.multiply(beta, elem);
        s = DComplex.multiply(alpha, s);
        expected[i][2 * j] = s[0] + elem[0];
        expected[i][2 * j + 1] = s[1] + elem[1];
      }
    }
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NROWS; c++) {
        expect(expected[r][2 * c], closeTo(C.getQuick(r, c)[0], TOL));
        expect(expected[r][2 * c + 1], closeTo(C.getQuick(r, c)[1], TOL));
      }
    }
    //---
    C = null;
    C = A.zMult2D(B, C, alpha, beta, false, true);
    expected = new List<Float64List>.generate(NROWS,
        (_) => new Float64List(2 * NROWS));
    for (int j = 0; j < NROWS; j++) {
      for (int i = 0; i < NROWS; i++) {
        Float64List s = new Float64List(2);
        for (int k = 0; k < NCOLUMNS; k++) {
          s = DComplex.plus(s, DComplex.multiply(A.getQuick(i, k), DComplex.conj(B.getQuick(j, k))));
        }
        s = DComplex.multiply(alpha, s);
        expected[i][2 * j] = s[0];
        expected[i][2 * j + 1] = s[1];
      }
    }
    for (int r = 0; r < NROWS; r++) {
      for (int c = 0; c < NROWS; c++) {
        expect(expected[r][2 * c], closeTo(C.getQuick(r, c)[0], TOL));
        expect(expected[r][2 * c + 1], closeTo(C.getQuick(r, c)[1], TOL));
      }
    }
    //transposeA and transposeB
    C = new DiagonalDComplexMatrix2D(NCOLUMNS, NCOLUMNS, 0);
    for (int i = 0; i < DLENGTH; i++) {
      C.setPartsQuick(i, i, random.nextDouble(), random.nextDouble());
    }
    expected = C.toArray();
    C = A.zMult2D(Bt, C, alpha, beta, true, true);
    for (int j = 0; j < NCOLUMNS; j++) {
      for (int i = 0; i < NCOLUMNS; i++) {
        Float64List s = new Float64List(2);
        for (int k = 0; k < NROWS; k++) {
          s = DComplex.plus(s, DComplex.multiply(DComplex.conj(A.getQuick(k, i)), DComplex.conj(Bt.getQuick(j, k))));
        }
        elem[0] = expected[i][2 * j];
        elem[1] = expected[i][2 * j + 1];
        elem = DComplex.multiply(beta, elem);
        s = DComplex.multiply(alpha, s);
        expected[i][2 * j] = s[0] + elem[0];
        expected[i][2 * j + 1] = s[1] + elem[1];
      }
    }
    for (int r = 0; r < NCOLUMNS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(expected[r][2 * c], closeTo(C.getQuick(r, c)[0], TOL));
        expect(expected[r][2 * c + 1], closeTo(C.getQuick(r, c)[1], TOL));
      }
    }
    //---
    C = null;
    C = A.zMult2D(Bt, C, alpha, beta, true, true);
    expected = new List<Float64List>.generate(NCOLUMNS,
        (_) => new Float64List(2 * NCOLUMNS));
    for (int j = 0; j < NCOLUMNS; j++) {
      for (int i = 0; i < NCOLUMNS; i++) {
        Float64List s = new Float64List(2);
        for (int k = 0; k < NROWS; k++) {
          s = DComplex.plus(s, DComplex.multiply(A.getQuick(k, i), Bt.getQuick(j, k)));
        }
        s = DComplex.multiply(alpha, s);
        expected[i][2 * j] = s[0];
        expected[i][2 * j + 1] = s[1];
      }
    }
    for (int r = 0; r < NCOLUMNS; r++) {
      for (int c = 0; c < NCOLUMNS; c++) {
        expect(expected[r][2 * c], closeTo(C.getQuick(r, c)[0], TOL));
        expect(expected[r][2 * c + 1], closeTo(C.getQuick(r, c)[1], TOL));
      }
    }

  }
}

class DiagonalDComplexMatrix2DViewTest extends DiagonalDComplexMatrix2DTest {
  void createMatrices() {
    DINDEX = 3;
    A = new DiagonalDComplexMatrix2D(NCOLUMNS, NROWS, -DINDEX);
    DLENGTH = (A as DiagonalDComplexMatrix2D).diagonalLength();
    A = A.viewDice();
    B = new DiagonalDComplexMatrix2D(NCOLUMNS, NROWS, -DINDEX).viewDice();
    Bt = new DiagonalDComplexMatrix2D(NROWS, NCOLUMNS, DINDEX).viewDice();
  }
}
