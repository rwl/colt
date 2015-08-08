part of cern.colt.matrix.complex.test;

const SIZE = 2 * 17 * 5;

const TOL = 1e-10;

testAbstractComplexVector(String kind, AbstractComplexVector make(int size)) {
  group("AbstractComplexVector ($kind)", () {
    AbstractComplexVector A, B;
    setUp(() {
      A = make(SIZE);
      B = make(SIZE);
      for (int i = 0; i < A.size; i++) {
        A.setParts(i, random.nextDouble(), random.nextDouble());
        B.setParts(i, random.nextDouble(), random.nextDouble());
      }
    });
    test('aggregate', () {
      var expected = Complex.ZERO;
      for (int i = 0; i < A.size; i++) {
        expected = expected + A.get(i).pow(2);
      }
      var result = A.aggregate(plus, square);
      assertEquals(expected, result, TOL);
    });

    test('apply', () {
      AbstractComplexVector Acopy = A.copy();
      A.apply(acos);
      for (int i = 0; i < A.size; i++) {
        var expected = Acopy.get(i).acos();
        assertEquals(expected, A.get(i), TOL);
      }
    });

    test('copyFrom', () {
      A.copyFrom(B);
      expect(A.size == B.size, isTrue);
      for (int i = 0; i < A.size; i++) {
        assertEquals(B.get(i), A.get(i), TOL);
      }
    });

    test('assign', () {
      AbstractComplexVector Acopy = A.copy();
      A.assign(B, div);
      for (int i = 0; i < A.size; i++) {
        assertEquals(Acopy.get(i) / B.get(i), A.get(i), TOL);
      }
    });

    test('applyReal', () {
      AbstractComplexVector Acopy = A.copy();
      A.applyReal(abs);
      for (int i = 0; i < A.size; i++) {
        var elem = A.get(i);
        expect(Acopy.get(i).abs(), closeTo(elem.real, TOL));
        expect(0, closeTo(elem.imaginary, TOL));
      }
    });

    test('setAll', () {
      var expected = new Float64List(2 * A.size);
      for (int i = 0; i < 2 * A.size; i++) {
        expected[i] = random.nextDouble();
      }
      A.setAll(expected);
      for (int i = 0; i < A.size; i++) {
        var elem = A.get(i);
        expect(expected[2 * i], closeTo(elem.real, TOL));
        expect(expected[2 * i + 1], closeTo(elem.imaginary, TOL));
      }
    });

    test('fill', () {
      double re = random.nextDouble();
      double im = random.nextDouble();
      A.fill(re, im);
      for (int i = 0; i < A.size; i++) {
        var elem = A.get(i);
        expect(re, closeTo(elem.real, TOL));
        expect(im, closeTo(elem.imaginary, TOL));
      }
    });

    test('setImaginary', () {
      AbstractComplexVector Acopy = A.copy();
      AbstractDoubleVector Im = new DoubleVector.random(A.size);
      A.setImaginary(Im);
      for (int i = 0; i < A.size; i++) {
        var elem = A.get(i);
        expect(Acopy.get(i).real, closeTo(elem.real, TOL));
        expect(Im.get(i), closeTo(elem.imaginary, TOL));
      }
    });

    test('setReal', () {
      AbstractComplexVector Acopy = A.copy();
      AbstractDoubleVector Re = new DoubleVector.random(A.size);
      A.setReal(Re);
      for (int i = 0; i < A.size; i++) {
        var elem = A.get(i);
        expect(Acopy.get(i).imaginary, closeTo(elem.imaginary, TOL));
        expect(Re.get(i), closeTo(elem.real, TOL));
      }
    });

    test('cardinality', () {
      int card = A.cardinality;
      expect(A.size, equals(card));
    });

    test('equals', () {
      expect(A.equals(A), isTrue);
      expect(A.equals(B), isFalse);
    });

    test('imaginary', () {
      AbstractDoubleVector Im = A.imaginary();
      for (int i = 0; i < A.size; i++) {
        expect(A.get(i).imaginary, closeTo(Im.get(i), TOL));
      }
    });

    test('real', () {
      AbstractDoubleVector Re = A.real();
      for (int i = 0; i < A.size; i++) {
        expect(A.get(i).real, closeTo(Re.get(i), TOL));
      }
    });

    test('nonzero', () {
      var indexList = <int>[];
      var valueList = <Complex>[];
      A.nonzero(indexList: indexList, valueList: valueList);
      expect(A.size, equals(indexList.length));
      expect(A.size, equals(valueList.length));
      for (int i = 0; i < A.size; i++) {
        assertEquals(A.get(indexList[i]), valueList[i], TOL);
        expect(valueList[i].real != 0 || valueList[i].imaginary != 0, isTrue);
      }
    });

    test('reshape', () {
      int rows = 10;
      int columns = 17;
      AbstractComplexMatrix B = A.reshape(rows, columns);
      int idx = 0;
      for (int c = 0; c < columns; c++) {
        for (int r = 0; r < rows; r++) {
          assertEquals(A.get(idx++), B.get(r, c), TOL);
        }
      }
    });

    test('toList', () {
      Float64List array = A.toList();
      for (int i = 0; i < A.size; i++) {
        var elem = A.get(i);
        expect(elem.real, closeTo(array[2 * i], TOL));
        expect(elem.imaginary, closeTo(array[2 * i + 1], TOL));
      }
    });

    test('fillList', () {
      var array = new Float64List(2 * A.size);
      A.fillList(array);
      for (int i = 0; i < A.size; i++) {
        var elem = A.get(i);
        expect(elem.real, closeTo(array[2 * i], TOL));
        expect(elem.imaginary, closeTo(array[2 * i + 1], TOL));
      }
    });

    test('flip', () {
      AbstractComplexVector B = A.flip();
      for (int i = 0; i < A.size; i++) {
        assertEquals(A.get(A.size - 1 - i), B.get(i), TOL);
      }
    });

    test('part', () {
      AbstractComplexVector B = A.part(A.size ~/ 2, A.size ~/ 3);
      for (int i = 0; i < A.size / 3; i++) {
        assertEquals(A.get(A.size ~/ 2 + i), B.get(i), TOL);
      }
    });

    test('select', () {
      var indexes = new Int32List.fromList(
          [A.size ~/ 6, A.size ~/ 5, A.size ~/ 4, A.size ~/ 3, A.size ~/ 2]);
      AbstractComplexVector B = A.select(indexes);
      for (int i = 0; i < indexes.length; i++) {
        assertEquals(A.get(indexes[i]), B.get(i), TOL);
      }
    });

    test('strides', () {
      int stride = 3;
      AbstractComplexVector B = A.strides(stride);
      for (int i = 0; i < B.size; i++) {
        assertEquals(A.get(i * stride), B.get(i), TOL);
      }
    });

    test('dot', () {
      var actual = A.dot(B);
      var expected = Complex.ZERO;
      for (int i = 0; i < A.size; i++) {
        expected += B.get(i).conjugate() * A.get(i);
      }
      assertEquals(expected, actual, TOL);
    });

    test('dot range', () {
      var actual = A.dot(B, 5, B.size - 10);
      var expected = Complex.ZERO;
      for (int i = 5; i < A.size - 5; i++) {
        expected += B.get(i).conjugate() * A.get(i);
      }
      assertEquals(expected, actual, TOL);
    });

    test('sum', () {
      var actual = A.sum();
      var expected = Complex.ZERO;
      for (int i = 0; i < A.size; i++) {
        expected += A.get(i);
      }
      assertEquals(expected, actual, TOL);
    });
  });
}

void assertEquals(Complex expected, Complex actual, double tol) {
  //for (int i = 0; i < actual.length; i++) {
    expect(cmath.isEqual(expected, actual, tol), isTrue);
  //}
}
