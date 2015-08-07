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
      Float64List expected = new Float64List(2);
      for (int i = 0; i < A.size; i++) {
        expected = cmath.plus(expected, cmath.square(A.get(i)));
      }
      Float64List result = A.aggregate(plus, square);
      assertEquals(expected, result, TOL);
    });

    test('apply', () {
      AbstractComplexVector Acopy = A.copy();
      A.apply(acos);
      for (int i = 0; i < A.size; i++) {
        Float64List expected = cmath.acos(Acopy.get(i));
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
        assertEquals(cmath.div_(Acopy.get(i), B.get(i)), A.get(i), TOL);
      }
    });

    test('applyReal', () {
      AbstractComplexVector Acopy = A.copy();
      A.applyReal(abs);
      for (int i = 0; i < A.size; i++) {
        Float64List elem = A.get(i);
        expect(cmath.abs(Acopy.get(i)), closeTo(elem[0], TOL));
        expect(0, closeTo(elem[1], TOL));
      }
    });

    test('setAll', () {
      Float64List expected = new Float64List(2 * A.size);
      for (int i = 0; i < 2 * A.size; i++) {
        expected[i] = random.nextDouble();
      }
      A.setAll(expected);
      for (int i = 0; i < A.size; i++) {
        Float64List elem = A.get(i);
        expect(expected[2 * i], closeTo(elem[0], TOL));
        expect(expected[2 * i + 1], closeTo(elem[1], TOL));
      }
    });

    test('fill', () {
      double re = random.nextDouble();
      double im = random.nextDouble();
      A.fill(re, im);
      for (int i = 0; i < A.size; i++) {
        Float64List elem = A.get(i);
        expect(re, closeTo(elem[0], TOL));
        expect(im, closeTo(elem[1], TOL));
      }
    });

    test('setImaginary', () {
      AbstractComplexVector Acopy = A.copy();
      AbstractDoubleVector Im = new DoubleVector.random(A.size);
      A.setImaginary(Im);
      for (int i = 0; i < A.size; i++) {
        Float64List elem = A.get(i);
        expect(Acopy.get(i)[0], closeTo(elem[0], TOL));
        expect(Im.get(i), closeTo(elem[1], TOL));
      }
    });

    test('setReal', () {
      AbstractComplexVector Acopy = A.copy();
      AbstractDoubleVector Re = new DoubleVector.random(A.size);
      A.setReal(Re);
      for (int i = 0; i < A.size; i++) {
        Float64List elem = A.get(i);
        expect(Acopy.get(i)[1], closeTo(elem[1], TOL));
        expect(Re.get(i), closeTo(elem[0], TOL));
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
        expect(A.get(i)[1], closeTo(Im.get(i), TOL));
      }
    });

    test('real', () {
      AbstractDoubleVector Re = A.real();
      for (int i = 0; i < A.size; i++) {
        expect(A.get(i)[0], closeTo(Re.get(i), TOL));
      }
    });

    test('nonzero', () {
      List<int> indexList = new List<int>();
      List<Float64List> valueList = new List<Float64List>();
      A.nonzero(indexList: indexList, valueList: valueList);
      expect(A.size, equals(indexList.length));
      expect(A.size, equals(valueList.length));
      for (int i = 0; i < A.size; i++) {
        assertEquals(A.get(indexList[i]), valueList[i], TOL);
        expect(valueList[i][0] != 0 || valueList[i][1] != 0, isTrue);
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
        Float64List elem = A.get(i);
        expect(elem[0], closeTo(array[2 * i], TOL));
        expect(elem[1], closeTo(array[2 * i + 1], TOL));
      }
    });

    test('fillList', () {
      Float64List array = new Float64List(2 * A.size);
      A.fillList(array);
      for (int i = 0; i < A.size; i++) {
        Float64List elem = A.get(i);
        expect(elem[0], closeTo(array[2 * i], TOL));
        expect(elem[1], closeTo(array[2 * i + 1], TOL));
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
      Int32List indexes = new Int32List.fromList(
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
      Float64List actual = A.dot(B);
      Float64List expected = new Float64List(2);
      for (int i = 0; i < A.size; i++) {
        expected = cmath.plus(
            expected, cmath.multiply(cmath.conj(B.get(i)), A.get(i)));
      }
      assertEquals(expected, actual, TOL);
    });

    test('dot range', () {
      Float64List actual = A.dot(B, 5, B.size - 10);
      Float64List expected = new Float64List(2);
      for (int i = 5; i < A.size - 5; i++) {
        expected = cmath.plus(
            expected, cmath.multiply(cmath.conj(B.get(i)), A.get(i)));
      }
      assertEquals(expected, actual, TOL);
    });

    test('sum', () {
      Float64List actual = A.sum();
      Float64List expected = new Float64List(2);
      for (int i = 0; i < A.size; i++) {
        expected = cmath.plus(expected, A.get(i));
      }
      assertEquals(expected, actual, TOL);
    });
  });
}

void assertEquals(List<double> expected, List<double> actual, double tol) {
  for (int i = 0; i < actual.length; i++) {
    expect(expected[i], closeTo(actual[i], tol));
  }
}
