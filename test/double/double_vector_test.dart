part of cern.colt.matrix.double.test;

const double TOL = 1e-10;

const int SIZE = 2 * 17 * 5;

testDoubleVector(String kind, DoubleVector make(int size)) {
  group("DoubleVector ($kind)", () {
    DoubleVector A, B;

    setUp(() {
      A = make(SIZE);
      B = make(SIZE);
      for (int i = 0; i < SIZE; i++) {
        A.set(i, _r.nextDouble());
        B.set(i, _r.nextDouble());
      }
    });

    test('aggregate', () {
      double expected = 0.0;
      for (int i = 0; i < A.size; i++) {
        double elem = A.get(i);
        expected += elem * elem;
      }
      double result = A.aggregate(plus, square);
      expect(result, closeTo(expected, TOL));
    });

    test('fill', () {
      double value = _r.nextDouble();
      A.fill(value);
      for (int i = 0; i < A.size; i++) {
        expect(value, closeTo(A.get(i), TOL));
      }
    });

    test('setAll', () {
      Float64List expected = new Float64List(A.size);
      for (int i = 0; i < A.size; i++) {
        expected[i] = _r.nextDouble();
      }
      A.setAll(expected);
      for (int i = 0; i < A.size; i++) {
        expect(expected[i], closeTo(A.get(i), TOL));
      }
    });

    test('apply', () {
      DoubleVector Acopy = A.copy();
      A.apply(acos);
      for (int i = 0; i < A.size; i++) {
        double expected = math.acos(Acopy.get(i));
        expect(expected, closeTo(A.get(i), TOL));
      }
    });

    test('copyFrom', () {
      A.copyFrom(B);
      expect(A.size, equals(B.size));
      for (int i = 0; i < A.size; i++) {
        expect(B.get(i), closeTo(A.get(i), TOL));
      }
    });

    test('assign', () {
      DoubleVector Acopy = A.copy();
      A.assign(B, div);
      for (int i = 0; i < A.size; i++) {
        expect(Acopy.get(i) / B.get(i), closeTo(A.get(i), TOL));
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

    test('max', () {
      A.fill(0.0);
      A.set(A.size ~/ 3, 0.7);
      A.set(A.size ~/ 2, 0.1);
      var maxAndLoc = A.max();
      expect(0.7, closeTo(maxAndLoc.value, TOL));
      expect(A.size ~/ 3, equals(maxAndLoc.location));
    });

    test('min', () {
      A.fill(0.0);
      A.set(A.size ~/ 3, -0.7);
      A.set(A.size ~/ 2, -0.1);
      var minAndLoc = A.min();
      expect(-0.7, closeTo(minAndLoc.value, TOL));
      expect(A.size ~/ 3, equals(minAndLoc.location));
    });

    test('negative', () {
      A.fill(0.0);
      A.set(A.size ~/ 3, -0.7);
      A.set(A.size ~/ 2, -0.1);
      List<int> indexList = new List<int>();
      List<double> valueList = new List<double>();
      A.negative(indexList: indexList, valueList: valueList);
      expect(2, indexList.length);
      expect(2, valueList.length);
      expect(indexList.contains(A.size ~/ 3), isTrue);
      expect(indexList.contains(A.size ~/ 2), isTrue);
      expect(valueList.contains(-0.7), isTrue);
      expect(valueList.contains(-0.1), isTrue);
    });

    test('nonzero', () {
      A.fill(0.0);
      A.set(A.size ~/ 3, 0.7);
      A.set(A.size ~/ 2, 0.1);
      List<int> indexList = new List<int>();
      List<double> valueList = new List<double>();
      A.nonzero(indexList: indexList, valueList: valueList);
      expect(2, indexList.length);
      expect(2, valueList.length);
      expect(indexList.contains(A.size ~/ 3), isTrue);
      expect(indexList.contains(A.size ~/ 2), isTrue);
      expect(valueList.contains(0.7), isTrue);
      expect(valueList.contains(0.1), isTrue);
    });

    test('positive', () {
      A.fill(0.0);
      A.set(A.size ~/ 3, 0.7);
      A.set(A.size ~/ 2, 0.1);
      List<int> indexList = new List<int>();
      List<double> valueList = new List<double>();
      A.positive(indexList: indexList, valueList: valueList);
      expect(2, indexList.length);
      expect(2, valueList.length);
      expect(indexList.contains(A.size ~/ 3), isTrue);
      expect(indexList.contains(A.size ~/ 2), isTrue);
      expect(valueList.contains(0.7), isTrue);
      expect(valueList.contains(0.1), isTrue);
    });

    test('toList', () {
      Float64List array = A.toList();
      expect(A.size == array.length, isTrue);
      for (int i = 0; i < A.size; i++) {
        expect(array[i], closeTo(A.get(i), TOL));
      }
    });

    test('fillList', () {
      var array = new Float64List(A.size);
      A.fillList(array);
      for (int i = 0; i < A.size; i++) {
        expect(A.get(i), closeTo(array[i], TOL));
      }
    });

    test('reshape', () {
      int rows = 10;
      int columns = 17;
      DoubleMatrix B = A.reshape(rows, columns);
      int idx = 0;
      for (int c = 0; c < columns; c++) {
        for (int r = 0; r < rows; r++) {
          expect(A.get(idx++), closeTo(B.get(r, c), TOL));
        }
      }
    });

    test('flip', () {
      DoubleVector b = A.flip();
      expect(A.size, equals(b.size));
      for (int i = 0; i < A.size; i++) {
        expect(A.get(i), closeTo(b.get(A.size - 1 - i), TOL));
      }
    });

    test('part', () {
      DoubleVector b = A.part(15, 11);
      for (int i = 0; i < 11; i++) {
        expect(A.get(15 + i), closeTo(b.get(i), TOL));
      }
    });

    test('select', () {
      var indexes = new Int32List.fromList([5, 11, 22, 37, 101]);
      DoubleVector b = A.select(indexes);
      for (int i = 0; i < indexes.length; i++) {
        expect(A.get(indexes[i]), closeTo(b.get(i), TOL));
      }
    });

    test('strides', () {
      int stride = 3;
      DoubleVector b = A.strides(stride);
      for (int i = 0; i < b.size; i++) {
        expect(A.get(i * stride), closeTo(b.get(i), TOL));
      }
    });

    test('dot', () {
      double product = A.dot(B);
      double expected = 0.0;
      for (int i = 0; i < A.size; i++) {
        expected += A.get(i) * B.get(i);
      }
      expect(expected, closeTo(product, TOL));
    });

    test('dot range', () {
      double product = A.dot(B, 5, B.size - 10);
      double expected = 0.0;
      for (int i = 5; i < A.size - 5; i++) {
        expected += A.get(i) * B.get(i);
      }
      expect(expected, closeTo(product, TOL));
    });

    test('sum', () {
      double sum = A.sum();
      double expected = 0.0;
      for (int i = 0; i < A.size; i++) {
        expected += A.get(i);
      }
      expect(expected, closeTo(sum, TOL));
    });
  });
}
