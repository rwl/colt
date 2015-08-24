library np;

import 'dart:typed_data';
import 'dart:math' as math;

import 'package:complex/complex.dart';
import 'package:quiver/iterables.dart' as iter;

import 'colt.dart';
import 'function/int.dart' as ifunc;
import 'function/double.dart' as dfunc;
import 'function/complex.dart' as cfunc;

const pi = math.PI;

DoubleVector array(Iterable<num> a) {
  var l = new Float64List.fromList(a.map((n) {
    return n.toDouble();
  }).toList(growable: false));
  return new DenseDoubleVector.fromList(l);
}

IntVector iarray(Iterable<num> a) {
  var l = new Int32List.fromList(a.map((n) {
    return n.toInt();
  }).toList(growable: false));
  return new DenseIntVector.fromList(l);
}

ComplexVector complex(List<double> re, [List<double> im]) {
  if (re != null && im != null) {
    var _re = new Float64List.fromList(re);
    var _im = new Float64List.fromList(im);
    _re = new DenseDoubleVector.fromList(_re);
    _im = new DenseDoubleVector.fromList(_im);
    return new DenseComplexVector.fromParts(_re, _im);
  } else if (im != null) {
    var _im = new Float64List.fromList(im);
    _im = new DenseDoubleVector.fromList(_im);
    return new DenseComplexVector.fromImaginary(_im);
  } else if (re != null) {
    var _re = new Float64List.fromList(re);
    _re = new DenseDoubleVector.fromList(_re);
    return new DenseComplexVector.fromReal(_re);
  }
  return new ComplexVector(0);
}

DoubleVector div(DoubleVector a, num b) {
  return a.copy()..apply(dfunc.divide(b.toDouble()));
}

ComplexVector cdiv(ComplexVector a, num b) {
  return a.copy()..apply(cfunc.divideBy(b.toDouble()));
}

DoubleVector mult(DoubleVector a, num b) {
  return a.copy()..apply(dfunc.multiply(b.toDouble()));
}

ComplexVector cmult(ComplexVector a, num re, [num im = 0]) {
  return a.copy()..apply(cfunc.multiply(new Complex(re, im)));
}

DoubleVector abs(ComplexVector a) => a.abs().real();

ComplexVector conj(ComplexVector a) => a.conj();

DoubleVector angle(ComplexVector a) => a.arg().real();

DoubleVector sin(DoubleVector a) {
  return a.copy()..apply(dfunc.sin);
}

DoubleVector cos(DoubleVector a) {
  return a.copy()..apply(dfunc.cos);
}

DoubleVector radians(Iterable<num> a) {
  return array(a).copy()..apply(dfunc.toRadians);
}

DoubleVector degrees(Iterable<num> a) {
  return array(a).copy()..apply(dfunc.toDegrees);
}

ComplexVector cinv(ComplexVector a) {
  return a.copy()..apply(cfunc.inv);
}

DoubleVector pow(DoubleVector a, num b) {
  return a.copy()..apply(dfunc.power(b.toDouble()));
}

double round(double a, num precision) {
  return double.parse(a.toStringAsFixed(precision));
}

double sum(List<num> a) => a.reduce((a, b) => a + b).toDouble();

double prod(List<num> a) => a.reduce((a, b) => a * b).toDouble();

DoubleVector real(ComplexVector a) => a.real();

DoubleVector imag(ComplexVector a) => a.imaginary();

DoubleVector fill(size, num a) {
  return new DoubleVector(size)..fill(a.toDouble());
}

DoubleVector zeros(int size) => new DoubleVector(size);

ComplexVector czeros(int size) => new ComplexVector(size);

IntVector ieq(IntVector a, int b) {
  return a.copy()..apply(ifunc.equalTo(b));
}

DoubleVector eq(DoubleVector a, num b) {
  return a.copy()..apply(dfunc.equalTo(b.toDouble()));
}

IntVector ineq(IntVector a, IntVector b) {
  return a.copy()..assign(b, ifunc.equals);
}

DoubleVector ones(int size) => fill(size, 1);

ComplexVector cones(int size) {
  return new ComplexVector(size)..fill(1.0, 0.0);
}

Int32List nonzero(List<num> a) {
  var ix = <int>[];
  array(a).nonzero(indexList: ix);
  return new Int32List.fromList(ix);
}

IntVector diff(IntVector a) {
  if (a.size == 0) {
    return a;
  }
  var b = new IntVector(a.size - 1);
  for (var i = 0; i < a.size - 1; i++) {
    b[i] = a[i + 1] - a[i];
  }
  return b;
}

Int32List range(int start_or_stop, [int stop, int step]) {
  var r = iter.range(start_or_stop, stop, step).map((r) {
    return r.toInt();
  }).toList(growable: false);
  return new Int32List.fromList(r);
}

ComplexVector polar(List<num> r, List<num> theta, [bool radians = true]) {
  return new DenseComplexVector.fromPolar(array(r), array(theta), radians);
}

DoubleVector concat(DoubleVector a, DoubleVector b) {
  return new DenseDoubleVector.append(a, b);
}

Int32List cat(List<int> a, List<int> b) {
  return new Int32List.fromList(iter.concat([a, b]).toList(growable: false));
}

ComplexMatrix csparse(int rows, int columns, List<int> rowIndexes,
    List<int> columnIndexes, List<Complex> values,
    {bool removeDuplicates: false, bool removeZeroes: false}) {
  var ridx = new Int32List.fromList(rowIndexes);
  var cidx = new Int32List.fromList(columnIndexes);
  var vals = new Float64List(values.length * 2);
  for (int i = 0; i < values.length; i++) {
    vals[2 * i] = values[i].real;
    vals[2 * i + 1] = values[i].imaginary;
  }
  return new SparseRCComplexMatrix.withValues(rows, columns, ridx, cidx, vals,
      removeDuplicates: removeDuplicates, removeZeroes: removeZeroes);
}

ComplexMatrix cspdiag(ComplexVector a) {
  int n = a.size;
  var ix = range(n);
  return csparse(n, n, ix, ix, a, removeDuplicates: true);
}

SparseComplexVector csparray(int size, Int32List indexes, Float64List values) {
  var v = new SparseComplexVector(size);
  v.select(indexes).setValues(values);
  return v;
}

DoubleMatrix sparse(int rows, int columns, List<int> rowIndexes,
    List<int> columnIndexes, List<double> values,
    {bool removeDuplicates: false, bool removeZeroes: false}) {
  var ridx = new Int32List.fromList(rowIndexes);
  var cidx = new Int32List.fromList(columnIndexes);
  var vals = new Float64List.fromList(values);
  return new SparseRCDoubleMatrix.withValues(rows, columns, ridx, cidx, vals,
      removeDuplicates: removeDuplicates, removeZeroes: removeZeroes);
}

DoubleMatrix spdiag(DoubleVector a) {
  int n = a.size;
  var ix = range(n);
  return sparse(n, n, ix, ix, a.toList(), removeDuplicates: true);
}
