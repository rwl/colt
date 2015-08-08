library cern.colt.matrix.function.complex;

import 'dart:math' as Math;
import 'package:complex/complex.dart';

typedef Complex ComplexComplexComplexFunction(Complex x, Complex y);
typedef Complex ComplexComplexFunction(Complex x);
typedef double ComplexComplexRealFunction(Complex x, Complex y);
typedef double ComplexComplexRealRealFunction(Complex x, Complex y, double tol);
typedef Complex ComplexRealComplexFunction(Complex x, double y);
typedef double ComplexRealFunction(Complex x);
typedef Complex IntIntComplexFunction(int x, int y, Complex z);
typedef Complex RealComplexComplexFunction(double x, Complex y);
typedef Complex RealComplexFunction(double x);

/// Unary functions

double abs(Complex x) => x.abs();

Complex acos(Complex x) => x.acos();

double arg(Complex x) => x.argument();

Complex asin(Complex x) => x.asin();

Complex atan(Complex x) => x.atan();

Complex conj(Complex x) => x.conjugate();

Complex cos(Complex x) => x.cos();

Complex exp(Complex x) => x.exp();

Complex identity(Complex x) => x;

Complex inv(Complex x) => x.reciprocal();

Complex log(Complex x) => x.log();

Complex neg(Complex x) => new Complex(-x.real, -x.imaginary);

Complex sin(Complex x) => x.sin();

Complex sqrt(Complex x) => x.sqrt();

Complex square(Complex x) => x.pow(2);

Complex tan(Complex x) => x.tan();

/// Binary functions

Complex div(Complex x, Complex y) => x / y;

double equals(Complex x, Complex y, double tol) {
  var diff = new Complex(x.real - y.real, x.imaginary - y.imaginary);
  if (diff.abs() <= tol.abs()) {
    return 1.0;
  } else {
    return 0.0;
  }
}

Complex minus(Complex x, Complex y) => x - y;

Complex mult(Complex x, Complex y) => x * y;

Complex multConjFirst(Complex x, Complex y) => x.conjugate() * y;

Complex multConjSecond(Complex x, Complex y) => x * y.conjugate();

Complex plus(Complex x, Complex y) => x + y;

ComplexComplexFunction bindArg(ComplexComplexComplexFunction fn, Complex c) {
  return (Complex v) => fn(c, v);
}

ComplexComplexFunction constant(Complex c) => (Complex x) => c;

ComplexComplexFunction divide(Complex b) {
  return multiply(inv(b));
}

ComplexComplexFunction divideBy(double b) {
  var tmp = new Complex(b);
  return multiply(inv(tmp));
}

ComplexRealFunction equalTo(Complex y) {
  return (Complex x) => (x == y ? 1.0 : 0.0);
}

ComplexComplexFunction subtract(Complex y) {
  return (Complex x) => x - y;
}

ComplexComplexComplexFunction minusMult(Complex constant) {
  return plusMultSecond(neg(constant));
}

ComplexComplexFunction multiply(Complex y) {
  return new ComplexMult(y);
}

ComplexComplexFunction scale(double y) {
  return new ComplexMult(new Complex(y));
}

ComplexComplexFunction add(Complex y) {
  return (Complex x) => x + y;
}

ComplexComplexComplexFunction plusMultSecond(Complex constant) {
  return new ComplexPlusMultSecond(constant);
}

ComplexComplexComplexFunction plusMultFirst(Complex constant) {
  return new ComplexPlusMultFirst(constant);
}

ComplexComplexFunction pow(double y) => (Complex x) => x.pow(y);

ComplexComplexFunction power(Complex y) => (Complex x) => x.power(y);

Math.Random _r;

Complex random(Complex _) {
  if (_r == null) {
    _r = new Math.Random();
  }
  return new Complex(_r.nextDouble(), _r.nextDouble());
}

ComplexComplexComplexFunction swap(ComplexComplexComplexFunction fn) {
  return (Complex x, Complex y) => fn(y, x);
}

/// Only for performance tuning of compute intensive linear algebraic
/// computations. Constructs functions that return one of:
/// - `a * constant`
/// - `a / constant`
/// `a` is variable, `constant` is fixed, but for performance
/// reasons publicly accessible. Intended to be passed to
/// `matrix.apply(function)` methods.
class ComplexMult {
  /// Public read/write access to avoid frequent object construction.
  Complex multiplicator;

  ComplexMult(this.multiplicator);

  /// Returns the result of the function evaluation.
  Complex call(Complex a) => a * multiplicator;

  /// `a / constant`.
  static ComplexMult div(Complex constant) {
    return mult(inv(constant));
  }

  /// `a * constant`.
  static ComplexMult mult(Complex constant) {
    return new ComplexMult(constant);
  }
}

/// Only for performance tuning of compute intensive linear algebraic
/// computations. Constructs functions that return one of:
/// - `a*constant + b`
/// - `a*constant - b`
/// - `a/constant + b`
/// - `a/constant - b`
/// `a` and `b` are variables, `constant` is fixed, but for
/// performance reasons publicly accessible. Intended to be passed to
/// `matrix.assign(otherMatrix,function)` methods.
class ComplexPlusMultFirst {
  /// Public read/write access to avoid frequent object construction.
  Complex multiplicator;

  ComplexPlusMultFirst(this.multiplicator);

  /// Returns the result of the function evaluation.
  Complex call(Complex a, Complex b) => (a * multiplicator) + b;

  /// `a - b/constant`.
  static ComplexPlusMultFirst minusDiv(Complex constant) {
    return new ComplexPlusMultFirst(neg(inv(constant)));
  }

  /// `a - b*constant`.
  static ComplexPlusMultFirst minusMult(Complex constant) {
    return new ComplexPlusMultFirst(neg(constant));
  }

  /// `a + b/constant`.
  static ComplexPlusMultFirst plusDiv(Complex constant) {
    return new ComplexPlusMultFirst(inv(constant));
  }

  /// `a + b*constant`.
  static ComplexPlusMultFirst plusMult(Complex constant) {
    return new ComplexPlusMultFirst(constant);
  }
}

/// Only for performance tuning of compute intensive linear algebraic
/// computations. Constructs functions that return one of:
/// - `a + b*constant`
/// - `a - b*constant`
/// - `a + b/constant`
/// - `a - b/constant`
/// `a` and `b` are variables, `constant` is fixed, but for
/// performance reasons publicly accessible. Intended to be passed to
/// `matrix.assign(otherMatrix,function)` methods.
class ComplexPlusMultSecond {
  /// Public read/write access to avoid frequent object construction.
  Complex multiplicator;

  ComplexPlusMultSecond(this.multiplicator);

  /// Returns the result of the function evaluation.
  Complex call(Complex a, Complex b) => a + (b * multiplicator);

  /// `a - b/constant`.
  static ComplexPlusMultSecond minusDiv(Complex constant) {
    return new ComplexPlusMultSecond(neg(inv(constant)));
  }

  /// `a - b*constant`.
  static ComplexPlusMultSecond minusMult(Complex constant) {
    return new ComplexPlusMultSecond(neg(constant));
  }

  /// `a + b/constant`.
  static ComplexPlusMultSecond plusDiv(Complex constant) {
    return new ComplexPlusMultSecond(inv(constant));
  }

  /// `a + b*constant`.
  static ComplexPlusMultSecond plusMult(Complex constant) {
    return new ComplexPlusMultSecond(constant);
  }
}
