library cern.colt.function.double;

import 'dart:math' as Math;

typedef double DoubleFunction(double x);
typedef double Double5Function(
    double a, double b, double c, double d, double e);
typedef double Double9Function(double a00, double a01, double a02, double a10,
    double a11, double a12, double a20, double a21, double a22);
typedef double DoubleDoubleFunction(double x, double y);
typedef double IntDoubleFunction(int first, double second);
typedef double IntIntDoubleFunction(int first, int second, double third);

double toRadians(double a) => a * (Math.PI / 180.0);

double toDegrees(double a) => a * (180.0 / Math.PI);

double or(double a, double b) => (a != 0.0 || b != 0.0) ? 1.0 : 0.0;

double and(double a, double b) => (a != 0.0 && b != 0.0) ? 1.0 : 0.0;

/// Unary functions

/// Function that returns `Math.abs(a)`.
double abs(double a) => a.abs();

/// Function that returns `Math.acos(a)`.
double acos(double a) => Math.acos(a);

/// Function that returns `acosh(a)`.
//double acosh(double a) => Math.acosh(a);

/// Function that returns `Math.asin(a)`.
double asin(double a) => Math.asin(a);

/// Function that returns `asinh(a)`.
//double asinh (double a) => Math.asinh(a);

/// Function that returns `Math.atan(a)`.
double atan(double a) => Math.atan(a);

/// Function that returns `atanh(a)`.
//double atanh(double a) => Math.atanh(a);

/// Function that returns `Math.ceil(a)`.
double ceil(double a) => a.ceilToDouble();

/// Function that returns `Math.cos(a)`.
double cos(double a) => Math.cos(a);

/// Function that returns `cosh(a)`.
//double cosh(double a) => Math.cosh(a);

/// Function that returns `cot(a)`.
//double cot(double a) => Math.cot(a);

/// Function that returns `erf(a)`.
//double erf(double a) => Math.erf(a);

/// Function that returns `erfc(a)`.
//double erfc(double a) => Math.erfc(a);

/// Function that returns `Math.exp(a)`.
double exp(double a) => Math.exp(a);

/// Function that returns `Math.floor(a)`.
double floor(double a) => a.floorToDouble();

/// Function that returns `gamma(a)`.
//double gamma(double a) => Math.gamma(a);

/// Function that returns its argument.
double identity(double a) => a;

/// Function that returns `1.0 / a`.
double inv(double a) => 1.0 / a;

/// Function that returns `Math.log(a)`.
double log(double a) => Math.log(a);

/// Function that returns `log10(a)`.
//double log10(double a) => Math.log10(a);

/// Function that returns `Math.log(a) / Math.log(2)`.
// 1.0 / Math.log(2) == 1.4426950408889634
double log2(double a) => Math.log(a) * 1.4426950408889634;

/// Function that returns `logGamma(a)`.
//double logGamma(double a) => Math.logGamma(a);

/// Function that returns `-a`.
double neg(double a) => -a;

/// Function that returns `Math.rint(a)`.
//double rint(double a) => Math.rint(a);

/// Function that returns `a < 0 ? -1 : a > 0 ? 1 : 0`.
double sign(double a) => a < 0 ? -1.0 : a > 0 ? 1.0 : 0.0;

/// Function that returns `Math.sin(a)`.
double sin(double a) => Math.sin(a);

/// Function that returns `sinh(a)`.
//double sinh(double a) => Math.sinh(a);

/// Function that returns `Math.sqrt(a)`.
double sqrt(double a) => Math.sqrt(a);

/// Function that returns `a * a`.
double square(double a) => a * a;

/// Function that returns `Math.tan(a)`.
double tan(double a) => Math.tan(a);

/// Function that returns `tanh(a)`.
//double tanh(double a) => Math.tanh(a);

/// Binary functions

/// Function that returns `Math.atan2(a,b)`.
double atan2(double a, double b) => Math.atan2(a, b);

/// Function that returns `logBeta(a,b)`.
//double logBeta(double a, double b) => Math.logBeta(a,b);

/// Function that returns `a < b ? -1 : a > b ? 1 : 0`.
double compare(double a, double b) => a < b ? -1.0 : a > b ? 1.0 : 0.0;

/// Function that returns `a / b`.
double div(double a, double b) => a / b;

/// Function that returns `-(a / b)`.
double divNeg(double a, double b) => -(a / b);

/// Function that returns `a == b ? 1 : 0`.
double equals(double a, double b) => a == b ? 1.0 : 0.0;

/// Function that returns `a > b ? 1 : 0`.
double greater(double a, double b) => a > b ? 1.0 : 0.0;

/// Function that returns `Math.IEEEremainder(a,b)`.
//double IEEEremainder(double a, double b) => Math.IEEEremainder(a, b);

/// Function that returns `a == b`.
bool isEqual(double a, double b) => a == b;

/// Function that returns `a < b`.
bool isLess(double a, double b) => a < b;

/// Function that returns `a > b`.
bool isGreater(double a, double b) => a > b;

/// Function that returns `a < b ? 1 : 0`.
double less(double a, double b) => a < b ? 1.0 : 0.0;

/// Function that returns `Math.log(a) / Math.log(b)`.
double lg(double a, double b) => Math.log(a) / Math.log(b);

/// Function that returns `Math.max(a,b)`.
double max(double a, double b) => Math.max(a, b);

/// Function that returns `Math.min(a,b)`.
double min(double a, double b) => Math.min(a, b);

/// Function that returns `a - b`.
final DoubleDoubleFunction minus = plusMultSecond(-1.0);

/// Function that returns `a % b`.
double mod(double a, double b) => a % b;

/// Function that returns `a * b`.
double mult(double a, double b) => a * b;

/// Function that returns `-(a * b)`.
double multNeg(double a, double b) => -(a * b);

/// Function that returns `a * b^2`.
double multSquare(double a, double b) => a * b * b;

/// Function that returns `a + b`.
final DoubleDoubleFunction plus = plusMultSecond(1.0);

/// Function that returns `Math.abs(a) + Math.abs(b)`.
double plusAbs(double a, double b) => a.abs() + b.abs();

/// Function that returns `Math.pow(a,b)`.
double pow(double a, double b) => Math.pow(a, b);

/// Constructs a function that returns `(from<=a && a<=to) ? 1 : 0`.
/// `a` is a variable, `from` and `to` are fixed.
DoubleFunction between(double from, double to) {
  return (double a) => (from <= a && a <= to) ? 1.0 : 0.0;
}

/// Constructs a unary function from a binary function with the first operand
/// (argument) fixed to the given constant `c`. The second operand is
/// variable (free).
DoubleFunction bindArg1(DoubleDoubleFunction fn, double c) {
  return (double v) => fn(c, v);
}

/// Constructs a unary function from a binary function with the second
/// operand (argument) fixed to the given constant `c`. The first
/// operand is variable (free).
DoubleFunction bindArg2(DoubleDoubleFunction fn, double c) {
  return (double v) => fn(v, c);
}

/// Constructs the function `f( g(a), h(b) )`.
DoubleDoubleFunction chainFGH(
    DoubleDoubleFunction f, DoubleFunction g, DoubleFunction h) {
  return (double a, double b) => f(g(a), h(b));
}

/// Constructs the function `g( h(a,b) )`.
DoubleDoubleFunction chain(DoubleFunction g, DoubleDoubleFunction h) {
  return (double a, double b) => g(h(a, b));
}

/// Constructs the function `g( h(a) )`.
DoubleFunction chainGH(DoubleFunction g, DoubleFunction h) {
  return (double a) => g(h(a));
}

/// Constructs a function that returns `a < b ? -1 : a > b ? 1 : 0`.
/// `a` is a variable, `b` is fixed.
DoubleFunction compareTo(double b) {
  return (double a) => a < b ? -1 : a > b ? 1.0 : 0.0;
}

/// Constructs a function that returns the constant `c`.
DoubleFunction constant(double c) {
  return (double a) => c;
}

/// Constructs a function that returns `a / b`. `a` is a
/// variable, `b` is fixed.
DoubleFunction divide(double b) {
  return multiply(1.0 / b);
}

/// Constructs a function that returns `a == b ? 1 : 0`. `a` is
/// a variable, `b` is fixed.
DoubleFunction equalTo(double b) {
  return (double a) {
    return a == b ? 1.0 : 0.0;
  };
}

/// Constructs a function that returns `a > b ? 1 : 0`. `a` is
/// a variable, `b` is fixed.
DoubleFunction greaterThan(double b) {
  return (double a) {
    return a > b ? 1.0 : 0.0;
  };
}

/// Constructs a function that returns `Math.IEEEremainder(a,b)`.
/// `a` is a variable, `b` is fixed.
/*DoubleFunction IEEEremainderOf(double b) {
    return (double a) {
            return Math.IEEEremainder(a, b);
    };
}*/

/// Constructs a function that returns `from<=a && a<=to`. `a`
/// is a variable, `from` and `to` are fixed.
//DoubleProcedure isBetween(double from, double to) {
//    return (double a) => from <= a && a <= to;
//}

/// Constructs a function that returns `a == b`. `a` is a
/// variable, `b` is fixed.
//DoubleProcedure isEqualTo(double b) {
//    return (double a) => a == b;
//}

/// Constructs a function that returns `a > b`. `a` is a
/// variable, `b` is fixed.
//DoubleProcedure isGreaterThan(double b) {
//    return (double a) => a > b;
//}

/**
 * Constructs a function that returns `a < b`. `a` is a
 * variable, `b` is fixed.
 */
//DoubleProcedure isLessThan(double b) {
//    return (double a) => a < b;
//}

/// Constructs a function that returns `a < b ? 1 : 0`. `a` is
/// a variable, `b` is fixed.
DoubleFunction lessThan(double b) {
  return (double a) => a < b ? 1.0 : 0.0;
}

/// Constructs a function that returns ``Math.log(a) / Math.log(b)`
/// `. `a` is a variable, `b` is fixed.
DoubleFunction lgFn(double b) {
  double logInv = 1 / Math.log(b); // cached for speed
  return (double a) => Math.log(a) * logInv;
}

/// Constructs a function that returns `Math.max(a,b)`. `a` is
/// a variable, `b` is fixed.
DoubleFunction maximum(double b) {
  return (double a) => Math.max(a, b);
}

/// Constructs a function that returns `Math.min(a,b)`. `a` is
/// a variable, `b` is fixed.
DoubleFunction minimum(double b) {
  return (double a) => Math.min(a, b);
}

/// Constructs a function that returns `a - b`. `a` is a
/// variable, `b` is fixed.
DoubleFunction subtract(num b) {
  return add(-b.toDouble());
}

/// Constructs a function that returns `a - b*constant`. `a`
/// and `b` are variables, `constant` is fixed.
DoubleDoubleFunction minusMult(double constant) {
  return plusMultSecond(-constant);
}

/// Constructs a function that returns `a % b`. `a` is a
/// variable, `b` is fixed.
DoubleFunction modulus(double b) {
  return (double a) {
    return a % b;
  };
}

/// Constructs a function that returns `a * b`. `a` is a
/// variable, `b` is fixed.
DoubleFunction multiply(double b) {
  return new DoubleMult(b);
}

/// Constructs a function that returns `a + b`. `a` is a
/// variable, `b` is fixed.
DoubleFunction add(double b) {
  return (double a) => a + b;
}

/// Constructs a function that returns `b*constant`.
DoubleDoubleFunction multSecond(double constant) {
  return (double a, double b) => b * constant;
}

/// Constructs a function that returns `a + b*constant`. `a`
/// and `b` are variables, `constant` is fixed.
DoubleDoubleFunction plusMultSecond(double constant) {
  return new DoublePlusMultSecond(constant);
}

/// Constructs a function that returns `a * constant + b`. `a`
/// and `b` are variables, `constant` is fixed.
DoubleDoubleFunction plusMultFirst(double constant) {
  return new DoublePlusMultFirst(constant);
}

/// Constructs a function that returns `Math.pow(a,b)`. `a` is
/// a variable, `b` is fixed.
DoubleFunction power(double b) {
  return (double a) => Math.pow(a, b);
}

/// Constructs a function that returns a new uniform random number in the
/// open unit interval <code>(0.0,1.0)</code> (excluding 0.0 and 1.0).
DoubleFunction random() {
  return RandomDoubleFunction;
}

final _r = new Math.Random();

double RandomDoubleFunction(double argument) {
  return _r.nextDouble();
}

/// Constructs a function that returns the number rounded to the given
/// precision; `Math.rint(a/precision)*precision`. Examples:
///     precision = 0.01 rounds 0.012 --> 0.01, 0.018 --> 0.02
///     precision = 10   rounds 123   --> 120 , 127   --> 130
DoubleFunction round(double precision) {
  return (double a) => (a / precision).round() * precision;
}

/// Constructs a function that returns `fn(b,a)`, i.e.
/// applies the function with the first operand as second operand and the
/// second operand as first operand.
DoubleDoubleFunction swapArgs(final DoubleDoubleFunction fn) {
  return (double a, double b) => fn(b, a);
}

/// Only for performance tuning of compute intensive linear algebraic
/// computations. Constructs functions that return one of:
/// - `a * constant`
/// - `a / constant`
/// `a` is variable, `constant` is fixed, but for performance
/// reasons publicly accessible. Intended to be passed to
/// `matrix.apply(function)` methods.
class DoubleMult {
    /// Public read/write access to avoid frequent object construction.
    double multiplicator;

    DoubleMult(this.multiplicator);

    /// Returns the result of the function evaluation.
    double call(double a) => a * multiplicator;

    /// `a / constant`.
    factory DoubleMult.div(final double constant) {
        return new DoubleMult.mult(1 / constant);
    }

    /// `a * constant`.
    factory DoubleMult.mult(final double constant) {
        return new DoubleMult(constant);
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
class DoublePlusMultFirst {
    /// Public read/write access to avoid frequent object construction.
    double multiplicator;

    DoublePlusMultFirst(this.multiplicator);

    /// Returns the result of the function evaluation.
    double call(double a, double b) => a * multiplicator + b;

    /// `a - b/constant`.
    factory DoublePlusMultFirst.minusDiv(final double constant) {
        return new DoublePlusMultFirst(-1 / constant);
    }

    /// `a - b*constant`.
    factory DoublePlusMultFirst.minusMult(final double constant) {
        return new DoublePlusMultFirst(-constant);
    }

    /// `a + b/constant`.
    factory DoublePlusMultFirst.plusDiv(final double constant) {
        return new DoublePlusMultFirst(1 / constant);
    }

    /// `a + b*constant`.
    factory DoublePlusMultFirst.plusMult(final double constant) {
        return new DoublePlusMultFirst(constant);
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
class DoublePlusMultSecond {
    /// Public read/write access to avoid frequent object construction.
    double multiplicator;

    DoublePlusMultSecond(this.multiplicator);

    /// Returns the result of the function evaluation.
    double call(double a, double b) => a + b * multiplicator;

    /// `a - b/constant`.
    factory DoublePlusMultSecond.minusDiv(final double constant) {
        return new DoublePlusMultSecond(-1 / constant);
    }

    /// `a - b*constant`.
    factory DoublePlusMultSecond.minusMult(final double constant) {
        return new DoublePlusMultSecond(-constant);
    }

    /// `a + b/constant`.
    factory DoublePlusMultSecond.plusDiv(final double constant) {
        return new DoublePlusMultSecond(1 / constant);
    }

    /// `a + b*constant`.
    factory DoublePlusMultSecond.plusMult(final double constant) {
        return new DoublePlusMultSecond(constant);
    }
}
