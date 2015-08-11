library cern.colt.function.int;

import 'dart:math' as math;
import '../math.dart' show MAX_INT;

typedef int IntFunction(int argument);
typedef int IntIntFunction(int x, int y);
typedef int IntIntIntFunction(int first, int second, int third);
//typedef bool IntIntIntProcedure(int first, int second, int third);
//typedef bool IntIntProcedure(int first, int second);
//typedef bool IntProcedure(int element);

/// Unary functions

/// Function that returns `Math.abs(a) == (a < 0) ? -a : a`.
int abs(int a) => (a < 0) ? -a : a;

/// Function that returns `a--`.
int dec(int a) => a--;

/// Function that returns `Arithmetic.factorial(a)`.
//int factorial(int a) => DoubleArithmetic.factorial(a).toInt();

/// Function that returns its argument.
int identity(int a) => a;

/// Function that returns `a++`.
int inc(int a) => a++;

/// Function that returns `-a`.
int neg(int a) => -a;

/// Function that returns `~a`.
int not(int a) => a == 0 ? 1 : 0;

/// Function that returns `a < 0 ? -1 : a > 0 ? 1 : 0`.
int sign(int a) => a < 0 ? -1 : a > 0 ? 1 : 0;

/// Function that returns `a * a`.
int square(int a) => a * a;

// Binary functions

/// Function that returns `a & b`.
int and(int a, int b) => a & b;

/// Function that returns `a < b ? -1 : a > b ? 1 : 0`.
int compare(int a, int b) => a < b ? -1 : a > b ? 1 : 0;

/// Function that returns `a / b`.
int div(int a, int b) => a ~/ b;

/// Function that returns `-(a / b)`.
int divNeg(int a, int b) => -(a ~/ b);

/// Function that returns `a == b ? 1 : 0`.
int equals(int a, int b) => a == b ? 1 : 0;

/// Function that returns `a == b`.
bool isEqual(int a, int b) => a == b;

/// Function that returns `a < b`.
bool isLess(int a, int b) => a < b;

/// Function that returns `a > b`.
bool isGreater(int a, int b) => a > b;

/// Function that returns `Math.max(a,b)`.
int max(int a, int b) => (a >= b) ? a : b;

/// Function that returns `Math.min(a,b)`.
int min(int a, int b) => (a <= b) ? a : b;

/// Function that returns `a - b`.
int minus(int a, int b) => a - b;

/// Function that returns `a % b`.
int mod(int a, int b) => a % b;

/// Function that returns `a * b`.
int mult(int a, int b) => a * b;

/// Function that returns `-(a * b)`.
int multNeg(int a, int b) => -(a * b);

/// Function that returns `a * b^2`.
int multSquare(int a, int b) => a * b * b;

/// Function that returns `a | b`.
int or(int a, int b) => a | b;

/// Function that returns `a + b`.
int plus(int a, int b) => a + b;

/// Function that returns `Math.abs(a) + Math.abs(b)`.
int plusAbs(int a, int b) => a.abs() + b.abs();

/// Function that returns `(int) Math.pow(a,b)`.
int pow(int a, int b) => math.pow(a, b);

/// Function that returns `a << b`.
int shiftLeft(int a, int b) => a << b;

/// Function that returns `a >> b`.
int shiftRightSigned(int a, int b) => a >> b;

/// Function that returns `a ^ b`.
int xor(int a, int b) => a ^ b;

/// Constructs a function that returns `a & b`. `a` is a
/// variable, `b` is fixed.
IntFunction ampersand(int b) {
  return (int a) => a & b;
}

/// Constructs a function that returns `(from<=a && a<=to) ? 1 : 0`.
/// `a` is a variable, `from` and `to` are fixed.
IntFunction between(int from, int to) {
  return (int a) => (from <= a && a <= to) ? 1 : 0;
}

/// Constructs a unary function from a binary function with the first operand
/// (argument) fixed to the given constant `c`. The second operand is
/// variable (free).
IntFunction bindArg1(IntIntFunction function, int c) {
  return (int v) => function(c, v);
}

/// Constructs a unary function from a binary function with the second
/// operand (argument) fixed to the given constant `c`. The first
/// operand is variable (free).
IntFunction bindArg2(IntIntFunction function, int c) {
  return (int v) => function(v, c);
}

/// Constructs the function `g( h(a) )`.
IntFunction chain2(IntFunction g, IntFunction h) {
  return (int a) => g(h(a));
}

/// Constructs the function `g( h(a,b) )`.
IntIntFunction chain3(IntFunction g, IntIntFunction h) {
  return (int a, int b) => g(h(a, b));
}

/// Constructs the function `f( g(a), h(b) )`.
IntIntFunction chain4(IntIntFunction f, IntFunction g, IntFunction h) {
  return (int a, int b) => f(g(a), h(b));
}

/// Constructs a function that returns `a < b ? -1 : a > b ? 1 : 0`.
/// `a` is a variable, `b` is fixed.
IntFunction compareTo(int b) {
  return (int a) => a < b ? -1 : a > b ? 1 : 0;
}

/// Constructs a function that returns the constant `c`.
IntFunction constant(int c) {
  return (int a) => c;
}

/// Constructs a function that returns `a / b`. `a` is a
/// variable, `b` is fixed.
IntFunction divide(int b) {
  return (int a) => a ~/ b;
}

/// Constructs a function that returns `a == b ? 1 : 0`. `a` is
/// a variable, `b` is fixed.
IntFunction equalTo(int b) {
  return (int a) => a == b ? 1 : 0;
}

/// Constructs a function that returns `from<=a && a<=to`. `a`
/// is a variable, `from` and `to` are fixed.
//IntProcedure isBetween(int from, int to) {
//  return (int a) => from <= a && a <= to;
//}

/// Constructs a function that returns `a == b`. `a` is a
/// variable, `b` is fixed.
//IntProcedure isEqualTo(int b) {
//  return (int a) => a == b;
//}

/// Constructs a function that returns `a > b`. `a` is a
/// variable, `b` is fixed.
//IntProcedure isGreaterThan(int b) {
//  return (int a) => a > b;
//}

/// Constructs a function that returns `a < b`. `a` is a
/// variable, `b` is fixed.
//IntProcedure isLessThan(int b) {
//  return (int a) => a < b;
//}

/// Constructs a function that returns `Math.max(a,b)`. `a` is
/// a variable, `b` is fixed.
IntFunction maximum(int b) {
  return (int a) => (a >= b) ? a : b;
}

/// Constructs a function that returns `Math.min(a,b)`. `a` is
/// a variable, `b` is fixed.
IntFunction minimum(int b) {
  return (int a) => (a <= b) ? a : b;
}

/// Constructs a function that returns `a - b`. `a` is a
/// variable, `b` is fixed.
IntFunction subtract(int b) {
  return (int a) => a - b;
}

/// Constructs a function that returns `a - b*constant`. `a`
/// and `b` are variables, `constant` is fixed.
IntIntFunction minusMult(int constant) {
  return plusMultSecond(-constant);
}

/// Constructs a function that returns `a % b`. `a` is a
/// variable, `b` is fixed.
IntFunction modulus(int b) {
  return (int a) => a % b;
}

/// Constructs a function that returns `a * b`. `a` is a
/// variable, `b` is fixed.
IntFunction multiply(int b) {
  return (int a) => a * b;
}

/// Constructs a function that returns `a | b`. `a` is a
/// variable, `b` is fixed.
IntFunction bar(int b) {
  return (int a) => a | b;
}

/// Constructs a function that returns `a + b`. `a` is a
/// variable, `b` is fixed.
IntFunction add(int b) {
  return (int a) => a + b;
}

/// Constructs a function that returns `b*constant`.
IntIntFunction multSecond(int constant) {
  return (int a, int b) => b * constant;
}

/// Constructs a function that returns `Math.pow(a,b)`.
/// `a` is a variable, `b` is fixed.
IntFunction power(int b) {
  return (int a) => math.pow(a, b);
}

/// Constructs a function that returns `a + b*constant`. `a`
/// and `b` are variables, `constant` is fixed.
IntIntFunction plusMultSecond(int constant) {
  return new IntPlusMultSecond(constant);
}

/// Constructs a function that returns `a * constant + b`. `a`
/// and `b` are variables, `constant` is fixed.
IntIntFunction plusMultFirst(int constant) {
  return new IntPlusMultFirst(constant);
}

final _r = new math.Random();

IntFunction random() {
  return (_) => _r.nextInt(MAX_INT);
}

/// Constructs a function that returns `a << b`. `a` is a
/// variable, `b` is fixed.
IntFunction leftShift(int b) {
  return (int a) => a << b;
}

/// Constructs a function that returns `a >> b`. `a` is a
/// variable, `b` is fixed.
IntFunction signedRightShift(int b) {
  return (int a) => a >> b;
}

/// Constructs a function that returns `function(b,a)`, i.e.
/// applies the function with the first operand as second operand and the
/// second operand as first operand.
IntIntFunction swapArgs(IntIntFunction function) {
  return (int a, int b) => function(b, a);
}

/// Constructs a function that returns `a | b`. `a` is a
/// variable, `b` is fixed.
IntFunction exclusiveOr(int b) {
  return (int a) => a ^ b;
}

/// Only for performance tuning of compute intensive linear algebraic
/// computations. Constructs functions that return one of:
/// - `a * constant`
/// - `a / constant`
/// `a` is variable, `constant` is fixed, but for performance
/// reasons publicly accessible. Intended to be passed to
/// `matrix.apply(function)` methods.
class IntMult {
  /// Public read/write access to avoid frequent object construction.
  int multiplicator;

  IntMult(this.multiplicator);

  /// Returns the result of the function evaluation.
  int call(int a) {
    return a * multiplicator;
  }

  /// `a / constant`.
  factory IntMult.div(int constant) {
    return new IntMult.mult(1 ~/ constant);
  }

  /// `a * constant`.
  factory IntMult.mult(int constant) {
    return new IntMult(constant);
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
class IntPlusMultFirst {
  /// Public read/write access to avoid frequent object construction.
  int multiplicator;

  IntPlusMultFirst(this.multiplicator);

  /// Returns the result of the function evaluation.
  int call(int a, int b) {
    return a * multiplicator + b;
  }

  /// `a - b*constant`.
  factory IntPlusMultFirst.minusMult(int constant) {
    return new IntPlusMultFirst(-constant);
  }

  /// `a + b*constant`.
  factory IntPlusMultFirst.plusMult(int constant) {
    return new IntPlusMultFirst(constant);
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
class IntPlusMultSecond {
  /// Public read/write access to avoid frequent object construction.
  int multiplicator;

  IntPlusMultSecond(this.multiplicator);

  /// Returns the result of the function evaluation.
  int call(int a, int b) {
    return a + b * multiplicator;
  }

  /// `a - b*constant`.
  factory IntPlusMultSecond.minusMult(int constant) {
    return new IntPlusMultSecond(-constant);
  }

  /// `a + b*constant`.
  factory IntPlusMultSecond.plusMult(int constant) {
    return new IntPlusMultSecond(constant);
  }
}
