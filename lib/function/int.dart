library cern.colt.function.int;

import 'dart:math' as math;
import '../math.dart' show MAX_INT;

typedef int IntFunction(int argument);
typedef int IntIntFunction(int x, int y);
typedef int IntIntIntFunction(int first, int second, int third);
typedef bool IntIntIntProcedure(int first, int second, int third);
typedef bool IntIntProcedure(int first, int second);
typedef bool IntProcedure(int element);


/***************************************************************************
 * <H3>Unary functions</H3>
 **************************************************************************/
/**
 * Function that returns <tt>Math.abs(a) == (a < 0) ? -a : a</tt>.
 */
int abs(int a) => (a < 0) ? -a : a;

/**
 * Function that returns <tt>a--</tt>.
 */
int dec(int a) => a--;

/**
 * Function that returns <tt>(int) Arithmetic.factorial(a)</tt>.
 */
//int factorial(int a) => DoubleArithmetic.factorial(a).toInt();

/**
 * Function that returns its argument.
 */
int identity(int a) => a;

/**
 * Function that returns <tt>a++</tt>.
 */
int inc(int a) => a++;

/**
 * Function that returns <tt>-a</tt>.
 */
int neg(int a) => -a;

/**
 * Function that returns <tt>~a</tt>.
 */
int not(int a) => a == 0 ? 1 : 0;

/**
 * Function that returns <tt>a < 0 ? -1 : a > 0 ? 1 : 0</tt>.
 */
int sign(int a) => a < 0 ? -1 : a > 0 ? 1 : 0;

/**
 * Function that returns <tt>a * a</tt>.
 */
int square(int a) => a * a;

/***************************************************************************
 * <H3>Binary functions</H3>
 **************************************************************************/

/**
 * Function that returns <tt>a & b</tt>.
 */
int and(int a, int b) => a & b;

/**
 * Function that returns <tt>a < b ? -1 : a > b ? 1 : 0</tt>.
 */
int compare(int a, int b) => a < b ? -1 : a > b ? 1 : 0;

/**
 * Function that returns <tt>a / b</tt>.
 */
int div(int a, int b) => a ~/ b;

/**
 * Function that returns <tt>-(a / b)</tt>.
 */
int divNeg(int a, int b) => -(a ~/ b);

/**
 * Function that returns <tt>a == b ? 1 : 0</tt>.
 */
int equals(int a, int b) => a == b ? 1 : 0;

/**
 * Function that returns <tt>a == b</tt>.
 */
bool isEqual(int a, int b) => a == b;

/**
 * Function that returns <tt>a < b</tt>.
 */
bool isLess(int a, int b) => a < b;

/**
 * Function that returns <tt>a > b</tt>.
 */
bool isGreater(int a, int b) => a > b;

/**
 * Function that returns <tt>Math.max(a,b)</tt>.
 */
int max(int a, int b) => (a >= b) ? a : b;

/**
 * Function that returns <tt>Math.min(a,b)</tt>.
 */
int min(int a, int b) => (a <= b) ? a : b;

/**
 * Function that returns <tt>a - b</tt>.
 */
int minus(int a, int b) => a - b;

/**
 * Function that returns <tt>a % b</tt>.
 */
int mod(int a, int b) => a % b;

/**
 * Function that returns <tt>a * b</tt>.
 */
int mult(int a, int b) => a * b;

/**
 * Function that returns <tt>-(a * b)</tt>.
 */
int multNeg(int a, int b) => -(a * b);

/**
 * Function that returns <tt>a * b^2</tt>.
 */
int multSquare(int a, int b) => a * b * b;

/**
 * Function that returns <tt>a | b</tt>.
 */
int or(int a, int b) => a | b;

/**
 * Function that returns <tt>a + b</tt>.
 */
int plus(int a, int b) => a + b;

/**
 * Function that returns <tt>Math.abs(a) + Math.abs(b)</tt>.
 */
int plusAbs(int a, int b) => a.abs() + b.abs();

/**
 * Function that returns <tt>(int) Math.pow(a,b)</tt>.
 */
int pow(int a, int b) => math.pow(a, b);

/**
 * Function that returns <tt>a << b</tt>.
 */
int shiftLeft(int a, int b) => a << b;

/**
 * Function that returns <tt>a >> b</tt>.
 */
int shiftRightSigned(int a, int b) => a >> b;

/**
 * Function that returns <tt>a >>> b</tt>.
 */
//int shiftRightUnsigned(int a, int b) => a >>> b;

/**
 * Function that returns <tt>a ^ b</tt>.
 */
int xor(int a, int b) => a ^ b;

/**
 * Constructs a function that returns <tt>a & b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
IntFunction ampersand(final int b) {
  return (int a) => a & b;
}

/**
 * Constructs a function that returns <tt>(from<=a && a<=to) ? 1 : 0</tt>.
 * <tt>a</tt> is a variable, <tt>from</tt> and <tt>to</tt> are fixed.
 */
IntFunction between(final int from, final int to) {
  return (int a) => (from <= a && a <= to) ? 1 : 0;
}

/**
 * Constructs a unary function from a binary function with the first operand
 * (argument) fixed to the given constant <tt>c</tt>. The second operand is
 * variable (free).
 *
 * @param function
 *            a binary function taking operands in the form
 *            <tt>function.apply(c,var)</tt>.
 * @return the unary function <tt>function(c,var)</tt>.
 */
IntFunction bindArg1(final IntIntFunction function, final int c) {
  return (int v) => function(c, v);
}

/**
 * Constructs a unary function from a binary function with the second
 * operand (argument) fixed to the given constant <tt>c</tt>. The first
 * operand is variable (free).
 *
 * @param function
 *            a binary function taking operands in the form
 *            <tt>function.apply(var,c)</tt>.
 * @return the unary function <tt>function(var,c)</tt>.
 */
IntFunction bindArg2(final IntIntFunction function, final int c) {
  return (int v) => function(v, c);
}

/**
 * Constructs the function <tt>g( h(a) )</tt>.
 *
 * @param g
 *            a unary function.
 * @param h
 *            a unary function.
 * @return the unary function <tt>g( h(a) )</tt>.
 */
IntFunction chain2(final IntFunction g, final IntFunction h) {
  return (int a) => g(h(a));
}

/**
 * Constructs the function <tt>g( h(a,b) )</tt>.
 *
 * @param g
 *            a unary function.
 * @param h
 *            a binary function.
 * @return the unary function <tt>g( h(a,b) )</tt>.
 */
IntIntFunction chain3(final IntFunction g, final IntIntFunction h) {
  return (int a, int b) => g(h(a, b));
}

/**
 * Constructs the function <tt>f( g(a), h(b) )</tt>.
 *
 * @param f
 *            a binary function.
 * @param g
 *            a unary function.
 * @param h
 *            a unary function.
 * @return the binary function <tt>f( g(a), h(b) )</tt>.
 */
IntIntFunction chain4(final IntIntFunction f, final IntFunction g, final IntFunction h) {
  return (int a, int b) => f(g(a), h(b));
}

/**
 * Constructs a function that returns <tt>a < b ? -1 : a > b ? 1 : 0</tt>.
 * <tt>a</tt> is a variable, <tt>b</tt> is fixed.
 */
IntFunction compareTo(final int b) {
  return (int a) => a < b ? -1 : a > b ? 1 : 0;
}

/**
 * Constructs a function that returns the constant <tt>c</tt>.
 */
IntFunction constant(final int c) {
  return (int a) => c;
}

/**
 * Constructs a function that returns <tt>a / b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
IntFunction divide(final int b) {
  return (int a) => a ~/ b;
}

/**
 * Constructs a function that returns <tt>a == b ? 1 : 0</tt>. <tt>a</tt> is
 * a variable, <tt>b</tt> is fixed.
 */
IntFunction equalTo(final int b) {
  return (int a) => a == b ? 1 : 0;
}

/**
 * Constructs a function that returns <tt>from<=a && a<=to</tt>. <tt>a</tt>
 * is a variable, <tt>from</tt> and <tt>to</tt> are fixed.
 */
IntProcedure isBetween(final int from, final int to) {
  return (int a) => from <= a && a <= to;
}

/**
 * Constructs a function that returns <tt>a == b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
IntProcedure isEqualTo(final int b) {
  return (int a) => a == b;
}

/**
 * Constructs a function that returns <tt>a > b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
IntProcedure isGreaterThan(final int b) {
  return (int a) => a > b;
}

/**
 * Constructs a function that returns <tt>a < b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
IntProcedure isLessThan(final int b) {
  return (int a) => a < b;
}

/**
 * Constructs a function that returns <tt>Math.max(a,b)</tt>. <tt>a</tt> is
 * a variable, <tt>b</tt> is fixed.
 */
IntFunction maximum(final int b) {
  return (int a) => (a >= b) ? a : b;
}

/**
 * Constructs a function that returns <tt>Math.min(a,b)</tt>. <tt>a</tt> is
 * a variable, <tt>b</tt> is fixed.
 */
IntFunction minimum(final int b) {
  return (int a) => (a <= b) ? a : b;
}

/**
 * Constructs a function that returns <tt>a - b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
IntFunction subtract(final int b) {
  return (int a) => a - b;
}

/**
 * Constructs a function that returns <tt>a - b*constant</tt>. <tt>a</tt>
 * and <tt>b</tt> are variables, <tt>constant</tt> is fixed.
 */
IntIntFunction minusMult(final int constant) {
  return plusMultSecond(-constant);
}

/**
 * Constructs a function that returns <tt>a % b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
IntFunction modulus(final int b) {
  return (int a) => a % b;
}

/**
 * Constructs a function that returns <tt>a * b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
IntFunction multiply(final int b) {
  return (int a) => a * b;
}

/**
 * Constructs a function that returns <tt>a | b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
IntFunction bar(final int b) {
  return (int a) => a | b;
}

/**
 * Constructs a function that returns <tt>a + b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
IntFunction add(final int b) {
  return (int a) => a + b;
}

/**
 * Constructs a function that returns <tt>b*constant</tt>.
 */
IntIntFunction multSecond(final int constant) {
  return (int a, int b) => b * constant;
}

/**
 * Constructs a function that returns <tt>(int) Math.pow(a,b)</tt>.
 * <tt>a</tt> is a variable, <tt>b</tt> is fixed.
 */
IntFunction power(final int b) {
  return (int a) => math.pow(a, b);
}

/**
 * Constructs a function that returns <tt>a + b*constant</tt>. <tt>a</tt>
 * and <tt>b</tt> are variables, <tt>constant</tt> is fixed.
 */
IntIntFunction plusMultSecond(final int constant) {
  return new IntPlusMultSecond(constant);
}

/**
 * Constructs a function that returns <tt>a * constant + b</tt>. <tt>a</tt>
 * and <tt>b</tt> are variables, <tt>constant</tt> is fixed.
 */
IntIntFunction plusMultFirst(final int constant) {
  return new IntPlusMultFirst(constant);
}

IntFunction random() {
  final r = new math.Random();
  return (_) => r.nextInt(MAX_INT);
    //return new cern.jet.random.tdouble.engine.DoubleMersenneTwister(new java.util.Date());
}

/**
 * Constructs a function that returns <tt>a << b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
IntFunction leftShift(final int b) {
  return (int a) => a << b;
}

/**
 * Constructs a function that returns <tt>a >> b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
IntFunction signedRightShift(final int b) {
  return (int a) => a >> b;
}

/**
 * Constructs a function that returns <tt>a >>> b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
/*IntFunction unsignedRightShift(final int b) {
    return (int a) => a >>> b;
}*/

/**
 * Constructs a function that returns <tt>function.apply(b,a)</tt>, i.e.
 * applies the function with the first operand as second operand and the
 * second operand as first operand.
 *
 * @param function
 *            a function taking operands in the form
 *            <tt>function.apply(a,b)</tt>.
 * @return the binary function <tt>function(b,a)</tt>.
 */
IntIntFunction swapArgs(final IntIntFunction function) {
  return (int a, int b) => function(b, a);
}

/**
 * Constructs a function that returns <tt>a | b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
IntFunction exclusiveOr(final int b) {
  return (int a) => a ^ b;
}

/**
 * Only for performance tuning of compute intensive linear algebraic
 * computations. Constructs functions that return one of
 * <ul>
 * <li><tt>a * constant</tt>
 * <li><tt>a / constant</tt>
 * </ul>
 * <tt>a</tt> is variable, <tt>constant</tt> is fixed, but for performance
 * reasons publicly accessible. Intended to be passed to
 * <tt>matrix.assign(function)</tt> methods.
 */
class IntMult {
  /**
   * Public read/write access to avoid frequent object construction.
   */
  int multiplicator;

  /**
   * Insert the method's description here. Creation date: (8/10/99 19:12:09)
   */
  IntMult(this.multiplicator);

  /**
   * Returns the result of the function evaluation.
   */
  int call(int a) {
    return a * multiplicator;
  }

  /**
   * <tt>a / constant</tt>.
   */
  factory IntMult.div(final int constant) {
    return new IntMult.mult(1 ~/ constant);
  }

  /**
   * <tt>a * constant</tt>.
   */
  factory IntMult.mult(final int constant) {
    return new IntMult(constant);
  }
}

/**
 * Only for performance tuning of compute intensive linear algebraic
 * computations. Constructs functions that return one of
 * <ul>
 * <li><tt>a + b*constant</tt>
 * <li><tt>a - b*constant</tt>
 * <li><tt>a + b/constant</tt>
 * <li><tt>a - b/constant</tt>
 * </ul>
 * <tt>a</tt> and <tt>b</tt> are variables, <tt>constant</tt> is fixed, but for
 * performance reasons publicly accessible. Intended to be passed to
 * <tt>matrix.assign(otherMatrix,function)</tt> methods.
 */
class IntPlusMultFirst {
  /**
   * Public read/write access to avoid frequent object construction.
   */
  int multiplicator;

  /**
   * Insert the method's description here. Creation date: (8/10/99 19:12:09)
   */
  IntPlusMultFirst(this.multiplicator);

  /**
   * Returns the result of the function evaluation.
   */
  int call(int a, int b) {
    return a * multiplicator + b;
  }

  /**
   * <tt>a - b*constant</tt>.
   */
  factory IntPlusMultFirst.minusMult(final int constant) {
    return new IntPlusMultFirst(-constant);
  }

  /**
   * <tt>a + b*constant</tt>.
   */
  factory IntPlusMultFirst.plusMult(final int constant) {
    return new IntPlusMultFirst(constant);
  }
}

/**
 * Only for performance tuning of compute intensive linear algebraic
 * computations. Constructs functions that return one of
 * <ul>
 * <li><tt>a + b*constant</tt>
 * <li><tt>a - b*constant</tt>
 * <li><tt>a + b/constant</tt>
 * <li><tt>a - b/constant</tt>
 * </ul>
 * <tt>a</tt> and <tt>b</tt> are variables, <tt>constant</tt> is fixed, but for
 * performance reasons publicly accessible. Intended to be passed to
 * <tt>matrix.assign(otherMatrix,function)</tt> methods.
 */
class IntPlusMultSecond {
  /**
   * Public read/write access to avoid frequent object construction.
   */
  int multiplicator;

  /**
   * Insert the method's description here. Creation date: (8/10/99 19:12:09)
   */
  IntPlusMultSecond(this.multiplicator);

  /**
   * Returns the result of the function evaluation.
   */
  int call(int a, int b) {
    return a + b * multiplicator;
  }

  /**
   * <tt>a - b*constant</tt>.
   */
  factory IntPlusMultSecond.minusMult(final int constant) {
    return new IntPlusMultSecond(-constant);
  }

  /**
   * <tt>a + b*constant</tt>.
   */
  factory IntPlusMultSecond.plusMult(final int constant) {
    return new IntPlusMultSecond(constant);
  }
}
