library cern.colt.function.double;

import 'dart:math' as Math;
import '../math.dart';

typedef double DoubleFunction(double x);
typedef double Double5Function(double a, double b, double c, double d, double e);
typedef double Double9Function(double a00, double a01, double a02, double a10, double a11, double a12, double a20,
    double a21, double a22);
typedef double DoubleDoubleFunction(double x, double y);
typedef bool DoubleDoubleProcedure(double first, double second);
typedef bool DoubleIntProcedure(double first, int second);
typedef bool DoubleProcedure(double element);
typedef double IntDoubleFunction(int first, double second);
typedef bool IntDoubleProcedure(int first, double second);
typedef double IntIntDoubleFunction(int first, int second, double third);
typedef bool IntIntDoubleProcedure(int first, int second, double third);


/***************************************************************************
 * <H3>Unary functions</H3>
 **************************************************************************/
/**
 * Function that returns <tt>Math.abs(a)</tt>.
 */
double abs(double a) => a.abs();

/**
 * Function that returns <tt>Math.acos(a)</tt>.
 */
double acos (double a) => Math.acos(a);

/**
 * Function that returns <tt>com.imsl.math.Sfun.acosh(a)</tt>.
 */
//double acosh(double a) => Math.acosh(a);

/**
 * Function that returns <tt>Math.asin(a)</tt>.
 */
double asin(double a) => Math.asin(a);

/**
 * Function that returns <tt>com.imsl.math.Sfun.asinh(a)</tt>.
 */
//double asinh (double a) => Math.asinh(a);

/**
 * Function that returns <tt>Math.atan(a)</tt>.
 */
double atan(double a) => Math.atan(a);

/**
 * Function that returns <tt>com.imsl.math.Sfun.atanh(a)</tt>.
 */
//double atanh(double a) => Math.atanh(a);

/**
 * Function that returns <tt>Math.ceil(a)</tt>.
 */
double ceil(double a) => a.ceilToDouble();

/**
 * Function that returns <tt>Math.cos(a)</tt>.
 */
double cos(double a) => Math.cos(a);

/**
 * Function that returns <tt>com.imsl.math.Sfun.cosh(a)</tt>.
 */
//double cosh(double a) => Math.cosh(a);

/**
 * Function that returns <tt>com.imsl.math.Sfun.cot(a)</tt>.
 */
//double cot(double a) => Math.cot(a);

/**
 * Function that returns <tt>com.imsl.math.Sfun.erf(a)</tt>.
 */
//double erf(double a) => Math.erf(a);

/**
 * Function that returns <tt>com.imsl.math.Sfun.erfc(a)</tt>.
 */
//double erfc(double a) => Math.erfc(a);

/**
 * Function that returns <tt>Math.exp(a)</tt>.
 */
double exp(double a) => Math.exp(a);

/**
 * Function that returns <tt>Math.floor(a)</tt>.
 */
double floor(double a) => a.floorToDouble();

/**
 * Function that returns <tt>com.imsl.math.Sfun.gamma(a)</tt>.
 */
//double gamma(double a) => Math.gamma(a);

/**
 * Function that returns its argument.
 */
double identity(double a) => a;

/**
 * Function that returns <tt>1.0 / a</tt>.
 */
double inv(double a) => 1.0 / a;

/**
 * Function that returns <tt>Math.log(a)</tt>.
 */
double log(double a) => Math.log(a);

/**
 * Function that returns <tt>com.imsl.math.Sfun.log10(a)</tt>.
 */
//double log10(double a) => Math.log10(a);

/**
 * Function that returns <tt>Math.log(a) / Math.log(2)</tt>.
 */
// 1.0 / Math.log(2) == 1.4426950408889634
double log2(double a) => Math.log(a) * 1.4426950408889634;

/**
 * Function that returns <tt>com.imsl.math.Sfun.logGamma(a)</tt>.
 */
//double logGamma(double a) => Math.logGamma(a);

/**
 * Function that returns <tt>-a</tt>.
 */
double neg(double a) => -a;

/**
 * Function that returns <tt>Math.rint(a)</tt>.
 */
//double rint(double a) => Math.rint(a);

/**
 * Function that returns <tt>a < 0 ? -1 : a > 0 ? 1 : 0</tt>.
 */
double sign(double a) => a < 0 ? -1.0 : a > 0 ? 1.0 : 0.0;

/**
 * Function that returns <tt>Math.sin(a)</tt>.
 */
double sin(double a) => Math.sin(a);

/**
 * Function that returns <tt>com.imsl.math.Sfun.sinh(a)</tt>.
 */
/*
 * DoubleFunction sinh = new DoubleFunction() { public
 * final double apply(double a) { return Sfun.sinh(a); } };
 */

/**
 * Function that returns <tt>Math.sqrt(a)</tt>.
 */
double sqrt(double a) => Math.sqrt(a);

/**
 * Function that returns <tt>a * a</tt>.
 */
double square(double a) => a * a;

/**
 * Function that returns <tt>Math.tan(a)</tt>.
 */
double tan(double a) => Math.tan(a);

/**
 * Function that returns <tt>com.imsl.math.Sfun.tanh(a)</tt>.
 */
//double tanh(double a) => Math.tanh(a);

/**
 * Function that returns <tt>Math.toDegrees(a)</tt>.
 */
//double toDegrees(double a) => Math.toDegrees(a);

/**
 * Function that returns <tt>Math.toRadians(a)</tt>.
 */
//double toRadians(double a) => Math.toRadians(a);

/***************************************************************************
 * <H3>Binary functions</H3>
 **************************************************************************/

/**
 * Function that returns <tt>Math.atan2(a,b)</tt>.
 */
double atan2(double a, double b) => Math.atan2(a, b);

/**
 * Function that returns <tt>com.imsl.math.Sfun.logBeta(a,b)</tt>.
 */
//double logBeta(double a, double b) => Math.logBeta(a,b);

/**
 * Function that returns <tt>a < b ? -1 : a > b ? 1 : 0</tt>.
 */
double compare(double a, double b) => a < b ? -1.0 : a > b ? 1.0 : 0.0;

/**
 * Function that returns <tt>a / b</tt>.
 */
double div(double a, double b) => a / b;

/**
 * Function that returns <tt>-(a / b)</tt>.
 */
double divNeg(double a, double b) => -(a / b);

/**
 * Function that returns <tt>a == b ? 1 : 0</tt>.
 */
double equals(double a, double b) => a == b ? 1.0 : 0.0;

/**
 * Function that returns <tt>a > b ? 1 : 0</tt>.
 */
double greater(double a, double b) => a > b ? 1.0 : 0.0;

/**
 * Function that returns <tt>Math.IEEEremainder(a,b)</tt>.
 */
//double IEEEremainder(double a, double b) => Math.IEEEremainder(a, b);

/**
 * Function that returns <tt>a == b</tt>.
 */
bool isEqual(double a, double b) => a == b;

/**
 * Function that returns <tt>a < b</tt>.
 */
bool isLess(double a, double b) => a < b;

/**
 * Function that returns <tt>a > b</tt>.
 */
bool isGreater(double a, double b) => a > b;

/**
 * Function that returns <tt>a < b ? 1 : 0</tt>.
 */
double less(double a, double b) => a < b ? 1.0 : 0.0;

/**
 * Function that returns <tt>Math.log(a) / Math.log(b)</tt>.
 */
double lg(double a, double b) => Math.log(a) / Math.log(b);

/**
 * Function that returns <tt>Math.max(a,b)</tt>.
 */
double max(double a, double b) => Math.max(a, b);

/**
 * Function that returns <tt>Math.min(a,b)</tt>.
 */
double min(double a, double b) => Math.min(a, b);

/**
 * Function that returns <tt>a - b</tt>.
 */
final DoubleDoubleFunction minus = plusMultSecond(-1.0);

/*
 * new DoubleDoubleFunction() { public final double apply(double a, double
 * b) { return a - b; } };
 */

/**
 * Function that returns <tt>a % b</tt>.
 */
double mod(double a, double b) => a % b;

/**
 * Function that returns <tt>a * b</tt>.
 */
double mult(double a, double b) => a * b;

/**
 * Function that returns <tt>-(a * b)</tt>.
 */
double multNeg(double a, double b) => -(a * b);

/**
 * Function that returns <tt>a * b^2</tt>.
 */
double multSquare(double a, double b) => a * b * b;

/**
 * Function that returns <tt>a + b</tt>.
 */
final DoubleDoubleFunction plus = plusMultSecond(1.0);
//        new DoubleDoubleFunction() {
//            public final double apply(double a, double b) {
//                return a + b;
//            }
//        };

/**
 * Function that returns <tt>Math.abs(a) + Math.abs(b)</tt>.
 */
double plusAbs(double a, double b) => a.abs() + b.abs();

/**
 * Function that returns <tt>Math.pow(a,b)</tt>.
 */
double pow(double a, double b) => Math.pow(a, b);

/**
 * Makes this class non instantiable, but still let's others inherit from
 * it.
 */

/**
 * Constructs a function that returns <tt>(from<=a && a<=to) ? 1 : 0</tt>.
 * <tt>a</tt> is a variable, <tt>from</tt> and <tt>to</tt> are fixed.
 */
DoubleFunction between(double from, double to) {
    return (double a) {
            return (from <= a && a <= to) ? 1 : 0;
    };
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
DoubleFunction bindArg1(DoubleDoubleFunction function, double c) {
    return (double v) {
            return function(c, v);
    };
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
DoubleFunction bindArg2(DoubleDoubleFunction function, double c) {
    return (double v) {
            return function(v, c);
    };
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
DoubleDoubleFunction chainFGH(DoubleDoubleFunction f, DoubleFunction g,
        DoubleFunction h) {
    return (double a, double b) {
            return f(g(a), h(b));
    };
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
DoubleDoubleFunction chain(DoubleFunction g, DoubleDoubleFunction h) {
    return (double a, double b) {
            return g(h(a, b));
    };
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
DoubleFunction chainGH(DoubleFunction g, DoubleFunction h) {
    return (double a) {
            return g(h(a));
    };
}

/**
 * Constructs a function that returns <tt>a < b ? -1 : a > b ? 1 : 0</tt>.
 * <tt>a</tt> is a variable, <tt>b</tt> is fixed.
 */
DoubleFunction compareTo(double b) {
    return (double a) {
            return a < b ? -1 : a > b ? 1 : 0;
    };
}

/**
 * Constructs a function that returns the constant <tt>c</tt>.
 */
DoubleFunction constant(double c) {
    return (double a) {
            return c;
    };
}

/**
 * Demonstrates usage of this class.
 */
void demo1() {
    //DoubleFunctions F = cern.jet.math.tdouble.DoubleFunctions.functions;
    double a = 0.5;
    double b = 0.2;
    double v = Math.sin(a) + Math.pow(Math.cos(b), 2);
    print(v);
    DoubleDoubleFunction f = chainFGH(plus, sin, chainGH(square, cos));
    // DoubleDoubleFunction f = F.chain(plus,sin,F.chain(square,cos));
    print(f(a, b));
    double g(double x, double y) {
            return Math.sin(x) + Math.pow(Math.cos(y), 2);
    };
    print(g(a, b));
    DoubleFunction m = add(3.0);
    DoubleFunction n = add(4.0);
    print(m(0.0));
    print(n(0.0));
}

/**
 * Benchmarks and demonstrates usage of trivial and complex functions.
 */
void demo2(int size) {
    //cern.jet.math.tdouble.DoubleFunctions F = cern.jet.math.tdouble.DoubleFunctions.functions;
    print("\n\n");
    double a = 0.0;
    double b = 0.0;
    double v = (Math.sin(a) + Math.pow(Math.cos(b), 2)).abs();
    // double v = Math.sin(a) + Math.pow(Math.cos(b),2);
    // double v = a + b;
    print(v);

    // DoubleDoubleFunction f = F.chain(F.plus,F.identity,F.identity);
    DoubleDoubleFunction f = chain(abs, chainFGH(plus, sin, chainGH(square, cos)));
    // DoubleDoubleFunction f =
    // F.chain(F.plus,F.sin,F.chain(F.square,F.cos));
    // DoubleDoubleFunction f = F.plus;

    print(f(a, b));
    double g(double x, double y) {
            return (Math.sin(x) + Math.pow(Math.cos(y), 2)).abs();
    };
    print(g(a, b));

    // emptyLoop
    //colt.Timer emptyLoop = new colt.Timer().start();
    a = 0.0;
    b = 0.0;
    double sum = 0.0;
    for (int i = size; --i >= 0;) {
        sum += a;
        a++;
        b++;
    }
    //emptyLoop.stop().display();
    print("empty sum=$sum");

    //colt.Timer timer = new colt.Timer().start();
    a = 0.0;
    b = 0.0;
    sum = 0.0;
    for (int i = size; --i >= 0;) {
        sum += (Math.sin(a) + Math.pow(Math.cos(b), 2)).abs();
        // sum += a + b;
        a++;
        b++;
    }
    //timer.stop().display();
    //print("evals / sec = ${size / timer.minus(emptyLoop).seconds()}");
    print("sum=$sum");

    //timer.reset().start();
    a = 0.0;
    b = 0.0;
    sum = 0.0;
    for (int i = size; --i >= 0;) {
        sum += f(a, b);
        a++;
        b++;
    }
    //timer.stop().display();
    //print("evals / sec = ${size / timer.minus(emptyLoop).seconds()}");
    print("sum=${sum}");

    //timer.reset().start();
    a = 0.0;
    b = 0.0;
    sum = 0.0;
    for (int i = size; --i >= 0;) {
        sum += g(a, b);
        a++;
        b++;
    }
    //timer.stop().display();
    //print("evals / sec = ${size / timer.minus(emptyLoop).seconds()}");
    print("sum=$sum");

}

/**
 * Constructs a function that returns <tt>a / b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
DoubleFunction divBy(double b) {
    return multiply(1.0 / b);
}

/**
 * Constructs a function that returns <tt>a == b ? 1 : 0</tt>. <tt>a</tt> is
 * a variable, <tt>b</tt> is fixed.
 */
DoubleFunction equalTo(double b) {
    return (double a) {
            return a == b ? 1 : 0;
    };
}

/**
 * Constructs a function that returns <tt>a > b ? 1 : 0</tt>. <tt>a</tt> is
 * a variable, <tt>b</tt> is fixed.
 */
DoubleFunction greaterThan(double b) {
    return (double a) {
            return a > b ? 1 : 0;
    };
}

/**
 * Constructs a function that returns <tt>Math.IEEEremainder(a,b)</tt>.
 * <tt>a</tt> is a variable, <tt>b</tt> is fixed.
 */
/*DoubleFunction IEEEremainderOf(double b) {
    return (double a) {
            return Math.IEEEremainder(a, b);
    };
}*/

/**
 * Constructs a function that returns <tt>from<=a && a<=to</tt>. <tt>a</tt>
 * is a variable, <tt>from</tt> and <tt>to</tt> are fixed.
 */
DoubleProcedure isBetween(double from, double to) {
    return (double a) {
            return from <= a && a <= to;
    };
}

/**
 * Constructs a function that returns <tt>a == b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
DoubleProcedure isEqualTo(double b) {
    return (double a) {
            return a == b;
    };
}

/**
 * Constructs a function that returns <tt>a > b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
DoubleProcedure isGreaterThan(double b) {
    return (double a) {
            return a > b;
    };
}

/**
 * Constructs a function that returns <tt>a < b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
DoubleProcedure isLessThan(double b) {
    return (double a) {
            return a < b;
    };
}

/**
 * Constructs a function that returns <tt>a < b ? 1 : 0</tt>. <tt>a</tt> is
 * a variable, <tt>b</tt> is fixed.
 */
DoubleFunction lessThan(double b) {
    return (double a) {
            return a < b ? 1 : 0;
    };
}

/**
 * Constructs a function that returns <tt><tt>Math.log(a) / Math.log(b)</tt>
 * </tt>. <tt>a</tt> is a variable, <tt>b</tt> is fixed.
 */
DoubleFunction lgFn(double b) {
    double logInv = 1 / Math.log(b); // cached for speed
    return (double a) {
            return Math.log(a) * logInv;
    };
}

/**
 * Tests various methods of this class.
 */
void main() {
    int size = 6;//Integer.parseInt(args[0]);
    demo2(size);
    // demo1();
}

/**
 * Constructs a function that returns <tt>Math.max(a,b)</tt>. <tt>a</tt> is
 * a variable, <tt>b</tt> is fixed.
 */
DoubleFunction maximum(double b) {
    return (double a) {
            return Math.max(a, b);
    };
}

/**
 * Constructs a function that returns <tt>Math.min(a,b)</tt>. <tt>a</tt> is
 * a variable, <tt>b</tt> is fixed.
 */
DoubleFunction minimum(double b) {
    return (double a) {
            return Math.min(a, b);
    };
}

/**
 * Constructs a function that returns <tt>a - b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
DoubleFunction subtract(num b) {
    return add(-b.toDouble());
}

/**
 * Constructs a function that returns <tt>a - b*constant</tt>. <tt>a</tt>
 * and <tt>b</tt> are variables, <tt>constant</tt> is fixed.
 */
DoubleDoubleFunction minusMult(double constant) {
    return plusMultSecond(-constant);
}

/**
 * Constructs a function that returns <tt>a % b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
DoubleFunction modulus(double b) {
    return (double a) {
            return a % b;
    };
}

/**
 * Constructs a function that returns <tt>a * b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
DoubleFunction multiply(double b) {
    return new DoubleMult(b);
    /*
     * return new DoubleFunction() { public final double apply(double a) {
     * return a * b; } };
     */
}

/**
 * Constructs a function that returns <tt>a + b</tt>. <tt>a</tt> is a
 * variable, <tt>b</tt> is fixed.
 */
DoubleFunction add(double b) {
    return (double a) {
            return a + b;
    };
}

/**
 * Constructs a function that returns <tt>b*constant</tt>.
 */
DoubleDoubleFunction multSecond(double constant) {
    return (double a, double b) {
            return b * constant;
    };
}

/**
 * Constructs a function that returns <tt>a + b*constant</tt>. <tt>a</tt>
 * and <tt>b</tt> are variables, <tt>constant</tt> is fixed.
 */
DoubleDoubleFunction plusMultSecond(double constant) {
    return new DoublePlusMultSecond(constant);
}

/**
 * Constructs a function that returns <tt>a * constant + b</tt>. <tt>a</tt>
 * and <tt>b</tt> are variables, <tt>constant</tt> is fixed.
 */
DoubleDoubleFunction plusMultFirst(double constant) {
    return new DoublePlusMultFirst(constant);
}

/**
 * Constructs a function that returns <tt>Math.pow(a,b)</tt>. <tt>a</tt> is
 * a variable, <tt>b</tt> is fixed.
 */
DoubleFunction power(double b) {
    return (double a) {
            return Math.pow(a, b);
    };
}

/**
 * Constructs a function that returns a new uniform random number in the
 * open unit interval <code>(0.0,1.0)</code> (excluding 0.0 and 1.0).
 * Currently the engine is
 * {@link cern.jet.random.tdouble.engine.DoubleMersenneTwister} and is
 * seeded with the current time.
 * <p>
 * Note that any random engine derived from
 * {@link cern.jet.random.tdouble.engine.DoubleRandomEngine} and any random
 * distribution derived from
 * {@link cern.jet.random.tdouble.AbstractDoubleDistribution} are function
 * objects, because they implement the proper interfaces. Thus, if you are
 * not happy with the default, just pass your favourite random generator to
 * function evaluating methods.
 */
DoubleFunction random() {
    return RandomDoubleFunction;
}

// TODO
double RandomDoubleFunction(double argument) {
        return new Math.Random().nextDouble();
}

/**
 * Constructs a function that returns the number rounded to the given
 * precision; <tt>Math.rint(a/precision)*precision</tt>. Examples:
 *
 * <pre>
 * precision = 0.01 rounds 0.012 --&gt; 0.01, 0.018 --&gt; 0.02
 * precision = 10   rounds 123   --&gt; 120 , 127   --&gt; 130
 * </pre>
 */
DoubleFunction round(double precision) {
    return (double a) {
            return (a / precision).round() * precision;
    };
}

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
DoubleDoubleFunction swapArgs(final DoubleDoubleFunction function) {
    return (double a, double b) {
            return function(b, a);
    };
}