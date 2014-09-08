/*
Copyright (C) 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
is hereby granted without fee, provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear in supporting documentation.
CERN makes no representations about the suitability of this software for any purpose.
It is provided "as is" without expressed or implied warranty.
 */
library cern.jet.math;

import 'dart:math' as Math;
import 'dart:typed_data';

const MAX_INT = 2147483647; // 2^31-1
//const MAX_INT = 576460752303423487;// 2^59-1

const EPSILON = 1e-9;

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
class DoubleMult {
  /**
   * Public read/write access to avoid frequent object construction.
   */
  double multiplicator;

  /**
   * Insert the method's description here. Creation date: (8/10/99 19:12:09)
   */
  DoubleMult(final double multiplicator) {
    this.multiplicator = multiplicator;
  }

  /**
   * Returns the result of the function evaluation.
   */
  double call(double a) {
    return a * multiplicator;
  }

  /**
   * <tt>a / constant</tt>.
   */
  factory DoubleMult.div(final double constant) {
    return new DoubleMult.mult(1 / constant);
  }

  /**
   * <tt>a * constant</tt>.
   */
  factory DoubleMult.mult(final double constant) {
    return new DoubleMult(constant);
  }
}

/**
 * Only for performance tuning of compute intensive linear algebraic
 * computations. Constructs functions that return one of
 * <ul>
 * <li><tt>a*constant + b</tt>
 * <li><tt>a*constant - b</tt>
 * <li><tt>a/constant + b</tt>
 * <li><tt>a/constant - b</tt>
 * </ul>
 * <tt>a</tt> and <tt>b</tt> are variables, <tt>constant</tt> is fixed, but for
 * performance reasons publicly accessible. Intended to be passed to
 * <tt>matrix.assign(otherMatrix,function)</tt> methods.
 */
class DoublePlusMultFirst {
  /**
   * Public read/write access to avoid frequent object construction.
   */
  double multiplicator;

  DoublePlusMultFirst(final double multiplicator) {
    this.multiplicator = multiplicator;
  }

  /**
   * Returns the result of the function evaluation.
   */
  double call(double a, double b) {
    return a * multiplicator + b;
  }

  /**
   * <tt>a - b/constant</tt>.
   */
  factory DoublePlusMultFirst.minusDiv(final double constant) {
    return new DoublePlusMultFirst(-1 / constant);
  }

  /**
   * <tt>a - b*constant</tt>.
   */
  factory DoublePlusMultFirst.minusMult(final double constant) {
    return new DoublePlusMultFirst(-constant);
  }

  /**
   * <tt>a + b/constant</tt>.
   */
  factory DoublePlusMultFirst.plusDiv(final double constant) {
    return new DoublePlusMultFirst(1 / constant);
  }

  /**
   * <tt>a + b*constant</tt>.
   */
  factory DoublePlusMultFirst.plusMult(final double constant) {
    return new DoublePlusMultFirst(constant);
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
class DoublePlusMultSecond {//implements cern.colt.function.tdouble.DoubleDoubleFunction {
  /**
   * Public read/write access to avoid frequent object construction.
   */
  double multiplicator;

  /**
   * Insert the method's description here. Creation date: (8/10/99 19:12:09)
   */
  DoublePlusMultSecond(final double multiplicator) {
    this.multiplicator = multiplicator;
  }

  /**
   * Returns the result of the function evaluation.
   */
  double call(double a, double b) {
    return a + b * multiplicator;
  }

  /**
   * <tt>a - b/constant</tt>.
   */
  factory DoublePlusMultSecond.minusDiv(final double constant) {
    return new DoublePlusMultSecond(-1 / constant);
  }

  /**
   * <tt>a - b*constant</tt>.
   */
  factory DoublePlusMultSecond.minusMult(final double constant) {
    return new DoublePlusMultSecond(-constant);
  }

  /**
   * <tt>a + b/constant</tt>.
   */
  factory DoublePlusMultSecond.plusDiv(final double constant) {
    return new DoublePlusMultSecond(1 / constant);
  }

  /**
   * <tt>a + b*constant</tt>.
   */
  factory DoublePlusMultSecond.plusMult(final double constant) {
    return new DoublePlusMultSecond(constant);
  }
}

/**
 * Defines some useful constants.
 */
class DoubleConstants {
    /*
     * machine constants
     */
    static final double MACHEP = 1.11022302462515654042E-16;

    static final double MAXLOG = 7.09782712893383996732E2;

    static final double MINLOG = -7.451332191019412076235E2;

    static final double MAXGAM = 171.624376956302725;

    static final double SQTPI = 2.50662827463100050242E0;

    static final double SQRTH = 7.07106781186547524401E-1;

    static final double LOGPI = 1.14472988584940017414;

    static final double big = 4.503599627370496e15;

    static final double biginv = 2.22044604925031308085e-16;

}

/**
 * Complex arithmetic
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
class Complex extends DoubleConstants {

    static double abs(Float64List x) {
        double absX = x[0].abs();
        double absY = x[1].abs();
        if (absX == 0.0 && absY == 0.0) {
            return 0.0;
        } else if (absX >= absY) {
            double d = x[1] / x[0];
            return absX * Math.sqrt(1.0 + d * d);
        } else {
            double d = x[0] / x[1];
            return absY * Math.sqrt(1.0 + d * d);
        }
    }

    static double abs_(double re, double im) {
        double absX = re.abs();
        double absY = im.abs();
        if (absX == 0.0 && absY == 0.0) {
            return 0.0;
        } else if (absX >= absY) {
            double d = im / re;
            return absX * Math.sqrt(1.0 + d * d);
        } else {
            double d = re / im;
            return absY * Math.sqrt(1.0 + d * d);
        }
    }

    static Float64List acos(Float64List x) {
        Float64List z = new Float64List(2);

        double re, im;

        re = 1.0 - ((x[0] * x[0]) - (x[1] * x[1]));
        im = -((x[0] * x[1]) + (x[1] * x[0]));

        z[0] = re;
        z[1] = im;
        z = sqrt(z);

        re = -z[1];
        im = z[0];

        z[0] = x[0] + re;
        z[1] = x[1] + im;

        re = Math.log(abs(z));
        im = Math.atan2(z[1], z[0]);

        z[0] = im;
        z[1] = -re;
        return z;
    }

    static double arg(Float64List x) {
        return Math.atan2(x[1], x[0]);
    }

    static double arg_(double re, double im) {
        return Math.atan2(im, re);
    }

    static Float64List asin(Float64List x) {
        Float64List z = new Float64List(2);

        double re, im;

        re = 1.0 - ((x[0] * x[0]) - (x[1] * x[1]));
        im = -((x[0] * x[1]) + (x[1] * x[0]));

        z[0] = re;
        z[1] = im;
        z = sqrt(z);

        re = -z[1];
        im = z[0];

        z[0] = z[0] + re;
        z[1] = z[1] + im;

        re = Math.log(abs(z));
        im = Math.atan2(z[1], z[0]);

        z[0] = im;
        z[1] = -re;
        return z;
    }

    static Float64List atan(Float64List x) {
        Float64List z = new Float64List(2);

        double re, im;

        z[0] = -x[0];
        z[1] = 1.0 - x[1];

        re = x[0];
        im = 1.0 + x[1];

        z = div(z, re, im);

        re = Math.log(abs(z));
        im = Math.atan2(z[1], z[0]);

        z[0] = 0.5 * im;
        z[1] = -0.5 * re;

        return z;
    }

    static Float64List conj(Float64List x) {
        Float64List z = new Float64List(2);
        z[0] = x[0];
        z[1] = -x[1];
        return z;
    }

    static Float64List cos(Float64List x) {
        Float64List z = new Float64List(2);

        double re1, im1, re2, im2;
        double scalar;
        double iz_re, iz_im;

        iz_re = -x[1];
        iz_im = x[0];

        scalar = Math.exp(iz_re);
        re1 = scalar * Math.cos(iz_im);
        im1 = scalar * Math.sin(iz_im);

        scalar = Math.exp(-iz_re);
        re2 = scalar * Math.cos(-iz_im);
        im2 = scalar * Math.sin(-iz_im);

        re1 = re1 + re2;
        im1 = im1 + im2;

        z[0] = 0.5 * re1;
        z[1] = 0.5 * im1;

        return z;
    }

    static Float64List div(Float64List x, double re, double im) {
        Float64List z = new Float64List(2);
        double scalar;

        if (re.abs() >= im.abs()) {
            scalar = 1.0 / (re + im * (im / re));

            z[0] = scalar * (x[0] + x[1] * (im / re));
            z[1] = scalar * (x[1] - x[0] * (im / re));

        } else {
            scalar = 1.0 / (re * (re / im) + im);

            z[0] = scalar * (x[0] * (re / im) + x[1]);
            z[1] = scalar * (x[1] * (re / im) - x[0]);
        }

        return z;
    }

    static Float64List div_(Float64List x, Float64List y) {
        return div(x, y[0], y[1]);
    }

    static double equals(Float64List x, Float64List y, double tol) {
        if (abs_(x[0] - y[0], x[1] - y[1]) <= tol.abs()) {
            return 1.0;
        } else {
            return 0.0;
        }
    }

    static bool isEqual(Float64List x, Float64List y, double tol) {
        if (abs_(x[0] - y[0], x[1] - y[1]) <= tol.abs()) {
            return true;
        } else {
            return false;
        }
    }

    static Float64List exp(Float64List x) {
        Float64List z = new Float64List(2);
        double scalar = Math.exp(x[0]);
        z[0] = scalar * Math.cos(x[1]);
        z[1] = scalar * Math.sin(x[1]);
        return z;
    }

    static Float64List inv(Float64List x) {
        Float64List z = new Float64List(2);
        if (x[1] != 0.0) {
            double scalar;
            if (x[0].abs() >= x[1].abs()) {
                scalar = 1.0 / (x[0] + x[1] * (x[1] / x[0]));
                z[0] = scalar;
                z[1] = scalar * (-x[1] / x[0]);
            } else {
                scalar = 1.0 / (x[0] * (x[0] / x[1]) + x[1]);
                z[0] = scalar * (x[0] / x[1]);
                z[1] = -scalar;
            }
        } else {
            z[0] = 1 / x[0];
            z[1] = 0.0;
        }
        return z;
    }

    static Float64List log(Float64List x) {
        Float64List z = new Float64List(2);
        z[0] = Math.log(abs(x));
        z[1] = arg(x);
        return z;
    }

    static Float64List minus(Float64List x, Float64List y) {
        Float64List z = new Float64List(2);
        z[0] = x[0] - y[0];
        z[1] = x[1] - y[1];
        return z;
    }

    static Float64List minusAbs(Float64List x, Float64List y) {
        Float64List z = new Float64List(2);
        z[0] = (x[0] - y[0]).abs();
        z[1] = (x[1] - y[1]).abs();
        return z;
    }

    static Float64List mult(Float64List x, double y) {
        Float64List z = new Float64List(2);
        z[0] = x[0] * y;
        z[1] = x[1] * y;
        return z;
    }

    static Float64List multiply(Float64List x, Float64List y) {
        Float64List z = new Float64List(2);
        z[0] = x[0] * y[0] - x[1] * y[1];
        z[1] = x[1] * y[0] + x[0] * y[1];
        return z;
    }

    static Float64List neg(Float64List x) {
        Float64List z = new Float64List(2);
        z[0] = -x[0];
        z[1] = -x[1];
        return z;
    }

    static Float64List plus(Float64List x, Float64List y) {
        Float64List z = new Float64List(2);
        z[0] = x[0] + y[0];
        z[1] = x[1] + y[1];
        return z;
    }

    static Float64List pow(Float64List x, double y) {
        Float64List z = new Float64List(2);
        double re = y * Math.log(abs(x));
        double im = y * arg(x);
        double scalar = Math.exp(re);
        z[0] = scalar * Math.cos(im);
        z[1] = scalar * Math.sin(im);
        return z;
    }

    static Float64List power_(double x, Float64List y) {
        Float64List z = new Float64List(2);
        double re = Math.log(x.abs());
        double im = Math.atan2(0.0, x);

        double re2 = (re * y[0]) - (im * y[1]);
        double im2 = (re * y[1]) + (im * y[0]);

        double scalar = Math.exp(re2);

        z[0] = scalar * Math.cos(im2);
        z[1] = scalar * Math.sin(im2);
        return z;
    }

    static Float64List power(Float64List x, Float64List y) {
        Float64List z = new Float64List(2);
        double re = Math.log(abs(x));
        double im = arg(x);

        double re2 = (re * y[0]) - (im * y[1]);
        double im2 = (re * y[1]) + (im * y[0]);

        double scalar = Math.exp(re2);

        z[0] = scalar * Math.cos(im2);
        z[1] = scalar * Math.sin(im2);
        return z;
    }

    static Float64List sin(Float64List x) {
        Float64List z = new Float64List(2);
        double re1, im1, re2, im2;
        double scalar;
        double iz_re, iz_im;

        iz_re = -x[1];
        iz_im = x[0];

        scalar = Math.exp(iz_re);
        re1 = scalar * Math.cos(iz_im);
        im1 = scalar * Math.sin(iz_im);

        scalar = Math.exp(-iz_re);
        re2 = scalar * Math.cos(-iz_im);
        im2 = scalar * Math.sin(-iz_im);

        re1 = re1 - re2;
        im1 = im1 - im2;

        z[0] = 0.5 * im1;
        z[1] = -0.5 * re1;

        return z;
    }

    static Float64List sqrt(Float64List x) {
        Float64List z = new Float64List(2);
        double absx = abs(x);
        double tmp;
        if (absx > 0.0) {
            if (x[0] > 0.0) {
                tmp = Math.sqrt(0.5 * (absx + x[0]));
                z[0] = tmp;
                z[1] = 0.5 * (x[1] / tmp);
            } else {
                tmp = Math.sqrt(0.5 * (absx - x[0]));
                if (x[1] < 0.0) {
                    tmp = -tmp;
                }
                z[0] = 0.5 * (x[1] / tmp);
                z[1] = tmp;
            }
        } else {
            z[0] = 0.0;
            z[1] = 0.0;
        }
        return z;
    }

    static Float64List square(Float64List x) {
        return multiply(x, x);
    }

    static Float64List tan(Float64List x) {
        Float64List z = new Float64List(2);
        double scalar;
        double iz_re, iz_im;
        double re1, im1, re2, im2, re3, im3;
        double cs_re, cs_im;

        iz_re = -x[1];
        iz_im = x[0];

        scalar = Math.exp(iz_re);
        re1 = scalar * Math.cos(iz_im);
        im1 = scalar * Math.sin(iz_im);

        scalar = Math.exp(-iz_re);
        re2 = scalar * Math.cos(-iz_im);
        im2 = scalar * Math.sin(-iz_im);

        re3 = re1 - re2;
        im3 = im1 - im2;

        z[0] = 0.5 * im3;
        z[1] = -0.5 * re3;

        re3 = re1 + re2;
        im3 = im1 + im2;

        cs_re = 0.5 * re3;
        cs_im = 0.5 * im3;

        z = div(z, cs_re, cs_im);

        return z;
    }

    /**
     * Makes this class non instantiable, but still let's others inherit from
     * it.
     */
    Complex._internal();

}
