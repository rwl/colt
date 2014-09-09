library cern.colt.matrix.function.complex;

import 'dart:math' as Math;
import 'dart:typed_data';

import '../math.dart';

typedef Float64List ComplexComplexComplexFunction(Float64List x, Float64List y);
typedef Float64List ComplexComplexFunc(double re, double im);
typedef Float64List ComplexComplexFunction(Float64List x);
typedef bool ComplexComplexProcedure(Float64List x, Float64List y);
typedef double ComplexComplexRealFunction(Float64List x, Float64List y);
typedef bool ComplexComplexRealProcedure(Float64List x, Float64List y, double tol);
typedef double ComplexComplexRealRealFunction(Float64List x, Float64List y, double tol);
typedef bool ComplexProcedure(Float64List x);
typedef Float64List ComplexRealComplexFunction(Float64List x, double y);
typedef double ComplexRealFunction(Float64List x);
typedef Float64List IntIntComplexFunction(int x, int y, Float64List z);
typedef Float64List RealComplexComplexFunction(double x, Float64List y);
typedef Float64List RealComplexFunction(double x);


/**
 * Complex function objects to be passed to generic methods.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */

/***************************************************************************
 * <H3>Unary functions</H3>
 **************************************************************************/

double abs(Float64List x) {
  double absX = x[0].abs();
  double absY = x[1].abs();
  if (absX == 0.0 && absY == 0.0) {
    return 0.0;
  } else if (absX >= absY) {
    double d = x[1] / x[0];
    return (absX * Math.sqrt(1.0 + d * d));
  } else {
    double d = x[0] / x[1];
    return (absY * Math.sqrt(1.0 + d * d));
  }
}

Float64List acos(Float64List x) {
  Float64List z = new Float64List(2);

  double re, im;

  re = (1.0 - ((x[0] * x[0]) - (x[1] * x[1])));
  im = -((x[0] * x[1]) + (x[1] * x[0]));

  z[0] = re;
  z[1] = im;
  z = Complex.sqrt(z);

  re = -z[1];
  im = z[0];

  z[0] = x[0] + re;
  z[1] = x[1] + im;

  re = Math.log(Complex.abs(z));
  im = Math.atan2(z[1], z[0]);

  z[0] = im;
  z[1] = -re;
  return z;
}

Float64List _acos(double re, double im) {
  Float64List z = new Float64List(2);

  double re2, im2;

  re2 = (1.0 - ((re * re) - (im * im)));
  im2 = -((re * im) + (im * re));

  z[0] = re2;
  z[1] = im2;
  z = Complex.sqrt(z);

  re2 = -z[1];
  im2 = z[0];

  z[0] = re + re2;
  z[1] = im + im2;

  re2 = Math.log(Complex.abs(z));
  im2 = Math.atan2(z[1], z[0]);

  z[0] = im2;
  z[1] = -re2;
  return z;
}

double arg(Float64List x) {
  return Math.atan2(x[1], x[0]);
}

Float64List asin(Float64List x) {
  Float64List z = new Float64List(2);

  double re, im;

  re = (1.0 - ((x[0] * x[0]) - (x[1] * x[1])));
  im = -((x[0] * x[1]) + (x[1] * x[0]));

  z[0] = re;
  z[1] = im;
  z = Complex.sqrt(z);

  re = -z[1];
  im = z[0];

  z[0] = z[0] + re;
  z[1] = z[1] + im;

  re = Math.log(Complex.abs(z));
  im = Math.atan2(z[1], z[0]);

  z[0] = im;
  z[1] = -re;
  return z;
}

Float64List _asin(double re, double im) {
  Float64List z = new Float64List(2);

  double re2, im2;

  re2 = (1.0 - ((re * re) - (im * im)));
  im2 = -((re * im) + (im * re));

  z[0] = re2;
  z[1] = im2;
  z = Complex.sqrt(z);

  re2 = -z[1];
  im2 = z[0];

  z[0] = z[0] + re2;
  z[1] = z[1] + im2;

  re2 = Math.log(Complex.abs(z));
  im2 = Math.atan2(z[1], z[0]);

  z[0] = im2;
  z[1] = -re2;
  return z;
}

Float64List atan(Float64List x) {
  Float64List z = new Float64List(2);

  double re, im;

  z[0] = -x[0];
  z[1] = 1.0 - x[1];

  re = x[0];
  im = 1.0 + x[1];

  z = Complex.div(z, re, im);

  re = Math.log(Complex.abs(z));
  im = Math.atan2(z[1], z[0]);

  z[0] = 0.5 * im;
  z[1] = -0.5 * re;

  return z;
}

Float64List _atan(double re, double im) {
  Float64List z = new Float64List(2);

  double re2, im2;

  z[0] = -re;
  z[1] = 1.0 - im;

  re2 = re;
  im2 = 1.0 + im;

  z = Complex.div(z, re2, im2);

  re2 = Math.log(Complex.abs(z));
  im2 = Math.atan2(z[1], z[0]);

  z[0] = 0.5 * im2;
  z[1] = -0.5 * re2;

  return z;
}

Float64List conj(Float64List x) {
  Float64List z = new Float64List(2);
  z[0] = x[0];
  z[1] = -x[1];
  return z;
}

Float64List _conj(double re, double im) {
  Float64List z = new Float64List(2);
  z[0] = re;
  z[1] = -im;
  return z;
}

Float64List cos(Float64List x) {
  Float64List z = new Float64List(2);

  double re1, im1, re2, im2;
  double scalar;
  double iz_re, iz_im;

  iz_re = -x[1];
  iz_im = x[0];

  scalar = Math.exp(iz_re);
  re1 = (scalar * Math.cos(iz_im));
  im1 = (scalar * Math.sin(iz_im));

  scalar = Math.exp(-iz_re);
  re2 = (scalar * Math.cos(-iz_im));
  im2 = (scalar * Math.sin(-iz_im));

  re1 = re1 + re2;
  im1 = im1 + im2;

  z[0] = 0.5 * re1;
  z[1] = 0.5 * im1;

  return z;
}

Float64List _cos(double re, double im) {
  Float64List z = new Float64List(2);

  double re1, im1, re2, im2;
  double scalar;
  double iz_re, iz_im;

  iz_re = -im;
  iz_im = re;

  scalar = Math.exp(iz_re);
  re1 = (scalar * Math.cos(iz_im));
  im1 = (scalar * Math.sin(iz_im));

  scalar = Math.exp(-iz_re);
  re2 = (scalar * Math.cos(-iz_im));
  im2 = (scalar * Math.sin(-iz_im));

  re1 = re1 + re2;
  im1 = im1 + im2;

  z[0] = 0.5 * re1;
  z[1] = 0.5 * im1;

  return z;
}

Float64List exp(Float64List x) {
  Float64List z = new Float64List(2);
  double scalar = Math.exp(x[0]);
  z[0] = (scalar * Math.cos(x[1]));
  z[1] = (scalar * Math.sin(x[1]));
  return z;
}

Float64List _exp(double re, double im) {
  Float64List z = new Float64List(2);
  double scalar = Math.exp(re);
  z[0] = (scalar * Math.cos(im));
  z[1] = (scalar * Math.sin(im));
  return z;
}

Float64List identity(Float64List x) {
  return x;
}

Float64List _identity(double re, double im) {
  return new Float64List.fromList([re, im]);
}

Float64List inv(Float64List x) {
  Float64List z = new Float64List(2);
  if (x[1] != 0.0) {
    double tmp = (x[0] * x[0]) + (x[1] * x[1]);
    z[0] = x[0] / tmp;
    z[1] = -x[1] / tmp;
  } else {
    z[0] = 1 / x[0];
    z[1] = 0.0;
  }
  return z;
}

Float64List _inv(double re, double im) {
  Float64List z = new Float64List(2);
  if (im != 0.0) {
    double scalar;
    if (re.abs() >= z[1].abs()) {
      scalar = (1.0 / (re + im * (im / re)));
      z[0] = scalar;
      z[1] = scalar * (-im / re);
    } else {
      scalar = (1.0 / (re * (re / im) + im));
      z[0] = scalar * (re / im);
      z[1] = -scalar;
    }
  } else {
    z[0] = 1 / re;
    z[1] = 0.0;
  }
  return z;
}

Float64List log(Float64List x) {
  Float64List z = new Float64List(2);
  z[0] = Math.log(Complex.abs(x));
  z[1] = Complex.arg(x);
  return z;
}

Float64List _log(double re, double im) {
  Float64List z = new Float64List(2);
  z[0] = Math.log(Complex.abs_(re, im));
  z[1] = Complex.arg_(re, im);
  return z;
}

Float64List neg(Float64List x) {
  return new Float64List.fromList([-x[0], -x[1]]);
}

Float64List _neg(double re, double im) {
  return new Float64List.fromList([-re, -im]);
}

Float64List sin(Float64List x) {
  Float64List z = new Float64List(2);
  double re1, im1, re2, im2;
  double scalar;
  double iz_re, iz_im;

  iz_re = -x[1];
  iz_im = x[0];

  scalar = Math.exp(iz_re);
  re1 = (scalar * Math.cos(iz_im));
  im1 = (scalar * Math.sin(iz_im));

  scalar = Math.exp(-iz_re);
  re2 = (scalar * Math.cos(-iz_im));
  im2 = (scalar * Math.sin(-iz_im));

  re1 = re1 - re2;
  im1 = im1 - im2;

  z[0] = 0.5 * im1;
  z[1] = -0.5 * re1;

  return z;
}

Float64List _sin(double re, double im) {
  Float64List z = new Float64List(2);
  double re1, im1, re2, im2;
  double scalar;
  double iz_re, iz_im;

  iz_re = -im;
  iz_im = re;

  scalar = Math.exp(iz_re);
  re1 = (scalar * Math.cos(iz_im));
  im1 = (scalar * Math.sin(iz_im));

  scalar = Math.exp(-iz_re);
  re2 = (scalar * Math.cos(-iz_im));
  im2 = (scalar * Math.sin(-iz_im));

  re1 = re1 - re2;
  im1 = im1 - im2;

  z[0] = 0.5 * im1;
  z[1] = -0.5 * re1;

  return z;
}

Float64List sqrt(Float64List x) {
  Float64List z = new Float64List(2);
  double absx = Complex.abs(x);
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

Float64List _sqrt(double re, double im) {
  Float64List z = new Float64List(2);
  double absx = Complex.abs_(re, im);
  double tmp;
  if (absx > 0.0) {
    if (re > 0.0) {
      tmp = Math.sqrt(0.5 * (absx + re));
      z[0] = tmp;
      z[1] = 0.5 * (im / tmp);
    } else {
      tmp = Math.sqrt(0.5 * (absx - re));
      if (im < 0.0) {
        tmp = -tmp;
      }
      z[0] = 0.5 * (im / tmp);
      z[1] = tmp;
    }
  } else {
    z[0] = 0.0;
    z[1] = 0.0;
  }
  return z;
}

Float64List square(Float64List x) {
  Float64List z = new Float64List(2);
  z[0] = x[0] * x[0] - x[1] * x[1];
  z[1] = x[1] * x[0] + x[0] * x[1];
  return z;
}

Float64List _square(double re, double im) {
  Float64List z = new Float64List(2);
  z[0] = re * re - im * im;
  z[1] = im * re + re * im;
  return z;
}

Float64List tan(Float64List x) {
  Float64List z = new Float64List(2);
  double scalar;
  double iz_re, iz_im;
  double re1, im1, re2, im2, re3, im3;
  double cs_re, cs_im;

  iz_re = -x[1];
  iz_im = x[0];

  scalar = Math.exp(iz_re);
  re1 = (scalar * Math.cos(iz_im));
  im1 = (scalar * Math.sin(iz_im));

  scalar = Math.exp(-iz_re);
  re2 = (scalar * Math.cos(-iz_im));
  im2 = (scalar * Math.sin(-iz_im));

  re3 = re1 - re2;
  im3 = im1 - im2;

  z[0] = 0.5 * im3;
  z[1] = -0.5 * re3;

  re3 = re1 + re2;
  im3 = im1 + im2;

  cs_re = 0.5 * re3;
  cs_im = 0.5 * im3;

  z = Complex.div(z, cs_re, cs_im);

  return z;
}

Float64List _tan(double re, double im) {
  Float64List z = new Float64List(2);
  double scalar;
  double iz_re, iz_im;
  double re1, im1, re2, im2, re3, im3;
  double cs_re, cs_im;

  iz_re = -im;
  iz_im = re;

  scalar = Math.exp(iz_re);
  re1 = (scalar * Math.cos(iz_im));
  im1 = (scalar * Math.sin(iz_im));

  scalar = Math.exp(-iz_re);
  re2 = (scalar * Math.cos(-iz_im));
  im2 = (scalar * Math.sin(-iz_im));

  re3 = re1 - re2;
  im3 = im1 - im2;

  z[0] = 0.5 * im3;
  z[1] = -0.5 * re3;

  re3 = re1 + re2;
  im3 = im1 + im2;

  cs_re = 0.5 * re3;
  cs_im = 0.5 * im3;

  z = Complex.div(z, cs_re, cs_im);

  return z;
}

/***************************************************************************
 * <H3>Binary functions</H3>
 **************************************************************************/

Float64List div(Float64List x, Float64List y) {
  double re = y[0];
  double im = y[1];

  Float64List z = new Float64List(2);
  double scalar;

  if (re.abs() >= im.abs()) {
    scalar = (1.0 / (re + im * (im / re)));

    z[0] = scalar * (x[0] + x[1] * (im / re));
    z[1] = scalar * (x[1] - x[0] * (im / re));

  } else {
    scalar = (1.0 / (re * (re / im) + im));

    z[0] = scalar * (x[0] * (re / im) + x[1]);
    z[1] = scalar * (x[1] * (re / im) - x[0]);
  }

  return z;
}

double equals(Float64List x, Float64List y, double tol) {
  if (Complex.abs_(x[0] - y[0], x[1] - y[1]) <= tol.abs()) {
    return 1.0;
  } else {
    return 0.0;
  }
}

bool isEqual(Float64List x, Float64List y, double tol) {
  if (Complex.abs_(x[0] - y[0], x[1] - y[1]) <= tol.abs()) {
    return true;
  } else {
    return false;
  }
}

Float64List minus(Float64List x, Float64List y) {
  Float64List z = new Float64List(2);
  z[0] = x[0] - y[0];
  z[1] = x[1] - y[1];
  return z;
}

Float64List mult(Float64List x, Float64List y) {
  Float64List z = new Float64List(2);
  z[0] = x[0] * y[0] - x[1] * y[1];
  z[1] = x[1] * y[0] + x[0] * y[1];
  return z;
}

Float64List multConjFirst(Float64List x, Float64List y) {
  Float64List z = new Float64List(2);
  z[0] = x[0] * y[0] + x[1] * y[1];
  z[1] = -x[1] * y[0] + x[0] * y[1];
  return z;
}

Float64List multConjSecond(Float64List x, Float64List y) {
  Float64List z = new Float64List(2);
  z[0] = x[0] * y[0] + x[1] * y[1];
  z[1] = x[1] * y[0] - x[0] * y[1];
  return z;
}

Float64List plus(Float64List x, Float64List y) {
  Float64List z = new Float64List(2);
  z[0] = x[0] + y[0];
  z[1] = x[1] + y[1];
  return z;
}

Float64List pow1(Float64List x, double y) {
  Float64List z = new Float64List(2);
  double re = (y * Math.log(Complex.abs(x)));
  double im = y * Complex.arg(x);
  double scalar = Math.exp(re);
  z[0] = (scalar * Math.cos(im));
  z[1] = (scalar * Math.sin(im));
  return z;
}

Float64List pow2(double x, Float64List y) {
  Float64List z = new Float64List(2);
  double re = Math.log(x.abs());
  double im = Math.atan2(0.0, x);

  double re2 = (re * y[0]) - (im * y[1]);
  double im2 = (re * y[1]) + (im * y[0]);

  double scalar = Math.exp(re2);

  z[0] = (scalar * Math.cos(im2));
  z[1] = (scalar * Math.sin(im2));
  return z;
}

Float64List pow3(Float64List x, Float64List y) {
  Float64List z = new Float64List(2);
  double re = Math.log(Complex.abs(x));
  double im = Complex.arg(x);

  double re2 = (re * y[0]) - (im * y[1]);
  double im2 = (re * y[1]) + (im * y[0]);

  double scalar = Math.exp(re2);

  z[0] = (scalar * Math.cos(im2));
  z[1] = (scalar * Math.sin(im2));
  return z;
}

ComplexComplexFunction bindArg1(final ComplexComplexComplexFunction function, final Float64List c) {
  return (Float64List v) {
    return function(c, v);
  };
}

//        final Float64List apply(double re, double im) {
//            return function.apply(c, new Float64List { re, im });
//        }
//    };
//}

ComplexComplexFunction bindArg2(final ComplexComplexComplexFunction function, final Float64List c) {
  return (Float64List v) {
    return function(v, c);
  };
}

//        final Float64List apply(double re, double im) {
//            return function.apply(new Float64List { re, im }, c);
//        }
//    };
//}

ComplexComplexComplexFunction chain(final ComplexComplexComplexFunction f, final ComplexComplexFunction g, final ComplexComplexFunction h) {
  return (Float64List x, Float64List y) {
    return f(g(x), h(y));
  };
}

ComplexComplexComplexFunction chain2(final ComplexComplexFunction g, final ComplexComplexComplexFunction h) {
  return (Float64List x, Float64List y) {
    return g(h(x, y));
  };
}

ComplexComplexFunction chain3(final ComplexComplexFunction g, final ComplexComplexFunction h) {
  return (Float64List x) {
    return g(h(x));
  };

  //        final Float64List apply(double re, double im) {
  //            return g.apply(h.apply(new Float64List { re, im }));
  //        }
  //    };
}

ComplexComplexFunction constant(final Float64List c) {
  return (Float64List x) {
    return c;
  };

  //        final Float64List apply(double re, double im) {
  //            return new Float64List { re, im };
  //        }
  //    };
}

ComplexComplexFunction divide(final Float64List b) {
  return multiply(Complex.inv(b));
}

ComplexComplexFunction divideBy(final double b) {
  Float64List tmp = new Float64List.fromList([b, 0]);
  return multiply(Complex.inv(tmp));
}

ComplexRealFunction equalTo(final Float64List y) {
  return (Float64List x) {
    if (x[0] == y[0] && x[1] == y[1]) {
      return 1;
    } else {
      return 0;
    }
  };
}

ComplexProcedure isEqualTo(final Float64List y) {
  return (Float64List x) {
    if (x[0] == y[0] && x[1] == y[1]) {
      return true;
    } else {
      return false;
    }
  };
}

ComplexComplexFunction subtract(final Float64List x) {
  Float64List negb = new Float64List(2);
  negb[0] = -x[0];
  negb[1] = -x[1];
  return add(negb);
}

ComplexComplexComplexFunction minusMult(final Float64List constant) {
  Float64List negconstant = new Float64List(2);
  negconstant[0] = -constant[0];
  negconstant[1] = -constant[1];
  return plusMultSecond(negconstant);
}

ComplexComplexFunction multiply(final Float64List x) {
  return new ComplexMult(x);
}

ComplexComplexFunction scale(final double x) {
  return new ComplexMult(new Float64List.fromList([x, 0.0]));
}

ComplexComplexFunction add(final Float64List y) {
  return (Float64List x) {
    Float64List z = new Float64List(2);
    z[0] = x[0] + y[0];
    z[1] = x[1] + y[1];
    return z;
  };

  //        final Float64List apply(double re, double im) {
  //            Float64List z = new Float64List(2);
  //            z[0] = re + y[0];
  //            z[1] = im + y[1];
  //            return z;
  //        }
  //    };
}

ComplexComplexComplexFunction plusMultSecond(Float64List constant) {
  return new ComplexPlusMultSecond(constant);
}

ComplexComplexComplexFunction plusMultFirst(Float64List constant) {
  return new ComplexPlusMultFirst(constant);
}

ComplexComplexFunction power1(final double y) {
  return (Float64List x) {
    Float64List z = new Float64List(2);
    double re = (y * Math.log(Complex.abs(x)));
    double im = y * Complex.arg(x);
    double scalar = Math.exp(re);
    z[0] = (scalar * Math.cos(im));
    z[1] = (scalar * Math.sin(im));
    return z;
  };

  //        final Float64List apply(double re, double im) {
  //            Float64List z = new Float64List(2);
  //            double re2 = (y * Math.log(Complex.abs(re, im)));
  //            double im2 = y * Complex.arg(re, im);
  //            double scalar = Math.exp(re2);
  //            z[0] = (scalar * Math.cos(im2));
  //            z[1] = (scalar * Math.sin(im2));
  //            return z;
  //        }
  //    };
}

RealComplexFunction power2(final Float64List y) {
  return (double x) {
    Float64List z = new Float64List(2);
    double re = Math.log(x.abs());
    double im = Math.atan2(0.0, x);

    double re2 = (re * y[0]) - (im * y[1]);
    double im2 = (re * y[1]) + (im * y[0]);

    double scalar = Math.exp(re2);

    z[0] = (scalar * Math.cos(im2));
    z[1] = (scalar * Math.sin(im2));
    return z;
  };
}

ComplexComplexFunction power3(final Float64List y) {
  return (Float64List x) {
    Float64List z = new Float64List(2);
    double re = Math.log(Complex.abs(x));
    double im = Complex.arg(x);

    double re2 = (re * y[0]) - (im * y[1]);
    double im2 = (re * y[1]) + (im * y[0]);

    double scalar = Math.exp(re2);

    z[0] = scalar * Math.cos(im2);
    z[1] = scalar * Math.sin(im2);
    return z;
  };

  //        final Float64List apply(double re, double im) {
  //            Float64List z = new Float64List(2);
  //            double re1 = (double) Math.log(Complex.abs(re, im));
  //            double im1 = Complex.arg(re, im);
  //
  //            double re2 = (re1 * y[0]) - (im1 * y[1]);
  //            double im2 = (re1 * y[1]) + (im1 * y[0]);
  //
  //            double scalar = (double) Math.exp(re2);
  //
  //            z[0] = (double) (scalar * Math.cos(im2));
  //            z[1] = (double) (scalar * Math.sin(im2));
  //            return z;
  //        }
  //    };
}

Float64List random(Float64List argument) {
  final r = new Math.Random();
  return new Float64List.fromList([r.nextDouble(), r.nextDouble()]);
}

Float64List _random(double re, double im) {
  final r = new Math.Random();
  return new Float64List.fromList([r.nextDouble(), r.nextDouble()]);
}

ComplexComplexComplexFunction swapArgs(final ComplexComplexComplexFunction function) {
  return (Float64List x, Float64List y) {
    return function(y, x);
  };
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
class ComplexMult {//implements cern.colt.function.tdcomplex.ComplexComplexFunction {
  /**
   * Public read/write access to avoid frequent object construction.
   */
  Float64List multiplicator;

  ComplexMult(final Float64List multiplicator) {
    this.multiplicator = multiplicator;
  }

  /**
   * Returns the result of the function evaluation.
   */
  Float64List call(Float64List a) {
    Float64List z = new Float64List(2);
    z[0] = a[0] * multiplicator[0] - a[1] * multiplicator[1];
    z[1] = a[1] * multiplicator[0] + a[0] * multiplicator[1];
    return z;
  }

  /**
   * Returns the result of the function evaluation.
   */
  Float64List apply(double re, double im) {
    Float64List z = new Float64List(2);
    z[0] = re * multiplicator[0] - im * multiplicator[1];
    z[1] = im * multiplicator[0] + re * multiplicator[1];
    return z;
  }

  /**
   * <tt>a / constant</tt>.
   */
  static ComplexMult div(final Float64List constant) {
    return mult(Complex.inv(constant));
  }

  /**
   * <tt>a * constant</tt>.
   */
  static ComplexMult mult(final Float64List constant) {
    return new ComplexMult(constant);
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
class ComplexPlusMultFirst {//implements cern.colt.function.tdcomplex.ComplexComplexComplexFunction {
  /**
   * Public read/write access to avoid frequent object construction.
   */
  Float64List multiplicator;

  /**
   * Insert the method's description here. Creation date: (8/10/99 19:12:09)
   */
  ComplexPlusMultFirst(final Float64List multiplicator) {
    this.multiplicator = multiplicator;
  }

  /**
   * Returns the result of the function evaluation.
   */
  Float64List call(Float64List a, Float64List b) {
    Float64List z = new Float64List(2);
    z[0] = a[0] * multiplicator[0] - a[1] * multiplicator[1];
    z[1] = a[1] * multiplicator[0] + a[0] * multiplicator[1];
    z[0] += b[0];
    z[1] += b[1];
    return z;
  }

  /**
   * <tt>a - b/constant</tt>.
   */
  static ComplexPlusMultFirst minusDiv(final Float64List constant) {
    return new ComplexPlusMultFirst(Complex.neg(Complex.inv(constant)));
  }

  /**
   * <tt>a - b*constant</tt>.
   */
  static ComplexPlusMultFirst minusMult(final Float64List constant) {
    return new ComplexPlusMultFirst(Complex.neg(constant));
  }

  /**
   * <tt>a + b/constant</tt>.
   */
  static ComplexPlusMultFirst plusDiv(final Float64List constant) {
    return new ComplexPlusMultFirst(Complex.inv(constant));
  }

  /**
   * <tt>a + b*constant</tt>.
   */
  static ComplexPlusMultFirst plusMult(final Float64List constant) {
    return new ComplexPlusMultFirst(constant);
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
class ComplexPlusMultSecond {//implements cern.colt.function.tdcomplex.ComplexComplexComplexFunction {
  /**
   * Public read/write access to avoid frequent object construction.
   */
  Float64List multiplicator;

  /**
   * Insert the method's description here. Creation date: (8/10/99 19:12:09)
   */
  ComplexPlusMultSecond(final Float64List multiplicator) {
    this.multiplicator = multiplicator;
  }

  /**
   * Returns the result of the function evaluation.
   */
  Float64List call(Float64List a, Float64List b) {
    Float64List z = new Float64List(2);
    z[0] = b[0] * multiplicator[0] - b[1] * multiplicator[1];
    z[1] = b[1] * multiplicator[0] + b[0] * multiplicator[1];
    z[0] += a[0];
    z[1] += a[1];
    return z;
  }

  /**
   * <tt>a - b/constant</tt>.
   */
  static ComplexPlusMultSecond minusDiv(final Float64List constant) {
    return new ComplexPlusMultSecond(Complex.neg(Complex.inv(constant)));
  }

  /**
   * <tt>a - b*constant</tt>.
   */
  static ComplexPlusMultSecond minusMult(final Float64List constant) {
    return new ComplexPlusMultSecond(Complex.neg(constant));
  }

  /**
   * <tt>a + b/constant</tt>.
   */
  static ComplexPlusMultSecond plusDiv(final Float64List constant) {
    return new ComplexPlusMultSecond(Complex.inv(constant));
  }

  /**
   * <tt>a + b*constant</tt>.
   */
  static ComplexPlusMultSecond plusMult(final Float64List constant) {
    return new ComplexPlusMultSecond(constant);
  }
}
