library cern.colt.matrix.function.complex;

import 'dart:math' as Math;
import 'dart:typed_data';

import '../math.dart' as cmath;

typedef List<double> ComplexComplexComplexFunction(
    List<double> x, List<double> y);
typedef List<double> ComplexComplexFunc(double re, double im);
typedef List<double> ComplexComplexFunction(List<double> x);
//typedef bool ComplexComplexProcedure(List<double> x, List<double> y);
typedef double ComplexComplexRealFunction(List<double> x, List<double> y);
//typedef bool ComplexComplexRealProcedure(List<double> x, List<double> y, double tol);
typedef double ComplexComplexRealRealFunction(
    List<double> x, List<double> y, double tol);
//typedef bool ComplexProcedure(List<double> x);
typedef List<double> ComplexRealComplexFunction(List<double> x, double y);
typedef double ComplexRealFunction(List<double> x);
typedef List<double> IntIntComplexFunction(int x, int y, List<double> z);
typedef List<double> RealComplexComplexFunction(double x, List<double> y);
typedef List<double> RealComplexFunction(double x);

/// Unary functions

double abs(List<double> x) {
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

List<double> acos(List<double> x) {
  Float64List z = new Float64List(2);

  double re, im;

  re = (1.0 - ((x[0] * x[0]) - (x[1] * x[1])));
  im = -((x[0] * x[1]) + (x[1] * x[0]));

  z[0] = re;
  z[1] = im;
  z = cmath.sqrt(z);

  re = -z[1];
  im = z[0];

  z[0] = x[0] + re;
  z[1] = x[1] + im;

  re = Math.log(cmath.abs(z));
  im = Math.atan2(z[1], z[0]);

  z[0] = im;
  z[1] = -re;
  return z;
}

List<double> acos_(double re, double im) {
  Float64List z = new Float64List(2);

  double re2, im2;

  re2 = (1.0 - ((re * re) - (im * im)));
  im2 = -((re * im) + (im * re));

  z[0] = re2;
  z[1] = im2;
  z = cmath.sqrt(z);

  re2 = -z[1];
  im2 = z[0];

  z[0] = re + re2;
  z[1] = im + im2;

  re2 = Math.log(cmath.abs(z));
  im2 = Math.atan2(z[1], z[0]);

  z[0] = im2;
  z[1] = -re2;
  return z;
}

double arg(List<double> x) {
  return Math.atan2(x[1], x[0]);
}

List<double> asin(List<double> x) {
  Float64List z = new Float64List(2);

  double re, im;

  re = (1.0 - ((x[0] * x[0]) - (x[1] * x[1])));
  im = -((x[0] * x[1]) + (x[1] * x[0]));

  z[0] = re;
  z[1] = im;
  z = cmath.sqrt(z);

  re = -z[1];
  im = z[0];

  z[0] = z[0] + re;
  z[1] = z[1] + im;

  re = Math.log(cmath.abs(z));
  im = Math.atan2(z[1], z[0]);

  z[0] = im;
  z[1] = -re;
  return z;
}

List<double> asin_(double re, double im) {
  Float64List z = new Float64List(2);

  double re2, im2;

  re2 = (1.0 - ((re * re) - (im * im)));
  im2 = -((re * im) + (im * re));

  z[0] = re2;
  z[1] = im2;
  z = cmath.sqrt(z);

  re2 = -z[1];
  im2 = z[0];

  z[0] = z[0] + re2;
  z[1] = z[1] + im2;

  re2 = Math.log(cmath.abs(z));
  im2 = Math.atan2(z[1], z[0]);

  z[0] = im2;
  z[1] = -re2;
  return z;
}

List<double> atan(List<double> x) {
  Float64List z = new Float64List(2);

  double re, im;

  z[0] = -x[0];
  z[1] = 1.0 - x[1];

  re = x[0];
  im = 1.0 + x[1];

  z = cmath.div(z, re, im);

  re = Math.log(cmath.abs(z));
  im = Math.atan2(z[1], z[0]);

  z[0] = 0.5 * im;
  z[1] = -0.5 * re;

  return z;
}

List<double> atan_(double re, double im) {
  Float64List z = new Float64List(2);

  double re2, im2;

  z[0] = -re;
  z[1] = 1.0 - im;

  re2 = re;
  im2 = 1.0 + im;

  z = cmath.div(z, re2, im2);

  re2 = Math.log(cmath.abs(z));
  im2 = Math.atan2(z[1], z[0]);

  z[0] = 0.5 * im2;
  z[1] = -0.5 * re2;

  return z;
}

List<double> conj(List<double> x) {
  Float64List z = new Float64List(2);
  z[0] = x[0];
  z[1] = -x[1];
  return z;
}

List<double> conj_(double re, double im) {
  Float64List z = new Float64List(2);
  z[0] = re;
  z[1] = -im;
  return z;
}

List<double> cos(List<double> x) {
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

List<double> cos_(double re, double im) {
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

List<double> exp(List<double> x) {
  Float64List z = new Float64List(2);
  double scalar = Math.exp(x[0]);
  z[0] = (scalar * Math.cos(x[1]));
  z[1] = (scalar * Math.sin(x[1]));
  return z;
}

List<double> exp_(double re, double im) {
  Float64List z = new Float64List(2);
  double scalar = Math.exp(re);
  z[0] = (scalar * Math.cos(im));
  z[1] = (scalar * Math.sin(im));
  return z;
}

List<double> identity(List<double> x) {
  return x;
}

List<double> identity_(double re, double im) {
  return new Float64List.fromList([re, im]);
}

List<double> inv(List<double> x) {
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

List<double> inv_(double re, double im) {
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

List<double> log(List<double> x) {
  Float64List z = new Float64List(2);
  z[0] = Math.log(cmath.abs(x));
  z[1] = cmath.arg(x);
  return z;
}

List<double> log_(double re, double im) {
  Float64List z = new Float64List(2);
  z[0] = Math.log(cmath.abs_(re, im));
  z[1] = cmath.arg_(re, im);
  return z;
}

List<double> neg(List<double> x) {
  return new Float64List.fromList([-x[0], -x[1]]);
}

List<double> neg_(double re, double im) {
  return new Float64List.fromList([-re, -im]);
}

List<double> sin(List<double> x) {
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

List<double> sin_(double re, double im) {
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

List<double> sqrt(List<double> x) {
  Float64List z = new Float64List(2);
  double absx = cmath.abs(x);
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

List<double> sqrt_(double re, double im) {
  Float64List z = new Float64List(2);
  double absx = cmath.abs_(re, im);
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

List<double> square(List<double> x) {
  Float64List z = new Float64List(2);
  z[0] = x[0] * x[0] - x[1] * x[1];
  z[1] = x[1] * x[0] + x[0] * x[1];
  return z;
}

List<double> square_(double re, double im) {
  Float64List z = new Float64List(2);
  z[0] = re * re - im * im;
  z[1] = im * re + re * im;
  return z;
}

List<double> tan(List<double> x) {
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

  z = cmath.div(z, cs_re, cs_im);

  return z;
}

List<double> tan_(double re, double im) {
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

  z = cmath.div(z, cs_re, cs_im);

  return z;
}

/// Binary functions

List<double> div(List<double> x, List<double> y) {
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

double equals(List<double> x, List<double> y, double tol) {
  if (cmath.abs_(x[0] - y[0], x[1] - y[1]) <= tol.abs()) {
    return 1.0;
  } else {
    return 0.0;
  }
}

bool isEqual(List<double> x, List<double> y, double tol) {
  if (cmath.abs_(x[0] - y[0], x[1] - y[1]) <= tol.abs()) {
    return true;
  } else {
    return false;
  }
}

List<double> minus(List<double> x, List<double> y) {
  Float64List z = new Float64List(2);
  z[0] = x[0] - y[0];
  z[1] = x[1] - y[1];
  return z;
}

List<double> mult(List<double> x, List<double> y) {
  Float64List z = new Float64List(2);
  z[0] = x[0] * y[0] - x[1] * y[1];
  z[1] = x[1] * y[0] + x[0] * y[1];
  return z;
}

List<double> multConjFirst(List<double> x, List<double> y) {
  Float64List z = new Float64List(2);
  z[0] = x[0] * y[0] + x[1] * y[1];
  z[1] = -x[1] * y[0] + x[0] * y[1];
  return z;
}

List<double> multConjSecond(List<double> x, List<double> y) {
  Float64List z = new Float64List(2);
  z[0] = x[0] * y[0] + x[1] * y[1];
  z[1] = x[1] * y[0] - x[0] * y[1];
  return z;
}

List<double> plus(List<double> x, List<double> y) {
  Float64List z = new Float64List(2);
  z[0] = x[0] + y[0];
  z[1] = x[1] + y[1];
  return z;
}

List<double> pow1(List<double> x, double y) {
  Float64List z = new Float64List(2);
  double re = (y * Math.log(cmath.abs(x)));
  double im = y * cmath.arg(x);
  double scalar = Math.exp(re);
  z[0] = (scalar * Math.cos(im));
  z[1] = (scalar * Math.sin(im));
  return z;
}

List<double> pow2(double x, List<double> y) {
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

List<double> pow3(List<double> x, List<double> y) {
  Float64List z = new Float64List(2);
  double re = Math.log(cmath.abs(x));
  double im = cmath.arg(x);

  double re2 = (re * y[0]) - (im * y[1]);
  double im2 = (re * y[1]) + (im * y[0]);

  double scalar = Math.exp(re2);

  z[0] = (scalar * Math.cos(im2));
  z[1] = (scalar * Math.sin(im2));
  return z;
}

ComplexComplexFunction bindArg1(
    final ComplexComplexComplexFunction function, final List<double> c) {
  return (List<double> v) {
    return function(c, v);
  };
}

ComplexComplexFunc bindArg1_(
    final ComplexComplexComplexFunction function, final List<double> c) {
  return (double re, double im) {
    return function(c, new Float64List.fromList([re, im]));
  };
}

ComplexComplexFunction bindArg2(
    final ComplexComplexComplexFunction function, final List<double> c) {
  return (List<double> v) {
    return function(v, c);
  };
}

ComplexComplexFunc bindArg2_(
    final ComplexComplexComplexFunction function, final List<double> c) {
  return (double re, double im) {
    return function(new Float64List.fromList([re, im]), c);
  };
}

ComplexComplexComplexFunction chain(final ComplexComplexComplexFunction f,
    final ComplexComplexFunction g, final ComplexComplexFunction h) {
  return (List<double> x, List<double> y) {
    return f(g(x), h(y));
  };
}

ComplexComplexComplexFunction chain2(
    final ComplexComplexFunction g, final ComplexComplexComplexFunction h) {
  return (List<double> x, List<double> y) {
    return g(h(x, y));
  };
}

ComplexComplexFunction chain3(
    final ComplexComplexFunction g, final ComplexComplexFunction h) {
  return (List<double> x) {
    return g(h(x));
  };
}

ComplexComplexFunc chain3_(
    final ComplexComplexFunction g, final ComplexComplexFunction h) {
  return (double re, double im) {
    return g(h(new Float64List.fromList([re, im])));
  };
}

ComplexComplexFunction constant(final List<double> c) {
  return (List<double> x) {
    return c;
  };
}

ComplexComplexFunc constant_(final List<double> c) {
  return (double re, double im) {
    return new Float64List.fromList([re, im]);
  };
}

ComplexComplexFunction divide(final List<double> b) {
  return multiply(cmath.inv(b));
}

ComplexComplexFunction divideBy(final double b) {
  Float64List tmp = new Float64List.fromList([b, 0.0]);
  return multiply(cmath.inv(tmp));
}

ComplexRealFunction equalTo(final List<double> y) {
  return (List<double> x) {
    if (x[0] == y[0] && x[1] == y[1]) {
      return 1;
    } else {
      return 0;
    }
  };
}

/*ComplexProcedure isEqualTo(final List<double> y) {
  return (List<double> x) {
    if (x[0] == y[0] && x[1] == y[1]) {
      return true;
    } else {
      return false;
    }
  };
}*/

ComplexComplexFunction subtract(final List<double> x) {
  Float64List negb = new Float64List(2);
  negb[0] = -x[0];
  negb[1] = -x[1];
  return add(negb);
}

ComplexComplexComplexFunction minusMult(final List<double> constant) {
  Float64List negconstant = new Float64List(2);
  negconstant[0] = -constant[0];
  negconstant[1] = -constant[1];
  return plusMultSecond(negconstant);
}

ComplexComplexFunction multiply(final List<double> x) {
  return new ComplexMult(x);
}

ComplexComplexFunction scale(final double x) {
  return new ComplexMult(new Float64List.fromList([x, 0.0]));
}

ComplexComplexFunction add(final List<double> y) {
  return (List<double> x) {
    Float64List z = new Float64List(2);
    z[0] = x[0] + y[0];
    z[1] = x[1] + y[1];
    return z;
  };
}

ComplexComplexFunc add_(final List<double> y) {
  return (double re, double im) {
    Float64List z = new Float64List(2);
    z[0] = re + y[0];
    z[1] = im + y[1];
    return z;
  };
}

ComplexComplexComplexFunction plusMultSecond(List<double> constant) {
  return new ComplexPlusMultSecond(constant);
}

ComplexComplexComplexFunction plusMultFirst(List<double> constant) {
  return new ComplexPlusMultFirst(constant);
}

ComplexComplexFunction power1(final double y) {
  return (List<double> x) {
    Float64List z = new Float64List(2);
    double re = (y * Math.log(cmath.abs(x)));
    double im = y * cmath.arg(x);
    double scalar = Math.exp(re);
    z[0] = (scalar * Math.cos(im));
    z[1] = (scalar * Math.sin(im));
    return z;
  };
}

ComplexComplexFunc power1_(final double y) {
  return (double re, double im) {
    Float64List z = new Float64List(2);
    double re2 = (y * Math.log(cmath.abs_(re, im)));
    double im2 = y * cmath.arg_(re, im);
    double scalar = Math.exp(re2);
    z[0] = (scalar * Math.cos(im2));
    z[1] = (scalar * Math.sin(im2));
    return z;
  };
}

RealComplexFunction power2(final List<double> y) {
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

ComplexComplexFunction power3(final List<double> y) {
  return (List<double> x) {
    Float64List z = new Float64List(2);
    double re = Math.log(cmath.abs(x));
    double im = cmath.arg(x);

    double re2 = (re * y[0]) - (im * y[1]);
    double im2 = (re * y[1]) + (im * y[0]);

    double scalar = Math.exp(re2);

    z[0] = scalar * Math.cos(im2);
    z[1] = scalar * Math.sin(im2);
    return z;
  };
}

ComplexComplexFunc power3_(final List<double> y) {
  return (double re, double im) {
    Float64List z = new Float64List(2);
    double re1 = Math.log(cmath.abs_(re, im));
    double im1 = cmath.arg_(re, im);

    double re2 = (re1 * y[0]) - (im1 * y[1]);
    double im2 = (re1 * y[1]) + (im1 * y[0]);

    double scalar = Math.exp(re2);

    z[0] = scalar * Math.cos(im2);
    z[1] = scalar * Math.sin(im2);
    return z;
  };
}

final _r = new Math.Random();

List<double> random(List<double> argument) {
  return new Float64List.fromList([_r.nextDouble(), _r.nextDouble()]);
}

List<double> random_(double re, double im) {
  return new Float64List.fromList([_r.nextDouble(), _r.nextDouble()]);
}

ComplexComplexComplexFunction swapArgs(
    final ComplexComplexComplexFunction function) {
  return (List<double> x, List<double> y) {
    return function(y, x);
  };
}

/// Only for performance tuning of compute intensive linear algebraic
/// computations. Constructs functions that return one of:
/// - `a * constant`
/// - `a / constant`
/// `a` is variable, `constant` is fixed, but for performance
/// reasons publicly accessible. Intended to be passed to
/// `matrix.assign(function)` methods.
class ComplexMult {
  /// Public read/write access to avoid frequent object construction.
  Float64List multiplicator;

  ComplexMult(final List<double> multiplicator) {
    this.multiplicator = new Float64List.fromList(multiplicator);
  }

  /// Returns the result of the function evaluation.
  List<double> call(List<double> a) {
    Float64List z = new Float64List(2);
    z[0] = a[0] * multiplicator[0] - a[1] * multiplicator[1];
    z[1] = a[1] * multiplicator[0] + a[0] * multiplicator[1];
    return z;
  }

  /// Returns the result of the function evaluation.
  List<double> apply(double re, double im) {
    Float64List z = new Float64List(2);
    z[0] = re * multiplicator[0] - im * multiplicator[1];
    z[1] = im * multiplicator[0] + re * multiplicator[1];
    return z;
  }

  /// `a / constant`.
  static ComplexMult div(final List<double> constant) {
    return mult(cmath.inv(constant));
  }

  /// `a * constant`.
  static ComplexMult mult(final List<double> constant) {
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
  Float64List multiplicator;

  ComplexPlusMultFirst(final List<double> multiplicator) {
    this.multiplicator = new Float64List.fromList(multiplicator);
  }

  /// Returns the result of the function evaluation.
  List<double> call(List<double> a, List<double> b) {
    Float64List z = new Float64List(2);
    z[0] = a[0] * multiplicator[0] - a[1] * multiplicator[1];
    z[1] = a[1] * multiplicator[0] + a[0] * multiplicator[1];
    z[0] += b[0];
    z[1] += b[1];
    return z;
  }

  /// `a - b/constant`.
  static ComplexPlusMultFirst minusDiv(final List<double> constant) {
    return new ComplexPlusMultFirst(cmath.neg(cmath.inv(constant)));
  }

  /// `a - b*constant`.
  static ComplexPlusMultFirst minusMult(final List<double> constant) {
    return new ComplexPlusMultFirst(cmath.neg(constant));
  }

  /// `a + b/constant`.
  static ComplexPlusMultFirst plusDiv(final List<double> constant) {
    return new ComplexPlusMultFirst(cmath.inv(constant));
  }

  /// `a + b*constant`.
  static ComplexPlusMultFirst plusMult(final List<double> constant) {
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
  Float64List multiplicator;

  ComplexPlusMultSecond(final List<double> multiplicator) {
    this.multiplicator = new Float64List.fromList(multiplicator);
  }

  /// Returns the result of the function evaluation.
  List<double> call(List<double> a, List<double> b) {
    Float64List z = new Float64List(2);
    z[0] = b[0] * multiplicator[0] - b[1] * multiplicator[1];
    z[1] = b[1] * multiplicator[0] + b[0] * multiplicator[1];
    z[0] += a[0];
    z[1] += a[1];
    return z;
  }

  /// `a - b/constant`.
  static ComplexPlusMultSecond minusDiv(final List<double> constant) {
    return new ComplexPlusMultSecond(cmath.neg(cmath.inv(constant)));
  }

  /// `a - b*constant`.
  static ComplexPlusMultSecond minusMult(final List<double> constant) {
    return new ComplexPlusMultSecond(cmath.neg(constant));
  }

  /// `a + b/constant`.
  static ComplexPlusMultSecond plusDiv(final List<double> constant) {
    return new ComplexPlusMultSecond(cmath.inv(constant));
  }

  /// `a + b*constant`.
  static ComplexPlusMultSecond plusMult(final List<double> constant) {
    return new ComplexPlusMultSecond(constant);
  }
}
