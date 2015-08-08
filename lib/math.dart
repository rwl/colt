// Copyright (C) 1999 CERN - European Organization for Nuclear Research.
//
// Permission to use, copy, modify, distribute and sell this software and
// its documentation for any purpose is hereby granted without fee, provided
// that the above copyright notice appear in all copies and that both that
// copyright notice and this permission notice appear in supporting
// documentation.
//
// CERN makes no representations about the suitability of this software for
// any purpose. It is provided "as is" without expressed or implied warranty.
library cern.jet.math;

import 'dart:math' as Math;
import 'dart:typed_data';
import 'package:complex/complex.dart';

const MAX_INT = 2147483647; // 2^31-1
//const MAX_INT = 576460752303423487;// 2^59-1
const MIN_INT = -MAX_INT;

const EPSILON = 1e-9;

const DEG_RAD = Math.PI / 180.0;

const double MACHEP = 1.11022302462515654042E-16;

const double MAXLOG = 7.09782712893383996732E2;

const double MINLOG = -7.451332191019412076235E2;

const double MAXGAM = 171.624376956302725;

const double SQTPI = 2.50662827463100050242E0;

const double SQRTH = 7.07106781186547524401E-1;

const double LOGPI = 1.14472988584940017414;

const double big = 4.503599627370496e15;

const double biginv = 2.22044604925031308085e-16;

bool isEqual(Complex x, Complex y, double tol) {
  return (x - y).abs() <= tol.abs();
}

/*
double abs(List<double> x) {
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

double abs_(double re, double im) {
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

List<double> acos(List<double> x) {
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

double arg(List<double> x) {
  return Math.atan2(x[1], x[0]);
}

double arg_(double re, double im) {
  return Math.atan2(im, re);
}

List<double> asin(List<double> x) {
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

List<double> atan(List<double> x) {
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

List<double> conj(List<double> x) {
  Float64List z = new Float64List(2);
  z[0] = x[0];
  z[1] = -x[1];
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

List<double> div(List<double> x, double re, double im) {
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

List<double> div_(List<double> x, List<double> y) {
  return div(x, y[0], y[1]);
}

double equals(List<double> x, List<double> y, double tol) {
  if (abs_(x[0] - y[0], x[1] - y[1]) <= tol.abs()) {
    return 1.0;
  } else {
    return 0.0;
  }
}

List<double> exp(List<double> x) {
  Float64List z = new Float64List(2);
  double scalar = Math.exp(x[0]);
  z[0] = scalar * Math.cos(x[1]);
  z[1] = scalar * Math.sin(x[1]);
  return z;
}

List<double> inv(List<double> x) {
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

List<double> log(List<double> x) {
  Float64List z = new Float64List(2);
  z[0] = Math.log(abs(x));
  z[1] = arg(x);
  return z;
}

List<double> minus(List<double> x, List<double> y) {
  Float64List z = new Float64List(2);
  z[0] = x[0] - y[0];
  z[1] = x[1] - y[1];
  return z;
}

List<double> minusAbs(List<double> x, List<double> y) {
  Float64List z = new Float64List(2);
  z[0] = (x[0] - y[0]).abs();
  z[1] = (x[1] - y[1]).abs();
  return z;
}

List<double> mult(List<double> x, double y) {
  Float64List z = new Float64List(2);
  z[0] = x[0] * y;
  z[1] = x[1] * y;
  return z;
}

List<double> multiply(List<double> x, List<double> y) {
  Float64List z = new Float64List(2);
  z[0] = x[0] * y[0] - x[1] * y[1];
  z[1] = x[1] * y[0] + x[0] * y[1];
  return z;
}

List<double> neg(List<double> x) {
  Float64List z = new Float64List(2);
  z[0] = -x[0];
  z[1] = -x[1];
  return z;
}

List<double> plus(List<double> x, List<double> y) {
  Float64List z = new Float64List(2);
  z[0] = x[0] + y[0];
  z[1] = x[1] + y[1];
  return z;
}

List<double> pow(List<double> x, double y) {
  Float64List z = new Float64List(2);
  double re = y * Math.log(abs(x));
  double im = y * arg(x);
  double scalar = Math.exp(re);
  z[0] = scalar * Math.cos(im);
  z[1] = scalar * Math.sin(im);
  return z;
}

List<double> power_(double x, List<double> y) {
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

List<double> power(List<double> x, List<double> y) {
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

List<double> sin(List<double> x) {
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

List<double> sqrt(List<double> x) {
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

List<double> square(List<double> x) {
  return multiply(x, x);
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
*/