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

final EPS = Math.pow(2, -52);

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

int searchRange(Int32List list, int key, int from, int to) {
  while (from <= to) {
    if (list[from] == key) {
      return from;
    } else {
      from++;
      continue;
    }
  }
  return -(from + 1); // not found
}

int cumsum(Int32List p, Int32List c) {
  int nz = 0;
  for (int k = 0; k < c.length; k++) {
    p[k] = nz;
    nz += c[k];
    c[k] = p[k];
  }
  p[c.length] = nz;
  return nz;
}

int toDense(Int32List columnPointersA, Int32List rowIndexesA, List<num> valuesA,
    int j, num beta, Int32List w, Float64List x, int mark,
    Int32List rowIndexesC, int nz) {
  for (var p = columnPointersA[j]; p < columnPointersA[j + 1]; p++) {
    var i = rowIndexesA[p];
    if (w[i] < mark) {
      w[i] = mark;
      rowIndexesC[nz++] = i;
      if (x != null) {
        x[i] = beta * valuesA[p];
      }
    } else if (x != null) {
      x[i] += beta * valuesA[p];
    }
  }
  return nz;
}
