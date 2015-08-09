library cern.colt.matrix.complex.algo;

import 'matrix.dart';

/// Infinity norm of vector `x`, which is:
///     max(abs(x[i]))
double normInfinity(ComplexVector x) {
  if (x.size == 0) {
    return 0.0;
  }
  var d = x.abs().real();
  return d.max().value;
}
