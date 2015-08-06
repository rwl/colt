library cern.colt.matrix.tdouble.algo;

import 'matrix.dart';
import '../../function/double.dart';

/// Infinity norm of vector `x`, which is:
///     max(abs(x[i]))
double normInfinity(AbstractDoubleVector x) {
  if (x.size == 0) {
    return 0.0;
  }
  return x.aggregate(max, abs);
}
