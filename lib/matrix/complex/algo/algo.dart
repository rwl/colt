library cern.colt.matrix.complex.algo;

import 'dart:typed_data';

import '../../matrix.dart';
import '../../former.dart';
import '../../../math.dart';

part 'property.dart';

/// Infinity norm of vector `x`, which is:
///     max(abs(x[i]))
double normInfinity(AbstractComplexVector x) {
  if (x.length == 0) {
    return 0.0;
  }
  var d = x.abs().real();
  return d.max().value;
}
