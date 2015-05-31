library cern.colt.matrix.tdouble.algo;

import 'dart:math' as Math;

import '../../matrix.dart';
import '../../former.dart';
import '../../../math.dart';
import '../../../function/double.dart';

part 'property.dart';
//part 'formatter.dart';

/// Infinity norm of vector `x`, which is:
///     max(abs(x[i]))
double normInfinity(AbstractDoubleVector x) {
  if (x.length == 0) {
    return 0.0;
  }
  return x.aggregate(max, abs);
}

typedef AbstractDoubleVector Solve(AbstractDoubleMatrix, AbstractDoubleVector);
typedef dynamic Factor(AbstractDoubleMatrix);
typedef AbstractDoubleVector SolveFactors(dynamic, AbstractDoubleVector);