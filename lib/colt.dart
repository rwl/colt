library colt;

import 'dart:typed_data';

export 'matrix/matrix.dart';
export 'matrix/double/matrix.dart';
export 'matrix/int/matrix.dart';
export 'matrix/complex/matrix.dart';

//export 'matrix/double/algo/algo.dart';
//export 'matrix/complex/algo/algo.dart';

export 'matrix/double/algo/decomposition.dart';
export 'matrix/complex/algo/decomposition.dart';

//Int32List irange(int stop, {int start: 0, int step: 1}) {
//  var r = new Int32List(stop - start);
//  int v = start;
//  for (int i = 0; i < r.length; i++) {
//    r[i] = v;
//    v += step;
//  }
//  return r;
//}
//
//Int32List icat(Int32List a, Int32List b) {
//  var c = new Int32List(a.length + b.length);
//  c.setRange(0, a.length, a);
//  c.setRange(a.length, a.length + b.length, b);
//  return c;
//}
