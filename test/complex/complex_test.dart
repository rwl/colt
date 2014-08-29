library cern.colt.matrix.complex.test;

import 'dart:math' as math;
import 'dart:typed_data';

import 'package:unittest/unittest.dart';
import 'package:colt/colt.dart';
import 'package:colt/math.dart';
import 'package:colt/function/double.dart' as F;
import 'package:colt/function/complex.dart' hide equals;

part 'complex_matrix_1d_test.dart';
part 'complex_matrix_2d_test.dart';

final math.Random random = new math.Random(0);