library cern.colt.matrix;

import 'dart:math' as Math;
import 'dart:async' show Future;
import 'dart:typed_data';

import 'package:csparse/csparse.dart';

import '../function/double.dart' show DoubleFunction, DoubleDoubleFunction, DoubleDoubleProcedure,
  DoubleIntProcedure, DoubleProcedure;
import '../function/double.dart' as func;
import '../function/complex.dart' as cfunc;
import '../math.dart';
import '../util.dart';
import 'double/algo/algo.dart';
import 'complex/algo/algo.dart';
import 'former.dart';

//part 'abstract_formatter.dart';
part 'abstract_matrix_1d.dart';
part 'abstract_matrix_2d.dart';
part 'abstract_matrix.dart';
part 'double/dense_double_matrix_1d.dart';
part 'double/dense_double_matrix_2d.dart';
part 'double/sparse_double_matrix_1d.dart';
part 'double/sparse_double_matrix_2d.dart';
part 'double/sparse_rc_double_matrix.dart';
part 'double/sparse_cc_double_matrix.dart';
part 'double/diagonal_double_matrix_2d.dart';
part 'double/wrapper_double_matrix_2d.dart';
part 'double/delegate_double_matrix.dart';
part 'double/double_matrix_1d.dart';
part 'double/double_matrix_2d.dart';
part 'double/factory.dart';
//part 'former.dart';
//part 'former_factory.dart';

part 'complex/complex_matrix_1d.dart';
part 'complex/complex_matrix_2d.dart';
part 'complex/sparse_complex_matrix_1d.dart';