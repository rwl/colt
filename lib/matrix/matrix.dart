library cern.colt.matrix;

import 'dart:math' as Math;
import 'dart:async' show Future;
import 'dart:typed_data';
import 'dart:collection' show ListMixin;

import 'package:intl/intl.dart';

import 'package:csparse/double/csparse.dart';
import 'package:csparse/complex/cxsparse.dart' show DZcs;
import 'package:csparse/complex/cxsparse.dart' as cxsparse;

import '../function/double.dart' show DoubleFunction, DoubleDoubleFunction, DoubleDoubleProcedure,
  DoubleIntProcedure, DoubleProcedure;
import '../function/double.dart' as func;
import '../function/complex.dart' as cfunc;
import '../function/int.dart' as ifunc;
import '../math.dart';
//import '../util.dart';
import 'double/algo/algo.dart' as dprop;
import 'complex/algo/algo.dart' as cprop;
import 'int/algo/algo.dart' as iprop;
import 'former.dart';

//part 'abstract_formatter.dart';
part 'abstract_vector.dart';
part 'abstract_matrix.dart';

part 'double/double_vector.dart';
part 'double/double_matrix.dart';
part 'double/large_double_matrix.dart';
part 'double/sparse_double_vector.dart';
part 'double/sparse_double_matrix.dart';
part 'double/sparse_rc_double_matrix.dart';
part 'double/sparse_cc_double_matrix.dart';
part 'double/diagonal_double_matrix.dart';
part 'double/wrapper_double_matrix.dart';
part 'double/delegate_double_vector.dart';
part 'double/abstract_double_vector.dart';
part 'double/abstract_double_matrix.dart';
//part 'former.dart';
//part 'former_factory.dart';

part 'complex/abstract_complex_vector.dart';
part 'complex/abstract_complex_matrix.dart';
part 'complex/sparse_complex_vector.dart';
part 'complex/sparse_complex_matrix.dart';
part 'complex/complex_vector.dart';
part 'complex/complex_matrix.dart';
part 'complex/large_complex_matrix.dart';
part 'complex/wrapper_complex_matrix.dart';
part 'complex/diagonal_complex_matrix.dart';
part 'complex/sparse_cc_complex_matrix.dart';
part 'complex/sparse_rc_complex_matrix.dart';
part 'complex/delegate_complex_vector.dart';

part 'int/abstract_int_vector.dart';
part 'int/abstract_int_matrix.dart';
part 'int/int_vector.dart';
part 'int/int_matrix.dart';
part 'int/sparse_int_vector.dart';
part 'int/sparse_int_matrix.dart';
part 'int/wrapper_int_matrix.dart';
part 'int/diagonal_int_matrix.dart';
part 'int/sparse_cc_int_matrix.dart';
part 'int/sparse_rc_int_matrix.dart';

part 'double/factory.dart';
part 'int/factory.dart';
part 'complex/factory.dart';