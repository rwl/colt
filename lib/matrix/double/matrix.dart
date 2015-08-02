library cern.colt.matrix.double;

import 'dart:math' as Math;
import 'dart:typed_data';
//import 'dart:collection' show ListMixin;

import '../matrix.dart'
    show
        AbstractMatrix,
        AbstractVector,
        checkSize,
        checkIndex,
        checkIndexes,
        vFlip,
        vPart,
        vStride,
        checkShape,
        checkColumn,
        checkRow,
        checkColumnIndexes,
        checkRowIndexes,
        vStrides,
        vColumnFlip,
        vRowFlip,
        vDice,
        vBox;

import '../../function/double.dart'
    show DoubleFunction, DoubleDoubleFunction; //,
//        DoubleDoubleProcedure,
//        DoubleIntProcedure,
//        DoubleProcedure;

import '../../function/double.dart' as func;

import 'factory.dart' as dfactory;

import '../../math.dart' show DoubleMult, DoublePlusMultSecond;

import 'algo/algo.dart' as dprop;

import '../former.dart';

part 'double_vector.dart';
part 'double_matrix.dart';
part 'large_double_matrix.dart';
part 'sparse_double_vector.dart';
part 'sparse_double_matrix.dart';
part 'sparse_rc_double_matrix.dart';
part 'sparse_cc_double_matrix.dart';
part 'diagonal_double_matrix.dart';
part 'wrapper_double_matrix.dart';
part 'delegate_double_vector.dart';
part 'abstract_double_vector.dart';
part 'abstract_double_matrix.dart';


class DoubleVectorLocation {
  final double value;
  final int location;
  DoubleVectorLocation._(this.value, this.location);
}

class DoubleMatrixLocation {
  final double value;
  final int row, column;
  DoubleMatrixLocation._(this.value, this.row, this.column);
}