library cern.colt.matrix.int;

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
        checkBox,
        checkColumn,
        checkRow,
        checkColumnIndexes,
        checkRowIndexes,
        vStrides,
        vColumnFlip,
        vRowFlip,
        vDice,
        vBox,
        setIsNoView,
        setRows,
        setColumns;

import '../../function/int.dart' as ifunc;

import 'factory.dart' as ifactory;

import '../../math.dart' show MAX_INT, cumsum, toDense, searchRange;

import 'property.dart' as iprop;

import '../formatter.dart';

import '../double/matrix.dart' show DenseDoubleVector;

part 'abstract_int_vector.dart';
part 'abstract_int_matrix.dart';
part 'dense_int_vector.dart';
part 'dense_int_matrix.dart';
part 'sparse_int_vector.dart';
part 'sparse_int_matrix.dart';
part 'wrapper_int_matrix.dart';
part 'diagonal_int_matrix.dart';
part 'sparse_cc_int_matrix.dart';
part 'sparse_rc_int_matrix.dart';
part 'delegate_int_vector.dart';
part 'large_int_matrix.dart';

class IntVectorLocation {
  final int value, location;
  IntVectorLocation._(this.value, this.location);
}

class IntMatrixLocation {
  final int value, row, column;
  IntMatrixLocation._(this.value, this.row, this.column);
}
