library cern.colt.matrix.complex;

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

import '../../function/complex.dart'
    show ComplexFunction, ComplexComplexFunction;

import '../../function/complex.dart' as cfunc;

import '../double/matrix.dart' show AbstractDoubleVector, AbstractDoubleMatrix;

import 'factory.dart' as cfactory;

import '../../math.dart' show EPSILON;

import 'algo/algo.dart' as cprop;

import '../former.dart';

part 'abstract_complex_vector.dart';
part 'abstract_complex_matrix.dart';
part 'sparse_complex_vector.dart';
part 'sparse_complex_matrix.dart';
part 'complex_vector.dart';
part 'complex_matrix.dart';
part 'large_complex_matrix.dart';
part 'wrapper_complex_matrix.dart';
part 'diagonal_complex_matrix.dart';
part 'sparse_cc_complex_matrix.dart';
part 'sparse_rc_complex_matrix.dart';
part 'delegate_complex_vector.dart';
