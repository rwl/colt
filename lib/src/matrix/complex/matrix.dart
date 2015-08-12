library cern.colt.matrix.complex;

import 'dart:math' as Math;
import 'dart:typed_data';
import 'package:intl/intl.dart' show NumberFormat;
import 'package:complex/complex.dart';

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

import '../../function/complex.dart' as cfunc;
import '../../function/double.dart' as func;

import '../double/matrix.dart'
    show
        DoubleVector,
        DoubleMatrix,
        DenseDoubleVector,
        DenseDoubleMatrix,
        LargeDoubleMatrix,
        SparseDoubleVector,
        SparseDoubleMatrix,
        SparseRCDoubleMatrix,
        makeSparseRCDoubleMatrix;

import '../../math.dart'
    show EPSILON, DEG_RAD, MAX_INT, isEqual, cumsum, searchRange, EPS;

import 'property.dart' as cprop;

part 'abstract_complex_vector.dart';
part 'abstract_complex_matrix.dart';
part 'sparse_complex_vector.dart';
part 'sparse_complex_matrix.dart';
part 'dense_complex_vector.dart';
part 'dense_complex_matrix.dart';
part 'large_complex_matrix.dart';
part 'wrapper_complex_matrix.dart';
part 'diagonal_complex_matrix.dart';
part 'sparse_cc_complex_matrix.dart';
part 'sparse_rc_complex_matrix.dart';
part 'delegate_complex_vector.dart';
