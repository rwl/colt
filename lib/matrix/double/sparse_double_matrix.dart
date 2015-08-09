// Copyright (C) 1999 CERN - European Organization for Nuclear Research.
//
// Permission to use, copy, modify, distribute and sell this software and
// its documentation for any purpose is hereby granted without fee, provided
// that the above copyright notice appear in all copies and that both that
// copyright notice and this permission notice appear in supporting
// documentation.
//
// CERN makes no representations about the suitability of this software for
// any purpose. It is provided "as is" without expressed or implied warranty.
part of cern.colt.matrix.double;

/// Sparse hashed 2-d matrix holding [double] elements in a [Map].
///
/// Cells that:
/// - are never set to non-zero values do not use any memory.
/// - switch from zero to non-zero state do use memory.
/// - switch back from non-zero to zero state also do use memory. However,
/// their memory is automatically reclaimed from time to time.
///
///
/// This class offers *expected* time complexity `O(1)` (i.e. constant time)
/// for the basic operations `get`, `set` and `size`. As such this sparse
/// class is expected to have no worse time complexity than its dense
/// counterpart [DenseDoubleMatrix]. However, constant factors are considerably
/// larger.
///
/// Cells are internally addressed in row-major order. Performance sensitive
/// applications can exploit this fact.
class SparseDoubleMatrix extends DoubleMatrix {
  Map<int, double> _elements;

  /// Constructs a matrix with a given number of rows and columns.
  factory SparseDoubleMatrix(int rows, int columns) {
    final elements = new Map<int, double>();
    return new SparseDoubleMatrix._internal(
        rows, columns, elements, 0, 0, columns, 1, false);
  }

  /// Constructs a matrix with a copy of the given indexes and a single value.
  factory SparseDoubleMatrix.withValue(int rows, int columns,
      Int32List rowIndexes, Int32List columnIndexes, double value) {
    return new SparseDoubleMatrix(rows, columns)
      .._insert(rowIndexes, columnIndexes, value);
  }

  /// Constructs a matrix with a copy of the given indexes and values.
  factory SparseDoubleMatrix.withValues(int rows, int columns,
      Int32List rowIndexes, Int32List columnIndexes, Float64List values) {
    return new SparseDoubleMatrix(rows, columns)
      .._insertValues(rowIndexes, columnIndexes, values);
  }

  SparseDoubleMatrix._internal(int rows, int columns, Map<int, double> elements,
      int rowZero, int columnZero, int rowStride, int columnStride, bool isView)
      : super(rows, columns, rowZero, columnZero, rowStride, columnStride,
          !isView) {
    _elements = elements;
  }

  static SparseDoubleMatrix create(int rows, int columns) {
    return new SparseDoubleMatrix(rows, columns);
  }

  factory SparseDoubleMatrix.compose(List<List<DoubleMatrix>> parts) {
    return dfactory.compose(
        parts, (rows, columns) => new SparseDoubleMatrix(rows, columns));
  }

  void fill(double value) {
    if (!isView && value == 0) {
      _elements.clear();
    } else {
      super.fill(value);
    }
  }

  void copyFrom(DoubleMatrix source) {
    if (!(source is SparseDoubleMatrix)) {
      super.copyFrom(source);
      return;
    }
    SparseDoubleMatrix other = source as SparseDoubleMatrix;
    if (other == this) {
      return; // nothing to do
    }
    checkShape(this, other);
    super.copyFrom(source);
  }

  void assign(final DoubleMatrix y, func.DoubleDoubleFunction fn) {
    if (isView) {
      super.assign(y, fn);
      return;
    }

    checkShape(this, y);

    if (fn is func.DoublePlusMultSecond) {
      // x[i] = x[i] + alpha*y[i]
      final double alpha = fn.multiplicator;
      if (alpha == 0) {
        return; // nothing to do
      }
      y.forEachNonZero((int i, int j, double value) {
        set(i, j, get(i, j) + alpha * value);
        return value;
      });
    } else if (fn == func.mult) {
      // x[i] = x[i] * y[i]
      this._elements.forEach((int key, double value) {
        int i = key ~/ columns;
        int j = key % columns;
        double r = value * y.get(i, j);
        if (r != value) _elements[key] = r;
        return true;
      });
    } else if (fn == func.div) {
      // x[i] = x[i] /  y[i]
      this._elements.forEach((int key, double value) {
        int i = key ~/ columns;
        int j = key % columns;
        double r = value / y.get(i, j);
        if (r != value) {
          _elements[key] = r;
        }
        return true;
      });
    } else {
      super.assign(y, fn);
    }
  }

  int get cardinality {
    if (!isView) {
      return _elements.length;
    } else {
      return super.cardinality;
    }
  }

  /// Returns a new matrix that has the same elements as this matrix, but is in
  /// a column-compressed form.
  SparseCCDoubleMatrix columnCompressed([bool sortRowIndexes = false]) {
    int nnz = cardinality;
    var keys = _elements.keys;
    var values = _elements.values;
    var rowIndexes = new Int32List(nnz);
    var columnIndexes = new Int32List(nnz);

    for (int k = 0; k < nnz; k++) {
      final key = keys[k];
      rowIndexes[k] = key ~/ columns;
      columnIndexes[k] = key % columns;
    }
    return new SparseCCDoubleMatrix.withValues(rows, columns, rowIndexes,
        columnIndexes, values, false, false, sortRowIndexes);
  }

  /// Returns a new matrix that has the same elements as this matrix, but is
  /// in a row-compressed form.
  SparseRCDoubleMatrix rowCompressed([bool sortColumnIndexes = false]) {
    int nnz = cardinality;
    var keys = _elements.keys;
    var values = _elements.values;
    var rowIndexes = new Int32List(nnz);
    var columnIndexes = new Int32List(nnz);
    for (int k = 0; k < nnz; k++) {
      int key = keys[k];
      rowIndexes[k] = key ~/ columns;
      columnIndexes[k] = key % columns;
    }
    return new SparseRCDoubleMatrix.withValues(
        rows, columns, rowIndexes, columnIndexes, values,
        removeDuplicates: false,
        removeZeroes: false,
        sortColumnIndexes: sortColumnIndexes);
  }

  Object get elements => _elements;

  void forEachNonZero(final func.IntIntDoubleFunction fn) {
    if (!isView) {
      _elements.forEach((int key, double value) {
        int i = key ~/ columns;
        int j = key % columns;
        double r = fn(i, j, value);
        if (r != value) _elements[key] = r;
        return true;
      });
    } else {
      super.forEachNonZero(fn);
    }
  }

  double get(int row, int column) {
    final i = rowZero + row * rowStride + columnZero + column * columnStride;
    if (_elements.containsKey(i)) {
      return _elements[i];
    }
    return 0.0;
  }

  int index(int row, int column) {
    return rowZero + row * rowStride + columnZero + column * columnStride;
  }

  DoubleMatrix like2D(int rows, int columns) {
    return new SparseDoubleMatrix(rows, columns);
  }

  DoubleVector like1D(int size) => new SparseDoubleVector(size);

  void set(int row, int column, double value) {
    int index = rowZero + row * rowStride + columnZero + column * columnStride;
    if (value == 0.0) {
      _elements.remove(index);
    } else {
      _elements[index] = value;
    }
  }

  String toString() {
    var buf = new StringBuffer();
    buf.write("$rows x $columns sparse matrix, nnz = $cardinality\n");
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        double elem = get(r, c);
        if (elem != 0) {
          buf.write('($r,c)    $elem\n');
        }
      }
    }
    return buf.toString();
  }

  DoubleVector mult(DoubleVector y,
      [DoubleVector z = null, final double alpha = 1.0,
      double beta = 0.0, final bool transposeA = false]) {
    int rowsA = rows;
    int columnsA = columns;
    if (transposeA) {
      rowsA = columns;
      columnsA = rows;
    }

    bool ignore = (z == null);
    if (z == null) {
      z = new DenseDoubleVector(rowsA);
    }

    if (!(!isView && y is DenseDoubleVector && z is DenseDoubleVector)) {
      return super.mult(y, z, alpha, beta, transposeA);
    }

    if (columnsA != y.size || rowsA > z.size) {
      throw new ArgumentError("Incompatible args: " +
          ((transposeA ? dice() : this).toStringShort()) +
          ", " +
          y.toStringShort() +
          ", " +
          z.toStringShort());
    }

    if (!ignore) {
      z.apply(func.multiply(beta));
    }

    DenseDoubleVector zz = z as DenseDoubleVector;
    Float64List elementsZ = zz._elements;
    int strideZ = zz.stride;
    int zeroZ = z.index(0);

    DenseDoubleVector yy = y as DenseDoubleVector;
    Float64List elementsY = yy._elements;
    int strideY = yy.stride;
    int zeroY = y.index(0);

    if (elementsY == null || elementsZ == null) {
      throw new Error();
    }

    _elements.forEach((int key, double value) {
      int i = key ~/ columns;
      int j = key % columns;
      if (transposeA) {
        int tmp = i;
        i = j;
        j = tmp;
      }
      elementsZ[zeroZ + strideZ * i] +=
          alpha * value * elementsY[zeroY + strideY * j];
      return true;
    });

    return z;
  }

  DoubleMatrix multiply(DoubleMatrix B,
      [DoubleMatrix C = null, final double alpha = 1.0,
      double beta = 0.0, final bool transposeA = false,
      bool transposeB = false]) {
    if (isView) {
      return super.multiply(B, C, alpha, beta, transposeA, transposeB);
    }
    if (transposeB) {
      B = B.dice();
    }
    int rowsA = rows;
    int columnsA = columns;
    if (transposeA) {
      rowsA = columns;
      columnsA = rows;
    }
    int p = B.columns;
    bool ignore = (C == null);
    if (C == null) {
      C = new DenseDoubleMatrix(rowsA, p);
    }

    if (B.rows != columnsA) {
      throw new ArgumentError("Matrix inner dimensions must agree:" +
          toStringShort() +
          ", " +
          (transposeB ? B.dice() : B).toStringShort());
    }
    if (C.rows != rowsA || C.columns != p) {
      throw new ArgumentError("Incompatibel result matrix: " +
          toStringShort() +
          ", " +
          (transposeB ? B.dice() : B).toStringShort() +
          ", " +
          C.toStringShort());
    }
    if (this == C || B == C) {
      throw new ArgumentError("Matrices must not be identical");
    }

    if (!ignore) {
      C.apply(func.multiply(beta));
    }

    // cache views
    var Brows = new List<DoubleVector>(columnsA);
    for (int i = columnsA; --i >= 0;) {
      Brows[i] = B.row(i);
    }
    var Crows = new List<DoubleVector>(rowsA);
    for (int i = rowsA; --i >= 0;) {
      Crows[i] = C.row(i);
    }

    var fun = new func.DoublePlusMultSecond.plusMult(0.0);

    _elements.forEach((int key, double value) {
      int i = key ~/ columns;
      int j = key % columns;
      fun.multiplicator = value * alpha;
      if (!transposeA) {
        Crows[i].assign(Brows[j], fun);
      } else {
        Crows[j].assign(Brows[i], fun);
      }
      return true;
    });

    return C;
  }

  void _insert(Int32List rowIndexes, Int32List columnIndexes, double value) {
    int size = rowIndexes.length;
    for (int i = 0; i < size; i++) {
      int row = rowIndexes[i];
      int column = columnIndexes[i];
      if (row >= rows || column >= columns) {
        throw new RangeError("row: $row, column: $column");
      }
      if (value != 0) {
        int index =
            rowZero + row * rowStride + columnZero + column * columnStride;
        double elem = _elements[index];
        if (elem == 0) {
          _elements[index] = value;
        } else {
          double sum = elem + value;
          if (sum == 0) {
            _elements.remove(index);
          } else {
            _elements[index] = sum;
          }
        }
      }
    }
  }

  void _insertValues(
      Int32List rowIndexes, Int32List columnIndexes, Float64List values) {
    int size = rowIndexes.length;
    for (int i = 0; i < size; i++) {
      double value = values[i];
      int row = rowIndexes[i];
      int column = columnIndexes[i];
      if (row >= rows || column >= columns) {
        throw new RangeError("row: $row, column: $column");
      }
      if (value != 0) {
        int index =
            rowZero + row * rowStride + columnZero + column * columnStride;
        double elem = _elements[index];
        if (elem == 0) {
          _elements[index] = value;
        } else {
          double sum = elem + value;
          if (sum == 0) {
            _elements.remove(index);
          } else {
            _elements[index] = sum;
          }
        }
      }
    }
  }

  bool _haveSharedCellsRaw(DoubleMatrix other) {
    if (other is SelectedSparseDoubleMatrix) {
      return _elements == other._elements;
    } else if (other is SparseDoubleMatrix) {
      return _elements == other._elements;
    }
    return false;
  }

  DoubleVector _like1D(int size, int offset, int stride) {
    return new SparseDoubleVector._internal(
        size, _elements, offset, stride, false);
  }

  DoubleMatrix _viewSelectionLike(
      Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedSparseDoubleMatrix(
        _elements, rowOffsets, columnOffsets, 0);
  }

  Object clone() {
    return new SparseDoubleMatrix._internal(rows, columns, _elements, rowZero,
        columnZero, rowStride, columnStride, isView);
  }
}

/// Selection view on sparse 2-d matrices holding <tt>double</tt> elements.
///
/// Objects of this class are typically constructed via `select`
/// methods on some source matrix. The interface introduced in abstract super
/// classes defines everything a user can do. From a user point of view there is
/// nothing special about this class; it presents the same functionality with the
/// same signatures and semantics as its abstract superclass(es) while
/// introducing no additional functionality. Thus, this class need not be visible
/// to users.
class SelectedSparseDoubleMatrix extends DoubleMatrix {
  Map<int, double> _elements;

  /// The offsets of the visible cells of this matrix.
  Int32List _rowOffsets;

  Int32List _columnOffsets;

  int _offset;

  factory SelectedSparseDoubleMatrix(Map<int, double> elements,
      Int32List rowOffsets, Int32List columnOffsets, int offset) {
    return new SelectedSparseDoubleMatrix._internal(rowOffsets.length,
        columnOffsets.length, elements, 0, 0, 1, 1, rowOffsets, columnOffsets,
        offset);
  }

  SelectedSparseDoubleMatrix._internal(int rows, int columns,
      Map<int, double> elements, int rowZero, int columnZero, int rowStride,
      int columnStride, Int32List rowOffsets, Int32List columnOffsets,
      int offset)
      : super(
          rows, columns, rowZero, columnZero, rowStride, columnStride, false) {
    _elements = elements;
    _rowOffsets = new Int32List.fromList(rowOffsets);
    _columnOffsets = new Int32List.fromList(columnOffsets);
    _offset = offset;
  }

  Object get elements => _elements;

  double get(int row, int column) {
    return _elements[_offset +
        _rowOffsets[rowZero + row * rowStride] +
        _columnOffsets[columnZero + column * columnStride]];
  }

  int index(int row, int column) {
    return _offset +
        _rowOffsets[rowZero + row * rowStride] +
        _columnOffsets[columnZero + column * columnStride];
  }

  DoubleMatrix like2D(int rows, int columns) {
    return new SparseDoubleMatrix(rows, columns);
  }

  DoubleVector like1D(int size) => new SparseDoubleVector(size);

  void set(int row, int column, double value) {
    int index = _offset +
        _rowOffsets[rowZero + row * rowStride] +
        _columnOffsets[columnZero + column * columnStride];

    if (value == 0) {
      _elements.remove(index);
    } else {
      _elements[index] = value;
    }
  }

  DoubleVector column(int column) {
    checkColumn(this, column);
    int viewSize = rows;
    int viewZero = rowZero;
    int viewStride = rowStride;
    Int32List viewOffsets = _rowOffsets;
    int viewOffset = _offset + columnOffset(columnRank(column));
    return new SelectedSparseDoubleVector._internal(
        viewSize, _elements, viewZero, viewStride, viewOffsets, viewOffset);
  }

  DoubleVector row(int row) {
    checkRow(this, row);
    int viewSize = columns;
    int viewZero = columnZero;
    int viewStride = columnStride;
    Int32List viewOffsets = _columnOffsets;
    int viewOffset = _offset + _rowOffset(rowRank(row));
    return new SelectedSparseDoubleVector._internal(viewSize, this._elements,
        viewZero, viewStride, viewOffsets, viewOffset);
  }

  int columnOffset(int absRank) => _columnOffsets[absRank];

  int _rowOffset(int absRank) => _rowOffsets[absRank];

  bool _haveSharedCellsRaw(DoubleMatrix other) {
    if (other is SelectedSparseDoubleMatrix) {
      return _elements == other._elements;
    } else if (other is SparseDoubleMatrix) {
      return _elements == other._elements;
    }
    return false;
  }

  DoubleVector _like1D(int size, int zero, int stride) {
    // this method is never called since [row] and [column] are overridden
    throw new Error();
  }

  DoubleMatrix dice() {
    var v = _view();
    vDice(v);
    Int32List tmp = _rowOffsets;
    _rowOffsets = _columnOffsets;
    _columnOffsets = tmp;
    setIsNoView(v, false);
    return v;
  }

  DoubleMatrix _viewSelectionLike(
      Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedSparseDoubleMatrix(
        _elements, rowOffsets, columnOffsets, this._offset);
  }

  Object clone() {
    return new SelectedSparseDoubleMatrix._internal(rows, columns, _elements,
        rowZero, columnZero, rowStride, columnStride, _rowOffsets,
        _columnOffsets, _offset);
  }
}
