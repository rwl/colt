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
part of cern.colt.matrix.int;

/// Sparse hashed 2-d matrix holding [int] elements in a [Map].
///
/// Cells that:
/// - are never set to non-zero values do not use any memory.
/// - switch from zero to non-zero state do use memory.
/// - switch back from non-zero to zero state also do use memory.
///
/// This class offers expected time complexity `O(1)` (i.e. constant time)
/// for the basic operations `get`, `set`, and `size`. As such this sparse
/// class is expected to have no worse time complexity than its dense
/// counterpart [IntMatrix]. However, constant factors are considerably
/// larger.
class SparseIntMatrix extends AbstractIntMatrix {
  Map<int, int> _elements;

  /// Constructs a matrix with a given number of rows and columns. All
  /// entries are initially `0`.
  factory SparseIntMatrix(int rows, int columns) {
    return new SparseIntMatrix._internal(rows, columns, new Map<int, int>(), 0, 0, null, 1);
  }

  /// Constructs a matrix with a copy of the given indexes and a single value.
  factory SparseIntMatrix.withValue(int rows, int columns, Int32List rowIndexes, Int32List columnIndexes, int value) {
    return new SparseIntMatrix(rows, columns).._insert(rowIndexes, columnIndexes, value);
  }

  /// Constructs a matrix with a copy of the given indexes and values.
  factory SparseIntMatrix.withValues(int rows, int columns, Int32List rowIndexes, Int32List columnIndexes, Int32List values) {
    return new SparseIntMatrix(rows, columns).._insertValues(rowIndexes, columnIndexes, values);
  }

  /**
   * Constructs a view with the given parameters.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @param elements
   *            the cells.
   * @param rowZero
   *            the position of the first element.
   * @param columnZero
   *            the position of the first element.
   * @param rowStride
   *            the number of elements between two rows, i.e.
   *            `index(i+1,j)-index(i,j)`.
   * @param columnStride
   *            the number of elements between two columns, i.e.
   *            `index(i,j+1)-index(i,j)`.
   * @throws ArgumentError
   *             if
   *             `rows<0 || columns<0 || (double)columns*rows > Integer.MAX_VALUE`
   *             or flip's are illegal.
   */
  SparseIntMatrix._internal(int rows, int columns, Map<int, int> elements, int rowZero, int columnZero, int rowStride, int columnStride)
      : super(rows, columns, rowZero, columnZero, rowStride, columnStride, false) {
    _elements = elements;
  }

  static SparseIntMatrix create(int rows, int columns) {
    return new SparseIntMatrix(rows, columns);
  }

  void fill(int value) {
    if (!isView && value == 0) {
      _elements.clear();
    } else {
      super.fill(value);
    }
    return;
  }

  void copyFrom(AbstractIntMatrix source) {
    // overriden for performance only
    if (source is! SparseIntMatrix) {
      super.copyFrom(source);
      return;
    }
    SparseIntMatrix other = source as SparseIntMatrix;
    if (other == this) {
      return; // nothing to do
    }
    checkShape(this, other);

    if (!isView && !other.isView) { // quickest
      _elements.addAll(other._elements);
      return;
    }
    super.copyFrom(source);
  }

  void assign(final AbstractIntMatrix y, ifunc.IntIntFunction fn) {
    if (isView) {
      super.assign(y, fn);
      return;
    }

    checkShape(this, y);

    if (fn is ifunc.IntPlusMultSecond) { // x[i] = x[i] + alpha*y[i]
      final int alpha = fn.multiplicator;
      if (alpha == 0) {
        return; // nothing to do
      }
      y.forEachNonZero((int i, int j, int value) {
        set(i, j, get(i, j) + alpha * value);
        return value;
      });
    } else if (fn == ifunc.mult) { // x[i] = x[i] * y[i]
      this._elements.forEach((int key, int value) {
        int i = key ~/ columns;
        int j = key % columns;
        int r = value * y.get(i, j);
        if (r != value) {
          _elements[key] = r;
        }
        return true;
      });
    } else if (fn == ifunc.div) { // x[i] = x[i] / y[i]
      this._elements.forEach((int key, int value) {
        int i = key ~/ columns;
        int j = key % columns;
        int r = value ~/ y.get(i, j);
        if (r != value) {
          _elements[key] = r;
        }
        return true;
      });
    } else {
      super.assign(y, fn);
    }
    return;
  }

  int get cardinality {
    if (!isView) {
      return _elements.length;
    } else {
      return super.cardinality;
    }
  }

  /// Returns a new matrix that has the same elements as this matrix, but is
  /// in a column-compressed form.
  SparseCCIntMatrix columnCompressed([bool sortRowIndexes=false]) {
    var nnz = cardinality;
    var keys = _elements.keys;
    var values = _elements.values;
    var rowIndexes = new Int32List(nnz);
    var columnIndexes = new Int32List(nnz);

    for (int k = 0; k < nnz; k++) {
      int key = keys[k];
      rowIndexes[k] = key ~/ columns;
      columnIndexes[k] = key % columns;
    }
    return new SparseCCIntMatrix.withValues(rows, columns, rowIndexes, columnIndexes, values, false, false, sortRowIndexes);
  }

  /// Returns a new matrix that has the same elements as this matrix, but is
  /// in a row-compressed form.
  SparseRCIntMatrix rowCompressed([bool sortColumnIndexes=false]) {
    var nnz = cardinality;
    var keys = _elements.keys;
    var values = _elements.values;
    var rowIndexes = new Int32List(nnz);
    var columnIndexes = new Int32List(nnz);
    for (int k = 0; k < nnz; k++) {
      int key = keys[k];
      rowIndexes[k] = key ~/ columns;
      columnIndexes[k] = key % columns;
    }
    return new SparseRCIntMatrix.withValues(rows, columns, rowIndexes, columnIndexes, values, false, false, sortColumnIndexes);
  }

  Object get elements => _elements;

  void forEachNonZero(final ifunc.IntIntIntFunction fn) {
    if (!isView) {
      _elements.forEach((int key, int value) {
        int i = key ~/ columns;
        int j = key % columns;
        int r = fn(i, j, value);
        if (r != value) {
          _elements[key] = r;
        }
        return true;
      });
    } else {
      super.forEachNonZero(fn);
    }
  }

  int get(int row, int column) {
    return _elements[rowZero + row * rowStride + columnZero + column * columnStride];
  }

  int index(int row, int column) {
    return rowZero + row * rowStride + columnZero + column * columnStride;
  }

  AbstractIntMatrix like2D(int rows, int columns) {
    return new SparseIntMatrix(rows, columns);
  }

  AbstractIntVector like1D(int size) => new SparseIntVector(size);

  void set(int row, int column, int value) {
    int index = rowZero + row * rowStride + columnZero + column * columnStride;
    if (value == 0) {
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
        int elem = get(r, c);
        if (elem != 0) {
          buf
              ..write('(')
              ..write(r)
              ..write(',')
              ..write(c)
              ..write(')')
              ..write('\t')
              ..write(elem)
              ..write('\n');
        }
      }
    }
    return buf.toString();
  }

  AbstractIntVector vectorize() {
    var v = new SparseIntVector(size);
    int idx = 0;
    for (int c = 0; c < columns; c++) {
      for (int r = 0; r < rows; r++) {
        int elem = get(r, c);
        v.set(idx++, elem);
      }
    }
    return v;
  }

  AbstractIntVector mult(AbstractIntVector y, [AbstractIntVector z = null, final int alpha = 1, int beta = null, final bool transposeA = false]) {
    if (beta == null) {
      beta = z == null ? 1 : 0;
    }
    int rowsA = rows;
    int columnsA = columns;
    if (transposeA) {
      rowsA = columns;
      columnsA = rows;
    }

    bool ignore = (z == null);
    if (z == null) {
      z = new IntVector(rowsA);
    }

    if (!(!isView && y is IntVector && z is IntVector)) {
      return super.mult(y, z, alpha, beta, transposeA);
    }

    if (columnsA != y.size || rowsA > z.size) {
      throw new ArgumentError("Incompatible args: " + ((transposeA ? dice() : this).toStringShort()) + ", " + y.toStringShort() + ", " + z.toStringShort());
    }

    if (!ignore) {
      z.apply(ifunc.multiply(beta));
    }

    IntVector zz = z as IntVector;
    Int32List elementsZ = zz._elements;
    int strideZ = zz.stride;
    int zeroZ = z.index(0);

    IntVector yy = y as IntVector;
    Int32List elementsY = yy._elements;
    int strideY = yy.stride;
    int zeroY = y.index(0);

    if (elementsY == null || elementsZ == null) {
      throw new Error();
    }

    _elements.forEach((int key, int value) {
      int i = key ~/ columns;
      int j = key % columns;
      if (transposeA) {
        int tmp = i;
        i = j;
        j = tmp;
      }
      elementsZ[zeroZ + strideZ * i] += alpha * value * elementsY[zeroY + strideY * j];
      return true;
    });

    return z;
  }

  AbstractIntMatrix multiply(AbstractIntMatrix B, [AbstractIntMatrix C = null, final int alpha = 1, int beta = null, final bool transposeA = false, final bool transposeB = false]) {
    if (beta == null) {
      beta = C == null ? 1 : 0;
    }
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
      C = new IntMatrix(rowsA, p);
    }

    if (B.rows != columnsA) {
      throw new ArgumentError("Matrix2D inner dimensions must agree:" + toStringShort() + ", " + (transposeB ? B.dice() : B).toStringShort());
    }
    if (C.rows != rowsA || C.columns != p) {
      throw new ArgumentError("Incompatibel result matrix: " + toStringShort() + ", " + (transposeB ? B.dice() : B).toStringShort() + ", " + C.toStringShort());
    }
    if (this == C || B == C) throw new ArgumentError("Matrices must not be identical");

    if (!ignore) {
      C.apply(ifunc.multiply(beta));
    }

    // cache views
    var Brows = new List<AbstractIntVector>(columnsA);
    for (int i = columnsA; --i >= 0; ) {
      Brows[i] = B.row(i);
    }
    var Crows = new List<AbstractIntVector>(rowsA);
    for (int i = rowsA; --i >= 0; ) {
      Crows[i] = C.row(i);
    }

    var fn = new ifunc.IntPlusMultSecond.plusMult(0);

    _elements.forEach((int key, int value) {
      int i = key ~/ columns;
      int j = key % columns;
      fn.multiplicator = value * alpha;
      if (!transposeA) {
        Crows[i].assign(Brows[j], fn);
      } else {
        Crows[j].assign(Brows[i], fn);
      }
      return true;
    });

    return C;
  }

  void _insert(Int32List rowIndexes, Int32List columnIndexes, int value) {
    int size = rowIndexes.length;
    for (int i = 0; i < size; i++) {
      int row = rowIndexes[i];
      int column = columnIndexes[i];
      if (row >= rows || column >= columns) {
        throw new RangeError("row: $row, column: $column");
      }
      if (value != 0) {
        int index = rowZero + row * rowStride + columnZero + column * columnStride;
        int elem = _elements[index];
        if (elem == 0) {
          _elements[index] = value;
        } else {
          int sum = elem + value;
          if (sum == 0) {
            _elements.remove(index);
          } else {
            _elements[index] = sum;
          }
        }
      }
    }
  }

  void _insertValues(Int32List rowIndexes, Int32List columnIndexes, Int32List values) {
    int size = rowIndexes.length;
    for (int i = 0; i < size; i++) {
      int value = values[i];
      int row = rowIndexes[i];
      int column = columnIndexes[i];
      if (row >= rows || column >= columns) {
        throw new RangeError("row: $row, column: $column");
      }
      if (value != 0) {
        int index = rowZero + row * rowStride + columnZero + column * columnStride;
        int elem = _elements[index];
        if (elem == 0) {
          _elements[index] = value;
        } else {
          int sum = elem + value;
          if (sum == 0) {
            _elements.remove(index);
          } else {
            _elements[index] = sum;
          }
        }
      }
    }
  }

  bool _haveSharedCellsRaw(AbstractIntMatrix other) {
    if (other is SelectedSparseIntMatrix) {
      return _elements == other._elements;
    } else if (other is SparseIntMatrix) {
      return _elements == other._elements;
    }
    return false;
  }

  AbstractIntVector _like1D(int size, int offset, int stride) {
    return new SparseIntVector._internal(size, _elements, offset, stride);
  }

  AbstractIntMatrix _viewSelectionLike(Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedSparseIntMatrix(_elements, rowOffsets, columnOffsets, 0);
  }

  Object clone() {
    return new SparseIntMatrix._internal(rows, columns, _elements, rowZero, columnZero, rowStride, columnStride);
  }

}

/// Selection view on sparse 2-d matrices holding [int] elements.
///
/// Objects of this class are typically constructed via `select` methods on
/// some source matrix. The interface introduced in abstract super classes
/// defines everything a user can do. From a user point of view there is
/// nothing special about this class; it presents the same functionality with
/// the same signatures and semantics as its abstract superclass(es) while
/// introducing no additional functionality. Thus, this class need not be
/// visible to users.
class SelectedSparseIntMatrix extends AbstractIntMatrix {
  Map<int, int> _elements;

  /// The offsets of the visible cells of this matrix.
  Int32List _rowOffsets;

  Int32List _columnOffsets;

  int _offset;

  factory SelectedSparseIntMatrix(Map<int, int> elements, Int32List rowOffsets, Int32List columnOffsets, int offset) {
    return new SelectedSparseIntMatrix._internal(rowOffsets.length, columnOffsets.length, elements, 0, 0, 1, 1, rowOffsets, columnOffsets, offset);
  }

  SelectedSparseIntMatrix._internal(int rows, int columns, Map<int, int> elements, int rowZero, int columnZero, int rowStride, int columnStride, Int32List rowOffsets, Int32List columnOffsets, int offset)
      : super(rows, columns, rowZero, columnZero, rowStride, columnStride, false) {
    _elements = elements;
    _rowOffsets = rowOffsets;
    _columnOffsets = columnOffsets;
    _offset = offset;
  }

  // This method is not supported.
  Object get elements {
    throw new ArgumentError("This method is not supported.");
  }

  int get(int row, int column) {
    return _elements[_offset + _rowOffsets[rowZero + row * rowStride] + _columnOffsets[columnZero + column * columnStride]];
  }

  int index(int row, int column) {
    return _offset + _rowOffsets[rowZero + row * rowStride] + _columnOffsets[columnZero + column * columnStride];
  }

  AbstractIntMatrix like2D(int rows, int columns) {
    return new SparseIntMatrix(rows, columns);
  }

  AbstractIntVector like1D(int size) => new SparseIntVector(size);

  void set(int row, int column, int value) {
    int index = _offset + _rowOffsets[rowZero + row * rowStride] + _columnOffsets[columnZero + column * columnStride];

    if (value == 0) {
      _elements.remove(index);
    } else {
      _elements[index] = value;
    }
  }

  AbstractIntVector column(int column) {
    checkColumn(this, column);
    int viewSize = rows;
    int viewZero = rowZero;
    int viewStride = rowStride;
    Int32List viewOffsets = _rowOffsets;
    int viewOffset = _offset + columnOffset(columnRank(column));
    return new SelectedSparseIntVector._internal(viewSize, _elements, viewZero, viewStride, viewOffsets, viewOffset);
  }

  AbstractIntVector row(int row) {
    checkRow(this, row);
    int viewSize = columns;
    int viewZero = columnZero;
    int viewStride = columnStride;
    Int32List viewOffsets = _columnOffsets;
    int viewOffset = _offset + rowOffset(rowRank(row));
    return new SelectedSparseIntVector._internal(viewSize, _elements, viewZero, viewStride, viewOffsets, viewOffset);
  }

  int columnOffset(int absRank) => _columnOffsets[absRank];

  int rowOffset(int absRank) => _rowOffsets[absRank];

  bool _haveSharedCellsRaw(AbstractIntMatrix other) {
    if (other is SelectedSparseIntMatrix) {
      return _elements == other._elements;
    } else if (other is SparseIntMatrix) {
      return _elements == other._elements;
    }
    return false;
  }

  AbstractIntVector _like1D(int size, int zero, int stride) {
    // never called since row() and column() are overridden.
    throw new Error();
  }

  AbstractIntMatrix dice() {
    var v = _view();
    vDice(v);
    Int32List tmp = v._rowOffsets;
    v._rowOffsets = v._columnOffsets;
    v._columnOffsets = tmp;
    setIsNoView(v, false);
    return v;
  }

  AbstractIntMatrix _viewSelectionLike(Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedSparseIntMatrix(this._elements, rowOffsets, columnOffsets, this._offset);
  }

  Object clone() {
    return new SelectedSparseIntMatrix._internal(rows, columns, _elements, rowZero, columnZero, rowStride, columnStride, _rowOffsets, _columnOffsets, _offset);
  }
}
