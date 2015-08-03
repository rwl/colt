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

/// Sparse hashed 1-d matrix (aka vector) holding [double] elements
/// in a [Map].
///
/// Cells that:
/// - are never set to non-zero values do not use any memory.
/// - switch from zero to non-zero state do use memory.
/// - switch back from non-zero to zero state also do use memory.
///
/// This class offers *expected* time complexity `O(1)` (i.e. constant
/// time) for the basic operations `get`, `set`, and `size`. As such
/// this sparse class is expected to have no worse time complexity than
/// its dense counterpart [DoubleVector]. However, constant factors are
/// considerably larger.
class SparseDoubleVector extends AbstractDoubleVector {
  Map<int, double> _elements;

  /// Constructs a matrix with a copy of the given values.
  factory SparseDoubleVector.fromList(Float64List values) {
    return new SparseDoubleVector(values.length)..setAll(values);
  }

  /// Constructs a matrix with a given number of parameters. All entries are
  /// initially `0`.
  factory SparseDoubleVector(int size) {
    final elements = new Map<int, double>();
    return new SparseDoubleVector._internal(size, elements, 0, 1, false);
  }

  SparseDoubleVector._internal(
      int size, Map<int, double> elements, int offset, int stride, bool isView)
      : super(size, offset, stride, !isView) {
    _elements = elements;
  }

  factory SparseDoubleVector.append(
      AbstractDoubleVector A, AbstractDoubleVector B) {
    return dfactory.append(A, B, (n) => new SparseDoubleVector(n));
  }

  /// Sets all cells to the state specified by [value].
  void fill(double value) {
    if (!isView && value == 0) {
      _elements.clear();
    } else {
      super.fill(value);
    }
  }

  /// Returns the number of cells having non-zero values.
  int get cardinality {
    if (!isView) {
      return _elements.length;
    } else {
      return super.cardinality;
    }
  }

  dynamic get elements => _elements;

  double get(int index) {
    final i = zero + index * stride;
    if (_elements.containsKey(i)) {
      return _elements[i];
    }
    return 0.0;
  }

  int index(int rank) => zero + rank * stride;

  AbstractDoubleVector like1D(int size) => new SparseDoubleVector(size);

  AbstractDoubleMatrix like2D(int rows, int columns) {
    return new SparseDoubleMatrix(rows, columns);
  }

  AbstractDoubleMatrix reshape(int rows, int columns) {
    if (rows * columns != size) {
      throw new ArgumentError("rows*columns != size");
    }
    AbstractDoubleMatrix M = new SparseDoubleMatrix(rows, columns);
    int idx = 0;
    for (int c = 0; c < columns; c++) {
      for (int r = 0; r < rows; r++) {
        double elem = get(idx++);
        if (elem != 0) {
          M.set(r, c, elem);
        }
      }
    }
    return M;
  }

  void set(int index, double value) {
    int i = zero + index * stride;
    if (value == 0) {
      _elements.remove(i);
    } else {
      _elements[i] = value;
    }
  }

  String toString() {
    var buf = new StringBuffer();
    buf.write("1 x $size sparse matrix, nnz = $cardinality\n");
    for (int i = 0; i < size; i++) {
      double elem = get(i);
      if (elem != 0) {
        buf
          ..write('(')
          ..write(i)
          ..write(')')
          ..write('\t')
          ..write(elem)
          ..write('\n');
      }
    }
    return buf.toString();
  }

  bool _haveSharedCellsRaw(AbstractDoubleVector other) {
    if (other is SelectedSparseDoubleVector) {
      return this._elements == other._elements;
    } else if (other is SparseDoubleVector) {
      return this._elements == other._elements;
    }
    return false;
  }

  AbstractDoubleVector _viewSelectionLike(Int32List offsets) {
    return new SelectedSparseDoubleVector(_elements, offsets);
  }

  Object clone() {
    return new SparseDoubleVector._internal(
        size, _elements, zero, stride, isView);
  }
}

/// Selection view on sparse 1-d matrices holding [double] elements.
///
/// Objects of this class are typically constructed via `select` methods on
/// some source matrix. The interface introduced in abstract super classes
/// defines everything a user can do. From a user point of view there is
/// nothing special about this class; it presents the same functionality with
/// the same signatures and semantics as its abstract superclass(es) while
/// introducing no additional functionality. Thus, this class need not be
/// visible to users.
class SelectedSparseDoubleVector extends AbstractDoubleVector {
  Map<int, double> _elements;

  /// The offsets of visible indexes of this matrix.
  Int32List _offsets;

  int __offset;

  factory SelectedSparseDoubleVector(
      Map<int, double> elements, Int32List offsets) {
    return new SelectedSparseDoubleVector._internal(
        offsets.length, elements, 0, 1, offsets, 0);
  }

  SelectedSparseDoubleVector._internal(int size, Map<int, double> elements,
      int zero, int stride, Int32List offsets, int offset)
      : super(size, zero, stride, false) {
    _elements = elements;
    _offsets = new Int32List.fromList(offsets);
    __offset = offset;
  }

  dynamic get elements => _elements;

  double get(int index) {
    final i = __offset + _offsets[zero + index * stride];
    if (_elements.containsKey(i)) {
      return _elements[i];
    }
    return 0.0;
  }

  int index(int rank) => __offset + _offsets[zero + rank * stride];

  AbstractDoubleVector like1D(int size) => new SparseDoubleVector(size);

  AbstractDoubleMatrix like2D(int rows, int columns) {
    return new SparseDoubleMatrix(rows, columns);
  }

  AbstractDoubleMatrix reshape(int rows, int columns) {
    if (rows * columns != size) {
      throw new ArgumentError("rows*columns != size");
    }
    AbstractDoubleMatrix M = new SparseDoubleMatrix(rows, columns);
    int idx = 0;
    for (int c = 0; c < columns; c++) {
      for (int r = 0; r < rows; r++) {
        M.set(r, c, get(idx++));
      }
    }
    return M;
  }

  void set(int index, double value) {
    int i = __offset + _offsets[zero + index * stride];
    if (value == 0) {
      _elements.remove(i);
    } else {
      _elements[i] = value;
    }
  }

  int offset(int absRank) => _offsets[absRank];

  bool _haveSharedCellsRaw(AbstractDoubleVector other) {
    if (other is SelectedSparseDoubleVector) {
      return _elements == other._elements;
    } else if (other is SparseDoubleVector) {
      return _elements == other._elements;
    }
    return false;
  }

  AbstractDoubleVector _viewSelectionLike(Int32List offsets) {
    return new SelectedSparseDoubleVector(_elements, offsets);
  }

  Object clone() {
    return new SelectedSparseDoubleVector._internal(
        size, _elements, zero, stride, _offsets, __offset);
  }
}
