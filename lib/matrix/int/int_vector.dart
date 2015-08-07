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

/// Dense 1-d matrix (aka vector) holding [int} elements.
///
/// Internally holds one single contigous one-dimensional array.
class IntVector extends AbstractIntVector {
  Int32List _elements;

  /// Constructs a matrix with a copy of the given [values].
  factory IntVector.fromList(Int32List values) {
    return new IntVector(values.length)..setAll(values);
  }

  /// Constructs a matrix with a given number of cells. All entries are
  /// initially `0`.
  factory IntVector(int size) {
    final elements = new Int32List(size);
    return new IntVector._internal(size, elements, 0, 1, false);
  }

  IntVector._internal(
      int size, Int32List elements, int zero, int stride, bool isView)
      : super(size, zero, stride, !isView) {
    if (elements == null) {
      elements = new Int32List(size);
    }
    _elements = elements;
  }

  static IntVector create(int size) => new IntVector(size);

  factory IntVector.nonzero(List<num> values) {
    var idx = <int>[];
    int i = 0;
    values.forEach((num value) {
      if (value != 0) {
        idx.add(i);
      }
      i++;
    });
    return new IntVector.fromList(idx);
  }

  factory IntVector.random(int size) => ifactory.random(size, create);

  factory IntVector.append(AbstractIntVector a, AbstractIntVector b) {
    return ifactory.append(a, b, create);
  }

  int aggregate(final ifunc.IntIntFunction aggr, final ifunc.IntFunction fn) {
    if (size == 0) {
      throw new ArgumentError("size == 0");
    }
    int a = fn(_elements[zero]);
    int idx = zero;
    for (int i = 1; i < size; i++) {
      idx += stride;
      a = aggr(a, fn(_elements[idx]));
    }
    return a;
  }

  void apply(final ifunc.IntFunction fn) {
    int multiplicator;
    if (fn is ifunc.IntMult) {
      // x[i] = mult*x[i]
      multiplicator = fn.multiplicator;
      if (multiplicator == 1) {
        return;
      }
    } else {
      multiplicator = 0;
    }
    int idx = zero - stride;
    // specialization for speed
    if (fn is ifunc.IntMult) {
      // x[i] = mult*x[i]
      for (int k = size; --k >= 0;) {
        _elements[idx += stride] *= multiplicator;
      }
    } else {
      // the general case x[i] = f(x[i])
      for (int k = size; --k >= 0;) {
        _elements[idx += stride] = fn(_elements[idx]);
      }
    }
  }

  void fill(final int value) {
    final Int32List elems = this._elements;
    int idx = zero;
    for (int i = 0; i < size; i++) {
      elems[idx] = value;
      idx += stride;
    }
  }

  void setAll(final Int32List values) {
    if (values.length != size) {
      throw new ArgumentError(
          "Must have same number of cells: length=${values.length} size()=$size");
    }
    if (!isView) {
      _elements.setAll(0, values);
    } else {
      int idx = zero;
      for (int i = 0; i < size; i++) {
        _elements[idx] = values[i];
        idx += stride;
      }
    }
  }

  void copyFrom(AbstractIntVector source) {
    // overriden for performance only
    if (source is! IntVector) {
      super.copyFrom(source);
      return;
    }
    IntVector other = source as IntVector;
    if (other == this) {
      return;
    }
    checkSize(this, other);
    if (!isView && !other.isView) {
      // quickest
      _elements.setAll(0, other._elements);
      return;
    }
    if (_haveSharedCells(other)) {
      AbstractIntVector c = other.copy();
      if (c is! IntVector) {
        // should not happen
        super.copyFrom(source);
        return;
      }
      other = c as IntVector;
    }

    Int32List elemsOther = other._elements;
    if (_elements == null || elemsOther == null) {
      throw new Error();
    }
    int zeroOther = other.index(0);
    int strideOther = other.stride;
    int idx = zero;
    int idxOther = zeroOther;
    for (int k = 0; k < size; k++) {
      _elements[idx] = elemsOther[idxOther];
      idx += stride;
      idxOther += strideOther;
    }
  }

  void assign(final AbstractIntVector y, final ifunc.IntIntFunction fn) {
    // overriden for performance only
    if (!(y is IntVector)) {
      super.assign(y, fn);
      return;
    }
    checkSize(this, y);
    int zeroOther = y.index(0);
    int strideOther = y.stride;
    Int32List elemsOther = y.elements as Int32List;
    // specialized for speed
    int idx = zero;
    int idxOther = zeroOther;
    if (fn == ifunc.mult) {
      // x[i] = x[i] * y[i]
      for (int k = 0; k < size; k++) {
        _elements[idx] *= elemsOther[idxOther];
        idx += stride;
        idxOther += strideOther;
      }
    } else if (fn == ifunc.div) {
      // x[i] = x[i] / y[i]
      for (int k = 0; k < size; k++) {
        _elements[idx] ~/= elemsOther[idxOther];
        idx += stride;
        idxOther += strideOther;
      }
    } else if (fn is ifunc.IntPlusMultSecond) {
      int multiplicator = fn.multiplicator;
      if (multiplicator == 0) {
        // x[i] = x[i] + 0*y[i]
        return;
      } else if (multiplicator == 1) {
        // x[i] = x[i] + y[i]
        for (int k = 0; k < size; k++) {
          _elements[idx] += elemsOther[idxOther];
          idx += stride;
          idxOther += strideOther;
        }
      } else if (multiplicator == -1) {
        // x[i] = x[i] - y[i]
        for (int k = 0; k < size; k++) {
          _elements[idx] -= elemsOther[idxOther];
          idx += stride;
          idxOther += strideOther;
        }
      } else {
        // the general case x[i] = x[i] + mult*y[i]
        for (int k = 0; k < size; k++) {
          _elements[idx] += multiplicator * elemsOther[idxOther];
          idx += stride;
          idxOther += strideOther;
        }
      }
    } else {
      // the general case x[i] = f(x[i],y[i])
      for (int k = 0; k < size; k++) {
        _elements[idx] = fn(_elements[idx], elemsOther[idxOther]);
        idx += stride;
        idxOther += strideOther;
      }
    }
  }

  int get cardinality {
    int cardinality = 0;
    int idx = zero;
    for (int i = 0; i < size; i++) {
      if (_elements[idx] != 0) {
        cardinality++;
      }
      idx += stride;
    }
    return cardinality;
  }

  Object get elements => _elements;

  void positive({List<int> indexList, final List<int> valueList}) {
    bool fillIndexList = indexList != null;
    bool fillValueList = valueList != null;
    if (fillIndexList) {
      indexList.clear();
    }
    if (fillValueList) {
      valueList.clear();
    }
    int idx = zero;
    int rem = size % 2;
    if (rem == 1) {
      int value = _elements[idx];
      if (value > 0) {
        if (fillIndexList) {
          indexList.add(0);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      idx += stride;
    }
    for (int i = rem; i < size; i += 2) {
      int value = _elements[idx];
      if (value > 0) {
        if (fillIndexList) {
          indexList.add(i);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      idx += stride;
      value = _elements[idx];
      if (value > 0) {
        if (fillIndexList) {
          indexList.add(i + 1);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      idx += stride;
    }
  }

  void negative({List<int> indexList, final List<int> valueList}) {
    bool fillIndexList = indexList != null;
    bool fillValueList = valueList != null;
    if (fillIndexList) {
      indexList.clear();
    }
    if (fillValueList) {
      valueList.clear();
    }
    int idx = zero;
    int rem = size % 2;
    if (rem == 1) {
      int value = _elements[idx];
      if (value < 0) {
        if (fillIndexList) {
          indexList.add(0);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      idx += stride;
    }
    for (int i = rem; i < size; i += 2) {
      int value = _elements[idx];
      if (value < 0) {
        if (fillIndexList) {
          indexList.add(i);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      idx += stride;
      value = _elements[idx];
      if (value < 0) {
        if (fillIndexList) {
          indexList.add(i + 1);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      idx += stride;
    }
  }

  IntVectorLocation max() {
    int maxValue = _elements[zero];
    int location = 0;
    int idx = zero;
    for (int i = 1; i < size; i++) {
      idx += stride;
      if (maxValue < _elements[idx]) {
        maxValue = _elements[idx];
        location = (idx - zero) ~/ stride;
      }
    }
    return new IntVectorLocation._(maxValue, location);
  }

  IntVectorLocation min() {
    int location = 0;
    int minValue = _elements[zero];
    int idx = zero;
    for (int i = 1; i < size; i++) {
      idx += stride;
      if (minValue > _elements[idx]) {
        minValue = _elements[idx];
        location = (idx - zero) ~/ stride;
      }
    }
    return new IntVectorLocation._(minValue, location);
  }

  int get(int index) => _elements[zero + index * stride];

  AbstractIntVector like1D(int size) => new IntVector(size);

  AbstractIntMatrix like2D(int rows, int columns) {
    return new IntMatrix(rows, columns);
  }

  AbstractIntMatrix reshape(final int rows, final int columns) {
    if (rows * columns != size) {
      throw new ArgumentError("rows*columns != size");
    }
    var M = new IntMatrix(rows, columns);
    Int32List elemsOther = M.elements as Int32List;
    int zeroOther = M.index(0, 0);
    int rowStrideOther = M.rowStride;
    int colStrideOther = M.columnStride;
    int idx = zero;
    for (int c = 0; c < columns; c++) {
      var idxOther = zeroOther + c * colStrideOther;
      for (int r = 0; r < rows; r++) {
        elemsOther[idxOther] = _elements[idx];
        idxOther += rowStrideOther;
        idx += stride;
      }
    }
    return M;
  }

  void set(int index, int value) {
    _elements[zero + index * stride] = value;
  }

  void fillList(Int32List values) {
    if (values.length < size) {
      throw new ArgumentError("values too small");
    }
    if (!isView) {
      values.setAll(0, _elements);
    } else {
      super.fillList(values);
    }
  }

  int dot(AbstractIntVector y, [final int from = 0, int length = null]) {
    if (length == null) {
      length = size;
    }
    if (!(y is IntVector)) {
      return super.dot(y, from, length);
    }
    IntVector yy = y as IntVector;

    int tail = from + length;
    if (from < 0 || length < 0) {
      return 0;
    }
    if (size < tail) {
      tail = size;
    }
    if (y.size < tail) {
      tail = y.size;
    }
    Int32List elementsOther = yy._elements;
    int zeroThis = index(from);
    int zeroOther = yy.index(from);
    int strideOther = yy.stride;
    if (_elements == null || elementsOther == null) {
      throw new Error();
    }
    int sum = 0;
    zeroThis -= stride;
    zeroOther -= strideOther;
    int min = tail - from;
    for (int k = min ~/ 4; --k >= 0;) {
      sum += _elements[zeroThis += stride] *
              elementsOther[zeroOther += strideOther] +
          _elements[zeroThis += stride] *
              elementsOther[zeroOther += strideOther] +
          _elements[zeroThis += stride] *
              elementsOther[zeroOther += strideOther] +
          _elements[zeroThis += stride] *
              elementsOther[zeroOther += strideOther];
    }
    for (int k = min % 4; --k >= 0;) {
      sum += _elements[zeroThis += stride] *
          elementsOther[zeroOther += strideOther];
    }
    return sum;
  }

  int sum() {
    int sum = 0;
    Int32List elems = this._elements;
    if (elems == null) {
      throw new Error();
    }
    int idx = zero;
    for (int k = 0; k < size; k++) {
      sum += elems[idx];
      idx += stride;
    }
    return sum;
  }

  bool _haveSharedCellsRaw(AbstractIntVector other) {
    if (other is SelectedDenseIntVector) {
      return _elements == other._elements;
    } else if (other is IntVector) {
      return _elements == other._elements;
    }
    return false;
  }

  int index(int rank) => zero + rank * stride;

  AbstractIntVector _viewSelectionLike(Int32List offsets) {
    return new SelectedDenseIntVector(_elements, offsets);
  }

  Object clone() {
    return new IntVector._internal(size, _elements, zero, stride, isView);
  }

  DoubleVector toDouble() {
    var l = new List<double>.generate(size, (int i) => get(i).toDouble());
    return new DoubleVector(size)..setAll(new Float64List.fromList(l));
  }
}

/// Selection view on dense 1-d matrices holding `int` elements.
///
/// Objects of this class are typically constructed via `select`
/// methods on some source matrix. The interface introduced in abstract super
/// classes defines everything a user can do. From a user point of view there
/// is nothing special about this class; it presents the same functionality
/// with the same signatures and semantics as its abstract superclass(es)
/// while introducing no additional functionality. Thus, this class need not
/// be visible to users.
class SelectedDenseIntVector extends AbstractIntVector {
  Int32List _elements;

  /// The offsets of visible indexes of this matrix.
  Int32List _offsets;

  int __offset;

  factory SelectedDenseIntVector(Int32List elements, Int32List offsets) {
    return new SelectedDenseIntVector._internal(
        offsets.length, elements, 0, 1, offsets, 0);
  }

  SelectedDenseIntVector._internal(int size, Int32List elements, int zero,
      int stride, Int32List offsets, int offset)
      : super(size, zero, stride, false) {
    _elements = elements;
    _offsets = offsets;
    __offset = offset;
  }

  Object get elements {
    throw new ArgumentError("This method is not supported.");
  }

  int get(int index) => _elements[__offset + _offsets[zero + index * stride]];

  int index(int rank) => __offset + _offsets[zero + rank * stride];

  AbstractIntVector like1D(int size) => new IntVector(size);

  AbstractIntMatrix like2D(int rows, int columns) {
    return new IntMatrix(rows, columns);
  }

  /// This method is not supported.
  AbstractIntMatrix reshape(int rows, int columns) {
    throw new ArgumentError("This method is not supported.");
  }

  void set(int index, int value) {
    _elements[__offset + _offsets[zero + index * stride]] = value;
  }

  int _offset(int absRank) => _offsets[absRank];

  bool _haveSharedCellsRaw(AbstractIntVector other) {
    if (other is SelectedDenseIntVector) {
      return _elements == other._elements;
    } else if (other is IntVector) {
      return _elements == other._elements;
    }
    return false;
  }

  AbstractIntVector _viewSelectionLike(Int32List offsets) {
    return new SelectedDenseIntVector(_elements, offsets);
  }

  Object clone() {
    return new SelectedDenseIntVector._internal(
        size, _elements, zero, stride, _offsets, __offset);
  }
}
