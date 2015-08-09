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

/// Dense 1-d matrix (aka vector) holding [double] elements.
///
/// Internally holds one single contigous one-dimensional array.
class DenseDoubleVector extends DoubleVector {
  Float64List _elements;

  /// Constructs a matrix with a copy of the given values.
  factory DenseDoubleVector.fromList(Float64List values) {
    return new DenseDoubleVector(values.length)..setAll(values);
  }

  /// Constructs a matrix with a given number of cells. All entries are
  /// initially `0`.
  factory DenseDoubleVector(int size) {
    final elements = new Float64List(size);
    return new DenseDoubleVector._internal(size, elements, 0, 1, false);
  }

  DenseDoubleVector._internal(
      int size, Float64List elements, int zero, int stride, bool isView)
      : super(size, zero, stride, !isView) {
    if (elements == null) {
      elements = new Float64List(size);
    }
    _elements = elements;
  }

  static DenseDoubleVector create(int size) => new DenseDoubleVector(size);

  factory DenseDoubleVector.random(int size) => dfactory.random(size, create);

  factory DenseDoubleVector.append(DoubleVector A, DoubleVector B) {
    return dfactory.append(A, B, create);
  }

  double aggregate(final DoubleDoubleFunction aggr, final DoubleFunction fn) {
    if (size == 0) {
      return double.NAN;
    }
    double a = 0.0;
    int idx = zero + (size - 1) * stride;
    a = fn(_elements[idx]);
    for (int i = size - 1; --i >= 0;) {
      a = aggr(a, fn(_elements[idx -= stride]));
    }
    return a;
  }

  void apply(final DoubleFunction fn) {
    double multiplicator = 0.0;
    if (fn is func.DoubleMult) {
      // x[i] = mult*x[i]
      multiplicator = fn.multiplicator;
      if (multiplicator == 1) {
        return;
      }
    }

    int idx = zero - stride;
    // specialization for speed
    if (fn is func.DoubleMult) {
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

  void fill(final double value) {
    int idx = zero;
    for (int i = 0; i < size; i++) {
      _elements[idx] = value;
      idx += stride;
    }
  }

  void setAll(final Float64List values) {
    if (values.length != size) {
      throw new ArgumentError(
          "Must have same number of cells: length=${values.length} size()=${size}");
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

  void copyFrom(DoubleVector source) {
    // overriden for performance only
    if (source is! DenseDoubleVector) {
      super.copyFrom(source);
      return;
    }

    DenseDoubleVector other = source as DenseDoubleVector;
    if (other == this) {
      return;
    }

    checkSize(this, other);

    if (!isView && !other.isView) {
      _elements.setAll(0, other._elements);
      return;
    }

    if (_haveSharedCells(other)) {
      DoubleVector c = other.copy();
      if (c is! DenseDoubleVector) {
        // should not happen
        super.copyFrom(source);
        return;
      }
      other = c as DenseDoubleVector;
    }

    if (_elements == null || other._elements == null) {
      throw new Error();
    }
    final int zeroOther = other.index(0);
    final int strideOther = other.stride;
    int idx = zero;
    int idxOther = zeroOther;
    for (int k = 0; k < size; k++) {
      _elements[idx] = other._elements[idxOther];
      idx += stride;
      idxOther += strideOther;
    }
  }

  void assign(final DoubleVector y, final DoubleDoubleFunction fn) {
    // overriden for performance only
    if (y is! DenseDoubleVector) {
      super.assign(y, fn);
      return;
    }

    checkSize(this, y);

    final int zeroOther = y.index(0);
    final int strideOther = y.stride;
    final Float64List elementsOther = y.elements as Float64List;
    // specialized for speed
    int idx = zero;
    int idxOther = zeroOther;
    if (fn == func.mult) {
      // x[i] = x[i] * y[i]
      for (int k = 0; k < size; k++) {
        _elements[idx] *= elementsOther[idxOther];
        idx += stride;
        idxOther += strideOther;
      }
    } else if (fn == func.div) {
      // x[i] = x[i] / y[i]
      for (int k = 0; k < size; k++) {
        _elements[idx] /= elementsOther[idxOther];
        idx += stride;
        idxOther += strideOther;
      }
    } else if (fn is func.DoublePlusMultSecond) {
      double multiplicator = fn.multiplicator;
      if (multiplicator == 0) {
        // x[i] = x[i] + 0*y[i]
        return;
      } else if (multiplicator == 1) {
        // x[i] = x[i] + y[i]
        for (int k = 0; k < size; k++) {
          _elements[idx] += elementsOther[idxOther];
          idx += stride;
          idxOther += strideOther;
        }
      } else if (multiplicator == -1) {
        // x[i] = x[i] - y[i]
        for (int k = 0; k < size; k++) {
          _elements[idx] -= elementsOther[idxOther];
          idx += stride;
          idxOther += strideOther;
        }
      } else {
        // the general case x[i] = x[i] + mult*y[i]
        for (int k = 0; k < size; k++) {
          _elements[idx] += multiplicator * elementsOther[idxOther];
          idx += stride;
          idxOther += strideOther;
        }
      }
    } else {
      // the general case x[i] = f(x[i],y[i])
      for (int k = 0; k < size; k++) {
        _elements[idx] = fn(_elements[idx], elementsOther[idxOther]);
        idx += stride;
        idxOther += strideOther;
      }
    }
  }

  int get cardinality {
    int cardinality = 0;
    int idx = zero;
    for (int i = 0; i < size; i++) {
      if (_elements[idx] != 0) cardinality++;
      idx += stride;
    }
    return cardinality;
  }

  Object get elements => _elements;

  void nonzero({List<int> indexList: null, List<double> valueList: null}) {
    bool fillIndexList = indexList != null;
    bool fillValueList = valueList != null;
    if (fillIndexList) {
      indexList.clear();
    }
    if (fillValueList) {
      valueList.clear();
    }
    int rem = size % 2;
    int idx = zero;
    if (rem == 1) {
      double value = _elements[idx];
      if (value != 0) {
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
      double value = _elements[idx];
      if (value != 0) {
        if (fillIndexList) {
          indexList.add(i);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      idx += stride;
      value = _elements[idx];
      if (value != 0) {
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

  void positive({List<int> indexList: null, List<double> valueList: null}) {
    bool fillIndexList = indexList != null;
    bool fillValueList = valueList != null;
    if (fillIndexList) {
      indexList.clear();
    }
    if (fillValueList) {
      valueList.clear();
    }
    int rem = size % 2;
    int idx = zero;
    if (rem == 1) {
      double value = _elements[idx];
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
      double value = _elements[idx];
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

  void negative({List<int> indexList: null, List<double> valueList: null}) {
    bool fillIndexList = indexList != null;
    bool fillValueList = valueList != null;
    if (fillIndexList) {
      indexList.clear();
    }
    if (fillValueList) {
      valueList.clear();
    }
    int rem = size % 2;
    int idx = zero;
    if (rem == 1) {
      double value = _elements[idx];
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
      double value = _elements[idx];
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

  DoubleVectorLocation max() {
    int location = 0;
    double maxValue = 0.0;
    maxValue = _elements[zero];
    location = 0;
    int idx = zero;
    for (int i = 1; i < size; i++) {
      idx += stride;
      if (maxValue < _elements[idx]) {
        maxValue = _elements[idx];
        location = (idx - zero) ~/ stride;
      }
    }
    //}
    return new DoubleVectorLocation._(maxValue, location);
  }

  DoubleVectorLocation min() {
    int location = 0;
    double minValue = 0.0;
    minValue = _elements[zero];
    location = 0;
    int idx = zero;
    for (int i = 1; i < size; i++) {
      idx += stride;
      if (minValue > _elements[idx]) {
        minValue = _elements[idx];
        location = (idx - zero) ~/ stride;
      }
    }
    return new DoubleVectorLocation._(minValue, location);
  }

  double get(int index) => _elements[zero + index * stride];

  DoubleVector like1D(int size) => new DenseDoubleVector(size);

  DoubleMatrix like2D(int rows, int columns) {
    return new DenseDoubleMatrix(rows, columns);
  }

  DoubleMatrix reshape(final int rows, final int columns) {
    if (rows * columns != size) {
      throw new ArgumentError("rows*columns != size");
    }
    var m = new DenseDoubleMatrix(rows, columns);
    final Float64List elementsOther = m.elements as Float64List;
    final int zeroOther = m.index(0, 0);
    final int rowStrideOther = m.rowStride;
    final int columnStrideOther = m.columnStride;
    int idx = zero;
    for (int c = 0; c < columns; c++) {
      int idxOther = zeroOther + c * columnStrideOther;
      for (int r = 0; r < rows; r++) {
        elementsOther[idxOther] = _elements[idx];
        idxOther += rowStrideOther;
        idx += stride;
      }
    }
    return m;
  }

  void set(int index, double value) {
    _elements[zero + index * stride] = value;
  }

  void fillList(Float64List values) {
    if (values.length < size) {
      throw new ArgumentError("values too small");
    }
    if (!isView) {
      values.setAll(0, _elements);
    } else {
      super.fillList(values);
    }
  }

  double dot(DoubleVector y, [final int from = 0, int length = null]) {
    if (length == null) {
      length = size;
    }
    if (y is! DenseDoubleVector) {
      return super.dot(y, from, length);
    }
    DenseDoubleVector yy = y as DenseDoubleVector;

    int tail = from + length;
    if (from < 0 || length < 0) {
      return 0.0;
    }
    if (size < tail) {
      tail = size;
    }
    if (y.size < tail) {
      tail = y.size;
    }
    final Float64List elementsOther = yy._elements;
    int zeroThis = index(from);
    int zeroOther = yy.index(from);
    int strideOther = yy.stride;
    if (_elements == null || elementsOther == null) {
      throw new Error();
    }
    double sum = 0.0;
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

  double sum() {
    double sum = 0.0;
    final Float64List elems = this._elements;
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

  bool _haveSharedCellsRaw(DoubleVector other) {
    if (other is SelectedDenseDoubleVector) {
      return _elements == other._elements;
    } else if (other is DenseDoubleVector) {
      return _elements == other._elements;
    }
    return false;
  }

  int index(int rank) => zero + rank * stride;

  DoubleVector _viewSelectionLike(Int32List offsets) {
    return new SelectedDenseDoubleVector._(_elements, offsets);
  }

  Object clone() {
    return new DenseDoubleVector._internal(size, _elements, zero, stride, isView);
  }
}

/// Selection view on dense 1-d matrices holding [double] elements.
///
/// Objects of this class are typically constructed via `viewIndexes`
/// methods on some source matrix. The interface introduced in abstract super
/// classes defines everything a user can do. From a user point of view there is
/// nothing special about this class; it presents the same functionality with the
/// same signatures and semantics as its abstract superclass(es) while
/// introducing no additional functionality.
///
/// This class uses no delegation. Its instances point directly to the data. Cell
/// addressing overhead is 1 additional array index access per get/set.
class SelectedDenseDoubleVector extends DoubleVector {
  Float64List _elements;

  /// The offsets of visible indexes of this matrix.
  Int32List _offsets;

  int __offset;

  factory SelectedDenseDoubleVector._(Float64List elements, Int32List offsets) {
    return new SelectedDenseDoubleVector._internal(
        offsets.length, elements, 0, 1, offsets, 0);
  }

  SelectedDenseDoubleVector._internal(int size, Float64List elements, int zero,
      int stride, Int32List offsets, int offset)
      : super(size, zero, stride, false) {
    _elements = elements;
    _offsets = offsets;
    __offset = offset;
  }

  Object get elements => _elements;

  double get(int i) => _elements[__offset + _offsets[zero + i * stride]];

  int index(int rank) => __offset + _offsets[zero + rank * stride];

  DoubleVector like1D(int size) => new DenseDoubleVector(size);

  DoubleMatrix like2D(int rows, int columns) {
    return new DenseDoubleMatrix(rows, columns);
  }

  DoubleMatrix reshape(int rows, int columns) {
    if (rows * columns != size) {
      throw new ArgumentError("rows*columns != size");
    }
    var m = new DenseDoubleMatrix(rows, columns);
    final Float64List elementsOther = m.elements as Float64List;
    final int zeroOther = m.index(0, 0);
    final int rowStrideOther = m.rowStride;
    final int colStrideOther = m.columnStride;
    int idxOther;
    int idx = 0;
    for (int c = 0; c < columns; c++) {
      idxOther = zeroOther + c * colStrideOther;
      for (int r = 0; r < rows; r++) {
        elementsOther[idxOther] = get(idx++);
        idxOther += rowStrideOther;
      }
    }
    return m;
  }

  void set(int index, double value) {
    _elements[__offset + _offsets[zero + index * stride]] = value;
  }

  int offset(int absRank) => _offsets[absRank];

  bool _haveSharedCellsRaw(DoubleVector other) {
    if (other is SelectedDenseDoubleVector) {
      return _elements == other._elements;
    } else if (other is DenseDoubleVector) {
      return _elements == other._elements;
    }
    return false;
  }

  DoubleVector _viewSelectionLike(Int32List offsets) {
    return new SelectedDenseDoubleVector._(_elements, offsets);
  }

  Object clone() {
    return new SelectedDenseDoubleVector._internal(
        size, _elements, zero, stride, _offsets, __offset);
  }
}
