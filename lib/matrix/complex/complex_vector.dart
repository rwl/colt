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
part of cern.colt.matrix.complex;

/// Dense 1-d matrix (aka vector) holding `complex` elements.
///
/// Internally holds one single contiguous one-dimensional array. Complex
/// data is represented by 2 double values in sequence, i.e.
/// `elements[zero + 2 * k * stride]` constitute real part and
/// `elements[zero + 2 * k * stride + 1]`
/// constitute imaginary part `(k=0,...,size()-1)`.
class ComplexVector extends AbstractComplexVector {
  Float64List _elements;

  /// Constructs a matrix with a copy of the given values. Due to the fact
  /// that complex data is represented by 2 double values in sequence: the
  /// real and imaginary parts, the size of new matrix will be equal to
  /// `values.length / 2`.
  factory ComplexVector.fromList(Float64List values) {
    return new ComplexVector(values.length ~/ 2)..setAll(values);
  }

  /// Constructs a complex matrix with the same size as [realPart] matrix
  /// and fills the real part of this matrix with elements of [realPart].
  factory ComplexVector.fromReal(AbstractDoubleVector realPart) {
    return new ComplexVector(realPart.size)..setReal(realPart);
  }

  factory ComplexVector.fromImaginary(AbstractDoubleVector imagPart) {
    return new ComplexVector(imagPart.size)..setImaginary(imagPart);
  }

  /// Constructs a complex matrix with the same size as [realPart] matrix
  /// and fills the real part of this matrix with elements of [realPart]
  /// and fills the imaginary part of this matrix with elements of
  /// [imaginaryPart].
  factory ComplexVector.fromParts(
      AbstractDoubleVector realPart, AbstractDoubleVector imaginaryPart) {
    return new ComplexVector(realPart.size)
      ..setReal(realPart)
      ..setImaginary(imaginaryPart);
  }

  /// Constructs a matrix with the given size.
  factory ComplexVector(int size) {
    final elements = new Float64List(2 * size);
    return new ComplexVector._internal(size, elements, 0, 2, true);
  }

  ComplexVector._internal(
      int size, Float64List elements, int zero, int stride, bool isNoView)
      : super(size, zero, stride, isNoView) {
    _elements = elements;
  }

  /// Constructs a matrix from the given polar representation. [r] is the
  /// polar radius. [theta] is the polar angle.
  factory ComplexVector.fromPolar(
      AbstractDoubleVector r, AbstractDoubleVector theta,
      [bool radians = true]) {
    final real = theta.copy();
    final imag = theta.copy();
    if (!radians) {
      real.apply(func.multiply(DEG_RAD));
      imag.apply(func.multiply(DEG_RAD));
    }
    real.apply(func.cos);
    real.assign(r, func.mult);
    imag.apply(func.sin);
    imag.assign(r, func.mult);

    return new ComplexVector.fromParts(real, imag);
  }

  static ComplexVector create(int size) {
    return new ComplexVector(size);
  }

  Complex aggregate(final cfunc.ComplexComplexComplexFunction aggr,
      final cfunc.ComplexComplexFunction fn) {
    if (size == 0) {
      return Complex.NAN;
    }
    var a = fn(new Complex(_elements[zero], _elements[zero + 1]));
    int idx = zero;
    for (int i = 1; i < size; i++) {
      idx += stride;
      a = aggr(a, fn(new Complex(_elements[idx], _elements[idx + 1])));
    }
    return a;
  }

  void apply(final cfunc.ComplexComplexFunction fn) {
    if (_elements == null) {
      throw new Error();
    }
    if (fn is cfunc.ComplexMult) {
      Complex multiplicator = fn.multiplicator;
      if (multiplicator.real == 1 && multiplicator.imaginary == 0) {
        return;
      }
    }
    int idx = zero;
    if (fn is cfunc.ComplexMult) {
      Complex multiplicator = fn.multiplicator;
      for (int k = 0; k < size; k++) {
        _elements[idx] = _elements[idx] * multiplicator.real -
            _elements[idx + 1] * multiplicator.imaginary;
        _elements[idx + 1] = _elements[idx + 1] * multiplicator.real +
            _elements[idx] * multiplicator.imaginary;
        idx += stride;
      }
    } else {
      for (int k = 0; k < size; k++) {
        var tmp = new Complex(_elements[idx], _elements[idx + 1]);
        tmp = fn(tmp);
        _elements[idx] = tmp.real;
        _elements[idx + 1] = tmp.imaginary;
        idx += stride;
      }
    }
  }

  void applyReal(final cfunc.ComplexRealFunction fn) {
    if (this._elements == null) {
      throw new Error();
    }
    int idx = zero;
    if (fn == cfunc.abs) {
      for (int k = 0; k < size; k++) {
        double absX = _elements[idx].abs();
        double absY = _elements[idx + 1].abs();
        if (absX == 0.0 && absY == 0.0) {
          _elements[idx] = 0.0;
        } else if (absX >= absY) {
          double d = _elements[idx + 1] / _elements[idx];
          _elements[idx] = absX * Math.sqrt(1.0 + d * d);
        } else {
          double d = _elements[idx] / _elements[idx + 1];
          _elements[idx] = absY * Math.sqrt(1.0 + d * d);
        }
        _elements[idx + 1] = 0.0;
        idx += stride;
      }
    } else {
      for (int k = 0; k < size; k++) {
        var tmp = new Complex(_elements[idx], _elements[idx + 1]);
        _elements[idx] = fn(tmp);
        _elements[idx + 1] = 0.0;
        idx += stride;
      }
    }
  }

  void copyFrom(AbstractComplexVector source) {
    if (source is! ComplexVector) {
      super.copyFrom(source);
      return;
    }
    ComplexVector other = source as ComplexVector;
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
      AbstractComplexVector c = other.copy();
      if (c is! ComplexVector) {
        // should not happen
        super.copyFrom(source);
        return;
      }
      other = c as ComplexVector;
    }

    Float64List elemsOther = other._elements;
    if (_elements == null || elemsOther == null) {
      throw new Error();
    }
    int strideOther = other.stride;
    int zeroOther = other.index(0);

    int idx = zero;
    int idxOther = zeroOther;
    for (int k = 0; k < size; k++) {
      _elements[idx] = elemsOther[idxOther];
      _elements[idx + 1] = elemsOther[idxOther + 1];
      idx += stride;
      idxOther += strideOther;
    }
  }

  void assign(
      AbstractComplexVector y, final cfunc.ComplexComplexComplexFunction fn) {
    if (y is! ComplexVector) {
      super.assign(y, fn);
      return;
    }
    checkSize(this, y);
    Float64List elemsOther = y.elements as Float64List;
    int zeroOther = y.index(0);
    int strideOther = y.stride;

    if (_elements == null || elemsOther == null) {
      throw new Error();
    }
    int idx = zero;
    int idxOther = zeroOther;
    if (fn == cfunc.plus) {
      for (int k = 0; k < size; k++) {
        _elements[idx] += elemsOther[idxOther];
        _elements[idx + 1] += elemsOther[idxOther + 1];
        idx += stride;
        idxOther += strideOther;
      }
    } else if (fn == cfunc.minus) {
      for (int k = 0; k < size; k++) {
        _elements[idx] -= elemsOther[idxOther];
        _elements[idx + 1] -= elemsOther[idxOther + 1];
        idx += stride;
        idxOther += strideOther;
      }
    } else if (fn == cfunc.div) {
      for (int k = 0; k < size; k++) {
        double re = elemsOther[idxOther];
        double im = elemsOther[idxOther + 1];
        Complex tmp;
        if (re.abs() >= im.abs()) {
          var scalar = (1.0 / (re + im * (im / re)));
          tmp = new Complex(
              scalar * (_elements[idx] + _elements[idx + 1] * (im / re)),
              scalar * (_elements[idx + 1] - _elements[idx] * (im / re)));
        } else {
          var scalar = (1.0 / (re * (re / im) + im));
          tmp = new Complex(
              scalar * (_elements[idx] * (re / im) + _elements[idx + 1]),
              scalar * (_elements[idx + 1] * (re / im) - _elements[idx]));
        }
        _elements[idx] = tmp.real;
        _elements[idx + 1] = tmp.imaginary;
        idx += stride;
        idxOther += strideOther;
      }
    } else if (fn == cfunc.mult) {
      for (int k = 0; k < size; k++) {
        var tmp = new Complex(_elements[idx] * elemsOther[idxOther] -
                _elements[idx + 1] * elemsOther[idxOther + 1],
            _elements[idx + 1] * elemsOther[idxOther] +
                _elements[idx] * elemsOther[idxOther + 1]);
        _elements[idx] = tmp.real;
        _elements[idx + 1] = tmp.imaginary;
        idx += stride;
        idxOther += strideOther;
      }
    } else if (fn == cfunc.multConjFirst) {
      for (int k = 0; k < size; k++) {
        var tmp = new Complex(_elements[idx] * elemsOther[idxOther] +
                _elements[idx + 1] * elemsOther[idxOther + 1],
            -_elements[idx + 1] * elemsOther[idxOther] +
                _elements[idx] * elemsOther[idxOther + 1]);
        _elements[idx] = tmp.real;
        _elements[idx + 1] = tmp.imaginary;
        idx += stride;
        idxOther += strideOther;
      }
    } else if (fn == cfunc.multConjSecond) {
      for (int k = 0; k < size; k++) {
        var tmp = new Complex(_elements[idx] * elemsOther[idxOther] +
                _elements[idx + 1] * elemsOther[idxOther + 1],
            _elements[idx + 1] * elemsOther[idxOther] -
                _elements[idx] * elemsOther[idxOther + 1]);
        _elements[idx] = tmp.real;
        _elements[idx + 1] = tmp.imaginary;
        idx += stride;
        idxOther += strideOther;
      }
    } else {
      for (int k = 0; k < size; k++) {
        var tmp1 = new Complex(_elements[idx], _elements[idx + 1]);
        var tmp2 = new Complex(elemsOther[idxOther], elemsOther[idxOther + 1]);
        tmp1 = fn(tmp1, tmp2);
        _elements[idx] = tmp1.real;
        _elements[idx + 1] = tmp1.imaginary;
        idx += stride;
        idxOther += strideOther;
      }
    }
  }

  void fill(final double re, final double im) {
    int idx = zero;
    for (int i = 0; i < size; i++) {
      _elements[idx] = re;
      _elements[idx + 1] = im;
      idx += stride;
    }
  }

  void setAll(Float64List values) {
    if (!isView) {
      if (values.length != 2 * size) {
        throw new ArgumentError(
            "The length of values[] must be equal to 2*size()=${2 * size}");
      }
      _elements.setAll(0, values);
    } else {
      super.setAll(values);
    }
  }

  void setImaginary(final AbstractDoubleVector other) {
    if (other is! DoubleVector) {
      super.setImaginary(other);
      return;
    }
    checkSize(this, other);
    int zeroOther = other.index(0);
    int strideOther = other.stride;
    Float64List elemsOther = other.elements as Float64List;
    int idx = zero;
    int idxOther = zeroOther;
    for (int i = 0; i < size; i++) {
      _elements[idx + 1] = elemsOther[idxOther];
      idx += stride;
      idxOther += strideOther;
    }
  }

  void setReal(final AbstractDoubleVector other) {
    if (other is! DoubleVector) {
      super.setReal(other);
      return;
    }
    checkSize(this, other);
    int zeroOther = other.index(0);
    int strideOther = other.stride;
    Float64List elemsOther = other.elements as Float64List;
    int idx = zero;
    int idxOther = zeroOther;
    for (int i = 0; i < size; i++) {
      _elements[idx] = elemsOther[idxOther];
      idx += stride;
      idxOther += strideOther;
    }
  }

  Object get elements => _elements;

  AbstractDoubleVector imaginary() {
    var Im = new DoubleVector(size);
    Float64List elemsOther = Im.elements;
    int zeroOther = Im.index(0);
    int strideOther = Im.stride;
    int idx = zero;
    int idxOther = zeroOther;
    for (int i = 0; i < size; i++) {
      elemsOther[idxOther] = _elements[idx + 1];
      idx += stride;
      idxOther += strideOther;
    }
    return Im;
  }

  void nonzero({List<int> indexList, List<Complex> valueList}) {
    bool fillIndexList = indexList != null;
    bool fillValueList = valueList != null;
    if (fillIndexList) {
      indexList.clear();
    }
    if (fillValueList) {
      valueList.clear();
    }

    int idx = zero;
    for (int k = 0; k < size; k++) {
      var value = new Complex(_elements[idx], _elements[idx + 1]);
      if (value.real != 0 || value.imaginary != 0) {
        if (fillIndexList) {
          indexList.add(k);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      idx += stride;
    }
  }

  Complex get(int index) {
    int idx = zero + index * stride;
    return new Complex(_elements[idx], _elements[idx + 1]);
  }

  AbstractDoubleVector real() {
    var R = new DoubleVector(size);
    Float64List elemsOther = R.elements;
    int zeroOther = R.index(0);
    int strideOther = R.stride;
    int idx = zero;
    int idxOther = zeroOther;
    for (int i = 0; i < size; i++) {
      elemsOther[idxOther] = _elements[idx];
      idx += stride;
      idxOther += strideOther;
    }
    return R;
  }

  AbstractComplexVector like1D(int size) => new ComplexVector(size);

  AbstractComplexMatrix like2D(int rows, int columns) {
    return new ComplexMatrix(rows, columns);
  }

  AbstractComplexMatrix reshape(final int rows, final int columns) {
    if (rows * columns != size) {
      throw new ArgumentError("rows*columns != size");
    }
    var M = new ComplexMatrix(rows, columns);
    Float64List elemsOther = M.elements as Float64List;
    int zeroOther = M.index(0, 0);
    int rowStrideOther = M.rowStride;
    int columnStrideOther = M.columnStride;
    int idx = zero;
    for (int c = 0; c < columns; c++) {
      var idxOther = zeroOther + c * columnStrideOther;
      for (int r = 0; r < rows; r++) {
        elemsOther[idxOther] = _elements[idx];
        elemsOther[idxOther + 1] = _elements[idx + 1];
        idxOther += rowStrideOther;
        idx += stride;
      }
    }
    return M;
  }

  void setParts(int index, double re, double im) {
    int idx = zero + index * stride;
    _elements[idx] = re;
    _elements[idx + 1] = im;
  }

  void set(int index, Complex value) {
    int idx = zero + index * stride;
    _elements[idx] = value.real;
    _elements[idx + 1] = value.imaginary;
  }

  void fillList(Float64List values) {
    if (values.length < 2 * size) {
      throw new ArgumentError("values too small");
    }
    if (!isView) {
      values.setAll(0, _elements);
    } else {
      super.fillList(values);
    }
  }

  Complex dot(final AbstractComplexVector y,
      [final int from = 0, int length = null]) {
    if (length == null) {
      length = this.size;
    }
    int size = this.size;
    if (from < 0 || length <= 0) {
      return Complex.ZERO;
    }

    int tail = from + length;
    if (size < tail) {
      tail = size;
    }
    if (y.size < tail) {
      tail = y.size;
    }
    length = tail - from;
    Float64List elemsOther = y.elements as Float64List;
    if (_elements == null || elemsOther == null) {
      throw new Error();
    }
    int strideOther = y.stride;
    int zero = index(from);
    int zeroOther = y.index(from);
    var sum = Complex.ZERO;

    int idx = zero;
    int idxOther = zeroOther;
    for (int k = 0; k < length; k++) {
      var re = _elements[idx] * elemsOther[idxOther] +
          _elements[idx + 1] * elemsOther[idxOther + 1];
      var im = _elements[idx + 1] * elemsOther[idxOther] -
          _elements[idx] * elemsOther[idxOther + 1];
      sum += new Complex(re, im);
      idx += stride;
      idxOther += strideOther;
    }
    return sum;
  }

  Complex sum() {
    var sum = Complex.ZERO;
    if (_elements == null) {
      throw new Error();
    }
    int idx = zero;
    for (int k = 0; k < size; k++) {
      sum[0] += _elements[idx];
      sum[1] += _elements[idx + 1];
      idx += stride;
    }
    return sum;
  }

  bool _haveSharedCellsRaw(AbstractComplexVector other) {
    if (other is SelectedDenseComplexVector) {
      return _elements == other._elements;
    } else if (other is ComplexVector) {
      return _elements == other._elements;
    }
    return false;
  }

  int index(int rank) => zero + rank * stride;

  AbstractComplexVector _viewSelectionLike(Int32List offsets) {
    return new SelectedDenseComplexVector(_elements, offsets);
  }

  Object clone() {
    return new ComplexVector._internal(size, _elements, zero, stride, !isView);
  }
}

/// Selection view on dense 1-d matrices holding `complex` elements.
///
/// Objects of this class are typically constructed via `select` methods on
/// some source matrix. The interface introduced in abstract super classes
/// defines everything a user can do. From a user point of view there is
/// nothing special about this class; it presents the same functionality with the
/// same signatures and semantics as its abstract superclass(es) while
/// introducing no additional functionality. Thus, this class need not be visible
/// to users.
class SelectedDenseComplexVector extends AbstractComplexVector {
  Float64List _elements;

  /// The offsets of visible indexes of this matrix.
  Int32List _offsets;

  int __offset;

  factory SelectedDenseComplexVector(Float64List elements, Int32List offsets) {
    return new SelectedDenseComplexVector._internal(
        offsets.length, elements, 0, 1, offsets, 0);
  }

  SelectedDenseComplexVector._internal(int size, Float64List elements, int zero,
      int stride, Int32List offsets, int offset)
      : super(size, zero, stride, false) {
    _elements = elements;
    _offsets = offsets;
    __offset = offset;
  }

  int _offset(int absRank) => _offsets[absRank];

  Complex get(int index) {
    int idx = zero + index * stride;
    return new Complex(_elements[__offset + _offsets[idx]],
        _elements[__offset + _offsets[idx] + 1]);
  }

  AbstractDoubleVector real() {
    var R = new DoubleVector(size);
    for (int i = 0; i < size; i++) {
      var tmp = get(i);
      R.set(i, tmp.real);
    }
    return R;
  }

  AbstractDoubleVector imaginary() {
    var Im = new DoubleVector(size);
    for (int i = 0; i < size; i++) {
      var tmp = get(i);
      Im.set(i, tmp.imaginary);
    }
    return Im;
  }

  /// This method is not supported.
  dynamic get elements {
    throw new UnsupportedError("This method is not supported.");
  }

  bool _haveSharedCellsRaw(AbstractComplexVector other) {
    if (other is SelectedDenseComplexVector) {
      return _elements == other._elements;
    } else if (other is ComplexVector) {
      return _elements == other._elements;
    }
    return false;
  }

  int index(int rank) => __offset + _offsets[zero + rank * stride];

  AbstractComplexVector like1D(int size) => new ComplexVector(size);

  AbstractComplexMatrix like2D(int rows, int columns) {
    return new ComplexMatrix(rows, columns);
  }

  /// This method is not supported.
  AbstractComplexMatrix reshape(int rows, int columns) {
    throw new UnsupportedError("This method is not supported.");
  }

  void set(int index, Complex value) {
    int idx = zero + index * stride;
    _elements[__offset + _offsets[idx]] = value.real;
    _elements[__offset + _offsets[idx] + 1] = value.imaginary;
  }

  void setParts(int index, double re, double im) {
    int idx = zero + index * stride;
    _elements[__offset + _offsets[idx]] = re;
    _elements[__offset + _offsets[idx] + 1] = im;
  }

  AbstractComplexVector _viewSelectionLike(Int32List offsets) {
    return new SelectedDenseComplexVector(_elements, offsets);
  }

  Object clone() {
    return new SelectedDenseComplexVector._internal(
        size, _elements, zero, stride, _offsets, __offset);
  }
}
