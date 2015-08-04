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

/// Abstract base class for 1-d matrices (aka vectors) holding `complex`
/// elements.
abstract class AbstractComplexVector extends AbstractVector {
  //with ListMixin<double> {
  AbstractComplexVector(int size,
      [int zero = 0, int stride = 1, bool isNoView = true])
      : super(size, zero, stride, isNoView);

  /// Applies a function to each cell and aggregates the results.
  Float64List aggregate(final cfunc.ComplexComplexComplexFunction aggr,
      final cfunc.ComplexComplexFunction f) {
    Float64List b = new Float64List(2);
    if (size == 0) {
      b[0] = double.NAN;
      b[1] = double.NAN;
      return b;
    }
    Float64List a = f(get(0));
    for (int i = 1; i < size; i++) {
      a = aggr(a, f(get(i)));
    }
    return a;
  }

  /// Applies a function to each cell.
  void apply(final cfunc.ComplexComplexFunction f) {
    for (int i = 0; i < size; i++) {
      set(i, f(get(i)));
    }
  }

  /// Applies a function to the real part of the receiver. The
  /// imaginary part of the receiver is reset to zero.
  void applyReal(final cfunc.ComplexRealFunction f) {
    for (int i = 0; i < size; i++) {
      setParts(i, f(get(i)), 0.0);
    }
  }

  /// Replaces all cell values of the receiver with the values of another
  /// matrix. Both matrices must have the same size. If both matrices share
  /// the same cells (as is the case if they are views derived from the same
  /// matrix) and intersect in an ambiguous way, then replaces as if using
  /// an intermediate auxiliary deep copy of [other].
  void copyFrom(AbstractComplexVector other) {
    if (other == this) {
      return;
    }
    checkSize(this, other);
    AbstractComplexVector otherLoc;
    if (_haveSharedCells(other)) {
      otherLoc = other.copy();
    } else {
      otherLoc = other;
    }
    for (int i = 0; i < size; i++) {
      set(i, otherLoc.get(i));
    }
  }

  /// Assigns the result of a function to each cell.
  void assign(final AbstractComplexVector y,
      final cfunc.ComplexComplexComplexFunction f) {
    checkSize(this, y);
    for (int i = 0; i < size; i++) {
      set(i, f(get(i), y.get(i)));
    }
  }

  /// Sets all cells to the state specified by [re] and [im].
  void fill(final double re, final double im) {
    for (int i = 0; i < size; i++) {
      setParts(i, re, im);
    }
  }

  /// Sets all cells to the state specified by [values]. [values] is required
  /// to have the same number of cells as the receiver. Complex data is
  /// represented by 2 double values in sequence: the real and imaginary
  /// parts, i.e. input array must be of size `2*size`.
  void setAll(final Float64List values) {
    if (values.length != 2 * size) {
      throw new ArgumentError(
          "The length of values[] must be equal to 2*size()=$size");
    }
    for (int i = 0; i < size; i++) {
      setParts(i, values[2 * i], values[2 * i + 1]);
    }
  }

  /// Replaces imaginary part of the receiver with the values of another
  /// real matrix. The real part remains unchanged. Both matrices must
  /// have the same size.
  void setImaginary(final AbstractDoubleVector other) {
    checkSize(this, other);
    for (int i = 0; i < size; i++) {
      double re = get(i)[0];
      double im = other.get(i);
      setParts(i, re, im);
    }
  }

  /// Replaces real part of the receiver with the values of another real
  /// matrix. The imaginary part remains unchanged. Both matrices must
  /// have the same size.
  void setReal(final AbstractDoubleVector other) {
    checkSize(this, other);
    for (int i = 0; i < size; i++) {
      double re = other.get(i);
      double im = get(i)[1];
      setParts(i, re, im);
    }
  }

  /// Returns the number of cells having non-zero values; ignores tolerance.
  int get cardinality {
    int cardinality = 0;
    for (int i = 0; i < size; i++) {
      var tmp = get(i);
      if ((tmp[0] != 0.0) || (tmp[1] != 0.0)) {
        cardinality++;
      }
    }
    return cardinality;
  }

  /// Constructs and returns a deep copy of the receiver.
  AbstractComplexVector copy() {
    return like()..copyFrom(this);
  }

  /// Compares this object against the specified object. The result is
  /// `true` if and only if the argument is not `null` and is at least
  /// a [AbstractComplexVector] object that has the same sizes as the
  /// receiver and has exactly the same values at the same indexes.
  bool equals(AbstractComplexVector obj) {
    if (identical(this, obj)) {
      return true;
    }
    if (obj == null) {
      return false;
    }

    return cprop.equalsVector(this, obj);
  }

  /// Returns the matrix cell value at coordinate `index`.
  Float64List operator [](int index) {
    checkIndex(this, index);
    return get(index);
  }

  /// Returns the elements of this matrix.
  Object get elements;

  /// Returns the imaginary part of this matrix
  AbstractDoubleVector imaginary();

  /// Fills the coordinates and values of cells having non-zero values into the
  /// specified lists. Fills into the lists, starting at index 0. After this
  /// call returns the specified lists all have a new size, the number of
  /// non-zero values.
  ///
  /// In general, fill order is unspecified.
  void nonZeros({List<int> indexList, final List<List<double>> valueList}) {
    bool fillIndexList = indexList != null;
    bool fillValueList = valueList != null;
    if (fillIndexList) {
      indexList.clear();
    }
    if (fillValueList) {
      valueList.clear();
    }
    indexList.clear();
    valueList.clear();
    for (int i = 0; i < size; i++) {
      Float64List value = get(i);
      if (value[0] != 0 || value[1] != 0) {
        if (fillIndexList) {
          indexList.add(i);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
    }
  }

  /// Returns the matrix cell value at coordinate [index].
  ///
  /// Provided with invalid parameters this method may return invalid objects
  /// without throwing any exception. You should only use this method when
  /// you are absolutely sure that the coordinate is within bounds.
  Float64List get(int index);

  /// Returns the real part of this matrix
  AbstractDoubleVector real();

  /// Construct and returns a new empty matrix of the same dynamic type
  /// as the receiver, having the same size.
  AbstractComplexVector like() {
    return like1D(size);
  }

  /// Construct and returns a new empty matrix of the same dynamic type
  /// as the receiver, having the specified size.
  AbstractComplexVector like1D(int size);

  /// Construct and returns a new 2-d matrix of the corresponding dynamic
  /// type, entirely independent of the receiver.
  AbstractComplexMatrix like2D(int rows, int columns);

  /// Returns new matrix of size rows x columns whose elements are
  /// taken column-wise from this matrix.
  AbstractComplexMatrix reshape(int rows, int columns);

  /// Sets the matrix cell at coordinate [index] to the specified value.
  void put(int index, double re, double im) {
    checkIndex(this, index);
    setParts(index, re, im);
  }

  /// Sets the matrix cell at coordinate [index] to the specified value.
  void operator []=(int index, Float64List value) {
    checkIndex(this, index);
    set(index, value);
  }

  /// Sets the matrix cell at coordinate [index] to the specified value.
  void setParts(int index, double re, double im);

  /// Sets the matrix cell at coordinate [index] to the specified value.
  void set(int index, Float64List value);

  /// Constructs and returns a 1-dimensional array containing the cell values.
  Float64List toList() {
    var values = new Float64List(2 * size);
    fillList(values);
    return values;
  }

  /// Fills the cell values into the specified 1-dimensional array.
  void fillList(final Float64List values) {
    if (values.length < 2 * size) {
      throw new ArgumentError("values too small");
    }
    for (int i = 0; i < size; i++) {
      var tmp = get(i);
      values[2 * i] = tmp[0];
      values[2 * i + 1] = tmp[1];
    }
  }

  /// Returns a string representation using default formatting ("%.4f").
  String toString() => toStringFormat('0.00E+00' /*"%.4f"*/);

  /// Returns a string representation using given [format].
  String toStringFormat(String format) {
    var f = new NumberFormat(format);
    var s = new StringBuffer("ComplexVector: ${size} elements\n\n");
    Float64List elem = new Float64List(2);
    for (int i = 0; i < size; i++) {
      elem = get(i);
      /*if (elem[1] == 0) {
        s.write(f.format(elem[0]) + "\n");
        continue;
      }
      if (elem[0] == 0) {
        s.write(f.format(elem[1]) + "i\n");
        continue;
      }*/
      if (elem[1] < 0) {
        s.write(f.format(elem[0]) + " - " + f.format(-elem[1]) + "i\n");
        continue;
      }
      s.write(f.format(elem[0]) + " + " + f.format(elem[1]) + "i\n");
    }
    return s.toString();
  }

  /// Constructs and returns a new *flip view*. What used to be index
  /// `0` is now index `size()-1`.
  AbstractComplexVector flip() {
    var v = _view();
    vFlip(v);
    return v;
  }

  /// Constructs and returns a new *sub-range view* that is a [width]
  /// sub matrix starting at [index].
  AbstractComplexVector part(int index, int width) {
    var v = _view();
    vPart(v, index, width);
    return v;
  }

  /// Constructs and returns a new *selection view* that is a matrix
  /// holding the indicated cells. Indexes can occur multiple times
  /// and can be in arbitrary order.
  AbstractComplexVector select(Int32List indexes) {
    // check for "all"
    if (indexes == null) {
      indexes = new Int32List(size);
      for (int i = size - 1; --i >= 0;) {
        indexes[i] = i;
      }
    }

    checkIndexes(this, indexes);
    var offsets = new Int32List(indexes.length);
    for (int i = 0; i < indexes.length; i++) {
      offsets[i] = index(indexes[i]);
    }
    return _viewSelectionLike(offsets);
  }

  /// Constructs and returns a new *stride view* which is a sub matrix
  /// consisting of every i-th cell.
  AbstractComplexVector strides(int stride) {
    var v = _view();
    vStride(v, stride);
    return v;
  }

  /// Returns the dot product of two vectors x and y. Operates on cells at
  /// indexes `from .. Min(size(),y.size(),from+length)-1`.
  Float64List dot(final AbstractComplexVector y,
      [final int from = 0, int length = null]) {
    if (length == null) {
      length = this.size;
    }
    if (from < 0 || length <= 0) {
      return new Float64List.fromList([0.0, 0.0]);
    }

    int tail = from + length;
    if (size < tail) {
      tail = size;
    }
    if (y.size < tail) {
      tail = y.size;
    }
    length = tail - from;
    var sum = new Float64List(2);

    for (int k = 0; k < length; k++) {
      var idx = k + from;
      var tmp = y.get(idx);
      tmp[1] = -tmp[1]; // complex conjugate
      sum = Complex.plus(sum, Complex.multiply(tmp, get(idx)));
    }
    return sum;
  }

  /// Returns the sum of all cells.
  Float64List sum() {
    var sum = new Float64List(2);
    for (int k = 0; k < size; k++) {
      sum = Complex.plus(sum, get(k));
    }
    return sum;
  }

  /// Returns the content of this matrix if it is a wrapper; or `this`
  /// otherwise. Override this method in wrappers.
  AbstractComplexVector _getContent() => this;

  /// Returns `true` if both matrices share at least one identical cell.
  bool _haveSharedCells(AbstractComplexVector other) {
    if (other == null) {
      return false;
    }
    if (this == other) {
      return true;
    }
    return _getContent()._haveSharedCellsRaw(other._getContent());
  }

  /// Always returns false
  bool _haveSharedCellsRaw(AbstractComplexVector other) => false;

  /// Constructs and returns a new view equal to the receiver.
  AbstractComplexVector _view() {
    return clone() as AbstractComplexVector;
  }

  /// Construct and returns a new selection view.
  AbstractComplexVector _viewSelectionLike(Int32List offsets);

  Object clone();

  AbstractComplexVector operator *(AbstractComplexVector y) {
    return copy()..assign(y, cfunc.mult);
  }

  AbstractComplexVector operator /(AbstractComplexVector y) {
    return copy()..assign(y, cfunc.div);
  }

  AbstractComplexVector operator +(AbstractComplexVector y) {
    return copy()..assign(y, cfunc.plus);
  }

  AbstractComplexVector operator -(AbstractComplexVector y) {
    return copy()..assign(y, cfunc.minus);
  }

  AbstractComplexVector operator -() => copy()..apply(cfunc.neg);

  AbstractComplexVector conj() => copy()..apply(cfunc.conj);

  AbstractComplexVector abs() => copy()..applyReal(cfunc.abs);

  AbstractComplexVector arg() => copy()..applyReal(cfunc.arg);
}
