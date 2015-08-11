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

/// Abstract base class for 1-d matrices (aka vectors) holding [int]
/// elements.
abstract class IntVector extends AbstractVector /*with ListMixin<int>*/ {
  IntVector._(int size, [int zero = 0, int stride = 1, bool isNoView = true])
      : super(size, zero, stride, isNoView);

  factory IntVector(int size) = DenseIntVector;

  factory IntVector.sparse(int size) = SparseIntVector;

  /// Applies a function to each cell and aggregates the results. Returns a
  /// value `v` such that `v==a(size())` where
  /// `a(i) == aggr( a(i-1), f(get(i)) )` and terminators are
  /// `a(1) == f(get(0)), a(0) == double.NaN`.
  ///
  /// Example:
  ///     matrix = 0 1 2 3
  ///
  ///     // Sum( x[i]*x[i] )
  ///     matrix.aggregate(F.plus,F.square);
  ///     --> 14
  int aggregate(ifunc.IntIntFunction aggr, ifunc.IntFunction fn) {
    if (size == 0) {
      throw new ArgumentError("size == 0");
    }
    int a = fn(get(0));
    for (int i = 1; i < size; i++) {
      a = aggr(a, fn(get(i)));
    }
    return a;
  }

  /// Applies the result of a function to each cell.
  ///
  /// Example:
  ///
  ///     // change each cell to its sine
  ///     matrix =   0.5      1.5      2.5       3.5
  ///     matrix.apply(F.sin);
  ///     -->
  ///     matrix ==  0.479426 0.997495 0.598472 -0.350783
  void apply(ifunc.IntFunction fn) {
    for (int i = 0; i < size; i++) {
      set(i, fn(get(i)));
    }
  }

  /// Sets all cells to the state specified by [value].
  void fill(int value) {
    for (int i = 0; i < size; i++) {
      set(i, value);
    }
  }

  /// Sets all cells to the state specified by [values]. [values]
  /// is required to have the same number of cells as the receiver.
  void setAll(Int32List values) {
    if (values.length != size) {
      throw new ArgumentError(
          "Must have same number of cells: length=${values.length} size()=$size");
    }
    for (int i = 0; i < size; i++) {
      set(i, values[i]);
    }
  }

  /// Replaces all cell values of the receiver with the values of another
  /// matrix. Both matrices must have the same size. If both matrices share the
  /// same cells (as is the case if they are views derived from the same
  /// matrix) and intersect in an ambiguous way, then replaces *as if*
  /// using an intermediate auxiliary deep copy of [other].
  void copyFrom(IntVector other) {
    if (other == this) {
      return;
    }
    checkSize(this, other);
    IntVector other_loc;
    if (_haveSharedCells(other)) {
      other_loc = other.copy();
    } else {
      other_loc = other;
    }
    for (int i = 0; i < size; i++) {
      set(i, other_loc.get(i));
    }
  }

  /// Assigns the result of a function to each cell.
  ///
  /// Example:
  ///
  ///     // assign x[i] = x[i]^y[i]
  ///     m1 = 0 1 2 3;
  ///     m2 = 0 2 4 6;
  ///     m1.assign(m2, F.pow);
  ///     -->
  ///     m1 == 1 1 16 729
  void assign(IntVector y, ifunc.IntIntFunction fn) {
    checkSize(this, y);
    for (int i = 0; i < size; i++) {
      set(i, fn(get(i), y.get(i)));
    }
  }

  /// Returns the number of cells having non-zero values; ignores tolerance.
  int get cardinality {
    int cardinality = 0;
    for (int i = 0; i < size; i++) {
      if (get(i) != 0) cardinality++;
    }
    return cardinality;
  }

  /// Constructs and returns a deep copy of the receiver.
  IntVector copy() {
    IntVector copy = like();
    copy.copyFrom(this);
    return copy;
  }

  dynamic get elements;

  /// Compares this object against the specified object. The result is
  /// `true` if and only if the argument is not `null` and is at least
  /// a `IntVector` object that has the same sizes as the receiver and
  /// has exactly the same values at the same indexes.
  bool equals(IntVector obj) {
    if (identical(this, obj)) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    if (!(obj is IntVector)) {
      return false;
    }

    return iprop.equalsVector(this, obj);
  }

  /// Returns the matrix cell value at coordinate `index`.
  int operator [](int index) {
    checkIndex(this, index);
    return get(index);
  }

  /// Returns the content of this matrix if it is a wrapper; or `this`
  /// otherwise. Override this method in wrappers.
  IntVector _getContent() => this;

  /// Fills the coordinates and values of cells having negative values into
  /// the specified lists. Fills into the lists, starting at index 0. After
  /// this call returns the specified lists all have a new size, the number
  /// of negative values.
  void negative({List<int> indexList, List<int> valueList}) {
    bool fillIndexList = indexList != null;
    bool fillValueList = valueList != null;
    if (fillIndexList) {
      indexList.clear();
    }
    if (fillValueList) {
      valueList.clear();
    }
    int rem = size % 2;
    if (rem == 1) {
      int value = get(0);
      if (value < 0) {
        if (fillIndexList) {
          indexList.add(0);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
    }

    for (int i = rem; i < size; i += 2) {
      int value = get(i);
      if (value < 0) {
        if (fillIndexList) {
          indexList.add(i);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      value = get(i + 1);
      if (value < 0) {
        if (fillIndexList) {
          indexList.add(i + 1);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
    }
  }

  /// Fills the coordinates and values of cells having non-zero values into
  /// the specified lists. Fills into the lists, starting at index 0. After
  /// this call returns the specified lists all have a new size, the number
  /// of non-zero values.
  ///
  /// In general, fill order is unspecified.
  ///
  /// Example:
  ///
  ///     0, 0, 8, 0, 7
  ///     -->
  ///     indexList  = (2, 4)
  ///     valueList  = (8, 7)
  ///
  /// In other words, `get(2) == 8, get(4) == 7`.
  void nonzero({List<int> indexList, List<int> valueList}) {
    bool fillIndexList = indexList != null;
    bool fillValueList = valueList != null;
    if (fillIndexList) {
      indexList.clear();
    }
    if (fillValueList) {
      valueList.clear();
    }
    int rem = size % 2;
    if (rem == 1) {
      int value = get(0);
      if (value != 0) {
        if (fillIndexList) {
          indexList.add(0);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
    }

    for (int i = rem; i < size; i += 2) {
      int value = get(i);
      if (value != 0) {
        if (fillIndexList) {
          indexList.add(i);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      value = get(i + 1);
      if (value != 0) {
        if (fillIndexList) {
          indexList.add(i + 1);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
    }
  }

  /// Fills the coordinates and values of cells having positive values into
  /// the specified lists. Fills into the lists, starting at index 0. After
  /// this call returns the specified lists all have a new size, the number
  /// of positive values.
  void positive({List<int> indexList, List<int> valueList}) {
    bool fillIndexList = indexList != null;
    bool fillValueList = valueList != null;
    if (fillIndexList) {
      indexList.clear();
    }
    if (fillValueList) {
      valueList.clear();
    }
    int rem = size % 2;
    if (rem == 1) {
      int value = get(0);
      if (value > 0) {
        if (fillIndexList) {
          indexList.add(0);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
    }

    for (int i = rem; i < size; i += 2) {
      int value = get(i);
      if (value > 0) {
        if (fillIndexList) {
          indexList.add(i);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      value = get(i + 1);
      if (value > 0) {
        if (fillIndexList) {
          indexList.add(i + 1);
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
  int get(int index);

  /// Returns `true` if both matrices share at least one identical cell.
  bool _haveSharedCells(IntVector other) {
    if (other == null) {
      return false;
    }
    if (this == other) {
      return true;
    }
    return _getContent()._haveSharedCellsRaw(other._getContent());
  }

  /// Returns `true` if both matrices share at least one identical cell.
  bool _haveSharedCellsRaw(IntVector other) => false;

  /// Construct and returns a new empty matrix of the same dynamic type
  /// as the receiver, having the same size.
  IntVector like() => like1D(size);

  /// Construct and returns a new empty matrix of the same dynamic type
  /// as the receiver, having the specified size.
  IntVector like1D(int size);

  /// Construct and returns a new 2-d matrix of the corresponding dynamic
  /// type, entirelly independent of the receiver.
  IntMatrix like2D(int rows, int columns);

  /// Return the maximum value of this matrix together with its location.
  IntVectorLocation max() {
    int location = 0;
    int maxValue = get(location);
    int elem;
    for (int i = 1; i < size; i++) {
      elem = get(i);
      if (maxValue < elem) {
        maxValue = elem;
        location = i;
      }
    }
    return new IntVectorLocation._(maxValue, location);
  }

  /// Return the minimum value of this matrix together with its location.
  IntVectorLocation min() {
    int location = 0;
    int minValue = get(location);
    int elem;
    for (int i = 1; i < size; i++) {
      elem = get(i);
      if (minValue > elem) {
        minValue = elem;
        location = i;
      }
    }
    return new IntVectorLocation._(minValue, location);
  }

  /// Returns new IntMatrix of size rows x columns whose elements are taken
  /// column-wise from this matrix.
  IntMatrix reshape(int rows, int columns);

  /// Sets the matrix cell at coordinate `index` to the specified value.
  void operator []=(int index, int value) {
    checkIndex(this, index);
    set(index, value);
  }

  /// Sets the matrix cell at coordinate [index] to the specified [value].
  ///
  /// Provided with invalid parameters this method may access illegal indexes
  /// without throwing any exception. You should only use this method when
  /// you are absolutely sure that the coordinate is within bounds.
  void set(int index, int value);

  /// Constructs and returns a 1-dimensional array containing the cell values.
  Int32List toList() {
    var result = new Int32List(size);
    fillList(result);
    return result;
  }

  /// Fills the cell values into the specified 1-dimensional array.
  void fillList(Int32List values) {
    if (values.length < size) {
      throw new ArgumentError("values too small");
    }
    for (int i = 0; i < size; i++) {
      values[i] = get(i);
    }
  }

  /// Returns a string representation using default formatting.
  String toString() {
    //return new IntFormatter().toString1D(this);
    return elements.toString();
  }

  /// Constructs and returns a new view equal to the receiver.
  IntVector _view() => clone() as IntVector;

  /// Constructs and returns a new *flip view*. What used to be index
  /// `0` is now index `size-1`.
  IntVector flip() {
    var v = _view();
    vFlip(v);
    return v;
  }

  /// Constructs and returns a new *sub-range view* that is a [width]
  /// sub matrix starting at [index].
  IntVector part(int index, int width) {
    var v = _view();
    vPart(v, index, width);
    return v;
  }

  /// Constructs and returns a new *selection view* that is a matrix
  /// holding the indicated cells. There holds
  /// `view.size() == indexes.length` and
  /// `view.get(i) == this.get(indexes[i])`. Indexes can occur multiple
  /// times and can be in arbitrary order.
  ///
  /// Example:
  ///
  ///     this     = (0,0,8,0,7)
  ///     indexes  = (0,2,4,2)
  ///     -->
  ///     view     = (0,8,7,8)
  IntVector select(Int32List indexes) {
    // check for "all"
    if (indexes == null) {
      indexes = new Int32List(size);
      for (int i = 0; i < size; i++) {
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

  /// Construct and returns a new selection view.
  IntVector _viewSelectionLike(Int32List offsets);

  /// Constructs and returns a new *stride view* which is a sub matrix
  /// consisting of every i-th cell.
  IntVector strides(int stride) {
    var v = _view();
    vStride(v, stride);
    return v;
  }

  /// Returns the dot product of two vectors x and y, which is
  /// `Sum(x[i]*y[i])`. Where `x == this`. Operates on cells at
  /// indexes `from .. Min(size(),y.size(),from+length)-1`.
  int dot(IntVector y, [int from = 0, int length = null]) {
    if (length == null) {
      length = size;
    }
    if (from < 0 || length <= 0) {
      return 0;
    }

    int tail = from + length;
    if (size < tail) {
      tail = size;
    }
    if (y.size < tail) {
      tail = y.size;
    }
    length = tail - from;

    int sum = 0;
    int i = tail - 1;
    for (int k = length; --k >= 0; i--) {
      sum += get(i) * y.get(i);
    }
    return sum;
  }

  /// Returns the sum of all cells; `Sum( x[i] )`.
  int sum() {
    if (size == 0) {
      return 0;
    }
    return aggregate(ifunc.plus, ifunc.identity);
  }

  Object clone();

  IntVector operator &(IntVector a) {
    return copy()..assign(a, ifunc.and);
  }

  IntVector operator |(IntVector a) {
    return copy()..assign(a, ifunc.or);
  }

  IntVector operator ~() {
    return copy()..apply(ifunc.not);
  }
}
