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

/// Abstract base class for 2-d matrices holding [double] elements.
///
/// A matrix has a number of rows and columns, which are assigned upon instance
/// construction - The matrix's size is then `rows*columns`. Elements
/// are accessed via `[row,column]` coordinates. Legal coordinates range
/// from `[0,0]` to `[rows-1,columns-1]`. Any attempt to access
/// an element at a coordinate outside this rangewill throw a [RangeError].
abstract class AbstractDoubleMatrix extends AbstractMatrix {
  AbstractDoubleMatrix(int rows, int columns, [int rowZero = 0,
      int columnZero = 0, int rowStride = null, int columnStride = 1,
      bool isNoView = true])
      : super(rows, columns, rowZero, columnZero, rowStride, columnStride,
          isNoView);

  /// Applies a function to each cell and aggregates the results. Returns a
  /// value `v` such that `v==a(size)` where
  /// `a(i) == aggr( a(i-1), f(get(row,column)) )` and terminators are
  /// `a(1) == f(get(0,0)), a(0)==Double.NaN`.
  ///
  /// Example:
  ///     2 x 2 matrix
  ///     0 1
  ///     2 3
  ///
  ///     // Sum( x[row,col]*x[row,col] )
  ///     matrix.aggregate(F.plus,F.square);
  ///     --> 14
  double aggregate(DoubleDoubleFunction aggr, DoubleFunction fn) {
    if (size == 0) {
      return double.NAN;
    }
    double a = fn(get(0, 0));
    int d = 1; // first cell already done
    for (int r = 0; r < rows; r++) {
      for (int c = d; c < columns; c++) {
        a = aggr(a, fn(get(r, c)));
      }
      d = 0;
    }
    return a;
  }

  /// Applies a function to each cell;
  /// `x[row,col] = function(x[row,col])`.
  ///
  /// Example:
  ///
  ///     matrix = 2 x 2 matrix
  ///     0.5 1.5
  ///     2.5 3.5
  ///
  ///     // change each cell to its sine
  ///     matrix.assign(cern.jet.math.Functions.sin);
  ///     -->
  ///     2 x 2 matrix
  ///     0.479426  0.997495
  ///     0.598472 -0.350783
  void apply(final DoubleFunction fn) {
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        set(r, c, fn(get(r, c)));
      }
    }
  }

  /// Sets all cells to the state specified by `value`.
  void fill(final double value) {
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        set(r, c, value.toDouble());
      }
    }
  }

  /// Sets all cells to the state specified by [values]. [values] is required
  /// to have the form `values[row*column]` and elements have to be stored in
  /// a row-wise order.
  void setAll(Float64List values) {
    if (values.length != rows * columns) {
      throw new ArgumentError(
          "Must have same length: length=${values.length}rows()*columns()=${rows * columns}");
    }
    int idx = 0;
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        set(r, c, values[idx++]);
      }
    }
  }

  /// Replaces all cell values of the receiver with the values of another
  /// matrix. Both matrices must have the same number of rows and columns.
  /// If both matrices share the same cells (as is the case if they are views
  /// derived from the same matrix) and intersect in an ambiguous way, then
  /// replaces as if using an intermediate auxiliary deep copy of [other].
  void copyFrom(AbstractDoubleMatrix other) {
    if (other == this) {
      return;
    }
    checkShape(this, other);
    AbstractDoubleMatrix source;
    if (_haveSharedCells(other)) {
      source = other.copy();
    } else {
      source = other;
    }
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        set(r, c, source.get(r, c));
      }
    }
  }

  /// Assigns the result of a function to each cell;
  /// `x[row,col] = function(x[row,col],y[row,col])`.
  ///
  /// Example:
  ///
  ///      // assign x[row,col] = x[row,col]^y[row,col];
  ///      m1 = 2 x 2 matrix
  ///      0 1
  ///      2 3
  ///
  ///      m2 = 2 x 2 matrix
  ///      0 2
  ///      4 6
  ///
  ///      m1.assign(m2, F.pow);
  ///      -->
  ///      m1 == 2 x 2 matrix
  ///      1   1
  ///      16 729
  void assign(final AbstractDoubleMatrix y, final DoubleDoubleFunction fn) {
    checkShape(this, y);
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        set(r, c, fn(get(r, c), y.get(r, c)));
      }
    }
  }

  /// Returns the number of cells having non-zero values; ignores tolerance.
  int get cardinality {
    int cardinality = 0;
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        if (get(r, c) != 0) cardinality++;
      }
    }
    return cardinality;
  }

  /// Constructs and returns a deep copy of the receiver.
  AbstractDoubleMatrix copy() => like()..copyFrom(this);

  /// Returns the elements of this matrix.
  dynamic get elements;

  /// Compares this object against the specified object. The result is
  /// `true` if and only if the argument is not `null` and is at least
  /// a [AbstractDoubleMatrix] object that has the same number of columns
  /// and rows as the receiver and has exactly the same values at the same
  /// coordinates.
  bool equals(AbstractDoubleMatrix obj) {
    if (identical(this, obj)) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    return dprop.equalsMatrix(this, obj);
  }

  /// Applies a function to each `non-zero` cell;
  void forEachNonZero(final func.IntIntDoubleFunction fn) {
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        double value = get(r, c);
        if (value != 0) {
          double a = fn(r, c, value);
          if (a != value) {
            set(r, c, a);
          }
        }
      }
    }
  }

  /// Returns the matrix cell value at coordinate `[row,column]`.
  double at(int row, int column) {
    if (column < 0 || column >= columns || row < 0 || row >= rows) {
      throw new RangeError("row:$row, column:$column");
    }
    return get(row, column);
  }

  /// Return the maximum value of this matrix together with its location.
  DoubleMatrixLocation max() {
    int rowLocation = 0;
    int columnLocation = 0;
    double maxValue = get(0, 0);
    double elem;
    int d = 1;
    for (int r = 0; r < rows; r++) {
      for (int c = d; c < columns; c++) {
        elem = get(r, c);
        if (maxValue < elem) {
          maxValue = elem;
          rowLocation = r;
          columnLocation = c;
        }
      }
      d = 0;
    }
    return new DoubleMatrixLocation._(maxValue, rowLocation, columnLocation);
  }

  /// Return the minimum value of this matrix together with its location
  DoubleMatrixLocation min() {
    int rowLocation = 0;
    int columnLocation = 0;
    double minValue = get(0, 0);
    double elem;
    int d = 1;
    for (int r = 0; r < rows; r++) {
      for (int c = d; c < columns; c++) {
        elem = get(r, c);
        if (minValue > elem) {
          minValue = elem;
          rowLocation = r;
          columnLocation = c;
        }
      }
      d = 0;
    }
    return new DoubleMatrixLocation._(minValue, rowLocation, columnLocation);
  }

  /// Fills the coordinates and values of cells having negative values into the
  /// specified lists. Fills into the lists, starting at index 0. After this
  /// call returns the specified lists all have a new size, the number of
  /// negative values.
  void negativeValues(final List<int> rowList, final List<int> columnList,
      final List<double> valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        double value = get(r, c);
        if (value < 0) {
          rowList.add(r);
          columnList.add(c);
          valueList.add(value);
        }
      }
    }
  }

  /// Fills the coordinates and values of cells having non-zero values into the
  /// specified lists. Fills into the lists, starting at index 0. After this
  /// call returns the specified lists all have a new size, the number of
  /// non-zero values.
  ///
  /// In general, fill order is unspecified. For example:
  ///
  ///     2 x 3 matrix:
  ///     0, 0, 8
  ///     0, 7, 0
  ///     -->
  ///     rowList    = (0,1)
  ///     columnList = (2,1)
  ///     valueList  = (8,7)
  ///
  /// In other words, `get(0,2) == 8, get(1,1) == 7`.
  void nonzero(final List<int> rowList, final List<int> columnList,
      final List<double> valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        double value = get(r, c);
        if (value != 0) {
          rowList.add(r);
          columnList.add(c);
          valueList.add(value);
        }
      }
    }
  }

  /// Fills the coordinates and values of cells having positive values into the
  /// specified lists. Fills into the lists, starting at index 0. After this
  /// call returns the specified lists all have a new size, the number of
  /// positive values.
  void positive(final List<int> rowList, final List<int> columnList,
      final List<double> valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        double value = get(r, c);
        if (value > 0) {
          rowList.add(r);
          columnList.add(c);
          valueList.add(value);
        }
      }
    }
  }

  /// Returns the matrix cell value at coordinate `[row,column]`.
  ///
  /// Provided with invalid parameters this method may return invalid objects
  /// without throwing any exception. You should only use this method when
  /// you are absolutely sure that the coordinate is within bounds.
  double get(int row, int column);

  /// Construct and returns a new empty matrix of the same dynamic type
  /// as the receiver, having the same number of rows and columns.
  AbstractDoubleMatrix like() => like2D(rows, columns);

  /// Construct and returns a new empty matrix of the same dynamic type
  /// as the receiver, having the specified number of rows and columns.
  AbstractDoubleMatrix like2D(int rows, int columns);

  /// Construct and returns a new 1-d matrix of the corresponding dynamic
  /// type, entirelly independent of the receiver.
  AbstractDoubleVector like1D(int size);

  /// Sets the matrix cell at coordinate `[row,column]` to the specified
  /// value.
  void put(int row, int column, double value) {
    if (column < 0 || column >= columns || row < 0 || row >= rows) {
      throw new RangeError("row:$row, column:$column");
    }
    set(row, column, value);
  }

  /// Sets the matrix cell at coordinate `[row,column]` to the specified
  /// value.
  ///
  /// Provided with invalid parameters this method may access illegal indexes
  /// without throwing any exception. You should only use this method when
  /// you are absolutely sure that the coordinate is within bounds.
  void set(int row, int column, double value);

  /// Returns a string representation using default formatting.
  String toString() {
    return new DoubleFormatter().toStringDouble2D(this);
  }

  /// Constructs and returns a new *slice view* representing the rows of
  /// the given column.
  AbstractDoubleVector column(int column) {
    checkColumn(this, column);
    int viewSize = rows;
    int viewZero = index(0, column);
    int viewStride = rowStride;
    return _like1D(viewSize, viewZero, viewStride);
  }

  /// Constructs and returns a new *flip view* along the column axis. What
  /// used to be column `0` is now column `columns-1`
  AbstractDoubleMatrix columnFlip() {
    var v = _view();
    vColumnFlip(v);
    return v;
  }

  /// Constructs and returns a new *dice (transposition) view*; Swaps
  /// axes; example: 3 x 4 matrix --> 4 x 3 matrix. This is a zero-copy
  /// transposition, taking O(1), i.e. constant time.
  ///
  /// Example:
  ///     2 x 3 matrix:
  ///     1, 2, 3
  ///     4, 5, 6
  ///     transpose ==>
  ///     3 x 2 matrix:
  ///     1, 4
  ///     2, 5
  ///     3, 6
  ///     transpose ==>
  ///     2 x 3 matrix:
  ///     1, 2, 3
  ///     4, 5, 6
  AbstractDoubleMatrix dice() {
    var v = _view();
    vDice(v);
    return v;
  }

  /// Synonym for [dice].
  AbstractDoubleMatrix get T => dice();

  /// Constructs and returns a new *sub-range view* that is a
  /// `height x width` sub matrix starting at `[row,column]`.
  AbstractDoubleMatrix part(int row, int column, int height, int width) {
    var v = _view();
    vBox(v, row, column, height, width);
    return v;
  }

  /// Constructs and returns a new *slice view* representing the columns
  /// of the given row.
  AbstractDoubleVector row(int r) {
    checkRow(this, r);
    int viewSize = columns;
    int viewZero = index(r, 0);
    int viewStride = columnStride;
    return _like1D(viewSize, viewZero, viewStride);
  }

  // Synonymous with [row].
  AbstractDoubleVector operator [](int r) => row(r);

  /// Constructs and returns a new *flip view* along the row axis. What
  /// used to be row `0` is now row `rows-1`.
  AbstractDoubleMatrix rowFlip() {
    var v = _view();
    vRowFlip(v);
    return v;
  }

  /// Constructs and returns a new *selection view* that is a matrix
  /// holding the indicated cells. There holds
  /// `view.rows == rowIndexes.length, view.columns == columnIndexes.length`
  /// and `view.get(i,j) == this.get(rowIndexes[i],columnIndexes[j])`.
  /// Indexes can occur multiple times and can be in arbitrary order.
  ///
  /// Example:
  ///
  ///     this = 2 x 3 matrix:
  ///     1, 2, 3
  ///     4, 5, 6
  ///     rowIndexes     = (0,1)
  ///     columnIndexes  = (1,0,1,0)
  ///     -->
  ///     view = 2 x 4 matrix:
  ///     2, 1, 2, 1
  ///     5, 4, 5, 4
  ///
  /// To indicate "all" rows or "all columns", simply set the respective
  /// parameter to `null`.
  AbstractDoubleMatrix select(Int32List rowIndexes, Int32List columnIndexes) {
    // check for "all"
    if (rowIndexes == null) {
      rowIndexes = new Int32List(rows);
      for (int i = 0; i < rows; i++) {
        rowIndexes[i] = i;
      }
    }
    if (columnIndexes == null) {
      columnIndexes = new Int32List(columns);
      for (int i = 0; i < columns; i++) {
        columnIndexes[i] = i;
      }
    }

    checkRowIndexes(this, rowIndexes);
    checkColumnIndexes(this, columnIndexes);
    var rowOffsets = new Int32List(rowIndexes.length);
    var columnOffsets = new Int32List(columnIndexes.length);
    for (int i = 0; i < rowIndexes.length; i++) {
      rowOffsets[i] = rowOffset(rowRank(rowIndexes[i]));
    }
    for (int i = 0; i < columnIndexes.length; i++) {
      columnOffsets[i] = columnOffset(columnRank(columnIndexes[i]));
    }
    return _viewSelectionLike(rowOffsets, columnOffsets);
  }

  /// Constructs and returns a new *stride view* which is a sub matrix
  /// consisting of every i-th cell. More specifically, the view has
  /// `this.rows/rowStride` rows and `this.columns/columnStride` columns
  /// holding cells `this.get(i*rowStride,j*columnStride)` for all
  /// `i = 0..rows/rowStride - 1, j = 0..columns/columnStride - 1`.
  AbstractDoubleMatrix strides(int rowStride, int columnStride) {
    var v = _view();
    vStrides(v, rowStride, columnStride);
    return v;
  }

  /// Linear algebraic matrix-vector multiplication;
  /// `z = alpha * A * y + beta*z`.
  /// `z[i] = alpha*Sum(A[i,j] * y[j]) + beta*z[i], i=0..A.rows-1, j=0..y.size-1`.
  AbstractDoubleVector mult(final AbstractDoubleVector y,
      [AbstractDoubleVector z = null, double alpha = 1.0, double beta = 0.0,
      bool transposeA = false]) {
    if (transposeA) {
      return dice().mult(y, z, alpha, beta, false);
    }
    AbstractDoubleVector zz;
    if (z == null) {
      zz = y.like1D(rows);
    } else {
      zz = z;
    }
    if (columns != y.size || rows > zz.size) {
      throw new ArgumentError("Incompatible args: " +
          toStringShort() +
          ", " +
          y.toStringShort() +
          ", " +
          zz.toStringShort());
    }

    for (int r = 0; r < rows; r++) {
      double s = 0.0;
      for (int c = 0; c < columns; c++) {
        s += get(r, c) * y.get(c);
      }
      zz.set(r, alpha * s + beta * zz.get(r));
    }
    return zz;
  }

  /// Linear algebraic matrix-matrix multiplication;
  /// `C = alpha * A x B + beta*C`.
  /// `C[i,j] = alpha*Sum(A[i,k] * B[k,j]) + beta*C[i,j], k=0..n-1`.
  ///
  /// Matrix shapes: `A(m x n), B(n x p), C(m x p)`.
  AbstractDoubleMatrix multiply(final AbstractDoubleMatrix B,
      [AbstractDoubleMatrix C = null, double alpha = 1.0, double beta = 0.0,
      bool transposeA = false, bool transposeB = false]) {
    if (transposeA) {
      return dice().multiply(B, C, alpha, beta, false, transposeB);
    }
    if (transposeB) {
      return this.multiply(B.dice(), C, alpha, beta, transposeA, false);
    }

    var m = rows;
    var n = columns;
    var p = B.columns;
    AbstractDoubleMatrix CC;
    if (C == null) {
      CC = like2D(m, p);
    } else {
      CC = C;
    }
    if (B.rows != n) {
      throw new ArgumentError("Matrix inner dimensions must agree:" +
          toStringShort() +
          ", " +
          B.toStringShort());
    }
    if (CC.rows != m || CC.columns != p) {
      throw new ArgumentError("Incompatibe result matrix: " +
          toStringShort() +
          ", " +
          B.toStringShort() +
          ", " +
          CC.toStringShort());
    }
    if (this == CC || B == CC) {
      throw new ArgumentError("Matrices must not be identical");
    }
    for (int a = 0; a < p; a++) {
      for (int b = 0; b < m; b++) {
        double s = 0.0;
        for (int c = 0; c < n; c++) {
          s += get(b, c) * B.get(c, a);
        }
        CC.set(b, a, alpha * s + beta * CC.get(b, a));
      }
    }
    return CC;
  }

  /// Returns the sum of all cells; `Sum( x[i,j] )`.
  double sum() {
    if (size == 0) {
      return 0.0;
    }
    return aggregate(func.plus, func.identity);
  }

  /// Returns the content of this matrix if it is a wrapper; or `this`
  /// otherwise. Override this method in wrappers.
  AbstractDoubleMatrix _getContent() => this;

  /// Returns `true` if both matrices share at least one identical cell.
  bool _haveSharedCells(AbstractDoubleMatrix other) {
    if (other == null) {
      return false;
    }
    if (this == other) {
      return true;
    }
    return _getContent()._haveSharedCellsRaw(other._getContent());
  }

  /// Returns `true` if both matrices share at least one identical cell.
  bool _haveSharedCellsRaw(AbstractDoubleMatrix other) => false;

  /// Construct and returns a new 1-d matrix of the corresponding dynamic
  /// type, sharing the same cells.
  AbstractDoubleVector _like1D(int size, int zero, int stride);

  /// Constructs and returns a new view equal to the receiver. The view is a
  /// shallow clone.
  AbstractDoubleMatrix _view() => clone() as AbstractDoubleMatrix;

  /// Construct and returns a new selection view.
  AbstractDoubleMatrix _viewSelectionLike(
      Int32List rowOffsets, Int32List columnOffsets);

  Object clone();

  AbstractDoubleMatrix operator *(AbstractDoubleMatrix y) {
    return copy()..assign(y, func.mult);
  }

  AbstractDoubleMatrix operator /(AbstractDoubleMatrix y) {
    return copy()..assign(y, func.div);
  }

  AbstractDoubleMatrix operator +(AbstractDoubleMatrix y) {
    return copy()..assign(y, func.plus);
  }

  AbstractDoubleMatrix operator -(AbstractDoubleMatrix y) {
    return copy()..assign(y, func.minus);
  }

  AbstractDoubleMatrix operator -() {
    return copy()..apply(func.neg);
  }
}
