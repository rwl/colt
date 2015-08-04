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

/// Abstract base class for 2-d matrices holding [int] elements.
abstract class AbstractIntMatrix extends AbstractMatrix {
  AbstractIntMatrix(int rows, int columns, [int rowZero = 0, int columnZero = 0,
      int rowStride = null, int columnStride = 1, bool isNoView = true])
      : super(rows, columns, rowZero, columnZero, rowStride, columnStride,
          isNoView);

  /// Applies a function to each cell and aggregates the results. Returns a
  /// value `v` such that `v==a(size())` where
  /// `a(i) == aggr( a(i-1), f(get(row,column)) )` and terminators are
  /// `a(1) == f(get(0,0)), a(0)==Int.NaN`.
  ///
  /// Example:
  ///
  ///     2 x 2 matrix
  ///     0 1
  ///     2 3
  ///
  ///     // Sum( x[row,col]*x[row,col] )
  ///     matrix.aggregate(F.plus, F.square);
  ///     --> 14
  int aggregate(final ifunc.IntIntFunction aggr, final ifunc.IntFunction fn) {
    if (size == 0) {
      throw new ArgumentError("size == 0");
    }
    int a = fn(get(0, 0));
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
  ///     matrix.assign(F.sin);
  ///     -->
  ///     2 x 2 matrix
  ///     0.479426  0.997495
  ///     0.598472 -0.350783
  void apply(final ifunc.IntFunction fn) {
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        set(r, c, fn(get(r, c)));
      }
    }
  }

  /// Sets all cells to the state specified by [value].
  void fill(final int value) {
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        set(r, c, value);
      }
    }
  }

  /// Sets all cells to the state specified by [values]. [values]
  /// is required to have the form `values[row*column]` and elements
  /// have to be stored in a row-wise order.
  void setAll(final Int32List values) {
    if (values.length != rows * columns) {
      throw new ArgumentError("Must have same length: length=${values.length} "
          "rows()*columns()=${rows * columns}");
    }
    int idx = 0;
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        set(r, c, values[idx++]);
      }
    }
  }

  /// Replaces all cell values of the receiver with the values of another
  /// matrix. Both matrices must have the same number of rows and columns. If
  /// both matrices share the same cells (as is the case if they are views
  /// derived from the same matrix) and intersect in an ambiguous way, then
  /// replaces as if using an intermediate auxiliary deep copy of [other].
  void copyFrom(AbstractIntMatrix other) {
    if (other == this) {
      return;
    }
    checkShape(this, other);
    AbstractIntMatrix other_loc;
    if (_haveSharedCells(other)) {
      other_loc = other.copy();
    } else {
      other_loc = other;
    }
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        set(r, c, other_loc.get(r, c));
      }
    }
  }

  /// Assigns the result of a function to each cell;
  /// `x[row,col] = function(x[row,col],y[row,col])`.
  ///
  /// Example:
  ///
  ///     // assign x[row,col] = x[row,col]^y[row,col]
  ///     m1 = 2 x 2 matrix
  ///     0 1
  ///     2 3
  ///
  ///     m2 = 2 x 2 matrix
  ///     0 2
  ///     4 6
  ///
  ///     m1.assign(m2, F.pow);
  ///     -->
  ///     m1 == 2 x 2 matrix
  ///     1   1
  ///     16 729
  void assign(final AbstractIntMatrix y, final ifunc.IntIntFunction fn) {
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
  AbstractIntMatrix copy() => like()..copyFrom(this);

  Object get elements;

  /// Compares this object against the specified object. The result is
  /// `true` if and only if the argument is not `null` and is at least
  /// a [AbstractIntMatrix] object that has the same number of columns
  /// and rows as the receiver and has exactly the same values at the
  /// same coordinates.
  bool equals(AbstractIntMatrix obj) {
    if (identical(this, obj)) {
      return true;
    }
    if (obj == null) {
      return false;
    }

    return iprop.equalsMatrix(this, obj);
  }

  /// Assigns the result of a function to each non-zero cell;
  /// `x[row,col] = function(x[row,col])`. Use this method for fast
  /// special-purpose iteration.
  void forEachNonZero(final ifunc.IntIntIntFunction fn) {
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        int value = get(r, c);
        if (value != 0) {
          int a = fn(r, c, value);
          if (a != value) set(r, c, a);
        }
      }
    }
  }

  /// Returns the matrix cell value at coordinate `[row,column]`.
  int at(int row, int column) {
    if (column < 0 || column >= columns || row < 0 || row >= rows) {
      throw new RangeError("row:$row, column:$column");
    }
    return get(row, column);
  }

  /// Returns the content of this matrix if it is a wrapper; or `this`
  /// otherwise. Override this method in wrappers.
  AbstractIntMatrix _getContent() => this;

  /// Fills the coordinates and values of cells having negative values into
  /// the specified lists. Fills into the lists, starting at index 0. After
  /// this call returns the specified lists all have a new size, the number
  /// of negative values.
  void negative(final Int32List rowList, final Int32List columnList,
      final Int32List valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        int value = get(r, c);
        if (value < 0) {
          rowList.add(r);
          columnList.add(c);
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
  ///     2 x 3 matrix:
  ///     0, 0, 8
  ///     0, 7, 0
  ///     -->
  ///     rowList    = (0,1)
  ///     columnList = (2,1)
  ///     valueList  = (8,7)
  ///
  /// In other words, `get(0,2) == 8, get(1,1) == 7`.
  void nonzero(final Int32List rowList, final Int32List columnList,
      final Int32List valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        int value = get(r, c);
        if (value != 0) {
          rowList.add(r);
          columnList.add(c);
          valueList.add(value);
        }
      }
    }
  }

  /// Fills the coordinates and values of cells having positive values into
  /// the specified lists. Fills into the lists, starting at index 0. After
  /// this call returns the specified lists all have a new size, the number
  /// of positive values.
  void positive(final Int32List rowList, final Int32List columnList,
      final Int32List valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        int value = get(r, c);
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
  int get(int row, int column);

  /// Returns `true` if both matrices share at least one identical cell.
  bool _haveSharedCells(AbstractIntMatrix other) {
    if (other == null) {
      return false;
    }
    if (this == other) {
      return true;
    }
    return _getContent()._haveSharedCellsRaw(other._getContent());
  }

  /// Returns `true` if both matrices share at least one identical cell.
  bool _haveSharedCellsRaw(AbstractIntMatrix other) => false;

  /// Construct and returns a new empty matrix of the same dynamic type
  /// as the receiver, having the same number of rows and columns.
  AbstractIntMatrix like() => like2D(rows, columns);

  /// Construct and returns a new empty matrix of the same dynamic type
  /// as the receiver, having the specified number of rows and columns.
  AbstractIntMatrix like2D(int rows, int columns);

  /// Construct and returns a new 1-d matrix of the corresponding dynamic
  /// type, entirelly independent of the receiver.
  AbstractIntVector like1D(int size);

  /// Construct and returns a new 1-d matrix of the corresponding dynamic
  /// type, sharing the same cells.
  AbstractIntVector _like1D(int size, int zero, int stride);

  /// Return the maximum value of this matrix together with its location.
  IntMatrixLocation max() {
    int rowLocation = 0;
    int columnLocation = 0;
    int maxValue = get(0, 0);
    int elem;
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
    return new IntMatrixLocation._(maxValue, rowLocation, columnLocation);
  }

  /// Return the minimum value of this matrix together with its location.
  IntMatrixLocation min() {
    int rowLocation = 0;
    int columnLocation = 0;
    int minValue = get(0, 0);
    int elem;
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
    return new IntMatrixLocation._(minValue, rowLocation, columnLocation);
  }

  /// Sets the matrix cell at coordinate `[row,column]` to the specified
  /// value.
  void put(int row, int column, int value) {
    if (column < 0 || column >= columns || row < 0 || row >= rows) {
      throw new RangeError("row:$row, column:$column");
    }
    set(row, column, value);
  }

  /// Sets the matrix cell at coordinate `[row,column]` to the specified
  /// value.
  void set(int row, int column, int value);

  /// Returns a string representation using default formatting.
  String toString() => new IntFormatter().toString2D(this);

  /// Constructs and returns a new view equal to the receiver. The view is a
  /// shallow clone.
  AbstractIntMatrix _view() => clone() as AbstractIntMatrix;

  /// Constructs and returns a new *slice view* representing the rows of
  /// the given column.
  ///
  /// Example:
  ///     2 x 3 matrix:
  ///     1, 2, 3
  ///     4, 5, 6
  ///
  ///     viewColumn(0) ==>
  ///     Vector of size 2:
  ///     1, 4
  AbstractIntVector column(int column) {
    checkColumn(this, column);
    int viewSize = rows;
    int viewZero = index(0, column);
    int viewStride = rowStride;
    return _like1D(viewSize, viewZero, viewStride);
  }

  /// Constructs and returns a new *flip view* along the column axis. What
  /// used to be column `0` is now column `columns()-1`
  ///
  /// Example:
  ///     2 x 3 matrix:
  ///     1, 2, 3
  ///     4, 5, 6
  ///
  ///     columnFlip ==>
  ///     2 x 3 matrix:
  ///     3, 2, 1
  ///     6, 5, 4
  ///
  ///     columnFlip ==>
  ///     2 x 3 matrix:
  ///     1, 2, 3
  ///     4, 5, 6
  AbstractIntMatrix columnFlip() {
    var v = _view();
    vColumnFlip(v);
    return v;
  }

  /// Constructs and returns a new *dice (transposition) view*; Swaps
  /// axes; example: 3 x 4 matrix --> 4 x 3 matrix.
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
  AbstractIntMatrix dice() {
    var v = _view();
    vDice(v);
    return v;
  }

  /// Constructs and returns a new *sub-range view* that is a
  /// `height x width` sub matrix starting at `[row,column]`.
  AbstractIntMatrix part(int row, int column, int height, int width) {
    var v = _view();
    vBox(v, row, column, height, width);
    return v;
  }

  /// Constructs and returns a new *slice view* representing the columns
  /// of the given row.
  ///
  /// Example:
  ///     2 x 3 matrix:
  ///     1, 2, 3
  ///     4, 5, 6
  ///     row(0) ==>
  ///     Vector of size 3:
  ///     1, 2, 3
  AbstractIntVector row(int row) {
    checkRow(this, row);
    int viewSize = columns;
    int viewZero = index(row, 0);
    int viewStride = columnStride;
    return _like1D(viewSize, viewZero, viewStride);
  }

  /// Constructs and returns a new *flip view* aint the row axis. What
  /// used to be row `0` is now row `rows()-1`.
  ///
  /// Example:
  ///     2 x 3 matrix:
  ///     1, 2, 3
  ///     4, 5, 6
  ///     rowFlip ==>
  ///     2 x 3 matrix:
  ///     4, 5, 6
  ///     1, 2, 3
  ///     rowFlip ==>
  ///     2 x 3 matrix:
  ///     1, 2, 3
  ///     4, 5, 6
  AbstractIntMatrix rowFlip() {
    var v = _view();
    vRowFlip(v);
    return v;
  }

  /// Constructs and returns a new *selection view* that is a matrix
  /// holding the indicated cells. Indexes can occur multiple times and
  /// can be in arbitrary order.
  ///
  /// Example:
  ///     2 x 3 matrix:
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
  /// parameter
  AbstractIntMatrix select(Int32List rowIndexes, Int32List columnIndexes) {
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

  /// Construct and returns a new selection view.
  AbstractIntMatrix _viewSelectionLike(
      Int32List rowOffsets, Int32List columnOffsets);

  /// Constructs and returns a new *stride view* which is a sub matrix
  /// consisting of every i-th cell.
  AbstractIntMatrix strides(int rowStride, int columnStride) {
    var v = _view();
    vStrides(v, rowStride, columnStride);
    return v;
  }

  /// Linear algebraic matrix-vector multiplication;
  /// `z = alpha * A * y + beta*z`.
  /// `z[i] = alpha*Sum(A[i,j] * y[j]) + beta*z[i], i=0..A.rows()-1, j=0..y.length-1`.
  /// Where `A == this`.
  AbstractIntVector mult(final AbstractIntVector y, [AbstractIntVector z = null,
      final int alpha = 1, int beta = null, final bool transposeA = false]) {
    if (beta == null) {
      beta = z == null ? 1 : 0;
    }
    if (transposeA) {
      return dice().mult(y, z, alpha, beta, false);
    }
    AbstractIntVector z_loc;
    if (z == null) {
      z_loc = new IntVector(rows);
    } else {
      z_loc = z;
    }
    if (columns != y.size || rows > z_loc.size) {
      throw new ArgumentError("Incompatible args: " +
          toStringShort() +
          ", " +
          y.toStringShort() +
          ", " +
          z_loc.toStringShort());
    }

    for (int r = 0; r < rows; r++) {
      int s = 0;
      for (int c = 0; c < columns; c++) {
        s += get(r, c) * y.get(c);
      }
      z_loc.set(r, alpha * s + beta * z_loc.get(r));
    }
    return z_loc;
  }

  /// Linear algebraic matrix-matrix multiplication;
  /// `C = alpha * A x B + beta*C`.
  /// `C[i,j] = alpha*Sum(A[i,k] * B[k,j]) + beta*C[i,j], k=0..n-1`.
  /// Matrix shapes: `A(m x n), B(n x p), C(m x p)`.
  AbstractIntMatrix multiply(final AbstractIntMatrix B,
      [AbstractIntMatrix C = null, final int alpha = 1, int beta = null,
      final bool transposeA = false, final bool transposeB = false]) {
    if (beta == null) {
      beta = C == null ? 1 : 0;
    }
    if (transposeA) {
      return dice().multiply(B, C, alpha, beta, false, transposeB);
    }
    if (transposeB) {
      return multiply(B.dice(), C, alpha, beta, transposeA, false);
    }

    int m = rows;
    int n = columns;
    int p = B.columns;
    AbstractIntMatrix C_loc;
    if (C == null) {
      C_loc = new IntMatrix(m, p);
    } else {
      C_loc = C;
    }
    if (B.rows != n) {
      throw new ArgumentError("Matrix inner dimensions must agree:" +
          toStringShort() +
          ", " +
          B.toStringShort());
    }
    if (C_loc.rows != m || C_loc.columns != p) {
      throw new ArgumentError("Incompatibe result matrix: " +
          toStringShort() +
          ", " +
          B.toStringShort() +
          ", " +
          C_loc.toStringShort());
    }
    if (this == C_loc || B == C_loc) {
      throw new ArgumentError("Matrices must not be identical");
    }

    for (int a = 0; a < p; a++) {
      for (int b = 0; b < m; b++) {
        int s = 0;
        for (int c = 0; c < n; c++) {
          s += get(b, c) * B.get(c, a);
        }
        C_loc.set(b, a, alpha * s + beta * C_loc.get(b, a));
      }
    }

    return C_loc;
  }

  /// Returns the sum of all cells; `Sum( x[i,j] )`.
  int sum() {
    if (size == 0) {
      return 0;
    }
    return aggregate(ifunc.plus, ifunc.identity);
  }

  Object clone();
}
