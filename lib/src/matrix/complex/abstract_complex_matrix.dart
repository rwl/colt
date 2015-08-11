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

/// Abstract base class for 2-d matrices holding [Complex] elements.
abstract class ComplexMatrix extends AbstractMatrix {
  ComplexMatrix._(int rows, int columns, [int rowZero = 0,
      int columnZero = 0, int rowStride = null, int columnStride = 1,
      bool isNoView = true])
      : super(rows, columns, rowZero, columnZero, rowStride, columnStride,
          isNoView);

  factory ComplexMatrix(int rows, int columns) = DenseComplexMatrix;

  factory ComplexMatrix.sparse(int rows, int columns) = SparseComplexMatrix;

  /// Applies a function to each cell and aggregates the results.
  Complex aggregate(cfunc.ComplexComplexComplexFunction aggr,
      cfunc.ComplexComplexFunction fn) {
    if (size == 0) {
      return Complex.NAN;
    }
    var a = fn(get(0, 0));
    int d = 1; // first cell already done
    for (int r = 0; r < rows; r++) {
      for (int c = d; c < columns; c++) {
        a = aggr(a, fn(get(r, c)));
      }
      d = 0;
    }
    return a;
  }

  /// Applies a function to each cell.
  void apply(cfunc.ComplexComplexFunction fn) {
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        set(r, c, fn(get(r, c)));
      }
    }
  }

  /// Assigns the result of a function to the real part of the receiver. The
  /// imaginary part of the receiver is reset to zero.
  void applyReal(cfunc.ComplexRealFunction fn) {
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        double re = fn(get(r, c));
        setParts(r, c, re, 0.0);
      }
    }
  }

  /// Replaces all cell values of the receiver with the values of another
  /// matrix. Both matrices must have the same number of rows and columns.
  /// If both matrices share the same cells (as is the case if they are views
  /// derived from the same matrix) and intersect in an ambiguous way, then
  /// replaces as if using an intermediate auxiliary deep copy of [other].
  void copyFrom(ComplexMatrix other) {
    if (other == this) {
      return;
    }
    checkShape(this, other);
    ComplexMatrix otherLoc;
    if (_haveSharedCells(other)) {
      otherLoc = other.copy();
    } else {
      otherLoc = other;
    }
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        set(r, c, otherLoc.get(r, c));
      }
    }
  }

  /// Assigns the result of a function to each cell.
  void assign(ComplexMatrix y,
      cfunc.ComplexComplexComplexFunction f) {
    checkShape(this, y);
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        set(r, c, f(get(r, c), y.get(r, c)));
      }
    }
  }

  /// Sets all cells to the state specified by [re] and [im].
  void fill(double re, double im) {
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        setParts(r, c, re, im);
      }
    }
  }

  /// Sets all cells to the state specified by [values]. [values] is required
  /// to have the form
  /// `re = values[row*rowStride+column*columnStride];`
  /// `im = values[row*rowStride+column*columnStride+1]`
  /// and have exactly the same number of rows and columns as the receiver.
  void setAll(Float64List values) {
    if (values.length != rows * 2 * columns) {
      throw new ArgumentError(
          "Must have same length: length=${values.length} rows()*2*columns()=${rows * 2 * columns}");
    }
    int idx = 0;
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        setParts(r, c, values[idx], values[idx + 1]);
        idx += 2;
      }
    }
  }

  /// Replaces imaginary part of the receiver with the values of another real
  /// matrix. The real part of the receiver remains unchanged. Both matrices
  /// must have the same size.
  void setImaginary(DoubleMatrix other) {
    checkShape(this, other);
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        double re = get(r, c).real;
        double im = other.get(r, c);
        setParts(r, c, re, im);
      }
    }
  }

  /// Replaces real part of the receiver with the values of another real
  /// matrix. The imaginary part of the receiver remains unchanged. Both
  /// matrices must have the same size.
  void setReal(DoubleMatrix other) {
    checkShape(this, other);
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        double re = other.get(r, c);
        double im = get(r, c).imaginary;
        setParts(r, c, re, im);
      }
    }
  }

  /// Returns the number of cells having non-zero values; ignores tolerance.
  int get cardinality {
    int cardinality = 0;
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        var tmp = get(r, c);
        if (tmp.real != 0 || tmp.imaginary != 0) {
          cardinality++;
        }
      }
    }
    return cardinality;
  }

  /// Constructs and returns a deep copy of the receiver.
  ComplexMatrix copy() => like()..copyFrom(this);

  /// Compares this object against the specified object. The result is
  /// `true` if and only if the argument is not `null` and is at least
  /// a [DoubleMatrix] object that has the same number of columns
  /// and rows as the receiver and has exactly the same values at the
  /// same coordinates.
  bool equals(ComplexMatrix obj) {
    if (identical(this, obj)) {
      return true;
    }
    if (obj == null) {
      return false;
    }

    return cprop.equalsMatrix(this, obj);
  }

  /// Assigns the result of a function to each non-zero cell. Use this
  /// method for fast special-purpose iteration.
  void forEachNonZero(cfunc.IntIntComplexFunction fn) {
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        var value = get(r, c);
        if (value.real != 0 || value.imaginary != 0) {
          var v = fn(r, c, value);
          set(r, c, v);
        }
      }
    }
  }

  /// Returns the matrix cell value at coordinate `[row,column]`.
  Complex at(int row, int column) {
    if (column < 0 || column >= columns || row < 0 || row >= rows) {
      throw new RangeError("row:$row, column:$column");
    }
    return get(row, column);
  }

  /// Returns a new matrix that is a complex conjugate of this matrix. If
  /// unconjugated complex transposition is needed, one should use [dice]
  /// method.
  ComplexMatrix conjugateTranspose() {
    ComplexMatrix transpose = dice().copy();
    for (int r = 0; r < columns; r++) {
      for (int c = 0; c < rows; c++) {
        var tmp = transpose.get(r, c);
        transpose.set(r, c, tmp.conjugate());
      }
    }
    return transpose;
  }

  // Synonym for [conjugateTranspose].
  ComplexMatrix get H => conjugateTranspose();

  Object get elements;

  /// Returns the imaginary part of this matrix
  DoubleMatrix imaginary();

  /// Fills the coordinates and values of cells having non-zero values into the
  /// specified lists. Fills into the lists, starting at index 0. After this
  /// call returns the specified lists all have a new size, the number of
  /// non-zero values.
  void nonzero(List<int> rowList, List<int> columnList,
      List<Complex> valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        var value = get(r, c);
        if (value.real != 0 || value.imaginary != 0) {
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
  Complex get(int row, int column);

  /// Returns the real part of this matrix.
  DoubleMatrix real();

  /// Construct and returns a new empty matrix of the same dynamic type
  /// as the receiver, having the same number of rows and columns.
  ComplexMatrix like() => like2D(rows, columns);

  /// Construct and returns a new empty matrix of the same dynamic type
  /// as the receiver, having the specified number of rows and columns.
  ComplexMatrix like2D(int rows, int columns);

  /// Construct and returns a new 1-d matrix of the corresponding dynamic
  /// type, entirelly independent of the receiver.
  ComplexVector like1D(int size);

  /// Sets the matrix cell at coordinate `[row,column]` to the specified
  /// value.
  void put(int row, int column, Complex value) {
    if (column < 0 || column >= columns || row < 0 || row >= rows) {
      throw new RangeError("row:$row, column:$column");
    }
    set(row, column, value);
  }

  /// Sets the matrix cell at coordinate `[row,column]` to the specified
  /// value.
  void putParts(int row, int column, double re, double im) {
    if (column < 0 || column >= columns || row < 0 || row >= rows) {
      throw new RangeError("row:$row, column:$column");
    }
    setParts(row, column, re, im);
  }

  /// Sets the matrix cell at coordinate `[row,column]` to the specified
  /// value.
  void setParts(int row, int column, double re, double im);

  /// Sets the matrix cell at coordinate `[row,column]` to the specified
  /// value.
  void set(int row, int column, Complex value);

  /// Returns a string representation using default formatting ("%.4f").
  String toString() {
    return toStringFormat("%.4f");
  }

  /// Returns a string representation using using given [format].
  String toStringFormat(String format) {
    var f = new NumberFormat(format);
    var s = new StringBuffer("ComplexMatrix: $rows rows, $columns columns\n\n");
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        var elem = get(r, c);
        if (elem.imaginary == 0) {
          s.write(f.format(elem.real) + "\t");
          continue;
        }
        if (elem.real == 0) {
          s.write(f.format(elem.imaginary) + "i\t");
          continue;
        }
        if (elem.imaginary < 0) {
          s.write(f.format(elem.real) + " - " + f.format(-elem.imaginary) + "i\t");
          continue;
        }
        s.write(f.format(elem.real) + " + " + f.format(elem.imaginary) + "i\t");
      }
      s.write("\n");
    }
    return s.toString();
  }

  /// Constructs and returns a new *slice view* representing the rows of
  /// the given column.
  ComplexVector column(int column) {
    checkColumn(this, column);
    int viewSize = rows;
    int viewZero = index(0, column);
    int viewStride = rowStride;
    return _like1D(viewSize, viewZero, viewStride);
  }

  /// Constructs and returns a new *flip view* along the column axis. What
  /// used to be column `0` is now column `columns()-1`.
  ComplexMatrix columnFlip() {
    var v = _view();
    vColumnFlip(v);
    return v;
  }

  /// Constructs and returns a new *dice (transposition) view*; Swaps
  /// axes; example: 3 x 4 matrix --> 4 x 3 matrix.
  ComplexMatrix dice() {
    var v = _view();
    vDice(v);
    return v;
  }

  // Synonym for [dice].
  ComplexMatrix get T => dice();

  /// Constructs and returns a new *sub-range view* that is a
  /// `height x width` sub matrix starting at `[row,column]`.
  ComplexMatrix part(int row, int column, int height, int width) {
    var v = _view();
    vBox(v, row, column, height, width);
    return v;
  }

  /// Constructs and returns a new *slice view* representing the columns
  /// of the given row.
  ComplexVector row(int row) {
    checkRow(this, row);
    int viewSize = columns;
    int viewZero = index(row, 0);
    int viewStride = columnStride;
    return _like1D(viewSize, viewZero, viewStride);
  }

  /// Constructs and returns a new *flip view* along the row axis. What
  /// used to be row `0` is now row `rows()-1`.
  ComplexMatrix rowFlip() {
    var v = _view();
    vRowFlip(v);
    return v;
  }

  /// Constructs and returns a new *selection view* that is a matrix
  /// holding the indicated cells. Indexes can occur multiple times and
  /// can be in arbitrary order.
  ///
  /// To indicate "all" rows or "all columns", simply set the respective
  /// parameter
  ComplexMatrix select(Int32List rowIndexes, Int32List columnIndexes) {
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
  /// consisting of every i-th cell.
  ComplexMatrix strides(int rowStride, int columnStride) {
    var v = _view();
    vStrides(v, rowStride, columnStride);
    return v;
  }

  /// Linear algebraic matrix-vector multiplication;
  /// `z = alpha * A * y + beta*z`. Where `A == this`.
  ComplexVector mult(ComplexVector y,
      [ComplexVector z = null, Complex alpha = null,
      Complex beta = null, bool transposeA = false]) {
    if (alpha == null) {
      alpha = Complex.ONE;
    }
    if (beta == null) {
      beta = (z == null ? Complex.ONE : Complex.ZERO);
    }
    if (transposeA) {
      return conjugateTranspose().mult(y, z, alpha, beta, false);
    }
    ComplexVector zz;
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
    var s = Complex.ZERO;
    for (int r = 0; r < rows; r++) {
      s[0] = 0.0;
      s[1] = 0.0;
      for (int c = 0; c < columns; c++) {
        s += get(r, c) * y.get(c);
      }
      zz.set(r, (s * alpha) + (zz.get(r) * beta));
    }
    return zz;
  }

  /// Linear algebraic matrix-matrix multiplication;
  /// `C = alpha * A x B + beta*C`. Matrix shapes:
  /// `A(m x n), B(n x p), C(m x p)`.
  ComplexMatrix multiply(ComplexMatrix B,
      [ComplexMatrix C = null, Complex alpha = null,
      Complex beta = null, bool transposeA = false, bool transposeB = false]) {
    if (alpha == null) {
      alpha = Complex.ONE;
    }
    if (beta == null) {
      beta = (C == null ? Complex.ONE : Complex.ZERO);
    }
    if (transposeA) {
      return conjugateTranspose().multiply(
          B, C, alpha, beta, false, transposeB);
    }
    if (transposeB) {
      return this.multiply(
          B.conjugateTranspose(), C, alpha, beta, transposeA, false);
    }
    int m = rows;
    int n = columns;
    int p = B.columns;
    ComplexMatrix CC;
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
        var s = Complex.ZERO;
        for (int c = 0; c < n; c++) {
          s += get(b, c) * B.get(c, a);
        }
        CC.set(b, a, (s * alpha) + (CC.get(b, a) * beta));
      }
    }
    return CC;
  }

  /// Returns the sum of all cells.
  Complex sum() {
    if (size == 0) {
      return Complex.ZERO;
    }
    return aggregate(cfunc.plus, cfunc.identity);
  }

  /// Returns the content of this matrix if it is a wrapper; or `this`
  /// otherwise. Override this method in wrappers.
  ComplexMatrix _getContent() => this;

  /// Returns `true` if both matrices share at least one identical cell.
  bool _haveSharedCells(ComplexMatrix other) {
    if (other == null) {
      return false;
    }
    if (this == other) {
      return true;
    }
    return _getContent()._haveSharedCellsRaw(other._getContent());
  }

  // Always returns false.
  bool _haveSharedCellsRaw(ComplexMatrix other) => false;

  /// Construct and returns a new 1-d matrix of the corresponding dynamic
  /// type, sharing the same cells.
  ComplexVector _like1D(int size, int zero, int stride);

  /// Constructs and returns a new view equal to the receiver.
  ComplexMatrix _view() => clone() as ComplexMatrix;

  /// Construct and returns a new selection view.
  ComplexMatrix _viewSelectionLike(
      Int32List rowOffsets, Int32List columnOffsets);

  Object clone();

  ComplexMatrix operator *(ComplexMatrix y) {
    return copy()..assign(y, cfunc.mult);
  }

  ComplexMatrix operator /(ComplexMatrix y) {
    return copy()..assign(y, cfunc.div);
  }

  ComplexMatrix operator +(ComplexMatrix y) {
    return copy()..assign(y, cfunc.plus);
  }

  ComplexMatrix operator -(ComplexMatrix y) {
    return copy()..assign(y, cfunc.minus);
  }

  ComplexMatrix conj() {
    return copy()..apply(cfunc.conj);
  }
}
