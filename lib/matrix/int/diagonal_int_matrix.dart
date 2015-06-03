/*
Copyright (C) 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
is hereby granted without fee, provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear in supporting documentation.
CERN makes no representations about the suitability of this software for any purpose.
It is provided "as is" without expressed or implied warranty.
 */
part of cern.colt.matrix;

/**
 * Diagonal 2-d matrix holding <tt>int</tt> elements. First see the <a
 * href="package-summary.html">package summary</a> and javadoc <a
 * href="package-tree.html">tree view</a> to get the broad picture.
 * <p>
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
class DiagonalIntMatrix extends WrapperIntMatrix {

  /*
   * The non zero elements of the matrix.
   */
  Int32List _elements;

  /*
   * Length of the diagonal
   */
  int _dlength;

  /*
   * An m-by-n matrix A has m+n-1 diagonals. Since the DiagonalIntMatrix2D can have only one
   * diagonal, dindex is a value from interval [-m+1, n-1] that denotes which diagonal is stored.
   */
  int _dindex;

  /**
   * Constructs a matrix with a copy of the given values. <tt>values</tt> is
   * required to have the form <tt>values[row][column]</tt> and have exactly
   * the same number of columns in every row. Only the values on the main
   * diagonal, i.e. values[i][i] are used.
   * <p>
   * The values are copied. So subsequent changes in <tt>values</tt> are not
   * reflected in the matrix, and vice-versa.
   *
   * @param values
   *            The values to be filled into the new matrix.
   * @param dindex
   *            index of the diagonal.
   * @throws ArgumentError
   *             if
   *
   *             <tt>for any 1 &lt;= row &lt; values.length: values[row].length != values[row-1].length || index < -rows+1 || index > columns - 1</tt>
   *             .
   */
  factory DiagonalIntMatrix.fromList(List<Int32List> values, int dindex) {
    return new DiagonalIntMatrix(values.length, values.length == 0 ? 0 : values[0].length, dindex)
      ..setAll2D(values);
  }

  /**
   * Constructs a matrix with a given number of rows and columns. All entries
   * are initially <tt>0</tt>.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @param dindex
   *            index of the diagonal.
   * @throws ArgumentError
   *             if <tt>size<0 (int)size > Integer.MAX_VALUE</tt>.
   */
  DiagonalIntMatrix(int rows, int columns, int dindex) : super(null) {
    try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if ("matrix too large" != exc.message) {
        throw exc;
      }
    }
    if ((dindex < -rows + 1) || (dindex > columns - 1)) {
      throw new ArgumentError("index is out of bounds");
    } else {
      this._dindex = dindex;
    }
    if (dindex == 0) {
      _dlength = Math.min(rows, columns);
    } else if (dindex > 0) {
      if (rows >= columns) {
        _dlength = columns - dindex;
      } else {
        int diff = columns - rows;
        if (dindex <= diff) {
          _dlength = rows;
        } else {
          _dlength = rows - (dindex - diff);
        }
      }
    } else {
      if (rows >= columns) {
        int diff = rows - columns;
        if (-dindex <= diff) {
          _dlength = columns;
        } else {
          _dlength = columns + dindex + diff;
        }
      } else {
        _dlength = rows + dindex;
      }
    }
    _elements = new Int32List(_dlength);
  }

  void forEach(final ifunc.IntFunction function) {
    if (function is ifunc.IntMult) { // x[i] = mult*x[i]
      final int alpha = (function as ifunc.IntMult).multiplicator;
      if (alpha == 1) {
        return;
      }
      if (alpha == 0) {
        fill(0);
        return;
      }
      if (alpha != alpha) {
        fill(alpha); // the funny definition of isNaN(). This should better not happen.
        return;
      }
      for (int j = _dlength; --j >= 0; ) {
        _elements[j] *= alpha;
      }
    } else {
      for (int j = _dlength; --j >= 0; ) {
        _elements[j] = function(_elements[j]);
      }
    }
    return;
  }

  void fill(int value) {
    for (int i = _dlength; --i >= 0; ) {
      _elements[i] = value;
    }
    return;
  }

  void setAll(final Int32List values) {
    if (values.length != _dlength) {
      throw new ArgumentError("Must have same length: length=${values.length} dlength=$_dlength");
    }
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_dlength >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _dlength);
      List<Future> futures = new List<Future>(nthreads);
      int k = _dlength ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _dlength : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int r = firstRow; r < lastRow; r++) {
            _elements[r] = values[r];
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int r = _dlength; --r >= 0; ) {
        _elements[r] = values[r];
      }
    //}
    return;
  }

  void setAll2D(final List<Int32List> values) {
    if (values.length != _rows) {
      throw new ArgumentError("Must have same number of rows: rows=${values.length} rows()=$rows");
    }
    int r, c;
    if (_dindex >= 0) {
      r = 0;
      c = _dindex;
    } else {
      r = -_dindex;
      c = 0;
    }
    for (int i = 0; i < _dlength; i++) {
      if (values[i].length != _columns) {
        throw new ArgumentError("Must have same number of columns in every row: columns=${values[r].length} columns()=$columns");
      }
      _elements[i] = values[r++][c++];
    }
    return;
  }

  void copyFrom(AbstractIntMatrix source) {
    // overriden for performance only
    if (source == this) {
      return; // nothing to do
    }
    checkShape(source);

    if (source is DiagonalIntMatrix) {
      DiagonalIntMatrix other = source;
      if ((_dindex != other._dindex) || (_dlength != other._dlength)) {
        throw new ArgumentError("source is DiagonalIntMatrix2D with different diagonal stored.");
      }
      // quickest
      //System.arraycopy(other._elements, 0, this._elements, 0, this._elements.length);
      this._elements.setAll(0, other._elements);
      return;
    } else {
      super.copyFrom(source);
      return;
    }
  }

  void forEachWith(final AbstractIntMatrix y, final ifunc.IntIntFunction function) {
    checkShape(y);
    if (y is DiagonalIntMatrix) {
      DiagonalIntMatrix other = y;
      if ((_dindex != other._dindex) || (_dlength != other._dlength)) {
        throw new ArgumentError("y is DiagonalIntMatrix2D with different diagonal stored.");
      }
      if (function is ifunc.IntPlusMultSecond) { // x[i] = x[i] + alpha*y[i]
        final int alpha = (function as ifunc.IntPlusMultSecond).multiplicator;
        if (alpha == 0) {
          return; // nothing to do
        }
      }
      final Int32List otherElements = other._elements;
      /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
      if ((nthreads > 1) && (_dlength >= ConcurrencyUtils.getThreadsBeginN_2D())) {
        nthreads = Math.min(nthreads, _dlength);
        List<Future> futures = new List<Future>(nthreads);
        int k = _dlength ~/ nthreads;
        for (int j = 0; j < nthreads; j++) {
          final int firstRow = j * k;
          final int lastRow = (j == nthreads - 1) ? _dlength : firstRow + k;
          futures[j] = ConcurrencyUtils.submit(() {
            if (function is ifunc.IntPlusMultSecond) { // x[i] = x[i] + alpha*y[i]
              final int alpha = (function as ifunc.IntPlusMultSecond).multiplicator;
              if (alpha == 1) {
                for (int j = firstRow; j < lastRow; j++) {
                  _elements[j] += otherElements[j];
                }
              } else {
                for (int j = firstRow; j < lastRow; j++) {
                  _elements[j] = _elements[j] + alpha * otherElements[j];
                }
              }
            } else if (function == ifunc.mult) { // x[i] = x[i] * y[i]
              for (int j = firstRow; j < lastRow; j++) {
                _elements[j] = _elements[j] * otherElements[j];
              }
            } else if (function == ifunc.div) { // x[i] = x[i] /  y[i]
              for (int j = firstRow; j < lastRow; j++) {
                _elements[j] = _elements[j] ~/ otherElements[j];
              }
            } else {
              for (int j = firstRow; j < lastRow; j++) {
                _elements[j] = function(_elements[j], otherElements[j]);
              }
            }
          });
        }
        ConcurrencyUtils.waitForCompletion(futures);
      } else {*/
        if (function is ifunc.IntPlusMultSecond) { // x[i] = x[i] + alpha*y[i]
          final int alpha = (function as ifunc.IntPlusMultSecond).multiplicator;
          if (alpha == 1) {
            for (int j = _dlength; --j >= 0; ) {
              _elements[j] += otherElements[j];
            }
          } else {
            for (int j = _dlength; --j >= 0; ) {
              _elements[j] = _elements[j] + alpha * otherElements[j];
            }
          }
        } else if (function == ifunc.mult) { // x[i] = x[i] * y[i]
          for (int j = _dlength; --j >= 0; ) {
            _elements[j] = _elements[j] * otherElements[j];
          }
        } else if (function == ifunc.div) { // x[i] = x[i] /  y[i]
          for (int j = _dlength; --j >= 0; ) {
            _elements[j] = _elements[j] ~/ otherElements[j];
          }
        } else {
          for (int j = _dlength; --j >= 0; ) {
            _elements[j] = function(_elements[j], otherElements[j]);
          }
        }
      //}
      return;
    } else {
      super.forEachWith(y, function);
      return;
    }
  }

  int get cardinality {
    int cardinality = 0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_dlength >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _dlength);
      List<Future> futures = new List<Future>(nthreads);
      Int32List results = new Int32List(nthreads);
      int k = _dlength ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _dlength : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int cardinality = 0;
          for (int r = firstRow; r < lastRow; r++) {
            if (_elements[r] != 0) cardinality++;
          }
          return cardinality;
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          results[j] = futures[j].get() as int;
        }
        cardinality = results[0];
        for (int j = 1; j < nthreads; j++) {
          cardinality += results[j];
        }
      } on ExecutionException catch (ex) {
        ex.printStackTrace();
      } on InterruptedException catch (e) {
        e.printStackTrace();
      }
    } else {*/
      for (int r = 0; r < _dlength; r++) {
        if (_elements[r] != 0) cardinality++;
      }
    //}
    return cardinality;
  }

  Object get elements => _elements;

  bool all(int value) {
    for (int r = 0; r < _dlength; r++) {
      int x = _elements[r];
      int diff = value - x;
      if (diff != 0) {
        return false;
      }
    }
    return true;
  }

  bool equals(AbstractIntMatrix obj) {
    if (obj is DiagonalIntMatrix) {
      DiagonalIntMatrix other = obj;
      if (this == obj) {
        return true;
      }
      if (!(this != null && obj != null)) {
        return false;
      }
      final int rows = this.rows;
      final int columns = this.columns;
      if (columns != other.columns || rows != other.rows) {
        return false;
      }
      if ((_dindex != other._dindex) || (_dlength != other._dlength)) {
        return false;
      }
      Int32List otherElements = other._elements;
      for (int r = 0; r < _dlength; r++) {
        int x = _elements[r];
        int value = otherElements[r];
        int diff = value - x;
        if (diff != 0) {
          return false;
        }
      }
      return true;
    } else {
      return super.equals(obj);
    }
  }

  void forEachNonZero(final ifunc.IntIntIntFunction function) {
    for (int j = _dlength; --j >= 0; ) {
      int value = _elements[j];
      if (value != 0) {
        _elements[j] = function(j, j, value);
      }
    }
    return;
  }

  /**
   * Returns the length of the diagonal
   *
   * @return the length of the diagonal
   */
  int get diagonalLength {
    return _dlength;
  }

  /**
   * Returns the index of the diagonal
   *
   * @return the index of the diagonal
   */
  int get diagonalIndex {
    return _dindex;
  }

  IntMatrixLocation max() {
    int location = 0;
    int maxValue = 0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_dlength >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _dlength);
      List<Future> futures = new List<Future>(nthreads);
      List<Int32List> results = new List<Int32List>(nthreads);//[2];
      int k = _dlength ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _dlength : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int location = firstRow;
          int maxValue = _elements[location];
          int elem;
          for (int r = firstRow + 1; r < lastRow; r++) {
            elem = _elements[r];
            if (maxValue < elem) {
              maxValue = elem;
              location = r;
            }
          }
          return [maxValue, location, location];
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          results[j] = futures[j].get() as Int32List;
        }
        maxValue = results[0][0];
        location = results[0][1] as int;
        for (int j = 1; j < nthreads; j++) {
          if (maxValue < results[j][0]) {
            maxValue = results[j][0];
            location = results[j][1] as int;
          }
        }
      } on ExecutionException catch (ex) {
        ex.printStackTrace();
      } on InterruptedException catch (e) {
        e.printStackTrace();
      }
    } else {*/
      maxValue = _elements[0];
      int elem;
      for (int r = 1; r < _dlength; r++) {
        elem = _elements[r];
        if (maxValue < elem) {
          maxValue = elem;
          location = r;
        }
      }
    //}
    int rowLocation;
    int columnLocation;
    if (_dindex > 0) {
      rowLocation = location;
      columnLocation = location + _dindex;
    } else if (_dindex < 0) {
      rowLocation = location - _dindex;
      columnLocation = location;
    } else {
      rowLocation = location;
      columnLocation = location;
    }
    return new IntMatrixLocation(maxValue, rowLocation, columnLocation);
  }

  IntMatrixLocation min() {
    int location = 0;
    int minValue = 0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_dlength >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _dlength);
      List<Future> futures = new List<Future>(nthreads);
      List<Int32List> results = new List<Int32List>(nthreads);//[2];
      int k = _dlength ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _dlength : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int location = firstRow;
          int minValue = _elements[location];
          int elem;
          for (int r = firstRow + 1; r < lastRow; r++) {
            elem = _elements[r];
            if (minValue > elem) {
              minValue = elem;
              location = r;
            }
          }
          return [minValue, location, location];
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          results[j] = futures[j].get() as Int32List;
        }
        minValue = results[0][0];
        location = results[0][1] as int;
        for (int j = 1; j < nthreads; j++) {
          if (minValue > results[j][0]) {
            minValue = results[j][0];
            location = results[j][1] as int;
          }
        }
      } on ExecutionException catch (ex) {
        ex.printStackTrace();
      } on InterruptedException catch (e) {
        e.printStackTrace();
      }
    } else {*/
      minValue = _elements[0];
      int elem;
      for (int r = 1; r < _dlength; r++) {
        elem = _elements[r];
        if (minValue > elem) {
          minValue = elem;
          location = r;
        }
      }
    //}
    int rowLocation;
    int columnLocation;
    if (_dindex > 0) {
      rowLocation = location;
      columnLocation = location + _dindex;
    } else if (_dindex < 0) {
      rowLocation = location - _dindex;
      columnLocation = location;
    } else {
      rowLocation = location;
      columnLocation = location;
    }
    return new IntMatrixLocation(minValue, rowLocation, columnLocation);
  }

  int get(int row, int column) {
    if (_dindex >= 0) {
      if (column < _dindex) {
        return 0;
      } else {
        if ((row < _dlength) && (row + _dindex == column)) {
          return _elements[row];
        } else {
          return 0;
        }
      }
    } else {
      if (row < -_dindex) {
        return 0;
      } else {
        if ((column < _dlength) && (row + _dindex == column)) {
          return _elements[column];
        } else {
          return 0;
        }
      }
    }
  }

  AbstractIntMatrix like2D(int rows, int columns) {
    return new SparseIntMatrix(rows, columns);
  }

  AbstractIntVector like1D(int size) {
    return new SparseIntVector(size);
  }

  void set(int row, int column, int value) {
    if (_dindex >= 0) {
      if (column < _dindex) {
        //do nothing
      } else {
        if ((row < _dlength) && (row + _dindex == column)) {
          _elements[row] = value;
        } else {
          // do nothing
        }
      }
    } else {
      if (row < -_dindex) {
        //do nothing
      } else {
        if ((column < _dlength) && (row + _dindex == column)) {
          _elements[column] = value;
        } else {
          //do nothing;
        }
      }
    }
  }

  AbstractIntVector mult(AbstractIntVector y, [AbstractIntVector z = null, final int alpha = 1, int beta = null, final bool transposeA = false]) {
    if (beta == null) {
      beta = z == null ? 1 : 0;
    }
    int rowsA = _rows;
    int columnsA = _columns;
    if (transposeA) {
      rowsA = _columns;
      columnsA = _rows;
    }

    bool ignore = (z == null);
    if (z == null) z = new IntVector(rowsA);

    if (!(this._isNoView && y is IntVector && z is IntVector)) {
      return super.mult(y, z, alpha, beta, transposeA);
    }

    if (columnsA != y.length || rowsA > z.length) {
      throw new ArgumentError("Incompatible args: " + ((transposeA ? dice() : this).toStringShort()) + ", " + y.toStringShort() + ", " + z.toStringShort());
    }

    if ((!ignore) && ((beta) != 1)) {
      z.forEach(ifunc.multiply(beta));
    }

    IntVector zz = z as IntVector;
    final Int32List elementsZ = zz._elements;
    final int strideZ = zz.stride();
    final int zeroZ = z.index(0);

    IntVector yy = y as IntVector;
    final Int32List elementsY = yy._elements;
    final int strideY = yy.stride();
    final int zeroY = y.index(0);

    if (elementsY == null || elementsZ == null) {
      throw new Error();
    }
    if (!transposeA) {
      if (_dindex >= 0) {
        for (int i = _dlength; --i >= 0; ) {
          elementsZ[zeroZ + strideZ * i] += alpha * _elements[i] * elementsY[_dindex + zeroY + strideY * i];
        }
      } else {
        for (int i = _dlength; --i >= 0; ) {
          elementsZ[-_dindex + zeroZ + strideZ * i] += alpha * _elements[i] * elementsY[zeroY + strideY * i];
        }
      }
    } else {
      if (_dindex >= 0) {
        for (int i = _dlength; --i >= 0; ) {
          elementsZ[_dindex + zeroZ + strideZ * i] += alpha * _elements[i] * elementsY[zeroY + strideY * i];
        }
      } else {
        for (int i = _dlength; --i >= 0; ) {
          elementsZ[zeroZ + strideZ * i] += alpha * _elements[i] * elementsY[-_dindex + zeroY + strideY * i];
        }
      }

    }
    return z;
  }

  AbstractIntMatrix _getContent() {
    return this;
  }

  Object clone() {
    return new DiagonalIntMatrix(rows, columns, _dindex);
  }
}
