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
 * Diagonal 2-d matrix holding <tt>double</tt> elements. First see the <a
 * href="package-summary.html">package summary</a> and javadoc <a
 * href="package-tree.html">tree view</a> to get the broad picture.
 * <p>
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
class DiagonalDoubleMatrix2D extends WrapperDoubleMatrix2D {

  /*
   * The non zero elements of the matrix.
   */
  Float64List _elements;

  /*
   * Length of the diagonal
   */
  int _dlength;

  /*
   * An m-by-n matrix A has m+n-1 diagonals. Since the DiagonalDoubleMatrix2D can have only one
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
  factory DiagonalDoubleMatrix2D.fromList(List<Float64List> values, int dindex) {
    return new DiagonalDoubleMatrix2D.index(values.length, values.length == 0 ? 0 : values[0].length, dindex)
      ..assignValues2D(values);
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
   *             if <tt>size<0 (double)size > Integer.MAX_VALUE</tt>.
   */
  factory DiagonalDoubleMatrix2D.index(int rows, int columns, int dindex) {
    if ((dindex < -rows + 1) || (dindex > columns - 1)) {
      throw new ArgumentError("index is out of bounds");
    }
    int dlength;
    if (dindex == 0) {
      dlength = Math.min(rows, columns);
    } else if (dindex > 0) {
      if (rows >= columns) {
        dlength = columns - dindex;
      } else {
        int diff = columns - rows;
        if (dindex <= diff) {
          dlength = rows;
        } else {
          dlength = rows - (dindex - diff);
        }
      }
    } else {
      if (rows >= columns) {
        int diff = rows - columns;
        if (-dindex <= diff) {
          dlength = columns;
        } else {
          dlength = columns + dindex + diff;
        }
      } else {
        dlength = rows + dindex;
      }
    }
    final elements = new Float64List(dlength);
    return new DiagonalDoubleMatrix2D._internal(rows, columns, elements, dlength, dindex); 
  }
  
  DiagonalDoubleMatrix2D._internal(int rows, int columns, this._elements, this._dlength, this._dindex) : super(null) {
    try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if ("matrix too large" != exc.message) {
        throw exc;
      }
    }    
  }

  DoubleMatrix2D assign(final func.DoubleFunction function) {
    if (function is DoubleMult) { // x[i] = mult*x[i]
      final double alpha = (function as DoubleMult).multiplicator;
      if (alpha == 1) {
        return this;
      }
      if (alpha == 0) {
        return assignValue(0.0);
      }
      if (alpha != alpha) {
        return assignValue(alpha); // the funny definition of isNaN(). This should better not happen.
      }
      for (int j = _dlength; --j >= 0; ) {
        _elements[j] *= alpha;
      }
    } else {
      for (int j = _dlength; --j >= 0; ) {
        _elements[j] = function(_elements[j]);
      }
    }
    return this;
  }

  DoubleMatrix2D assignValue(double value) {
    for (int i = _dlength; --i >= 0; ) _elements[i] = value;
    return this;
  }

  DoubleMatrix2D assignValues(final Float64List values) {
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
    return this;
  }

  DoubleMatrix2D assignValues2D(final List<Float64List> values) {
    if (values.length != _rows) {
      throw new ArgumentError("Must have same number of rows: rows=${values.length} rows()=${rows()}");
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
        throw new ArgumentError("Must have same number of columns in every row: columns=${values[r].length} columns()=${columns()}");
      }
      _elements[i] = values[r++][c++];
    }
    return this;
  }

  DoubleMatrix2D assignMatrix(DoubleMatrix2D source) {
    // overriden for performance only
    if (source == this) return this; // nothing to do
    checkShape(source);

    if (source is DiagonalDoubleMatrix2D) {
      DiagonalDoubleMatrix2D other = source;
      if ((_dindex != other._dindex) || (_dlength != other._dlength)) {
        throw new ArgumentError("source is DiagonalDoubleMatrix2D with different diagonal stored.");
      }
      // quickest
      this._elements.setAll(0, other._elements);
      //System.arraycopy(other._elements, 0, this._elements, 0, this._elements.length);
      return this;
    } else {
      return super.assignMatrix(source);
    }
  }

  DoubleMatrix2D assignFunc(final DoubleMatrix2D y, final func.DoubleDoubleFunction function) {
    checkShape(y);
    if (y is DiagonalDoubleMatrix2D) {
      DiagonalDoubleMatrix2D other = y;
      if ((_dindex != other._dindex) || (_dlength != other._dlength)) {
        throw new ArgumentError("y is DiagonalDoubleMatrix2D with different diagonal stored.");
      }
      if (function is DoublePlusMultSecond) { // x[i] = x[i] + alpha*y[i]
        final double alpha = (function as DoublePlusMultSecond).multiplicator;
        if (alpha == 0) {
          return this; // nothing to do
        }
      }
      final Float64List otherElements = other._elements;
      /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
      if ((nthreads > 1) && (_dlength >= ConcurrencyUtils.getThreadsBeginN_2D())) {
        nthreads = Math.min(nthreads, _dlength);
        List<Future> futures = new List<Future>(nthreads);
        int k = _dlength ~/ nthreads;
        for (int j = 0; j < nthreads; j++) {
          final int firstRow = j * k;
          final int lastRow = (j == nthreads - 1) ? _dlength : firstRow + k;
          futures[j] = ConcurrencyUtils.submit(() {
            if (function is DoublePlusMultSecond) { // x[i] = x[i] + alpha*y[i]
              final double alpha = (function as DoublePlusMultSecond).multiplicator;
              if (alpha == 1) {
                for (int j = firstRow; j < lastRow; j++) {
                  _elements[j] += otherElements[j];
                }
              } else {
                for (int j = firstRow; j < lastRow; j++) {
                  _elements[j] = _elements[j] + alpha * otherElements[j];
                }
              }
            } else if (function == func.mult) { // x[i] = x[i] * y[i]
              for (int j = firstRow; j < lastRow; j++) {
                _elements[j] = _elements[j] * otherElements[j];
              }
            } else if (function == func.div) { // x[i] = x[i] /  y[i]
              for (int j = firstRow; j < lastRow; j++) {
                _elements[j] = _elements[j] / otherElements[j];
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
        if (function is DoublePlusMultSecond) { // x[i] = x[i] + alpha*y[i]
          final double alpha = (function as DoublePlusMultSecond).multiplicator;
          if (alpha == 1) {
            for (int j = _dlength; --j >= 0; ) {
              _elements[j] += otherElements[j];
            }
          } else {
            for (int j = _dlength; --j >= 0; ) {
              _elements[j] = _elements[j] + alpha * otherElements[j];
            }
          }
        } else if (function == func.mult) { // x[i] = x[i] * y[i]
          for (int j = _dlength; --j >= 0; ) {
            _elements[j] = _elements[j] * otherElements[j];
          }
        } else if (function == func.div) { // x[i] = x[i] /  y[i]
          for (int j = _dlength; --j >= 0; ) {
            _elements[j] = _elements[j] / otherElements[j];
          }
        } else {
          for (int j = _dlength; --j >= 0; ) {
            _elements[j] = function(_elements[j], otherElements[j]);
          }
        }
      //}
      return this;
    } else {
      return super.assignFunc(y, function);
    }
  }

  int cardinality() {
    int cardinality = 0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_dlength >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _dlength);
      List<Future> futures = new List<Future>(nthreads);
      List<int> results = new List<int>(nthreads);
      int k = _dlength / nthreads;
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
          results[j] = futures[j].get();
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

  Float64List elements() {
    return _elements;
  }

  bool equalsValues(double value) {
    double epsilon = DoubleProperty.DEFAULT.tolerance();
    for (int r = 0; r < _dlength; r++) {
      double x = _elements[r];
      double diff = (value - x).abs();
      if ((diff != diff) && ((value != value && x != x) || value == x)) {
        diff = 0.0;
      }
      if (!(diff <= epsilon)) {
        return false;
      }
    }
    return true;
  }

  bool equals(Object obj) {
    if (obj is DiagonalDoubleMatrix2D) {
      DiagonalDoubleMatrix2D other = obj;
      double epsilon = DoubleProperty.DEFAULT.tolerance();
      if (this == obj) return true;
      if (!(this != null && obj != null)) return false;
      final int rows = this.rows();
      final int columns = this.columns();
      if (columns != other.columns() || rows != other.rows()) return false;
      if ((_dindex != other._dindex) || (_dlength != other._dlength)) {
        return false;
      }
      Float64List otherElements = other._elements;
      for (int r = 0; r < _dlength; r++) {
        double x = _elements[r];
        double value = otherElements[r];
        double diff = (value - x).abs();
        if ((diff != diff) && ((value != value && x != x) || value == x)) {
          diff = 0.0;
        }
        if (!(diff <= epsilon)) {
          return false;
        }
      }
      return true;
    } else {
      return super.equals(obj);
    }
  }

  DoubleMatrix2D forEachNonZero(final func.IntIntDoubleFunction function) {
    for (int j = _dlength; --j >= 0; ) {
      double value = _elements[j];
      if (value != 0) {
        _elements[j] = function(j, j, value);
      }
    }
    return this;
  }

  /**
   * Returns the length of the diagonal
   *
   * @return the length of the diagonal
   */
  int diagonalLength() {
    return _dlength;
  }

  /**
   * Returns the index of the diagonal
   *
   * @return the index of the diagonal
   */
  int diagonalIndex() {
    return _dindex;
  }

  Float64List getMaxLocation() {
    int location = 0;
    double maxValue = 0.0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_dlength >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _dlength);
      List<Future> futures = new List<Future>(nthreads);
      List<Float64List> results = new List<double>(nthreads);//[2];
      int k = _dlength ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _dlength : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int location = firstRow;
          double maxValue = _elements[location];
          double elem;
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
          results[j] = futures[j].get() as Float64List;
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
      double elem;
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
    return [maxValue, rowLocation, columnLocation];
  }

  Float64List getMinLocation() {
    int location = 0;
    double minValue = 0.0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_dlength >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _dlength);
      List<Future> futures = new List<Future>(nthreads);
      List<Float64List> results = new List<double>(nthreads);//[2];
      int k = _dlength ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _dlength : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int location = firstRow;
          double minValue = _elements[location];
          double elem;
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
          results[j] = futures[j].get() as Float64List;
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
      double elem;
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
    return [minValue, rowLocation, columnLocation];
  }

  double getQuick(int row, int column) {
    if (_dindex >= 0) {
      if (column < _dindex) {
        return 0.0;
      } else {
        if ((row < _dlength) && (row + _dindex == column)) {
          return _elements[row];
        } else {
          return 0.0;
        }
      }
    } else {
      if (row < -_dindex) {
        return 0.0;
      } else {
        if ((column < _dlength) && (row + _dindex == column)) {
          return _elements[column];
        } else {
          return 0.0;
        }
      }
    }
  }

  DoubleMatrix2D like2D(int rows, int columns) {
    return new SparseDoubleMatrix2D(rows, columns);
  }

  DoubleMatrix1D like1D(int size) {
    return new SparseDoubleMatrix1D(size);
  }

  void setQuick(int row, int column, double value) {
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

  DoubleMatrix1D zMult(DoubleMatrix1D y, DoubleMatrix1D z, [double alpha=1.0, double beta=0.0, final bool transposeA=false]) {
    int rowsA = _rows;
    int columnsA = _columns;
    if (transposeA) {
      rowsA = _columns;
      columnsA = _rows;
    }

    bool ignore = (z == null);
    if (z == null) z = new DenseDoubleMatrix1D(rowsA);

    if (!(this._isNoView && y is DenseDoubleMatrix1D && z is DenseDoubleMatrix1D)) {
      return super.zMult(y, z, alpha, beta, transposeA);
    }

    if (columnsA != y.size() || rowsA > z.size()) {
      throw new ArgumentError("Incompatible args: " + ((transposeA ? viewDice() : this).toStringShort()) + ", " + y.toStringShort() + ", " + z.toStringShort());
    }

    if ((!ignore) && ((beta) != 1)) {
      z.assign(func.multiply(beta));
    }

    DenseDoubleMatrix1D zz = z as DenseDoubleMatrix1D;
    final Float64List elementsZ = zz._elements;
    final int strideZ = zz.stride();
    final int zeroZ = z.index(0);

    DenseDoubleMatrix1D yy = y as DenseDoubleMatrix1D;
    final Float64List elementsY = yy._elements;
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

  DoubleMatrix2D _getContent() {
    return this;
  }
  
  Object clone() {
    return new DiagonalDoubleMatrix2D._internal(_rows, _columns, _elements, _dlength, _dindex);
  }
}
