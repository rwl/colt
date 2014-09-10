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
 * Diagonal 2-d matrix holding <tt>complex</tt> elements. First see the <a
 * href="package-summary.html">package summary</a> and javadoc <a
 * href="package-tree.html">tree view</a> to get the broad picture.
 * <p>
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
class DiagonalComplexMatrix extends WrapperComplexMatrix {

  /*
   * The non zero elements of the matrix.
   */
  Float64List _elements;

  /*
   * Length of the diagonal
   */
  int _dlength;

  /*
   * An m-by-n matrix A has m+n-1 diagonals. Since the DiagonalComplexMatrix can have only one
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
  factory DiagonalComplexMatrix.fromList(List<Float64List> values, int dindex) {
    return new DiagonalComplexMatrix(values.length, values.length == 0 ? 0 : values[0].length, dindex)
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
   *             if <tt>size<0 (double)size > Integer.MAX_VALUE</tt>.
   */
  DiagonalComplexMatrix(int rows, int columns, this._dindex) : super(null) {
    try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if ("matrix too large" != exc.message) {
        throw exc;
      }
    }
    if ((_dindex < -rows + 1) || (_dindex > columns - 1)) {
      throw new ArgumentError("index is out of bounds");
    }
    if (_dindex == 0) {
      _dlength = Math.min(rows, columns);
    } else if (_dindex > 0) {
      if (rows >= columns) {
        _dlength = columns - _dindex;
      } else {
        int diff = columns - rows;
        if (_dindex <= diff) {
          _dlength = rows;
        } else {
          _dlength = rows - (_dindex - diff);
        }
      }
    } else {
      if (rows >= columns) {
        int diff = rows - columns;
        if (-_dindex <= diff) {
          _dlength = columns;
        } else {
          _dlength = columns + _dindex + diff;
        }
      } else {
        _dlength = rows + _dindex;
      }
    }
    _elements = new Float64List(2 * _dlength);
  }

  DiagonalComplexMatrix._internal(int rows, int columns, this._elements, this._dlength, this._dindex) : super(null) {
    try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if ("matrix too large" != exc.message) {
        throw exc;
      }
    }
  }

  void forEach(final cfunc.ComplexComplexFunction function) {
    if (function is cfunc.ComplexMult) { // x[i] = mult*x[i]
      final Float64List alpha = (function as cfunc.ComplexMult).multiplicator;
      if (alpha[0] == 1 && alpha[1] == 0) {
        return;
      }
      if (alpha[0] == 0 && alpha[1] == 0) {
        fill(alpha[0], alpha[1]);
        return;
      }
      if (alpha[0] != alpha[0] || alpha[1] != alpha[1]) {
        fill(alpha[0], alpha[1]); // the funny definition of isNaN(). This should better not happen.
        return;
      }
      Float64List elem = new Float64List(2);
      for (int j = 0; j < _dlength; j++) {
        elem[0] = _elements[2 * j];
        elem[1] = _elements[2 * j + 1];
        elem = Complex.multiply(elem, alpha);
        _elements[2 * j] = elem[0];
        _elements[2 * j + 1] = elem[1];
      }
    } else {
      Float64List elem = new Float64List(2);
      for (int j = 0; j < _dlength; j++) {
        elem[0] = _elements[2 * j];
        elem[1] = _elements[2 * j + 1];
        elem = function(elem);
        _elements[2 * j] = elem[0];
        _elements[2 * j + 1] = elem[1];
      }
    }
  }

  void fill(double re, double im) {
    for (int j = 0; j < _dlength; j++) {
      _elements[2 * j] = re;
      _elements[2 * j + 1] = im;
    }
  }

  void setAll(final Float64List values) {
    if (values.length != 2 * _dlength) {
      throw new ArgumentError("Must have same length: length=${values.length} 2*dlength=${2 * _dlength}");
    }
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_dlength >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _dlength);
      list<Future> futures = new List<Future>(nthreads);
      int k = _dlength / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _dlength : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int i = firstIdx; i < lastIdx; i++) {
            _elements[2 * i] = values[2 * i];
            _elements[2 * i + 1] = values[2 * i + 1];
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int i = 0; i < _dlength; i++) {
        _elements[2 * i] = values[2 * i];
        _elements[2 * i + 1] = values[2 * i + 1];
      }
    //}
  }

  void setAll2D(final List<Float64List> values) {
    if (values.length != _rows) {
      throw new ArgumentError("Must have same number of rows: rows=${values.length} rows()=${rows}");
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
      if (values[i].length != 2 * _columns) {
        throw new ArgumentError("Must have same number of columns in every row: columns=${values[r].length} 2 * columns()=${2 * columns}");
      }
      _elements[2 * i] = values[r][2 * c];
      _elements[2 * i + 1] = values[r][2 * c + 1];
      c++;
      r++;
    }
  }

  void copyFrom(AbstractComplexMatrix source) {
    // overriden for performance only
    if (source == this) {
      return ; // nothing to do
    }
    checkShape(source);

    if (source is DiagonalComplexMatrix) {
      DiagonalComplexMatrix other = source;
      if ((_dindex != other._dindex) || (_dlength != other._dlength)) {
        throw new ArgumentError("source is DiagonalComplexMatrix with different diagonal stored.");
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

  void forEachMatrix(final AbstractComplexMatrix y, final cfunc.ComplexComplexComplexFunction function) {
    checkShape(y);
    if (y is DiagonalComplexMatrix) {
      DiagonalComplexMatrix other = y;
      if ((_dindex != other._dindex) || (_dlength != other._dlength)) {
        throw new ArgumentError("y is DiagonalComplexMatrix with different diagonal stored.");
      }
      if (function is cfunc.ComplexPlusMultSecond) { // x[i] = x[i] + alpha*y[i]
        final Float64List alpha = (function as cfunc.ComplexPlusMultSecond).multiplicator;
        if (alpha[0] == 0 && alpha[1] == 0) {
          return; // nothing to do
        }
      }
      final Float64List otherElements = other._elements;
      /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
      if ((nthreads > 1) && (_dlength >= ConcurrencyUtils.getThreadsBeginN_2D())) {
        nthreads = Math.min(nthreads, _dlength);
        list<Future> futures = new List<Future>(nthreads);
        int k = _dlength / nthreads;
        for (int j = 0; j < nthreads; j++) {
          final int firstIdx = j * k;
          final int lastIdx = (j == nthreads - 1) ? _dlength : firstIdx + k;
          futures[j] = ConcurrencyUtils.submit(() {
            if (function is ComplexPlusMultSecond) { // x[i] = x[i] + alpha*y[i]
              final Float64List alpha = (function as ComplexPlusMultSecond).multiplicator;
              if (alpha[0] == 1 && alpha[1] == 0) {
                for (int j = firstIdx; j < lastIdx; j++) {
                  _elements[2 * j] += otherElements[2 * j];
                  _elements[2 * j + 1] += otherElements[2 * j + 1];
                }
              } else {
                Float64List elem = new Float64List(2);
                for (int j = firstIdx; j < lastIdx; j++) {
                  elem[0] = otherElements[2 * j];
                  elem[1] = otherElements[2 * j + 1];
                  elem = Complex.mult(alpha, elem);
                  _elements[2 * j] += elem[0];
                  _elements[2 * j + 1] += elem[1];
                }
              }
            } else if (function == ComplexFunctions.mult) { // x[i] = x[i] * y[i]
              Float64List elem = new Float64List(2);
              Float64List otherElem = new Float64List(2);
              for (int j = firstIdx; j < lastIdx; j++) {
                otherElem[0] = otherElements[2 * j];
                otherElem[1] = otherElements[2 * j + 1];
                elem[0] = _elements[2 * j];
                elem[1] = _elements[2 * j + 1];
                elem = Complex.mult(elem, otherElem);
                _elements[2 * j] = elem[0];
                _elements[2 * j + 1] = elem[1];
              }
            } else if (function == ComplexFunctions.div) { // x[i] = x[i] /  y[i]
              Float64List elem = new Float64List(2);
              Float64List otherElem = new Float64List(2);
              for (int j = firstIdx; j < lastIdx; j++) {
                otherElem[0] = otherElements[2 * j];
                otherElem[1] = otherElements[2 * j + 1];
                elem[0] = _elements[2 * j];
                elem[1] = _elements[2 * j + 1];
                elem = Complex.div(elem, otherElem);
                _elements[2 * j] = elem[0];
                _elements[2 * j + 1] = elem[1];
              }
            } else {
              Float64List elem = new Float64List(2);
              Float64List otherElem = new Float64List(2);
              for (int j = firstIdx; j < lastIdx; j++) {
                otherElem[0] = otherElements[2 * j];
                otherElem[1] = otherElements[2 * j + 1];
                elem[0] = _elements[2 * j];
                elem[1] = _elements[2 * j + 1];
                elem = function.apply(elem, otherElem);
                _elements[2 * j] = elem[0];
                _elements[2 * j + 1] = elem[1];
              }
            }
          });
        }
        ConcurrencyUtils.waitForCompletion(futures);
      } else {*/
        if (function is cfunc.ComplexPlusMultSecond) { // x[i] = x[i] + alpha*y[i]
          final Float64List alpha = (function as cfunc.ComplexPlusMultSecond).multiplicator;
          if (alpha[0] == 1 && alpha[1] == 0) {
            for (int j = 0; j < _dlength; j++) {
              _elements[2 * j] += otherElements[2 * j];
              _elements[2 * j + 1] += otherElements[2 * j + 1];
            }
          } else {
            Float64List elem = new Float64List(2);
            for (int j = 0; j < _dlength; j++) {
              elem[0] = otherElements[2 * j];
              elem[1] = otherElements[2 * j + 1];
              elem = Complex.multiply(alpha, elem);
              _elements[2 * j] += elem[0];
              _elements[2 * j + 1] += elem[1];
            }
          }
        } else if (function == cfunc.mult) { // x[i] = x[i] * y[i]
          Float64List elem = new Float64List(2);
          Float64List otherElem = new Float64List(2);
          for (int j = 0; j < _dlength; j++) {
            otherElem[0] = otherElements[2 * j];
            otherElem[1] = otherElements[2 * j + 1];
            elem[0] = _elements[2 * j];
            elem[1] = _elements[2 * j + 1];
            elem = Complex.multiply(elem, otherElem);
            _elements[2 * j] = elem[0];
            _elements[2 * j + 1] = elem[1];
          }
        } else if (function == cfunc.div) { // x[i] = x[i] /  y[i]
          Float64List elem = new Float64List(2);
          Float64List otherElem = new Float64List(2);
          for (int j = 0; j < _dlength; j++) {
            otherElem[0] = otherElements[2 * j];
            otherElem[1] = otherElements[2 * j + 1];
            elem[0] = _elements[2 * j];
            elem[1] = _elements[2 * j + 1];
            elem = Complex.div_(elem, otherElem);
            _elements[2 * j] = elem[0];
            _elements[2 * j + 1] = elem[1];
          }
        } else {
          Float64List elem = new Float64List(2);
          Float64List otherElem = new Float64List(2);
          for (int j = 0; j < _dlength; j++) {
            otherElem[0] = otherElements[2 * j];
            otherElem[1] = otherElements[2 * j + 1];
            elem[0] = _elements[2 * j];
            elem[1] = _elements[2 * j + 1];
            elem = function(elem, otherElem);
            _elements[2 * j] = elem[0];
            _elements[2 * j + 1] = elem[1];
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
      list<Future> futures = new List<Future>(nthreads);
      List<int> results = new List<int>(nthreads);
      int k = _dlength / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _dlength : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int cardinality = 0;
          for (int i = firstIdx; i < lastIdx; i++) {
            if (_elements[2 * i] != 0 || _elements[2 * i + 1] != 0) cardinality++;
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
      for (int i = 0; i < _dlength; i++) {
        if (_elements[2 * i] != 0 || _elements[2 * i + 1] != 0) {
          cardinality++;
        }
      }
    //}
    return cardinality;
  }

  Float64List elements() {
    return _elements;
  }

  /*bool equalsValue(Float64List value) {
    double epsilon = ComplexProperty.DEFAULT.tolerance();
    Float64List x = new Float64List(2);
    Float64List diff = new Float64List(2);
    for (int i = 0; i < _dlength; i++) {
      x[0] = _elements[2 * i];
      x[1] = _elements[2 * i + 1];
      diff[0] = (value[0] - x[0]).abs();
      diff[1] = (value[1] - x[1]).abs();
      if (((diff[0] != diff[0]) || (diff[1] != diff[1])) && ((((value[0] != value[0]) || (value[1] != value[1])) && ((x[0] != x[0]) || (x[1] != x[1])))) || (Complex.isEqual(value, x, epsilon))) {
        diff[0] = 0.0;
        diff[1] = 0.0;
      }
      if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
        return false;
      }
    }
    return true;
  }*/

  bool operator ==(var obj) {
    if (obj is Float64List) {
      final value = obj;
      double epsilon = EPSILON;
      Float64List x = new Float64List(2);
      Float64List diff = new Float64List(2);
      for (int i = 0; i < _dlength; i++) {
        x[0] = _elements[2 * i];
        x[1] = _elements[2 * i + 1];
        diff[0] = (value[0] - x[0]).abs();
        diff[1] = (value[1] - x[1]).abs();
        if (((diff[0] != diff[0]) || (diff[1] != diff[1])) && ((((value[0] != value[0]) || (value[1] != value[1])) && ((x[0] != x[0]) || (x[1] != x[1])))) || (Complex.isEqual(value, x, epsilon))) {
          diff[0] = 0.0;
          diff[1] = 0.0;
        }
        if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
          return false;
        }
      }
      return true;
    }
    if (obj is DiagonalComplexMatrix) {
      DiagonalComplexMatrix other = obj;
      double epsilon = EPSILON;
      if (identical(this, obj)) {
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
      Float64List otherElements = other._elements;
      Float64List x = new Float64List(2);
      Float64List value = new Float64List(2);
      Float64List diff = new Float64List(2);
      for (int i = 0; i < _dlength; i++) {
        x[0] = _elements[2 * i];
        x[1] = _elements[2 * i + 1];
        value[0] = otherElements[2 * i];
        value[1] = otherElements[2 * i + 1];
        diff[0] = (value[0] - x[0]).abs();
        diff[1] = (value[1] - x[1]).abs();
        if (((diff[0] != diff[0]) || (diff[1] != diff[1])) && ((((value[0] != value[0]) || (value[1] != value[1])) && ((x[0] != x[0]) || (x[1] != x[1])))) || (Complex.isEqual(value, x, epsilon))) {
          diff[0] = 0.0;
          diff[1] = 0.0;
        }
        if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
          return false;
        }
      }
      return true;
    } else {
      return super ==(obj);
    }
  }

  void forEachNonZero(final cfunc.IntIntComplexFunction function) {
    Float64List value = new Float64List(2);
    for (int i = 0; i < _dlength; i++) {
      value[0] = _elements[2 * i];
      value[1] = _elements[2 * i + 1];
      if (value[0] != 0 || value[1] != 0) {
        value = function(i, i, value);
        _elements[2 * i] = value[0];
        _elements[2 * i + 1] = value[1];
      }
    }
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

  Float64List get(int row, int column) {
    if (_dindex >= 0) {
      if (column < _dindex) {
        return new Float64List(2);
      } else {
        if ((row < _dlength) && (row + _dindex == column)) {
          return new Float64List.fromList([_elements[2 * row], _elements[2 * row + 1]]);
        } else {
          return new Float64List(2);
        }
      }
    } else {
      if (row < -_dindex) {
        return new Float64List(2);
      } else {
        if ((column < _dlength) && (row + _dindex == column)) {
          return new Float64List.fromList([_elements[2 * column], _elements[2 * column + 1]]);
        } else {
          return new Float64List(2);
        }
      }
    }
  }

  AbstractComplexMatrix like2D(int rows, int columns) {
    return new SparseComplexMatrix(rows, columns);
  }

  AbstractComplexVector like1D(int size) {
    return new SparseComplexVector(size);
  }

  void set(int row, int column, Float64List value) {
    if (_dindex >= 0) {
      if (column < _dindex) {
        //do nothing
      } else {
        if ((row < _dlength) && (row + _dindex == column)) {
          _elements[2 * row] = value[0];
          _elements[2 * row + 1] = value[1];
        } else {
          // do nothing
        }
      }
    } else {
      if (row < -_dindex) {
        //do nothing
      } else {
        if ((column < _dlength) && (row + _dindex == column)) {
          _elements[2 * column] = value[0];
          _elements[2 * column + 1] = value[1];
        } else {
          //do nothing;
        }
      }
    }
  }

  void setParts(int row, int column, double re, double im) {
    if (_dindex >= 0) {
      if (column < _dindex) {
        //do nothing
      } else {
        if ((row < _dlength) && (row + _dindex == column)) {
          _elements[2 * row] = re;
          _elements[2 * row + 1] = im;
        } else {
          // do nothing
        }
      }
    } else {
      if (row < -_dindex) {
        //do nothing
      } else {
        if ((column < _dlength) && (row + _dindex == column)) {
          _elements[2 * column] = re;
          _elements[2 * column + 1] = im;
        } else {
          //do nothing;
        }
      }
    }
  }

  AbstractComplexVector mult(AbstractComplexVector y, [AbstractComplexVector z = null, Float64List alpha = null, Float64List beta = null, bool transposeA = false]) {
    if (alpha == null) {
      alpha = new Float64List.fromList([1.0, 0.0]);
    }
    if (beta == null) {
      beta = (z == null ? new Float64List.fromList([1.0, 0.0]) : new Float64List.fromList([0.0, 0.0]));
    }
    int rowsA = _rows;
    int columnsA = _columns;
    if (transposeA) {
      rowsA = _columns;
      columnsA = _rows;
    }

    bool ignore = (z == null);
    if (z == null) z = new ComplexVector(rowsA);

    if (!(this._isNoView && y is ComplexVector && z is ComplexVector)) {
      return super.mult(y, z, alpha, beta, transposeA);
    }

    if (columnsA != y.length || rowsA > z.length) {
      throw new ArgumentError("Incompatible args: " + ((transposeA ? dice() : this).toStringShort()) + ", " + y.toStringShort() + ", " + z.toStringShort());
    }

    if ((!ignore) && !((beta[0] == 1) && (beta[1] == 0))) {
      z.forEach(cfunc.multiply(beta));
    }

    ComplexVector zz = z as ComplexVector;
    final Float64List elementsZ = zz._elements;
    final int strideZ = zz.stride();
    final int zeroZ = z.index(0);

    ComplexVector yy = y as ComplexVector;
    final Float64List elementsY = yy._elements;
    final int strideY = yy.stride();
    final int zeroY = y.index(0);

    if (elementsY == null || elementsZ == null) {
      throw new Error();
    }
    Float64List elemA = new Float64List(2);
    Float64List elemY = new Float64List(2);
    if (!transposeA) {
      if (_dindex >= 0) {
        for (int i = 0; i < _dlength; i++) {
          elemA[0] = _elements[2 * i];
          elemA[1] = _elements[2 * i + 1];
          elemY[0] = elementsY[2 * _dindex + zeroY + strideY * i];
          elemY[1] = elementsY[2 * _dindex + zeroY + strideY * i + 1];
          elemA = Complex.multiply(elemA, elemY);
          elemA = Complex.multiply(alpha, elemA);
          elementsZ[zeroZ + strideZ * i] += elemA[0];
          elementsZ[zeroZ + strideZ * i + 1] += elemA[1];
        }
      } else {
        for (int i = 0; i < _dlength; i++) {
          elemA[0] = _elements[2 * i];
          elemA[1] = _elements[2 * i + 1];
          elemY[0] = elementsY[zeroY + strideY * i];
          elemY[1] = elementsY[zeroY + strideY * i + 1];
          elemA = Complex.multiply(elemA, elemY);
          elemA = Complex.multiply(alpha, elemA);
          elementsZ[-2 * _dindex + zeroZ + strideZ * i] += elemA[0];
          elementsZ[-2 * _dindex + zeroZ + strideZ * i + 1] += elemA[1];
        }
      }
    } else {
      if (_dindex >= 0) {
        for (int i = 0; i < _dlength; i++) {
          elemA[0] = _elements[2 * i];
          elemA[1] = -_elements[2 * i + 1];
          elemY[0] = elementsY[zeroY + strideY * i];
          elemY[1] = elementsY[zeroY + strideY * i + 1];
          elemA = Complex.multiply(elemA, elemY);
          elemA = Complex.multiply(alpha, elemA);
          elementsZ[2 * _dindex + zeroZ + strideZ * i] += elemA[0];
          elementsZ[2 * _dindex + zeroZ + strideZ * i + 1] += elemA[1];
        }
      } else {
        for (int i = 0; i < _dlength; i++) {
          elemA[0] = _elements[2 * i];
          elemA[1] = -_elements[2 * i + 1];
          elemY[0] = elementsY[-2 * _dindex + zeroY + strideY * i];
          elemY[1] = elementsY[-2 * _dindex + zeroY + strideY * i + 1];
          elemA = Complex.multiply(elemA, elemY);
          elemA = Complex.multiply(alpha, elemA);
          elementsZ[zeroZ + strideZ * i] += elemA[0];
          elementsZ[zeroZ + strideZ * i + 1] += elemA[1];
        }
      }

    }
    return z;
  }

  AbstractComplexMatrix _getContent() {
    return this;
  }

  Object clone() {
    return new DiagonalComplexMatrix._internal(_rows, _columns, _elements, _dlength, _dindex);
  }
}
