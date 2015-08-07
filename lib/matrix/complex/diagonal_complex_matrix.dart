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

/// Diagonal 2-d matrix holding `complex` elements.
class DiagonalComplexMatrix extends WrapperComplexMatrix {
  /// The non zero elements of the matrix.
  Float64List _elements;

  /// Length of the diagonal
  int _dlength;

  /// An m-by-n matrix A has m+n-1 diagonals. Since the DiagonalComplexMatrix
  /// can have only one diagonal, dindex is a value from interval [-m+1, n-1]
  /// that denotes which diagonal is stored.
  int _dindex;

  /// Constructs a matrix with a given number of rows and columns. All entries
  /// are initially `0`.
  DiagonalComplexMatrix(int rows, int columns, [this._dindex = 0]) : super._(rows, columns) {
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

  DiagonalComplexMatrix._internal(int rows, int columns, this._elements, this._dlength, this._dindex) : super._(rows, columns);

  static DiagonalComplexMatrix create(int rows, int columns) {
    return new DiagonalComplexMatrix(rows, columns);
  }

  void apply(final cfunc.ComplexComplexFunction fn) {
    if (fn is cfunc.ComplexMult) { // x[i] = mult*x[i]
      final Float64List alpha = fn.multiplicator;
      if (alpha[0] == 1 && alpha[1] == 0) {
        return;
      }
      if (alpha[0] == 0 && alpha[1] == 0) {
        fill(alpha[0], alpha[1]);
        return;
      }
      if (alpha[0] != alpha[0] || alpha[1] != alpha[1]) {
        fill(alpha[0], alpha[1]); // isNaN. This should not happen.
        return;
      }
      Float64List elem = new Float64List(2);
      for (int j = 0; j < _dlength; j++) {
        elem[0] = _elements[2 * j];
        elem[1] = _elements[2 * j + 1];
        elem = cmath.multiply(elem, alpha);
        _elements[2 * j] = elem[0];
        _elements[2 * j + 1] = elem[1];
      }
    } else {
      Float64List elem = new Float64List(2);
      for (int j = 0; j < _dlength; j++) {
        elem[0] = _elements[2 * j];
        elem[1] = _elements[2 * j + 1];
        elem = fn(elem);
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
    for (int i = 0; i < _dlength; i++) {
      _elements[2 * i] = values[2 * i];
      _elements[2 * i + 1] = values[2 * i + 1];
    }
  }

  void copyFrom(AbstractComplexMatrix source) {
    // overriden for performance only
    if (source == this) {
      return ; // nothing to do
    }
    checkShape(this, source);

    if (source is DiagonalComplexMatrix) {
      DiagonalComplexMatrix other = source;
      if ((_dindex != other._dindex) || (_dlength != other._dlength)) {
        throw new ArgumentError("source is DiagonalComplexMatrix with different diagonal stored.");
      }
      // quickest
      _elements.setAll(0, other._elements);
      return;
    } else {
      super.copyFrom(source);
      return;
    }
  }

  void assign(final AbstractComplexMatrix y, final cfunc.ComplexComplexComplexFunction fn) {
    checkShape(this, y);
    if (y is DiagonalComplexMatrix) {
      DiagonalComplexMatrix other = y;
      if ((_dindex != other._dindex) || (_dlength != other._dlength)) {
        throw new ArgumentError("y is DiagonalComplexMatrix with different diagonal stored.");
      }
      if (fn is cfunc.ComplexPlusMultSecond) { // x[i] = x[i] + alpha*y[i]
        var alpha = fn.multiplicator;
        if (alpha[0] == 0 && alpha[1] == 0) {
          return; // nothing to do
        }
      }
      Float64List otherElements = other._elements;
        if (fn is cfunc.ComplexPlusMultSecond) { // x[i] = x[i] + alpha*y[i]
          Float64List alpha = fn.multiplicator;
          if (alpha[0] == 1 && alpha[1] == 0) {
            for (int j = 0; j < _dlength; j++) {
              _elements[2 * j] += otherElements[2 * j];
              _elements[2 * j + 1] += otherElements[2 * j + 1];
            }
          } else {
            var elem = new Float64List(2);
            for (int j = 0; j < _dlength; j++) {
              elem[0] = otherElements[2 * j];
              elem[1] = otherElements[2 * j + 1];
              elem = cmath.multiply(alpha, elem);
              _elements[2 * j] += elem[0];
              _elements[2 * j + 1] += elem[1];
            }
          }
        } else if (fn == cfunc.mult) { // x[i] = x[i] * y[i]
          var elem = new Float64List(2);
          var otherElem = new Float64List(2);
          for (int j = 0; j < _dlength; j++) {
            otherElem[0] = otherElements[2 * j];
            otherElem[1] = otherElements[2 * j + 1];
            elem[0] = _elements[2 * j];
            elem[1] = _elements[2 * j + 1];
            elem = cmath.multiply(elem, otherElem);
            _elements[2 * j] = elem[0];
            _elements[2 * j + 1] = elem[1];
          }
        } else if (fn == cfunc.div) { // x[i] = x[i] /  y[i]
          var elem = new Float64List(2);
          var otherElem = new Float64List(2);
          for (int j = 0; j < _dlength; j++) {
            otherElem[0] = otherElements[2 * j];
            otherElem[1] = otherElements[2 * j + 1];
            elem[0] = _elements[2 * j];
            elem[1] = _elements[2 * j + 1];
            elem = cmath.div_(elem, otherElem);
            _elements[2 * j] = elem[0];
            _elements[2 * j + 1] = elem[1];
          }
        } else {
          var elem = new Float64List(2);
          var otherElem = new Float64List(2);
          for (int j = 0; j < _dlength; j++) {
            otherElem[0] = otherElements[2 * j];
            otherElem[1] = otherElements[2 * j + 1];
            elem[0] = _elements[2 * j];
            elem[1] = _elements[2 * j + 1];
            elem = fn(elem, otherElem);
            _elements[2 * j] = elem[0];
            _elements[2 * j + 1] = elem[1];
          }
        }
      return;
    } else {
      super.assign(y, fn);
      return;
    }
  }

  int get cardinality {
    int cardinality = 0;
      for (int i = 0; i < _dlength; i++) {
        if (_elements[2 * i] != 0 || _elements[2 * i + 1] != 0) {
          cardinality++;
        }
      }
    return cardinality;
  }

  Object get elements => _elements;

  bool all(Float64List value) {
    double epsilon = EPSILON;
    var x = new Float64List(2);
    var diff = new Float64List(2);
    for (int i = 0; i < _dlength; i++) {
      x[0] = _elements[2 * i];
      x[1] = _elements[2 * i + 1];
      diff[0] = (value[0] - x[0]).abs();
      diff[1] = (value[1] - x[1]).abs();
      if (((diff[0] != diff[0]) || (diff[1] != diff[1])) && ((((value[0] != value[0]) || (value[1] != value[1])) && ((x[0] != x[0]) || (x[1] != x[1])))) || (cmath.isEqual(value, x, epsilon))) {
        diff[0] = 0.0;
        diff[1] = 0.0;
      }
      if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
        return false;
      }
    }
    return true;
  }

  bool equals(AbstractComplexMatrix obj) {
    if (obj is DiagonalComplexMatrix) {
      DiagonalComplexMatrix other = obj;
      double epsilon = EPSILON;
      if (identical(this, obj)) {
        return true;
      }
      if (!(this != null && obj != null)) {
        return false;
      }
      if (columns != other.columns || rows != other.rows) {
        return false;
      }
      if ((_dindex != other._dindex) || (_dlength != other._dlength)) {
        return false;
      }
      Float64List otherElements = other._elements;
      var x = new Float64List(2);
      var value = new Float64List(2);
      var diff = new Float64List(2);
      for (int i = 0; i < _dlength; i++) {
        x[0] = _elements[2 * i];
        x[1] = _elements[2 * i + 1];
        value[0] = otherElements[2 * i];
        value[1] = otherElements[2 * i + 1];
        diff[0] = (value[0] - x[0]).abs();
        diff[1] = (value[1] - x[1]).abs();
        if (((diff[0] != diff[0]) || (diff[1] != diff[1])) && ((((value[0] != value[0]) || (value[1] != value[1])) && ((x[0] != x[0]) || (x[1] != x[1])))) || (cmath.isEqual(value, x, epsilon))) {
          diff[0] = 0.0;
          diff[1] = 0.0;
        }
        if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
          return false;
        }
      }
      return true;
    } else {
      return super.equals(obj);
    }
  }

  void forEachNonZero(final cfunc.IntIntComplexFunction fn) {
    var value = new Float64List(2);
    for (int i = 0; i < _dlength; i++) {
      value[0] = _elements[2 * i];
      value[1] = _elements[2 * i + 1];
      if (value[0] != 0 || value[1] != 0) {
        value = fn(i, i, value);
        _elements[2 * i] = value[0];
        _elements[2 * i + 1] = value[1];
      }
    }
  }

  int get diagonalLength => _dlength;

  int get diagonalIndex => _dindex;

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

  AbstractComplexVector like1D(int size) => new SparseComplexVector(size);

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
    int rowsA = rows;
    int columnsA = columns;
    if (transposeA) {
      rowsA = columns;
      columnsA = rows;
    }

    bool ignore = (z == null);
    if (z == null) {
      z = new ComplexVector(rowsA);
    }

    if (!(!isView && y is ComplexVector && z is ComplexVector)) {
      return super.mult(y, z, alpha, beta, transposeA);
    }

    if (columnsA != y.size || rowsA > z.size) {
      throw new ArgumentError("Incompatible args: " + ((transposeA ? dice() : this).toStringShort()) + ", " + y.toStringShort() + ", " + z.toStringShort());
    }

    if ((!ignore) && !((beta[0] == 1) && (beta[1] == 0))) {
      z.apply(cfunc.multiply(beta));
    }

    ComplexVector zz = z as ComplexVector;
    Float64List elementsZ = zz._elements;
    int strideZ = zz.stride;
    int zeroZ = z.index(0);

    ComplexVector yy = y as ComplexVector;
    Float64List elementsY = yy._elements;
    int strideY = yy.stride;
    int zeroY = y.index(0);

    if (elementsY == null || elementsZ == null) {
      throw new Error();
    }
    var elemA = new Float64List(2);
    var elemY = new Float64List(2);
    if (!transposeA) {
      if (_dindex >= 0) {
        for (int i = 0; i < _dlength; i++) {
          elemA[0] = _elements[2 * i];
          elemA[1] = _elements[2 * i + 1];
          elemY[0] = elementsY[2 * _dindex + zeroY + strideY * i];
          elemY[1] = elementsY[2 * _dindex + zeroY + strideY * i + 1];
          elemA = cmath.multiply(elemA, elemY);
          elemA = cmath.multiply(alpha, elemA);
          elementsZ[zeroZ + strideZ * i] += elemA[0];
          elementsZ[zeroZ + strideZ * i + 1] += elemA[1];
        }
      } else {
        for (int i = 0; i < _dlength; i++) {
          elemA[0] = _elements[2 * i];
          elemA[1] = _elements[2 * i + 1];
          elemY[0] = elementsY[zeroY + strideY * i];
          elemY[1] = elementsY[zeroY + strideY * i + 1];
          elemA = cmath.multiply(elemA, elemY);
          elemA = cmath.multiply(alpha, elemA);
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
          elemA = cmath.multiply(elemA, elemY);
          elemA = cmath.multiply(alpha, elemA);
          elementsZ[2 * _dindex + zeroZ + strideZ * i] += elemA[0];
          elementsZ[2 * _dindex + zeroZ + strideZ * i + 1] += elemA[1];
        }
      } else {
        for (int i = 0; i < _dlength; i++) {
          elemA[0] = _elements[2 * i];
          elemA[1] = -_elements[2 * i + 1];
          elemY[0] = elementsY[-2 * _dindex + zeroY + strideY * i];
          elemY[1] = elementsY[-2 * _dindex + zeroY + strideY * i + 1];
          elemA = cmath.multiply(elemA, elemY);
          elemA = cmath.multiply(alpha, elemA);
          elementsZ[zeroZ + strideZ * i] += elemA[0];
          elementsZ[zeroZ + strideZ * i + 1] += elemA[1];
        }
      }

    }
    return z;
  }

  AbstractComplexMatrix _getContent() => this;

  Object clone() {
    return new DiagonalComplexMatrix._internal(rows, columns, _elements, _dlength, _dindex);
  }
}
