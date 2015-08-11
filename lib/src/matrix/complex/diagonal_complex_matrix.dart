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

/// Diagonal 2-d matrix holding [Complex] elements.
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
  DiagonalComplexMatrix(int rows, int columns, [this._dindex = 0])
      : super._(rows, columns) {
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

  DiagonalComplexMatrix._internal(
      int rows, int columns, this._elements, this._dlength, this._dindex)
      : super._(rows, columns);

  static DiagonalComplexMatrix create(int rows, int columns) {
    return new DiagonalComplexMatrix(rows, columns);
  }

  void apply(cfunc.ComplexComplexFunction fn) {
    if (fn is cfunc.ComplexMult) {
      // x[i] = mult*x[i]
      Complex alpha = fn.multiplicator;
      if (alpha.real == 1 && alpha.imaginary == 0) {
        return;
      }
      if (alpha.real == 0 && alpha.imaginary == 0) {
        fill(alpha.real, alpha.imaginary);
        return;
      }
      if (alpha.real != alpha.real || alpha.imaginary != alpha.imaginary) {
        fill(alpha.real, alpha.imaginary); // isNaN. This should not happen.
        return;
      }
      for (int j = 0; j < _dlength; j++) {
        var elem = new Complex(_elements[2 * j], _elements[2 * j + 1]);
        elem = elem * alpha;
        _elements[2 * j] = elem.real;
        _elements[2 * j + 1] = elem.imaginary;
      }
    } else {
      for (int j = 0; j < _dlength; j++) {
        var elem = new Complex(_elements[2 * j], _elements[2 * j + 1]);
        elem = fn(elem);
        _elements[2 * j] = elem.real;
        _elements[2 * j + 1] = elem.imaginary;
      }
    }
  }

  void fill(double re, double im) {
    for (int j = 0; j < _dlength; j++) {
      _elements[2 * j] = re;
      _elements[2 * j + 1] = im;
    }
  }

  void setAll(Float64List values) {
    if (values.length != 2 * _dlength) {
      throw new ArgumentError(
          "Must have same length: length=${values.length} 2*dlength=${2 * _dlength}");
    }
    for (int i = 0; i < _dlength; i++) {
      _elements[2 * i] = values[2 * i];
      _elements[2 * i + 1] = values[2 * i + 1];
    }
  }

  void copyFrom(ComplexMatrix source) {
    // overriden for performance only
    if (source == this) {
      return; // nothing to do
    }
    checkShape(this, source);

    if (source is DiagonalComplexMatrix) {
      DiagonalComplexMatrix other = source;
      if ((_dindex != other._dindex) || (_dlength != other._dlength)) {
        throw new ArgumentError(
            "source is DiagonalComplexMatrix with different diagonal stored.");
      }
      // quickest
      _elements.setAll(0, other._elements);
      return;
    } else {
      super.copyFrom(source);
      return;
    }
  }

  void assign(ComplexMatrix y,
      cfunc.ComplexComplexComplexFunction fn) {
    checkShape(this, y);
    if (y is DiagonalComplexMatrix) {
      DiagonalComplexMatrix other = y;
      if ((_dindex != other._dindex) || (_dlength != other._dlength)) {
        throw new ArgumentError(
            "y is DiagonalComplexMatrix with different diagonal stored.");
      }
      if (fn is cfunc.ComplexPlusMultSecond) {
        // x[i] = x[i] + alpha*y[i]
        Complex alpha = fn.multiplicator;
        if (alpha.real == 0 && alpha.imaginary == 0) {
          return; // nothing to do
        }
      }
      Float64List otherElements = other._elements;
      if (fn is cfunc.ComplexPlusMultSecond) {
        // x[i] = x[i] + alpha*y[i]
        Complex alpha = fn.multiplicator;
        if (alpha.real == 1 && alpha.imaginary == 0) {
          for (int j = 0; j < _dlength; j++) {
            _elements[2 * j] += otherElements[2 * j];
            _elements[2 * j + 1] += otherElements[2 * j + 1];
          }
        } else {
          for (int j = 0; j < _dlength; j++) {
            var elem =
                new Complex(otherElements[2 * j], otherElements[2 * j + 1]);
            elem = alpha * elem;
            _elements[2 * j] += elem.real;
            _elements[2 * j + 1] += elem.imaginary;
          }
        }
      } else if (fn == cfunc.mult) {
        // x[i] = x[i] * y[i]
        for (int j = 0; j < _dlength; j++) {
          var otherElem =
              new Complex(otherElements[2 * j], otherElements[2 * j + 1]);
          var elem = new Complex(_elements[2 * j], _elements[2 * j + 1]);
          elem = elem * otherElem;
          _elements[2 * j] = elem.real;
          _elements[2 * j + 1] = elem.imaginary;
        }
      } else if (fn == cfunc.div) {
        // x[i] = x[i] /  y[i]
        for (int j = 0; j < _dlength; j++) {
          var otherElem =
              new Complex(otherElements[2 * j], otherElements[2 * j + 1]);
          var elem = new Complex(_elements[2 * j], _elements[2 * j + 1]);
          elem = elem / otherElem;
          _elements[2 * j] = elem.real;
          _elements[2 * j + 1] = elem.imaginary;
        }
      } else {
        for (int j = 0; j < _dlength; j++) {
          var otherElem =
              new Complex(otherElements[2 * j], otherElements[2 * j + 1]);
          var elem = new Complex(_elements[2 * j], _elements[2 * j + 1]);
          elem = fn(elem, otherElem);
          _elements[2 * j] = elem.real;
          _elements[2 * j + 1] = elem.imaginary;
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

  bool equals(ComplexMatrix obj) {
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
      for (int i = 0; i < _dlength; i++) {
        var x = new Complex(_elements[2 * i], _elements[2 * i + 1]);
        var value = new Complex(otherElements[2 * i], otherElements[2 * i + 1]);
        var diff =
            new Complex((value.real - x.real).abs(), (value.imaginary - x.imaginary).abs());
        if (((diff.real != diff.real) || (diff.imaginary != diff.imaginary)) &&
                ((((value.real != value.real) ||
                        (value.imaginary != value.imaginary)) &&
                    ((x.real != x.real) || (x.imaginary != x.imaginary)))) ||
            (isEqual(value, x, epsilon))) {
          diff = Complex.ZERO;
        }
        if ((diff.real > epsilon) || (diff.imaginary > epsilon)) {
          return false;
        }
      }
      return true;
    } else {
      return super.equals(obj);
    }
  }

  void forEachNonZero(cfunc.IntIntComplexFunction fn) {
    for (int i = 0; i < _dlength; i++) {
      var value = new Complex(_elements[2 * i], _elements[2 * i + 1]);
      if (value.real != 0 || value.imaginary != 0) {
        value = fn(i, i, value);
        _elements[2 * i] = value.real;
        _elements[2 * i + 1] = value.imaginary;
      }
    }
  }

  int get diagonalLength => _dlength;

  int get diagonalIndex => _dindex;

  Complex get(int row, int column) {
    if (_dindex >= 0) {
      if (column < _dindex) {
        return Complex.ZERO;
      } else {
        if ((row < _dlength) && (row + _dindex == column)) {
          return new Complex(_elements[2 * row], _elements[2 * row + 1]);
        } else {
          return Complex.ZERO;
        }
      }
    } else {
      if (row < -_dindex) {
        return Complex.ZERO;
      } else {
        if ((column < _dlength) && (row + _dindex == column)) {
          return new Complex(_elements[2 * column], _elements[2 * column + 1]);
        } else {
          return Complex.ZERO;
        }
      }
    }
  }

  ComplexMatrix like2D(int rows, int columns) {
    return new SparseComplexMatrix(rows, columns);
  }

  ComplexVector like1D(int size) => new SparseComplexVector(size);

  void set(int row, int column, Complex value) {
    if (_dindex >= 0) {
      if (column < _dindex) {
        //do nothing
      } else {
        if ((row < _dlength) && (row + _dindex == column)) {
          _elements[2 * row] = value.real;
          _elements[2 * row + 1] = value.imaginary;
        } else {
          // do nothing
        }
      }
    } else {
      if (row < -_dindex) {
        //do nothing
      } else {
        if ((column < _dlength) && (row + _dindex == column)) {
          _elements[2 * column] = value.real;
          _elements[2 * column + 1] = value.imaginary;
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

  ComplexVector mult(ComplexVector y,
      [ComplexVector z = null, Complex alpha = null,
      Complex beta = null, bool transposeA = false]) {
    if (alpha == null) {
      alpha = Complex.ONE;
    }
    if (beta == null) {
      beta = (z == null ? Complex.ONE : Complex.ZERO);
    }
    int rowsA = rows;
    int columnsA = columns;
    if (transposeA) {
      rowsA = columns;
      columnsA = rows;
    }

    bool ignore = (z == null);
    if (z == null) {
      z = new DenseComplexVector(rowsA);
    }

    if (!(!isView && y is DenseComplexVector && z is DenseComplexVector)) {
      return super.mult(y, z, alpha, beta, transposeA);
    }

    if (columnsA != y.size || rowsA > z.size) {
      throw new ArgumentError("Incompatible args: " +
          ((transposeA ? dice() : this).toStringShort()) +
          ", " +
          y.toStringShort() +
          ", " +
          z.toStringShort());
    }

    if ((!ignore) && !((beta.real == 1) && (beta.imaginary == 0))) {
      z.apply(cfunc.multiply(beta));
    }

    DenseComplexVector zz = z as DenseComplexVector;
    Float64List elementsZ = zz._elements;
    int strideZ = zz.stride;
    int zeroZ = z.index(0);

    DenseComplexVector yy = y as DenseComplexVector;
    Float64List elementsY = yy._elements;
    int strideY = yy.stride;
    int zeroY = y.index(0);

    if (elementsY == null || elementsZ == null) {
      throw new Error();
    }
    if (!transposeA) {
      if (_dindex >= 0) {
        for (int i = 0; i < _dlength; i++) {
          var elemA = new Complex(_elements[2 * i], _elements[2 * i + 1]);
          var elemY = new Complex(elementsY[2 * _dindex + zeroY + strideY * i],
              elementsY[2 * _dindex + zeroY + strideY * i + 1]);
          elemA = elemA * elemY;
          elemA = alpha * elemA;
          elementsZ[zeroZ + strideZ * i] += elemA.real;
          elementsZ[zeroZ + strideZ * i + 1] += elemA.imaginary;
        }
      } else {
        for (int i = 0; i < _dlength; i++) {
          var elemA = new Complex(_elements[2 * i], _elements[2 * i + 1]);
          var elemY = new Complex(elementsY[zeroY + strideY * i],
              elementsY[zeroY + strideY * i + 1]);
          elemA = elemA * elemY;
          elemA = alpha * elemA;
          elementsZ[-2 * _dindex + zeroZ + strideZ * i] += elemA.real;
          elementsZ[-2 * _dindex + zeroZ + strideZ * i + 1] += elemA.imaginary;
        }
      }
    } else {
      if (_dindex >= 0) {
        for (int i = 0; i < _dlength; i++) {
          var elemA = new Complex(_elements[2 * i], -_elements[2 * i + 1]);
          var elemY = new Complex(elementsY[zeroY + strideY * i],
              elementsY[zeroY + strideY * i + 1]);
          elemA = elemA * elemY;
          elemA = alpha * elemA;
          elementsZ[2 * _dindex + zeroZ + strideZ * i] += elemA.real;
          elementsZ[2 * _dindex + zeroZ + strideZ * i + 1] += elemA.imaginary;
        }
      } else {
        for (int i = 0; i < _dlength; i++) {
          var elemA = new Complex(_elements[2 * i], -_elements[2 * i + 1]);
          var elemY = new Complex(elementsY[-2 * _dindex + zeroY + strideY * i],
              elementsY[-2 * _dindex + zeroY + strideY * i + 1]);
          elemA = elemA * elemY;
          elemA = alpha * elemA;
          elementsZ[zeroZ + strideZ * i] += elemA.real;
          elementsZ[zeroZ + strideZ * i + 1] += elemA.imaginary;
        }
      }
    }
    return z;
  }

  ComplexMatrix _getContent() => this;

  Object clone() {
    return new DiagonalComplexMatrix._internal(
        rows, columns, _elements, _dlength, _dindex);
  }
}
