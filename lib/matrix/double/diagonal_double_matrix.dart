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

/// Diagonal 2-d matrix holding [double] elements.
class DiagonalDoubleMatrix extends WrapperDoubleMatrix {

  /// The non zero elements of the matrix.
  Float64List _elements;

  /// Length of the diagonal.
  int _dlength;

  /// An m-by-n matrix A has m+n-1 diagonals. Since the [DiagonalDoubleMatrix]
  /// can have only one diagonal, [_dindex] is a value from interval
  /// `[-m+1, n-1]` that denotes which diagonal is stored.
  int _dindex;

  /// Constructs a matrix with a given number of rows and columns. All entries
  /// are initially `0`.
  DiagonalDoubleMatrix(int rows, int columns, [this._dindex = 0])
      : super._(rows, columns) {
    if ((_dindex < -rows + 1) || (_dindex > columns - 1)) {
      throw new ArgumentError("index is out of bounds");
    }
    _dlength;
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
    _elements = new Float64List(_dlength);
  }

  DiagonalDoubleMatrix._internal(
      int rows, int columns, this._elements, this._dlength, this._dindex)
      : super._(rows, columns);

  static DiagonalDoubleMatrix create(int rows, int columns) {
    return new DiagonalDoubleMatrix(rows, columns);
  }

  void apply(final func.DoubleFunction fn) {
    if (fn is DoubleMult) {
      // x[i] = mult*x[i]
      final double alpha = fn.multiplicator;
      if (alpha == 1) {
        return;
      }
      if (alpha == 0) {
        fill(0.0);
        return;
      }
      if (alpha != alpha) {
        fill(
            alpha); // the funny definition of isNaN(). This should better not happen.
        return;
      }
      for (int j = _dlength; --j >= 0;) {
        _elements[j] *= alpha;
      }
    } else {
      for (int j = _dlength; --j >= 0;) {
        _elements[j] = fn(_elements[j]);
      }
    }
  }

  void fill(double value) {
    for (int i = _dlength; --i >= 0;) {
      _elements[i] = value;
    }
  }

  void setAll(final Float64List values) {
    if (values.length != _dlength) {
      throw new ArgumentError(
          "Must have same length: length=${values.length} dlength=$_dlength");
    }
    for (int r = _dlength; --r >= 0;) {
      _elements[r] = values[r];
    }
  }

  void copyFrom(AbstractDoubleMatrix source) {
    // overriden for performance only
    if (source == this) {
      return; // nothing to do
    }
    checkShape(this, source);

    if (source is DiagonalDoubleMatrix) {
      DiagonalDoubleMatrix other = source;
      if ((_dindex != other._dindex) || (_dlength != other._dlength)) {
        throw new ArgumentError(
            "source is DiagonalDoubleMatrix with different diagonal stored.");
      }
      // quickest
      _elements.setAll(0, other._elements);
      return;
    } else {
      super.copyFrom(source);
    }
  }

  void assign(
      final AbstractDoubleMatrix y, final func.DoubleDoubleFunction fn) {
    checkShape(this, y);
    if (y is DiagonalDoubleMatrix) {
      DiagonalDoubleMatrix other = y;
      if ((_dindex != other._dindex) || (_dlength != other._dlength)) {
        throw new ArgumentError(
            "y is DiagonalDoubleMatrix with different diagonal stored.");
      }
      if (fn is DoublePlusMultSecond) {
        // x[i] = x[i] + alpha*y[i]
        double alpha = fn.multiplicator;
        if (alpha == 0) {
          return; // nothing to do
        }
      }
      final Float64List otherElements = other._elements;

      if (fn is DoublePlusMultSecond) {
        // x[i] = x[i] + alpha*y[i]
        var alpha = fn.multiplicator;
        if (alpha == 1) {
          for (int j = _dlength; --j >= 0;) {
            _elements[j] += otherElements[j];
          }
        } else {
          for (int j = _dlength; --j >= 0;) {
            _elements[j] = _elements[j] + alpha * otherElements[j];
          }
        }
      } else if (fn == func.mult) {
        // x[i] = x[i] * y[i]
        for (int j = _dlength; --j >= 0;) {
          _elements[j] = _elements[j] * otherElements[j];
        }
      } else if (fn == func.div) {
        // x[i] = x[i] /  y[i]
        for (int j = _dlength; --j >= 0;) {
          _elements[j] = _elements[j] / otherElements[j];
        }
      } else {
        for (int j = _dlength; --j >= 0;) {
          _elements[j] = fn(_elements[j], otherElements[j]);
        }
      }
      return;
    } else {
      super.assign(y, fn);
    }
  }

  int get cardinality {
    int cardinality = 0;
    for (int r = 0; r < _dlength; r++) {
      if (_elements[r] != 0) {
        cardinality++;
      }
    }
    return cardinality;
  }

  Object get elements => _elements;

  bool all(double value) {
    double epsilon = EPSILON;
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

  bool equals(AbstractDoubleMatrix obj) {
    if (obj is DiagonalDoubleMatrix) {
      DiagonalDoubleMatrix other = obj;
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

  void forEachNonZero(final func.IntIntDoubleFunction fn) {
    for (int j = _dlength; --j >= 0;) {
      double value = _elements[j];
      if (value != 0) {
        _elements[j] = fn(j, j, value);
      }
    }
  }

  /// The length of the diagonal.
  int get diagonalLength => _dlength;

  /// The index of the diagonal.
  int get diagonalIndex => _dindex;

  DoubleMatrixLocation max() {
    int location = 0;
    double maxValue = _elements[0];
    double elem;
    for (int r = 1; r < _dlength; r++) {
      elem = _elements[r];
      if (maxValue < elem) {
        maxValue = elem;
        location = r;
      }
    }

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
    return new DoubleMatrixLocation._(maxValue, rowLocation, columnLocation);
  }

  DoubleMatrixLocation min() {
    int location = 0;
    double minValue = _elements[0];
    double elem;
    for (int r = 1; r < _dlength; r++) {
      elem = _elements[r];
      if (minValue > elem) {
        minValue = elem;
        location = r;
      }
    }

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
    return new DoubleMatrixLocation._(minValue, rowLocation, columnLocation);
  }

  double get(int row, int column) {
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

  AbstractDoubleMatrix like2D(int rows, int columns) {
    return new SparseDoubleMatrix(rows, columns);
  }

  AbstractDoubleVector like1D(int size) {
    return new SparseDoubleVector(size);
  }

  void set(int row, int column, double value) {
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

  AbstractDoubleVector mult(AbstractDoubleVector y,
      [AbstractDoubleVector z = null, double alpha = 1.0, double beta = 0.0,
      final bool transposeA = false]) {
    int rowsA = rows;
    int columnsA = columns;
    if (transposeA) {
      rowsA = columns;
      columnsA = rows;
    }

    bool ignore = (z == null);
    if (z == null) {
      z = new DoubleVector(rowsA);
    }

    if (!(!isView && y is DoubleVector && z is DoubleVector)) {
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

    if ((!ignore) && ((beta) != 1)) {
      z.apply(func.multiply(beta));
    }

    DoubleVector zz = z as DoubleVector;
    final Float64List elementsZ = zz._elements;
    final int strideZ = zz.stride;
    final int zeroZ = z.index(0);

    DoubleVector yy = y as DoubleVector;
    final Float64List elementsY = yy._elements;
    final int strideY = yy.stride;
    final int zeroY = y.index(0);

    if (elementsY == null || elementsZ == null) {
      throw new Error();
    }
    if (!transposeA) {
      if (_dindex >= 0) {
        for (int i = _dlength; --i >= 0;) {
          elementsZ[zeroZ + strideZ * i] +=
              alpha * _elements[i] * elementsY[_dindex + zeroY + strideY * i];
        }
      } else {
        for (int i = _dlength; --i >= 0;) {
          elementsZ[-_dindex + zeroZ + strideZ * i] +=
              alpha * _elements[i] * elementsY[zeroY + strideY * i];
        }
      }
    } else {
      if (_dindex >= 0) {
        for (int i = _dlength; --i >= 0;) {
          elementsZ[_dindex + zeroZ + strideZ * i] +=
              alpha * _elements[i] * elementsY[zeroY + strideY * i];
        }
      } else {
        for (int i = _dlength; --i >= 0;) {
          elementsZ[zeroZ + strideZ * i] +=
              alpha * _elements[i] * elementsY[-_dindex + zeroY + strideY * i];
        }
      }
    }
    return z;
  }

  AbstractDoubleMatrix _getContent() => this;

  Object clone() {
    return new DiagonalDoubleMatrix._internal(
        rows, columns, _elements, _dlength, _dindex);
  }
}
