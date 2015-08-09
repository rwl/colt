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

/// Diagonal 2-d matrix holding [int] elements.
class DiagonalIntMatrix extends WrapperIntMatrix {
  Int32List _elements;

  /// Length of the diagonal
  int _dlength;

  /// An m-by-n matrix A has m+n-1 diagonals. Since the [DiagonalIntMatrix]
  /// can have only one diagonal, dindex is a value from interval `[-m+1, n-1]`
  /// that denotes which diagonal is stored.
  int _dindex;

  /// Constructs a matrix with a given number of rows and columns. All entries
  /// are initially `0`.
  DiagonalIntMatrix(int rows, int columns, int dindex)
      : super._(rows, columns) {
    if ((dindex < -rows + 1) || (dindex > columns - 1)) {
      throw new ArgumentError("index is out of bounds");
    } else {
      _dindex = dindex;
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

  void apply(final ifunc.IntFunction fn) {
    if (fn is ifunc.IntMult) {
      // x[i] = mult*x[i]
      int alpha = fn.multiplicator;
      if (alpha == 1) {
        return;
      }
      if (alpha == 0) {
        fill(0);
        return;
      }
      if (alpha != alpha) {
        // The definition of isNaN(). This should not happen.
        fill(alpha);
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
    return;
  }

  void fill(int value) {
    for (int i = _dlength; --i >= 0;) {
      _elements[i] = value;
    }
    return;
  }

  void setAll(final Int32List values) {
    if (values.length != _dlength) {
      throw new ArgumentError(
          "Must have same length: length=${values.length} dlength=$_dlength");
    }
    for (int r = _dlength; --r >= 0;) {
      _elements[r] = values[r];
    }
    return;
  }

  void copyFrom(IntMatrix source) {
    // overriden for performance only
    if (source == this) {
      return; // nothing to do
    }
    checkShape(this, source);

    if (source is DiagonalIntMatrix) {
      DiagonalIntMatrix other = source;
      if ((_dindex != other._dindex) || (_dlength != other._dlength)) {
        throw new ArgumentError(
            "source is DiagonalIntMatrix with different diagonal stored.");
      }
      _elements.setAll(0, other._elements);
      return;
    } else {
      super.copyFrom(source);
      return;
    }
  }

  void assign(final IntMatrix y, final ifunc.IntIntFunction fn) {
    checkShape(this, y);
    if (y is DiagonalIntMatrix) {
      DiagonalIntMatrix other = y;
      if ((_dindex != other._dindex) || (_dlength != other._dlength)) {
        throw new ArgumentError(
            "y is DiagonalIntMatrix2D with different diagonal stored.");
      }
      if (fn is ifunc.IntPlusMultSecond) {
        // x[i] = x[i] + alpha*y[i]
        final int alpha = fn.multiplicator;
        if (alpha == 0) {
          return; // nothing to do
        }
      }
      Int32List otherElements = other._elements;
      if (fn is ifunc.IntPlusMultSecond) {
        // x[i] = x[i] + alpha*y[i]
        final int alpha = fn.multiplicator;
        if (alpha == 1) {
          for (int j = _dlength; --j >= 0;) {
            _elements[j] += otherElements[j];
          }
        } else {
          for (int j = _dlength; --j >= 0;) {
            _elements[j] = _elements[j] + alpha * otherElements[j];
          }
        }
      } else if (fn == ifunc.mult) {
        // x[i] = x[i] * y[i]
        for (int j = _dlength; --j >= 0;) {
          _elements[j] = _elements[j] * otherElements[j];
        }
      } else if (fn == ifunc.div) {
        // x[i] = x[i] /  y[i]
        for (int j = _dlength; --j >= 0;) {
          _elements[j] = _elements[j] ~/ otherElements[j];
        }
      } else {
        for (int j = _dlength; --j >= 0;) {
          _elements[j] = fn(_elements[j], otherElements[j]);
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
    for (int r = 0; r < _dlength; r++) {
      if (_elements[r] != 0) cardinality++;
    }
    return cardinality;
  }

  Object get elements => _elements;

  bool equals(IntMatrix obj) {
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
    for (int j = _dlength; --j >= 0;) {
      int value = _elements[j];
      if (value != 0) {
        _elements[j] = function(j, j, value);
      }
    }
    return;
  }

  int get diagonalLength => _dlength;

  int get diagonalIndex => _dindex;

  IntMatrixLocation max() {
    int location = 0;
    int maxValue = _elements[0];
    for (int r = 1; r < _dlength; r++) {
      var elem = _elements[r];
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
    return new IntMatrixLocation._(maxValue, rowLocation, columnLocation);
  }

  IntMatrixLocation min() {
    int location = 0;
    int minValue = _elements[0];
    for (int r = 1; r < _dlength; r++) {
      var elem = _elements[r];
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
    return new IntMatrixLocation._(minValue, rowLocation, columnLocation);
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

  IntMatrix like2D(int rows, int columns) {
    return new SparseIntMatrix(rows, columns);
  }

  IntVector like1D(int size) => new SparseIntVector(size);

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

  IntVector mult(IntVector y, [IntVector z = null,
      final int alpha = 1, int beta = null, final bool transposeA = false]) {
    if (beta == null) {
      beta = z == null ? 1 : 0;
    }
    int rowsA = rows;
    int columnsA = columns;
    if (transposeA) {
      rowsA = columns;
      columnsA = rows;
    }

    bool ignore = (z == null);
    if (z == null) {
      z = new DenseIntVector(rowsA);
    }

    if (!(!isView && y is DenseIntVector && z is DenseIntVector)) {
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
      z.apply(ifunc.multiply(beta));
    }

    DenseIntVector zz = z as DenseIntVector;
    Int32List elementsZ = zz._elements;
    int strideZ = zz.stride;
    int zeroZ = z.index(0);

    DenseIntVector yy = y as DenseIntVector;
    Int32List elementsY = yy._elements;
    int strideY = yy.stride;
    int zeroY = y.index(0);

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

  IntMatrix _getContent() => this;

  Object clone() => new DiagonalIntMatrix(rows, columns, _dindex);
}
