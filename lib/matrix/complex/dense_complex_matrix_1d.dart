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
 * Dense 1-d matrix (aka <i>vector</i>) holding <tt>complex</tt> elements.
 * <p>
 * Internally holds one single contiguous one-dimensional array. Complex data is
 * represented by 2 double values in sequence, i.e. elements[zero + 2 * k *
 * stride] constitute real part and elements[zero + 2 * k * stride + 1]
 * constitute imaginary part (k=0,...,size()-1).
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
class DenseDComplexMatrix1D extends DComplexMatrix1D {

  /**
   * The elements of this matrix. Complex data is represented by 2 double
   * values in sequence, i.e. elements[zero + 2 * k * stride] constitute real
   * part and elements[zero + 2 * k * stride] constitute imaginary part
   * (k=0,...,size()-1).
   */
  Float64List _elements;

  /**
   * Constructs a matrix with a copy of the given values. The values are
   * copied. So subsequent changes in <tt>values</tt> are not reflected in the
   * matrix, and vice-versa. Due to the fact that complex data is represented
   * by 2 double values in sequence: the real and imaginary parts, the size of
   * new matrix will be equal to values.length / 2.
   *
   * @param values
   *            The values to be filled into the new matrix.
   */
  factory DenseDComplexMatrix1D.fromValues(Float64List values) {
    return new DenseDComplexMatrix1D(values.length ~/ 2)..assignValues(values);
  }

  /**
   * Constructs a complex matrix with the same size as <tt>realPart</tt>
   * matrix and fills the real part of this matrix with elements of
   * <tt>realPart</tt>.
   *
   * @param realPart
   *            a real matrix whose elements become a real part of this matrix
   * @throws IllegalArgumentException
   *             if <tt>size<0</tt>.
   */
  factory DenseDComplexMatrix1D.fromRealPart(DoubleMatrix1D realPart) {
    return new DenseDComplexMatrix1D(realPart.size())..assignReal(realPart);
  }

  /**
   * Constructs a complex matrix with the same size as <tt>realPart</tt>
   * matrix and fills the real part of this matrix with elements of
   * <tt>realPart</tt> and fills the imaginary part of this matrix with
   * elements of <tt>imaginaryPart</tt>.
   *
   * @param realPart
   *            a real matrix whose elements become a real part of this matrix
   * @param imaginaryPart
   *            a imaginary matrix whose elements become a imaginary part of
   *            this matrix
   * @throws IllegalArgumentException
   *             if <tt>size<0</tt>.
   */
  factory DenseDComplexMatrix1D.fromParts(DoubleMatrix1D realPart, DoubleMatrix1D imaginaryPart) {
    return new DenseDComplexMatrix1D(realPart.size())
      ..assignReal(realPart)
      ..assignImaginary(imaginaryPart);
  }

  /**
   * Constructs a matrix with the given parameters.
   *
   * @param size
   *            the number of cells the matrix shall have.
   * @param elements
   *            the cells.
   * @param zero
   *            the index of the first element.
   * @param stride
   *            the number of indexes between any two elements, i.e.
   *            <tt>index(i+1)-index(i)</tt>.
   * @param isNoView
   *            if false then the view is constructed
   * @throws IllegalArgumentException
   *             if <tt>size<0</tt>.
   */
  DenseDComplexMatrix1D(int size, [Float64List elements = null, int zero = 0, int stride = 1, bool isNoView = true]) {
    if (elements == null) {
      elements = new Float64List(2 * size);
    }
    _setUp(size, zero, stride);
    this._elements = elements;
    this._isNoView = isNoView;
  }

  Float64List aggregate(final cfunc.DComplexDComplexDComplexFunction aggr, final cfunc.DComplexDComplexFunction f) {
    Float64List b = new Float64List(2);
    if (_size == 0) {
      b[0] = double.NAN;
      b[1] = double.NAN;
      return b;
    }
    Float64List a = null;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = _zero + firstIdx * _stride;
          Float64List a = f([_elements[idx], _elements[idx + 1]]);
          for (int i = firstIdx + 1; i < lastIdx; i++) {
            idx += _stride;
            a = aggr(a, f([_elements[idx], _elements[idx + 1]]));
          }
          return a;
        });
      }
      a = ConcurrencyUtils.waitForCompletion(futures, aggr);
    } else {*/
      a = f([_elements[_zero], _elements[_zero + 1]]);
      int idx = _zero;
      for (int i = 1; i < _size; i++) {
        idx += _stride;
        a = aggr(a, f([_elements[idx], _elements[idx + 1]]));
      }
    //}
    return a;
  }

  Float64List aggregateMatrix(final DComplexMatrix1D other, final cfunc.DComplexDComplexDComplexFunction aggr, final cfunc.DComplexDComplexDComplexFunction f) {
    if (!(other is DenseDComplexMatrix1D)) {
      return super.aggregateMatrix(other, aggr, f);
    }
    checkSize(other);
    if (_size == 0) {
      Float64List b = new Float64List(2);
      b[0] = double.NAN;
      b[1] = double.NAN;
      return b;
    }
    final int zeroOther = other.index(0);
    final int strideOther = other.stride();
    final Float64List elemsOther = other.elements() as Float64List;
    Float64List a = null;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = _zero + firstIdx * _stride;
          int idxOther = zeroOther + firstIdx * strideOther;
          Float64List a = f([_elements[idx], _elements[idx + 1]], [elemsOther[idxOther], elemsOther[idxOther + 1]]);
          for (int i = firstIdx + 1; i < lastIdx; i++) {
            idx += _stride;
            idxOther += strideOther;
            a = aggr(a, f([_elements[idx], _elements[idx + 1]], [elemsOther[idxOther], elemsOther[idxOther + 1]]));
          }
          return a;
        });
      }
      a = ConcurrencyUtils.waitForCompletion(futures, aggr);
    } else {*/
      int idx = _zero;
      int idxOther = zeroOther;
      a = f([_elements[_zero], _elements[_zero + 1]], [elemsOther[zeroOther], elemsOther[zeroOther + 1]]);
      for (int i = 1; i < _size; i++) {
        idx += _stride;
        idxOther += strideOther;
        a = aggr(a, f([_elements[idx], _elements[idx + 1]], [elemsOther[idxOther], elemsOther[idxOther + 1]]));
      }
    //}
    return a;
  }

  DComplexMatrix1D assign(final cfunc.DComplexDComplexFunction function) {
    if (this._elements == null) {
      throw new Error();
    }
    if (function is DComplexMult) {
      Float64List multiplicator = (function as DComplexMult).multiplicator;
      if (multiplicator[0] == 1 && multiplicator[1] == 0) {
        return this;
      }
    }
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          Float64List tmp = new Float64List(2);
          int idx = _zero + firstIdx * _stride;
          if (function is DComplexMult) {
            Float64List multiplicator = (function as DComplexMult).multiplicator;
            for (int k = firstIdx; k < lastIdx; k++) {
              _elements[idx] = _elements[idx] * multiplicator[0] - _elements[idx + 1] * multiplicator[1];
              _elements[idx + 1] = _elements[idx + 1] * multiplicator[0] + _elements[idx] * multiplicator[1];
              idx += _stride;
            }
          } else {
            for (int k = firstIdx; k < lastIdx; k++) {
              tmp[0] = _elements[idx];
              tmp[1] = _elements[idx + 1];
              tmp = function(tmp);
              _elements[idx] = tmp[0];
              _elements[idx + 1] = tmp[1];
              idx += _stride;
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      Float64List tmp = new Float64List(2);
      int idx = _zero;
      if (function is DComplexMult) {
        Float64List multiplicator = (function as DComplexMult).multiplicator;
        for (int k = 0; k < _size; k++) {
          _elements[idx] = _elements[idx] * multiplicator[0] - _elements[idx + 1] * multiplicator[1];
          _elements[idx + 1] = _elements[idx + 1] * multiplicator[0] + _elements[idx] * multiplicator[1];
          idx += _stride;
        }
      } else {
        for (int k = 0; k < _size; k++) {
          tmp[0] = _elements[idx];
          tmp[1] = _elements[idx + 1];
          tmp = function(tmp);
          _elements[idx] = tmp[0];
          _elements[idx + 1] = tmp[1];
          idx += _stride;
        }
      }
    //}
    return this;
  }

  DComplexMatrix1D assignProc(final cfunc.DComplexProcedure cond, final cfunc.DComplexDComplexFunction function) {
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          Float64List elem = new Float64List(2);
          int idx = _zero + firstIdx * _stride;
          for (int i = firstIdx; i < lastIdx; i++) {
            elem[0] = _elements[idx];
            elem[1] = _elements[idx + 1];
            if (cond(elem) == true) {
              elem = function(elem);
              _elements[idx] = elem[0];
              _elements[idx + 1] = elem[1];
            }
            idx += _stride;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      Float64List elem = new Float64List(2);
      int idx = _zero;
      for (int i = 0; i < _size; i++) {
        elem[0] = _elements[idx];
        elem[1] = _elements[idx + 1];
        if (cond(elem) == true) {
          elem = function(elem);
          _elements[idx] = elem[0];
          _elements[idx + 1] = elem[1];
        }
        idx += _stride;
      }
    //}
    return this;
  }

  DComplexMatrix1D assignProcValue(final cfunc.DComplexProcedure cond, final Float64List value) {
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          Float64List elem = new Float64List(2);
          int idx = _zero + firstIdx * _stride;
          for (int i = firstIdx; i < lastIdx; i++) {
            elem[0] = _elements[idx];
            elem[1] = _elements[idx + 1];
            if (cond(elem) == true) {
              _elements[idx] = value[0];
              _elements[idx + 1] = value[1];
            }
            idx += _stride;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      Float64List elem = new Float64List(2);
      int idx = _zero;
      for (int i = 0; i < _size; i++) {
        elem[0] = _elements[idx];
        elem[1] = _elements[idx + 1];
        if (cond(elem) == true) {
          _elements[idx] = value[0];
          _elements[idx + 1] = value[1];
        }
        idx += _stride;
      }
    //}
    return this;
  }

  DComplexMatrix1D assignRealFunc(final cfunc.DComplexRealFunction function) {
    if (this._elements == null) {
      throw new Error();
    }
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = _zero + firstIdx * _stride;
          if (function == DComplexFunctions.abs) {
            for (int k = firstIdx; k < lastIdx; k++) {
              double absX = Math.abs(_elements[idx]);
              double absY = Math.abs(_elements[idx + 1]);
              if (absX == 0.0 && absY == 0.0) {
                _elements[idx] = 0;
              } else if (absX >= absY) {
                double d = _elements[idx + 1] / _elements[idx];
                _elements[idx] = absX * Math.sqrt(1.0 + d * d);
              } else {
                double d = _elements[idx] / _elements[idx + 1];
                _elements[idx] = absY * Math.sqrt(1.0 + d * d);
              }
              _elements[idx + 1] = 0;
              idx += _stride;
            }
          } else {
            Float64List tmp = new Float64List(2);
            for (int k = firstIdx; k < lastIdx; k++) {
              tmp[0] = _elements[idx];
              tmp[1] = _elements[idx + 1];
              tmp[0] = function(tmp);
              _elements[idx] = tmp[0];
              _elements[idx + 1] = 0;
              idx += _stride;
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      int idx = _zero;
      if (function == cfunc.abs) {
        for (int k = 0; k < _size; k++) {
          double absX = _elements[idx].abs();
          double absY = _elements[idx + 1].abs();
          if (absX == 0.0 && absY == 0.0) {
            _elements[idx] = 0.0;
          } else if (absX >= absY) {
            double d = _elements[idx + 1] / _elements[idx];
            _elements[idx] = absX * Math.sqrt(1.0 + d * d);
          } else {
            double d = _elements[idx] / _elements[idx + 1];
            _elements[idx] = absY * Math.sqrt(1.0 + d * d);
          }
          _elements[idx + 1] = 0.0;
          idx += _stride;
        }
      } else {
        Float64List tmp = new Float64List(2);
        for (int k = 0; k < _size; k++) {
          tmp[0] = _elements[idx];
          tmp[1] = _elements[idx + 1];
          tmp[0] = function(tmp);
          _elements[idx] = tmp[0];
          _elements[idx + 1] = 0.0;
          idx += _stride;
        }
      }
    //}
    return this;
  }

  DComplexMatrix1D assignMatrix(DComplexMatrix1D source) {
    if (!(source is DenseDComplexMatrix1D)) {
      return super.assignMatrix(source);
    }
    DenseDComplexMatrix1D other = source as DenseDComplexMatrix1D;
    if (other == this) return this;
    checkSize(other);
    if (_isNoView && other._isNoView) { // quickest
      //System.arraycopy(other._elements, 0, this._elements, 0, this._elements.length);
      this._elements.setAll(0, other._elements);
      return this;
    }
    if (_haveSharedCells(other)) {
      DComplexMatrix1D c = other.copy();
      if (!(c is DenseDComplexMatrix1D)) { // should not happen
        return super.assignMatrix(source);
      }
      other = c as DenseDComplexMatrix1D;
    }

    final Float64List elemsOther = other._elements;
    if (_elements == null || elemsOther == null) {
      throw new Error();
    }
    final int strideOther = other._stride;
    final int zeroOther = other.index(0);

    int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = _zero + firstIdx * _stride;
          int idxOther = zeroOther + firstIdx * strideOther;
          for (int k = firstIdx; k < lastIdx; k++) {
            _elements[idx] = elemsOther[idxOther];
            _elements[idx + 1] = elemsOther[idxOther + 1];
            idx += _stride;
            idxOther += strideOther;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {
      int idx = _zero;
      int idxOther = zeroOther;
      for (int k = 0; k < _size; k++) {
        _elements[idx] = elemsOther[idxOther];
        _elements[idx + 1] = elemsOther[idxOther + 1];
        idx += _stride;
        idxOther += strideOther;
      }
    }
    return this;
  }

  DComplexMatrix1D assignFunc(DComplexMatrix1D y, final cfunc.DComplexDComplexDComplexFunction function) {
    if (!(y is DenseDComplexMatrix1D)) {
      return super.assignMatrixFunc(y, function);
    }
    checkSize(y);
    final Float64List elemsOther = y.elements() as Float64List;
    final int zeroOther = y.index(0);
    final int strideOther = y.stride();

    if (_elements == null || elemsOther == null) {
      throw new Error();
    }
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = _zero + firstIdx * _stride;
          int idxOther = zeroOther + firstIdx * strideOther;
          if (function == cfunc.plus) {
            for (int k = firstIdx; k < lastIdx; k++) {
              _elements[idx] += elemsOther[idxOther];
              _elements[idx + 1] += elemsOther[idxOther + 1];
              idx += _stride;
              idxOther += strideOther;
            }
          } else if (function == cfunc.minus) {
            for (int k = firstIdx; k < lastIdx; k++) {
              _elements[idx] -= elemsOther[idxOther];
              _elements[idx + 1] -= elemsOther[idxOther + 1];
              idx += _stride;
              idxOther += strideOther;
            }
          } else if (function == cfunc.div) {
            Float64List tmp = new Float64List(2);
            for (int k = firstIdx; k < lastIdx; k++) {
              double re = elemsOther[idxOther];
              double im = elemsOther[idxOther + 1];
              double scalar;
              if (re.abs() >= im.abs()) {
                scalar = (1.0 / (re + im * (im / re)));
                tmp[0] = scalar * (_elements[idx] + _elements[idx + 1] * (im / re));
                tmp[1] = scalar * (_elements[idx + 1] - _elements[idx] * (im / re));
              } else {
                scalar = (1.0 / (re * (re / im) + im));
                tmp[0] = scalar * (_elements[idx] * (re / im) + _elements[idx + 1]);
                tmp[1] = scalar * (_elements[idx + 1] * (re / im) - _elements[idx]);
              }
              _elements[idx] = tmp[0];
              _elements[idx + 1] = tmp[1];
              idx += _stride;
              idxOther += strideOther;
            }
          } else if (function == cfunc.mult) {
            Float64List tmp = new Float64List(2);
            for (int k = firstIdx; k < lastIdx; k++) {
              tmp[0] = _elements[idx] * elemsOther[idxOther] - _elements[idx + 1] * elemsOther[idxOther + 1];
              tmp[1] = _elements[idx + 1] * elemsOther[idxOther] + _elements[idx] * elemsOther[idxOther + 1];
              _elements[idx] = tmp[0];
              _elements[idx + 1] = tmp[1];
              idx += _stride;
              idxOther += strideOther;
            }
          } else if (function == cfunc.multConjFirst) {
            Float64List tmp = new Float64List(2);
            for (int k = firstIdx; k < lastIdx; k++) {
              tmp[0] = _elements[idx] * elemsOther[idxOther] + _elements[idx + 1] * elemsOther[idxOther + 1];
              tmp[1] = -_elements[idx + 1] * elemsOther[idxOther] + _elements[idx] * elemsOther[idxOther + 1];
              _elements[idx] = tmp[0];
              _elements[idx + 1] = tmp[1];
              idx += _stride;
              idxOther += strideOther;
            }
          } else if (function == cfunc.multConjSecond) {
            Float64List tmp = new Float64List(2);
            for (int k = firstIdx; k < lastIdx; k++) {
              tmp[0] = _elements[idx] * elemsOther[idxOther] + _elements[idx + 1] * elemsOther[idxOther + 1];
              tmp[1] = _elements[idx + 1] * elemsOther[idxOther] - _elements[idx] * elemsOther[idxOther + 1];
              _elements[idx] = tmp[0];
              _elements[idx + 1] = tmp[1];
              idx += _stride;
              idxOther += strideOther;
            }
          } else {
            Float64List tmp1 = new Float64List(2);
            Float64List tmp2 = new Float64List(2);
            for (int k = firstIdx; k < lastIdx; k++) {
              tmp1[0] = _elements[idx];
              tmp1[1] = _elements[idx + 1];
              tmp2[0] = elemsOther[idxOther];
              tmp2[1] = elemsOther[idxOther + 1];
              tmp1 = function(tmp1, tmp2);
              _elements[idx] = tmp1[0];
              _elements[idx + 1] = tmp1[1];
              idx += _stride;
              idxOther += strideOther;
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      int idx = _zero;
      int idxOther = zeroOther;
      if (function == cfunc.plus) {
        for (int k = 0; k < _size; k++) {
          _elements[idx] += elemsOther[idxOther];
          _elements[idx + 1] += elemsOther[idxOther + 1];
          idx += _stride;
          idxOther += strideOther;
        }
      } else if (function == cfunc.minus) {
        for (int k = 0; k < _size; k++) {
          _elements[idx] -= elemsOther[idxOther];
          _elements[idx + 1] -= elemsOther[idxOther + 1];
          idx += _stride;
          idxOther += strideOther;
        }
      } else if (function == cfunc.div) {
        Float64List tmp = new Float64List(2);
        for (int k = 0; k < _size; k++) {
          double re = elemsOther[idxOther];
          double im = elemsOther[idxOther + 1];
          double scalar;
          if (re.abs() >= im.abs()) {
            scalar = (1.0 / (re + im * (im / re)));
            tmp[0] = scalar * (_elements[idx] + _elements[idx + 1] * (im / re));
            tmp[1] = scalar * (_elements[idx + 1] - _elements[idx] * (im / re));
          } else {
            scalar = (1.0 / (re * (re / im) + im));
            tmp[0] = scalar * (_elements[idx] * (re / im) + _elements[idx + 1]);
            tmp[1] = scalar * (_elements[idx + 1] * (re / im) - _elements[idx]);
          }
          _elements[idx] = tmp[0];
          _elements[idx + 1] = tmp[1];
          idx += _stride;
          idxOther += strideOther;
        }
      } else if (function == cfunc.mult) {
        Float64List tmp = new Float64List(2);
        for (int k = 0; k < _size; k++) {
          tmp[0] = _elements[idx] * elemsOther[idxOther] - _elements[idx + 1] * elemsOther[idxOther + 1];
          tmp[1] = _elements[idx + 1] * elemsOther[idxOther] + _elements[idx] * elemsOther[idxOther + 1];
          _elements[idx] = tmp[0];
          _elements[idx + 1] = tmp[1];
          idx += _stride;
          idxOther += strideOther;
        }
      } else if (function == cfunc.multConjFirst) {
        Float64List tmp = new Float64List(2);
        for (int k = 0; k < _size; k++) {
          tmp[0] = _elements[idx] * elemsOther[idxOther] + _elements[idx + 1] * elemsOther[idxOther + 1];
          tmp[1] = -_elements[idx + 1] * elemsOther[idxOther] + _elements[idx] * elemsOther[idxOther + 1];
          _elements[idx] = tmp[0];
          _elements[idx + 1] = tmp[1];
          idx += _stride;
          idxOther += strideOther;
        }
      } else if (function == cfunc.multConjSecond) {
        Float64List tmp = new Float64List(2);
        for (int k = 0; k < _size; k++) {
          tmp[0] = _elements[idx] * elemsOther[idxOther] + _elements[idx + 1] * elemsOther[idxOther + 1];
          tmp[1] = _elements[idx + 1] * elemsOther[idxOther] - _elements[idx] * elemsOther[idxOther + 1];
          _elements[idx] = tmp[0];
          _elements[idx + 1] = tmp[1];
          idx += _stride;
          idxOther += strideOther;
        }
      } else {
        Float64List tmp1 = new Float64List(2);
        Float64List tmp2 = new Float64List(2);
        for (int k = 0; k < _size; k++) {
          tmp1[0] = _elements[idx];
          tmp1[1] = _elements[idx + 1];
          tmp2[0] = elemsOther[idxOther];
          tmp2[1] = elemsOther[idxOther + 1];
          tmp1 = function(tmp1, tmp2);
          _elements[idx] = tmp1[0];
          _elements[idx + 1] = tmp1[1];
          idx += _stride;
          idxOther += strideOther;
        }
      }
    //}
    return this;
  }

  DComplexMatrix1D assignValue(final double re, final double im) {
    int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = _zero + firstIdx * _stride;
          for (int k = firstIdx; k < lastIdx; k++) {
            _elements[idx] = re;
            _elements[idx + 1] = im;
            idx += _stride;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {
      int idx = _zero;
      for (int i = 0; i < _size; i++) {
        this._elements[idx] = re;
        this._elements[idx + 1] = im;
        idx += _stride;
      }
    }
    return this;
  }

  DComplexMatrix1D assignValues(Float64List values) {
    if (_isNoView) {
      if (values.length != 2 * _size) {
        throw new ArgumentError("The length of values[] must be equal to 2*size()=${2 * size()}");
      }
      //System.arraycopy(values, 0, this._elements, 0, values.length);
      this._elements.setAll(0, values);
    } else {
      super.assignValues(values);
    }
    return this;
  }

  DComplexMatrix1D assignImaginary(final DoubleMatrix1D other) {
    if (!(other is DenseDoubleMatrix1D)) {
      return super.assignImaginary(other);
    }
    checkSize(other);
    final int zeroOther = other.index(0);
    final int strideOther = other.stride();
    final Float64List elemsOther = other.elements() as Float64List;
    int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = _zero + firstIdx * _stride;
          int idxOther = zeroOther + firstIdx * strideOther;
          for (int i = firstIdx; i < lastIdx; i++) {
            _elements[idx + 1] = elemsOther[idxOther];
            idx += _stride;
            idxOther += strideOther;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {
      int idx = _zero;
      int idxOther = zeroOther;
      for (int i = 0; i < _size; i++) {
        _elements[idx + 1] = elemsOther[idxOther];
        idx += _stride;
        idxOther += strideOther;
      }
    }
    return this;
  }

  DComplexMatrix1D assignReal(final DoubleMatrix1D other) {
    if (!(other is DenseDoubleMatrix1D)) {
      return super.assignReal(other);
    }
    checkSize(other);
    final int zeroOther = other.index(0);
    final int strideOther = other.stride();
    final Float64List elemsOther = other.elements() as Float64List;
    int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = _zero + firstIdx * _stride;
          int idxOther = zeroOther + firstIdx * strideOther;
          for (int i = firstIdx; i < lastIdx; i++) {
            _elements[idx] = elemsOther[idxOther];
            idx += _stride;
            idxOther += strideOther;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {
      int idx = _zero;
      int idxOther = zeroOther;
      for (int i = 0; i < _size; i++) {
        _elements[idx] = elemsOther[idxOther];
        idx += _stride;
        idxOther += strideOther;
      }
    }
    return this;
  }

  Float64List elements() {
    return _elements;
  }

  DoubleMatrix1D getImaginaryPart() {
    final DenseDoubleMatrix1D Im = new DenseDoubleMatrix1D(_size);
    final Float64List elemsOther = Im.elements();
    final int zeroOther = Im.index(0);
    final int strideOther = Im.stride();
    int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = _zero + firstIdx * _stride;
          int idxOther = zeroOther + firstIdx * strideOther;
          for (int k = firstIdx; k < lastIdx; k++) {
            elemsOther[idxOther] = _elements[idx + 1];
            idx += _stride;
            idxOther += strideOther;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {
      int idx = _zero;
      int idxOther = zeroOther;
      for (int i = 0; i < _size; i++) {
        elemsOther[idxOther] = _elements[idx + 1];
        idx += _stride;
        idxOther += strideOther;
      }
    }
    return Im;
  }

  void getNonZeros(final List<int> indexList, final List<Float64List> valueList) {
    indexList.clear();
    valueList.clear();
    int s = size();

    int idx = _zero;
    for (int k = 0; k < s; k++) {
      Float64List value = new Float64List(2);
      value[0] = _elements[idx];
      value[1] = _elements[idx + 1];
      if (value[0] != 0 || value[1] != 0) {
        indexList.add(k);
        valueList.add(value);
      }
      idx += _stride;
    }
  }

  Float64List getQuick(int index) {
    int idx = _zero + index * _stride;
    return new Float64List.fromList([_elements[idx], _elements[idx + 1]]);
  }

  DoubleMatrix1D getRealPart() {
    final DenseDoubleMatrix1D R = new DenseDoubleMatrix1D(_size);
    final Float64List elemsOther = R.elements();
    final int zeroOther = R.index(0);
    final int strideOther = R.stride();
    int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = _zero + firstIdx * _stride;
          int idxOther = zeroOther + firstIdx * strideOther;
          for (int k = firstIdx; k < lastIdx; k++) {
            elemsOther[idxOther] = _elements[idx];
            idx += _stride;
            idxOther += strideOther;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {
      int idx = _zero;
      int idxOther = zeroOther;
      for (int i = 0; i < _size; i++) {
        elemsOther[idxOther] = _elements[idx];
        idx += _stride;
        idxOther += strideOther;
      }
    }
    return R;
  }

  DComplexMatrix1D like1D(int size) {
    return new DenseDComplexMatrix1D(size);
  }

  DComplexMatrix2D like2D(int rows, int columns) {
    return new DenseDComplexMatrix2D(rows, columns);
  }

  DComplexMatrix2D reshape(final int rows, final int columns) {
    if (rows * columns != _size) {
      throw new ArgumentError("rows*columns != size");
    }
    DComplexMatrix2D M = new DenseDComplexMatrix2D(rows, columns);
    final Float64List elemsOther = M.elements() as Float64List;
    final int zeroOther = M.index(0, 0);
    final int rowStrideOther = M.rowStride();
    final int columnStrideOther = M.columnStride();
    int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = columns ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstColumn = j * k;
        final int lastColumn = (j == nthreads - 1) ? columns : firstColumn + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx;
          int idxOther;
          for (int c = firstColumn; c < lastColumn; c++) {
            idxOther = zeroOther + c * columnStrideOther;
            idx = _zero + (c * rows) * _stride;
            for (int r = 0; r < rows; r++) {
              elemsOther[idxOther] = _elements[idx];
              elemsOther[idxOther + 1] = _elements[idx + 1];
              idxOther += rowStrideOther;
              idx += _stride;
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {
      int idxOther;
      int idx = _zero;
      for (int c = 0; c < columns; c++) {
        idxOther = zeroOther + c * columnStrideOther;
        for (int r = 0; r < rows; r++) {
          elemsOther[idxOther] = _elements[idx];
          elemsOther[idxOther + 1] = _elements[idx + 1];
          idxOther += rowStrideOther;
          idx += _stride;
        }
      }
    }
    return M;
  }

  /*DComplexMatrix3D reshape3D(final int slices, final int rows, final int columns) {
    if (slices * rows * columns != _size) {
      throw new IllegalArgumentException("slices*rows*columns != size");
    }
    DComplexMatrix3D M = new DenseDComplexMatrix3D(slices, rows, columns);
    final Float64List elemsOther = M.elements() as Float64List;
    final int zeroOther = M.index(0, 0, 0);
    final int sliceStrideOther = M.sliceStride();
    final int rowStrideOther = M.rowStride();
    final int columnStrideOther = M.columnStride();
    int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = slices ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstSlice = j * k;
        final int lastSlice = (j == nthreads - 1) ? slices : firstSlice + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx;
          int idxOther;
          for (int s = firstSlice; s < lastSlice; s++) {
            for (int c = 0; c < columns; c++) {
              idxOther = zeroOther + s * sliceStrideOther + c * columnStrideOther;
              idx = _zero + (s * rows * columns + c * rows) * _stride;
              for (int r = 0; r < rows; r++) {
                elemsOther[idxOther] = _elements[idx];
                elemsOther[idxOther + 1] = _elements[idx + 1];
                idxOther += rowStrideOther;
                idx += _stride;
              }
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {
      int idxOther;
      int idx = _zero;
      for (int s = 0; s < slices; s++) {
        for (int c = 0; c < columns; c++) {
          idxOther = zeroOther + s * sliceStrideOther + c * columnStrideOther;
          for (int r = 0; r < rows; r++) {
            elemsOther[idxOther] = _elements[idx];
            elemsOther[idxOther + 1] = _elements[idx + 1];
            idxOther += rowStrideOther;
            idx += _stride;
          }
        }
      }
    }
    return M;
  }*/

  void setPartsQuick(int index, double re, double im) {
    int idx = _zero + index * _stride;
    this._elements[idx] = re;
    this._elements[idx + 1] = im;
  }

  void setQuick(int index, Float64List value) {
    int idx = _zero + index * _stride;
    this._elements[idx] = value[0];
    this._elements[idx + 1] = value[1];
  }

  void swap(DComplexMatrix1D other) {
    if (!(other is DenseDComplexMatrix1D)) {
      super.swap(other);
    }
    DenseDComplexMatrix1D y = other as DenseDComplexMatrix1D;
    if (y == this) {
      return;
    }
    checkSize(y);

    final Float64List elemsOther = y._elements;
    if (_elements == null || elemsOther == null) {
      throw new Error();
    }
    final int strideOther = y._stride;
    final int zeroOther = y.index(0);

    int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = _zero + firstIdx * _stride;
          int idxOther = zeroOther + firstIdx * strideOther;
          double tmp;
          for (int k = firstIdx; k < lastIdx; k++) {
            tmp = _elements[idx];
            _elements[idx] = elemsOther[idxOther];
            elemsOther[idxOther] = tmp;
            tmp = _elements[idx + 1];
            _elements[idx + 1] = elemsOther[idxOther + 1];
            elemsOther[idxOther + 1] = tmp;
            idx += _stride;
            idxOther += strideOther;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {
      int idx = _zero;
      int idxOther = zeroOther;
      double tmp;
      for (int k = 0; k < _size; k++) {
        tmp = _elements[idx];
        _elements[idx] = elemsOther[idxOther];
        elemsOther[idxOther] = tmp;
        tmp = _elements[idx + 1];
        _elements[idx + 1] = elemsOther[idxOther + 1];
        elemsOther[idxOther + 1] = tmp;
        idx += _stride;
        idxOther += strideOther;
      }
    }
  }

  void toArrayFill(Float64List values) {
    if (values.length < 2 * _size) {
      throw new ArgumentError("values too small");
    }
    if (_isNoView) {
      //System.arraycopy(this._elements, 0, values, 0, this._elements.length);
      values.setAll(0, this._elements);
    } else {
      super.toArrayFill(values);
    }
  }

  Float64List zDotProduct(final DComplexMatrix1D y, [final int from=0, int length=null]) {
    if (length == null) {
      length = this.size();
    }
    int size = this.size();
    if (from < 0 || length <= 0) {
      return new Float64List.fromList([0, 0]);
    }

    int tail = from + length;
    if (size < tail) {
      tail = size;
    }
    if (y.size() < tail) {
      tail = y.size();
    }
    length = tail - from;
    final Float64List elemsOther = y.elements() as Float64List;
    if (_elements == null || elemsOther == null) {
      throw new Error();
    }
    final int strideOther = y.stride();
    final int zero = index(from);
    final int zeroOther = y.index(from);
    Float64List sum = new Float64List(2);

    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (length >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, length);
      List<Future> futures = new List<Future>(nthreads);
      List<Float64List> results = new List<Float64List>(nthreads);//[2];
      int k = length ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? length : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          Float64List sum = new Float64List(2);
          int idx = zero + firstIdx * _stride;
          int idxOther = zeroOther + firstIdx * strideOther;
          for (int k = firstIdx; k < lastIdx; k++) {
            sum[0] += _elements[idx] * elemsOther[idxOther] + _elements[idx + 1] * elemsOther[idxOther + 1];
            sum[1] += _elements[idx + 1] * elemsOther[idxOther] - _elements[idx] * elemsOther[idxOther + 1];
            idx += _stride;
            idxOther += strideOther;
          }
          return sum;
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          results[j] = futures[j].get() as Float64List;
        }
        sum = results[0];
        for (int j = 1; j < nthreads; j++) {
          sum = DComplex.plus(sum, results[j]);
        }
      } on ExecutionException catch (ex) {
        ex.printStackTrace();
      } on InterruptedException catch (e) {
        e.printStackTrace();
      }
    } else {*/
      int idx = zero;
      int idxOther = zeroOther;
      for (int k = 0; k < length; k++) {
        sum[0] += _elements[idx] * elemsOther[idxOther] + _elements[idx + 1] * elemsOther[idxOther + 1];
        sum[1] += _elements[idx + 1] * elemsOther[idxOther] - _elements[idx] * elemsOther[idxOther + 1];
        idx += _stride;
        idxOther += strideOther;
      }
    //}
    return sum;
  }

  Float64List zSum() {
    Float64List sum = new Float64List(2);
    if (this._elements == null) {
      throw new Error();
    }
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      List<Float64List> results = new List<Float64List>(nthreads);//[2];
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          Float64List sum = new Float64List(2);
          int idx = _zero + firstIdx * _stride;
          for (int k = firstIdx; k < lastIdx; k++) {
            sum[0] += _elements[idx];
            sum[1] += _elements[idx + 1];
            idx += _stride;
          }
          return sum;
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          results[j] = futures[j].get() as Float64List;
        }
        sum = results[0];
        for (int j = 1; j < nthreads; j++) {
          sum[0] = sum[0] + results[j][0];
          sum[1] = sum[1] + results[j][1];
        }
      } on ExecutionException catch (ex) {
        ex.printStackTrace();
      } on InterruptedException catch (e) {
        e.printStackTrace();
      }
    } else {*/
      int idx = _zero;
      for (int k = 0; k < _size; k++) {
        sum[0] += _elements[idx];
        sum[1] += _elements[idx + 1];
        idx += _stride;
      }
    //}
    return sum;
  }

  int _cardinality(int maxCardinality) {
    int cardinality = 0;
    int idx = _zero;
    int i = 0;
    while (i < _size && cardinality < maxCardinality) {
      if (_elements[idx] != 0 || _elements[idx + 1] != 0) {
        cardinality++;
      }
      idx += _stride;
      i++;
    }
    return cardinality;
  }

  bool _haveSharedCellsRaw(DComplexMatrix1D other) {
    if (other is SelectedDenseDComplexMatrix1D) {
      return this._elements == other._elements;
    } else if (other is DenseDComplexMatrix1D) {
      return this._elements == other._elements;
    }
    return false;
  }

  int index(int rank) {
    return _zero + rank * _stride;
  }

  DComplexMatrix1D _viewSelectionLike(Int32List offsets) {
    return new SelectedDenseDComplexMatrix1D(this._elements, offsets);
  }
}
