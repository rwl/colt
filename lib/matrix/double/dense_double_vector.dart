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
 * Dense 1-d matrix (aka <i>vector</i>) holding <tt>double</tt> elements. First
 * see the <a href="package-summary.html">package summary</a> and javadoc <a
 * href="package-tree.html">tree view</a> to get the broad picture.
 * <p>
 * <b>Implementation:</b>
 * <p>
 * Internally holds one single contigous one-dimensional array. Note that this
 * implementation is not synchronized.
 * <p>
 * <b>Time complexity:</b>
 * <p>
 * <tt>O(1)</tt> (i.e. constant time) for the basic operations <tt>get</tt>,
 * <tt>getQuick</tt>, <tt>set</tt>, <tt>setQuick</tt> and <tt>size</tt>,
 * <p>
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 09/24/99
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
class DenseDoubleVector extends DoubleVector {

  /*DoubleFFT_1D _fft;

    DoubleDCT_1D _dct;

    DoubleDST_1D _dst;

    DoubleDHT_1D _dht;*/

  /**
   * The elements of this matrix.
   */
  Float64List _elements;

  /**
   * Constructs a matrix with a copy of the given values. The values are
   * copied. So subsequent changes in <tt>values</tt> are not reflected in the
   * matrix, and vice-versa.
   *
   * @param values
   *            The values to be filled into the new matrix.
   */
  factory DenseDoubleVector.fromList(Float64List values) {
    final m = new DenseDoubleVector(values.length);
    m.setValues(values);
    return m;
  }

  /**
   * Constructs a matrix with a given number of cells. All entries are
   * initially <tt>0</tt>.
   *
   * @param size
   *            the number of cells the matrix shall have.
   * @throws ArgumentError
   *             if <tt>size<0</tt>.
   */
  //    DenseDoubleVector(int size) {
  //        _setUp(size);
  //        this._elements = new Float64List(size);
  //    }

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
   * @param isView
   *            if true then a matrix view is constructed
   * @throws ArgumentError
   *             if <tt>size<0</tt>.
   */
  DenseDoubleVector(int size, [Float64List elements = null, int zero = 0, int stride = 1, bool isView = false]) {
    if (elements == null) {
      elements = new Float64List(size);
    }
    _setUp(size, zero, stride);
    this._elements = elements;
    this._isNoView = !isView;
  }

  double reduce(final DoubleDoubleFunction aggr, final DoubleFunction f) {
    if (_size == 0) return double.NAN;
    double a = 0.0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = _size - j * k;
        final int lastIdx = (j == (nthreads - 1)) ? 0 : firstIdx - k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = _zero + (firstIdx - 1) * _stride;
          double a = f(_elements[idx]);
          for (int i = firstIdx - 1; --i >= lastIdx; ) {
            a = aggr(a, f(_elements[idx -= _stride]));
          }
          return a;
        });
      }
      a = ConcurrencyUtils.waitForCompletion(futures, aggr);
    } else {*/
    int idx = _zero + (_size - 1) * _stride;
    a = f(_elements[idx]);
    for (int i = _size - 1; --i >= 0; ) {
      a = aggr(a, f(_elements[idx -= _stride]));
    }
    //}
    return a;
  }

  double reduceRange(final DoubleDoubleFunction aggr, final DoubleFunction f, final /*IntArrayList*/List<int> indexList) {
    if (this.length == 0) return double.NAN;
    final int size = indexList.length;
    final List<int> indexElements = indexList;//.elements();
    double a = 0.0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, size);
      List<Future> futures = new List<Future>(nthreads);
      int k = size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = _zero + indexElements[firstIdx] * _stride;
          double a = f(_elements[idx]);
          double elem;
          for (int i = firstIdx + 1; i < lastIdx; i++) {
            idx = _zero + indexElements[i] * _stride;
            elem = _elements[idx];
            a = aggr(a, f(elem));
          }
          return a;
        });
      }
      a = ConcurrencyUtils.waitForCompletion(futures, aggr);
    } else {*/
    double elem;
    int idx = _zero + indexElements[0] * _stride;
    a = f(_elements[idx]);
    for (int i = 1; i < size; i++) {
      idx = _zero + indexElements[i] * _stride;
      elem = _elements[idx];
      a = aggr(a, f(elem));
    }
    //}
    return a;
  }

  double reduceVector(final DoubleVector other, final DoubleDoubleFunction aggr, final DoubleDoubleFunction f) {
    if (!(other is DenseDoubleVector)) {
      return super.reduceVector(other, aggr, f);
    }
    checkSize(other);
    if (_size == 0) return double.NAN;
    final int zeroOther = other.index(0);
    final int strideOther = other.stride();
    final Float64List elementsOther = other.elements() as Float64List;
    double a = 0.0;
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
          double a = f(_elements[idx], elementsOther[idxOther]);
          for (int i = firstIdx + 1; i < lastIdx; i++) {
            idx += _stride;
            idxOther += strideOther;
            a = aggr(a, f(_elements[idx], elementsOther[idxOther]));
          }
          return a;
        });
      }
      a = ConcurrencyUtils.waitForCompletion(futures, aggr);
    } else {*/
    a = f(_elements[_zero], elementsOther[zeroOther]);
    int idx = _zero;
    int idxOther = zeroOther;
    for (int i = 1; i < _size; i++) {
      idx += _stride;
      idxOther += strideOther;
      a = aggr(a, f(_elements[idx], elementsOther[idxOther]));
    }
    //}
    return a;
  }

  DoubleVector forEach(final DoubleFunction function) {
    double multiplicator;
    if (function is DoubleMult) {
      // x[i] = mult*x[i]
      multiplicator = (function as DoubleMult).multiplicator;
      if (multiplicator == 1) {
        return this;
      }
    } else {
      multiplicator = 0.0;
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
          // specialization for speed
          if (function is DoubleMult) {
            // x[i] = mult*x[i]
            for (int k = firstIdx; k < lastIdx; k++) {
              _elements[idx] *= multiplicator;
              idx += _stride;
            }
          } else {
            // the general case x[i] = f(x[i])
            for (int k = firstIdx; k < lastIdx; k++) {
              _elements[idx] = function(_elements[idx]);
              idx += _stride;
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idx = _zero - _stride;
    // specialization for speed
    if (function is DoubleMult) {
      // x[i] = mult*x[i]
      for (int k = _size; --k >= 0; ) {
        _elements[idx += _stride] *= multiplicator;
      }
    } else {
      // the general case x[i] = f(x[i])
      for (int k = _size; --k >= 0; ) {
        _elements[idx += _stride] = function(_elements[idx]);
      }
    }
    //}
    return this;
  }

  DoubleVector forEachWhere(final DoubleProcedure cond, final DoubleFunction function) {
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
          for (int i = firstIdx; i < lastIdx; i++) {
            if (cond(_elements[idx]) == true) {
              _elements[idx] = function(_elements[idx]);
            }
            idx += _stride;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idx = _zero;
    for (int i = 0; i < _size; i++) {
      if (cond(_elements[idx]) == true) {
        _elements[idx] = function(_elements[idx]);
      }
      idx += _stride;
    }
    //}
    return this;
  }

  DoubleVector fillWhere(final DoubleProcedure cond, final double value) {
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
          for (int i = firstIdx; i < lastIdx; i++) {
            if (cond(_elements[idx]) == true) {
              _elements[idx] = value;
            }
            idx += _stride;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idx = _zero;
    for (int i = 0; i < _size; i++) {
      if (cond(_elements[idx]) == true) {
        _elements[idx] = value;
      }
      idx += _stride;
    }
    //}
    return this;
  }

  DoubleVector fill(final double value) {
    final Float64List elems = this._elements;
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
          for (int k = firstIdx; k < lastIdx; k++) {
            elems[idx] = value;
            idx += _stride;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      int idx = _zero;
      for (int i = 0; i < _size; i++) {
        elems[idx] = value;
        idx += _stride;
      }
    //}
    return this;
  }

  DoubleVector setValues(final Float64List values) {
    if (values.length != _size) {
      throw new ArgumentError("Must have same number of cells: length=${values.length} size()=${length}");
    }
    //int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if (_isNoView) {
      //System.arraycopy(values, 0, this._elements, 0, values.length);
      this._elements.setAll(0, values);
    } else {
      /*if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
        nthreads = Math.min(nthreads, _size);
        List<Future> futures = new List<Future>(nthreads);
        int k = _size ~/ nthreads;
        for (int j = 0; j < nthreads; j++) {
          final int firstIdx = j * k;
          final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
          futures[j] = ConcurrencyUtils.submit(() {
            int idx = _zero + firstIdx * _stride;
            for (int i = firstIdx; i < lastIdx; i++) {
              _elements[idx] = values[i];
              idx += _stride;
            }
          });
        }
        ConcurrencyUtils.waitForCompletion(futures);
      } else {*/
      int idx = _zero;
      for (int i = 0; i < _size; i++) {
        _elements[idx] = values[i];
        idx += _stride;
      }
      //}
    }
    return this;
  }

  DoubleVector setAll(DoubleVector source) {
    // overriden for performance only
    if (!(source is DenseDoubleVector)) {
      super.setAll(source);
      return this;
    }
    DenseDoubleVector other = source as DenseDoubleVector;
    if (other == this) return this;
    checkSize(other);
    if (_isNoView && other._isNoView) {
      // quickest
      this._elements.setAll(0, other._elements);
      //System.arraycopy(other._elements, 0, this._elements, 0, this._elements.length);
      return this;
    }
    if (_haveSharedCells(other)) {
      DoubleVector c = other.copy();
      if (!(c is DenseDoubleVector)) {
        // should not happen
        super.setAll(source);
        return this;
      }
      other = c as DenseDoubleVector;
    }

    final Float64List elementsOther = other._elements;
    if (_elements == null || elementsOther == null) {
      throw new Error();
    }
    final int zeroOther = other.index(0);
    final int strideOther = other._stride;
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
          for (int k = firstIdx; k < lastIdx; k++) {
            _elements[idx] = elementsOther[idxOther];
            idx += _stride;
            idxOther += strideOther;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idx = _zero;
    int idxOther = zeroOther;
    for (int k = 0; k < _size; k++) {
      _elements[idx] = elementsOther[idxOther];
      idx += _stride;
      idxOther += strideOther;
    }
    //}
    return this;
  }

  DoubleVector forEachVector(final DoubleVector y, final DoubleDoubleFunction function) {
    // overriden for performance only
    if (!(y is DenseDoubleVector)) {
      super.forEachVector(y, function);
      return this;
    }
    checkSize(y);
    final int zeroOther = y.index(0);
    final int strideOther = y.stride();
    final Float64List elementsOther = y.elements() as Float64List;
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
          // specialized for speed
          if (function == func.mult) {
            // x[i] = x[i] * y[i]
            for (int k = firstIdx; k < lastIdx; k++) {
              _elements[idx] *= elementsOther[idxOther];
              idx += _stride;
              idxOther += strideOther;
            }
          } else if (function == func.div) {
            // x[i] = x[i] / y[i]
            for (int k = firstIdx; k < lastIdx; k++) {
              _elements[idx] /= elementsOther[idxOther];
              idx += _stride;
              idxOther += strideOther;

            }
          } else if (function is DoublePlusMultFirst) {
            double multiplicator = function.multiplicator;
            if (multiplicator == 0) {
              // x[i] = 0*x[i] + y[i]
              for (int k = firstIdx; k < lastIdx; k++) {
                _elements[idx] = elementsOther[idxOther];
                idx += _stride;
                idxOther += strideOther;
              }
            } else if (multiplicator == 1) {
              // x[i] = x[i] + y[i]
              for (int k = firstIdx; k < lastIdx; k++) {
                _elements[idx] += elementsOther[idxOther];
                idx += _stride;
                idxOther += strideOther;
              }
            } else if (multiplicator == -1) {
              // x[i] = -x[i] + y[i]
              for (int k = firstIdx; k < lastIdx; k++) {
                _elements[idx] = elementsOther[idxOther] - _elements[idx];
                idx += _stride;
                idxOther += strideOther;
              }
            } else {
              // the general case x[i] = mult*x[i] + y[i]
              for (int k = firstIdx; k < lastIdx; k++) {
                _elements[idx] = multiplicator * _elements[idx] + elementsOther[idxOther];
                idx += _stride;
                idxOther += strideOther;
              }
            }
          } else if (function is DoublePlusMultSecond) {
            double multiplicator = function.multiplicator;
            if (multiplicator == 0) {
              // x[i] = x[i] + 0*y[i]
              return;
            } else if (multiplicator == 1) {
              // x[i] = x[i] + y[i]
              for (int k = firstIdx; k < lastIdx; k++) {
                _elements[idx] += elementsOther[idxOther];
                idx += _stride;
                idxOther += strideOther;
              }
            } else if (multiplicator == -1) {
              // x[i] = x[i] - y[i]
              for (int k = firstIdx; k < lastIdx; k++) {
                _elements[idx] -= elementsOther[idxOther];
                idx += _stride;
                idxOther += strideOther;
              }
            } else {
              // the general case x[i] = x[i] + mult*y[i]
              for (int k = firstIdx; k < lastIdx; k++) {
                _elements[idx] += multiplicator * elementsOther[idxOther];
                idx += _stride;
                idxOther += strideOther;
              }

            }
          } else {
            // the general case x[i] = f(x[i],y[i])
            for (int k = firstIdx; k < lastIdx; k++) {
              _elements[idx] = function(_elements[idx], elementsOther[idxOther]);
              idx += _stride;
              idxOther += strideOther;
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    // specialized for speed
    int idx = _zero;
    int idxOther = zeroOther;
    if (function == func.mult) {
      // x[i] = x[i] * y[i]
      for (int k = 0; k < _size; k++) {
        _elements[idx] *= elementsOther[idxOther];
        idx += _stride;
        idxOther += strideOther;
      }
    } else if (function == func.div) {
      // x[i] = x[i] / y[i]
      for (int k = 0; k < _size; k++) {
        _elements[idx] /= elementsOther[idxOther];
        idx += _stride;
        idxOther += strideOther;
      }
    } else if (function is DoublePlusMultSecond) {
      double multiplicator = (function as DoublePlusMultSecond).multiplicator;
      if (multiplicator == 0) {
        // x[i] = x[i] + 0*y[i]
        return this;
      } else if (multiplicator == 1) {
        // x[i] = x[i] + y[i]
        for (int k = 0; k < _size; k++) {
          _elements[idx] += elementsOther[idxOther];
          idx += _stride;
          idxOther += strideOther;
        }
      } else if (multiplicator == -1) {
        // x[i] = x[i] - y[i]
        for (int k = 0; k < _size; k++) {
          _elements[idx] -= elementsOther[idxOther];
          idx += _stride;
          idxOther += strideOther;
        }
      } else {
        // the general case x[i] = x[i] + mult*y[i]
        for (int k = 0; k < _size; k++) {
          _elements[idx] += multiplicator * elementsOther[idxOther];
          idx += _stride;
          idxOther += strideOther;
        }
      }
    } else {
      // the general case x[i] = f(x[i],y[i])
      for (int k = 0; k < _size; k++) {
        _elements[idx] = function(_elements[idx], elementsOther[idxOther]);
        idx += _stride;
        idxOther += strideOther;
      }
    }
    //}
    return this;
  }

  int get cardinality {
    int cardinality = 0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      List<int> results = new List<int>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int cardinality = 0;
          int idx = _zero + firstIdx * _stride;
          for (int i = firstIdx; i < lastIdx; i++) {
            if (_elements[idx] != 0) cardinality++;
            idx += _stride;
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
    int idx = _zero;
    for (int i = 0; i < _size; i++) {
      if (_elements[idx] != 0) cardinality++;
      idx += _stride;
    }
    //}
    return cardinality;
  }

  Float64List elements() {
    return _elements;
  }

  void nonZeros(final /*IntArrayList*/List<int> indexList, final /*DoubleArrayList*/List<double> valueList) {
    bool fillIndexList = indexList != null;
    bool fillValueList = valueList != null;
    if (fillIndexList) indexList.clear();
    if (fillValueList) valueList.clear();
    int rem = _size % 2;
    int idx = _zero;
    if (rem == 1) {
      double value = _elements[idx];
      if (value != 0) {
        if (fillIndexList) {
          indexList.add(0);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      idx += _stride;

    }
    for (int i = rem; i < _size; i += 2) {
      double value = _elements[idx];
      if (value != 0) {
        if (fillIndexList) {
          indexList.add(i);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      idx += _stride;
      value = _elements[idx];
      if (value != 0) {
        if (fillIndexList) {
          indexList.add(i + 1);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      idx += _stride;
    }
  }

  void positiveValues(final /*IntArrayList*/List<int> indexList, final /*DoubleArrayList*/List<double> valueList) {
    bool fillIndexList = indexList != null;
    bool fillValueList = valueList != null;
    if (fillIndexList) indexList.clear();
    if (fillValueList) valueList.clear();
    int rem = _size % 2;
    int idx = _zero;
    if (rem == 1) {
      double value = _elements[idx];
      if (value > 0) {
        if (fillIndexList) {
          indexList.add(0);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      idx += _stride;

    }
    for (int i = rem; i < _size; i += 2) {
      double value = _elements[idx];
      if (value > 0) {
        if (fillIndexList) {
          indexList.add(i);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      idx += _stride;
      value = _elements[idx];
      if (value > 0) {
        if (fillIndexList) {
          indexList.add(i + 1);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      idx += _stride;
    }
  }

  void negativeValues(final /*IntArrayList*/List<int> indexList, final /*DoubleArrayList*/List<double> valueList) {
    bool fillIndexList = indexList != null;
    bool fillValueList = valueList != null;
    if (fillIndexList) indexList.clear();
    if (fillValueList) valueList.clear();
    int rem = _size % 2;
    int idx = _zero;
    if (rem == 1) {
      double value = _elements[idx];
      if (value < 0) {
        if (fillIndexList) {
          indexList.add(0);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      idx += _stride;

    }
    for (int i = rem; i < _size; i += 2) {
      double value = _elements[idx];
      if (value < 0) {
        if (fillIndexList) {
          indexList.add(i);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      idx += _stride;
      value = _elements[idx];
      if (value < 0) {
        if (fillIndexList) {
          indexList.add(i + 1);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      idx += _stride;
    }
  }

  DoubleVectorLocation max() {
    int location = 0;
    double maxValue = 0.0;
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
          int idx = _zero + firstIdx * _stride;
          double maxValue = _elements[idx];
          int location = (idx - _zero) ~/ _stride;
          for (int i = firstIdx + 1; i < lastIdx; i++) {
            idx += _stride;
            if (maxValue < _elements[idx]) {
              maxValue = _elements[idx];
              location = (idx - _zero) ~/ _stride;
            }
          }
          return [maxValue, location];
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
    maxValue = _elements[_zero];
    location = 0;
    int idx = _zero;
    for (int i = 1; i < _size; i++) {
      idx += _stride;
      if (maxValue < _elements[idx]) {
        maxValue = _elements[idx];
        location = (idx - _zero) ~/ _stride;
      }
    }
    //}
    return new DoubleVectorLocation(maxValue, location);
  }

  DoubleVectorLocation min() {
    int location = 0;
    double minValue = 0.0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      List<List<double>> results = new List<Float64List>(nthreads);//[2];
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = _zero + firstIdx * _stride;
          double minValue = _elements[idx];
          int location = (idx - _zero) ~/ _stride;
          for (int i = firstIdx + 1; i < lastIdx; i++) {
            idx += _stride;
            if (minValue > _elements[idx]) {
              minValue = _elements[idx];
              location = (idx - _zero) ~/ _stride;
            }
          }
          return [minValue, location];
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
    minValue = _elements[_zero];
    location = 0;
    int idx = _zero;
    for (int i = 1; i < _size; i++) {
      idx += _stride;
      if (minValue > _elements[idx]) {
        minValue = _elements[idx];
        location = (idx - _zero) ~/ _stride;
      }
    }
    //}
    return new DoubleVectorLocation(minValue, location);
  }

  double get(int index) {
    return _elements[_zero + index * _stride];
  }

  DoubleVector like1D(int size) {
    return new DenseDoubleVector(size);
  }

  DoubleMatrix like2D(int rows, int columns) {
    return new DenseDoubleMatrix(rows, columns);
  }

  DoubleMatrix reshape(final int rows, final int columns) {
    if (rows * columns != _size) {
      throw new ArgumentError("rows*columns != size");
    }
    DoubleMatrix M = new DenseDoubleMatrix(rows, columns);
    final Float64List elementsOther = M.elements() as Float64List;
    final int zeroOther = M.index(0, 0);
    final int rowStrideOther = M.rowStride;
    final int columnStrideOther = M.columnStride;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
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
              elementsOther[idxOther] = _elements[idx];
              idxOther += rowStrideOther;
              idx += _stride;
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      int idx = _zero;
      for (int c = 0; c < columns; c++) {
        int idxOther = zeroOther + c * columnStrideOther;
        for (int r = 0; r < rows; r++) {
          elementsOther[idxOther] = _elements[idx];
          idxOther += rowStrideOther;
          idx += _stride;
        }
      }
    //}
    return M;
  }

  /*DoubleMatrix3D reshape3D(final int slices, final int rows, final int columns) {
    if (slices * rows * columns != _size) {
      throw new ArgumentError("slices*rows*columns != size");
    }
    DoubleMatrix3D M = new DenseDoubleMatrix3D(slices, rows, columns);
    final Float64List elementsOther = M.elements() as Float64List;
    final int zeroOther = M.index(0, 0, 0) as int;
    final int sliceStrideOther = M.sliceStride();
    final int rowStrideOther = M.rowStride();
    final int columnStrideOther = M.columnStride();
    int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = slices / nthreads;
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
                elementsOther[idxOther] = _elements[idx];
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
            elementsOther[idxOther] = _elements[idx];
            idxOther += rowStrideOther;
            idx += _stride;
          }
        }
      }
    }
    return M;
  }*/

  void set(int index, double value) {
    _elements[_zero + index * _stride] = value;
  }

  void swap(final DoubleVector other) {
    // overriden for performance only
    if (!(other is DenseDoubleVector)) {
      super.swap(other);
    }
    DenseDoubleVector y = other as DenseDoubleVector;
    if (y == this) return;
    checkSize(y);
    final Float64List elementsOther = y._elements;
    if (_elements == null || elementsOther == null) throw new Error();
    final int zeroOther = other.index(0);
    final int strideOther = other.stride();
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
          for (int k = firstIdx; k < lastIdx; k++) {
            double tmp = _elements[idx];
            _elements[idx] = elementsOther[idxOther];
            elementsOther[idxOther] = tmp;
            idx += _stride;
            idxOther += strideOther;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idx = _zero;
    int idxOther = zeroOther;
    for (int k = 0; k < _size; k++) {
      double tmp = _elements[idx];
      _elements[idx] = elementsOther[idxOther];
      elementsOther[idxOther] = tmp;
      idx += _stride;
      idxOther += strideOther;
    }
    //}
  }

  void fillList(Float64List values) {
    if (values.length < _size) throw new ArgumentError("values too small");
    if (this._isNoView) {
      values.setAll(0, this._elements);
      //System.arraycopy(this._elements, 0, values, 0, this._elements.length);
    } else {
      super.fillList(values);
    }
  }

  double dot(DoubleVector y, [final int from = 0, int length = null]) {
    if (length == null) {
      length = _size;
    }
    if (!(y is DenseDoubleVector)) {
      return super.dot(y, from, length);
    }
    DenseDoubleVector yy = y as DenseDoubleVector;

    int tail = from + length;
    if (from < 0 || length < 0) return 0.0;
    if (_size < tail) tail = _size;
    if (y.length < tail) tail = y.length;
    final Float64List elementsOther = yy._elements;
    int zeroThis = index(from);
    int zeroOther = yy.index(from);
    int strideOther = yy._stride;
    if (_elements == null || elementsOther == null) throw new Error();
    double sum = 0.0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (length >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      final int zeroThisF = zeroThis;
      final int zeroOtherF = zeroOther;
      final int strideOtherF = strideOther;
      nthreads = Math.min(nthreads, length);
      List<Future> futures = new List<Future>(nthreads);
      Float64List results = new Float64List(nthreads);
      int k = length ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? length : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = zeroThisF + firstIdx * _stride;
          int idxOther = zeroOtherF + firstIdx * strideOtherF;
          idx -= _stride;
          idxOther -= strideOtherF;
          double sum = 0.0;
          int min = lastIdx - firstIdx;
          for (int k = min ~/ 4; --k >= 0; ) {
            sum += _elements[idx += _stride] * elementsOther[idxOther += strideOtherF] + _elements[idx += _stride] * elementsOther[idxOther += strideOtherF] + _elements[idx += _stride] * elementsOther[idxOther += strideOtherF] + _elements[idx += _stride] * elementsOther[idxOther += strideOtherF];
          }
          for (int k = min % 4; --k >= 0; ) {
            sum += _elements[idx += _stride] * elementsOther[idxOther += strideOtherF];
          }
          return sum;
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          results[j] = futures[j].get() as double;
        }
        sum = results[0];
        for (int j = 1; j < nthreads; j++) {
          sum += results[j];
        }
      } on ExecutionException catch (ex) {
        ex.printStackTrace();
      } on InterruptedException catch (e) {
        e.printStackTrace();
      }
    } else {*/
    zeroThis -= _stride;
    zeroOther -= strideOther;
    int min = tail - from;
    for (int k = min ~/ 4; --k >= 0; ) {
      sum += _elements[zeroThis += _stride] * elementsOther[zeroOther += strideOther] + _elements[zeroThis += _stride] * elementsOther[zeroOther += strideOther] + _elements[zeroThis += _stride] * elementsOther[zeroOther += strideOther] + _elements[zeroThis += _stride] * elementsOther[zeroOther += strideOther];
    }
    for (int k = min % 4; --k >= 0; ) {
      sum += _elements[zeroThis += _stride] * elementsOther[zeroOther += strideOther];
    }
    //}
    return sum;
  }

  double sum() {
    double sum = 0.0;
    final Float64List elems = this._elements;
    if (elems == null) throw new Error();
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      Float64List results = new Float64List(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          double sum = 0.0;
          int idx = _zero + firstIdx * _stride;
          for (int i = firstIdx; i < lastIdx; i++) {
            sum += elems[idx];
            idx += _stride;
          }
          return sum;
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          results[j] = futures[j].get() as double;
        }
        sum = results[0];
        for (int j = 1; j < nthreads; j++) {
          sum += results[j];
        }
      } on ExecutionException catch (ex) {
        ex.printStackTrace();
      } on InterruptedException catch (e) {
        e.printStackTrace();
      }
    } else {*/
    int idx = _zero;
    for (int k = 0; k < _size; k++) {
      sum += elems[idx];
      idx += _stride;
    }
    //}
    return sum;
  }

  int _cardinality(int maxCardinality) {
    int cardinality = 0;
    int index = _zero;
    Float64List elems = this._elements;
    int i = _size;
    while (--i >= 0 && cardinality < maxCardinality) {
      if (elems[index] != 0) cardinality++;
      index += _stride;
    }
    return cardinality;
  }

  bool _haveSharedCellsRaw(DoubleVector other) {
    if (other is SelectedDenseDoubleVector) {
      //SelectedDenseDoubleVector otherMatrix = other as SelectedDenseDoubleVector;
      return this._elements == other._elements;
    } else if (other is DenseDoubleVector) {
      return this._elements == other._elements;
    }
    return false;
  }

  int index(int rank) {
    return _zero + rank * _stride;
  }

  DoubleVector _viewSelectionLike(List<int> offsets) {
    return new SelectedDenseDoubleVector.offset(this._elements, offsets);
  }

  Object clone() {
    return new DenseDoubleVector(_size, _elements, _zero, _stride, !_isNoView);
  }
}

/**
 * Selection view on dense 1-d matrices holding <tt>double</tt> elements. First
 * see the <a href="package-summary.html">package summary</a> and javadoc <a
 * href="package-tree.html">tree view</a> to get the broad picture.
 * <p>
 * <b>Implementation:</b>
 * <p>
 * Objects of this class are typically constructed via <tt>viewIndexes</tt>
 * methods on some source matrix. The interface introduced in abstract super
 * classes defines everything a user can do. From a user point of view there is
 * nothing special about this class; it presents the same functionality with the
 * same signatures and semantics as its abstract superclass(es) while
 * introducing no additional functionality. Thus, this class need not be visible
 * to users. By the way, the same principle applies to concrete DenseXXX,
 * SparseXXX classes: they presents the same functionality with the same
 * signatures and semantics as abstract superclass(es) while introducing no
 * additional functionality. Thus, they need not be visible to users, either.
 * Factory methods could hide all these concrete types.
 * <p>
 * This class uses no delegation. Its instances point directly to the data. Cell
 * addressing overhead is 1 additional array index access per get/set.
 * <p>
 * Note that this implementation is not synchronized.
 * <p>
 * <b>Memory requirements:</b>
 * <p>
 * <tt>memory [bytes] = 4*indexes.length</tt>. Thus, an index view with 1000
 * indexes additionally uses 4 KB.
 * <p>
 * <b>Time complexity:</b>
 * <p>
 * Depends on the parent view holding cells.
 * <p>
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 09/24/99
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @version 1.1, 08/22/2007
 */
class SelectedDenseDoubleVector extends DoubleVector {

  /**
   * The elements of this matrix.
   */
  Float64List _elements;

  /**
   * The offsets of visible indexes of this matrix.
   */
  Int32List _offsets;

  /**
   * The offset.
   */
  int __offset;

  /**
   * Constructs a matrix view with the given parameters.
   *
   * @param elements
   *            the cells.
   * @param indexes
   *            The indexes of the cells that shall be visible.
   */
  factory SelectedDenseDoubleVector.offset(Float64List elements, Int32List offsets) {
    return new SelectedDenseDoubleVector(offsets.length, elements, 0, 1, offsets, 0);
  }

  /**
   * Constructs a matrix view with the given parameters.
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
   * @param offsets
   *            the offsets of the cells that shall be visible.
   * @param offset
   */
  SelectedDenseDoubleVector(int size, Float64List elements, int zero, int stride, Int32List offsets, int offset) {
    _setUp(size, zero, stride);

    this._elements = elements;
    this._offsets = offsets;
    this.__offset = offset;
    this._isNoView = false;
  }

  Float64List elements() {
    return _elements;
  }

  /**
   * Returns the matrix cell value at coordinate <tt>index</tt>.
   *
   * <p>
   * Provided with invalid parameters this method may return invalid objects
   * without throwing any exception. <b>You should only use this method when
   * you are absolutely sure that the coordinate is within bounds.</b>
   * Precondition (unchecked): <tt>index&lt;0 || index&gt;=size()</tt>.
   *
   * @param index
   *            the index of the cell.
   * @return the value of the specified cell.
   */
  double get(int index) {
    // if (debug) if (index<0 || index>=size) checkIndex(index);
    // return elements[index(index)];
    // manually inlined:
    return _elements[__offset + _offsets[_zero + index * _stride]];
  }

  /**
   * Returns the position of the element with the given relative rank within
   * the (virtual or non-virtual) internal 1-dimensional array. You may want
   * to override this method for performance.
   *
   * @param rank
   *            the rank of the element.
   */
  int index(int rank) {
    // return this.offset + super.index(rank);
    // manually inlined:
    return __offset + _offsets[_zero + rank * _stride];
  }

  /**
   * Construct and returns a new empty matrix <i>of the same dynamic type</i>
   * as the receiver, having the specified size. For example, if the receiver
   * is an instance of type <tt>DenseDoubleVector</tt> the new matrix must
   * also be of type <tt>DenseDoubleVector</tt>, if the receiver is an
   * instance of type <tt>SparseDoubleVector</tt> the new matrix must also
   * be of type <tt>SparseDoubleVector</tt>, etc. In general, the new matrix
   * should have internal parametrization as similar as possible.
   *
   * @param size
   *            the number of cell the matrix shall have.
   * @return a new empty matrix of the same dynamic type.
   */
  DoubleVector like1D(int size) {
    return new DenseDoubleVector(size);
  }

  /**
   * Construct and returns a new 2-d matrix <i>of the corresponding dynamic
   * type</i>, entirelly independent of the receiver. For example, if the
   * receiver is an instance of type <tt>DenseDoubleVector</tt> the new
   * matrix must be of type <tt>DenseDoubleMatrix</tt>, if the receiver is
   * an instance of type <tt>SparseDoubleVector</tt> the new matrix must be
   * of type <tt>SparseDoubleMatrix</tt>, etc.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @return a new matrix of the corresponding dynamic type.
   */
  DoubleMatrix like2D(int rows, int columns) {
    return new DenseDoubleMatrix(rows, columns);
  }

  DoubleMatrix reshape(int rows, int columns) {
    if (rows * columns != _size) {
      throw new ArgumentError("rows*columns != size");
    }
    DoubleMatrix M = new DenseDoubleMatrix(rows, columns);
    final Float64List elementsOther = M.elements() as Float64List;
    final int zeroOther = M.index(0, 0);
    final int rowStrideOther = M.rowStride;
    final int colStrideOther = M.columnStride;
    int idxOther;
    int idx = 0;
    for (int c = 0; c < columns; c++) {
      idxOther = zeroOther + c * colStrideOther;
      for (int r = 0; r < rows; r++) {
        elementsOther[idxOther] = get(idx++);
        idxOther += rowStrideOther;
      }
    }
    return M;
  }

  /*DoubleMatrix3D reshape3D(int slices, int rows, int columns) {
    if (slices * rows * columns != _size) {
      throw new ArgumentError("slices*rows*columns != size");
    }
    DoubleMatrix3D M = new DenseDoubleMatrix3D(slices, rows, columns);
    final Float64List elementsOther = M.elements() as Float64List;
    final int zeroOther = M.index(0, 0, 0);
    final int sliceStrideOther = M.sliceStride();
    final int rowStrideOther = M.rowStride();
    final int colStrideOther = M.columnStride();
    int idxOther;
    int idx = 0;
    for (int s = 0; s < slices; s++) {
      for (int c = 0; c < columns; c++) {
        idxOther = zeroOther + s * sliceStrideOther + c * colStrideOther;
        for (int r = 0; r < rows; r++) {
          elementsOther[idxOther] = getQuick(idx++);
          idxOther += rowStrideOther;
        }
      }
    }
    return M;
  }*/

  /**
   * Sets the matrix cell at coordinate <tt>index</tt> to the specified value.
   *
   * <p>
   * Provided with invalid parameters this method may access illegal indexes
   * without throwing any exception. <b>You should only use this method when
   * you are absolutely sure that the coordinate is within bounds.</b>
   * Precondition (unchecked): <tt>index&lt;0 || index&gt;=size()</tt>.
   *
   * @param index
   *            the index of the cell.
   * @param value
   *            the value to be filled into the specified cell.
   */
  void set(int index, double value) {
    // if (debug) if (index<0 || index>=size) checkIndex(index);
    // elements[index(index)] = value;
    // manually inlined:
    _elements[__offset + _offsets[_zero + index * _stride]] = value;
  }

  /**
   * Returns the position of the given absolute rank within the (virtual or
   * non-virtual) internal 1-dimensional array. Default implementation.
   * Override, if necessary.
   *
   * @param rank
   *            the absolute rank of the element.
   * @return the position.
   */
  int _offset(int absRank) {
    return _offsets[absRank];
  }

  /**
   * Returns <tt>true</tt> if both matrices share at least one identical cell.
   */
  bool _haveSharedCellsRaw(DoubleVector other) {
    if (other is SelectedDenseDoubleVector) {
      return this._elements == other._elements;
    } else if (other is DenseDoubleVector) {
      return this._elements == other._elements;
    }
    return false;
  }

  /**
   * Sets up a matrix with a given number of cells.
   *
   * @param size
   *            the number of cells the matrix shall have.
   */
  void _setUp1D(int size) {
    super._setUp(size);
    this._stride = 1;
    this.__offset = 0;
  }

  /**
   * Construct and returns a new selection view.
   *
   * @param offsets
   *            the offsets of the visible elements.
   * @return a new view.
   */
  DoubleVector _viewSelectionLike(Int32List offsets) {
    return new SelectedDenseDoubleVector.offset(this._elements, offsets);
  }

  Object clone() {
    return new SelectedDenseDoubleVector(_size, _elements, _zero, _stride, _offsets, __offset);
  }
}
