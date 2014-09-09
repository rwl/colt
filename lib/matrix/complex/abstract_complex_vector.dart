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
 * Abstract base class for 1-d matrices (aka <i>vectors</i>) holding
 * <tt>complex</tt> elements. A matrix has a number of cells (its <i>size</i>),
 * which are assigned upon instance construction. Elements are accessed via zero
 * based indexes. Legal indexes are of the form <tt>[0..size()-1]</tt>. Any
 * attempt to access an element at a coordinate
 * <tt>index&lt;0 || index&gt;=size()</tt> will throw an
 * <tt>IndexOutOfBoundsException</tt>.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
abstract class AbstractComplexVector extends AbstractVector {

  /**
   * Makes this class non instantiable, but still let's others inherit from
   * it.
   */
  /*ComplexVector() {
    }*/

  /**
   * Applies a function to each cell and aggregates the results.
   *
   * @param aggr
   *            an aggregation function taking as first argument the current
   *            aggregation and as second argument the transformed current
   *            cell value.
   * @param f
   *            a function transforming the current cell value.
   * @return the aggregated measure.
   * @see cern.jet.math.tdcomplex.ComplexFunctions
   */
  Float64List reduce(final cfunc.ComplexComplexComplexFunction aggr, final cfunc.ComplexComplexFunction f) {
    Float64List b = new Float64List(2);
    int size = this.length;
    if (size == 0) {
      b[0] = double.NAN;
      b[1] = double.NAN;
      return b;
    }
    Float64List a = f(get(0));
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, size);
      List<Future> futures = new List<Future>(nthreads);
      int k = size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          Float64List a = f(get(firstIdx));
          for (int i = firstIdx + 1; i < lastIdx; i++) {
            a = aggr(a, f(get(i)));
          }
          return a;
        });
      }
      a = ConcurrencyUtils.waitForCompletion(futures, aggr);
    } else {*/
    for (int i = 1; i < size; i++) {
      a = aggr(a, f(get(i)));
    }
    //}
    return a;
  }

  /**
   * Applies a function to each corresponding cell of two matrices and
   * aggregates the results.
   *
   * @param other
   *            the secondary matrix to operate on.
   *
   * @param aggr
   *            an aggregation function taking as first argument the current
   *            aggregation and as second argument the transformed current
   *            cell values.
   * @param f
   *            a function transforming the current cell values.
   * @return the aggregated measure.
   * @throws ArgumentError
   *             if <tt>size() != other.size()</tt>.
   * @see cern.jet.math.tdcomplex.ComplexFunctions
   */
  Float64List reduceVector(final AbstractComplexVector other, final cfunc.ComplexComplexComplexFunction aggr, final cfunc.ComplexComplexComplexFunction f) {
    checkSize(other);
    int size = this.length;
    if (size == 0) {
      Float64List b = new Float64List(2);
      b[0] = double.NAN;
      b[1] = double.NAN;
      return b;
    }
    Float64List a = f(get(0), other.get(0));
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, size);
      List<Future> futures = new List<Future>(nthreads);
      int k = size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          Float64List a = f(get(firstIdx), other.get(firstIdx));
          for (int i = firstIdx + 1; i < lastIdx; i++) {
            a = aggr(a, f(get(i), other.get(i)));
          }
          return a;
        });
      }
      a = ConcurrencyUtils.waitForCompletion(futures, aggr);
    } else {*/
    for (int i = 1; i < size; i++) {
      a = aggr(a, f(get(i), other.get(i)));
    }
    //}
    return a;
  }

  /**
   * Assigns the result of a function to each cell;
   *
   * @param f
   *            a function object taking as argument the current cell's value.
   * @return <tt>this</tt> (for convenience only).
   * @see cern.jet.math.tdcomplex.ComplexFunctions
   */
  void forEach(final cfunc.ComplexComplexFunction f) {
    int size = this.length;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, size);
      List<Future> futures = new List<Future>(nthreads);
      int k = size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int i = firstIdx; i < lastIdx; i++) {
            setQuick(i, f(get(i)));
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    for (int i = 0; i < size; i++) {
      set(i, f(get(i)));
    }
    //}
  }

  /**
   * Assigns the result of a function to all cells that satisfy a condition.
   *
   * @param cond
   *            a condition.
   *
   * @param f
   *            a function object.
   * @return <tt>this</tt> (for convenience only).
   * @see cern.jet.math.tdcomplex.ComplexFunctions
   */
  void forEachWhere(final cfunc.ComplexProcedure cond, final cfunc.ComplexComplexFunction f) {
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int i = firstIdx; i < lastIdx; i++) {
            Float64List elem = get(i);
            if (cond(elem) == true) {
              setQuick(i, f(elem));
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    Float64List elem;
    for (int i = 0; i < _size; i++) {
      elem = get(i);
      if (cond(elem) == true) {
        set(i, f(elem));
      }
    }
    //}
  }

  /**
   * Assigns a value to all cells that satisfy a condition.
   *
   * @param cond
   *            a condition.
   *
   * @param value
   *            a value (re=value[0], im=value[1]).
   * @return <tt>this</tt> (for convenience only).
   *
   */
  void fillWhere(final cfunc.ComplexProcedure cond, final Float64List value) {
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int i = firstIdx; i < lastIdx; i++) {
            final elem = get(i);
            if (cond(elem) == true) {
              setQuick(i, value);
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    Float64List elem;
    for (int i = 0; i < _size; i++) {
      elem = get(i);
      if (cond(elem) == true) {
        set(i, value);
      }
    }
    //}
  }

  /**
   * Assigns the result of a function to the real part of the receiver. The
   * imaginary part of the receiver is reset to zero.
   *
   * @param f
   *            a function object taking as argument the current cell's value.
   * @return <tt>this</tt> (for convenience only).
   * @see cern.jet.math.tdcomplex.ComplexFunctions
   */
  void forEachReal(final cfunc.ComplexRealFunction f) {
    int size = this.length;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, size);
      List<Future> futures = new List<Future>(nthreads);
      int k = size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int i = firstIdx; i < lastIdx; i++) {
            setPartsQuick(i, f(get(i)), 0.0);
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    for (int i = 0; i < size; i++) {
      setParts(i, f(get(i)), 0.0);
    }
    //}
  }

  /**
   * Replaces all cell values of the receiver with the values of another
   * matrix. Both matrices must have the same size. If both matrices share the
   * same cells (as is the case if they are views derived from the same
   * matrix) and intersect in an ambiguous way, then replaces <i>as if</i>
   * using an intermediate auxiliary deep copy of <tt>other</tt>.
   *
   * @param other
   *            the source matrix to copy from (may be identical to the
   *            receiver).
   * @return <tt>this</tt> (for convenience only).
   * @throws ArgumentError
   *             if <tt>size() != other.size()</tt>.
   */
  void copyFrom(AbstractComplexVector other) {
    if (other == this) {
      return;
    }
    checkSize(other);
    AbstractComplexVector otherLoc;
    if (_haveSharedCells(other)) {
      otherLoc = other.copy();
    } else {
      otherLoc = other;
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
          for (int i = firstIdx; i < lastIdx; i++) {
            setQuick(i, otherLoc.get(i));
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int i = 0; i < _size; i++) {
        set(i, otherLoc.get(i));
      }
    //}
  }

  /**
   * Assigns the result of a function to each cell;
   *
   * @param y
   *            the secondary matrix to operate on.
   * @param f
   *            a function object taking as first argument the current cell's
   *            value of <tt>this</tt>, and as second argument the current
   *            cell's value of <tt>y</tt>,
   * @return <tt>this</tt> (for convenience only).
   * @throws ArgumentError
   *             if <tt>size() != y.size()</tt>.
   * @see cern.jet.math.tdcomplex.ComplexFunctions
   */
  void forEachWith(final AbstractComplexVector y, final cfunc.ComplexComplexComplexFunction f) {
    int size = this.length;
    checkSize(y);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, size);
      List<Future> futures = new List<Future>(nthreads);
      int k = size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int i = firstIdx; i < lastIdx; i++) {
            setQuick(i, f(get(i), y.get(i)));
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int i = 0; i < size; i++) {
        set(i, f(get(i), y.get(i)));
      }
    //}
  }

  /**
   * Sets all cells to the state specified by <tt>re</tt> and <tt>im</tt>.
   *
   * @param re
   *            the real part of the value to be filled into the cells.
   * @param im
   *            the imaginary part of the value to be filled into the cells.
   *
   * @return <tt>this</tt> (for convenience only).
   */
  void fill(final double re, final double im) {
    int size = this.length;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, size);
      List<Future> futures = new List<Future>(nthreads);
      int k = size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int i = firstIdx; i < lastIdx; i++) {
            setPartsQuick(i, re, im);
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int i = 0; i < size; i++) {
        setParts(i, re, im);
      }
    //}
  }

  /**
   * Sets all cells to the state specified by <tt>values</tt>. <tt>values</tt>
   * is required to have the same number of cells as the receiver. Complex
   * data is represented by 2 double values in sequence: the real and
   * imaginary parts, i.e. input array must be of size 2*size().
   * <p>
   * The values are copied. So subsequent changes in <tt>values</tt> are not
   * reflected in the matrix, and vice-versa.
   *
   * @param values
   *            the values to be filled into the cells.
   * @return <tt>this</tt> (for convenience only).
   * @throws ArgumentError
   *             if <tt>values.length != 2*size()</tt>.
   */
  void setAll(final Float64List values) {
    int size = this.length;
    if (values.length != 2 * size) {
      throw new ArgumentError("The length of values[] must be equal to 2*size()=$size");
    }
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, size);
      List<Future> futures = new List<Future>(nthreads);
      int k = size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int i = firstIdx; i < lastIdx; i++) {
            setPartsQuick(i, values[2 * i], values[2 * i + 1]);
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int i = 0; i < size; i++) {
        setParts(i, values[2 * i], values[2 * i + 1]);
      }
    //}
  }

  /**
   * Replaces imaginary part of the receiver with the values of another real
   * matrix. The real part remains unchanged. Both matrices must have the same
   * size.
   *
   * @param other
   *            the source matrix to copy from
   * @return <tt>this</tt> (for convenience only).
   * @throws ArgumentError
   *             if <tt>size() != other.size()</tt>.
   */
  void setImaginary(final AbstractDoubleVector other) {
    checkSize(other);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int i = firstIdx; i < lastIdx; i++) {
            double re = get(i)[0];
            double im = other.get(i);
            setPartsQuick(i, re, im);
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int i = 0; i < _size; i++) {
        double re = get(i)[0];
        double im = other.get(i);
        setParts(i, re, im);
      }
    //}
  }

  /**
   * Replaces real part of the receiver with the values of another real
   * matrix. The imaginary part remains unchanged. Both matrices must have the
   * same size.
   *
   * @param other
   *            the source matrix to copy from
   * @return <tt>this</tt> (for convenience only).
   * @throws ArgumentError
   *             if <tt>size() != other.size()</tt>.
   */
  void setReal(final AbstractDoubleVector other) {
    checkSize(other);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int i = firstIdx; i < lastIdx; i++) {
            double re = other.get(i);
            double im = get(i)[1];
            setPartsQuick(i, re, im);
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int i = 0; i < _size; i++) {
        double re = other.get(i);
        double im = get(i)[1];
        setParts(i, re, im);
      }
    //}
  }

  /**
   * Returns the number of cells having non-zero values; ignores tolerance.
   *
   * @return the number of cells having non-zero values.
   */
  int get cardinality {
    int size = this.length;
    int cardinality = 0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, size);
      List<Future> futures = new List<Future>(nthreads);
      List<int> results = new List<int>(nthreads);
      int k = size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int cardinality = 0;
          Float64List tmp = new Float64List(2);
          for (int i = firstIdx; i < lastIdx; i++) {
            tmp = get(i);
            if ((tmp[0] != 0.0) || (tmp[1] != 0.0)) cardinality++;
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
    Float64List tmp = new Float64List(2);
    for (int i = 0; i < size; i++) {
      tmp = get(i);
      if ((tmp[0] != 0.0) || (tmp[1] != 0.0)) cardinality++;
    }
    //}
    return cardinality;
  }

  /**
   * Constructs and returns a deep copy of the receiver.
   * <p>
   * <b>Note that the returned matrix is an independent deep copy.</b> The
   * returned matrix is not backed by this matrix, so changes in the returned
   * matrix are not reflected in this matrix, and vice-versa.
   *
   * @return a deep copy of the receiver.
   */
  AbstractComplexVector copy() {
    AbstractComplexVector copy = like();
    copy.copyFrom(this);
    return copy;
  }

  /**
   * Returns whether all cells are equal to the given value.
   *
   * @param value
   *            the value to test against (re=value[0], im=value[1]).
   *
   * @return <tt>true</tt> if all cells are equal to the given value,
   *         <tt>false</tt> otherwise.
   */
  /*bool equalsValue(Float64List value) {
    return ComplexProperty.DEFAULT.equalsValue1D(this, value);
  }*/

  /**
   * Compares this object against the specified object. The result is
   * <code>true</code> if and only if the argument is not <code>null</code>
   * and is at least a <code>ComplexVector</code> object that has the same
   * sizes as the receiver and has exactly the same values at the same
   * indexes.
   *
   * @param obj
   *            the object to compare with.
   * @return <code>true</code> if the objects are the same; <code>false</code>
   *         otherwise.
   */
  bool operator ==(var obj) {
    if (obj is num) {
      final c = new Float64List.fromList([obj.toDouble, 0.0]);
      return cprop.equalsValue1D(this, c);
    }
    if (obj is Float64List) {
      return cprop.equalsValue1D(this, obj);
    }
    if (identical(this, obj)) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    if (obj is! AbstractComplexVector) {
      return false;
    }

    return cprop.equalsVector(this, obj as AbstractComplexVector);
  }

  /**
   * Returns the matrix cell value at coordinate <tt>index</tt>.
   *
   * @param index
   *            the index of the cell.
   * @return the value of the specified cell.
   * @throws IndexOutOfBoundsException
   *             if <tt>index&lt;0 || index&gt;=size()</tt>.
   */
  Float64List operator [](int index) {
    int size = this.length;
    if (index < 0 || index >= size) {
      _checkIndex(index);
    }
    return get(index);
  }

  /**
   * Returns the elements of this matrix.
   *
   * @return the elements
   */
  Object elements();

  /**
   * Returns the imaginary part of this matrix
   *
   * @return the imaginary part
   */
  AbstractDoubleVector imaginary();

  /**
   * Fills the coordinates and values of cells having non-zero values into the
   * specified lists. Fills into the lists, starting at index 0. After this
   * call returns the specified lists all have a new size, the number of
   * non-zero values.
   * <p>
   * In general, fill order is <i>unspecified</i>. This implementation fills
   * like: <tt>for (index = 0..size()-1)  do ... </tt>. However, subclasses
   * are free to us any other order, even an order that may change over time
   * as cell values are changed. (Of course, result lists indexes are
   * guaranteed to correspond to the same cell).
   *
   * @param indexList
   *            the list to be filled with indexes, can have any size.
   * @param valueList
   *            the list to be filled with values, can have any size.
   */
  void nonZeros(final List<int> indexList, final List<Float64List> valueList) {
    indexList.clear();
    valueList.clear();
    int s = length;
    for (int i = 0; i < s; i++) {
      Float64List value = get(i);
      if (value[0] != 0 || value[1] != 0) {
        indexList.add(i);
        valueList.add(value);
      }
    }

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
  Float64List get(int index);

  /**
   * Returns the real part of this matrix
   *
   * @return the real part
   */
  AbstractDoubleVector real();

  /**
   * Construct and returns a new empty matrix <i>of the same dynamic type</i>
   * as the receiver, having the same size. For example, if the receiver is an
   * instance of type <tt>DenseComplexVector</tt> the new matrix must also
   * be of type <tt>DenseComplexVector</tt>. In general, the new matrix
   * should have internal parametrization as similar as possible.
   *
   * @return a new empty matrix of the same dynamic type.
   */
  AbstractComplexVector like() {
    int size = this.length;
    return like1D(size);
  }

  /**
   * Construct and returns a new empty matrix <i>of the same dynamic type</i>
   * as the receiver, having the specified size. For example, if the receiver
   * is an instance of type <tt>DenseComplexVector</tt> the new matrix must
   * also be of type <tt>DenseComplexVector</tt>. In general, the new
   * matrix should have internal parametrization as similar as possible.
   *
   * @param size
   *            the number of cell the matrix shall have.
   * @return a new empty matrix of the same dynamic type.
   */
  AbstractComplexVector like1D(int size);

  /**
   * Construct and returns a new 2-d matrix <i>of the corresponding dynamic
   * type</i>, entirely independent of the receiver. For example, if the
   * receiver is an instance of type <tt>DenseComplexVector</tt> the new
   * matrix must be of type <tt>DenseComplexMatrix</tt>.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @return a new matrix of the corresponding dynamic type.
   */
  AbstractComplexMatrix like2D(int rows, int columns);

  /**
   * Returns new DoubleMatrix of size rows x columns whose elements are
   * taken column-wise from this matrix.
   *
   * @param rows
   *            number of rows
   * @param columns
   *            number of columns
   * @return new 2D matrix with columns being the elements of this matrix.
   */
  AbstractComplexMatrix reshape(int rows, int columns);

  /**
   * Returns new DoubleMatrix3D of size slices x rows x columns, whose
   * elements are taken column-wise from this matrix.
   *
   * @param rows
   *            number of rows
   * @param columns
   *            number of columns
   * @return new 2D matrix with columns being the elements of this matrix.
   */
  //ComplexMatrix3D reshape(int slices, int rows, int columns);

  /**
   * Sets the matrix cell at coordinate <tt>index</tt> to the specified value.
   *
   * @param index
   *            the index of the cell.
   * @param re
   *            the real part of the value to be filled into the specified
   *            cell.
   * @param im
   *            the imaginary part of the value to be filled into the
   *            specified cell.
   *
   * @throws IndexOutOfBoundsException
   *             if <tt>index&lt;0 || index&gt;=size()</tt>.
   */
  void put(int index, double re, double im) {
    int size = this.length;
    if (index < 0 || index >= size) _checkIndex(index);
    setParts(index, re, im);
  }

  /**
   * Sets the matrix cell at coordinate <tt>index</tt> to the specified value.
   *
   * @param index
   *            the index of the cell.
   * @param value
   *            the value to be filled into the specified cell (re=value[0],
   *            im=value[1]).
   *
   * @throws IndexOutOfBoundsException
   *             if <tt>index&lt;0 || index&gt;=size()</tt>.
   */
  void operator []=(int index, Float64List value) {
    int size = this.length;
    if (index < 0 || index >= size) _checkIndex(index);
    set(index, value);
  }

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
   * @param re
   *            the real part of the value to be filled into the specified
   *            cell.
   * @param im
   *            the imaginary part of the value to be filled into the
   *            specified cell.
   */
  void setParts(int index, double re, double im);

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
   *            the value to be filled into the specified cell (re=value[0],
   *            im=value[1]).
   */
  void set(int index, Float64List value);

  /**
   * Swaps each element <tt>this[i]</tt> with <tt>other[i]</tt>.
   *
   * @throws ArgumentError
   *             if <tt>size() != other.size()</tt>.
   */
  void swap(final AbstractComplexVector other) {
    int size = this.length;
    checkSize(other);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, size);
      List<Future> futures = new List<Future>(nthreads);
      int k = size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          Float64List tmp;
          for (int i = firstIdx; i < lastIdx; i++) {
            tmp = get(i);
            setQuick(i, other.get(i));
            other.setQuick(i, tmp);
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      Float64List tmp;
      for (int i = 0; i < size; i++) {
        tmp = get(i);
        set(i, other.get(i));
        other.set(i, tmp);
      }
    //}
  }

  /**
   * Constructs and returns a 1-dimensional array containing the cell values.
   * The values are copied. So subsequent changes in <tt>values</tt> are not
   * reflected in the matrix, and vice-versa. The returned array
   * <tt>values</tt> has the form <br>
   * <tt>for (int i = 0; i < size; i++) {
   *      tmp = get(i);
   *      values[2 * i] = tmp[0]; //real part
   *      values[2 * i + 1] = tmp[1]; //imaginary part
   *     }</tt>
   *
   * @return an array filled with the values of the cells.
   */
  Float64List toList() {
    int size = this.length;
    Float64List values = new Float64List(2 * size);
    fillList(values);
    return values;
  }

  /**
   * Fills the cell values into the specified 1-dimensional array. The values
   * are copied. So subsequent changes in <tt>values</tt> are not reflected in
   * the matrix, and vice-versa. After this call returns the array
   * <tt>values</tt> has the form <br>
   * <tt>for (int i = 0; i < size; i++) {
   *      tmp = get(i);
   *      values[2 * i] = tmp[0]; //real part
   *      values[2 * i + 1] = tmp[1]; //imaginary part
   *     }</tt>
   *
   * @throws ArgumentError
   *             if <tt>values.length < 2*size()</tt>.
   */
  void fillList(final Float64List values) {
    int size = this.length;
    if (values.length < 2 * size) throw new ArgumentError("values too small");
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, size);
      List<Future> futures = new List<Future>(nthreads);
      int k = size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          Float64List tmp;
          for (int i = firstIdx; i < lastIdx; i++) {
            tmp = get(i);
            values[2 * i] = tmp[0];
            values[2 * i + 1] = tmp[1];
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      Float64List tmp;
      for (int i = 0; i < size; i++) {
        tmp = get(i);
        values[2 * i] = tmp[0];
        values[2 * i + 1] = tmp[1];
      }
    //}
  }

  /**
   * Returns a string representation using default formatting ("%.4f").
   *
   * @return a string representation of the matrix.
   */
  String toString() {
    return toStringFormat("%.4f");
  }

  /**
   * Returns a string representation using given <tt>format</tt>
   *
   * @param format
   *            a format for java.lang.String.format().
   * @return a string representation of the matrix.
   */
  String toStringFormat(String format) {
    final f = new NumberFormat(format);
    StringBuffer s = new StringBuffer("ComplexVector: ${length} elements\n\n");
    Float64List elem = new Float64List(2);
    for (int i = 0; i < length; i++) {
      elem = get(i);
      if (elem[1] == 0) {
        s.write(f.format(elem[0]) + "\n");
        continue;
      }
      if (elem[0] == 0) {
        s.write(f.format(elem[1]) + "i\n");
        continue;
      }
      if (elem[1] < 0) {
        s.write(f.format(elem[0]) + " - " + f.format(-elem[1]) + "i\n");
        continue;
      }
      s.write(f.format(elem[0]) + " + " + f.format(elem[1]) + "i\n");
    }
    return s.toString();
  }

  /**
   * Constructs and returns a new <i>flip view</i>. What used to be index
   * <tt>0</tt> is now index <tt>size()-1</tt>, ..., what used to be index
   * <tt>size()-1</tt> is now index <tt>0</tt>. The returned view is backed by
   * this matrix, so changes in the returned view are reflected in this
   * matrix, and vice-versa.
   *
   * @return a new flip view.
   */
  AbstractComplexVector flip() {
    return _view().._vFlip();
  }

  /**
   * Constructs and returns a new <i>sub-range view</i> that is a
   * <tt>width</tt> sub matrix starting at <tt>index</tt>.
   *
   * Operations on the returned view can only be applied to the restricted
   * range. Any attempt to access coordinates not contained in the view will
   * throw an <tt>IndexOutOfBoundsException</tt>.
   * <p>
   * <b>Note that the view is really just a range restriction:</b> The
   * returned matrix is backed by this matrix, so changes in the returned
   * matrix are reflected in this matrix, and vice-versa.
   * <p>
   * The view contains the cells from <tt>index..index+width-1</tt>. and has
   * <tt>view.size() == width</tt>. A view's legal coordinates are again zero
   * based, as usual. In other words, legal coordinates of the view are
   * <tt>0 .. view.size()-1==width-1</tt>. As usual, any attempt to access a
   * cell at other coordinates will throw an
   * <tt>IndexOutOfBoundsException</tt>.
   *
   * @param index
   *            The index of the first cell.
   * @param width
   *            The width of the range.
   * @throws IndexOutOfBoundsException
   *             if <tt>index<0 || width<0 || index+width>size()</tt>.
   * @return the new view.
   *
   */
  AbstractComplexVector part(int index, int width) {
    return _view().._vPart(index, width);
  }

  /**
   * Constructs and returns a new <i>selection view</i> that is a matrix
   * holding the cells matching the given condition. Applies the condition to
   * each cell and takes only those cells where
   * <tt>condition(get(i))</tt> yields <tt>true</tt>.
   *
   * The returned view is backed by this matrix, so changes in the returned
   * view are reflected in this matrix, and vice-versa.
   *
   * @param condition
   *            The condition to be matched.
   * @return the new view.
   */
  AbstractComplexVector where(cfunc.ComplexProcedure condition) {
    int size = this.length;
    List<int> matches = new List<int>();
    for (int i = 0; i < size; i++) {
      if (condition(get(i))) {
        matches.add(i);
      }
    }
    //matches.trimToSize();
    return select(new Int32List.fromList(matches));
  }

  /**
   * Constructs and returns a new <i>selection view</i> that is a matrix
   * holding the indicated cells. There holds
   * <tt>view.size() == indexes.length</tt> and
   * <tt>view.get(i) == this.get(indexes[i])</tt>. Indexes can occur multiple
   * times and can be in arbitrary order.
   *
   * Note that modifying <tt>indexes</tt> after this call has returned has no
   * effect on the view. The returned view is backed by this matrix, so
   * changes in the returned view are reflected in this matrix, and
   * vice-versa.
   *
   * @param indexes
   *            The indexes of the cells that shall be visible in the new
   *            view. To indicate that <i>all</i> cells shall be visible,
   *            simply set this parameter to <tt>null</tt>.
   * @return the new view.
   * @throws IndexOutOfBoundsException
   *             if <tt>!(0 <= indexes[i] < size())</tt> for any
   *             <tt>i=0..indexes.length()-1</tt>.
   */
  AbstractComplexVector select(Int32List indexes) {
    // check for "all"
    int size = this.length;
    if (indexes == null) {
      indexes = new Int32List(size);
      for (int i = size - 1; --i >= 0; ) indexes[i] = i;
    }

    _checkIndexes(indexes);
    Int32List offsets = new Int32List(indexes.length);
    for (int i = 0; i < indexes.length; i++) {
      offsets[i] = index(indexes[i]);
    }
    return _viewSelectionLike(offsets);
  }

  /**
   * Constructs and returns a new <i>stride view</i> which is a sub matrix
   * consisting of every i-th cell. More specifically, the view has size
   * <tt>this.size()/stride</tt> holding cells <tt>this.get(i*stride)</tt> for
   * all <tt>i = 0..size()/stride - 1</tt>.
   *
   * @param stride
   *            the step factor.
   * @throws IndexOutOfBoundsException
   *             if <tt>stride <= 0</tt>.
   * @return the new view.
   *
   */
  AbstractComplexVector strides(int stride) {
    return _view().._vStrides(stride);
  }

  /**
   * Returns the dot product of two vectors x and y. Operates on cells at
   * indexes <tt>0 .. Math.min(size(),y.size())</tt>.
   *
   * @param y
   *            the second vector.
   * @return the sum of products.
   */
  /*Float64List zDotProduct(ComplexVector y) {
        int size = this.size();
        return zDotProduct(y, 0, size);
    }*/

  /**
   * Returns the dot product of two vectors x and y. Operates on cells at
   * indexes <tt>from .. Min(size(),y.size(),from+length)-1</tt>.
   *
   * @param y
   *            the second vector.
   * @param from
   *            the first index to be considered.
   * @param length
   *            the number of cells to be considered.
   * @return the sum of products; zero if <tt>from<0 || length<0</tt>.
   */
  Float64List dot(final AbstractComplexVector y, [final int from = 0, int length = null]) {
    if (length == null) {
      length = this.length;
    }
    int size = this.length;
    if (from < 0 || length <= 0) return new Float64List.fromList([0.0, 0.0]);

    int tail = from + length;
    if (size < tail) tail = size;
    if (y._size < tail) tail = y._size;
    length = tail - from;
    Float64List sum = new Float64List(2);

    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, length);
      List<Future> futures = new List<Future>(nthreads);
      List<Float64List> results = new List<Float64List>(nthreads);//[2];
      int k = length / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? length : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          Float64List sum = new Float64List(2);
          Float64List tmp;
          int idx;
          for (int k = firstIdx; k < lastIdx; k++) {
            idx = k + from;
            tmp = y.get(idx);
            tmp[1] = -tmp[1]; // complex conjugate
            sum = Complex.plus(sum, Complex.mult(tmp, get(idx)));
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
          sum = Complex.plus(sum, results[j]);
        }
      } on ExecutionException catch (ex) {
        ex.printStackTrace();
      } on InterruptedException catch (e) {
        e.printStackTrace();
      }
    } else {*/
    Float64List tmp;
    int idx;
    for (int k = 0; k < length; k++) {
      idx = k + from;
      tmp = y.get(idx);
      tmp[1] = -tmp[1]; // complex conjugate
      sum = Complex.plus(sum, Complex.multiply(tmp, get(idx)));
    }
    //}
    return sum;
  }

  /**
   * Returns the dot product of two vectors x and y.
   *
   * @param y
   *            the second vector.
   * @param nonZeroIndexes
   *            the indexes of cells in <tt>y</tt>having a non-zero value.
   * @return the sum of products.
   */
  Float64List dotNonZero(AbstractComplexVector y, Int32List nonZeroIndexes, [int from = 0, int length = null]) {
    if (length == null) {
      length = this.length;
    }
    int size = this.length;
    if (from < 0 || length <= 0) {
      return new Float64List.fromList([0.0, 0.0]);
    }

    int tail = from + length;
    if (size < tail) tail = size;
    if (y._size < tail) tail = y._size;
    length = tail - from;
    if (length <= 0) {
      return new Float64List.fromList([0.0, 0.0]);
    }

    // setup
    //List<int> indexesCopy = new List<int>.from(nonZeroIndexes);
    Int32List indexesCopy = new Int32List.fromList(nonZeroIndexes);
    //indexesCopy.trimToSize();
    indexesCopy.sort();
    Int32List nonZeroIndexElements = indexesCopy;//.elements();
    int index = 0;
    int s = indexesCopy.length;

    // skip to start
    while ((index < s) && nonZeroIndexElements[index] < from) {
      index++;
    }

    // now the sparse dot product
    int i;
    Float64List sum = new Float64List(2);
    Float64List tmp;
    while ((--length >= 0) && (index < s) && ((i = nonZeroIndexElements[index]) < tail)) {
      tmp = y.get(i);
      tmp[1] = -tmp[1]; // complex conjugate
      sum = Complex.plus(sum, Complex.multiply(tmp, get(i)));
      index++;
    }

    return sum;
  }

  /**
   * Returns the sum of all cells.
   *
   * @return the sum.
   */
  Float64List sum() {
    Float64List sum = new Float64List(2);
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
          for (int k = firstIdx; k < lastIdx; k++) {
            sum = Complex.plus(sum, get(k));
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
          sum[0] += results[j][0];
          sum[1] += results[j][1];
        }
      } on ExecutionException catch (ex) {
        ex.printStackTrace();
      } on InterruptedException catch (e) {
        e.printStackTrace();
      }
    } else {*/
    for (int k = 0; k < _size; k++) {
      sum = Complex.plus(sum, get(k));
    }
    //}
    return sum;
  }

  /**
   * Returns the number of cells having non-zero values, but at most
   * maxCardinality; ignores tolerance.
   *
   * @param maxCardinality
   *            maximal cardinality
   * @return number of cells having non-zero values, but at most
   *         maxCardinality.
   */
  int _cardinality(int maxCardinality) {
    int size = this.length;
    int cardinality = 0;
    int i = 0;
    Float64List tmp = new Float64List(2);
    while (i++ < size && cardinality < maxCardinality) {
      tmp = get(i);
      if ((tmp[0] != 0.0) || (tmp[1] != 0.0)) cardinality++;
    }
    return cardinality;
  }

  /**
   * Returns the content of this matrix if it is a wrapper; or <tt>this</tt>
   * otherwise. Override this method in wrappers.
   */
  AbstractComplexVector _getContent() {
    return this;
  }

  /**
   * Returns <tt>true</tt> if both matrices share at least one identical cell.
   *
   * @param other
   *            matrix
   * @return <tt>true</tt> if both matrices share at least one identical cell
   */
  bool _haveSharedCells(AbstractComplexVector other) {
    if (other == null) return false;
    if (this == other) return true;
    return _getContent()._haveSharedCellsRaw(other._getContent());
  }

  /**
   * Always returns false
   *
   * @param other
   *            matrix
   * @return false
   */
  bool _haveSharedCellsRaw(AbstractComplexVector other) {
    return false;
  }

  /**
   * Constructs and returns a new view equal to the receiver. The view is a
   * shallow clone. Calls <code>clone()</code> and casts the result.
   * <p>
   * <b>Note that the view is not a deep copy.</b> The returned matrix is
   * backed by this matrix, so changes in the returned matrix are reflected in
   * this matrix, and vice-versa.
   * <p>
   * Use {@link #copy()} to construct an independent deep copy rather than a
   * new view.
   *
   * @return a new view of the receiver.
   */
  AbstractComplexVector _view() {
    return clone() as AbstractComplexVector;
  }

  /**
   * Construct and returns a new selection view.
   *
   * @param offsets
   *            the offsets of the visible elements.
   * @return a new view.
   */
  AbstractComplexVector _viewSelectionLike(Int32List offsets);

  /**
   * Returns the dot product of two vectors x and y.
   *
   * @param y
   *            the second vector.
   * @param nonZeroIndexes
   *            the indexes of cells in <tt>y</tt>having a non-zero value.
   * @return the sum of products.
   */
  /*Float64List _zDotProduct(ComplexVector y, IntArrayList nonZeroIndexes) {
    return zDotProduct(y, 0, size(), nonZeroIndexes);
  }*/

  Object clone();

  AbstractComplexVector operator *(AbstractComplexVector y) {
    return this.copy()..forEachWith(y, cfunc.mult);
  }

  AbstractComplexVector operator /(AbstractComplexVector y) {
    return this.copy()..forEachWith(y, cfunc.div);
  }

  AbstractComplexVector operator +(AbstractComplexVector y) {
    return this.copy()..forEachWith(y, cfunc.plus);
  }

  AbstractComplexVector operator -(AbstractComplexVector y) {
    return this.copy()..forEachWith(y, cfunc.minus);
  }

  AbstractComplexVector conj() {
    return this.copy()..forEach(cfunc.conj);
  }
}
