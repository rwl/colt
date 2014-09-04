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
abstract class DComplexMatrix1D extends AbstractMatrix1D {

  /**
   * Makes this class non instantiable, but still let's others inherit from
   * it.
   */
  /*DComplexMatrix1D() {
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
   * @see cern.jet.math.tdcomplex.DComplexFunctions
   */
  Float64List aggregate(final cfunc.DComplexDComplexDComplexFunction aggr, final cfunc.DComplexDComplexFunction f) {
    Float64List b = new Float64List(2);
    int size = this.length;
    if (size == 0) {
      b[0] = double.NAN;
      b[1] = double.NAN;
      return b;
    }
    Float64List a = f(getQuick(0));
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, size);
      List<Future> futures = new List<Future>(nthreads);
      int k = size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          Float64List a = f(getQuick(firstIdx));
          for (int i = firstIdx + 1; i < lastIdx; i++) {
            a = aggr(a, f(getQuick(i)));
          }
          return a;
        });
      }
      a = ConcurrencyUtils.waitForCompletion(futures, aggr);
    } else {*/
    for (int i = 1; i < size; i++) {
      a = aggr(a, f(getQuick(i)));
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
   * @see cern.jet.math.tdcomplex.DComplexFunctions
   */
  Float64List aggregateMatrix(final DComplexMatrix1D other, final cfunc.DComplexDComplexDComplexFunction aggr, final cfunc.DComplexDComplexDComplexFunction f) {
    checkSize(other);
    int size = this.length;
    if (size == 0) {
      Float64List b = new Float64List(2);
      b[0] = double.NAN;
      b[1] = double.NAN;
      return b;
    }
    Float64List a = f(getQuick(0), other.getQuick(0));
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, size);
      List<Future> futures = new List<Future>(nthreads);
      int k = size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          Float64List a = f(getQuick(firstIdx), other.getQuick(firstIdx));
          for (int i = firstIdx + 1; i < lastIdx; i++) {
            a = aggr(a, f(getQuick(i), other.getQuick(i)));
          }
          return a;
        });
      }
      a = ConcurrencyUtils.waitForCompletion(futures, aggr);
    } else {*/
    for (int i = 1; i < size; i++) {
      a = aggr(a, f(getQuick(i), other.getQuick(i)));
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
   * @see cern.jet.math.tdcomplex.DComplexFunctions
   */
  DComplexMatrix1D assign(final cfunc.DComplexDComplexFunction f) {
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
            setQuick(i, f(getQuick(i)));
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    for (int i = 0; i < size; i++) {
      setQuick(i, f(getQuick(i)));
    }
    //}
    return this;
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
   * @see cern.jet.math.tdcomplex.DComplexFunctions
   */
  DComplexMatrix1D assignProc(final cfunc.DComplexProcedure cond, final cfunc.DComplexDComplexFunction f) {
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
            Float64List elem = getQuick(i);
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
      elem = getQuick(i);
      if (cond(elem) == true) {
        setQuick(i, f(elem));
      }
    }
    //}
    return this;
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
          for (int i = firstIdx; i < lastIdx; i++) {
            final elem = getQuick(i);
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
      elem = getQuick(i);
      if (cond(elem) == true) {
        setQuick(i, value);
      }
    }
    //}
    return this;
  }

  /**
   * Assigns the result of a function to the real part of the receiver. The
   * imaginary part of the receiver is reset to zero.
   *
   * @param f
   *            a function object taking as argument the current cell's value.
   * @return <tt>this</tt> (for convenience only).
   * @see cern.jet.math.tdcomplex.DComplexFunctions
   */
  DComplexMatrix1D assignRealFunc(final cfunc.DComplexRealFunction f) {
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
            setPartsQuick(i, f(getQuick(i)), 0.0);
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    for (int i = 0; i < size; i++) {
      setPartsQuick(i, f(getQuick(i)), 0.0);
    }
    //}
    return this;
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
  DComplexMatrix1D assignMatrix(DComplexMatrix1D other) {
    if (other == this) return this;
    checkSize(other);
    DComplexMatrix1D otherLoc;
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
            setQuick(i, otherLoc.getQuick(i));
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int i = 0; i < _size; i++) {
        setQuick(i, otherLoc.getQuick(i));
      }
    //}
    return this;
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
   * @see cern.jet.math.tdcomplex.DComplexFunctions
   */
  DComplexMatrix1D assignMatrixFunc(final DComplexMatrix1D y, final cfunc.DComplexDComplexDComplexFunction f) {
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
            setQuick(i, f(getQuick(i), y.getQuick(i)));
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int i = 0; i < size; i++) {
        setQuick(i, f(getQuick(i), y.getQuick(i)));
      }
    //}
    return this;
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
  DComplexMatrix1D assignValue(final double re, final double im) {
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
        setPartsQuick(i, re, im);
      }
    //}
    return this;
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
  DComplexMatrix1D assignValues(final Float64List values) {
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
        setPartsQuick(i, values[2 * i], values[2 * i + 1]);
      }
    //}
    return this;
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
  DComplexMatrix1D assignImaginary(final DoubleMatrix1D other) {
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
            double re = getQuick(i)[0];
            double im = other.getQuick(i);
            setPartsQuick(i, re, im);
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int i = 0; i < _size; i++) {
        double re = getQuick(i)[0];
        double im = other.get(i);
        setPartsQuick(i, re, im);
      }
    //}
    return this;
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
  DComplexMatrix1D assignReal(final DoubleMatrix1D other) {
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
            double re = other.getQuick(i);
            double im = getQuick(i)[1];
            setPartsQuick(i, re, im);
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int i = 0; i < _size; i++) {
        double re = other.get(i);
        double im = getQuick(i)[1];
        setPartsQuick(i, re, im);
      }
    //}
    return this;
  }

  /**
   * Returns the number of cells having non-zero values; ignores tolerance.
   *
   * @return the number of cells having non-zero values.
   */
  int cardinality() {
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
            tmp = getQuick(i);
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
      tmp = getQuick(i);
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
  DComplexMatrix1D copy() {
    DComplexMatrix1D copy = like();
    copy.assignMatrix(this);
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
  bool equalsValue(Float64List value) {
    return DComplexProperty.DEFAULT.equalsValue1D(this, value);
  }

  /**
   * Compares this object against the specified object. The result is
   * <code>true</code> if and only if the argument is not <code>null</code>
   * and is at least a <code>ComplexMatrix1D</code> object that has the same
   * sizes as the receiver and has exactly the same values at the same
   * indexes.
   *
   * @param obj
   *            the object to compare with.
   * @return <code>true</code> if the objects are the same; <code>false</code>
   *         otherwise.
   */

  bool equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (!(obj is DComplexMatrix1D)) return false;

    return DComplexProperty.DEFAULT.equalsMatrix1D(this, obj as DComplexMatrix1D);
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
  Float64List get(int index) {
    int size = this.length;
    if (index < 0 || index >= size) _checkIndex(index);
    return getQuick(index);
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
  DoubleMatrix1D getImaginaryPart();

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
  void getNonZeros(final List<int> indexList, final List<Float64List> valueList) {
    indexList.clear();
    valueList.clear();
    int s = length;
    for (int i = 0; i < s; i++) {
      Float64List value = getQuick(i);
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
  Float64List getQuick(int index);

  /**
   * Returns the real part of this matrix
   *
   * @return the real part
   */
  DoubleMatrix1D getRealPart();

  /**
   * Construct and returns a new empty matrix <i>of the same dynamic type</i>
   * as the receiver, having the same size. For example, if the receiver is an
   * instance of type <tt>DenseComplexMatrix1D</tt> the new matrix must also
   * be of type <tt>DenseComplexMatrix1D</tt>. In general, the new matrix
   * should have internal parametrization as similar as possible.
   *
   * @return a new empty matrix of the same dynamic type.
   */
  DComplexMatrix1D like() {
    int size = this.length;
    return like1D(size);
  }

  /**
   * Construct and returns a new empty matrix <i>of the same dynamic type</i>
   * as the receiver, having the specified size. For example, if the receiver
   * is an instance of type <tt>DenseDComplexMatrix1D</tt> the new matrix must
   * also be of type <tt>DenseDComplexMatrix1D</tt>. In general, the new
   * matrix should have internal parametrization as similar as possible.
   *
   * @param size
   *            the number of cell the matrix shall have.
   * @return a new empty matrix of the same dynamic type.
   */
  DComplexMatrix1D like1D(int size);

  /**
   * Construct and returns a new 2-d matrix <i>of the corresponding dynamic
   * type</i>, entirely independent of the receiver. For example, if the
   * receiver is an instance of type <tt>DenseDComplexMatrix1D</tt> the new
   * matrix must be of type <tt>DenseDComplexMatrix2D</tt>.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @return a new matrix of the corresponding dynamic type.
   */
  DComplexMatrix2D like2D(int rows, int columns);

  /**
   * Returns new DoubleMatrix2D of size rows x columns whose elements are
   * taken column-wise from this matrix.
   *
   * @param rows
   *            number of rows
   * @param columns
   *            number of columns
   * @return new 2D matrix with columns being the elements of this matrix.
   */
  DComplexMatrix2D reshape(int rows, int columns);

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
  //DComplexMatrix3D reshape(int slices, int rows, int columns);

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
  void setParts(int index, double re, double im) {
    int size = this.length;
    if (index < 0 || index >= size) _checkIndex(index);
    setPartsQuick(index, re, im);
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
  void set(int index, Float64List value) {
    int size = this.length;
    if (index < 0 || index >= size) _checkIndex(index);
    setQuick(index, value);
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
  void setPartsQuick(int index, double re, double im);

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
  void setQuick(int index, Float64List value);

  /**
   * Swaps each element <tt>this[i]</tt> with <tt>other[i]</tt>.
   *
   * @throws ArgumentError
   *             if <tt>size() != other.size()</tt>.
   */
  void swap(final DComplexMatrix1D other) {
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
            tmp = getQuick(i);
            setQuick(i, other.getQuick(i));
            other.setQuick(i, tmp);
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      Float64List tmp;
      for (int i = 0; i < size; i++) {
        tmp = getQuick(i);
        setQuick(i, other.getQuick(i));
        other.setQuick(i, tmp);
      }
    //}
  }

  /**
   * Constructs and returns a 1-dimensional array containing the cell values.
   * The values are copied. So subsequent changes in <tt>values</tt> are not
   * reflected in the matrix, and vice-versa. The returned array
   * <tt>values</tt> has the form <br>
   * <tt>for (int i = 0; i < size; i++) {
   *      tmp = getQuick(i);
   *      values[2 * i] = tmp[0]; //real part
   *      values[2 * i + 1] = tmp[1]; //imaginary part
   *     }</tt>
   *
   * @return an array filled with the values of the cells.
   */
  Float64List toArray() {
    int size = this.length;
    Float64List values = new Float64List(2 * size);
    toArrayFill(values);
    return values;
  }

  /**
   * Fills the cell values into the specified 1-dimensional array. The values
   * are copied. So subsequent changes in <tt>values</tt> are not reflected in
   * the matrix, and vice-versa. After this call returns the array
   * <tt>values</tt> has the form <br>
   * <tt>for (int i = 0; i < size; i++) {
   *      tmp = getQuick(i);
   *      values[2 * i] = tmp[0]; //real part
   *      values[2 * i + 1] = tmp[1]; //imaginary part
   *     }</tt>
   *
   * @throws ArgumentError
   *             if <tt>values.length < 2*size()</tt>.
   */
  void toArrayFill(final Float64List values) {
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
            tmp = getQuick(i);
            values[2 * i] = tmp[0];
            values[2 * i + 1] = tmp[1];
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      Float64List tmp;
      for (int i = 0; i < size; i++) {
        tmp = getQuick(i);
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
    StringBuffer s = new StringBuffer("ComplexMatrix1D: ${length} elements\n\n");
    Float64List elem = new Float64List(2);
    for (int i = 0; i < length; i++) {
      elem = getQuick(i);
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
  DComplexMatrix1D viewFlip() {
    return _view()._vFlip() as DComplexMatrix1D;
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
  DComplexMatrix1D viewPart(int index, int width) {
    return _view()._vPart(index, width) as DComplexMatrix1D;
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
  DComplexMatrix1D viewSelectionProc(cfunc.DComplexProcedure condition) {
    int size = this.length;
    List<int> matches = new List<int>();
    for (int i = 0; i < size; i++) {
      if (condition(getQuick(i))) {
        matches.add(i);
      }
    }
    //matches.trimToSize();
    return viewSelection(new Int32List.fromList(matches));
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
  DComplexMatrix1D viewSelection(Int32List indexes) {
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
  DComplexMatrix1D viewStrides(int stride) {
    return _view()._vStrides(stride) as DComplexMatrix1D;
  }

  /**
   * Returns the dot product of two vectors x and y. Operates on cells at
   * indexes <tt>0 .. Math.min(size(),y.size())</tt>.
   *
   * @param y
   *            the second vector.
   * @return the sum of products.
   */
  /*Float64List zDotProduct(DComplexMatrix1D y) {
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
  Float64List zDotProduct(final DComplexMatrix1D y, [final int from = 0, int length = null]) {
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
            tmp = y.getQuick(idx);
            tmp[1] = -tmp[1]; // complex conjugate
            sum = DComplex.plus(sum, DComplex.mult(tmp, getQuick(idx)));
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
    Float64List tmp;
    int idx;
    for (int k = 0; k < length; k++) {
      idx = k + from;
      tmp = y.getQuick(idx);
      tmp[1] = -tmp[1]; // complex conjugate
      sum = DComplex.plus(sum, DComplex.multiply(tmp, getQuick(idx)));
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
  Float64List zDotProductIndex(DComplexMatrix1D y, Int32List nonZeroIndexes, [int from = 0, int length = null]) {
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
      tmp = y.getQuick(i);
      tmp[1] = -tmp[1]; // complex conjugate
      sum = DComplex.plus(sum, DComplex.multiply(tmp, getQuick(i)));
      index++;
    }

    return sum;
  }

  /**
   * Returns the sum of all cells.
   *
   * @return the sum.
   */
  Float64List zSum() {
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
            sum = DComplex.plus(sum, getQuick(k));
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
      sum = DComplex.plus(sum, getQuick(k));
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
      tmp = getQuick(i);
      if ((tmp[0] != 0.0) || (tmp[1] != 0.0)) cardinality++;
    }
    return cardinality;
  }

  /**
   * Returns the content of this matrix if it is a wrapper; or <tt>this</tt>
   * otherwise. Override this method in wrappers.
   */
  DComplexMatrix1D _getContent() {
    return this;
  }

  /**
   * Returns <tt>true</tt> if both matrices share at least one identical cell.
   *
   * @param other
   *            matrix
   * @return <tt>true</tt> if both matrices share at least one identical cell
   */
  bool _haveSharedCells(DComplexMatrix1D other) {
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
  bool _haveSharedCellsRaw(DComplexMatrix1D other) {
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
  DComplexMatrix1D _view() {
    return clone() as DComplexMatrix1D;
  }

  /**
   * Construct and returns a new selection view.
   *
   * @param offsets
   *            the offsets of the visible elements.
   * @return a new view.
   */
  DComplexMatrix1D _viewSelectionLike(Int32List offsets);

  /**
   * Returns the dot product of two vectors x and y.
   *
   * @param y
   *            the second vector.
   * @param nonZeroIndexes
   *            the indexes of cells in <tt>y</tt>having a non-zero value.
   * @return the sum of products.
   */
  /*Float64List _zDotProduct(DComplexMatrix1D y, IntArrayList nonZeroIndexes) {
    return zDotProduct(y, 0, size(), nonZeroIndexes);
  }*/

  Object clone();
}
