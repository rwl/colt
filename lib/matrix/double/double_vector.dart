/*
Copyright (C) 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
is hereby granted without fee, provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear in supporting documentation.
CERN makes no representations about the suitability of this software for any purpose.
It is provided "as is" without expressed or implied warranty.
 */
part of cern.colt.matrix;

class DoubleVectorLocation {
  final double value;
  final int location;
  DoubleVectorLocation(this.value, this.location);
}

/**
 * Abstract base class for 1-d matrices (aka <i>vectors</i>) holding
 * <tt>double</tt> elements. First see the <a
 * href="package-summary.html">package summary</a> and javadoc <a
 * href="package-tree.html">tree view</a> to get the broad picture.
 * <p>
 * A matrix has a number of cells (its <i>size</i>), which are assigned upon
 * instance construction. Elements are accessed via zero based indexes. Legal
 * indexes are of the form <tt>[0..size()-1]</tt>. Any attempt to access an
 * element at a coordinate <tt>index&lt;0 || index&gt;=size()</tt> will throw an
 * <tt>IndexOutOfBoundsException</tt>.
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 09/24/99
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
abstract class DoubleVector extends AbstractVector {

  /**
   * Makes this class non instantiable, but still let's others inherit from
   * it.
   */
  DoubleVector();

  /**
   * Applies a function to each cell and aggregates the results. Returns a
   * value <tt>v</tt> such that <tt>v==a(size())</tt> where
   * <tt>a(i) == aggr( a(i-1), f(get(i)) )</tt> and terminators are
   * <tt>a(1) == f(get(0)), a(0)==double.NAN</tt>.
   * <p>
   * <b>Example:</b>
   *
   * <pre>
   * 	 cern.jet.math.Functions F = cern.jet.math.Functions.functions;
   * 	 matrix = 0 1 2 3
   *
   * 	 // Sum( x[i]*x[i] )
   * 	 matrix.aggregate(F.plus,F.square);
   * 	 --&gt; 14
   *
   * </pre>
   *
   * For further examples, see the <a
   * href="package-summary.html#FunctionObjects">package doc</a>.
   *
   * @param aggr
   *            an aggregation function taking as first argument the current
   *            aggregation and as second argument the transformed current
   *            cell value.
   * @param f
   *            a function transforming the current cell value.
   * @return the aggregated measure.
   * @see cern.jet.math.tdouble.DoubleFunctions
   */
  double reduce(final DoubleDoubleFunction aggr, DoubleFunction f) {
    if (_size == 0) return double.NAN;
    double a = f(get(0));
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
        nthreads = Math.min(nthreads, _size);
        List<Future<double>> futures = new List<Future<double>>(nthreads);
        int k = _size ~/ nthreads;
        for (int j = 0; j < nthreads; j++) {
            final int firstIdx = j * k;
            final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
            futures[j] = ConcurrencyUtils.submitDouble(() {
                double a = f(getQuick(firstIdx));
                for (int i = firstIdx + 1; i < lastIdx; i++) {
                    a = aggr(a, f(getQuick(i)));
                }
                return a;
            });
        }
        a = ConcurrencyUtils.waitForCompletionAggr(futures, aggr);
    } else {*/
    for (int i = 1; i < _size; i++) {
      a = aggr(a, f(get(i)));
    }
    //}
    return a;
  }

  /**
   * Applies a function to all cells with a given indexes and aggregates the
   * results.
   *
   * @param aggr
   *            an aggregation function taking as first argument the current
   *            aggregation and as second argument the transformed current
   *            cell value.
   * @param f
   *            a function transforming the current cell value.
   * @param indexList
   *            indexes.
   *
   * @return the aggregated measure.
   * @see cern.jet.math.tdouble.DoubleFunctions
   */
  double reduceRange(final DoubleDoubleFunction aggr, DoubleFunction f,  /*IntArrayList*/Int32List indexList) {
    if (this.length == 0) return double.NAN;
    final int size = indexList.length;
    final Int32List indexElements = indexList;
    double a = 0.0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
        nthreads = Math.min(nthreads, size);
        List<Future> futures = new List<Future>(nthreads);
        int k = size ~/ nthreads;
        for (int j = 0; j < nthreads; j++) {
            final int firstIdx = j * k;
            final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
            futures[j] = ConcurrencyUtils.submitDouble(() {
                  double a = f(getQuick(indexElements[firstIdx]));
                  double elem;
                  for (int i = firstIdx + 1; i < lastIdx; i++) {
                      elem = getQuick(indexElements[i]);
                      a = aggr(a, f(elem));
                  }
                  return a;
            });
        }
        a = ConcurrencyUtils.waitForCompletionAggr(futures, aggr);
    } else {*/
    double elem;
    a = f(get(indexElements[0]));
    for (int i = 1; i < size; i++) {
      elem = get(indexElements[i]);
      a = aggr(a, f(elem));
    }
    //}
    return a;
  }

  /**
   * Applies a function to each corresponding cell of two matrices and
   * aggregates the results. Returns a value <tt>v</tt> such that
   * <tt>v==a(size())</tt> where
   * <tt>a(i) == aggr( a(i-1), f(get(i),other.get(i)) )</tt> and terminators
   * are <tt>a(1) == f(get(0),other.get(0)), a(0)==double.NAN</tt>.
   * <p>
   * <b>Example:</b>
   *
   * <pre>
   * 	 cern.jet.math.Functions F = cern.jet.math.Functions.functions;
   * 	 x = 0 1 2 3
   * 	 y = 0 1 2 3
   *
   * 	 // Sum( x[i]*y[i] )
   * 	 x.aggregate(y, F.plus, F.mult);
   * 	 --&gt; 14
   *
   * 	 // Sum( (x[i]+y[i])&circ;2 )
   * 	 x.aggregate(y, F.plus, F.chain(F.square,F.plus));
   * 	 --&gt; 56
   *
   * </pre>
   *
   * For further examples, see the <a
   * href="package-summary.html#FunctionObjects">package doc</a>.
   *
   * @param aggr
   *            an aggregation function taking as first argument the current
   *            aggregation and as second argument the transformed current
   *            cell values.
   * @param f
   *            a function transforming the current cell values.
   * @return the aggregated measure.
   * @throws IllegalArgumentException
   *             if <tt>size() != other.size()</tt>.
   * @see cern.jet.math.tdouble.DoubleFunctions
   */
  double reduceVector(final DoubleVector other, DoubleDoubleFunction aggr, DoubleDoubleFunction f) {
    checkSize(other);
    if (_size == 0) return double.NAN;
    double a = f(get(0), other.get(0));
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
        nthreads = Math.min(nthreads, _size);
        List<Future> futures = new List<Future>(nthreads);
        int k = _size ~/ nthreads;
        for (int j = 0; j < nthreads; j++) {
            final int firstIdx = j * k;
            final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
            futures[j] = ConcurrencyUtils.submitDouble(() {
                double a = f(getQuick(firstIdx), other.getQuick(firstIdx));
                for (int i = firstIdx + 1; i < lastIdx; i++) {
                    a = aggr(a, f(getQuick(i), other.getQuick(i)));
                }
                return a;
            });
        }
        a = ConcurrencyUtils.waitForCompletionAggr(futures, aggr);
    } else {*/
    for (int i = 1; i < _size; i++) {
      a = aggr(a, f(get(i), other.get(i)));
    }
    //}
    return a;
  }

  /**
   * Assigns the result of a function to each cell;
   * <tt>x[i] = function(x[i])</tt>. (Iterates downwards from
   * <tt>[size()-1]</tt> to <tt>[0]</tt>).
   * <p>
   * <b>Example:</b>
   *
   * <pre>
   * 	 // change each cell to its sine
   * 	 matrix =   0.5      1.5      2.5       3.5
   * 	 matrix.assign(cern.jet.math.Functions.sin);
   * 	 --&gt;
   * 	 matrix ==  0.479426 0.997495 0.598472 -0.350783
   *
   * </pre>
   *
   * For further examples, see the <a
   * href="package-summary.html#FunctionObjects">package doc</a>.
   *
   * @param f
   *            a function object taking as argument the current cell's value.
   * @return <tt>this</tt> (for convenience only).
   * @see cern.jet.math.tdouble.DoubleFunctions
   */
  DoubleVector forEach(final DoubleFunction f) {
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
                    setQuick(i, f(getQuick(i)));
                }
            });
        }
        ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    for (int i = 0; i < _size; i++) {
      set(i, f(get(i)));
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
   * @see cern.jet.math.tdouble.DoubleFunctions
   */
  DoubleVector forEachWhere(final DoubleProcedure cond, DoubleFunction f) {
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          double elem;
          for (int i = firstIdx; i < lastIdx; i++) {
            elem = getQuick(i);
            if (cond(elem) == true) {
              setQuick(i, f(elem));
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      double elem;
      for (int i = 0; i < _size; i++) {
        elem = get(i);
        if (cond(elem) == true) {
          set(i, f(elem));
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
   *            a value.
   * @return <tt>this</tt> (for convenience only).
   *
   */
  DoubleVector fillWhere(DoubleProcedure cond, double value) {
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          double elem;
          for (int i = firstIdx; i < lastIdx; i++) {
            elem = getQuick(i);
            if (cond(elem) == true) {
              setQuick(i, value);
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      double elem;
      for (int i = 0; i < _size; i++) {
        elem = get(i);
        if (cond(elem) == true) {
          set(i, value);
        }
      }
    //}
    return this;
  }

  /**
   * Sets all cells to the state specified by <tt>value</tt>.
   *
   * @param value
   *            the value to be filled into the cells.
   * @return <tt>this</tt> (for convenience only).
   */
  DoubleVector fill(final double value) {
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
            setQuick(i, value);
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int i = 0; i < _size; i++) {
        set(i, value);
      }
    //}
    return this;
  }

  /**
   * Sets all cells to the state specified by <tt>values</tt>. <tt>values</tt>
   * is required to have the same number of cells as the receiver.
   * <p>
   * The values are copied. So subsequent changes in <tt>values</tt> are not
   * reflected in the matrix, and vice-versa.
   *
   * @param values
   *            the values to be filled into the cells.
   * @return <tt>this</tt> (for convenience only).
   * @throws IllegalArgumentException
   *             if <tt>values.length != size()</tt>.
   */
  DoubleVector setValues(List<double> values) {
    if (values.length != _size) {
      throw new ArgumentError("Must have same number of cells: length=${values.length} size()=${length}");
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
            setQuick(i, values[i]);
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int i = 0; i < _size; i++) {
        set(i, values[i]);
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
   * @throws IllegalArgumentException
   *             if <tt>size() != other.size()</tt>.
   */
  DoubleVector setAll(DoubleVector other) {
    if (other == this) return this;
    checkSize(other);
    DoubleVector source;
    if (_haveSharedCells(other)) {
      source = other.copy();
    } else {
      source = other;
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
            setQuick(i, source.getQuick(i));
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int i = 0; i < _size; i++) {
        set(i, source.get(i));
      }
    //}
    return this;
  }

  /**
   * Assigns the result of a function to each cell;
   * <tt>x[i] = function(x[i],y[i])</tt>.
   * <p>
   * <b>Example:</b>
   *
   * <pre>
   * 	 // assign x[i] = x[i]&lt;sup&gt;y[i]&lt;/sup&gt;
   * 	 m1 = 0 1 2 3;
   * 	 m2 = 0 2 4 6;
   * 	 m1.assign(m2, cern.jet.math.Functions.pow);
   * 	 --&gt;
   * 	 m1 == 1 1 16 729
   *
   * </pre>
   *
   * For further examples, see the <a
   * href="package-summary.html#FunctionObjects">package doc</a>.
   *
   * @param y
   *            the secondary matrix to operate on.
   * @param function
   *            a function object taking as first argument the current cell's
   *            value of <tt>this</tt>, and as second argument the current
   *            cell's value of <tt>y</tt>,
   * @return <tt>this</tt> (for convenience only).
   * @throws IllegalArgumentException
   *             if <tt>size() != y.size()</tt>.
   * @see cern.jet.math.tdouble.DoubleFunctions
   */
  DoubleVector forEachVector(DoubleVector y, DoubleDoubleFunction function) {
    checkSize(y);
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
            setQuick(i, function(getQuick(i), y.getQuick(i)));
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int i = 0; i < _size; i++) {
        set(i, function(get(i), y.get(i)));
      }
    //}
    return this;
  }

  /**
   * Assigns the result of a function to each cell;
   * <tt>x[i] = function(x[i],y[i])</tt>. (Iterates downwards from
   * <tt>[size()-1]</tt> to <tt>[0]</tt>).
   * <p>
   * <b>Example:</b>
   *
   * <pre>
   * 	 // assign x[i] = x[i]&lt;sup&gt;y[i]&lt;/sup&gt;
   * 	 m1 = 0 1 2 3;
   * 	 m2 = 0 2 4 6;
   * 	 m1.assign(m2, cern.jet.math.Functions.pow);
   * 	 --&gt;
   * 	 m1 == 1 1 16 729
   *
   * 	 // for non-standard functions there is no shortcut:
   * 	 m1.assign(m2,
   * 	    new DoubleDoubleFunction() {
   * 	       double apply(double x, double y) { return Math.pow(x,y); }
   * 	    }
   * 	 );
   *
   * </pre>
   *
   * For further examples, see the <a
   * href="package-summary.html#FunctionObjects">package doc</a>.
   *
   * @param y
   *            the secondary matrix to operate on.
   * @param function
   *            a function object taking as first argument the current cell's
   *            value of <tt>this</tt>, and as second argument the current
   *            cell's value of <tt>y</tt>.
   * @param nonZeroIndexes
   *            list of indexes of non-zero values
   * @return <tt>this</tt> (for convenience only).
   * @throws IllegalArgumentException
   *             if <tt>size() != y.size()</tt>.
   * @see cern.jet.math.tdouble.DoubleFunctions
   */
  DoubleVector forEachVectorNonZero(DoubleVector y, DoubleDoubleFunction function, /*IntArrayList*/Int32List nonZeroIndexes) {
    checkSize(y);
    Int32List nonZeroElements = nonZeroIndexes;//.elements();

    // specialized for speed
    if (function == func.mult) { // x[i] = x[i] * y[i]
      int j = 0;
      for (int index = nonZeroIndexes.length; --index >= 0; ) {
        int i = nonZeroElements[index];
        for ( ; j < i; j++) set(j, 0.0); // x[i] = 0 for all zeros
        set(i, get(i) * y.get(i)); // x[i] * y[i] for all nonZeros
        j++;
      }
    } else if (function is DoublePlusMultSecond) {
      double multiplicator = (function as DoublePlusMultSecond).multiplicator;
      if (multiplicator == 0) { // x[i] = x[i] + 0*y[i]
        return this;
      } else if (multiplicator == 1) { // x[i] = x[i] + y[i]
        for (int index = nonZeroIndexes.length; --index >= 0; ) {
          int i = nonZeroElements[index];
          set(i, get(i) + y.get(i));
        }
      } else if (multiplicator == -1) { // x[i] = x[i] - y[i]
        for (int index = nonZeroIndexes.length; --index >= 0; ) {
          int i = nonZeroElements[index];
          set(i, get(i) - y.get(i));
        }
      } else { // the general case x[i] = x[i] + mult*y[i]
        for (int index = nonZeroIndexes.length; --index >= 0; ) {
          int i = nonZeroElements[index];
          set(i, get(i) + multiplicator * y.get(i));
        }
      }
    } else { // the general case x[i] = f(x[i],y[i])
      return forEachVector(y, function);
    }
    return this;
  }

  /**
   * Returns the number of cells having non-zero values; ignores tolerance.
   *
   * @return the number of cells having non-zero values.
   */
  int get cardinality {
    int cardinality = 0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      Int32List results = new Int32List(nthreads);
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
          final int firstIdx = j * k;
          final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
          futures[j] = ConcurrencyUtils.submit(() {
              int cardinality = 0;
              for (int i = firstIdx; i < lastIdx; i++) {
                  if (getQuick(i) != 0)
                      cardinality++;
              }
              return cardinality;
          });
      }
      try {
          for (int j = 0; j < nthreads; j++) {
              results[j] = futures[j].get() as int;
          }
          Future.wait(futures).then((List results) {
          });
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
    for (int i = 0; i < _size; i++) {
      if (get(i) != 0) cardinality++;
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
  DoubleVector copy() {
    DoubleVector copy = like();
    copy.setAll(this);
    return copy;
  }

  /**
   * Returns the elements of this matrix.
   *
   * @return the elements
   */
  Object elements();

  /**
   * Returns whether all cells are equal to the given value.
   *
   * @param value
   *            the value to test against.
   * @return <tt>true</tt> if all cells are equal to the given value,
   *         <tt>false</tt> otherwise.
   */
  /*bool equalsValue(double value) {
    return DoubleProperty.DEFAULT.equals(this, value);
  }*/

  /**
   * Compares this object against the specified object. The result is
   * <code>true</code> if and only if the argument is not <code>null</code>
   * and is at least a <code>DoubleVector</code> object that has the same
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
      return DoubleProperty.DEFAULT.equals(this, obj.toDouble());      
    }
    if (identical(this, obj)) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    if (obj is! DoubleVector) {
      return false;
    }

    return DoubleProperty.DEFAULT.equalsVector(this, obj as DoubleVector);
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
  double operator [](int index) {
    if (index < 0 || index >= _size) {
      _checkIndex(index);
    }
    return get(index);
  }

  /**
   * Return the maximum value of this matrix together with its location
   *
   * @return { maximum_value, location };
   */
  DoubleVectorLocation max() {
    int location = 0;
    double maxValue = 0.0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      List<List<double>> results = new List<List<double>>(nthreads);//[2];
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int location = firstIdx;
          double maxValue = getQuick(location);
          double elem;
          for (int i = firstIdx + 1; i < lastIdx; i++) {
            elem = getQuick(i);
            if (maxValue < elem) {
              maxValue = elem;
              location = i;
            }
          }
          return [maxValue, location];
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          results[j] = futures[j].get() as List<double>;
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
      maxValue = get(location);
      double elem;
      for (int i = 1; i < length; i++) {
        elem = get(i);
        if (maxValue < elem) {
          maxValue = elem;
          location = i;
        }
      }
    //}
    return new DoubleVectorLocation(maxValue, location);
  }

  /**
   * Return the minimum value of this matrix together with its location
   *
   * @return { minimum_value, location };
   */
  DoubleVectorLocation min() {
    int location = 0;
    double minValue = 0.0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, _size);
      List<Future> futures = new List<Future>(nthreads);
      List<List<double>> results = new List<List<double>>(nthreads);//[2];
      int k = _size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int location = firstIdx;
          double minValue = getQuick(location);
          double elem;
          for (int i = firstIdx + 1; i < lastIdx; i++) {
            elem = getQuick(i);
            if (minValue > elem) {
              minValue = elem;
              location = i;
            }
          }
          return [minValue, location];
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          results[j] = futures[j].get() as List<double>;
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
      minValue = get(location);
      double elem;
      for (int i = 1; i < length; i++) {
        elem = get(i);
        if (minValue > elem) {
          minValue = elem;
          location = i;
        }
      }
    //}
    return new DoubleVectorLocation(minValue, location);
  }

  /**
   * Fills the coordinates and values of cells having negative values into the
   * specified lists. Fills into the lists, starting at index 0. After this
   * call returns the specified lists all have a new size, the number of
   * non-zero values.
   *
   * @param indexList
   *            the list to be filled with indexes, can have any size.
   * @param valueList
   *            the list to be filled with values, can have any size.
   */
  void negativeValues(final /*IntArrayList*/List<int> indexList, final /*DoubleArrayList*/List<double> valueList) {
    bool fillIndexList = indexList != null;
    bool fillValueList = valueList != null;
    if (fillIndexList) {
      indexList.clear();
    }
    if (fillValueList) {
      valueList.clear();
    }
    int rem = _size % 2;
    if (rem == 1) {
      double value = get(0);
      if (value < 0) {
        if (fillIndexList) {
          indexList.add(0);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
    }

    for (int i = rem; i < _size; i += 2) {
      double value = get(i);
      if (value < 0) {
        if (fillIndexList) {
          indexList.add(i);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      value = get(i + 1);
      if (value < 0) {
        if (fillIndexList) {
          indexList.add(i + 1);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
    }
  }

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
   * <p>
   * <b>Example:</b> <br>
   *
   * <pre>
   * 	 0, 0, 8, 0, 7
   * 	 --&gt;
   * 	 indexList  = (2,4)
   * 	 valueList  = (8,7)
   *
   * </pre>
   *
   * In other words, <tt>get(2)==8, get(4)==7</tt>.
   *
   * @param indexList
   *            the list to be filled with indexes, can have any size.
   * @param valueList
   *            the list to be filled with values, can have any size.
   */
  void nonZeros(final /*IntArrayList*/List<int> indexList, final /*DoubleArrayList*/List<double> valueList) {
    bool fillIndexList = indexList != null;
    bool fillValueList = valueList != null;
    if (fillIndexList) {
      indexList.clear();
    }
    if (fillValueList) {
      valueList.clear();
    }
    int rem = _size % 2;
    if (rem == 1) {
      double value = get(0);
      if (value != 0) {
        if (fillIndexList) {
          indexList.add(0);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
    }

    for (int i = rem; i < _size; i += 2) {
      double value = get(i);
      if (value != 0) {
        if (fillIndexList) {
          indexList.add(i);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      value = get(i + 1);
      if (value != 0) {
        if (fillIndexList) {
          indexList.add(i + 1);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
    }
  }

  /**
   * Fills the coordinates and values of the first <tt>maxCardinality</tt>
   * cells having non-zero values into the specified lists. Fills into the
   * lists, starting at index 0. After this call returns the specified lists
   * all have a new size, the number of non-zero values.
   * <p>
   * In general, fill order is <i>unspecified</i>. This implementation fills
   * like: <tt>for (index = 0..size()-1)  do ... </tt>. However, subclasses
   * are free to us any other order, even an order that may change over time
   * as cell values are changed. (Of course, result lists indexes are
   * guaranteed to correspond to the same cell).
   * <p>
   * <b>Example:</b> <br>
   *
   * <pre>
   * 	 0, 0, 8, 0, 7
   * 	 --&gt;
   * 	 indexList  = (2,4)
   * 	 valueList  = (8,7)
   *
   * </pre>
   *
   * In other words, <tt>get(2)==8, get(4)==7</tt>.
   *
   * @param indexList
   *            the list to be filled with indexes, can have any size.
   * @param valueList
   *            the list to be filled with values, can have any size.
   * @param maxCardinality
   *            maximal cardinality
   */
  void nonZerosCardinality(/*IntArrayList*/List<int> indexList, /*DoubleArrayList*/List<double> valueList, int maxCardinality) {
    bool fillIndexList = indexList != null;
    bool fillValueList = valueList != null;
    if (fillIndexList) {
      indexList.clear();
    }
    if (fillValueList) {
      valueList.clear();
    }
    int s = _size;
    int currentSize = 0;
    for (int i = 0; i < s; i++) {
      double value = get(i);
      if (value != 0) {
        if (fillIndexList) {
          indexList.add(i);
        }
        if (fillValueList) {
          valueList.add(value);
        }
        currentSize++;
      }
      if (currentSize >= maxCardinality) {
        break;
      }
    }
  }

  /**
   * Fills the coordinates and values of cells having positive values into the
   * specified lists. Fills into the lists, starting at index 0. After this
   * call returns the specified lists all have a new size, the number of
   * non-zero values.
   *
   * @param indexList
   *            the list to be filled with indexes, can have any size.
   * @param valueList
   *            the list to be filled with values, can have any size.
   */
  void positiveValues(final /*IntArrayList*/List<int> indexList, final /*DoubleArrayList*/List<double> valueList) {
    bool fillIndexList = indexList != null;
    bool fillValueList = valueList != null;
    if (fillIndexList) indexList.clear();
    if (fillValueList) valueList.clear();
    int rem = _size % 2;
    if (rem == 1) {
      double value = get(0);
      if (value > 0) {
        if (fillIndexList) {
          indexList.add(0);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
    }

    for (int i = rem; i < _size; i += 2) {
      double value = get(i);
      if (value > 0) {
        if (fillIndexList) {
          indexList.add(i);
        }
        if (fillValueList) {
          valueList.add(value);
        }
      }
      value = get(i + 1);
      if (value > 0) {
        if (fillIndexList) {
          indexList.add(i + 1);
        }
        if (fillValueList) {
          valueList.add(value);
        }
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
  double get(int index);

  /**
   * Construct and returns a new empty matrix <i>of the same dynamic type</i>
   * as the receiver, having the same size. For example, if the receiver is an
   * instance of type <tt>DenseDoubleVector</tt> the new matrix must also be
   * of type <tt>DenseDoubleVector</tt>, if the receiver is an instance of
   * type <tt>SparseDoubleVector</tt> the new matrix must also be of type
   * <tt>SparseDoubleVector</tt>, etc. In general, the new matrix should
   * have internal parametrization as similar as possible.
   *
   * @return a new empty matrix of the same dynamic type.
   */
  DoubleVector like() {
    return like1D(_size);
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
  DoubleVector like1D(int size);

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
  DoubleMatrix like2D(int rows, int columns);

  /**
   * Normalizes this matrix, i.e. makes the sum of all elements equal to 1.0
   * If the matrix contains negative elements then all the values are shifted
   * to ensure non-negativity.
   */
  void normalize() {
    double min = this.min().value;
    if (min < 0) {
      forEach(func.subtract(min));
    }
    if (max().value == 0) {
      fill(1.0 / length);
    } else {
      double sumScaleFactor = sum();
      sumScaleFactor = 1.0 / sumScaleFactor;
      forEach(func.multiply(sumScaleFactor));
    }
  }

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
  DoubleMatrix reshape(int rows, int columns);

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
  //DoubleMatrix3D reshape3D(int slices, int rows, int columns);

  /**
   * Sets the matrix cell at coordinate <tt>index</tt> to the specified value.
   *
   * @param index
   *            the index of the cell.
   * @param value
   *            the value to be filled into the specified cell.
   * @throws IndexOutOfBoundsException
   *             if <tt>index&lt;0 || index&gt;=size()</tt>.
   */
  void operator []=(int index, double value) {
    if (index < 0 || index >= _size) _checkIndex(index);
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
   * @param value
   *            the value to be filled into the specified cell.
   */
  void set(int index, double value);

  /**
   * Swaps each element <tt>this[i]</tt> with <tt>other[i]</tt>.
   *
   * @throws IllegalArgumentException
   *             if <tt>size() != other.size()</tt>.
   */
  void swap(final DoubleVector other) {
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
            double tmp = getQuick(i);
            setQuick(i, other.getQuick(i));
            other.setQuick(i, tmp);
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int i = 0; i < _size; i++) {
        double tmp = get(i);
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
   * <tt>for (int i=0; i < size(); i++) values[i] = get(i);</tt>
   *
   * @return an array filled with the values of the cells.
   */
  List<double> toList() {
    Float64List values = new Float64List(_size);
    fillList(values);
    return values;
  }

  /**
   * Fills the cell values into the specified 1-dimensional array. The values
   * are copied. So subsequent changes in <tt>values</tt> are not reflected in
   * the matrix, and vice-versa. After this call returns the array
   * <tt>values</tt> has the form <br>
   * <tt>for (int i=0; i < size(); i++) values[i] = get(i);</tt>
   *
   * @throws IllegalArgumentException
   *             if <tt>values.length < size()</tt>.
   */
  void fillList(List<double> values) {
    if (values.length < _size) throw new ArgumentError("values too small");
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
            values[i] = getQuick(i);
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int i = 0; i < _size; i++) {
        values[i] = get(i);
      }
    //}
  }

  /**
   * Returns a string representation using default formatting.
   *
   * @see cern.colt.matrix.tdouble.algo.DoubleFormatter
   */
  String toString() {
    return new DoubleFormatter().toStringDouble1D(this);
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
  DoubleVector flip() {
    return _view()._vFlip() as DoubleVector;
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
  DoubleVector part(int index, int width) {
    return _view()._vPart(index, width) as DoubleVector;
  }

  /**
   * Constructs and returns a new <i>selection view</i> that is a matrix
   * holding the cells matching the given condition. Applies the condition to
   * each cell and takes only those cells where
   * <tt>condition(get(i))</tt> yields <tt>true</tt>.
   * <p>
   * <b>Example:</b> <br>
   *
   * <pre>
   * 	 // extract and view all cells with even value
   * 	 matrix = 0 1 2 3
   * 	 matrix.viewSelection(
   * 	    new DoubleProcedure() {
   * 	       final bool apply(double a) { return a % 2 == 0; }
   * 	    }
   * 	 );
   * 	 --&gt;
   * 	 matrix ==  0 2
   *
   * </pre>
   *
   * For further examples, see the <a
   * href="package-summary.html#FunctionObjects">package doc</a>. The returned
   * view is backed by this matrix, so changes in the returned view are
   * reflected in this matrix, and vice-versa.
   *
   * @param condition
   *            The condition to be matched.
   * @return the new view.
   */
  DoubleVector where(DoubleProcedure condition) {
    /*IntArrayList*/List<int> matches = new /*IntArrayList*/List<int>();
    for (int i = 0; i < _size; i++) {
      if (condition(get(i))) matches.add(i);
    }
    //matches.trimToSize();
    return select(new Int32List.fromList(matches)/*.elements()*/);
  }

  /**
   * Constructs and returns a new <i>selection view</i> that is a matrix
   * holding the indicated cells. There holds
   * <tt>view.size() == indexes.length</tt> and
   * <tt>view.get(i) == this.get(indexes[i])</tt>. Indexes can occur multiple
   * times and can be in arbitrary order.
   * <p>
   * <b>Example:</b> <br>
   *
   * <pre>
   * 	 this     = (0,0,8,0,7)
   * 	 indexes  = (0,2,4,2)
   * 	 --&gt;
   * 	 view     = (0,8,7,8)
   *
   * </pre>
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
  DoubleVector select(Int32List indexes) {
    // check for "all"
    if (indexes == null) {
      indexes = new Int32List(_size);
      for (int i = 0; i < _size; i++) indexes[i] = i;
    }

    _checkIndexes(indexes);
    Int32List offsets = new Int32List(indexes.length);
    for (int i = 0; i < indexes.length; i++) {
      offsets[i] = index(indexes[i]);
    }
    return _viewSelectionLike(offsets);
  }

  /**
   * Sorts the vector into ascending order, according to the <i>natural
   * ordering</i>. This sort is guaranteed to be <i>stable</i>. For further
   * information, see
   * {@link cern.colt.matrix.tdouble.algo.DoubleSorting#sort(DoubleVector)}.
   * For more advanced sorting functionality, see
   * {@link cern.colt.matrix.tdouble.algo.DoubleSorting}.
   *
   * @return a new sorted vector (matrix) view.
   */
  /*DoubleVector viewSorted() {
    return DoubleSorting.mergeSort.sort(this);
  }*/

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
  DoubleVector strides(int stride) {
    return _view()._vStrides(stride) as DoubleVector;
  }

  /**
   * Returns the dot product of two vectors x and y, which is
   * <tt>Sum(x[i]*y[i])</tt>. Where <tt>x == this</tt>. Operates on cells at
   * indexes <tt>0 .. Math.min(size(),y.size())</tt>.
   *
   * @param y
   *            the second vector.
   * @return the sum of products.
   */
  /*double zDotProduct(DoubleVector y) {
    return zDotProductRange(y, 0, _size);
  }*/

  /**
   * Returns the dot product of two vectors x and y, which is
   * <tt>Sum(x[i]*y[i])</tt>. Where <tt>x == this</tt>. Operates on cells at
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
  double dot(final DoubleVector y, [final int from = 0, int length = null]) {
    if (length == null) {
      length = _size;
    }
    if (from < 0 || length <= 0) {
      return 0.0;
    }

    int tail = from + length;
    if (_size < tail) {
      tail = _size;
    }
    if (y._size < tail) {
      tail = y._size;
    }
    length = tail - from;

    double sum = 0.0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, length);
      List<Future> futures = new List<Future>(nthreads);
      List<double> results = new List<double>(nthreads);
      int k = length ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? length : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          double sum = 0;
          int idx;
          for (int k = firstIdx; k < lastIdx; k++) {
            idx = k + from;
            sum += getQuick(idx) * y.getQuick(idx);
          }
          return sum;
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          results[j] = futures[j].get() as double;
        }
        sum = results[0].doubleValue();
        for (int j = 1; j < nthreads; j++) {
          sum += results[j].doubleValue();
        }
      } on ExecutionException catch (ex) {
        ex.printStackTrace();
      } on InterruptedException catch (e) {
        e.printStackTrace();
      }
    } else {*/
      int i = tail - 1;
      for (int k = length; --k >= 0; i--) {
        sum += get(i) * y.get(i);
      }
    //}
    return sum;
  }

  /**
   * Returns the dot product of two vectors x and y, which is
   * <tt>Sum(x[i]*y[i])</tt>. Where <tt>x == this</tt>.
   *
   * @param y
   *            the second vector.
   * @param nonZeroIndexes
   *            the indexes of cells in <tt>y</tt>having a non-zero value.
   * @return the sum of products.
   */
  double dotNonZero(DoubleVector y,  /*IntArrayList*/Int32List nonZeroIndexes, [int from=0, int length=null]) {
    if (length == null) {
      length = _size;
    }
    // determine minimum length
    if (from < 0 || length <= 0) {
      return 0.0;
    }

    int tail = from + length;
    if (_size < tail) {
      tail = _size;
    }
    if (y._size < tail) {
      tail = y._size;
    }
    length = tail - from;
    if (length <= 0) {
      return 0.0;
    }
    /*IntArrayList*/Int32List indexesCopy = new Int32List.fromList(nonZeroIndexes);
    //indexesCopy.trimToSize();
    indexesCopy.sort();
    Int32List nonZeroIndexElements = indexesCopy;//.elements();
    int index = 0;
    int s = indexesCopy.length;
    // skip to start
    while ((index < s) && nonZeroIndexElements[index] < from) index++;
    // now the sparse dot product
    int i;
    double sum = 0.0;
    while ((--length >= 0) && (index < s) && ((i = nonZeroIndexElements[index]) < tail)) {
      sum += get(i) * y.get(i);
      index++;
    }
    return sum;
  }

  /**
   * Returns the sum of all cells; <tt>Sum( x[i] )</tt>.
   *
   * @return the sum.
   */
  double sum() {
    if (length == 0) {
      return 0.0;
    }
    return reduce(func.plus, func.identity);
  }

  /**
   * Returns the number of cells having non-zero values, but at most
   * maxCardinality; ignores tolerance.
   */
  int _cardinality(int maxCardinality) {
    int cardinality = 0;
    int i = _size;
    while (--i >= 0 && cardinality < maxCardinality) {
      if (get(i) != 0) cardinality++;
    }
    return cardinality;
  }

  /**
   * Returns the content of this matrix if it is a wrapper; or <tt>this</tt>
   * otherwise. Override this method in wrappers.
   */
  DoubleVector _getContent() {
    return this;
  }

  /**
   * Returns <tt>true</tt> if both matrices share at least one identical cell.
   */
  bool _haveSharedCells(DoubleVector other) {
    if (other == null) return false;
    if (this == other) return true;
    return _getContent()._haveSharedCellsRaw(other._getContent());
  }

  /**
   * Returns <tt>true</tt> if both matrices share at least one identical cell.
   */
  bool _haveSharedCellsRaw(DoubleVector other) {
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
  DoubleVector _view() {
    return clone() as DoubleVector;
  }

  /**
   * Construct and returns a new selection view.
   *
   * @param offsets
   *            the offsets of the visible elements.
   * @return a new view.
   */
  DoubleVector _viewSelectionLike(Int32List offsets);

  /**
   * Returns the dot product of two vectors x and y, which is
   * <tt>Sum(x[i]*y[i])</tt>. Where <tt>x == this</tt>.
   *
   * @param y
   *            the second vector.
   * @param nonZeroIndexes
   *            the indexes of cells in <tt>y</tt>having a non-zero value.
   * @return the sum of products.
   */
  /*double _zDotProduct(DoubleVector y,  /*IntArrayList*/Int32List nonZeroIndexes) {
    return dotNonZero(y, 0, _size, nonZeroIndexes);
  }*/

  Object clone();
}
