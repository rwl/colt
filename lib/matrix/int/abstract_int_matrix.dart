/*
Copyright (C) 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
is hereby granted without fee, provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear in supporting documentation.
CERN makes no representations about the suitability of this software for any purpose.
It is provided "as is" without expressed or implied warranty.
 */
part of cern.colt.matrix;

typedef bool IntVectorProcedure(AbstractIntVector element);

class IntMatrixLocation {
  final int value, row, column;
  IntMatrixLocation(this.value, this.row, this.column);
}

/**
 * Abstract base class for 2-d matrices holding <tt>int</tt> elements. First see
 * the <a href="package-summary.html">package summary</a> and javadoc <a
 * href="package-tree.html">tree view</a> to get the broad picture.
 * <p>
 * A matrix has a number of rows and columns, which are assigned upon instance
 * construction - The matrix's size is then <tt>rows()*columns()</tt>. Elements
 * are accessed via <tt>[row,column]</tt> coordinates. Legal coordinates range
 * from <tt>[0,0]</tt> to <tt>[rows()-1,columns()-1]</tt>. Any attempt to access
 * an element at a coordinate
 * <tt>column&lt;0 || column&gt;=columns() || row&lt;0 || row&gt;=rows()</tt>
 * will throw an <tt>IndexOutOfBoundsException</tt>.
 * <p>
 * <b>Note</b> that this implementation is not synchronized.
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 09/24/99
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
abstract class AbstractIntMatrix extends AbstractMatrix {

  /**
   * Applies a function to each cell and aggregates the results. Returns a
   * value <tt>v</tt> such that <tt>v==a(size())</tt> where
   * <tt>a(i) == aggr( a(i-1), f(get(row,column)) )</tt> and terminators are
   * <tt>a(1) == f(get(0,0)), a(0)==Int.NaN</tt>.
   * <p>
   * <b>Example:</b>
   *
   * <pre>
   *   cern.jet.math.Functions F = cern.jet.math.Functions.functions;
   *   2 x 2 matrix
   *   0 1
   *   2 3
   *
   *   // Sum( x[row,col]*x[row,col] )
   *   matrix.aggregate(F.plus,F.square);
   *   --&gt; 14
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
   * @see cern.jet.math.tint.IntFunctions
   */
  int reduce(final ifunc.IntIntFunction aggr, final ifunc.IntFunction f) {
    if (length == 0) {
      throw new ArgumentError("size == 0");
    }
    int a = 0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_rows * _columns >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int a = f(get(firstRow, 0));
          int d = 1;
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = d; c < _columns; c++) {
              a = aggr(a, f(get(r, c)));
            }
            d = 0;
          }
          return a;
        });
      }
      a = ConcurrencyUtils.waitForCompletion(futures, aggr);
    } else {*/
      a = f(get(0, 0));
      int d = 1; // first cell already done
      for (int r = 0; r < _rows; r++) {
        for (int c = d; c < _columns; c++) {
          a = aggr(a, f(get(r, c)));
        }
        d = 0;
      }
    //}
    return a;
  }

  /**
   * Applies a function to each cell that satisfies a condition and aggregates
   * the results.
   *
   * @param aggr
   *            an aggregation function taking as first argument the current
   *            aggregation and as second argument the transformed current
   *            cell value.
   * @param f
   *            a function transforming the current cell value.
   * @param cond
   *            a condition.
   * @return the aggregated measure.
   * @see cern.jet.math.tint.IntFunctions
   */
  int reduceWhere(final ifunc.IntIntFunction aggr, final ifunc.IntFunction f, final ifunc.IntProcedure cond) {
    if (length == 0) {
      throw new ArgumentError("size == 0");
    }
    int a = 0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_rows * _columns >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int elem = get(firstRow, 0);
          int a = 0;
          if (cond(elem) == true) {
            a = aggr(a, f(elem));
          }
          int d = 1;
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = d; c < _columns; c++) {
              elem = get(r, c);
              if (cond(elem) == true) {
                a = aggr(a, f(elem));
              }
            }
            d = 0;
          }
          return a;
        });
      }
      a = ConcurrencyUtils.waitForCompletion(futures, aggr);
    } else {*/
      int elem = get(0, 0);
      if (cond(elem) == true) {
        a = aggr(a, f(elem));
      }
      int d = 1; // first cell already done
      for (int r = 0; r < _rows; r++) {
        for (int c = d; c < _columns; c++) {
          elem = get(r, c);
          if (cond(elem) == true) {
            a = aggr(a, f(elem));
          }
        }
        d = 0;
      }
    //}
    return a;
  }

  /**
   *
   * Applies a function to all cells with a given indexes and aggregates the
   * results.
   *
   * @param aggr
   *            an aggregation function taking as first argument the current
   *            aggregation and as second argument the transformed current
   *            cell value.
   * @param f
   *            a function transforming the current cell value.
   * @param rowList
   *            row indexes.
   * @param columnList
   *            column indexes.
   *
   * @return the aggregated measure.
   * @see cern.jet.math.tint.IntFunctions
   */
  int reduceRange(final ifunc.IntIntFunction aggr, final ifunc.IntFunction f, final Int32List rowList, final Int32List columnList) {
    if (length == 0) {
      throw new ArgumentError("size == 0");
    }
    final int size = rowList.length;
    final Int32List rowElements = rowList;//.elements();
    final Int32List columnElements = columnList;//.elements();
    int a = 0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, size);
      List<Future> futures = new List<Future>(nthreads);
      int k = size / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int a = f(get(rowElements[firstIdx], columnElements[firstIdx]));
          int elem;
          for (int i = firstIdx + 1; i < lastIdx; i++) {
            elem = get(rowElements[i], columnElements[i]);
            a = aggr(a, f(elem));
          }
          return a;
        });
      }
      a = ConcurrencyUtils.waitForCompletion(futures, aggr);
    } else {*/
      int elem;
      a = f(get(rowElements[0], columnElements[0]));
      for (int i = 1; i < size; i++) {
        elem = get(rowElements[i], columnElements[i]);
        a = aggr(a, f(elem));
      }
    //}
    return a;
  }

  /**
   * Applies a function to each corresponding cell of two matrices and
   * aggregates the results. Returns a value <tt>v</tt> such that
   * <tt>v==a(size())</tt> where
   * <tt>a(i) == aggr( a(i-1), f(get(row,column),other.get(row,column)) )</tt>
   * and terminators are
   * <tt>a(1) == f(get(0,0),other.get(0,0)), a(0)==Int.NaN</tt>.
   * <p>
   * <b>Example:</b>
   *
   * <pre>
   *   cern.jet.math.Functions F = cern.jet.math.Functions.functions;
   *   x == 2 x 2 matrix
   *   0 1
   *   2 3
   *
   *   y == 2 x 2 matrix
   *   0 1
   *   2 3
   *
   *   // Sum( x[row,col] * y[row,col] )
   *   x.aggregate(y, F.plus, F.mult);
   *   --&gt; 14
   *
   *   // Sum( (x[row,col] + y[row,col])&circ;2 )
   *   x.aggregate(y, F.plus, F.chain(F.square,F.plus));
   *   --&gt; 56
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
   * @throws ArgumentError
   *             if
   *             <tt>columns() != other.columns() || rows() != other.rows()</tt>
   * @see cern.jet.math.tint.IntFunctions
   */
  int reduceWith(final AbstractIntMatrix other, final ifunc.IntIntFunction aggr, final ifunc.IntIntFunction f) {
    checkShape(other);
    if (length == 0) {
      throw new ArgumentError("size == 0");
    }
    int a = 0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_rows * _columns >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int a = f(get(firstRow, 0), other.get(firstRow, 0));
          int d = 1;
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = d; c < _columns; c++) {
              a = aggr(a, f(get(r, c), other.get(r, c)));
            }
            d = 0;
          }
          return a;
        });
      }
      a = ConcurrencyUtils.waitForCompletion(futures, aggr);
    } else {*/
      a = f(get(0, 0), other.get(0, 0));
      int d = 1; // first cell already done
      for (int r = 0; r < _rows; r++) {
        for (int c = d; c < _columns; c++) {
          a = aggr(a, f(get(r, c), other.get(r, c)));
        }
        d = 0;
      }
    //}
    return a;
  }

  /**
   * Assigns the result of a function to each cell;
   * <tt>x[row,col] = function(x[row,col])</tt>.
   * <p>
   * <b>Example:</b>
   *
   * <pre>
   *   matrix = 2 x 2 matrix
   *   0.5 1.5
   *   2.5 3.5
   *
   *   // change each cell to its sine
   *   matrix.assign(cern.jet.math.Functions.sin);
   *   --&gt;
   *   2 x 2 matrix
   *   0.479426  0.997495
   *   0.598472 -0.350783
   *
   * </pre>
   *
   * For further examples, see the <a
   * href="package-summary.html#FunctionObjects">package doc</a>.
   *
   * @param f
   *            a function object taking as argument the current cell's value.
   * @return <tt>this</tt> (for convenience only).
   * @see cern.jet.math.tint.IntFunctions
   */
  void forEach(final ifunc.IntFunction f) {
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_rows * _columns >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

        futures[j] = ConcurrencyUtils.submit(() {
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = 0; c < _columns; c++) {
              set(r, c, f(get(r, c)));
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int r = 0; r < _rows; r++) {
        for (int c = 0; c < _columns; c++) {
          set(r, c, f(get(r, c)));
        }
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
   * @see cern.jet.math.tint.IntFunctions
   */
  void forEachWhere(final ifunc.IntProcedure cond, final ifunc.IntFunction f) {
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_rows * _columns >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

        futures[j] = ConcurrencyUtils.submit(() {
          int elem;
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = 0; c < _columns; c++) {
              elem = get(r, c);
              if (cond(elem) == true) {
                set(r, c, f(elem));
              }
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      int elem;
      for (int r = 0; r < _rows; r++) {
        for (int c = 0; c < _columns; c++) {
          elem = get(r, c);
          if (cond(elem) == true) {
            set(r, c, f(elem));
          }
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
   *            a value.
   * @return <tt>this</tt> (for convenience only).
   *
   */
  void fillWhere(final ifunc.IntProcedure cond, final int value) {
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_rows * _columns >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

        futures[j] = ConcurrencyUtils.submit(() {
          int elem;
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = 0; c < _columns; c++) {
              elem = get(r, c);
              if (cond(elem) == true) {
                set(r, c, value);
              }
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      int elem;
      for (int r = 0; r < _rows; r++) {
        for (int c = 0; c < _columns; c++) {
          elem = get(r, c);
          if (cond(elem) == true) {
            set(r, c, value);
          }
        }
      }
    //}
  }

  /**
   * Sets all cells to the state specified by <tt>value</tt>.
   *
   * @param value
   *            the value to be filled into the cells.
   * @return <tt>this</tt> (for convenience only).
   */
  void fill(final int value) {
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_rows * _columns >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

        futures[j] = ConcurrencyUtils.submit(() {
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = 0; c < _columns; c++) {
              set(r, c, value);
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int r = 0; r < _rows; r++) {
        for (int c = 0; c < _columns; c++) {
          set(r, c, value);
        }
      }
    //}
  }

  /**
   * Sets all cells to the state specified by <tt>values</tt>. <tt>values</tt>
   * is required to have the form <tt>values[row*column]</tt> and elements
   * have to be stored in a row-wise order.
   * <p>
   * The values are copied. So subsequent changes in <tt>values</tt> are not
   * reflected in the matrix, and vice-versa.
   *
   * @param values
   *            the values to be filled into the cells.
   * @return <tt>this</tt> (for convenience only).
   * @throws ArgumentError
   *             if <tt>values.length != rows()*columns()</tt>.
   */
  void setAll(final Int32List values) {
    if (values.length != _rows * _columns) {
      throw new ArgumentError("Must have same length: length=${values.length} rows()*columns()=${rows * columns}");
    }
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_rows * _columns >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

        futures[j] = ConcurrencyUtils.submit(() {
          int idx = firstRow * _columns;
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = 0; c < _columns; c++) {
              set(r, c, values[idx++]);
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      int idx = 0;
      for (int r = 0; r < _rows; r++) {
        for (int c = 0; c < _columns; c++) {
          set(r, c, values[idx++]);
        }
      }
    //}

  }

  /**
   * Sets all cells to the state specified by <tt>values</tt>. <tt>values</tt>
   * is required to have the form <tt>values[row][column]</tt> and have
   * exactly the same number of rows and columns as the receiver.
   * <p>
   * The values are copied. So subsequent changes in <tt>values</tt> are not
   * reflected in the matrix, and vice-versa.
   *
   * @param values
   *            the values to be filled into the cells.
   * @return <tt>this</tt> (for convenience only).
   * @throws ArgumentError
   *             if
   *             <tt>values.length != rows() || for any 0 &lt;= row &lt; rows(): values[row].length != columns()</tt>
   *             .
   */
  void setAll2D(final List<Int32List> values) {
    if (values.length != _rows) {
      throw new ArgumentError("Must have same number of rows: rows=${values.length} rows()=${rows}");
    }
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_rows * _columns >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

        futures[j] = ConcurrencyUtils.submit(() {
          for (int r = firstRow; r < lastRow; r++) {
            Int32List currentRow = values[r];
            if (currentRow.length != _columns) throw new ArgumentError("Must have same number of columns in every row: columns=" + currentRow.length + "columns()=" + columns());
            for (int c = 0; c < _columns; c++) {
              set(r, c, currentRow[c]);
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int r = 0; r < _rows; r++) {
        Int32List currentRow = values[r];
        if (currentRow.length != _columns) {
          throw new ArgumentError("Must have same number of columns in every row: columns=${currentRow.length} columns()=$columns");
        }
        for (int c = 0; c < _columns; c++) {
          set(r, c, currentRow[c]);
        }
      }
    //}
  }

  /**
   * Replaces all cell values of the receiver with the values of another
   * matrix. Both matrices must have the same number of rows and columns. If
   * both matrices share the same cells (as is the case if they are views
   * derived from the same matrix) and intersect in an ambiguous way, then
   * replaces <i>as if</i> using an intermediate auxiliary deep copy of
   * <tt>other</tt>.
   *
   * @param other
   *            the source matrix to copy from (may be identical to the
   *            receiver).
   * @return <tt>this</tt> (for convenience only).
   * @throws ArgumentError
   *             if
   *             <tt>columns() != other.columns() || rows() != other.rows()</tt>
   */
  void copyFrom(AbstractIntMatrix other) {
    if (other == this) {
      return;
    }
    checkShape(other);
    AbstractIntMatrix other_loc;
    if (_haveSharedCells(other)) {
      other_loc = other.copy();
    } else {
      other_loc = other;
    }
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_rows * _columns >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

        futures[j] = ConcurrencyUtils.submit(() {
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = 0; c < _columns; c++) {
              set(r, c, other_loc.get(r, c));
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int r = 0; r < _rows; r++) {
        for (int c = 0; c < _columns; c++) {
          set(r, c, other_loc.get(r, c));
        }
      }
    //}
  }

  /**
   * Assigns the result of a function to each cell;
   * <tt>x[row,col] = function(x[row,col],y[row,col])</tt>.
   * <p>
   * <b>Example:</b>
   *
   * <pre>
   *   // assign x[row,col] = x[row,col]&lt;sup&gt;y[row,col]&lt;/sup&gt;
   *   m1 = 2 x 2 matrix
   *   0 1
   *   2 3
   *
   *   m2 = 2 x 2 matrix
   *   0 2
   *   4 6
   *
   *   m1.assign(m2, cern.jet.math.Functions.pow);
   *   --&gt;
   *   m1 == 2 x 2 matrix
   *   1   1
   *   16 729
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
   * @throws ArgumentError
   *             if
   *             <tt>columns() != other.columns() || rows() != other.rows()</tt>
   * @see cern.jet.math.tint.IntFunctions
   */
  void forEachWith(final AbstractIntMatrix y, final ifunc.IntIntFunction function) {
    checkShape(y);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_rows * _columns >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

        futures[j] = ConcurrencyUtils.submit(() {
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = 0; c < _columns; c++) {
              set(r, c, function(get(r, c), y.get(r, c)));
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int r = 0; r < _rows; r++) {
        for (int c = 0; c < _columns; c++) {
          set(r, c, function(get(r, c), y.get(r, c)));
        }
      }
    //}
  }

  /**
   * Assigns the result of a function to all cells with a given indexes
   *
   * @param y
   *            the secondary matrix to operate on.
   * @param function
   *            a function object taking as first argument the current cell's
   *            value of <tt>this</tt>, and as second argument the current
   *            cell's value of <tt>y</tt>,
   * @param rowList
   *            row indexes.
   * @param columnList
   *            column indexes.
   *
   * @return <tt>this</tt> (for convenience only).
   * @throws ArgumentError
   *             if
   *             <tt>columns() != other.columns() || rows() != other.rows()</tt>
   * @see cern.jet.math.tint.IntFunctions
   */
  void forEachWithRange(final AbstractIntMatrix y, final ifunc.IntIntFunction function, Int32List rowList, Int32List columnList) {
    checkShape(y);
    final int size = rowList.length;
    final Int32List rowElements = rowList;//.elements();
    final Int32List columnElements = columnList;//.elements();
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, size);
      List<Future> futures = new List<Future>(nthreads);
      int k = size / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int i = firstIdx; i < lastIdx; i++) {
            set(rowElements[i], columnElements[i], function(get(rowElements[i], columnElements[i]), y.get(rowElements[i], columnElements[i])));
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int i = 0; i < size; i++) {
        set(rowElements[i], columnElements[i], function(get(rowElements[i], columnElements[i]), y.get(rowElements[i], columnElements[i])));
      }
    //}
  }

  /**
   * Returns the number of cells having non-zero values; ignores tolerance.
   *
   * @return cardinality
   */
  int get cardinality {
    int cardinality = 0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_rows * _columns >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      List<int> results = new List<int>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

        futures[j] = ConcurrencyUtils.submit(() {
          int cardinality = 0;
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = 0; c < _columns; c++) {
              if (get(r, c) != 0) cardinality++;
            }
          }
          return cardinality;
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          results[j] = futures[j].get() as int;
        }
        cardinality = results[0].intValue();
        for (int j = 1; j < nthreads; j++) {
          cardinality += results[j].intValue();
        }
      } on ExecutionException catch (ex) {
        ex.printStackTrace();
      } on InterruptedException catch (e) {
        e.printStackTrace();
      }
    } else {*/
      for (int r = 0; r < _rows; r++) {
        for (int c = 0; c < _columns; c++) {
          if (get(r, c) != 0) cardinality++;
        }
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
  AbstractIntMatrix copy() {
    return like()..copyFrom(this);
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
  /*bool equals(int value) {
    return IntProperty.DEFAULT.equalsMatrixValue(this, value);
  }*/

  /**
   * Compares this object against the specified object. The result is
   * <code>true</code> if and only if the argument is not <code>null</code>
   * and is at least a <code>IntMatrix</code> object that has the same
   * number of columns and rows as the receiver and has exactly the same
   * values at the same coordinates.
   *
   * @param obj
   *            the object to compare with.
   * @return <code>true</code> if the objects are the same; <code>false</code>
   *         otherwise.
   */

  bool operator ==(var obj) {
    if (obj is int) {
      return iprop.equalsMatrixValue(this, obj);
    }
    if (identical(this, obj)) {
      return true;
    }
    if (obj == null) {
      return false;
    }
    if (!(obj is AbstractIntMatrix)) {
      return false;
    }

    return iprop.equalsMatrix(this, obj as AbstractIntMatrix);
  }

  /**
   * Assigns the result of a function to each <i>non-zero</i> cell;
   * <tt>x[row,col] = function(x[row,col])</tt>. Use this method for fast
   * special-purpose iteration. If you want to modify another matrix instead
   * of <tt>this</tt> (i.e. work in read-only mode), simply return the input
   * value unchanged.
   *
   * Parameters to function are as follows: <tt>first==row</tt>,
   * <tt>second==column</tt>, <tt>third==nonZeroValue</tt>.
   *
   * @param function
   *            a function object taking as argument the current non-zero
   *            cell's row, column and value.
   * @return <tt>this</tt> (for convenience only).
   */
  void forEachNonZero(final ifunc.IntIntIntFunction function) {
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

        futures[j] = ConcurrencyUtils.submit(() {
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = 0; c < _columns; c++) {
              int value = get(r, c);
              if (value != 0) {
                int a = function(r, c, value);
                if (a != value) set(r, c, a);
              }
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int r = 0; r < _rows; r++) {
        for (int c = 0; c < _columns; c++) {
          int value = get(r, c);
          if (value != 0) {
            int a = function(r, c, value);
            if (a != value) set(r, c, a);
          }
        }
      }
    //}
  }

  /**
   * Returns the matrix cell value at coordinate <tt>[row,column]</tt>.
   *
   * @param row
   *            the index of the row-coordinate.
   * @param column
   *            the index of the column-coordinate.
   * @return the value of the specified cell.
   * @throws IndexOutOfBoundsException
   *             if
   *             <tt>column&lt;0 || column&gt;=columns() || row&lt;0 || row&gt;=rows()</tt>
   */
  int at(int row, int column) {
    if (column < 0 || column >= _columns || row < 0 || row >= _rows) {
      throw new RangeError("row:$row, column:$column");
    }
    return get(row, column);
  }

  /**
   * Returns the content of this matrix if it is a wrapper; or <tt>this</tt>
   * otherwise. Override this method in wrappers.
   */
  AbstractIntMatrix _getContent() {
    return this;
  }

  /**
   * Fills the coordinates and values of cells having negative values into the
   * specified lists. Fills into the lists, starting at index 0. After this
   * call returns the specified lists all have a new size, the number of
   * non-zero values.
   *
   * @param rowList
   *            the list to be filled with row indexes, can have any size.
   * @param columnList
   *            the list to be filled with column indexes, can have any size.
   * @param valueList
   *            the list to be filled with values, can have any size.
   */
  void negativeValues(final Int32List rowList, final Int32List columnList, final Int32List valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    for (int r = 0; r < _rows; r++) {
      for (int c = 0; c < _columns; c++) {
        int value = get(r, c);
        if (value < 0) {
          rowList.add(r);
          columnList.add(c);
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
   * like <tt>for (row = 0..rows-1) for (column = 0..columns-1) do ... </tt>.
   * However, subclasses are free to us any other order, even an order that
   * may change over time as cell values are changed. (Of course, result lists
   * indexes are guaranteed to correspond to the same cell).
   * <p>
   * <b>Example:</b> <br>
   *
   * <pre>
   *   2 x 3 matrix:
   *   0, 0, 8
   *   0, 7, 0
   *   --&gt;
   *   rowList    = (0,1)
   *   columnList = (2,1)
   *   valueList  = (8,7)
   *
   * </pre>
   *
   * In other words, <tt>get(0,2)==8, get(1,1)==7</tt>.
   *
   * @param rowList
   *            the list to be filled with row indexes, can have any size.
   * @param columnList
   *            the list to be filled with column indexes, can have any size.
   * @param valueList
   *            the list to be filled with values, can have any size.
   */
  void nonZeros(final Int32List rowList, final Int32List columnList, final Int32List valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    for (int r = 0; r < _rows; r++) {
      for (int c = 0; c < _columns; c++) {
        int value = get(r, c);
        if (value != 0) {
          rowList.add(r);
          columnList.add(c);
          valueList.add(value);
        }
      }
    }
  }

  /**
   * Fills the coordinates and values of cells having positive values into the
   * specified lists. Fills into the lists, starting at index 0. After this
   * call returns the specified lists all have a new size, the number of
   * non-zero values.
   *
   * @param rowList
   *            the list to be filled with row indexes, can have any size.
   * @param columnList
   *            the list to be filled with column indexes, can have any size.
   * @param valueList
   *            the list to be filled with values, can have any size.
   */
  void positiveValues(final Int32List rowList, final Int32List columnList, final Int32List valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    for (int r = 0; r < _rows; r++) {
      for (int c = 0; c < _columns; c++) {
        int value = get(r, c);
        if (value > 0) {
          rowList.add(r);
          columnList.add(c);
          valueList.add(value);
        }
      }
    }
  }

  /**
   * Returns the matrix cell value at coordinate <tt>[row,column]</tt>.
   *
   * <p>
   * Provided with invalid parameters this method may return invalid objects
   * without throwing any exception. <b>You should only use this method when
   * you are absolutely sure that the coordinate is within bounds.</b>
   * Precondition (unchecked):
   * <tt>0 &lt;= column &lt; columns() && 0 &lt;= row &lt; rows()</tt>.
   *
   * @param row
   *            the index of the row-coordinate.
   * @param column
   *            the index of the column-coordinate.
   * @return the value at the specified coordinate.
   */
  int get(int row, int column);

  /**
   * Returns <tt>true</tt> if both matrices share at least one identical cell.
   */
  bool _haveSharedCells(AbstractIntMatrix other) {
    if (other == null) return false;
    if (this == other) return true;
    return _getContent()._haveSharedCellsRaw(other._getContent());
  }

  /**
   * Returns <tt>true</tt> if both matrices share at least one identical cell.
   */
  bool _haveSharedCellsRaw(AbstractIntMatrix other) {
    return false;
  }

  /**
   * Construct and returns a new empty matrix <i>of the same dynamic type</i>
   * as the receiver, having the same number of rows and columns. For example,
   * if the receiver is an instance of type <tt>DenseIntMatrix</tt> the new
   * matrix must also be of type <tt>DenseIntMatrix</tt>, if the receiver is
   * an instance of type <tt>SparseIntMatrix</tt> the new matrix must also
   * be of type <tt>SparseIntMatrix</tt>, etc. In general, the new matrix
   * should have internal parametrization as similar as possible.
   *
   * @return a new empty matrix of the same dynamic type.
   */
  AbstractIntMatrix like() {
    return like2D(_rows, _columns);
  }

  /**
   * Construct and returns a new empty matrix <i>of the same dynamic type</i>
   * as the receiver, having the specified number of rows and columns. For
   * example, if the receiver is an instance of type <tt>DenseIntMatrix</tt>
   * the new matrix must also be of type <tt>DenseIntMatrix</tt>, if the
   * receiver is an instance of type <tt>SparseIntMatrix</tt> the new matrix
   * must also be of type <tt>SparseIntMatrix</tt>, etc. In general, the new
   * matrix should have internal parametrization as similar as possible.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @return a new empty matrix of the same dynamic type.
   */
  AbstractIntMatrix like2D(int rows, int columns);

  /**
   * Construct and returns a new 1-d matrix <i>of the corresponding dynamic
   * type</i>, entirelly independent of the receiver. For example, if the
   * receiver is an instance of type <tt>DenseIntMatrix</tt> the new matrix
   * must be of type <tt>DenseIntVector</tt>, if the receiver is an instance
   * of type <tt>SparseIntMatrix</tt> the new matrix must be of type
   * <tt>SparseIntVector</tt>, etc.
   *
   * @param size
   *            the number of cells the matrix shall have.
   * @return a new matrix of the corresponding dynamic type.
   */
  AbstractIntVector like1D(int size);

  /**
   * Construct and returns a new 1-d matrix <i>of the corresponding dynamic
   * type</i>, sharing the same cells. For example, if the receiver is an
   * instance of type <tt>DenseIntMatrix</tt> the new matrix must be of type
   * <tt>DenseIntVector</tt>, if the receiver is an instance of type
   * <tt>SparseIntMatrix</tt> the new matrix must be of type
   * <tt>SparseIntVector</tt>, etc.
   *
   * @param size
   *            the number of cells the matrix shall have.
   * @param zero
   *            the index of the first element.
   * @param stride
   *            the number of indexes between any two elements, i.e.
   *            <tt>index(i+1)-index(i)</tt>.
   * @return a new matrix of the corresponding dynamic type.
   */
  AbstractIntVector _like1D(int size, int zero, int stride);

  /**
   * Return the maximum value of this matrix together with its location
   *
   * @return maximum_value, row_location, column_location };
   */
  IntMatrixLocation max() {
    int rowLocation = 0;
    int columnLocation = 0;
    int maxValue = 0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      List<Int32List> results = new List<Int32List>(nthreads);//[2];
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

        futures[j] = ConcurrencyUtils.submit(() {
          int rowLocation = firstRow;
          int columnLocation = 0;
          int maxValue = get(rowLocation, 0);
          int d = 1;
          int elem;
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = d; c < _columns; c++) {
              elem = get(r, c);
              if (maxValue < elem) {
                maxValue = elem;
                rowLocation = r;
                columnLocation = c;
              }
            }
            d = 0;
          }
          return [maxValue, rowLocation, columnLocation];
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          results[j] = futures[j].get() as Int32List;
        }
        maxValue = results[0][0];
        rowLocation = results[0][1] as int;
        columnLocation = results[0][2] as int;
        for (int j = 1; j < nthreads; j++) {
          if (maxValue < results[j][0]) {
            maxValue = results[j][0];
            rowLocation = results[j][1] as int;
            columnLocation = results[j][2] as int;
          }
        }
      } on ExecutionException catch (ex) {
        ex.printStackTrace();
      } on InterruptedException catch (e) {
        e.printStackTrace();
      }
    } else {*/
      maxValue = get(0, 0);
      int elem;
      int d = 1;
      for (int r = 0; r < _rows; r++) {
        for (int c = d; c < _columns; c++) {
          elem = get(r, c);
          if (maxValue < elem) {
            maxValue = elem;
            rowLocation = r;
            columnLocation = c;
          }
        }
        d = 0;
      }
    //}
    return new IntMatrixLocation(maxValue, rowLocation, columnLocation);
  }

  /**
   * Return the minimum value of this matrix together with its location
   *
   * @return minimum_value, row_location, column_location};
   */
  IntMatrixLocation min() {
    int rowLocation = 0;
    int columnLocation = 0;
    int minValue = 0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      List<Int32List> results = new List<Int32List>(nthreads);//[2];
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

        futures[j] = ConcurrencyUtils.submit(() {
          int rowLocation = firstRow;
          int columnLocation = 0;
          int minValue = get(rowLocation, 0);
          int d = 1;
          int elem;
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = d; c < _columns; c++) {
              elem = get(r, c);
              if (minValue > elem) {
                minValue = elem;
                rowLocation = r;
                columnLocation = c;
              }
            }
            d = 0;
          }
          return [minValue, rowLocation, columnLocation];
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          results[j] = futures[j].get() as Int32List;
        }
        minValue = results[0][0];
        rowLocation = results[0][1] as int;
        columnLocation = results[0][2] as int;
        for (int j = 1; j < nthreads; j++) {
          if (minValue > results[j][0]) {
            minValue = results[j][0];
            rowLocation = results[j][1] as int;
            columnLocation = results[j][2] as int;
          }
        }
      } on ExecutionException catch (ex) {
        ex.printStackTrace();
      } on InterruptedException catch (e) {
        e.printStackTrace();
      }
    } else {*/
      minValue = get(0, 0);
      int elem;
      int d = 1;
      for (int r = 0; r < _rows; r++) {
        for (int c = d; c < _columns; c++) {
          elem = get(r, c);
          if (minValue > elem) {
            minValue = elem;
            rowLocation = r;
            columnLocation = c;
          }
        }
        d = 0;
      }
    //}
    return new IntMatrixLocation(minValue, rowLocation, columnLocation);
  }

  /**
   * Sets the matrix cell at coordinate <tt>[row,column]</tt> to the specified
   * value.
   *
   * @param row
   *            the index of the row-coordinate.
   * @param column
   *            the index of the column-coordinate.
   * @param value
   *            the value to be filled into the specified cell.
   * @throws IndexOutOfBoundsException
   *             if
   *             <tt>column&lt;0 || column&gt;=columns() || row&lt;0 || row&gt;=rows()</tt>
   */
  void put(int row, int column, int value) {
    if (column < 0 || column >= _columns || row < 0 || row >= _rows) {
      throw new RangeError("row:$row, column:$column");
    }
    set(row, column, value);
  }

  /**
   * Sets the matrix cell at coordinate <tt>[row,column]</tt> to the specified
   * value.
   *
   * <p>
   * Provided with invalid parameters this method may access illegal indexes
   * without throwing any exception. <b>You should only use this method when
   * you are absolutely sure that the coordinate is within bounds.</b>
   * Precondition (unchecked):
   * <tt>0 &lt;= column &lt; columns() && 0 &lt;= row &lt; rows()</tt>.
   *
   * @param row
   *            the index of the row-coordinate.
   * @param column
   *            the index of the column-coordinate.
   * @param value
   *            the value to be filled into the specified cell.
   */
  void set(int row, int column, int value);

  /**
   * Constructs and returns a 2-dimensional array containing the cell values.
   * The returned array <tt>values</tt> has the form
   * <tt>values[row][column]</tt> and has the same number of rows and columns
   * as the receiver.
   * <p>
   * The values are copied. So subsequent changes in <tt>values</tt> are not
   * reflected in the matrix, and vice-versa.
   *
   * @return an array filled with the values of the cells.
   */
  List<Int32List> toList() {
    final List<Int32List> values = new List<Int32List>.generate(_rows, (_) => new Int32List(_columns));
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

        futures[j] = ConcurrencyUtils.submit(() {
          for (int r = firstRow; r < lastRow; r++) {
            Int32List currentRow = values[r];
            for (int c = 0; c < _columns; c++) {
              currentRow[c] = get(r, c);
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int r = 0; r < _rows; r++) {
        Int32List currentRow = values[r];
        for (int c = 0; c < _columns; c++) {
          currentRow[c] = get(r, c);
        }
      }
    //}
    return values;
  }

  /**
   * Returns a string representation using default formatting.
   *
   * @see IntFormatter
   */

  String toString() {
    return new IntFormatter().toString2D(this);
  }

  /**
   * Returns a vector obtained by stacking the columns of the matrix on top of
   * one another.
   *
   * @return a vector of columns of this matrix.
   */
  AbstractIntVector vectorize();

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
  AbstractIntMatrix _view() {
    return clone() as AbstractIntMatrix;
  }

  /**
   * Constructs and returns a new <i>slice view</i> representing the rows of
   * the given column. The returned view is backed by this matrix, so changes
   * in the returned view are reflected in this matrix, and vice-versa. To
   * obtain a slice view on subranges, construct a sub-ranging view (
   * <tt>viewPart(...)</tt>), then apply this method to the sub-range view.
   * <p>
   * <b>Example:</b>
   * <table border="0">
   * <tr nowrap>
   * <td valign="top">2 x 3 matrix: <br>
   * 1, 2, 3<br>
   * 4, 5, 6</td>
   * <td>viewColumn(0) ==></td>
   * <td valign="top">Vector of size 2:<br>
   * 1, 4</td>
   * </tr>
   * </table>
   *
   * @param column
   *            the column to fix.
   * @return a new slice view.
   * @throws IndexOutOfBoundsException
   *             if <tt>column < 0 || column >= columns()</tt>.
   * @see #row(int)
   */
  AbstractIntVector column(int column) {
    _checkColumn(column);
    int viewSize = this._rows;
    int viewZero = index(0, column);
    int viewStride = this._rowStride;
    return _like1D(viewSize, viewZero, viewStride);
  }

  /**
   * Constructs and returns a new <i>flip view</i> aint the column axis. What
   * used to be column <tt>0</tt> is now column <tt>columns()-1</tt>, ...,
   * what used to be column <tt>columns()-1</tt> is now column <tt>0</tt>. The
   * returned view is backed by this matrix, so changes in the returned view
   * are reflected in this matrix, and vice-versa.
   * <p>
   * <b>Example:</b>
   * <table border="0">
   * <tr nowrap>
   * <td valign="top">2 x 3 matrix: <br>
   * 1, 2, 3<br>
   * 4, 5, 6</td>
   * <td>columnFlip ==></td>
   * <td valign="top">2 x 3 matrix:<br>
   * 3, 2, 1 <br>
   * 6, 5, 4</td>
   * <td>columnFlip ==></td>
   * <td valign="top">2 x 3 matrix: <br>
   * 1, 2, 3<br>
   * 4, 5, 6</td>
   * </tr>
   * </table>
   *
   * @return a new flip view.
   * @see #rowFlip()
   */
  AbstractIntMatrix columnFlip() {
    return (_view().._vColumnFlip());
  }

  /**
   * Constructs and returns a new <i>dice (transposition) view</i>; Swaps
   * axes; example: 3 x 4 matrix --> 4 x 3 matrix. The view has both
   * dimensions exchanged; what used to be columns become rows, what used to
   * be rows become columns. In other words:
   * <tt>view.get(row,column)==this.get(column,row)</tt>. This is a zero-copy
   * transposition, taking O(1), i.e. constant time. The returned view is
   * backed by this matrix, so changes in the returned view are reflected in
   * this matrix, and vice-versa. Use idioms like
   * <tt>result = viewDice(A).copy()</tt> to generate an independent
   * transposed matrix.
   * <p>
   * <b>Example:</b>
   * <table border="0">
   * <tr nowrap>
   * <td valign="top">2 x 3 matrix: <br>
   * 1, 2, 3<br>
   * 4, 5, 6</td>
   * <td>transpose ==></td>
   * <td valign="top">3 x 2 matrix:<br>
   * 1, 4 <br>
   * 2, 5 <br>
   * 3, 6</td>
   * <td>transpose ==></td>
   * <td valign="top">2 x 3 matrix: <br>
   * 1, 2, 3<br>
   * 4, 5, 6</td>
   * </tr>
   * </table>
   *
   * @return a new dice view.
   */
  AbstractIntMatrix dice() {
    return _view().._vDice();
  }

  /**
   * Constructs and returns a new <i>sub-range view</i> that is a
   * <tt>height x width</tt> sub matrix starting at <tt>[row,column]</tt>.
   *
   * Operations on the returned view can only be applied to the restricted
   * range. Any attempt to access coordinates not contained in the view will
   * throw an <tt>IndexOutOfBoundsException</tt>.
   * <p>
   * <b>Note that the view is really just a range restriction:</b> The
   * returned matrix is backed by this matrix, so changes in the returned
   * matrix are reflected in this matrix, and vice-versa.
   * <p>
   * The view contains the cells from <tt>[row,column]</tt> to
   * <tt>[row+height-1,column+width-1]</tt>, all inclusive. and has
   * <tt>view.rows() == height; view.columns() == width;</tt>. A view's legal
   * coordinates are again zero based, as usual. In other words, legal
   * coordinates of the view range from <tt>[0,0]</tt> to
   * <tt>[view.rows()-1==height-1,view.columns()-1==width-1]</tt>. As usual,
   * any attempt to access a cell at a coordinate
   * <tt>column&lt;0 || column&gt;=view.columns() || row&lt;0 || row&gt;=view.rows()</tt>
   * will throw an <tt>IndexOutOfBoundsException</tt>.
   *
   * @param row
   *            The index of the row-coordinate.
   * @param column
   *            The index of the column-coordinate.
   * @param height
   *            The height of the box.
   * @param width
   *            The width of the box.
   * @throws IndexOutOfBoundsException
   *             if
   *             <tt>column<0 || width<0 || column+width>columns() || row<0 || height<0 || row+height>rows()</tt>
   * @return the new view.
   *
   */
  AbstractIntMatrix part(int row, int column, int height, int width) {
    return _view().._vPart(row, column, height, width);
  }

  /**
   * Constructs and returns a new <i>slice view</i> representing the columns
   * of the given row. The returned view is backed by this matrix, so changes
   * in the returned view are reflected in this matrix, and vice-versa. To
   * obtain a slice view on subranges, construct a sub-ranging view (
   * <tt>viewPart(...)</tt>), then apply this method to the sub-range view.
   * <p>
   * <b>Example:</b>
   * <table border="0">
   * <tr nowrap>
   * <td valign="top">2 x 3 matrix: <br>
   * 1, 2, 3<br>
   * 4, 5, 6</td>
   * <td>viewRow(0) ==></td>
   * <td valign="top">Vector of size 3:<br>
   * 1, 2, 3</td>
   * </tr>
   * </table>
   *
   * @param row
   *            the row to fix.
   * @return a new slice view.
   * @throws IndexOutOfBoundsException
   *             if <tt>row < 0 || row >= rows()</tt>.
   * @see #column(int)
   */
  AbstractIntVector row(int row) {
    _checkRow(row);
    int viewSize = this._columns;
    int viewZero = index(row, 0);
    int viewStride = this._columnStride;
    return _like1D(viewSize, viewZero, viewStride);
  }

  /**
   * Constructs and returns a new <i>flip view</i> aint the row axis. What
   * used to be row <tt>0</tt> is now row <tt>rows()-1</tt>, ..., what used to
   * be row <tt>rows()-1</tt> is now row <tt>0</tt>. The returned view is
   * backed by this matrix, so changes in the returned view are reflected in
   * this matrix, and vice-versa.
   * <p>
   * <b>Example:</b>
   * <table border="0">
   * <tr nowrap>
   * <td valign="top">2 x 3 matrix: <br>
   * 1, 2, 3<br>
   * 4, 5, 6</td>
   * <td>rowFlip ==></td>
   * <td valign="top">2 x 3 matrix:<br>
   * 4, 5, 6 <br>
   * 1, 2, 3</td>
   * <td>rowFlip ==></td>
   * <td valign="top">2 x 3 matrix: <br>
   * 1, 2, 3<br>
   * 4, 5, 6</td>
   * </tr>
   * </table>
   *
   * @return a new flip view.
   * @see #columnFlip()
   */
  AbstractIntMatrix rowFlip() {
    return _view().._vRowFlip();
  }

  /**
   * Constructs and returns a new <i>selection view</i> that is a matrix
   * holding all <b>rows</b> matching the given condition. Applies the
   * condition to each row and takes only those row where
   * <tt>condition(viewRow(i))</tt> yields <tt>true</tt>. To match
   * columns, use a dice view.
   * <p>
   * <b>Example:</b> <br>
   *
   * <pre>
   *   // extract and view all rows which have a value &lt; threshold in the first column (representing &quot;age&quot;)
   *   final int threshold = 16;
   *   matrix.viewSelection(
   *      new IntVectorProcedure() {
   *         final bool apply(IntVector m) { return m.get(0) &lt; threshold; }
   *      }
   *   );
   *
   *   // extract and view all rows with RMS &lt; threshold
   *   // The RMS (Root-Mean-Square) is a measure of the average &quot;size&quot; of the elements of a data sequence.
   *   matrix = 0 1 2 3
   *   final int threshold = 0.5;
   *   matrix.viewSelection(
   *      new IntVectorProcedure() {
   *         final bool apply(IntVector m) { return Math.sqrt(m.aggregate(F.plus,F.square) / m.length) &lt; threshold; }
   *      }
   *   );
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
  AbstractIntMatrix where(IntVectorProcedure condition) {
    List<int> matches = new List<int>();
    for (int i = 0; i < _rows; i++) {
      if (condition(row(i))) {
        matches.add(i);
      }
    }

    //matches.trimToSize();
    return select(new Int32List.fromList(matches), null); // take all columns
  }

  /**
   * Constructs and returns a new <i>selection view</i> that is a matrix
   * holding the indicated cells. There holds
   * <tt>view.rows() == rowIndexes.length, view.columns() == columnIndexes.length</tt>
   * and <tt>view.get(i,j) == this.get(rowIndexes[i],columnIndexes[j])</tt>.
   * Indexes can occur multiple times and can be in arbitrary order.
   * <p>
   * <b>Example:</b>
   *
   * <pre>
   *   this = 2 x 3 matrix:
   *   1, 2, 3
   *   4, 5, 6
   *   rowIndexes     = (0,1)
   *   columnIndexes  = (1,0,1,0)
   *   --&gt;
   *   view = 2 x 4 matrix:
   *   2, 1, 2, 1
   *   5, 4, 5, 4
   *
   * </pre>
   *
   * Note that modifying the index arguments after this call has returned has
   * no effect on the view. The returned view is backed by this matrix, so
   * changes in the returned view are reflected in this matrix, and
   * vice-versa.
   * <p>
   * To indicate "all" rows or "all columns", simply set the respective
   * parameter
   *
   * @param rowIndexes
   *            The rows of the cells that shall be visible in the new view.
   *            To indicate that <i>all</i> rows shall be visible, simply set
   *            this parameter to <tt>null</tt>.
   * @param columnIndexes
   *            The columns of the cells that shall be visible in the new
   *            view. To indicate that <i>all</i> columns shall be visible,
   *            simply set this parameter to <tt>null</tt>.
   * @return the new view.
   * @throws IndexOutOfBoundsException
   *             if <tt>!(0 <= rowIndexes[i] < rows())</tt> for any
   *             <tt>i=0..rowIndexes.length()-1</tt>.
   * @throws IndexOutOfBoundsException
   *             if <tt>!(0 <= columnIndexes[i] < columns())</tt> for any
   *             <tt>i=0..columnIndexes.length()-1</tt>.
   */
  AbstractIntMatrix select(Int32List rowIndexes, Int32List columnIndexes) {
    // check for "all"
    if (rowIndexes == null) {
      rowIndexes = new Int32List(_rows);
      for (int i = 0; i < _rows; i++) {
        rowIndexes[i] = i;
      }
    }
    if (columnIndexes == null) {
      columnIndexes = new Int32List(_columns);
      for (int i = 0; i < _columns; i++) {
        columnIndexes[i] = i;
      }
    }

    _checkRowIndexes(rowIndexes);
    _checkColumnIndexes(columnIndexes);
    Int32List rowOffsets = new Int32List(rowIndexes.length);
    Int32List columnOffsets = new Int32List(columnIndexes.length);
    for (int i = 0; i < rowIndexes.length; i++) {
      rowOffsets[i] = _rowOffset(_rowRank(rowIndexes[i]));
    }
    for (int i = 0; i < columnIndexes.length; i++) {
      columnOffsets[i] = _columnOffset(_columnRank(columnIndexes[i]));
    }
    return _viewSelectionLike(rowOffsets, columnOffsets);
  }

  AbstractIntMatrix selectIterable(Set<Int32List> indexes) {
    int n = indexes.length;
    Int32List rowIndexes = new Int32List(n);
    Int32List columnIndexes = new Int32List(n);
    int idx = 0;
    for (Iterator<Int32List> iterator = indexes.iterator; iterator.current != null; ) {
      iterator.moveNext();
      Int32List is_ = iterator.current;
      rowIndexes[idx] = is_[0];
      columnIndexes[idx] = is_[1];
      idx++;
    }
    _checkRowIndexes(rowIndexes);
    _checkColumnIndexes(columnIndexes);
    Int32List rowOffsets = new Int32List(rowIndexes.length);
    Int32List columnOffsets = new Int32List(columnIndexes.length);
    for (int i = 0; i < rowIndexes.length; i++) {
      rowOffsets[i] = _rowOffset(_rowRank(rowIndexes[i]));
    }
    for (int i = 0; i < columnIndexes.length; i++) {
      columnOffsets[i] = _columnOffset(_columnRank(columnIndexes[i]));
    }
    return _viewSelectionLike(rowOffsets, columnOffsets);
  }

  /**
   * Construct and returns a new selection view.
   *
   * @param rowOffsets
   *            the offsets of the visible elements.
   * @param columnOffsets
   *            the offsets of the visible elements.
   * @return a new view.
   */
  AbstractIntMatrix _viewSelectionLike(Int32List rowOffsets, Int32List columnOffsets);

  /**
   * Sorts the matrix rows into ascending order, according to the <i>natural
   * ordering</i> of the matrix values in the given column. This sort is
   * guaranteed to be <i>stable</i>. For further information, see
   * {@link IntSorting#sort(IntMatrix,int)}. For
   * more advanced sorting functionality, see
   * {@link IntSorting}.
   *
   * @return a new sorted vector (matrix) view.
   * @throws IndexOutOfBoundsException
   *             if <tt>column < 0 || column >= columns()</tt>.
   */
  /*IntMatrix sorted(int column) {
    return IntSorting.mergeSort.sort(this, column);
  }*/

  /**
   * Constructs and returns a new <i>stride view</i> which is a sub matrix
   * consisting of every i-th cell. More specifically, the view has
   * <tt>this.rows()/rowStride</tt> rows and
   * <tt>this.columns()/columnStride</tt> columns holding cells
   * <tt>this.get(i*rowStride,j*columnStride)</tt> for all
   * <tt>i = 0..rows()/rowStride - 1, j = 0..columns()/columnStride - 1</tt>.
   * The returned view is backed by this matrix, so changes in the returned
   * view are reflected in this matrix, and vice-versa.
   *
   * @param rowStride
   *            the row step factor.
   * @param columnStride
   *            the column step factor.
   * @return a new view.
   * @throws IndexOutOfBoundsException
   *             if <tt>rowStride<=0 || columnStride<=0</tt>.
   */
  AbstractIntMatrix strides(int rowStride, int columnStride) {
    return _view().._vStrides(rowStride, columnStride);
  }

  /**
   * Linear algebraic matrix-vector multiplication; <tt>z = A * y</tt>;
   * Equivalent to <tt>return A.zMult(y,z,1,0);</tt>
   */
  /*IntVector mult(IntVector y, IntVector z) {
      return mult(y, z, 1, (z == null ? 1 : 0), false);
  }*/

  /**
   * Linear algebraic matrix-vector multiplication;
   * <tt>z = alpha * A * y + beta*z</tt>.
   * <tt>z[i] = alpha*Sum(A[i,j] * y[j]) + beta*z[i], i=0..A.rows()-1, j=0..y.length-1</tt>
   * . Where <tt>A == this</tt>. <br>
   * Note: Matrix shape conformance is checked <i>after</i> potential
   * transpositions.
   *
   * @param y
   *            the source vector.
   * @param z
   *            the vector where results are to be stored. Set this parameter
   *            to <tt>null</tt> to indicate that a new result vector shall be
   *            constructed.
   * @return z (for convenience only).
   *
   * @throws ArgumentError
   *             if <tt>A.columns() != y.length || A.rows() > z.length)</tt>.
   */
  AbstractIntVector mult(final AbstractIntVector y, [AbstractIntVector z = null, final int alpha = 1, int beta = null, final bool transposeA = false]) {
    if (beta == null) {
      beta = z == null ? 1 : 0;
    }
    if (transposeA) {
      return dice().mult(y, z, alpha, beta, false);
    }
    AbstractIntVector z_loc;
    if (z == null) {
      z_loc = new IntVector(this._rows);
    } else {
      z_loc = z;
    }
    if (_columns != y.length || _rows > z_loc.length) {
      throw new ArgumentError("Incompatible args: " + toStringShort() + ", " + y.toStringShort() + ", " + z_loc.toStringShort());
    }

    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

        futures[j] = ConcurrencyUtils.submit(() {
          for (int r = firstRow; r < lastRow; r++) {
            int s = 0;
            for (int c = 0; c < _columns; c++) {
              s += get(r, c) * y.get(c);
            }
            z_loc.set(r, alpha * s + beta * z_loc.get(r));
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int r = 0; r < _rows; r++) {
        int s = 0;
        for (int c = 0; c < _columns; c++) {
          s += get(r, c) * y.get(c);
        }
        z_loc.set(r, alpha * s + beta * z_loc.get(r));
      }
    //}
    return z_loc;
  }

  /**
   * Linear algebraic matrix-matrix multiplication; <tt>C = A x B</tt>;
   * Equivalent to <tt>A.zMult(B,C,1,0,false,false)</tt>.
   */
  /*IntMatrix multiply(IntMatrix B, IntMatrix C) {
        return multiply(B, C, 1, (C == null ? 1 : 0), false, false);
    }*/

  /**
   * Linear algebraic matrix-matrix multiplication;
   * <tt>C = alpha * A x B + beta*C</tt>.
   * <tt>C[i,j] = alpha*Sum(A[i,k] * B[k,j]) + beta*C[i,j], k=0..n-1</tt>. <br>
   * Matrix shapes: <tt>A(m x n), B(n x p), C(m x p)</tt>. <br>
   * Note: Matrix shape conformance is checked <i>after</i> potential
   * transpositions.
   *
   * @param B
   *            the second source matrix.
   * @param C
   *            the matrix where results are to be stored. Set this parameter
   *            to <tt>null</tt> to indicate that a new result matrix shall be
   *            constructed.
   * @return C (for convenience only).
   *
   * @throws ArgumentError
   *             if <tt>B.rows() != A.columns()</tt>.
   * @throws ArgumentError
   *             if
   *             <tt>C.rows() != A.rows() || C.columns() != B.columns()</tt>.
   * @throws ArgumentError
   *             if <tt>A == C || B == C</tt>.
   */
  AbstractIntMatrix multiply(final AbstractIntMatrix B, [AbstractIntMatrix C = null, final int alpha = 1, int beta = null, final bool transposeA = false, final bool transposeB = false]) {
    if (beta == null) {
      beta = C == null ? 1 : 0;
    }
    if (transposeA) {
      return dice().multiply(B, C, alpha, beta, false, transposeB);
    }
    if (transposeB) {
      return this.multiply(B.dice(), C, alpha, beta, transposeA, false);
    }

    final int m = _rows;
    final int n = _columns;
    final int p = B._columns;
    AbstractIntMatrix C_loc;
    if (C == null) {
      C_loc = new IntMatrix(m, p);
    } else {
      C_loc = C;
    }
    if (B._rows != n) {
      throw new ArgumentError("Matrix inner dimensions must agree:" + toStringShort() + ", " + B.toStringShort());
    }
    if (C_loc._rows != m || C_loc._columns != p) {
      throw new ArgumentError("Incompatibe result matrix: " + toStringShort() + ", " + B.toStringShort() + ", " + C_loc.toStringShort());
    }
    if (this == C_loc || B == C_loc) {
      throw new ArgumentError("Matrices must not be identical");
    }
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, p);
      List<Future> futures = new List<Future>(nthreads);
      int k = p / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? p : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int a = firstIdx; a < lastIdx; a++) {
            for (int b = 0; b < m; b++) {
              int s = 0;
              for (int c = 0; c < n; c++) {
                s += get(b, c) * B.get(c, a);
              }
              C_loc.set(b, a, alpha * s + beta * C_loc.get(b, a));
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      for (int a = 0; a < p; a++) {
        for (int b = 0; b < m; b++) {
          int s = 0;
          for (int c = 0; c < n; c++) {
            s += get(b, c) * B.get(c, a);
          }
          C_loc.set(b, a, alpha * s + beta * C_loc.get(b, a));
        }
      }
    //}
    return C_loc;
  }

  /**
   * Returns the sum of all cells; <tt>Sum( x[i,j] )</tt>.
   *
   * @return the sum.
   */
  int sum() {
    if (length == 0) {
      return 0;
    }
    return reduce(ifunc.plus, ifunc.identity);
  }

  Object clone();
}
