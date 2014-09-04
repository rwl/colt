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
 * Dense 2-d matrix holding <tt>double</tt> elements. First see the <a
 * href="package-summary.html">package summary</a> and javadoc <a
 * href="package-tree.html">tree view</a> to get the broad picture.
 * <p>
 * <b>Implementation:</b>
 * <p>
 * Internally holds one single contigous one-dimensional array, addressed in row
 * major. Note that this implementation is not synchronized.
 * <p>
 * <b>Time complexity:</b>
 * <p>
 * <tt>O(1)</tt> (i.e. constant time) for the basic operations <tt>get</tt>,
 * <tt>getQuick</tt>, <tt>set</tt>, <tt>setQuick</tt> and <tt>size</tt>,
 * <p>
 * Cells are internally addressed in row-major. Applications demanding utmost
 * speed can exploit this fact. Setting/getting values in a loop row-by-row is
 * quicker than column-by-column. Thus
 *
 * <pre>
 * for (int row = 0; row &lt; rows; row++) {
 *     for (int column = 0; column &lt; columns; column++) {
 *         matrix.setQuick(row, column, someValue);
 *     }
 * }
 * </pre>
 *
 * is quicker than
 *
 * <pre>
 * for (int column = 0; column &lt; columns; column++) {
 *     for (int row = 0; row &lt; rows; row++) {
 *         matrix.setQuick(row, column, someValue);
 *     }
 * }
 * </pre>
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 09/24/99
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
class DenseDoubleMatrix2D extends DoubleMatrix2D {

  Float64List _elements;

  /**
   * Constructs a matrix with a copy of the given values. <tt>values</tt> is
   * required to have the form <tt>values[row][column]</tt> and have exactly
   * the same number of columns in every row.
   * <p>
   * The values are copied. So subsequent changes in <tt>values</tt> are not
   * reflected in the matrix, and vice-versa.
   *
   * @param values
   *            The values to be filled into the new matrix.
   * @throws ArgumentError
   *             if
   *             <tt>for any 1 &lt;= row &lt; values.length: values[row].length != values[row-1].length</tt>
   *             .
   */
  factory DenseDoubleMatrix2D.from(List<Float64List> values) {
    return new DenseDoubleMatrix2D(values.length, values.length == 0 ? 0 : values[0].length)
      ..setAll2D(values);
  }

  /**
   * Constructs a matrix with a given number of rows and columns. All entries
   * are initially <tt>0</tt>.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @throws ArgumentError
   *             if
   *             <tt>rows<0 || columns<0 || (double)columns*rows > Integer.MAX_VALUE</tt>
   *             .
   */
  //    DenseDoubleMatrix2D(int rows, int columns) {
  //        _setUp(rows, columns);
  //        this._elements = new Float64List(rows * columns);
  //    }

  /**
   * Constructs a matrix with the given parameters.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @param elements
   *            the cells.
   * @param rowZero
   *            the position of the first element.
   * @param columnZero
   *            the position of the first element.
   * @param rowStride
   *            the number of elements between two rows, i.e.
   *            <tt>index(i+1,j)-index(i,j)</tt>.
   * @param columnStride
   *            the number of elements between two columns, i.e.
   *            <tt>index(i,j+1)-index(i,j)</tt>.
   * @param isView
   *            if true then a matrix view is constructed
   * @throws ArgumentError
   *             if
   *             <tt>rows<0 || columns<0 || (double)columns*rows > Integer.MAX_VALUE</tt>
   *             or flip's are illegal.
   */
  DenseDoubleMatrix2D(int rows, int columns, [Float64List elements = null, int rowZero = 0, int columnZero = 0, int rowStride = null, int columnStride = 1, bool isView = false]) {
    if (elements == null) {
      elements = new Float64List(rows * columns);
    }
    _setUp(rows, columns, rowZero, columnZero, rowStride, columnStride);
    this._elements = elements;
    this._isNoView = !isView;
  }

  /**
   * Constructs a matrix from MatrixVectorReader.
   *
   * @param reader
   *            matrix reader
   * @throws IOException
   */
  /*factory DenseDoubleMatrix2D.read(MatrixVectorReader reader) {//throws IOException {
    MatrixInfo info;
    if (reader.hasInfo()) {
      info = reader.readMatrixInfo();
    } else {
      info = new MatrixInfo(true, MatrixInfo.MatrixField.Real, MatrixInfo.MatrixSymmetry.General);
    }

    if (info.isPattern()) throw new UnsupportedOperationException("Pattern matrices are not supported");
    if (info.isDense()) throw new UnsupportedOperationException("Dense matrices are not supported");
    if (info.isComplex()) throw new UnsupportedOperationException("Complex matrices are not supported");

    MatrixSize size = reader.readMatrixSize(info);
    final m = new DenseDoubleMatrix2D(size.numRows(), size.numColumns());
    m._elements = new Float64List(_rows * _columns);
    int numEntries = size.numEntries();
    List<int> columnIndexes = new List<int>(numEntries);
    List<int> rowIndexes = new List<int>(numEntries);
    Float64List values = new Float64List(numEntries);
    reader.readCoordinate(rowIndexes, columnIndexes, values);
    for (int i = 0; i < numEntries; i++) {
      m.setQuick(rowIndexes[i], columnIndexes[i], values[i]);
    }
    if (info.isSymmetric()) {
      for (int i = 0; i < numEntries; i++) {
        if (rowIndexes[i] != columnIndexes[i]) {
          m.setQuick(columnIndexes[i], rowIndexes[i], values[i]);
        }
      }
    } else if (info.isSkewSymmetric()) {
      for (int i = 0; i < numEntries; i++) {
        if (rowIndexes[i] != columnIndexes[i]) {
          m.setQuick(columnIndexes[i], rowIndexes[i], -values[i]);
        }
      }
    }
    return m;
  }*/

  double reduce(DoubleDoubleFunction aggr, DoubleFunction f) {
    if (length == 0) return double.NAN;
    final int zero = index(0, 0);
    double a = 0.0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = _rows - j * k;
        final int lastRow = (j == (nthreads - 1)) ? 0 : firstRow - k;
        futures[j] = ConcurrencyUtils.submit(() {
          double a = f(_elements[zero + (firstRow - 1) * _rowStride + (_columns - 1) * _columnStride]);
          int d = 1;
          for (int r = firstRow; --r >= lastRow; ) {
            int ridx = zero + r * _rowStride;
            for (int c = _columns - d; --c >= 0; ) {
              a = aggr(a, f(_elements[ridx + c * _columnStride]));
            }
            d = 0;
          }
          return a;
        });
      }
      a = ConcurrencyUtils.waitForCompletion(futures, aggr);
    } else {*/
    a = f(_elements[zero + (_rows - 1) * _rowStride + (_columns - 1) * _columnStride]);
    int d = 1;
    for (int r = _rows; --r >= 0; ) {
      int ridx = zero + r * _rowStride;
      for (int c = _columns - d; --c >= 0; ) {
        a = aggr(a, f(_elements[ridx + c * _columnStride]));
      }
      d = 0;
    }
    //}
    return a;
  }

  double reduceWhere(DoubleDoubleFunction aggr, DoubleFunction f, DoubleProcedure cond) {
    if (length == 0) return double.NAN;
    final int zero = index(0, 0);
    double a = 0.0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          double elem = _elements[zero + firstRow * _rowStride];
          double a = 0;
          if (cond(elem) == true) {
            a = f(elem);
          }
          int d = 1;
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = d; c < _columns; c++) {
              elem = _elements[zero + r * _rowStride + c * _columnStride];
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
    double elem = _elements[zero];
    if (cond(elem) == true) {
      a = f(_elements[zero]);
    }
    int d = 1; // first cell already done
    for (int r = 0; r < _rows; r++) {
      for (int c = d; c < _columns; c++) {
        elem = _elements[zero + r * _rowStride + c * _columnStride];
        if (cond(elem) == true) {
          a = aggr(a, f(elem));
        }
      }
      d = 0;
    }
    //}
    return a;
  }

  double reduceRange(DoubleDoubleFunction aggr, DoubleFunction f, final /*IntArrayList*/List<int> rowList, final /*IntArrayList*/List<int> columnList) {
    if (this.length == 0) return double.NAN;
    final int zero = index(0, 0);
    final int size = rowList.length;
    final List<int> rowElements = rowList;//.elements();
    final List<int> columnElements = columnList;//.elements();
    double a = 0.0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, size);
      List<Future> futures = new List<Future>(nthreads);
      int k = size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          double a = f(_elements[zero + rowElements[firstIdx] * _rowStride + columnElements[firstIdx] * _columnStride]);
          double elem;
          for (int i = firstIdx + 1; i < lastIdx; i++) {
            elem = _elements[zero + rowElements[i] * _rowStride + columnElements[i] * _columnStride];
            a = aggr(a, f(elem));
          }
          return a;
        });
      }
      a = ConcurrencyUtils.waitForCompletion(futures, aggr);
    } else {*/
    double elem;
    a = f(_elements[zero + rowElements[0] * _rowStride + columnElements[0] * _columnStride]);
    for (int i = 1; i < size; i++) {
      elem = _elements[zero + rowElements[i] * _rowStride + columnElements[i] * _columnStride];
      a = aggr(a, f(elem));
    }
    //}
    return a;
  }

  double reduceMatrix(final DoubleMatrix2D other, DoubleDoubleFunction aggr, DoubleDoubleFunction f) {
    if (!(other is DenseDoubleMatrix2D)) {
      return super.reduceMatrix(other, aggr, f);
    }
    checkShape(other);
    if (length == 0) return double.NAN;
    final int zero = index(0, 0);
    final int zeroOther = other.index(0, 0);
    final int rowStrideOther = other.rowStride;
    final int colStrideOther = other.columnStride;
    final Float64List elementsOther = other.elements() as Float64List;
    double a = 0.0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          double a = f(_elements[zero + firstRow * _rowStride], elementsOther[zeroOther + firstRow * rowStrideOther]);
          int d = 1;
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = d; c < _columns; c++) {
              a = aggr(a, f(_elements[zero + r * _rowStride + c * _columnStride], elementsOther[zeroOther + r * rowStrideOther + c * colStrideOther]));
            }
            d = 0;
          }
          return a;
        });
      }
      a = ConcurrencyUtils.waitForCompletion(futures, aggr);
    } else {*/
    int d = 1; // first cell already done
    a = f(_elements[zero], elementsOther[zeroOther]);
    for (int r = 0; r < _rows; r++) {
      for (int c = d; c < _columns; c++) {
        a = aggr(a, f(_elements[zero + r * _rowStride + c * _columnStride], elementsOther[zeroOther + r * rowStrideOther + c * colStrideOther]));
      }
      d = 0;
    }
    //}
    return a;
  }

  DoubleMatrix2D forEach(DoubleFunction function) {
    final Float64List elems = this._elements;
    if (elems == null) throw new Error();
    final int zero = index(0, 0);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      if (function is DoubleMult) { // x[i] =
        // mult*x[i]
        double multiplicator = function.multiplicator;
        if (multiplicator == 1) return this;
        if (multiplicator == 0) return assign(0);
      }
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = zero + firstRow * _rowStride;
          // specialization for speed
          if (function is DoubleMult) {
            // x[i] = mult*x[i]
            double multiplicator = function.multiplicator;
            if (multiplicator == 1) return;
            for (int r = firstRow; r < lastRow; r++) {
              for (int i = idx,
                  c = 0; c < _columns; c++) {
                elems[i] *= multiplicator;
                i += _columnStride;
              }
              idx += _rowStride;
            }
          } else {
            // the general case x[i] = f(x[i])
            for (int r = firstRow; r < lastRow; r++) {
              for (int i = idx,
                  c = 0; c < _columns; c++) {
                elems[i] = function(elems[i]);
                i += _columnStride;
              }
              idx += _rowStride;
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idx = zero + (_rows - 1) * _rowStride + (_columns - 1) * _columnStride;
    // specialization for speed
    if (function is DoubleMult) { // x[i] =
      // mult*x[i]
      double multiplicator = (function as DoubleMult).multiplicator;
      if (multiplicator == 1) return this;
      if (multiplicator == 0) return fill(0.0);
      for (int r = _rows; --r >= 0; ) { // the general case
        for (int i = idx,
            c = _columns; --c >= 0; ) {
          elems[i] *= multiplicator;
          i -= _columnStride;
        }
        idx -= _rowStride;
      }
    } else { // the general case x[i] = f(x[i])
      for (int r = _rows; --r >= 0; ) {
        for (int i = idx,
            c = _columns; --c >= 0; ) {
          elems[i] = function(elems[i]);
          i -= _columnStride;
        }
        idx -= _rowStride;
      }
    }
    //}
    return this;
  }

  DoubleMatrix2D forEachWhere(DoubleProcedure cond, func.DoubleFunction function) {
    final int zero = index(0, 0);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          double elem;
          int idx = zero + firstRow * _rowStride;
          for (int r = firstRow; r < lastRow; r++) {
            for (int i = idx,
                c = 0; c < _columns; c++) {
              elem = _elements[i];
              if (cond(elem) == true) {
                _elements[i] = function(elem);
              }
              i += _columnStride;
            }
            idx += _rowStride;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    double elem;
    int idx = zero;
    for (int r = 0; r < _rows; r++) {
      for (int i = idx,
          c = 0; c < _columns; c++) {
        elem = _elements[i];
        if (cond(elem) == true) {
          _elements[i] = function(elem);
        }
        i += _columnStride;
      }
      idx += _rowStride;
    }
    //}
    return this;
  }

  DoubleMatrix2D fillWhere(DoubleProcedure cond, final double value) {
    final int zero = index(0, 0);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          double elem;
          int idx = zero + firstRow * _rowStride;
          for (int r = firstRow; r < lastRow; r++) {
            for (int i = idx,
                c = 0; c < _columns; c++) {
              elem = _elements[i];
              if (cond(elem) == true) {
                _elements[i] = value;
              }
              i += _columnStride;
            }
            idx += _rowStride;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    double elem;
    int idx = zero;
    for (int r = 0; r < _rows; r++) {
      for (int i = idx,
          c = 0; c < _columns; c++) {
        elem = _elements[i];
        if (cond(elem) == true) {
          _elements[i] = value;
        }
        i += _columnStride;
      }
      idx += _rowStride;
    }
    //}
    return this;
  }

  DoubleMatrix2D fill(final double value) {
    final Float64List elems = this._elements;
    final int zero = index(0, 0);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = zero + firstRow * _rowStride;
          for (int r = firstRow; r < lastRow; r++) {
            for (int i = idx,
                c = 0; c < _columns; c++) {
              elems[i] = value;
              i += _columnStride;
            }
            idx += _rowStride;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idx = zero;
    for (int r = 0; r < _rows; r++) {
      for (int i = idx,
          c = 0; c < _columns; c++) {
        elems[i] = value;
        i += _columnStride;
      }
      idx += _rowStride;
    }
    //}
    return this;
  }

  DoubleMatrix2D setAll(final Float64List values) {
    if (values.length != length) {
      throw new ArgumentError("Must have same length: length=${values.length} rows()*columns()=${rows * columns}");
    }
    //int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if (this._isNoView) {
      this._elements.setAll(0, values);
      //System.arraycopy(values, 0, this._elements, 0, values.length);
    } else {
      final int zero = index(0, 0);
      /*if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
        nthreads = Math.min(nthreads, _rows);
        List<Future> futures = new List<Future>(nthreads);
        int k = _rows ~/ nthreads;
        for (int j = 0; j < nthreads; j++) {
          final int firstRow = j * k;
          final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
          futures[j] = ConcurrencyUtils.submit(() {
            int idxOther = firstRow * _columns;
            int idx = zero + firstRow * _rowStride;
            for (int r = firstRow; r < lastRow; r++) {
              for (int i = idx,
                  c = 0; c < _columns; c++) {
                _elements[i] = values[idxOther++];
                i += _columnStride;
              }
              idx += _rowStride;
            }
          });
        }
        ConcurrencyUtils.waitForCompletion(futures);
      } else {*/
      int idxOther = 0;
      int idx = zero;
      for (int r = 0; r < _rows; r++) {
        for (int i = idx,
            c = 0; c < _columns; c++) {
          _elements[i] = values[idxOther++];
          i += _columnStride;
        }
        idx += _rowStride;
      }
      //}
    }
    return this;
  }

  DoubleMatrix2D setAll2D(final List<Float64List> values) {
    if (values.length != _rows) {
      throw new ArgumentError("Must have same number of rows: rows=${values.length} rows()=${rows}");
    }
    //int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if (this._isNoView) {
      /*if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
        nthreads = Math.min(nthreads, _rows);
        List<Future> futures = new List<Future>(nthreads);
        int k = _rows ~/ nthreads;
        for (int j = 0; j < nthreads; j++) {
          final int firstRow = j * k;
          final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
          futures[j] = ConcurrencyUtils.submit(() {
            int i = firstRow * _rowStride;
            for (int r = firstRow; r < lastRow; r++) {
              Float64List currentRow = values[r];
              if (currentRow.length != _columns) throw new ArgumentError("Must have same number of columns in every row: columns=" + currentRow.length + "columns()=" + columns());
              System.arraycopy(currentRow, 0, _elements, i, _columns);
              i += _columns;
            }
          });
        }
        ConcurrencyUtils.waitForCompletion(futures);
      } else {*/
      int i = 0;
      for (int r = 0; r < _rows; r++) {
        Float64List currentRow = values[r];
        if (currentRow.length != _columns) {
          throw new ArgumentError("Must have same number of columns in every row: columns=${currentRow.length} columns()=${columns}");
        }
        this._elements.setAll(i, currentRow);
        //System.arraycopy(currentRow, 0, this._elements, i, _columns);
        i += _columns;
      }
      //}
    } else {
      final int zero = this.index(0, 0);
      /*if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
        nthreads = Math.min(nthreads, _rows);
        List<Future> futures = new List<Future>(nthreads);
        int k = _rows ~/ nthreads;
        for (int j = 0; j < nthreads; j++) {
          final int firstRow = j * k;
          final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
          futures[j] = ConcurrencyUtils.submit(() {
            int idx = zero + firstRow * _rowStride;
            for (int r = firstRow; r < lastRow; r++) {
              Float64List currentRow = values[r];
              if (currentRow.length != _columns) throw new ArgumentError("Must have same number of columns in every row: columns=" + currentRow.length + "columns()=" + columns());
              for (int i = idx,
                  c = 0; c < _columns; c++) {
                _elements[i] = currentRow[c];
                i += _columnStride;
              }
              idx += _rowStride;
            }
          });
        }
        ConcurrencyUtils.waitForCompletion(futures);
      } else {*/
      int idx = zero;
      for (int r = 0; r < _rows; r++) {
        Float64List currentRow = values[r];
        if (currentRow.length != _columns) {
          throw new ArgumentError("Must have same number of columns in every row: columns=${currentRow.length} columns()=${columns}");
        }
        for (int i = idx,
            c = 0; c < _columns; c++) {
          _elements[i] = currentRow[c];
          i += _columnStride;
        }
        idx += _rowStride;
      }
      //}
      return this;
    }
    return this;
  }

  DoubleMatrix2D copyFrom(final DoubleMatrix2D source) {
    // overriden for performance only
    if (!(source is DenseDoubleMatrix2D)) {
      super.copyFrom(source);
      return this;
    }
    DenseDoubleMatrix2D other = source as DenseDoubleMatrix2D;
    if (other == this) return this; // nothing to do
    checkShape(other);
    //int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if (this._isNoView && other._isNoView) { // quickest
      this._elements.setAll(0, other._elements);
      //System.arraycopy(other._elements, 0, this._elements, 0, this._elements.length);
      return this;
    }
    if (_haveSharedCells(other)) {
      DoubleMatrix2D c = other.copy();
      if (!(c is DenseDoubleMatrix2D)) { // should not happen
        super.copyFrom(other);
        return this;
      }
      other = c as DenseDoubleMatrix2D;
    }

    final Float64List elementsOther = other._elements;
    if (_elements == null || elementsOther == null) {
      throw new Error();
    }
    final int zeroOther = other.index(0, 0);
    final int zero = this.index(0, 0);
    final int columnStrideOther = other._columnStride;
    final int rowStrideOther = other._rowStride;
    /*if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = zero + firstRow * _rowStride;
          int idxOther = zeroOther + firstRow * rowStrideOther;
          for (int r = firstRow; r < lastRow; r++) {
            for (int i = idx,
                j = idxOther,
                c = 0; c < _columns; c++) {
              _elements[i] = elementsOther[j];
              i += _columnStride;
              j += columnStrideOther;
            }
            idx += _rowStride;
            idxOther += rowStrideOther;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idx = zero;
    int idxOther = zeroOther;
    for (int r = 0; r < _rows; r++) {
      for (int i = idx,
          j = idxOther,
          c = 0; c < _columns; c++) {
        _elements[i] = elementsOther[j];
        i += _columnStride;
        j += columnStrideOther;
      }
      idx += _rowStride;
      idxOther += rowStrideOther;
    }
    //}
    return this;
  }

  DoubleMatrix2D forEachMatrix(final DoubleMatrix2D y, DoubleDoubleFunction function) {
    // overriden for performance only
    if (!(y is DenseDoubleMatrix2D)) {
      super.forEachMatrix(y, function);
      return this;
    }
    DenseDoubleMatrix2D other = y as DenseDoubleMatrix2D;
    checkShape(y);
    final Float64List elementsOther = other._elements;
    if (_elements == null || elementsOther == null) throw new Error();
    final int zeroOther = other.index(0, 0);
    final int zero = this.index(0, 0);
    final int columnStrideOther = other._columnStride;
    final int rowStrideOther = other._rowStride;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      if (function is DoublePlusMultSecond) {
        double multiplicator = function.multiplicator;
        if (multiplicator == 0) { // x[i] = x[i] + 0*y[i]
          return this;
        }
      }
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx;
          int idxOther;
          // specialized for speed
          if (function == cern.jet.math.tdouble.DoubleFunctions.mult) {
            // x[i] = x[i]*y[i]
            idx = zero + firstRow * _rowStride;
            idxOther = zeroOther + firstRow * rowStrideOther;
            for (int r = firstRow; r < lastRow; r++) {
              for (int i = idx,
                  j = idxOther,
                  c = 0; c < _columns; c++) {
                _elements[i] *= elementsOther[j];
                i += _columnStride;
                j += columnStrideOther;
              }
              idx += _rowStride;
              idxOther += rowStrideOther;
            }
          } else if (function == cern.jet.math.tdouble.DoubleFunctions.div) {
            // x[i] = x[i] / y[i]
            idx = zero + firstRow * _rowStride;
            idxOther = zeroOther + firstRow * rowStrideOther;
            for (int r = firstRow; r < lastRow; r++) {
              for (int i = idx,
                  j = idxOther,
                  c = 0; c < _columns; c++) {
                _elements[i] /= elementsOther[j];
                i += _columnStride;
                j += columnStrideOther;
              }
              idx += _rowStride;
              idxOther += rowStrideOther;
            }
          } else if (function is DoublePlusMultSecond) {
            double multiplicator = function.multiplicator;
            if (multiplicator == 1) {
              // x[i] = x[i] + y[i]
              idx = zero + firstRow * _rowStride;
              idxOther = zeroOther + firstRow * rowStrideOther;
              for (int r = firstRow; r < lastRow; r++) {
                for (int i = idx,
                    j = idxOther,
                    c = 0; c < _columns; c++) {
                  _elements[i] += elementsOther[j];
                  i += _columnStride;
                  j += columnStrideOther;
                }
                idx += _rowStride;
                idxOther += rowStrideOther;
              }
            } else if (multiplicator == -1) {
              // x[i] = x[i] - y[i]
              idx = zero + firstRow * _rowStride;
              idxOther = zeroOther + firstRow * rowStrideOther;
              for (int r = firstRow; r < lastRow; r++) {
                for (int i = idx,
                    j = idxOther,
                    c = 0; c < _columns; c++) {
                  _elements[i] -= elementsOther[j];
                  i += _columnStride;
                  j += columnStrideOther;
                }
                idx += _rowStride;
                idxOther += rowStrideOther;
              }
            } else { // the general case
              // x[i] = x[i] + mult*y[i]
              idx = zero + firstRow * _rowStride;
              idxOther = zeroOther + firstRow * rowStrideOther;
              for (int r = firstRow; r < lastRow; r++) {
                for (int i = idx,
                    j = idxOther,
                    c = 0; c < _columns; c++) {
                  _elements[i] += multiplicator * elementsOther[j];
                  i += _columnStride;
                  j += columnStrideOther;
                }
                idx += _rowStride;
                idxOther += rowStrideOther;
              }
            }
          } else { // the general case x[i] = f(x[i],y[i])
            idx = zero + firstRow * _rowStride;
            idxOther = zeroOther + firstRow * rowStrideOther;
            for (int r = firstRow; r < lastRow; r++) {
              for (int i = idx,
                  j = idxOther,
                  c = 0; c < _columns; c++) {
                _elements[i] = function(_elements[i], elementsOther[j]);
                i += _columnStride;
                j += columnStrideOther;
              }
              idx += _rowStride;
              idxOther += rowStrideOther;
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idx;
    int idxOther;
    // specialized for speed
    if (function == func.mult) {
      // x[i] = x[i] * y[i]
      idx = zero;
      idxOther = zeroOther;
      for (int r = 0; r < _rows; r++) {
        for (int i = idx,
            j = idxOther,
            c = 0; c < _columns; c++) {
          _elements[i] *= elementsOther[j];
          i += _columnStride;
          j += columnStrideOther;
        }
        idx += _rowStride;
        idxOther += rowStrideOther;
      }
    } else if (function == func.div) {
      // x[i] = x[i] / y[i]
      idx = zero;
      idxOther = zeroOther;
      for (int r = 0; r < _rows; r++) {
        for (int i = idx,
            j = idxOther,
            c = 0; c < _columns; c++) {
          _elements[i] /= elementsOther[j];
          i += _columnStride;
          j += columnStrideOther;
        }
        idx += _rowStride;
        idxOther += rowStrideOther;
      }
    } else if (function is DoublePlusMultSecond) {
      double multiplicator = (function as DoublePlusMultSecond).multiplicator;
      if (multiplicator == 0) { // x[i] = x[i] + 0*y[i]
        return this;
      } else if (multiplicator == 1) { // x[i] = x[i] + y[i]
        idx = zero;
        idxOther = zeroOther;
        for (int r = 0; r < _rows; r++) {
          for (int i = idx,
              j = idxOther,
              c = 0; c < _columns; c++) {
            _elements[i] += elementsOther[j];
            i += _columnStride;
            j += columnStrideOther;
          }
          idx += _rowStride;
          idxOther += rowStrideOther;
        }

      } else if (multiplicator == -1) { // x[i] = x[i] - y[i]
        idx = zero;
        idxOther = zeroOther;
        for (int r = 0; r < _rows; r++) {
          for (int i = idx,
              j = idxOther,
              c = 0; c < _columns; c++) {
            _elements[i] -= elementsOther[j];
            i += _columnStride;
            j += columnStrideOther;
          }
          idx += _rowStride;
          idxOther += rowStrideOther;
        }
      } else { // the general case
        // x[i] = x[i] + mult*y[i]
        idx = zero;
        idxOther = zeroOther;
        for (int r = 0; r < _rows; r++) {
          for (int i = idx,
              j = idxOther,
              c = 0; c < _columns; c++) {
            _elements[i] += multiplicator * elementsOther[j];
            i += _columnStride;
            j += columnStrideOther;
          }
          idx += _rowStride;
          idxOther += rowStrideOther;
        }
      }
    } else { // the general case x[i] = f(x[i],y[i])
      idx = zero;
      idxOther = zeroOther;
      for (int r = 0; r < _rows; r++) {
        for (int i = idx,
            j = idxOther,
            c = 0; c < _columns; c++) {
          _elements[i] = function(_elements[i], elementsOther[j]);
          i += _columnStride;
          j += columnStrideOther;
        }
        idx += _rowStride;
        idxOther += rowStrideOther;
      }
    }
    //}
    return this;
  }

  DoubleMatrix2D forEachMatrixRange(final DoubleMatrix2D y, DoubleDoubleFunction function,  /*IntArrayList*/List<int> rowList,  /*IntArrayList*/List<int> columnList) {
    checkShape(y);
    final int size = rowList.length;
    final List<int> rowElements = rowList;//.elements();
    final List<int> columnElements = columnList;//.elements();
    final Float64List elementsOther = y.elements() as Float64List;
    final int zeroOther = y.index(0, 0);
    final int zero = this.index(0, 0);
    final int columnStrideOther = y.columnStride;
    final int rowStrideOther = y.rowStride;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx;
          int idxOther;
          for (int i = firstIdx; i < lastIdx; i++) {
            idx = zero + rowElements[i] * _rowStride + columnElements[i] * _columnStride;
            idxOther = zeroOther + rowElements[i] * rowStrideOther + columnElements[i] * columnStrideOther;
            _elements[idx] = function(_elements[idx], elementsOther[idxOther]);
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idx;
    int idxOther;
    for (int i = 0; i < size; i++) {
      idx = zero + rowElements[i] * _rowStride + columnElements[i] * _columnStride;
      idxOther = zeroOther + rowElements[i] * rowStrideOther + columnElements[i] * columnStrideOther;
      _elements[idx] = function(_elements[idx], elementsOther[idxOther]);
    }
    //}
    return this;
  }

  int get cardinality {
    int cardinality = 0;
    final int zero = this.index(0, 0);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      List<int> results = new List<int>(nthreads);
      int k = _rows ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int cardinality = 0;
          int idx = zero + firstRow * _rowStride;
          for (int r = firstRow; r < lastRow; r++) {
            for (int i = idx,
                c = 0; c < _columns; c++) {
              if (_elements[i] != 0) cardinality++;
              i += _columnStride;
            }
            idx += _rowStride;
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
    int idx = zero;
    for (int r = 0; r < _rows; r++) {
      for (int i = idx,
          c = 0; c < _columns; c++) {
        if (_elements[i] != 0) cardinality++;
        i += _columnStride;
      }
      idx += _rowStride;
    }
    //}
    return cardinality;
  }

  Float64List elements() {
    return _elements;
  }

  DoubleMatrix2D forEachNonZero(func.IntIntDoubleFunction function) {
    final int zero = this.index(0, 0);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = zero + firstRow * _rowStride;
          for (int r = firstRow; r < lastRow; r++) {
            for (int i = idx,
                c = 0; c < _columns; c++) {
              double value = _elements[i];
              if (value != 0) {
                _elements[i] = function(r, c, value);
              }
              i += _columnStride;
            }
            idx += _rowStride;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idx = zero;
    for (int r = 0; r < _rows; r++) {
      for (int i = idx,
          c = 0; c < _columns; c++) {
        double value = _elements[i];
        if (value != 0) {
          _elements[i] = function(r, c, value);
        }
        i += _columnStride;
      }
      idx += _rowStride;
    }
    //}
    return this;
  }

  /**
   * Returns a new matrix that has the same elements as this matrix, but they
   * are addressed internally in column major. This method creates a new
   * object (not a view), so changes in the returned matrix are NOT reflected
   * in this matrix.
   *
   * @return this matrix with elements addressed internally in column major
   */
  /*DenseColumnDoubleMatrix2D getColumnMajor() {
    DenseColumnDoubleMatrix2D R = new DenseColumnDoubleMatrix2D(_rows, _columns);
    final int zeroR = R.index(0, 0) as int;
    final int rowStrideR = R.rowStride();
    final int columnStrideR = R.columnStride();
    final Float64List elementsR = R.elements();
    final int zero = this.index(0, 0);
    int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == (nthreads - 1)) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = zero + (firstRow - 1) * _rowStride;
          int idxR = zeroR + (firstRow - 1) * rowStrideR;
          for (int r = firstRow; r < lastRow; r++) {
            for (int i = idx,
                j = idxR,
                c = 0; r < _columns; c++) {
              elementsR[j] = _elements[i];
              i += _rowStride;
              j += rowStrideR;
            }
            idx += _columnStride;
            idxR += columnStrideR;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {
      int idx = zero;
      int idxR = zeroR;
      for (int r = 0; r < _rows; r++) {
        for (int i = idx,
            j = idxR,
            c = 0; r < _columns; c++) {
          elementsR[j] = _elements[i];
          i += _rowStride;
          j += rowStrideR;
        }
        idx += _columnStride;
        idxR += columnStrideR;
      }
    }
    return R;
  }*/

  DoubleMatrixLocation max() {
    int rowLocation = 0;
    int columnLocation = 0;
    final int zero = this.index(0, 0);
    double maxValue = 0.0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      List<Float64List> results = new List<Float64List>(nthreads);//[2];
      int k = _rows ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          double maxValue = _elements[zero + firstRow * _rowStride];
          int rowLocation = firstRow;
          int colLocation = 0;
          double elem;
          int d = 1;
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = d; c < _columns; c++) {
              elem = _elements[zero + r * _rowStride + c * _columnStride];
              if (maxValue < elem) {
                maxValue = elem;
                rowLocation = r;
                colLocation = c;
              }
            }
            d = 0;
          }
          return [maxValue, rowLocation, colLocation];
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          results[j] = futures[j].get() as Float64List;
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
    maxValue = _elements[zero];
    int d = 1;
    double elem;
    for (int r = 0; r < _rows; r++) {
      for (int c = d; c < _columns; c++) {
        elem = _elements[zero + r * _rowStride + c * _columnStride];
        if (maxValue < elem) {
          maxValue = elem;
          rowLocation = r;
          columnLocation = c;
        }
      }
      d = 0;
    }
    //}
    return new DoubleMatrixLocation(maxValue, rowLocation, columnLocation);
  }

  DoubleMatrixLocation min() {
    int rowLocation = 0;
    int columnLocation = 0;
    final int zero = this.index(0, 0);
    double minValue = 0.0;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      List<Float64List> results = new List<Float64List>(nthreads);//[2];
      int k = _rows ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int rowLocation = firstRow;
          int columnLocation = 0;
          double minValue = _elements[zero + firstRow * _rowStride];
          double elem;
          int d = 1;
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = d; c < _columns; c++) {
              elem = _elements[zero + r * _rowStride + c * _columnStride];
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
          results[j] = futures[j].get() as Float64List;
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
    minValue = _elements[zero];
    int d = 1;
    double elem;
    for (int r = 0; r < _rows; r++) {
      for (int c = d; c < _columns; c++) {
        elem = _elements[zero + r * _rowStride + c * _columnStride];
        if (minValue > elem) {
          minValue = elem;
          rowLocation = r;
          columnLocation = c;
        }
      }
      d = 0;
    }
    //}
    return new DoubleMatrixLocation(minValue, rowLocation, columnLocation);
  }

  void negativeValues(final /*IntArrayList*/List<int> rowList, final /*IntArrayList*/List<int> columnList, final /*DoubleArrayList*/List<double> valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    int idx = this.index(0, 0);
    for (int r = 0; r < _rows; r++) {
      for (int i = idx,
          c = 0; c < _columns; c++) {
        double value = _elements[i];
        if (value < 0) {
          rowList.add(r);
          columnList.add(c);
          valueList.add(value);
        }
        i += _columnStride;
      }
      idx += _rowStride;
    }
  }

  void nonZeros(final /*IntArrayList*/List<int> rowList, final /*IntArrayList*/List<int> columnList, final /*DoubleArrayList*/List<double> valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    int idx = this.index(0, 0);
    for (int r = 0; r < _rows; r++) {
      for (int i = idx,
          c = 0; c < _columns; c++) {
        double value = _elements[i];
        if (value != 0) {
          rowList.add(r);
          columnList.add(c);
          valueList.add(value);
        }
        i += _columnStride;
      }
      idx += _rowStride;
    }
  }

  void positiveValues(final /*IntArrayList*/List<int> rowList, final /*IntArrayList*/List<int> columnList, final /*DoubleArrayList*/List<double> valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    int idx = this.index(0, 0);
    for (int r = 0; r < _rows; r++) {
      for (int i = idx,
          c = 0; c < _columns; c++) {
        double value = _elements[i];
        if (value > 0) {
          rowList.add(r);
          columnList.add(c);
          valueList.add(value);
        }
        i += _columnStride;
      }
      idx += _rowStride;
    }
  }

  double get(int row, int column) {
    return _elements[_rowZero + row * _rowStride + _columnZero + column * _columnStride];
  }

  int index(int row, int column) {
    return _rowZero + row * _rowStride + _columnZero + column * _columnStride;
  }

  DoubleMatrix2D like2D(int rows, int columns) {
    return new DenseDoubleMatrix2D(rows, columns);
  }

  DoubleMatrix1D like1D(int size) {
    return new DenseDoubleMatrix1D(size);
  }

  void set(int row, int column, double value) {
    _elements[_rowZero + row * _rowStride + _columnZero + column * _columnStride] = value;
  }

  List<Float64List> toList() {
    final List<Float64List> values = new List<Float64List>(_rows);//[_columns];
    final int zero = this.index(0, 0);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = zero + firstRow * _rowStride;
          for (int r = firstRow; r < lastRow; r++) {
            Float64List currentRow = values[r];
            for (int i = idx,
                c = 0; c < _columns; c++) {
              currentRow[c] = _elements[i];
              i += _columnStride;
            }
            idx += _rowStride;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idx = zero;
    for (int r = 0; r < _rows; r++) {
      values[r] = new Float64List(_columns);
      Float64List currentRow = values[r];
      for (int i = idx,
          c = 0; c < _columns; c++) {
        currentRow[c] = _elements[i];
        i += _columnStride;
      }
      idx += _rowStride;
    }
    //}
    return values;
  }

  DoubleMatrix1D vectorize() {
    final DenseDoubleMatrix1D v = new DenseDoubleMatrix1D(length);
    final int zero = index(0, 0);
    final int zeroOther = v.index(0);
    final int strideOther = v.stride();
    final Float64List elementsOther = v.elements();
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _columns);
      List<Future> futures = new List<Future>(nthreads);
      int k = _columns ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstColumn = j * k;
        final int lastColumn = (j == nthreads - 1) ? _columns : firstColumn + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = 0;
          int idxOther = zeroOther + firstColumn * _rows;
          for (int c = firstColumn; c < lastColumn; c++) {
            idx = zero + c * _columnStride;
            for (int r = 0; r < _rows; r++) {
              elementsOther[idxOther] = _elements[idx];
              idx += _rowStride;
              idxOther += strideOther;
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idx = zero;
    int idxOther = zeroOther;
    for (int c = 0; c < _columns; c++) {
      idx = zero + c * _columnStride;
      for (int r = 0; r < _rows; r++) {
        elementsOther[idxOther] = _elements[idx];
        idx += _rowStride;
        idxOther += strideOther;
      }
    }
    //}
    return v;
  }

  /*void zAssign8Neighbors(DoubleMatrix2D B, func.Double9Function function) {
    // 1. using only 4-5 out of the 9 cells in "function" is *not* the
    // limiting factor for performance.

    // 2. if the "function" would be hardwired into the innermost loop, a
    // speedup of 1.5-2.0 would be seen
    // but then the multi-purpose interface is gone...

    if (!(B is DenseDoubleMatrix2D)) {
      super.zAssign8Neighbors(B, function);
      return;
    }
    if (function == null) throw new ArgumentError("function must not be null.");
    checkShape(B);
    int r = _rows - 1;
    int c = _columns - 1;
    if (_rows < 3 || _columns < 3) return; // nothing to do

    DenseDoubleMatrix2D BB = B as DenseDoubleMatrix2D;
    int A_rs = _rowStride;
    int B_rs = BB._rowStride;
    int A_cs = _columnStride;
    int B_cs = BB._columnStride;
    Float64List elems = this._elements;
    Float64List B_elems = BB._elements;
    if (elems == null || B_elems == null) throw new Error();

    int A_index = index(1, 1);
    int B_index = BB.index(1, 1);
    for (int i = 1; i < r; i++) {
      double a00, a01, a02;
      double a10, a11, a12;
      double a20, a21, a22;

      int B11 = B_index;

      int A02 = A_index - A_rs - A_cs;
      int A12 = A02 + A_rs;
      int A22 = A12 + A_rs;

      // in each step six cells can be remembered in registers - they
      // don't need to be reread from slow memory
      a00 = elems[A02];
      A02 += A_cs;
      a01 = elems[A02]; // A02+=A_cs;
      a10 = elems[A12];
      A12 += A_cs;
      a11 = elems[A12]; // A12+=A_cs;
      a20 = elems[A22];
      A22 += A_cs;
      a21 = elems[A22]; // A22+=A_cs;

      for (int j = 1; j < c; j++) {
        // in each step 3 instead of 9 cells need to be read from
        // memory.
        a02 = elems[A02 += A_cs];
        a12 = elems[A12 += A_cs];
        a22 = elems[A22 += A_cs];

        B_elems[B11] = function(a00, a01, a02, a10, a11, a12, a20, a21, a22);
        B11 += B_cs;

        // move remembered cells
        a00 = a01;
        a01 = a02;
        a10 = a11;
        a11 = a12;
        a20 = a21;
        a21 = a22;
      }
      A_index += A_rs;
      B_index += B_rs;
    }
  }*/

  DoubleMatrix1D mult(final DoubleMatrix1D y, DoubleMatrix1D z, [final double alpha = 1.0, final double beta = 0.0, final bool transposeA = false]) {
    if (transposeA) {
      return dice().mult(y, z, alpha, beta, false);
    }
    if (z == null) {
      z = new DenseDoubleMatrix1D(_rows);
    }
    if (!(y is DenseDoubleMatrix1D && z is DenseDoubleMatrix1D)) return super.mult(y, z, alpha, beta, transposeA);

    if (_columns != y.length || _rows > z.length) throw new ArgumentError("Incompatible args: " + toStringShort() + ", " + y.toStringShort() + ", " + z.toStringShort());

    final Float64List elemsY = y.elements() as Float64List;
    final Float64List elemsZ = z.elements() as Float64List;
    if (_elements == null || elemsY == null || elemsZ == null) {
      throw new Error();
    }
    final int strideY = y.stride();
    final int strideZ = z.stride();
    final int zero = index(0, 0);
    final int zeroY = y.index(0);
    final int zeroZ = z.index(0);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idxZero = zero + firstRow * _rowStride;
          int idxZeroZ = zeroZ + firstRow * strideZ;
          for (int r = firstRow; r < lastRow; r++) {
            double sum = 0;
            int idx = idxZero;
            int idxY = zeroY;
            for (int c = 0; c < _columns; c++) {
              sum += _elements[idx] * elemsY[idxY];
              idx += _columnStride;
              idxY += strideY;
            }
            elemsZ[idxZeroZ] = alpha * sum + beta * elemsZ[idxZeroZ];
            idxZero += _rowStride;
            idxZeroZ += strideZ;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idxZero = zero;
    int idxZeroZ = zeroZ;
    for (int r = 0; r < _rows; r++) {
      double sum = 0.0;
      int idx = idxZero;
      int idxY = zeroY;
      for (int c = 0; c < _columns; c++) {
        sum += _elements[idx] * elemsY[idxY];
        idx += _columnStride;
        idxY += strideY;
      }
      elemsZ[idxZeroZ] = alpha * sum + beta * elemsZ[idxZeroZ];
      idxZero += _rowStride;
      idxZeroZ += strideZ;
    }
    //}
    return z;
  }

  DoubleMatrix2D multiply(final DoubleMatrix2D B, DoubleMatrix2D C, [final double alpha = 1.0, final double beta = 0.0, final bool transposeA = false, bool transposeB = false]) {
    final int rowsA = _rows;
    final int columnsA = _columns;
    final int rowsB = B.rows;
    final int columnsB = B.columns;
    final int rowsC = transposeA ? columnsA : rowsA;
    final int columnsC = transposeB ? rowsB : columnsB;

    if (C == null) {
      C = new DenseDoubleMatrix2D(rowsC, columnsC);
    }

    /*
    * determine how to split and parallelize best into blocks if more
    * B.columns than tasks --> split B.columns, as follows:
    *
    * xx|xx|xxx B xx|xx|xxx xx|xx|xxx A xxx xx|xx|xxx C xxx xx|xx|xxx xxx
    * xx|xx|xxx xxx xx|xx|xxx xxx xx|xx|xxx
    *
    * if less B.columns than tasks --> split A.rows, as follows:
    *
    * xxxxxxx B xxxxxxx xxxxxxx A xxx xxxxxxx C xxx xxxxxxx --- ------- xxx
    * xxxxxxx xxx xxxxxxx --- ------- xxx xxxxxxx
    */
    if (transposeA) {
      return dice().multiply(B, C, alpha, beta, false, transposeB);
    }
    if (B is SparseDoubleMatrix2D || B is SparseRCDoubleMatrix2D) {
      // exploit quick sparse mult
      // A*B = (B' * A')'
      if (C == null) {
        return B.multiply(this, null, alpha, beta, !transposeB, true).dice();
      } else {
        B.multiply(this, C.dice(), alpha, beta, !transposeB, true);
        return C;
      }
    }
    if (transposeB) {
      return this.multiply(B.dice(), C, alpha, beta, transposeA, false);
    }

    if (!(C is DenseDoubleMatrix2D)) {
      return super.multiply(B, C, alpha, beta, transposeA, transposeB);
    }

    if (B.rows != columnsA) {
      throw new ArgumentError("Matrix2D inner dimensions must agree:" + this.toStringShort() + ", " + B.toStringShort());
    }
    if (C.rows != rowsA || C.columns != columnsB) {
      throw new ArgumentError("Incompatibe result matrix: " + this.toStringShort() + ", " + B.toStringShort() + ", " + C.toStringShort());
    }
    if (this == C || B == C) {
      throw new ArgumentError("Matrices must not be identical");
    }

    int flops = 2 * rowsA * columnsA * columnsB;
    int noOfTasks = 1;//Math.min(flops / 30000, ConcurrencyUtils.getNumberOfThreads()) as int; // each
    /* thread should process at least 30000 flops */
    bool splitB = (columnsB >= noOfTasks);
    int width = splitB ? columnsB : rowsA;
    noOfTasks = Math.min(width, noOfTasks);

    if (noOfTasks < 2) { //parallelization doesn't pay off (too much start up overhead)
      return this.multiplySequential(B, C, alpha, beta, transposeA, transposeB);
    }

    // set up concurrent tasks
    int span = width ~/ noOfTasks;
    final List<Future> subTasks = new List<Future>(noOfTasks);
    for (int i = 0; i < noOfTasks; i++) {
      final int offset = i * span;
      if (i == noOfTasks - 1) span = width - span * i; // last span may be a bit larger

      DoubleMatrix2D AA, BB, CC;
      if (splitB) {
        // split B along columns into blocks
        AA = this;
        BB = B.part(0, offset, columnsA, span);
        CC = C.part(0, offset, rowsA, span);
      } else {
        // split A along rows into blocks
        AA = this.part(offset, 0, span, columnsA);
        BB = B;
        CC = C.part(offset, 0, span, columnsB);
      }

      /*subTasks[i] = ConcurrencyUtils.submit(() {
        (AA as DenseDoubleMatrix2D).zMultSequential(BB, CC, alpha, beta, transposeA, transposeB);
      });*/
      (AA as DenseDoubleMatrix2D).multiplySequential(BB, CC, alpha, beta, transposeA, transposeB);
    }

    //ConcurrencyUtils.waitForCompletion(subTasks);
    return C;
  }

  double sum() {
    double sum = 0.0;
    if (_elements == null) {
      throw new Error();
    }
    final int zero = index(0, 0);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          double sum = 0;
          int idx = zero + firstRow * _rowStride;
          for (int r = firstRow; r < lastRow; r++) {
            for (int i = idx,
                c = 0; c < _columns; c++) {
              sum += _elements[i];
              i += _columnStride;
            }
            idx += _rowStride;
          }
          return sum;
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          sum += futures[j].get() as double;
        }
      } on ExecutionException catch (ex) {
        ex.printStackTrace();
      } on InterruptedException catch (e) {
        e.printStackTrace();
      }
    } else {*/
    int idx = zero;
    for (int r = 0; r < _rows; r++) {
      for (int i = idx,
          c = 0; c < _columns; c++) {
        sum += _elements[i];
        i += _columnStride;
      }
      idx += _rowStride;
    }
    //}
    return sum;
  }

  DoubleMatrix2D multiplySequential(DoubleMatrix2D B, DoubleMatrix2D C, double alpha, double beta, bool transposeA, bool transposeB) {
    if (transposeA) {
      return dice().multiply(B, C, alpha, beta, false, transposeB);
    }
    if (B is SparseDoubleMatrix2D || B is SparseRCDoubleMatrix2D || B is SparseCCDoubleMatrix2D) {
      // exploit quick sparse mult
      // A*B = (B' * A')'
      if (C == null) {
        return B.multiply(this, null, alpha, beta, !transposeB, true).dice();
      } else {
        B.multiply(this, C.dice(), alpha, beta, !transposeB, true);
        return C;
      }
    }
    if (transposeB) {
      return this.multiply(B.dice(), C, alpha, beta, transposeA, false);
    }

    int rowsA = _rows;
    int columnsA = _columns;
    int p = B.columns;
    if (C == null) {
      C = new DenseDoubleMatrix2D(rowsA, p);
    }
    if (!(B is DenseDoubleMatrix2D) || !(C is DenseDoubleMatrix2D)) {
      return super.multiply(B, C, alpha, beta, transposeA, transposeB);
    }
    if (B.rows != columnsA) {
      throw new ArgumentError("Matrix2D inner dimensions must agree:" + toStringShort() + ", " + B.toStringShort());
    }
    if (C.rows != rowsA || C.columns != p) {
      throw new ArgumentError("Incompatibel result matrix: " + toStringShort() + ", " + B.toStringShort() + ", " + C.toStringShort());
    }
    if (this == C || B == C) {
      throw new ArgumentError("Matrices must not be identical");
    }

    DenseDoubleMatrix2D BB = B as DenseDoubleMatrix2D;
    DenseDoubleMatrix2D CC = C as DenseDoubleMatrix2D;
    final Float64List AElems = this._elements;
    final Float64List BElems = BB._elements;
    final Float64List CElems = CC._elements;
    if (AElems == null || BElems == null || CElems == null) {
      throw new Error();
    }

    int cA = this._columnStride;
    int cB = BB._columnStride;
    int cC = CC._columnStride;

    int rA = this._rowStride;
    int rB = BB._rowStride;
    int rC = CC._rowStride;

    /*
   * A is blocked to hide memory latency xxxxxxx B xxxxxxx xxxxxxx A xxx
   * xxxxxxx C xxx xxxxxxx --- ------- xxx xxxxxxx xxx xxxxxxx --- -------
   * xxx xxxxxxx
   */
    final int BLOCK_SIZE = 30000; // * 8 == Level 2 cache in bytes
    int m_optimal = (BLOCK_SIZE - columnsA) ~/ (columnsA + 1);
    if (m_optimal <= 0) m_optimal = 1;
    int blocks = rowsA ~/ m_optimal;
    int rr = 0;
    if (rowsA % m_optimal != 0) blocks++;
    for ( ; --blocks >= 0; ) {
      int jB = BB.index(0, 0);
      int indexA = index(rr, 0);
      int jC = CC.index(rr, 0);
      rr += m_optimal;
      if (blocks == 0) m_optimal += rowsA - rr;

      for (int j = p; --j >= 0; ) {
        int iA = indexA;
        int iC = jC;
        for (int i = m_optimal; --i >= 0; ) {
          int kA = iA;
          int kB = jB;
          double s = 0.0;

          // loop unrolled
          kA -= cA;
          kB -= rB;

          for (int k = columnsA % 4; --k >= 0; ) {
            s += AElems[kA += cA] * BElems[kB += rB];
          }
          for (int k = columnsA ~/ 4; --k >= 0; ) {
            s += AElems[kA += cA] * BElems[kB += rB] + AElems[kA += cA] * BElems[kB += rB] + AElems[kA += cA] * BElems[kB += rB] + AElems[kA += cA] * BElems[kB += rB];
          }

          CElems[iC] = alpha * s + beta * CElems[iC];
          iA += rA;
          iC += rC;
        }
        jB += cB;
        jC += cC;
      }
    }
    return C;
  }

  bool _haveSharedCellsRaw(DoubleMatrix2D other) {
    if (other is SelectedDenseDoubleMatrix2D) {
      return this._elements == other._elements;
    } else if (other is DenseDoubleMatrix2D) {
      return this._elements == other._elements;
    }
    return false;
  }

  DoubleMatrix1D _like1D(int size, int zero, int stride) {
    return new DenseDoubleMatrix1D(size, this._elements, zero, stride, true);
  }

  DoubleMatrix2D _viewSelectionLike(Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedDenseDoubleMatrix2D.offset(this._elements, rowOffsets, columnOffsets, 0);
  }

  Object clone() {
    return new DenseDoubleMatrix2D(_rows, _columns, _elements, _rowZero, _columnZero, _rowStride, _columnStride, !_isNoView);
  }
}

/**
 * Selection view on dense 2-d matrices holding <tt>double</tt> elements. First
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
 * to users. By the way, the same principle applies to concrete DenseXXX and
 * SparseXXX classes: they presents the same functionality with the same
 * signatures and semantics as abstract superclass(es) while introducing no
 * additional functionality. Thus, they need not be visible to users, either.
 * Factory methods could hide all these concrete types.
 * <p>
 * This class uses no delegation. Its instances point directly to the data. Cell
 * addressing overhead is 1 additional int addition and 2 additional array index
 * accesses per get/set.
 * <p>
 * Note that this implementation is not synchronized.
 * <p>
 * <b>Memory requirements:</b>
 * <p>
 * <tt>memory [bytes] = 4*(rowIndexes.length+columnIndexes.length)</tt>. Thus,
 * an index view with 1000 x 1000 indexes additionally uses 8 KB.
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
class SelectedDenseDoubleMatrix2D extends DoubleMatrix2D {

  /**
   * The elements of this matrix.
   */
  Float64List _elements;

  /**
   * The offsets of the visible cells of this matrix.
   */
  Int32List _rowOffsets;

  Int32List _columnOffsets;

  /**
   * The offset.
   */
  int _offset;

  /**
   * Constructs a matrix view with the given parameters.
   *
   * @param elements
   *            the cells.
   * @param rowOffsets
   *            The row offsets of the cells that shall be visible.
   * @param columnOffsets
   *            The column offsets of the cells that shall be visible.
   * @param offset
   */
  factory SelectedDenseDoubleMatrix2D.offset(Float64List elements, Int32List rowOffsets, Int32List columnOffsets, int offset) {
    return new SelectedDenseDoubleMatrix2D(rowOffsets.length, columnOffsets.length, elements, 0, 0, 1, 1, rowOffsets, columnOffsets, offset, true);
  }

  /**
   * Constructs a matrix view with the given parameters.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @param elements
   *            the cells.
   * @param rowZero
   *            the position of the first element.
   * @param columnZero
   *            the position of the first element.
   * @param rowStride
   *            the number of elements between two rows, i.e.
   *            <tt>index(i+1,j)-index(i,j)</tt>.
   * @param columnStride
   *            the number of elements between two columns, i.e.
   *            <tt>index(i,j+1)-index(i,j)</tt>.
   * @param rowOffsets
   *            The row offsets of the cells that shall be visible.
   * @param columnOffsets
   *            The column offsets of the cells that shall be visible.
   * @param offset
   */
  SelectedDenseDoubleMatrix2D(int rows, int columns, Float64List elements, int rowZero, int columnZero, int rowStride, int columnStride, Int32List rowOffsets, Int32List columnOffsets, int offset, bool isView) {
    // be sure parameters are valid, we do not check...
    _setUp(rows, columns, rowZero, columnZero, rowStride, columnStride);

    this._elements = elements;
    this._rowOffsets = rowOffsets;
    this._columnOffsets = columnOffsets;
    this._offset = offset;

    this._isNoView = !isView;
  }

  Float64List elements() {
    return _elements;
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
  double get(int row, int column) {
    // if (debug) if (column<0 || column>=columns || row<0 || row>=rows)
    // throw new IndexOutOfBoundsException("row:"+row+", column:"+column);
    // return elements[index(row,column)];
    // manually inlined:
    return _elements[_offset + _rowOffsets[_rowZero + row * _rowStride] + _columnOffsets[_columnZero + column * _columnStride]];
  }

  /**
   * Returns the position of the given coordinate within the (virtual or
   * non-virtual) internal 1-dimensional array.
   *
   * @param row
   *            the index of the row-coordinate.
   * @param column
   *            the index of the column-coordinate.
   */
  int index(int row, int column) {
    // return this.offset + super.index(row,column);
    // manually inlined:
    return this._offset + _rowOffsets[_rowZero + row * _rowStride] + _columnOffsets[_columnZero + column * _columnStride];
  }

  /**
   * Construct and returns a new empty matrix <i>of the same dynamic type</i>
   * as the receiver, having the specified number of rows and columns. For
   * example, if the receiver is an instance of type
   * <tt>DenseDoubleMatrix2D</tt> the new matrix must also be of type
   * <tt>DenseDoubleMatrix2D</tt>, if the receiver is an instance of type
   * <tt>SparseDoubleMatrix2D</tt> the new matrix must also be of type
   * <tt>SparseDoubleMatrix2D</tt>, etc. In general, the new matrix should
   * have internal parametrization as similar as possible.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @return a new empty matrix of the same dynamic type.
   */
  DoubleMatrix2D like2D(int rows, int columns) {
    return new DenseDoubleMatrix2D(rows, columns);
  }

  /**
   * Construct and returns a new 1-d matrix <i>of the corresponding dynamic
   * type</i>, entirelly independent of the receiver. For example, if the
   * receiver is an instance of type <tt>DenseDoubleMatrix2D</tt> the new
   * matrix must be of type <tt>DenseDoubleMatrix1D</tt>, if the receiver is
   * an instance of type <tt>SparseDoubleMatrix2D</tt> the new matrix must be
   * of type <tt>SparseDoubleMatrix1D</tt>, etc.
   *
   * @param size
   *            the number of cells the matrix shall have.
   * @return a new matrix of the corresponding dynamic type.
   */
  DoubleMatrix1D like1D(int size) {
    return new DenseDoubleMatrix1D(size);
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
  void set(int row, int column, double value) {
    // if (debug) if (column<0 || column>=columns || row<0 || row>=rows)
    // throw new IndexOutOfBoundsException("row:"+row+", column:"+column);
    // elements[index(row,column)] = value;
    // manually inlined:
    _elements[_offset + _rowOffsets[_rowZero + row * _rowStride] + _columnOffsets[_columnZero + column * _columnStride]] = value;
  }

  /**
   * Returns a vector obtained by stacking the columns of the matrix on top of
   * one another.
   *
   * @return
   */
  DoubleMatrix1D vectorize() {
    DenseDoubleMatrix1D v = new DenseDoubleMatrix1D(length);
    int idx = 0;
    for (int c = 0; c < _columns; c++) {
      for (int r = 0; r < _rows; r++) {
        v.set(idx++, get(c, r));
      }
    }
    return v;
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
   * <td valign="top">Matrix1D of size 2:<br>
   * 1, 4</td>
   * </tr>
   * </table>
   *
   * @param the
   *            column to fix.
   * @return a new slice view.
   * @throws IllegalArgumentException
   *             if <tt>column < 0 || column >= columns()</tt>.
   * @see #viewRow(int)
   */
  DoubleMatrix1D column(int column) {
    _checkColumn(column);
    int viewSize = this._rows;
    int viewZero = this._rowZero;
    int viewStride = this._rowStride;
    Int32List viewOffsets = this._rowOffsets;
    int viewOffset = this._offset + _columnOffset(_columnRank(column));
    return new SelectedDenseDoubleMatrix1D(viewSize, this._elements, viewZero, viewStride, viewOffsets, viewOffset);
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
   * <td valign="top">Matrix1D of size 3:<br>
   * 1, 2, 3</td>
   * </tr>
   * </table>
   *
   * @param the
   *            row to fix.
   * @return a new slice view.
   * @throws IndexOutOfBoundsException
   *             if <tt>row < 0 || row >= rows()</tt>.
   * @see #viewColumn(int)
   */
  DoubleMatrix1D row(int row) {
    _checkRow(row);
    int viewSize = this._columns;
    int viewZero = _columnZero;
    int viewStride = this._columnStride;
    Int32List viewOffsets = this._columnOffsets;
    int viewOffset = this._offset + _rowOffset(_rowRank(row));
    return new SelectedDenseDoubleMatrix1D(viewSize, this._elements, viewZero, viewStride, viewOffsets, viewOffset);
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
  int _columnOffset(int absRank) {
    return _columnOffsets[absRank];
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
  int _rowOffset(int absRank) {
    return _rowOffsets[absRank];
  }

  /**
   * Returns <tt>true</tt> if both matrices share common cells. More formally,
   * returns <tt>true</tt> if <tt>other != null</tt> and at least one of the
   * following conditions is met
   * <ul>
   * <li>the receiver is a view of the other matrix
   * <li>the other matrix is a view of the receiver
   * <li><tt>this == other</tt>
   * </ul>
   */
  bool _haveSharedCellsRaw(DoubleMatrix2D other) {
    if (other is SelectedDenseDoubleMatrix2D) {
      return this._elements == other._elements;
    } else if (other is DenseDoubleMatrix2D) {
      return this._elements == other._elements;
    }
    return false;
  }

  /**
   * Construct and returns a new 1-d matrix <i>of the corresponding dynamic
   * type</i>, sharing the same cells. For example, if the receiver is an
   * instance of type <tt>DenseDoubleMatrix2D</tt> the new matrix must be of
   * type <tt>DenseDoubleMatrix1D</tt>, if the receiver is an instance of type
   * <tt>SparseDoubleMatrix2D</tt> the new matrix must be of type
   * <tt>SparseDoubleMatrix1D</tt>, etc.
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
  DoubleMatrix1D _like1D(int size, int zero, int stride) {
    throw new Error(); // this method is never called since
    // viewRow() and viewColumn are overridden
    // properly.
  }

  /**
   * Sets up a matrix with a given number of rows and columns.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @throws IllegalArgumentException
   *             if <tt>(double)columns*rows > Integer.MAX_VALUE</tt>.
   */
  void _setUp2D(int rows, int columns) {
    super._setUp(rows, columns);
    this._rowStride = 1;
    this._columnStride = 1;
    this._offset = 0;
  }

  /**
   * Self modifying version of viewDice().
   */
  AbstractMatrix2D _vDice() {
    super._vDice();
    // swap
    Int32List tmp = _rowOffsets;
    _rowOffsets = _columnOffsets;
    _columnOffsets = tmp;

    // flips stay unaffected

    this._isNoView = false;
    return this;
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
  DoubleMatrix2D _viewSelectionLike(Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedDenseDoubleMatrix2D.offset(this._elements, rowOffsets, columnOffsets, this._offset);
  }

  Object clone() {
    return new SelectedDenseDoubleMatrix2D(_rows, _columns, _elements, _rowZero, _columnZero, _rowStride, _columnStride, _rowOffsets, _columnOffsets, _offset, !_isNoView);
  }
}
