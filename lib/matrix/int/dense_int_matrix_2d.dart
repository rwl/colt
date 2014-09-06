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
 * Dense 2-d matrix holding <tt>int</tt> elements. First see the <a
 * href="package-summary.html">package summary</a> and javadoc <a
 * href="package-tree.html">tree view</a> to get the broad picture.
 * <p>
 * <b>Implementation:</b>
 * <p>
 * Internally holds one single contigous one-dimensional array, addressed in row
 * major. Note that this implementation is not synchronized.
 * <p>
 * <b>Memory requirements:</b>
 * <p>
 * <tt>memory [bytes] = 8*rows()*columns()</tt>. Thus, a 1000*1000 matrix uses 8
 * MB.
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
class DenseIntMatrix2D extends IntMatrix2D {

    /**
     * The elements of this matrix. elements are stored in row major, i.e.
     * index==row*columns + column columnOf(index)==index%columns
     * rowOf(index)==index/columns i.e. {row0 column0..m}, {row1 column0..m},
     * ..., {rown column0..m}
     */
    Int32List _elements;

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
    factory DenseIntMatrix2D.fromList(List<Int32List> values) {
        return new DenseIntMatrix2D(values.length, values.length == 0 ? 0 : values[0].length)
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
     *             <tt>rows<0 || columns<0 || (int)columns*rows > Int.MAX_VALUE</tt>
     *             .
     */
    /*DenseIntMatrix2D(int rows, int columns) {
        _setUp(rows, columns);
        this._elements = new Int32List(rows * columns);
    }*/

    /**
     * Constructs a view with the given parameters.
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
     *             <tt>rows<0 || columns<0 || (int)columns*rows > Int.MAX_VALUE</tt>
     *             or flip's are illegal.
     */
    DenseIntMatrix2D(int rows, int columns, [Int32List elements=null, int rowZero=0, int columnZero=0, int rowStride=null,
            int columnStride=1, bool isView=false]) {
        _setUp(rows, columns, rowZero, columnZero, rowStride, columnStride);
        this._elements = elements;
        this._isNoView = !isView;
    }

    int reduce(final ifunc.IntIntFunction aggr, final ifunc.IntFunction f) {
        if (length == 0) {
            throw new ArgumentError("size == 0");
        }
        final int zero = index(0, 0);
        int a = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int a = f(_elements[zero + firstRow * _rowStride]);
                        int d = 1;
                        for (int r = firstRow; r < lastRow; r++) {
                            for (int c = d; c < _columns; c++) {
                                a = aggr(a, f(_elements[zero + r * _rowStride + c * _columnStride]));
                            }
                            d = 0;
                        }
                        return a;
                });
            }
            a = ConcurrencyUtils.waitForCompletion(futures, aggr);
        } else {
            a = f(_elements[zero]);
            int d = 1; // first cell already done
            for (int r = 0; r < _rows; r++) {
                for (int c = d; c < _columns; c++) {
                    a = aggr(a, f(_elements[zero + r * _rowStride + c * _columnStride]));
                }
                d = 0;
            }
        }
        return a;
    }

    int reduceWhere(final ifunc.IntIntFunction aggr,
            final ifunc.IntFunction f, final ifunc.IntProcedure cond) {
        if (size() == 0)
            throw new ArgumentError("size == 0");
        final int zero = index(0, 0);
        int a = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int elem = _elements[zero + firstRow * _rowStride];
                        int a = 0;
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
        } else {
            int elem = _elements[zero];
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
        }
        return a;
    }

    int reduceRange(final ifunc.IntIntFunction aggr,
            final ifunc.IntFunction f, final IntArrayList rowList, final IntArrayList columnList) {
        if (length == 0) {
            throw new ArgumentError("length == 0");
        }
        final int zero = index(0, 0);
        final int size = rowList.size();
        final Int32List rowElements = rowList.elements();
        final Int32List columnElements = columnList.elements();
        int a = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, size);
            List<Future> futures = new List<Future>(nthreads);
            int k = size / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstIdx = j * k;
                final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int a = f(_elements[zero + rowElements[firstIdx] * _rowStride + columnElements[firstIdx]
                                * _columnStride]);
                        int elem;
                        for (int i = firstIdx + 1; i < lastIdx; i++) {
                            elem = _elements[zero + rowElements[i] * _rowStride + columnElements[i] * _columnStride];
                            a = aggr(a, f(elem));
                        }
                        return a;
                });
            }
            a = ConcurrencyUtils.waitForCompletion(futures, aggr);
        } else {
            int elem;
            a = f(_elements[zero + rowElements[0] * _rowStride + columnElements[0] * _columnStride]);
            for (int i = 1; i < size; i++) {
                elem = _elements[zero + rowElements[i] * _rowStride + columnElements[i] * _columnStride];
                a = aggr(a, f(elem));
            }
        }
        return a;
    }

    int reduceWith(final IntMatrix2D other, final ifunc.IntIntFunction aggr,
            final ifunc.IntIntFunction f) {
        if (!(other is DenseIntMatrix2D)) {
            return super.reduceWith(other, aggr, f);
        }
        checkShape(other);
        if (size() == 0)
            throw new ArgumentError("size == 0");
        final int zero = index(0, 0);
        final int zeroOther = other.index(0, 0);
        final int rowStrideOther = other.rowStride();
        final int colStrideOther = other.columnStride();
        final Int32List elemsOther = other.elements() as Int32List;
        int a = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int a = f(_elements[zero + firstRow * _rowStride], elemsOther[zeroOther + firstRow
                                * rowStrideOther]);
                        int d = 1;
                        for (int r = firstRow; r < lastRow; r++) {
                            for (int c = d; c < _columns; c++) {
                                a = aggr(a, f(_elements[zero + r * _rowStride + c * _columnStride],
                                        elemsOther[zeroOther + r * rowStrideOther + c * colStrideOther]));
                            }
                            d = 0;
                        }
                        return a;
                });
            }
            a = ConcurrencyUtils.waitForCompletion(futures, aggr);
        } else {
            int d = 1; // first cell already done
            a = f(_elements[zero], elemsOther[zeroOther]);
            for (int r = 0; r < _rows; r++) {
                for (int c = d; c < _columns; c++) {
                    a = aggr(a, f(_elements[zero + r * _rowStride + c * _columnStride], elemsOther[zeroOther
                            + r * rowStrideOther + c * colStrideOther]));
                }
                d = 0;
            }
        }
        return a;
    }

    IntMatrix2D forEach(final ifunc.IntFunction function) {
        final Int32List elems = this._elements;
        if (elems == null) {
            throw new InternalError();
        }
        final int zero = index(0, 0);
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            if (function is ifunc.IntMult) { // x[i] =
                // mult*x[i]
                int multiplicator = (function as ifunc.IntMult).multiplicator;
                if (multiplicator == 1) {
                    return this;
                }
                if (multiplicator == 0) {
                    return fill(0);
                }
            }
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int idx = zero + firstRow * _rowStride;
                        // specialization for speed
                        if (function is ifunc.IntMult) {
                            // x[i] = mult*x[i]
                            int multiplicator = (function as ifunc.IntMult).multiplicator;
                            if (multiplicator == 1)
                                return;
                            for (int r = firstRow; r < lastRow; r++) {
                                for (int i = idx, c = 0; c < _columns; c++) {
                                    elems[i] *= multiplicator;
                                    i += _columnStride;
                                }
                                idx += _rowStride;
                            }
                        } else {
                            // the general case x[i] = f(x[i])
                            for (int r = firstRow; r < lastRow; r++) {
                                for (int i = idx, c = 0; c < _columns; c++) {
                                    elems[i] = function(elems[i]);
                                    i += _columnStride;
                                }
                                idx += _rowStride;
                            }
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            int idx = zero;
            // specialization for speed
            if (function is ifunc.IntMult) { // x[i] =
                // mult*x[i]
                int multiplicator = (function as ifunc.IntMult).multiplicator;
                if (multiplicator == 1) {
                    return this;
                }
                if (multiplicator == 0) {
                    return fill(0);
                }
                for (int r = 0; r < _rows; r++) { // the general case
                    for (int i = idx, c = 0; c < _columns; c++) {
                        elems[i] *= multiplicator;
                        i += _columnStride;
                    }
                    idx += _rowStride;
                }
            } else { // the general case x[i] = f(x[i])
                for (int r = 0; r < _rows; r++) {
                    for (int i = idx, c = 0; c < _columns; c++) {
                        elems[i] = function(elems[i]);
                        i += _columnStride;
                    }
                    idx += _rowStride;
                }
            }
        }
        return this;
    }

    IntMatrix2D forEachWhere(final ifunc.IntProcedure cond,
            final ifunc.IntFunction function) {
        final int zero = index(0, 0);
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int elem;
                        int idx = zero + firstRow * _rowStride;
                        for (int r = firstRow; r < lastRow; r++) {
                            for (int i = idx, c = 0; c < _columns; c++) {
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
        } else {
            int elem;
            int idx = zero;
            for (int r = 0; r < _rows; r++) {
                for (int i = idx, c = 0; c < _columns; c++) {
                    elem = _elements[i];
                    if (cond(elem) == true) {
                        _elements[i] = function(elem);
                    }
                    i += _columnStride;
                }
                idx += _rowStride;
            }
        }
        return this;
    }

    IntMatrix2D fillWhere(final ifunc.IntProcedure cond, final int value) {
        final int zero = index(0, 0);
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int elem;
                        int idx = zero + firstRow * _rowStride;
                        for (int r = firstRow; r < lastRow; r++) {
                            for (int i = idx, c = 0; c < _columns; c++) {
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
        } else {
            int elem;
            int idx = zero;
            for (int r = 0; r < _rows; r++) {
                for (int i = idx, c = 0; c < _columns; c++) {
                    elem = _elements[i];
                    if (cond(elem) == true) {
                        _elements[i] = value;
                    }
                    i += _columnStride;
                }
                idx += _rowStride;
            }
        }
        return this;
    }

    IntMatrix2D fill(final int value) {
        final Int32List elems = this._elements;
        final int zero = index(0, 0);
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int idx = zero + firstRow * _rowStride;
                        for (int r = firstRow; r < lastRow; r++) {
                            for (int i = idx, c = 0; c < _columns; c++) {
                                elems[i] = value;
                                i += _columnStride;
                            }
                            idx += _rowStride;
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            int idx = zero;
            for (int r = 0; r < _rows; r++) {
                for (int i = idx, c = 0; c < _columns; c++) {
                    elems[i] = value;
                    i += _columnStride;
                }
                idx += _rowStride;
            }
        }
        return this;
    }

    IntMatrix2D setAll(final Int32List values) {
        if (values.length != size())
            throw new ArgumentError("Must have same length: length=" + values.length + " rows()*columns()="
                    + rows() * columns());
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if (this._isNoView) {
            System.arraycopy(values, 0, this._elements, 0, values.length);
        } else {
            final int zero = index(0, 0);
            if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
                nthreads = Math.min(nthreads, _rows);
                List<Future> futures = new List<Future>(nthreads);
                int k = _rows / nthreads;
                for (int j = 0; j < nthreads; j++) {
                    final int firstRow = j * k;
                    final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

                    futures[j] = ConcurrencyUtils.submit(() {
                            int idxOther = firstRow * _columns;
                            int idx = zero + firstRow * _rowStride;
                            for (int r = firstRow; r < lastRow; r++) {
                                for (int i = idx, c = 0; c < _columns; c++) {
                                    _elements[i] = values[idxOther++];
                                    i += _columnStride;
                                }
                                idx += _rowStride;
                            }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {

                int idxOther = 0;
                int idx = zero;
                for (int r = 0; r < _rows; r++) {
                    for (int i = idx, c = 0; c < _columns; c++) {
                        _elements[i] = values[idxOther++];
                        i += _columnStride;
                    }
                    idx += _rowStride;
                }
            }
        }
        return this;
    }

    IntMatrix2D setAll2D(final List<Int32List> values) {
        if (values.length != _rows) {
            throw new ArgumentError("Must have same number of rows: rows=${values.length} rows()=${rows()}");
        }
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if (this._isNoView) {
            if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
                nthreads = Math.min(nthreads, _rows);
                List<Future> futures = new List<Future>(nthreads);
                int k = _rows / nthreads;
                for (int j = 0; j < nthreads; j++) {
                    final int firstRow = j * k;
                    final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

                    futures[j] = ConcurrencyUtils.submit(() {
                            int i = firstRow * _rowStride;
                            for (int r = firstRow; r < lastRow; r++) {
                                Int32List currentRow = values[r];
                                if (currentRow.length != _columns)
                                    throw new ArgumentError(
                                            "Must have same number of columns in every row: columns="
                                                    + currentRow.length + "columns()=" + columns());
                                System.arraycopy(currentRow, 0, _elements, i, _columns);
                                i += _columns;
                            }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {
                int i = 0;
                for (int r = 0; r < _rows; r++) {
                    Int32List currentRow = values[r];
                    if (currentRow.length != _columns)
                        throw new ArgumentError("Must have same number of columns in every row: columns="
                                + currentRow.length + "columns()=" + columns());
                    System.arraycopy(currentRow, 0, this._elements, i, _columns);
                    i += _columns;
                }
            }
        } else {
            final int zero = index(0, 0);
            if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
                nthreads = Math.min(nthreads, _rows);
                List<Future> futures = new List<Future>(nthreads);
                int k = _rows / nthreads;
                for (int j = 0; j < nthreads; j++) {
                    final int firstRow = j * k;
                    final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

                    futures[j] = ConcurrencyUtils.submit(() {
                            int idx = zero + firstRow * _rowStride;
                            for (int r = firstRow; r < lastRow; r++) {
                                Int32List currentRow = values[r];
                                if (currentRow.length != _columns)
                                    throw new ArgumentError(
                                            "Must have same number of columns in every row: columns="
                                                    + currentRow.length + "columns()=" + columns());
                                for (int i = idx, c = 0; c < _columns; c++) {
                                    _elements[i] = currentRow[c];
                                    i += _columnStride;
                                }
                                idx += _rowStride;
                            }
                    });
                }
                ConcurrencyUtils.waitForCompletion(futures);
            } else {
                int idx = zero;
                for (int r = 0; r < _rows; r++) {
                    Int32List currentRow = values[r];
                    if (currentRow.length != _columns)
                        throw new ArgumentError("Must have same number of columns in every row: columns="
                                + currentRow.length + "columns()=" + columns());
                    for (int i = idx, c = 0; c < _columns; c++) {
                        _elements[i] = currentRow[c];
                        i += _columnStride;
                    }
                    idx += _rowStride;
                }
            }
            return this;
        }
        return this;
    }

    IntMatrix2D copyFrom(final IntMatrix2D source) {
        // overriden for performance only
        if (!(source is DenseIntMatrix2D)) {
            super.copyFrom(source);
            return this;
        }
        final DenseIntMatrix2D other_final = source as DenseIntMatrix2D;
        if (other_final == this)
            return this; // nothing to do
        checkShape(other_final);
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if (this._isNoView && other_final._isNoView) { // quickest
            System.arraycopy(other_final._elements, 0, this._elements, 0, this._elements.length);
            return this;
        }
        DenseIntMatrix2D other = source as DenseIntMatrix2D;
        if (_haveSharedCells(other)) {
            IntMatrix2D c = other.copy();
            if (!(c is DenseIntMatrix2D)) { // should not happen
                super.copyFrom(other);
                return this;
            }
            other = c as DenseIntMatrix2D;
        }

        final Int32List elemsOther = other._elements;
        if (_elements == null || elemsOther == null) {
            throw new Error();
        }
        final int zeroOther = other.index(0, 0);
        final int zero = index(0, 0);
        final int columnStrideOther = other._columnStride;
        final int rowStrideOther = other._rowStride;
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

                futures[j] = ConcurrencyUtils.submit(() {
                        int idx = zero + firstRow * _rowStride;
                        int idxOther = zeroOther + firstRow * rowStrideOther;
                        for (int r = firstRow; r < lastRow; r++) {
                            for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
                                _elements[i] = elemsOther[j];
                                i += _columnStride;
                                j += columnStrideOther;
                            }
                            idx += _rowStride;
                            idxOther += rowStrideOther;
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            int idx = zero;
            int idxOther = zeroOther;
            for (int r = 0; r < _rows; r++) {
                for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
                    _elements[i] = elemsOther[j];
                    i += _columnStride;
                    j += columnStrideOther;
                }
                idx += _rowStride;
                idxOther += rowStrideOther;
            }
        }
        return this;
    }

    IntMatrix2D forEachWith(final IntMatrix2D y, final ifunc.IntIntFunction function) {
        // overriden for performance only
        if (!(y is DenseIntMatrix2D)) {
            super.forEachWith(y, function);
            return this;
        }
        DenseIntMatrix2D other = y as DenseIntMatrix2D;
        checkShape(y);
        final Int32List elemsOther = other._elements;
        if (_elements == null || elemsOther == null) {
            throw new Error();
        }
        final int zeroOther = other.index(0, 0);
        final int zero = index(0, 0);
        final int columnStrideOther = other._columnStride;
        final int rowStrideOther = other._rowStride;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            if (function is ifunc.IntPlusMultSecond) {
                int multiplicator = (function as ifunc.IntPlusMultSecond).multiplicator;
                if (multiplicator == 0) { // x[i] = x[i] + 0*y[i]
                    return this;
                }
            }
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

                futures[j] = ConcurrencyUtils.submit(() {
                        int idx;
                        int idxOther;
                        // specialized for speed
                        if (function == ifunc.IntFunctions.mult) {
                            // x[i] = x[i]*y[i]
                            idx = zero + firstRow * _rowStride;
                            idxOther = zeroOther + firstRow * rowStrideOther;
                            for (int r = firstRow; r < lastRow; r++) {
                                for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
                                    _elements[i] *= elemsOther[j];
                                    i += _columnStride;
                                    j += columnStrideOther;
                                }
                                idx += _rowStride;
                                idxOther += rowStrideOther;
                            }
                        } else if (function == ifunc.IntFunctions.div) {
                            // x[i] = x[i] / y[i]
                            idx = zero + firstRow * _rowStride;
                            idxOther = zeroOther + firstRow * rowStrideOther;
                            for (int r = firstRow; r < lastRow; r++) {
                                for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
                                    _elements[i] /= elemsOther[j];
                                    i += _columnStride;
                                    j += columnStrideOther;
                                }
                                idx += _rowStride;
                                idxOther += rowStrideOther;
                            }
                        } else if (function is ifunc.IntPlusMultSecond) {
                            int multiplicator = (function as ifunc.IntPlusMultSecond).multiplicator;
                            if (multiplicator == 1) {
                                // x[i] = x[i] + y[i]
                                idx = zero + firstRow * _rowStride;
                                idxOther = zeroOther + firstRow * rowStrideOther;
                                for (int r = firstRow; r < lastRow; r++) {
                                    for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
                                        _elements[i] += elemsOther[j];
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
                                    for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
                                        _elements[i] -= elemsOther[j];
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
                                    for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
                                        _elements[i] += multiplicator * elemsOther[j];
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
                                for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
                                    _elements[i] = function(_elements[i], elemsOther[j]);
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
        } else {
            int idx;
            int idxOther;
            // specialized for speed
            if (function == ifunc.IntFunctions.mult) {
                // x[i] = x[i] * y[i]
                idx = zero;
                idxOther = zeroOther;
                for (int r = 0; r < _rows; r++) {
                    for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
                        _elements[i] *= elemsOther[j];
                        i += _columnStride;
                        j += columnStrideOther;
                    }
                    idx += _rowStride;
                    idxOther += rowStrideOther;
                }
            } else if (function == ifunc.IntFunctions.div) {
                // x[i] = x[i] / y[i]
                idx = zero;
                idxOther = zeroOther;
                for (int r = 0; r < _rows; r++) {
                    for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
                        _elements[i] /= elemsOther[j];
                        i += _columnStride;
                        j += columnStrideOther;
                    }
                    idx += _rowStride;
                    idxOther += rowStrideOther;
                }
            } else if (function is ifunc.IntPlusMultSecond) {
                int multiplicator = (function as ifunc.IntPlusMultSecond).multiplicator;
                if (multiplicator == 0) { // x[i] = x[i] + 0*y[i]
                    return this;
                } else if (multiplicator == 1) { // x[i] = x[i] + y[i]
                    idx = zero;
                    idxOther = zeroOther;
                    for (int r = 0; r < _rows; r++) {
                        for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
                            _elements[i] += elemsOther[j];
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
                        for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
                            _elements[i] -= elemsOther[j];
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
                        for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
                            _elements[i] += multiplicator * elemsOther[j];
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
                    for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
                        _elements[i] = function(_elements[i], elemsOther[j]);
                        i += _columnStride;
                        j += columnStrideOther;
                    }
                    idx += _rowStride;
                    idxOther += rowStrideOther;
                }
            }
        }
        return this;
    }

    IntMatrix2D forEachWithRange(final IntMatrix2D y, final ifunc.IntIntFunction function,
            IntArrayList rowList, IntArrayList columnList) {
        checkShape(y);
        final int size = rowList.size();
        final Int32List rowElements = rowList.elements();
        final Int32List columnElements = columnList.elements();
        final Int32List elemsOther = y.elements() as Int32List;
        final int zeroOther = y.index(0, 0);
        final int zero = index(0, 0);
        final int columnStrideOther = y.columnStride();
        final int rowStrideOther = y.rowStride();
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, size);
            List<Future> futures = new List<Future>(nthreads);
            int k = size / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstIdx = j * k;
                final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int idx;
                        int idxOther;
                        for (int i = firstIdx; i < lastIdx; i++) {
                            idx = zero + rowElements[i] * _rowStride + columnElements[i] * _columnStride;
                            idxOther = zeroOther + rowElements[i] * rowStrideOther + columnElements[i]
                                    * columnStrideOther;
                            _elements[idx] = function(_elements[idx], elemsOther[idxOther]);
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            int idx;
            int idxOther;
            for (int i = 0; i < size; i++) {
                idx = zero + rowElements[i] * _rowStride + columnElements[i] * _columnStride;
                idxOther = zeroOther + rowElements[i] * rowStrideOther + columnElements[i] * columnStrideOther;
                _elements[idx] = function(_elements[idx], elemsOther[idxOther]);
            }
        }
        return this;
    }

    int cardinality() {
        int cardinality = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        final int zero = index(0, 0);
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            List<int> results = new List<int>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

                futures[j] = ConcurrencyUtils.submit(() {
                        int cardinality = 0;
                        int idx = zero + firstRow * _rowStride;
                        for (int r = firstRow; r < lastRow; r++) {
                            for (int i = idx, c = 0; c < _columns; c++) {
                                if (_elements[i] != 0)
                                    cardinality++;
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
        } else {
            int idx = zero;
            for (int r = 0; r < _rows; r++) {
                for (int i = idx, c = 0; c < _columns; c++) {
                    if (_elements[i] != 0)
                        cardinality++;
                    i += _columnStride;
                }
                idx += _rowStride;
            }
        }
        return cardinality;
    }

    Int32List elements() {
        return _elements;
    }

    IntMatrix2D forEachNonZero(final ifunc.IntIntIntFunction function) {
        final int zero = index(0, 0);
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

                futures[j] = ConcurrencyUtils.submit(() {
                        int idx = zero + firstRow * _rowStride;
                        for (int r = firstRow; r < lastRow; r++) {
                            for (int i = idx, c = 0; c < _columns; c++) {
                                int value = _elements[i];
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
        } else {
            int idx = zero;
            for (int r = 0; r < _rows; r++) {
                for (int i = idx, c = 0; c < _columns; c++) {
                    int value = _elements[i];
                    if (value != 0) {
                        _elements[i] = function(r, c, value);
                    }
                    i += _columnStride;
                }
                idx += _rowStride;
            }
        }
        return this;
    }

    void negativeValues(final IntArrayList rowList, final IntArrayList columnList,
            final IntArrayList valueList) {
        rowList.clear();
        columnList.clear();
        valueList.clear();
        int idx = index(0, 0);
        for (int r = 0; r < _rows; r++) {
            for (int i = idx, c = 0; c < _columns; c++) {
                int value = _elements[i];
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

    void nonZeros(final IntArrayList rowList, final IntArrayList columnList, final IntArrayList valueList) {
        rowList.clear();
        columnList.clear();
        valueList.clear();
        int idx = index(0, 0);
        for (int r = 0; r < _rows; r++) {
            for (int i = idx, c = 0; c < _columns; c++) {
                int value = _elements[i];
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

    void positiveValues(final IntArrayList rowList, final IntArrayList columnList,
            final IntArrayList valueList) {
        rowList.clear();
        columnList.clear();
        valueList.clear();
        int idx = index(0, 0);
        for (int r = 0; r < _rows; r++) {
            for (int i = idx, c = 0; c < _columns; c++) {
                int value = _elements[i];
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

    int get(int row, int column) {
        return _elements[_rowZero + row * _rowStride + _columnZero + column * _columnStride];
    }

    long index(int row, int column) {
        return _rowZero + row * _rowStride + _columnZero + column * _columnStride;
    }

    IntMatrix2D like(int rows, int columns) {
        return new DenseIntMatrix2D(rows, columns);
    }

    IntMatrix1D like1D(int size) {
        return new DenseIntMatrix1D(size);
    }

    Int32List max() {
        int rowLocation = 0;
        int columnLocation = 0;
        final int zero = index(0, 0);
        int maxValue = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            List<Int32List> results = new List<Int32List>(nthreads);//[2];
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

                futures[j] = ConcurrencyUtils.submit(() {
                        int maxValue = _elements[zero + firstRow * _rowStride];
                        int rowLocation = firstRow;
                        int colLocation = 0;
                        int elem;
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
                        return [ maxValue, rowLocation, colLocation ];
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
        } else {
            maxValue = _elements[zero];
            int d = 1;
            int elem;
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
        }
        return [ maxValue, rowLocation, columnLocation ];
    }

    Int32List min() {
        int rowLocation = 0;
        int columnLocation = 0;
        final int zero = index(0, 0);
        int minValue = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
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
                        int minValue = _elements[zero + firstRow * _rowStride];
                        int elem;
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
                        return [ minValue, rowLocation, columnLocation ];
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
        } else {
            minValue = _elements[zero];
            int d = 1;
            int elem;
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
        }
        return [ minValue, rowLocation, columnLocation ];
    }

    void set(int row, int column, int value) {
        _elements[_rowZero + row * _rowStride + _columnZero + column * _columnStride] = value;
    }

    List<Int32List> toList() {
        final List<Int32List> values = new List<Int32List>.generate(_rows,
            (_) => new Int32List(_columns));
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        final int zero = index(0, 0);
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;

                futures[j] = ConcurrencyUtils.submit(() {
                        int idx = zero + firstRow * _rowStride;
                        for (int r = firstRow; r < lastRow; r++) {
                            Int32List currentRow = values[r];
                            for (int i = idx, c = 0; c < _columns; c++) {
                                currentRow[c] = _elements[i];
                                i += _columnStride;
                            }
                            idx += _rowStride;
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            int idx = zero;
            for (int r = 0; r < _rows; r++) {
                Int32List currentRow = values[r];
                for (int i = idx, c = 0; c < _columns; c++) {
                    currentRow[c] = _elements[i];
                    i += _columnStride;
                }
                idx += _rowStride;
            }
        }
        return values;
    }

    IntMatrix1D vectorize() {
        final DenseIntMatrix1D v = new DenseIntMatrix1D(length);
        final int zero = index(0, 0);
        final int zeroOther = v.index(0);
        final int strideOther = v.stride();
        final Int32List elemsOther = v.elements();
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _columns);
            List<Future> futures = new List<Future>(nthreads);
            int k = _columns / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstColumn = j * k;
                final int lastColumn = (j == nthreads - 1) ? _columns : firstColumn + k;
                final int startidx = j * k * _rows;

                futures[j] = ConcurrencyUtils.submit(() {
                        int idx = 0;
                        int idxOther = zeroOther + startidx * strideOther;
                        for (int c = firstColumn; c < lastColumn; c++) {
                            idx = zero + c * _columnStride;
                            for (int r = 0; r < _rows; r++) {
                                elemsOther[idxOther] = _elements[idx];
                                idx += _rowStride;
                                idxOther += strideOther;
                            }
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            int idx = zero;
            int idxOther = zeroOther;
            for (int c = 0; c < _columns; c++) {
                idx = zero + c * _columnStride;
                for (int r = 0; r < _rows; r++) {
                    elemsOther[idxOther] = _elements[idx];
                    idx += _rowStride;
                    idxOther += strideOther;
                }
            }
        }
        return v;
    }

    IntMatrix1D mult(final IntMatrix1D y, IntMatrix1D z, final int alpha, final int beta,
            final bool transposeA) {
        if (transposeA)
            return dice().mult(y, z, alpha, beta, false);
        if (z == null) {
            z = new DenseIntMatrix1D(_rows);
        }
        if (!(y is DenseIntMatrix1D && z is DenseIntMatrix1D))
            return super.mult(y, z, alpha, beta, transposeA);

        if (_columns != y.size() || _rows > z.size())
            throw new ArgumentError("Incompatible args: " + toStringShort() + ", " + y.toStringShort()
                    + ", " + z.toStringShort());

        final Int32List elemsY = y.elements() as Int32List;
        final Int32List elemsZ = z.elements() as Int32List;
        if (_elements == null || elemsY == null || elemsZ == null) {
            throw new Error();
        }
        final int strideY = y.stride();
        final int strideZ = z.stride();
        final int zero = index(0, 0);
        final int zeroY = y.index(0);
        final int zeroZ = z.index(0);
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int idxZero = zero + firstRow * _rowStride;
                        int idxZeroZ = zeroZ + firstRow * strideZ;
                        for (int r = firstRow; r < lastRow; r++) {
                            int sum = 0;
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
        } else {
            int idxZero = zero;
            int idxZeroZ = zeroZ;
            for (int r = 0; r < _rows; r++) {
                int sum = 0;
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
        }
        return z;
    }

    IntMatrix2D multiply(final IntMatrix2D B, IntMatrix2D C, final int alpha, final int beta,
            final bool transposeA, final bool transposeB) {
        final int rowsA = _rows;
        final int columnsA = _columns;
        final int rowsB = B.rows();
        final int columnsB = B.columns();
        final int rowsC = transposeA ? columnsA : rowsA;
        final int columnsC = transposeB ? rowsB : columnsB;

        if (C == null) {
            C = new DenseIntMatrix2D(rowsC, columnsC);
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
        if (transposeA)
            return dice().multiply(B, C, alpha, beta, false, transposeB);
        if (B is SparseIntMatrix2D || B is SparseRCIntMatrix2D) {
            // exploit quick sparse mult
            // A*B = (B' * A')'
            if (C == null) {
                return B.multiply(this, null, alpha, beta, !transposeB, true).dice();
            } else {
                B.multiply(this, C.dice(), alpha, beta, !transposeB, true);
                return C;
            }
        }
        if (transposeB)
            return this.multiply(B.dice(), C, alpha, beta, transposeA, false);

        if (!(C is DenseIntMatrix2D))
            return super.multiply(B, C, alpha, beta, transposeA, transposeB);

        if (B.rows() != columnsA)
            throw new ArgumentError("Matrix2D inner dimensions must agree:" + this.toStringShort() + ", "
                    + B.toStringShort());
        if (C.rows() != rowsA || C.columns() != columnsB)
            throw new ArgumentError("Incompatibe result matrix: " + this.toStringShort() + ", "
                    + B.toStringShort() + ", " + C.toStringShort());
        if (this == C || B == C)
            throw new ArgumentError("Matrices must not be identical");

        int flops = 2 * rowsA * columnsA * columnsB;
        int noOfTasks = Math.min(flops / 30000, ConcurrencyUtils.getNumberOfThreads()); // each
        /* thread should process at least 30000 flops */
        bool splitB = (columnsB >= noOfTasks);
        int width = splitB ? columnsB : rowsA;
        noOfTasks = Math.min(width, noOfTasks);

        if (noOfTasks < 2) { //parallelization doesn't pay off (too much start up overhead)
            return this._zMultSequential(B, C, alpha, beta, transposeA, transposeB);
        }

        // set up concurrent tasks
        int span = width / noOfTasks;
        final List<Future> subTasks = new List<Future>(noOfTasks);
        for (int i = 0; i < noOfTasks; i++) {
            final int offset = i * span;
            if (i == noOfTasks - 1)
                span = width - span * i; // last span may be a bit larger

            final IntMatrix2D AA, BB, CC;
            if (splitB) {
                // split B aint columns into blocks
                AA = this;
                BB = B.part(0, offset, columnsA, span);
                CC = C.part(0, offset, rowsA, span);
            } else {
                // split A aint rows into blocks
                AA = this.part(offset, 0, span, columnsA);
                BB = B;
                CC = C.part(offset, 0, span, columnsB);
            }

            subTasks[i] = ConcurrencyUtils.submit(() {
              (AA as DenseIntMatrix2D)._zMultSequential(BB, CC, alpha, beta, transposeA, transposeB);
            });
        }

        ConcurrencyUtils.waitForCompletion(subTasks);
        return C;
    }

    int sum() {
        int sum = 0;
        if (_elements == null)
            throw new InternalError();
        final int zero = index(0, 0);
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int sum = 0;
                        int idx = zero + firstRow * _rowStride;
                        for (int r = firstRow; r < lastRow; r++) {
                            for (int i = idx, c = 0; c < _columns; c++) {
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
                    sum += futures[j].get() as int;
                }
            } on ExecutionException catch (ex) {
                ex.printStackTrace();
            } on InterruptedException catch (e) {
                e.printStackTrace();
            }
        } else {
            int idx = zero;
            for (int r = 0; r < _rows; r++) {
                for (int i = idx, c = 0; c < _columns; c++) {
                    sum += _elements[i];
                    i += _columnStride;
                }
                idx += _rowStride;
            }
        }
        return sum;
    }

    IntMatrix2D _zMultSequential(IntMatrix2D B, IntMatrix2D C, int alpha, int beta, bool transposeA,
            bool transposeB) {
        if (transposeA) {
            return dice().multiply(B, C, alpha, beta, false, transposeB);
        }
        if (B is SparseIntMatrix2D || B is SparseRCIntMatrix2D || B is SparseCCIntMatrix2D) {
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
        int p = B.columns();
        if (C == null) {
            C = new DenseIntMatrix2D(rowsA, p);
        }
        if (!(B is DenseIntMatrix2D) || !(C is DenseIntMatrix2D)) {
            return super.multiply(B, C, alpha, beta, transposeA, transposeB);
        }
        if (B.rows() != columnsA) {
            throw new ArgumentError("Matrix2D inner dimensions must agree:" + toStringShort() + ", "
                    + B.toStringShort());
        }
        if (C.rows() != rowsA || C.columns() != p) {
            throw new ArgumentError("Incompatibel result matrix: " + toStringShort() + ", "
                    + B.toStringShort() + ", " + C.toStringShort());
        }
        if (this == C || B == C) {
            throw new ArgumentError("Matrices must not be identical");
        }

        DenseIntMatrix2D BB = B as DenseIntMatrix2D;
        DenseIntMatrix2D CC = C as DenseIntMatrix2D;
        final Int32List AElems = this._elements;
        final Int32List BElems = BB._elements;
        final Int32List CElems = CC._elements;
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
        int m_optimal = (BLOCK_SIZE - columnsA) / (columnsA + 1);
        if (m_optimal <= 0) {
            m_optimal = 1;
        }
        int blocks = rowsA / m_optimal;
        int rr = 0;
        if (rowsA % m_optimal != 0) {
            blocks++;
        }
        for (; --blocks >= 0;) {
            int jB = BB.index(0, 0);
            int indexA = index(rr, 0);
            int jC = CC.index(rr, 0);
            rr += m_optimal;
            if (blocks == 0) {
                m_optimal += rowsA - rr;
            }

            for (int j = p; --j >= 0;) {
                int iA = indexA;
                int iC = jC;
                for (int i = m_optimal; --i >= 0;) {
                    int kA = iA;
                    int kB = jB;
                    int s = 0;

                    // loop unrolled
                    kA -= cA;
                    kB -= rB;

                    for (int k = columnsA % 4; --k >= 0;) {
                        s += AElems[kA += cA] * BElems[kB += rB];
                    }
                    for (int k = columnsA / 4; --k >= 0;) {
                        s += AElems[kA += cA] * BElems[kB += rB] + AElems[kA += cA] * BElems[kB += rB]
                                + AElems[kA += cA] * BElems[kB += rB] + AElems[kA += cA] * BElems[kB += rB];
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

    bool _haveSharedCellsRaw(IntMatrix2D other) {
        if (other is SelectedDenseIntMatrix2D) {
            return this._elements == other._elements;
        } else if (other is DenseIntMatrix2D) {
            return this._elements == other._elements;
        }
        return false;
    }

    IntMatrix1D _like1D(int size, int zero, int stride) {
        return new DenseIntMatrix1D(size, this._elements, zero, stride, true);
    }

    IntMatrix2D _viewSelectionLike(Int32List rowOffsets, Int32List columnOffsets) {
        return new SelectedDenseIntMatrix2D(this._elements, rowOffsets, columnOffsets, 0);
    }
}
