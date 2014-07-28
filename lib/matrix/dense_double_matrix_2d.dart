/*
Copyright (C) 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
is hereby granted without fee, provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear in supporting documentation.
CERN makes no representations about the suitability of this software for any purpose.
It is provided "as is" without expressed or implied warranty.
 */
part of cern.colt.matrix;

/*import java.io.IOException;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import cern.colt.list.tdouble.DoubleArrayList;
import cern.colt.list.tint.IntArrayList;
import cern.colt.matrix.io.MatrixInfo;
import cern.colt.matrix.io.MatrixSize;
import cern.colt.matrix.io.MatrixVectorReader;
import cern.colt.matrix.tdcomplex.impl.DenseDComplexMatrix2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import edu.emory.mathcs.jtransforms.dct.DoubleDCT_2D;
import edu.emory.mathcs.jtransforms.dht.DoubleDHT_2D;
import edu.emory.mathcs.jtransforms.dst.DoubleDST_2D;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_2D;
import edu.emory.mathcs.util.ConcurrencyUtils;*/

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

    DoubleFFT_2D _fft2;

    DoubleDCT_2D _dct2;

    DoubleDST_2D _dst2;

    DoubleDHT_2D _dht2;

    List<double> _elements;

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
     * @throws IllegalArgumentException
     *             if
     *             <tt>for any 1 &lt;= row &lt; values.length: values[row].length != values[row-1].length</tt>
     *             .
     */
    DenseDoubleMatrix2D(List<List<double>> values) {
        this(values.length, values.length == 0 ? 0 : values[0].length);
        assign(values);
    }

    /**
     * Constructs a matrix with a given number of rows and columns. All entries
     * are initially <tt>0</tt>.
     *
     * @param rows
     *            the number of rows the matrix shall have.
     * @param columns
     *            the number of columns the matrix shall have.
     * @throws IllegalArgumentException
     *             if
     *             <tt>rows<0 || columns<0 || (double)columns*rows > Integer.MAX_VALUE</tt>
     *             .
     */
    DenseDoubleMatrix2D(int rows, int columns) {
        _setUp(rows, columns);
        this._elements = new List<double>(rows * columns);
    }

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
     * @throws IllegalArgumentException
     *             if
     *             <tt>rows<0 || columns<0 || (double)columns*rows > Integer.MAX_VALUE</tt>
     *             or flip's are illegal.
     */
    DenseDoubleMatrix2D(int rows, int columns, List<double> elements, int rowZero, int columnZero, int rowStride,
            int columnStride, bool isView) {
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
    DenseDoubleMatrix2D(MatrixVectorReader reader) {//throws IOException {
        MatrixInfo info;
        if (reader.hasInfo())
            info = reader.readMatrixInfo();
        else
            info = new MatrixInfo(true, MatrixInfo.MatrixField.Real, MatrixInfo.MatrixSymmetry.General);

        if (info.isPattern())
            throw new UnsupportedOperationException("Pattern matrices are not supported");
        if (info.isDense())
            throw new UnsupportedOperationException("Dense matrices are not supported");
        if (info.isComplex())
            throw new UnsupportedOperationException("Complex matrices are not supported");

        MatrixSize size = reader.readMatrixSize(info);
        _setUp(size.numRows(), size.numColumns());
        this._elements = new List<double>(_rows * _columns);
        int numEntries = size.numEntries();
        List<int> columnIndexes = new List<int>(numEntries);
        List<int> rowIndexes = new List<int>(numEntries);
        List<double> values = new List<double>(numEntries);
        reader.readCoordinate(rowIndexes, columnIndexes, values);
        for (int i = 0; i < numEntries; i++) {
            setQuick(rowIndexes[i], columnIndexes[i], values[i]);
        }
        if (info.isSymmetric()) {
            for (int i = 0; i < numEntries; i++) {
                if (rowIndexes[i] != columnIndexes[i]) {
                    setQuick(columnIndexes[i], rowIndexes[i], values[i]);
                }
            }
        } else if (info.isSkewSymmetric()) {
            for (int i = 0; i < numEntries; i++) {
                if (rowIndexes[i] != columnIndexes[i]) {
                    setQuick(columnIndexes[i], rowIndexes[i], -values[i]);
                }
            }
        }
    }

    double aggregate(DoubleDoubleFunction aggr,
            DoubleFunction f) {
        if (size() == 0)
            return Double.NaN;
        final int zero = index(0, 0) as int;
        double a = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = _rows - j * k;
                final int lastRow = (j == (nthreads - 1)) ? 0 : firstRow - k;
                futures[j] = ConcurrencyUtils.submit(() {
                        double a = f.apply(_elements[zero + (firstRow - 1) * _rowStride + (_columns - 1) * _columnStride]);
                        int d = 1;
                        for (int r = firstRow; --r >= lastRow;) {
                            int ridx = zero + r * _rowStride;
                            for (int c = _columns - d; --c >= 0;) {
                                a = aggr.apply(a, f.apply(_elements[ridx + c * _columnStride]));
                            }
                            d = 0;
                        }
                        return a;
                });
            }
            a = ConcurrencyUtils.waitForCompletion(futures, aggr);
        } else {
            a = f.apply(_elements[zero + (_rows - 1) * _rowStride + (_columns - 1) * _columnStride]);
            int d = 1;
            for (int r = _rows; --r >= 0;) {
                int ridx = zero + r * _rowStride;
                for (int c = _columns - d; --c >= 0;) {
                    a = aggr.apply(a, f.apply(_elements[ridx + c * _columnStride]));
                }
                d = 0;
            }
        }
        return a;
    }

    double aggregate(DoubleDoubleFunction aggr,
            DoubleFunction f, DoubleProcedure cond) {
        if (size() == 0)
            return Double.NaN;
        final int zero = index(0, 0) as int;
        double a = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        double elem = _elements[zero + firstRow * _rowStride];
                        double a = 0;
                        if (cond.apply(elem) == true) {
                            a = f.apply(elem);
                        }
                        int d = 1;
                        for (int r = firstRow; r < lastRow; r++) {
                            for (int c = d; c < _columns; c++) {
                                elem = _elements[zero + r * _rowStride + c * _columnStride];
                                if (cond.apply(elem) == true) {
                                    a = aggr.apply(a, f.apply(elem));
                                }
                            }
                            d = 0;
                        }
                        return a;
                });
            }
            a = ConcurrencyUtils.waitForCompletion(futures, aggr);
        } else {
            double elem = _elements[zero];
            if (cond.apply(elem) == true) {
                a = f.apply(_elements[zero]);
            }
            int d = 1; // first cell already done
            for (int r = 0; r < _rows; r++) {
                for (int c = d; c < _columns; c++) {
                    elem = _elements[zero + r * _rowStride + c * _columnStride];
                    if (cond.apply(elem) == true) {
                        a = aggr.apply(a, f.apply(elem));
                    }
                }
                d = 0;
            }
        }
        return a;
    }

    double aggregate(DoubleDoubleFunction aggr,
            DoubleFunction f, final IntArrayList rowList, final IntArrayList columnList) {
        if (size() == 0)
            return Double.NaN;
        final int zero = index(0, 0) as int;
        final int size = rowList.size();
        final List<int> rowElements = rowList.elements();
        final List<int> columnElements = columnList.elements();
        double a = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, size);
            List<Future> futures = new List<Future>(nthreads);
            int k = size / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstIdx = j * k;
                final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        double a = f.apply(_elements[zero + rowElements[firstIdx] * _rowStride + columnElements[firstIdx]
                                * _columnStride]);
                        double elem;
                        for (int i = firstIdx + 1; i < lastIdx; i++) {
                            elem = _elements[zero + rowElements[i] * _rowStride + columnElements[i] * _columnStride];
                            a = aggr.apply(a, f.apply(elem));
                        }
                        return a;
                });
            }
            a = ConcurrencyUtils.waitForCompletion(futures, aggr);
        } else {
            double elem;
            a = f.apply(_elements[zero + rowElements[0] * _rowStride + columnElements[0] * _columnStride]);
            for (int i = 1; i < size; i++) {
                elem = _elements[zero + rowElements[i] * _rowStride + columnElements[i] * _columnStride];
                a = aggr.apply(a, f.apply(elem));
            }
        }
        return a;
    }

    double aggregate(final DoubleMatrix2D other, DoubleDoubleFunction aggr,
            DoubleDoubleFunction f) {
        if (!(other is DenseDoubleMatrix2D)) {
            return super.aggregate(other, aggr, f);
        }
        checkShape(other);
        if (size() == 0)
            return Double.NaN;
        final int zero = index(0, 0) as int;
        final int zeroOther = other.index(0, 0) as int;
        final int rowStrideOther = other.rowStride();
        final int colStrideOther = other.columnStride();
        final List<double> elementsOther = other.elements() as List<double>;
        double a = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        double a = f.apply(_elements[zero + firstRow * _rowStride], elementsOther[zeroOther + firstRow
                                * rowStrideOther]);
                        int d = 1;
                        for (int r = firstRow; r < lastRow; r++) {
                            for (int c = d; c < _columns; c++) {
                                a = aggr.apply(a, f.apply(_elements[zero + r * _rowStride + c * _columnStride],
                                        elementsOther[zeroOther + r * rowStrideOther + c * colStrideOther]));
                            }
                            d = 0;
                        }
                        return a;
                });
            }
            a = ConcurrencyUtils.waitForCompletion(futures, aggr);
        } else {
            int d = 1; // first cell already done
            a = f.apply(_elements[zero], elementsOther[zeroOther]);
            for (int r = 0; r < _rows; r++) {
                for (int c = d; c < _columns; c++) {
                    a = aggr.apply(a, f.apply(_elements[zero + r * _rowStride + c * _columnStride],
                            elementsOther[zeroOther + r * rowStrideOther + c * colStrideOther]));
                }
                d = 0;
            }
        }
        return a;
    }

    DoubleMatrix2D assign(DoubleFunction function) {
        final List<double> elems = this._elements;
        if (elems == null)
            throw new InternalError();
        final int zero = index(0, 0) as int;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            if (function is DoubleMult) { // x[i] =
                // mult*x[i]
                double multiplicator = function.multiplicator;
                if (multiplicator == 1)
                    return this;
                if (multiplicator == 0)
                    return assign(0);
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
                        if (function is DoubleMult) {
                            // x[i] = mult*x[i]
                            double multiplicator = function.multiplicator;
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
                                    elems[i] = function.apply(elems[i]);
                                    i += _columnStride;
                                }
                                idx += _rowStride;
                            }
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            int idx = zero + (_rows - 1) * _rowStride + (_columns - 1) * _columnStride;
            // specialization for speed
            if (function is DoubleMult) { // x[i] =
                // mult*x[i]
                double multiplicator = function.multiplicator;
                if (multiplicator == 1)
                    return this;
                if (multiplicator == 0)
                    return assign(0);
                for (int r = _rows; --r >= 0;) { // the general case
                    for (int i = idx, c = _columns; --c >= 0;) {
                        elems[i] *= multiplicator;
                        i -= _columnStride;
                    }
                    idx -= _rowStride;
                }
            } else { // the general case x[i] = f(x[i])
                for (int r = _rows; --r >= 0;) {
                    for (int i = idx, c = _columns; --c >= 0;) {
                        elems[i] = function.apply(elems[i]);
                        i -= _columnStride;
                    }
                    idx -= _rowStride;
                }
            }
        }
        return this;
    }

    DoubleMatrix2D assign(DoubleProcedure cond,
            DoubleFunction function) {
        final int zero = index(0, 0) as int;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        double elem;
                        int idx = zero + firstRow * _rowStride;
                        for (int r = firstRow; r < lastRow; r++) {
                            for (int i = idx, c = 0; c < _columns; c++) {
                                elem = _elements[i];
                                if (cond.apply(elem) == true) {
                                    _elements[i] = function.apply(elem);
                                }
                                i += _columnStride;
                            }
                            idx += _rowStride;
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            double elem;
            int idx = zero;
            for (int r = 0; r < _rows; r++) {
                for (int i = idx, c = 0; c < _columns; c++) {
                    elem = _elements[i];
                    if (cond.apply(elem) == true) {
                        _elements[i] = function.apply(elem);
                    }
                    i += _columnStride;
                }
                idx += _rowStride;
            }
        }
        return this;
    }

    DoubleMatrix2D assign(DoubleProcedure cond, final double value) {
        final int zero = index(0, 0) as int;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        double elem;
                        int idx = zero + firstRow * _rowStride;
                        for (int r = firstRow; r < lastRow; r++) {
                            for (int i = idx, c = 0; c < _columns; c++) {
                                elem = _elements[i];
                                if (cond.apply(elem) == true) {
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
            double elem;
            int idx = zero;
            for (int r = 0; r < _rows; r++) {
                for (int i = idx, c = 0; c < _columns; c++) {
                    elem = _elements[i];
                    if (cond.apply(elem) == true) {
                        _elements[i] = value;
                    }
                    i += _columnStride;
                }
                idx += _rowStride;
            }
        }
        return this;
    }

    DoubleMatrix2D assign(final double value) {
        final List<double> elems = this._elements;
        final int zero = index(0, 0) as int;
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

    DoubleMatrix2D assign(final List<double> values) {
        if (values.length != size())
            throw new IllegalArgumentException("Must have same length: length=" + values.length + " rows()*columns()="
                    + rows() * columns());
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if (this._isNoView) {
            System.arraycopy(values, 0, this._elements, 0, values.length);
        } else {
            final int zero = index(0, 0) as int;
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

    DoubleMatrix2D assign(final List<List<double>> values) {
        if (values.length != _rows)
            throw new ArgumentError("Must have same number of rows: rows=" + values.length + "rows()="
                    + rows());
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
                                List<double> currentRow = values[r];
                                if (currentRow.length != _columns)
                                    throw new IllegalArgumentException(
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
                    List<double> currentRow = values[r];
                    if (currentRow.length != _columns)
                        throw new IllegalArgumentException("Must have same number of columns in every row: columns="
                                + currentRow.length + "columns()=" + columns());
                    System.arraycopy(currentRow, 0, this._elements, i, _columns);
                    i += _columns;
                }
            }
        } else {
            final int zero = index(0, 0) as int;
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
                                List<double> currentRow = values[r];
                                if (currentRow.length != _columns)
                                    throw new IllegalArgumentException(
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
                    List<double> currentRow = values[r];
                    if (currentRow.length != _columns)
                        throw new IllegalArgumentException("Must have same number of columns in every row: columns="
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

    DoubleMatrix2D assign(final DoubleMatrix2D source) {
        // overriden for performance only
        if (!(source is DenseDoubleMatrix2D)) {
            super.assign(source);
            return this;
        }
        DenseDoubleMatrix2D other = double as DenseDoubleMatrix2D;
        if (other == this)
            return this; // nothing to do
        checkShape(other);
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if (this._isNoView && other._isNoView) { // quickest
            System.arraycopy(other._elements, 0, this._elements, 0, this._elements.length);
            return this;
        }
        if (_haveSharedCells(other)) {
            DoubleMatrix2D c = other.copy();
            if (!(c is DenseDoubleMatrix2D)) { // should not happen
                super.assign(other);
                return this;
            }
            other = c as DenseDoubleMatrix2D;
        }

        final List<double> elementsOther = other._elements;
        if (_elements == null || elementsOther == null)
            throw new InternalError();
        final int zeroOther = other.index(0, 0) as int;
        final int zero = index(0, 0) as int;
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
        } else {
            int idx = zero;
            int idxOther = zeroOther;
            for (int r = 0; r < _rows; r++) {
                for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
                    _elements[i] = elementsOther[j];
                    i += _columnStride;
                    j += columnStrideOther;
                }
                idx += _rowStride;
                idxOther += rowStrideOther;
            }
        }
        return this;
    }

    DoubleMatrix2D assign(final DoubleMatrix2D y, DoubleDoubleFunction function) {
        // overriden for performance only
        if (!(y is DenseDoubleMatrix2D)) {
            super.assign(y, function);
            return this;
        }
        DenseDoubleMatrix2D other = y as DenseDoubleMatrix2D;
        checkShape(y);
        final List<double> elementsOther = other._elements;
        if (_elements == null || elementsOther == null)
            throw new InternalError();
        final int zeroOther = other.index(0, 0) as int;
        final int zero = index(0, 0) as int;
        final int columnStrideOther = other._columnStride;
        final int rowStrideOther = other._rowStride;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            if (function is DoublePlusMultSecond) {
                double multiplicator = function.multiplicator;
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
                        if (function == cern.jet.math.tdouble.DoubleFunctions.mult) {
                            // x[i] = x[i]*y[i]
                            idx = zero + firstRow * _rowStride;
                            idxOther = zeroOther + firstRow * rowStrideOther;
                            for (int r = firstRow; r < lastRow; r++) {
                                for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
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
                                for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
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
                                    for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
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
                                    for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
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
                                    for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
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
                                for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
                                    _elements[i] = function.apply(_elements[i], elementsOther[j]);
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
            if (function == cern.jet.math.tdouble.DoubleFunctions.mult) {
                // x[i] = x[i] * y[i]
                idx = zero;
                idxOther = zeroOther;
                for (int r = 0; r < _rows; r++) {
                    for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
                        _elements[i] *= elementsOther[j];
                        i += _columnStride;
                        j += columnStrideOther;
                    }
                    idx += _rowStride;
                    idxOther += rowStrideOther;
                }
            } else if (function == cern.jet.math.tdouble.DoubleFunctions.div) {
                // x[i] = x[i] / y[i]
                idx = zero;
                idxOther = zeroOther;
                for (int r = 0; r < _rows; r++) {
                    for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
                        _elements[i] /= elementsOther[j];
                        i += _columnStride;
                        j += columnStrideOther;
                    }
                    idx += _rowStride;
                    idxOther += rowStrideOther;
                }
            } else if (function is DoublePlusMultSecond) {
                double multiplicator = function.multiplicator;
                if (multiplicator == 0) { // x[i] = x[i] + 0*y[i]
                    return this;
                } else if (multiplicator == 1) { // x[i] = x[i] + y[i]
                    idx = zero;
                    idxOther = zeroOther;
                    for (int r = 0; r < _rows; r++) {
                        for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
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
                        for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
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
                        for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
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
                    for (int i = idx, j = idxOther, c = 0; c < _columns; c++) {
                        _elements[i] = function.apply(_elements[i], elementsOther[j]);
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

    DoubleMatrix2D assign(final DoubleMatrix2D y,
            DoubleDoubleFunction function, IntArrayList rowList,
            IntArrayList columnList) {
        checkShape(y);
        final int size = rowList.size();
        final List<int> rowElements = rowList.elements();
        final List<int> columnElements = columnList.elements();
        final List<double> elementsOther = y.elements() as List<double>;
        final int zeroOther = y.index(0, 0) as int;
        final int zero = index(0, 0) as int;
        final int columnStrideOther = y.columnStride();
        final int rowStrideOther = y.rowStride();
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
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
                            _elements[idx] = function.apply(_elements[idx], elementsOther[idxOther]);
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
                _elements[idx] = function.apply(_elements[idx], elementsOther[idxOther]);
            }
        }
        return this;
    }

    DoubleMatrix2D assign(List<float> values) {
        if (values.length != size())
            throw new ArgumentError("Must have same length: length=" + values.length + "rows()*columns()="
                    + rows() * columns());
        final int zero = index(0, 0) as int;
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
                        int idxOther = firstRow * _columns;
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
        return this;
    }

    int cardinality() {
        int cardinality = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        final int zero = index(0, 0) as int;
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
                        return Integer.valueOf(cardinality);
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

    /**
     * Computes the 2D discrete cosine transform (DCT-II) of this matrix.
     *
     * @param scale
     *            if true then scaling is performed
     *
     */
    void dct2(bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        if (_dct2 == null) {
            _dct2 = new DoubleDCT_2D(_rows, _columns);
        }
        if (_isNoView == true) {
            _dct2.forward(_elements, scale);
        } else {
            DoubleMatrix2D copy = this.copy();
            _dct2.forward(copy.elements() as List<double>, scale);
            this.assign(copy.elements() as List<double>);
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the discrete cosine transform (DCT-II) of each column of this
     * matrix.
     *
     * @param scale
     *            if true then scaling is performed
     *
     */
    void dctColumns(final bool scale) {
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_2Threads(Integer.MAX_VALUE);
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_4Threads(Integer.MAX_VALUE);
            nthreads = Math.min(nthreads, _columns);
            List<Future> futures = new List<Future>(nthreads);
            int k = _columns / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstColumn = j * k;
                final int lastColumn = (j == nthreads - 1) ? _columns : firstColumn + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        for (int c = firstColumn; c < lastColumn; c++) {
                            (viewColumn(c) as DenseDoubleMatrix1D).dct(scale);
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
            ConcurrencyUtils.resetThreadsBeginN_FFT();
        } else {
            for (int c = 0; c < _columns; c++) {
                (viewColumn(c) as DenseDoubleMatrix1D).dct(scale);
            }
        }
    }

    /**
     * Computes the discrete cosine transform (DCT-II) of each row of this
     * matrix.
     *
     * @param scale
     *            if true then scaling is performed
     *
     */
    void dctRows(final bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_2Threads(Integer.MAX_VALUE);
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_4Threads(Integer.MAX_VALUE);
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        for (int r = firstRow; r < lastRow; r++) {
                            (viewRow(r) as DenseDoubleMatrix1D).dct(scale);
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
            ConcurrencyUtils.resetThreadsBeginN_FFT();
        } else {
            for (int r = 0; r < _rows; r++) {
                (viewRow(r) as DenseDoubleMatrix1D).dct(scale);
            }
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the 2D discrete Hartley transform (DHT) of this matrix.
     *
     */
    void dht2() {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        if (_dht2 == null) {
            _dht2 = new DoubleDHT_2D(_rows, _columns);
        }
        if (_isNoView == true) {
            _dht2.forward(_elements);
        } else {
            DoubleMatrix2D copy = this.copy();
            _dht2.forward(copy.elements() as List<double>);
            this.assign(copy.elements() as List<double>);
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the discrete Hartley transform (DHT) of each column of this
     * matrix.
     *
     */
    void dhtColumns() {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_2Threads(Integer.MAX_VALUE);
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_4Threads(Integer.MAX_VALUE);
            nthreads = Math.min(nthreads, _columns);
            List<Future> futures = new List<Future>(nthreads);
            int k = _columns / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstColumn = j * k;
                final int lastColumn = (j == nthreads - 1) ? _columns : firstColumn + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        for (int c = firstColumn; c < lastColumn; c++) {
                            (viewColumn(c) as DenseDoubleMatrix1D).dht();
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
            ConcurrencyUtils.resetThreadsBeginN_FFT();
        } else {
            for (int c = 0; c < _columns; c++) {
                (viewColumn(c) as DenseDoubleMatrix1D).dht();
            }
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the discrete Hartley transform (DHT) of each row of this matrix.
     *
     */
    void dhtRows() {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_2Threads(Integer.MAX_VALUE);
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_4Threads(Integer.MAX_VALUE);
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        for (int r = firstRow; r < lastRow; r++) {
                            (viewRow(r) as DenseDoubleMatrix1D).dht();
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
            ConcurrencyUtils.resetThreadsBeginN_FFT();
        } else {
            for (int r = 0; r < _rows; r++) {
                (viewRow(r) as DenseDoubleMatrix1D).dht();
            }
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the 2D discrete sine transform (DST-II) of this matrix.
     *
     * @param scale
     *            if true then scaling is performed
     *
     */
    void dst2(bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        if (_dst2 == null) {
            _dst2 = new DoubleDST_2D(_rows, _columns);
        }
        if (_isNoView == true) {
            _dst2.forward(_elements, scale);
        } else {
            DoubleMatrix2D copy = this.copy();
            _dst2.forward(copy.elements() as List<double>, scale);
            this.assign(copy.elements() as List<double>);
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the discrete sine transform (DST-II) of each column of this
     * matrix.
     *
     * @param scale
     *            if true then scaling is performed
     *
     */
    void dstColumns(final bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_2Threads(Integer.MAX_VALUE);
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_4Threads(Integer.MAX_VALUE);
            nthreads = Math.min(nthreads, _columns);
            List<Future> futures = new List<Future>(nthreads);
            int k = _columns / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstColumn = j * k;
                final int lastColumn = (j == nthreads - 1) ? _columns : firstColumn + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        for (int c = firstColumn; c < lastColumn; c++) {
                            (viewColumn(c) as DenseDoubleMatrix1D).dst(scale);
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
            ConcurrencyUtils.resetThreadsBeginN_FFT();
        } else {
            for (int c = 0; c < _columns; c++) {
                (viewColumn(c) as DenseDoubleMatrix1D).dst(scale);
            }
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the discrete sine transform (DST-II) of each row of this matrix.
     *
     * @param scale
     *            if true then scaling is performed
     *
     */
    void dstRows(final bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_2Threads(Integer.MAX_VALUE);
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_4Threads(Integer.MAX_VALUE);
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        for (int r = firstRow; r < lastRow; r++) {
                            (viewRow(r) as DenseDoubleMatrix1D).dst(scale);
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
            ConcurrencyUtils.resetThreadsBeginN_FFT();
        } else {
            for (int r = 0; r < _rows; r++) {
                (viewRow(r) as DenseDoubleMatrix1D).dst(scale);
            }
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    List<double> elements() {
        return _elements;
    }

    /**
     * Computes the 2D discrete Fourier transform (DFT) of this matrix. The
     * physical layout of the output data is as follows:
     *
     * <pre>
     * this[k1][2*k2] = Re[k1][k2] = Re[rows-k1][columns-k2],
     * this[k1][2*k2+1] = Im[k1][k2] = -Im[rows-k1][columns-k2],
     *       0&lt;k1&lt;rows, 0&lt;k2&lt;columns/2,
     * this[0][2*k2] = Re[0][k2] = Re[0][columns-k2],
     * this[0][2*k2+1] = Im[0][k2] = -Im[0][columns-k2],
     *       0&lt;k2&lt;columns/2,
     * this[k1][0] = Re[k1][0] = Re[rows-k1][0],
     * this[k1][1] = Im[k1][0] = -Im[rows-k1][0],
     * this[rows-k1][1] = Re[k1][columns/2] = Re[rows-k1][columns/2],
     * this[rows-k1][0] = -Im[k1][columns/2] = Im[rows-k1][columns/2],
     *       0&lt;k1&lt;rows/2,
     * this[0][0] = Re[0][0],
     * this[0][1] = Re[0][columns/2],
     * this[rows/2][0] = Re[rows/2][0],
     * this[rows/2][1] = Re[rows/2][columns/2]
     * </pre>
     *
     * This method computes only half of the elements of the real transform. The
     * other half satisfies the symmetry condition. If you want the full real
     * forward transform, use <code>getFft2</code>. To get back the original
     * data, use <code>ifft2</code>.
     *
     * @throws IllegalArgumentException
     *             if the row size or the column size of this matrix is not a
     *             power of 2 number.
     *
     */
    void fft2() {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        if (_fft2 == null) {
            _fft2 = new DoubleFFT_2D(_rows, _columns);
        }
        if (_isNoView == true) {
            _fft2.realForward(_elements);
        } else {
            DoubleMatrix2D copy = this.copy();
            _fft2.realForward(copy.elements() as List<double>);
            this.assign(copy.elements() as List<double>);
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    DoubleMatrix2D forEachNonZero(IntIntDoubleFunction function) {
        final int zero = index(0, 0) as int;
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
                                double value = _elements[i];
                                if (value != 0) {
                                    _elements[i] = function.apply(r, c, value);
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
                    double value = _elements[i];
                    if (value != 0) {
                        _elements[i] = function.apply(r, c, value);
                    }
                    i += _columnStride;
                }
                idx += _rowStride;
            }
        }
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
    DenseColumnDoubleMatrix2D getColumnMajor() {
        DenseColumnDoubleMatrix2D R = new DenseColumnDoubleMatrix2D(_rows, _columns);
        final int zeroR = R.index(0, 0) as int;
        final int rowStrideR = R.rowStride();
        final int columnStrideR = R.columnStride();
        final List<double> elementsR = R.elements();
        final int zero = index(0, 0) as int;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == (nthreads - 1)) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int idx = zero + (firstRow - 1) * _rowStride;
                        int idxR = zeroR + (firstRow - 1) * rowStrideR;
                        for (int r = firstRow; r < lastRow; r++) {
                            for (int i = idx, j = idxR, c = 0; r < _columns; c++) {
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
                for (int i = idx, j = idxR, c = 0; r < _columns; c++) {
                    elementsR[j] = _elements[i];
                    i += _rowStride;
                    j += rowStrideR;
                }
                idx += _columnStride;
                idxR += columnStrideR;
            }
        }
        return R;
    }

    /**
     * Returns new complex matrix which is the 2D discrete Fourier transform
     * (DFT) of this matrix.
     *
     * @return the 2D discrete Fourier transform (DFT) of this matrix.
     *
     */
    DenseDComplexMatrix2D getFft2() {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        if (_fft2 == null) {
            _fft2 = new DoubleFFT_2D(_rows, _columns);
        }
        final List<double> elementsA;
        if (_isNoView == true) {
            elementsA = _elements;
        } else {
            elementsA = this.copy().elements() as List<double>;
        }
        DenseDComplexMatrix2D C = new DenseDComplexMatrix2D(_rows, _columns);
        final List<double> elementsC = (C).elements();
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        for (int r = firstRow; r < lastRow; r++) {
                            System.arraycopy(elementsA, r * _columns, elementsC, r * _columns, _columns);
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            for (int r = 0; r < _rows; r++) {
                System.arraycopy(elementsA, r * _columns, elementsC, r * _columns, _columns);
            }
        }
        _fft2.realForwardFull(elementsC);
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
        return C;
    }

    /**
     * Returns new complex matrix which is the discrete Fourier transform (DFT)
     * of each column of this matrix.
     *
     * @return the discrete Fourier transform (DFT) of each column of this
     *         matrix.
     */
    DenseDComplexMatrix2D getFftColumns() {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        final DenseDComplexMatrix2D C = new DenseDComplexMatrix2D(_rows, _columns);
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_2Threads(Integer.MAX_VALUE);
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_4Threads(Integer.MAX_VALUE);
            nthreads = Math.min(nthreads, _columns);
            List<Future> futures = new List<Future>(nthreads);
            int k = _columns / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstColumn = j * k;
                final int lastColumn = (j == nthreads - 1) ? _columns : firstColumn + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        for (int c = firstColumn; c < lastColumn; c++) {
                            C.viewColumn(c).assign((viewColumn(c) as DenseDoubleMatrix1D).getFft());
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
            ConcurrencyUtils.resetThreadsBeginN_FFT();
        } else {
            for (int c = 0; c < _columns; c++) {
                C.viewColumn(c).assign((viewColumn(c) as DenseDoubleMatrix1D).getFft());
            }
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
        return C;
    }

    /**
     * Returns new complex matrix which is the discrete Fourier transform (DFT)
     * of each row of this matrix.
     *
     * @return the discrete Fourier transform (DFT) of each row of this matrix.
     */
    DenseDComplexMatrix2D getFftRows() {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        final DenseDComplexMatrix2D C = new DenseDComplexMatrix2D(_rows, _columns);
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_2Threads(Integer.MAX_VALUE);
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_4Threads(Integer.MAX_VALUE);
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        for (int r = firstRow; r < lastRow; r++) {
                            C.viewRow(r).assign((viewRow(r) as DenseDoubleMatrix1D).getFft());
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
            ConcurrencyUtils.resetThreadsBeginN_FFT();
        } else {
            for (int r = 0; r < _rows; r++) {
                C.viewRow(r).assign((viewRow(r) as DenseDoubleMatrix1D).getFft());
            }
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
        return C;
    }

    /**
     * Returns new complex matrix which is the 2D inverse of the discrete
     * Fourier transform (IDFT) of this matrix.
     *
     * @return the 2D inverse of the discrete Fourier transform (IDFT) of this
     *         matrix.
     */
    DenseDComplexMatrix2D getIfft2(bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        DenseDComplexMatrix2D C = new DenseDComplexMatrix2D(_rows, _columns);
        final List<double> elementsC = (C).elements();
        final List<double> elementsA;
        if (_isNoView == true) {
            elementsA = _elements;
        } else {
            elementsA = this.copy().elements() as List<double>;
        }
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow;
                if (j == nthreads - 1) {
                    lastRow = _rows;
                } else {
                    lastRow = firstRow + k;
                }
                futures[j] = ConcurrencyUtils.submit(() {
                        for (int r = firstRow; r < lastRow; r++) {
                            System.arraycopy(elementsA, r * _columns, elementsC, r * _columns, _columns);
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            for (int r = 0; r < _rows; r++) {
                System.arraycopy(elementsA, r * _columns, elementsC, r * _columns, _columns);
            }
        }
        if (_fft2 == null) {
            _fft2 = new DoubleFFT_2D(_rows, _columns);
        }
        _fft2.realInverseFull(elementsC, scale);
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
        return C;
    }

    /**
     * Returns new complex matrix which is the inverse of the discrete Fourier
     * transform (IDFT) of each column of this matrix.
     *
     * @return the inverse of the discrete Fourier transform (IDFT) of each
     *         column of this matrix.
     */
    DenseDComplexMatrix2D getIfftColumns(final bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        final DenseDComplexMatrix2D C = new DenseDComplexMatrix2D(_rows, _columns);
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_2Threads(Integer.MAX_VALUE);
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_4Threads(Integer.MAX_VALUE);
            nthreads = Math.min(nthreads, _columns);
            List<Future> futures = new List<Future>(nthreads);
            int k = _columns / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstColumn = j * k;
                final int lastColumn = (j == nthreads - 1) ? _columns : firstColumn + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        for (int c = firstColumn; c < lastColumn; c++) {
                            C.viewColumn(c).assign((viewColumn(c) as DenseDoubleMatrix1D).getIfft(scale));
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
            ConcurrencyUtils.resetThreadsBeginN_FFT();
        } else {
            for (int c = 0; c < _columns; c++) {
                C.viewColumn(c).assign((viewColumn(c) as DenseDoubleMatrix1D).getIfft(scale));
            }
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
        return C;
    }

    /**
     * Returns new complex matrix which is the inverse of the discrete Fourier
     * transform (IDFT) of each row of this matrix.
     *
     * @return the inverse of the discrete Fourier transform (IDFT) of each row
     *         of this matrix.
     */
    DenseDComplexMatrix2D getIfftRows(final bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        final DenseDComplexMatrix2D C = new DenseDComplexMatrix2D(_rows, _columns);
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_2Threads(Integer.MAX_VALUE);
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_4Threads(Integer.MAX_VALUE);
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        for (int r = firstRow; r < lastRow; r++) {
                            C.viewRow(r).assign((viewRow(r) as DenseDoubleMatrix1D).getIfft(scale));
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
            ConcurrencyUtils.resetThreadsBeginN_FFT();
        } else {
            for (int r = 0; r < _rows; r++) {
                C.viewRow(r).assign((viewRow(r) as DenseDoubleMatrix1D).getIfft(scale));
            }
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
        return C;
    }

    List<double> getMaxLocation() {
        int rowLocation = 0;
        int columnLocation = 0;
        final int zero = index(0, 0) as int;
        double maxValue = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            List<List<double>> results = new List<List<double>>(nthreads);//[2];
            int k = _rows / nthreads;
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
                        return [ maxValue, rowLocation, colLocation ];
                });
            }
            try {
                for (int j = 0; j < nthreads; j++) {
                    results[j] = futures[j].get() as List<double>;
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
        }
        return [ maxValue, rowLocation, columnLocation ];
    }

    List<double> getMinLocation() {
        int rowLocation = 0;
        int columnLocation = 0;
        final int zero = index(0, 0) as int;
        double minValue = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            List<List<double>> results = new List<List<double>>(nthreads);//[2];
            int k = _rows / nthreads;
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
                        return [ minValue, rowLocation, columnLocation ];
                });
            }
            try {
                for (int j = 0; j < nthreads; j++) {
                    results[j] = futures[j].get() as List<double>;
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
        }
        return [ minValue, rowLocation, columnLocation ];
    }

    void getNegativeValues(final IntArrayList rowList, final IntArrayList columnList,
            final DoubleArrayList valueList) {
        rowList.clear();
        columnList.clear();
        valueList.clear();
        int idx = index(0, 0) as int;
        for (int r = 0; r < _rows; r++) {
            for (int i = idx, c = 0; c < _columns; c++) {
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

    void getNonZeros(final IntArrayList rowList, final IntArrayList columnList, final DoubleArrayList valueList) {
        rowList.clear();
        columnList.clear();
        valueList.clear();
        int idx = index(0, 0) as int;
        for (int r = 0; r < _rows; r++) {
            for (int i = idx, c = 0; c < _columns; c++) {
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

    void getPositiveValues(final IntArrayList rowList, final IntArrayList columnList,
            final DoubleArrayList valueList) {
        rowList.clear();
        columnList.clear();
        valueList.clear();
        int idx = index(0, 0) as int;
        for (int r = 0; r < _rows; r++) {
            for (int i = idx, c = 0; c < _columns; c++) {
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

    double getQuick(int row, int column) {
        return _elements[_rowZero + row * _rowStride + _columnZero + column * _columnStride];
    }

    /**
     * Computes the 2D inverse of the discrete cosine transform (DCT-III) of
     * this matrix.
     *
     * @param scale
     *            if true then scaling is performed
     *
     */
    void idct2(bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        if (_dct2 == null) {
            _dct2 = new DoubleDCT_2D(_rows, _columns);
        }
        if (_isNoView == true) {
            _dct2.inverse(_elements, scale);
        } else {
            DoubleMatrix2D copy = this.copy();
            _dct2.inverse(copy.elements() as List<double>, scale);
            this.assign(copy.elements() as List<double>);
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the inverse of the discrete cosine transform (DCT-III) of each
     * column of this matrix.
     *
     * @param scale
     *            if true then scaling is performed
     *
     */
    void idctColumns(final bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_2Threads(Integer.MAX_VALUE);
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_4Threads(Integer.MAX_VALUE);
            nthreads = Math.min(nthreads, _columns);
            List<Future> futures = new List<Future>(nthreads);
            int k = _columns / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstColumn = j * k;
                final int lastColumn = (j == nthreads - 1) ? _columns : firstColumn + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        for (int c = firstColumn; c < lastColumn; c++) {
                            (viewColumn(c) as DenseDoubleMatrix1D).idct(scale);
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
            ConcurrencyUtils.resetThreadsBeginN_FFT();
        } else {
            for (int c = 0; c < _columns; c++) {
                (viewColumn(c) as DenseDoubleMatrix1D).idct(scale);
            }
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the inverse of the discrete cosine transform (DCT-III) of each
     * row of this matrix.
     *
     * @param scale
     *            if true then scaling is performed
     *
     */
    void idctRows(final bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_2Threads(Integer.MAX_VALUE);
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_4Threads(Integer.MAX_VALUE);
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        for (int r = firstRow; r < lastRow; r++) {
                            (viewRow(r) as DenseDoubleMatrix1D).idct(scale);
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
            ConcurrencyUtils.resetThreadsBeginN_FFT();
        } else {
            for (int r = 0; r < _rows; r++) {
                (viewRow(r) as DenseDoubleMatrix1D).idct(scale);
            }
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the 2D inverse of the discrete Hartley transform (IDHT) of this
     * matrix.
     *
     * @param scale
     *            if true then scaling is performed
     *
     */
    void idht2(bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        if (_dht2 == null) {
            _dht2 = new DoubleDHT_2D(_rows, _columns);
        }
        if (_isNoView == true) {
            _dht2.inverse(_elements, scale);
        } else {
            DoubleMatrix2D copy = this.copy();
            _dht2.inverse(copy.elements() as List<double>, scale);
            this.assign(copy.elements() as List<double>);
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the inverse of the discrete Hartley transform (IDHT) of each
     * column of this matrix.
     *
     * @param scale
     *            if true then scaling is performed
     *
     */
    void idhtColumns(final bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_2Threads(Integer.MAX_VALUE);
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_4Threads(Integer.MAX_VALUE);
            nthreads = Math.min(nthreads, _columns);
            List<Future> futures = new List<Future>(nthreads);
            int k = _columns / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstColumn = j * k;
                final int lastColumn = (j == nthreads - 1) ? _columns : firstColumn + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        for (int c = firstColumn; c < lastColumn; c++) {
                            (viewColumn(c) as DenseDoubleMatrix1D).idht(scale);
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
            ConcurrencyUtils.resetThreadsBeginN_FFT();
        } else {
            for (int c = 0; c < _columns; c++) {
                (viewColumn(c) as DenseDoubleMatrix1D).idht(scale);
            }
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the inverse of the discrete Hartley transform (IDHT) of each row
     * of this matrix.
     *
     * @param scale
     *            if true then scaling is performed
     *
     */
    void idhtRows(final bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_2Threads(Integer.MAX_VALUE);
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_4Threads(Integer.MAX_VALUE);
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        for (int r = firstRow; r < lastRow; r++) {
                            (viewRow(r) as DenseDoubleMatrix1D).idht(scale);
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
            ConcurrencyUtils.resetThreadsBeginN_FFT();
        } else {
            for (int r = 0; r < _rows; r++) {
                (viewRow(r) as DenseDoubleMatrix1D).idht(scale);
            }
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the 2D inverse of the discrete sine transform (DST-III) of this
     * matrix.
     *
     * @param scale
     *            if true then scaling is performed
     *
     */
    void idst2(bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        if (_dst2 == null) {
            _dst2 = new DoubleDST_2D(_rows, _columns);
        }
        if (_isNoView == true) {
            _dst2.inverse(_elements, scale);
        } else {
            DoubleMatrix2D copy = this.copy();
            _dst2.inverse(copy.elements() as List<double>, scale);
            this.assign(copy.elements() as List<double>);
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the inverse of the discrete sine transform (DST-III) of each
     * column of this matrix.
     *
     * @param scale
     *            if true then scaling is performed
     *
     */
    void idstColumns(final bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_2Threads(Integer.MAX_VALUE);
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_4Threads(Integer.MAX_VALUE);
            nthreads = Math.min(nthreads, _columns);
            List<Future> futures = new List<Future>(nthreads);
            int k = _columns / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstColumn = j * k;
                final int lastColumn = (j == nthreads - 1) ? _columns : firstColumn + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        for (int c = firstColumn; c < lastColumn; c++) {
                            (viewColumn(c) as DenseDoubleMatrix1D).idst(scale);
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
            ConcurrencyUtils.resetThreadsBeginN_FFT();
        } else {
            for (int c = 0; c < _columns; c++) {
                (viewColumn(c) as DenseDoubleMatrix1D).idst(scale);
            }
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the inverse of the discrete sine transform (DST-III) of each row
     * of this matrix.
     *
     * @param scale
     *            if true then scaling is performed
     *
     */
    void idstRows(final bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_2Threads(Integer.MAX_VALUE);
            ConcurrencyUtils.setThreadsBeginN_1D_FFT_4Threads(Integer.MAX_VALUE);
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        for (int r = firstRow; r < lastRow; r++) {
                            (viewRow(r) as DenseDoubleMatrix1D).idst(scale);
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
            ConcurrencyUtils.resetThreadsBeginN_FFT();
        } else {
            for (int r = 0; r < _rows; r++) {
                (viewRow(r) as DenseDoubleMatrix1D).idst(scale);
            }
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the 2D inverse of the discrete Fourier transform (IDFT) of this
     * matrix. The physical layout of the input data has to be as follows:
     *
     * <pre>
     * this[k1][2*k2] = Re[k1][k2] = Re[rows-k1][columns-k2],
     * this[k1][2*k2+1] = Im[k1][k2] = -Im[rows-k1][columns-k2],
     *       0&lt;k1&lt;rows, 0&lt;k2&lt;columns/2,
     * this[0][2*k2] = Re[0][k2] = Re[0][columns-k2],
     * this[0][2*k2+1] = Im[0][k2] = -Im[0][columns-k2],
     *       0&lt;k2&lt;columns/2,
     * this[k1][0] = Re[k1][0] = Re[rows-k1][0],
     * this[k1][1] = Im[k1][0] = -Im[rows-k1][0],
     * this[rows-k1][1] = Re[k1][columns/2] = Re[rows-k1][columns/2],
     * this[rows-k1][0] = -Im[k1][columns/2] = Im[rows-k1][columns/2],
     *       0&lt;k1&lt;rows/2,
     * this[0][0] = Re[0][0],
     * this[0][1] = Re[0][columns/2],
     * this[rows/2][0] = Re[rows/2][0],
     * this[rows/2][1] = Re[rows/2][columns/2]
     * </pre>
     *
     * This method computes only half of the elements of the real transform. The
     * other half satisfies the symmetry condition. If you want the full real
     * inverse transform, use <code>getIfft2</code>.
     *
     * @throws IllegalArgumentException
     *             if the row size or the column size of this matrix is not a
     *             power of 2 number.
     *
     * @param scale
     *            if true then scaling is performed
     *
     */
    void ifft2(bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        if (_fft2 == null) {
            _fft2 = new DoubleFFT_2D(_rows, _columns);
        }
        if (_isNoView == true) {
            _fft2.realInverse(_elements, scale);
        } else {
            DoubleMatrix2D copy = this.copy();
            _fft2.realInverse(copy.elements() as List<double>, scale);
            this.assign(copy.elements() as List<double>);
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    int index(int row, int column) {
        return _rowZero + row * _rowStride + _columnZero + column * _columnStride;
    }

    DoubleMatrix2D like(int rows, int columns) {
        return new DenseDoubleMatrix2D(rows, columns);
    }

    DoubleMatrix1D like1D(int size) {
        return new DenseDoubleMatrix1D(size);
    }

    void setQuick(int row, int column, double value) {
        _elements[_rowZero + row * _rowStride + _columnZero + column * _columnStride] = value;
    }

    List<List<double>> toArray() {
        final List<List<double>> values = new List<List<double>>(_rows);//[_columns];
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        final int zero = index(0, 0) as int;
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
                            List<double> currentRow = values[r];
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
                List<double> currentRow = values[r];
                for (int i = idx, c = 0; c < _columns; c++) {
                    currentRow[c] = _elements[i];
                    i += _columnStride;
                }
                idx += _rowStride;
            }
        }
        return values;
    }

    DoubleMatrix1D vectorize() {
        final DenseDoubleMatrix1D v = new DenseDoubleMatrix1D(size() as int);
        final int zero = index(0, 0) as int;
        final int zeroOther = v.index(0) as int;
        final int strideOther = v.stride();
        final List<double> elementsOther = v.elements();
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _columns);
            List<Future> futures = new List<Future>(nthreads);
            int k = _columns / nthreads;
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
        } else {
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
        }
        return v;
    }

    void zAssign8Neighbors(DoubleMatrix2D B, Double9Function function) {
        // 1. using only 4-5 out of the 9 cells in "function" is *not* the
        // limiting factor for performance.

        // 2. if the "function" would be hardwired into the innermost loop, a
        // speedup of 1.5-2.0 would be seen
        // but then the multi-purpose interface is gone...

        if (!(B is DenseDoubleMatrix2D)) {
            super.zAssign8Neighbors(B, function);
            return;
        }
        if (function == null)
            throw new NullPointerException("function must not be null.");
        checkShape(B);
        int r = _rows - 1;
        int c = _columns - 1;
        if (_rows < 3 || _columns < 3)
            return; // nothing to do

        DenseDoubleMatrix2D BB = B as DenseDoubleMatrix2D;
        int A_rs = _rowStride;
        int B_rs = BB._rowStride;
        int A_cs = _columnStride;
        int B_cs = BB._columnStride;
        List<double> elems = this._elements;
        List<double> B_elems = BB._elements;
        if (elems == null || B_elems == null)
            throw new InternalError();

        int A_index = index(1, 1) as int;
        int B_index = BB.index(1, 1) as int;
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

                B_elems[B11] = function.apply(a00, a01, a02, a10, a11, a12, a20, a21, a22);
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

    }

    DoubleMatrix1D zMult(final DoubleMatrix1D y, DoubleMatrix1D z, final double alpha, final double beta,
            final bool transposeA) {
        if (transposeA)
            return viewDice().zMult(y, z, alpha, beta, false);
        if (z == null) {
            z = new DenseDoubleMatrix1D(_rows);
        }
        if (!(y is DenseDoubleMatrix1D && z is DenseDoubleMatrix1D))
            return super.zMult(y, z, alpha, beta, transposeA);

        if (_columns != y.size() || _rows > z.size())
            throw new IllegalArgumentException("Incompatible args: " + toStringShort() + ", " + y.toStringShort()
                    + ", " + z.toStringShort());

        final List<double> elemsY = y.elements() as List<double>;
        final List<double> elemsZ = z.elements() as List<double>;
        if (_elements == null || elemsY == null || elemsZ == null)
            throw new InternalError();
        final int strideY = y.stride();
        final int strideZ = z.stride();
        final int zero = index(0, 0) as int;
        final int zeroY = y.index(0) as int;
        final int zeroZ = z.index(0) as int;
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
        } else {
            int idxZero = zero;
            int idxZeroZ = zeroZ;
            for (int r = 0; r < _rows; r++) {
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
        }
        return z;
    }

    DoubleMatrix2D zMult(final DoubleMatrix2D B, DoubleMatrix2D C, final double alpha, final double beta,
            final bool transposeA, final bool transposeB) {
        final int rowsA = _rows;
        final int columnsA = _columns;
        final int rowsB = B.rows();
        final int columnsB = B.columns();
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
        if (transposeA)
            return viewDice().zMult(B, C, alpha, beta, false, transposeB);
        if (B is SparseDoubleMatrix2D || B is SparseRCDoubleMatrix2D) {
            // exploit quick sparse mult
            // A*B = (B' * A')'
            if (C == null) {
                return B.zMult(this, null, alpha, beta, !transposeB, true).viewDice();
            } else {
                B.zMult(this, C.viewDice(), alpha, beta, !transposeB, true);
                return C;
            }
        }
        if (transposeB)
            return this.zMult(B.viewDice(), C, alpha, beta, transposeA, false);

        if (!(C is DenseDoubleMatrix2D))
            return super.zMult(B, C, alpha, beta, transposeA, transposeB);

        if (B.rows() != columnsA)
            throw new IllegalArgumentException("Matrix2D inner dimensions must agree:" + this.toStringShort() + ", "
                    + B.toStringShort());
        if (C.rows() != rowsA || C.columns() != columnsB)
            throw new IllegalArgumentException("Incompatibe result matrix: " + this.toStringShort() + ", "
                    + B.toStringShort() + ", " + C.toStringShort());
        if (this == C || B == C)
            throw new IllegalArgumentException("Matrices must not be identical");

        int flops = 2 * rowsA * columnsA * columnsB;
        int noOfTasks = Math.min(flops / 30000, ConcurrencyUtils.getNumberOfThreads()) as int; // each
        /* thread should process at least 30000 flops */
        bool splitB = (columnsB >= noOfTasks);
        int width = splitB ? columnsB : rowsA;
        noOfTasks = Math.min(width, noOfTasks);

        if (noOfTasks < 2) { //parallelization doesn't pay off (too much start up overhead)
            return this.zMultSequential(B, C, alpha, beta, transposeA, transposeB);
        }

        // set up concurrent tasks
        int span = width / noOfTasks;
        final List<Future> subTasks = new List<Future>(noOfTasks);
        for (int i = 0; i < noOfTasks; i++) {
            final int offset = i * span;
            if (i == noOfTasks - 1)
                span = width - span * i; // last span may be a bit larger

            final DoubleMatrix2D AA, BB, CC;
            if (splitB) {
                // split B along columns into blocks
                AA = this;
                BB = B.viewPart(0, offset, columnsA, span);
                CC = C.viewPart(0, offset, rowsA, span);
            } else {
                // split A along rows into blocks
                AA = this.viewPart(offset, 0, span, columnsA);
                BB = B;
                CC = C.viewPart(offset, 0, span, columnsB);
            }

            subTasks[i] = ConcurrencyUtils.submit(() {
                    (AA as DenseDoubleMatrix2D).zMultSequential(BB, CC, alpha, beta, transposeA, transposeB);
            });
        }

        ConcurrencyUtils.waitForCompletion(subTasks);
        return C;
    }

    double zSum() {
        double sum = 0;
        if (_elements == null)
            throw new InternalError();
        final int zero = index(0, 0) as int;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        double sum = 0;
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
                    sum += futures[j].get() as double;
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

    DoubleMatrix2D zMultSequential(DoubleMatrix2D B, DoubleMatrix2D C, double alpha, double beta,
            bool transposeA, bool transposeB) {
        if (transposeA)
            return viewDice().zMult(B, C, alpha, beta, false, transposeB);
        if (B is SparseDoubleMatrix2D || B is SparseRCDoubleMatrix2D
                || B is SparseCCDoubleMatrix2D) {
            // exploit quick sparse mult
            // A*B = (B' * A')'
            if (C == null) {
                return B.zMult(this, null, alpha, beta, !transposeB, true).viewDice();
            } else {
                B.zMult(this, C.viewDice(), alpha, beta, !transposeB, true);
                return C;
            }
        }
        if (transposeB)
            return this.zMult(B.viewDice(), C, alpha, beta, transposeA, false);

        int rowsA = _rows;
        int columnsA = _columns;
        int p = B.columns();
        if (C == null) {
            C = new DenseDoubleMatrix2D(rowsA, p);
        }
        if (!(B is DenseDoubleMatrix2D) || !(C is DenseDoubleMatrix2D))
            return super.zMult(B, C, alpha, beta, transposeA, transposeB);
        if (B.rows() != columnsA)
            throw new IllegalArgumentException("Matrix2D inner dimensions must agree:" + toStringShort() + ", "
                    + B.toStringShort());
        if (C.rows() != rowsA || C.columns() != p)
            throw new IllegalArgumentException("Incompatibel result matrix: " + toStringShort() + ", "
                    + B.toStringShort() + ", " + C.toStringShort());
        if (this == C || B == C)
            throw new IllegalArgumentException("Matrices must not be identical");

        DenseDoubleMatrix2D BB = B as DenseDoubleMatrix2D;
        DenseDoubleMatrix2D CC = C as DenseDoubleMatrix2D;
        final List<double> AElems = this._elements;
        final List<double> BElems = BB._elements;
        final List<double> CElems = CC._elements;
        if (AElems == null || BElems == null || CElems == null)
            throw new InternalError();

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
        if (m_optimal <= 0)
            m_optimal = 1;
        int blocks = rowsA / m_optimal;
        int rr = 0;
        if (rowsA % m_optimal != 0)
            blocks++;
        for (; --blocks >= 0;) {
            int jB = BB.index(0, 0) as int;
            int indexA = index(rr, 0) as int;
            int jC = CC.index(rr, 0) as int;
            rr += m_optimal;
            if (blocks == 0)
                m_optimal += rowsA - rr;

            for (int j = p; --j >= 0;) {
                int iA = indexA;
                int iC = jC;
                for (int i = m_optimal; --i >= 0;) {
                    int kA = iA;
                    int kB = jB;
                    double s = 0;

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

    bool _haveSharedCellsRaw(DoubleMatrix2D other) {
        if (other is SelectedDenseDoubleMatrix2D) {
            SelectedDenseDoubleMatrix2D otherMatrix = other as SelectedDenseDoubleMatrix2D;
            return this._elements == otherMatrix._elements;
        } else if (other is DenseDoubleMatrix2D) {
            DenseDoubleMatrix2D otherMatrix = other as DenseDoubleMatrix2D;
            return this._elements == otherMatrix._elements;
        }
        return false;
    }

    DoubleMatrix1D _like1D(int size, int zero, int stride) {
        return new DenseDoubleMatrix1D(size, this._elements, zero, stride, true);
    }

    DoubleMatrix2D _viewSelectionLike(List<int> rowOffsets, List<int> columnOffsets) {
        return new SelectedDenseDoubleMatrix2D(this._elements, rowOffsets, columnOffsets, 0);
    }
}
