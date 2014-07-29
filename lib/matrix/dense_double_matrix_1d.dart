/*
Copyright (C) 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
is hereby granted without fee, provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear in supporting documentation.
CERN makes no representations about the suitability of this software for any purpose.
It is provided "as is" without expressed or implied warranty.
 */
part of cern.colt.matrix;

/*import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import cern.colt.list.tdouble.DoubleArrayList;
import cern.colt.list.tint.IntArrayList;
import cern.colt.matrix.tdcomplex.impl.DenseDComplexMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tdouble.DoubleMatrix3D;
import edu.emory.mathcs.jtransforms.dct.DoubleDCT_1D;
import edu.emory.mathcs.jtransforms.dht.DoubleDHT_1D;
import edu.emory.mathcs.jtransforms.dst.DoubleDST_1D;
import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;
import edu.emory.mathcs.util.ConcurrencyUtils;*/

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
class DenseDoubleMatrix1D extends DoubleMatrix1D {

    DoubleFFT_1D _fft;

    DoubleDCT_1D _dct;

    DoubleDST_1D _dst;

    DoubleDHT_1D _dht;

    /**
     * The elements of this matrix.
     */
    List<double> _elements;

    /**
     * Constructs a matrix with a copy of the given values. The values are
     * copied. So subsequent changes in <tt>values</tt> are not reflected in the
     * matrix, and vice-versa.
     *
     * @param values
     *            The values to be filled into the new matrix.
     */
    factory DenseDoubleMatrix1D.from(List<double> values) {
        final m = new DenseDoubleMatrix1D(values.length);
        m.assign(values);
        return m;
    }

    /**
     * Constructs a matrix with a given number of cells. All entries are
     * initially <tt>0</tt>.
     *
     * @param size
     *            the number of cells the matrix shall have.
     * @throws IllegalArgumentException
     *             if <tt>size<0</tt>.
     */
//    DenseDoubleMatrix1D(int size) {
//        _setUp(size);
//        this._elements = new List<double>(size);
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
     * @throws IllegalArgumentException
     *             if <tt>size<0</tt>.
     */
    DenseDoubleMatrix1D(int size, [List<double> elements=null, int zero=0, int stride=1, bool isView=false]) {
        if (elements == null) {
          elements = new List<double>(size);
        }
        _setUp(size, zero, stride);
        this._elements = elements;
        this._isNoView = !isView;
    }

    double aggregate(final DoubleDoubleFunction aggr,
            final DoubleFunction f) {
        if (_size == 0)
            return Double.NaN;
        double a = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
            nthreads = Math.min(nthreads, _size);
            List<Future> futures = new List<Future>(nthreads);
            int k = _size / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstIdx = _size - j * k;
                final int lastIdx = (j == (nthreads - 1)) ? 0 : firstIdx - k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int idx = _zero + (firstIdx - 1) * _stride;
                        double a = f.apply(_elements[idx]);
                        for (int i = firstIdx - 1; --i >= lastIdx;) {
                            a = aggr.apply(a, f.apply(_elements[idx -= _stride]));
                        }
                        return a;
                });
            }
            a = ConcurrencyUtils.waitForCompletion(futures, aggr);
        } else {
            int idx = _zero + (_size - 1) * _stride;
            a = f.apply(_elements[idx]);
            for (int i = _size - 1; --i >= 0;) {
                a = aggr.apply(a, f.apply(_elements[idx -= _stride]));
            }
        }
        return a;
    }

    double aggregateIndex(final DoubleDoubleFunction aggr,
            final DoubleFunction f, final IntArrayList indexList) {
        if (this.size() == 0)
            return Double.NaN;
        final int size = indexList.size();
        final List<int> indexElements = indexList.elements();
        double a = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
            nthreads = Math.min(nthreads, size);
            List<Future> futures = new List<Future>(nthreads);
            int k = size / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstIdx = j * k;
                final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int idx = _zero + indexElements[firstIdx] * _stride;
                        double a = f.apply(_elements[idx]);
                        double elem;
                        for (int i = firstIdx + 1; i < lastIdx; i++) {
                            idx = _zero + indexElements[i] * _stride;
                            elem = _elements[idx];
                            a = aggr.apply(a, f.apply(elem));
                        }
                        return a;
                });
            }
            a = ConcurrencyUtils.waitForCompletion(futures, aggr);
        } else {
            double elem;
            int idx = _zero + indexElements[0] * _stride;
            a = f.apply(_elements[idx]);
            for (int i = 1; i < size; i++) {
                idx = _zero + indexElements[i] * _stride;
                elem = _elements[idx];
                a = aggr.apply(a, f.apply(elem));
            }
        }
        return a;
    }

    double aggregateMatrix(final DoubleMatrix1D other, final DoubleDoubleFunction aggr,
            final DoubleDoubleFunction f) {
        if (!(other is DenseDoubleMatrix1D)) {
            return super.aggregate(other, aggr, f);
        }
        checkSize(other);
        if (_size == 0)
            return Double.NaN;
        final int zeroOther = other.index(0) as int;
        final int strideOther = other.stride();
        final List<double> elementsOther = other.elements() as List<double>;
        double a = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
            nthreads = Math.min(nthreads, _size);
            List<Future> futures = new List<Future>(nthreads);
            int k = _size / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstIdx = j * k;
                final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int idx = _zero + firstIdx * _stride;
                        int idxOther = zeroOther + firstIdx * strideOther;
                        double a = f.apply(_elements[idx], elementsOther[idxOther]);
                        for (int i = firstIdx + 1; i < lastIdx; i++) {
                            idx += _stride;
                            idxOther += strideOther;
                            a = aggr.apply(a, f.apply(_elements[idx], elementsOther[idxOther]));
                        }
                        return a;
                });
            }
            a = ConcurrencyUtils.waitForCompletion(futures, aggr);
        } else {
            a = f.apply(_elements[_zero], elementsOther[zeroOther]);
            int idx = _zero;
            int idxOther = zeroOther;
            for (int i = 1; i < _size; i++) {
                idx += _stride;
                idxOther += strideOther;
                a = aggr.apply(a, f.apply(_elements[idx], elementsOther[idxOther]));
            }
        }
        return a;
    }

    DoubleMatrix1D assignFunc(final DoubleFunction function) {
        final double multiplicator;
        if (function is DoubleMult) {
            // x[i] = mult*x[i]
            multiplicator = function.multiplicator;
            if (multiplicator == 1) {
                return this;
            }
        } else {
            multiplicator = 0;
        }
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
            nthreads = Math.min(nthreads, _size);
            List<Future> futures = new List<Future>(nthreads);
            int k = _size / nthreads;
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
                                _elements[idx] = function.apply(_elements[idx]);
                                idx += _stride;
                            }
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            int idx = _zero - _stride;
            // specialization for speed
            if (function is DoubleMult) {
                // x[i] = mult*x[i]
                for (int k = _size; --k >= 0;) {
                    _elements[idx += _stride] *= multiplicator;
                }
            } else {
                // the general case x[i] = f(x[i])
                for (int k = _size; --k >= 0;) {
                    _elements[idx += _stride] = function.apply(_elements[idx]);
                }
            }
        }
        return this;
    }

    DoubleMatrix1D assignProc(final DoubleProcedure cond,
            final DoubleFunction function) {
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
            nthreads = Math.min(nthreads, _size);
            List<Future> futures = new List<Future>(nthreads);
            int k = _size / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstIdx = j * k;
                final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int idx = _zero + firstIdx * _stride;
                        for (int i = firstIdx; i < lastIdx; i++) {
                            if (cond.apply(_elements[idx]) == true) {
                                _elements[idx] = function.apply(_elements[idx]);
                            }
                            idx += _stride;
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            int idx = _zero;
            for (int i = 0; i < _size; i++) {
                if (cond.apply(_elements[idx]) == true) {
                    _elements[idx] = function.apply(_elements[idx]);
                }
                idx += _stride;
            }
        }
        return this;
    }

    DoubleMatrix1D assignProcValue(final DoubleProcedure cond, final double value) {
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
            nthreads = Math.min(nthreads, _size);
            List<Future> futures = new List<Future>(nthreads);
            int k = _size / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstIdx = j * k;
                final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int idx = _zero + firstIdx * _stride;
                        for (int i = firstIdx; i < lastIdx; i++) {
                            if (cond.apply(_elements[idx]) == true) {
                                _elements[idx] = value;
                            }
                            idx += _stride;
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            int idx = _zero;
            for (int i = 0; i < _size; i++) {
                if (cond.apply(_elements[idx]) == true) {
                    _elements[idx] = value;
                }
                idx += _stride;
            }
        }
        return this;
    }

    DoubleMatrix1D assign(final double value) {
        final List<double> elems = this._elements;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
            nthreads = Math.min(nthreads, _size);
            List<Future> futures = new List<Future>(nthreads);
            int k = _size / nthreads;
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
        } else {
            int idx = _zero;
            for (int i = 0; i < _size; i++) {
                elems[idx] = value;
                idx += _stride;
            }
        }
        return this;
    }

    DoubleMatrix1D assignValues(final List<double> values) {
        if (values.length != _size)
            throw new IllegalArgumentException("Must have same number of cells: length=" + values.length + "size()="
                    + size());
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if (_isNoView) {
            System.arraycopy(values, 0, this._elements, 0, values.length);
        } else {
            if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
                nthreads = Math.min(nthreads, _size);
                List<Future> futures = new List<Future>(nthreads);
                int k = _size / nthreads;
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
            } else {
                int idx = _zero;
                for (int i = 0; i < _size; i++) {
                    _elements[idx] = values[i];
                    idx += _stride;
                }
            }
        }
        return this;
    }

    DoubleMatrix1D assignMatrix(DoubleMatrix1D source) {
        // overriden for performance only
        if (!(source is DenseDoubleMatrix1D)) {
            super.assign(source);
            return this;
        }
        DenseDoubleMatrix1D other = source as DenseDoubleMatrix1D;
        if (other == this)
            return this;
        checkSize(other);
        if (_isNoView && other._isNoView) {
            // quickest
            System.arraycopy(other._elements, 0, this._elements, 0, this._elements.length);
            return this;
        }
        if (_haveSharedCells(other)) {
            DoubleMatrix1D c = other.copy();
            if (!(c is DenseDoubleMatrix1D)) {
                // should not happen
                super.assign(source);
                return this;
            }
            other = c as DenseDoubleMatrix1D;
        }

        final List<double> elementsOther = other._elements;
        if (_elements == null || elementsOther == null)
            throw new InternalError();
        final int zeroOther = other.index(0) as int;
        final int strideOther = other._stride;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
            nthreads = Math.min(nthreads, _size);
            List<Future> futures = new List<Future>(nthreads);
            int k = _size / nthreads;
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
        } else {
            int idx = _zero;
            int idxOther = zeroOther;
            for (int k = 0; k < _size; k++) {
                _elements[idx] = elementsOther[idxOther];
                idx += _stride;
                idxOther += strideOther;
            }
        }
        return this;
    }

    DoubleMatrix1D assignMatrixFunc(final DoubleMatrix1D y, final DoubleDoubleFunction function) {
        // overriden for performance only
        if (!(y is DenseDoubleMatrix1D)) {
            super.assign(y, function);
            return this;
        }
        checkSize(y);
        final int zeroOther = y.index(0) as int;
        final int strideOther = y.stride();
        final List<double> elementsOther = y.elements() as List<double>;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
            nthreads = Math.min(nthreads, _size);
            List<Future> futures = new List<Future>(nthreads);
            int k = _size / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstIdx = j * k;
                final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int idx = _zero + firstIdx * _stride;
                        int idxOther = zeroOther + firstIdx * strideOther;
                        // specialized for speed
                        if (function == cern.jet.math.tdouble.DoubleFunctions.mult) {
                            // x[i] = x[i] * y[i]
                            for (int k = firstIdx; k < lastIdx; k++) {
                                _elements[idx] *= elementsOther[idxOther];
                                idx += _stride;
                                idxOther += strideOther;
                            }
                        } else if (function == cern.jet.math.tdouble.DoubleFunctions.div) {
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
                                _elements[idx] = function.apply(_elements[idx], elementsOther[idxOther]);
                                idx += _stride;
                                idxOther += strideOther;
                            }
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {
            // specialized for speed
            int idx = _zero;
            int idxOther = zeroOther;
            if (function == cern.jet.math.tdouble.DoubleFunctions.mult) {
                // x[i] = x[i] * y[i]
                for (int k = 0; k < _size; k++) {
                    _elements[idx] *= elementsOther[idxOther];
                    idx += _stride;
                    idxOther += strideOther;
                }
            } else if (function == cern.jet.math.tdouble.DoubleFunctions.div) {
                // x[i] = x[i] / y[i]
                for (int k = 0; k < _size; k++) {
                    _elements[idx] /= elementsOther[idxOther];
                    idx += _stride;
                    idxOther += strideOther;
                }
            } else if (function is DoublePlusMultSecond) {
                double multiplicator = function.multiplicator;
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
                    _elements[idx] = function.apply(_elements[idx], elementsOther[idxOther]);
                    idx += _stride;
                    idxOther += strideOther;
                }
            }
        }
        return this;
    }

    int cardinality() {
        int cardinality = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
            nthreads = Math.min(nthreads, _size);
            List<Future> futures = new List<Future>(nthreads);
            List<int> results = new List<int>(nthreads);
            int k = _size / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstIdx = j * k;
                final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int cardinality = 0;
                        int idx = _zero + firstIdx * _stride;
                        for (int i = firstIdx; i < lastIdx; i++) {
                            if (_elements[idx] != 0)
                                cardinality++;
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
        } else {
            int idx = _zero;
            for (int i = 0; i < _size; i++) {
                if (_elements[idx] != 0)
                    cardinality++;
                idx += _stride;
            }
        }
        return cardinality;
    }

    /**
     * Computes the discrete cosine transform (DCT-II) of this matrix.
     *
     * @param scale
     *            if true then scaling is performed
     */
    void dct(bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        if (_dct == null) {
            _dct = new DoubleDCT_1D(_size);
        }
        if (_isNoView) {
            _dct.forward(_elements, scale);
        } else {
            DoubleMatrix1D copy = this.copy();
            _dct.forward(copy.elements() as List<double>, scale);
            this.assign(copy.elements() as List<double>);
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the discrete Hartley transform (DHT) of this matrix.
     *
     */
    void dht() {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        if (_dht == null) {
            _dht = new DoubleDHT_1D(_size);
        }
        if (_isNoView) {
            _dht.forward(_elements);
        } else {
            DoubleMatrix1D copy = this.copy();
            _dht.forward(copy.elements() as List<double>);
            this.assign(copy.elements() as List<double>);
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the discrete sine transform (DST-II) of this matrix.
     *
     * @param scale
     *            if true then scaling is performed
     */
    void dst(bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        if (_dst == null) {
            _dst = new DoubleDST_1D(_size);
        }
        if (_isNoView) {
            _dst.forward(_elements, scale);
        } else {
            DoubleMatrix1D copy = this.copy();
            _dst.forward(copy.elements() as List<double>, scale);
            this.assign(copy.elements() as List<double>);
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    List<double> elements() {
        return _elements;
    }

    /**
     * Computes the discrete Fourier transform (DFT) of this matrix. The
     * physical layout of the output data is as follows:
     *
     * <pre>
     * this[2*k] = Re[k], 0&lt;=k&lt;size/2
     * this[2*k+1] = Im[k], 0&lt;k&lt;size/2
     * this[1] = Re[size/2]
     * </pre>
     *
     * This method computes only half of the elements of the real transform. The
     * other half satisfies the symmetry condition. If you want the full real
     * forward transform, use <code>getFft</code>. To get back the original
     * data, use <code>ifft</code>.
     *
     */
    void fft() {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        if (_fft == null) {
            _fft = new DoubleFFT_1D(_size);
        }
        if (_isNoView) {
            _fft.realForward(_elements);
        } else {
            DoubleMatrix1D copy = this.copy();
            _fft.realForward(copy.elements() as List<double>);
            this.assign(copy.elements() as List<double>);
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Returns new complex matrix which is the discrete Fourier transform (DFT)
     * of this matrix.
     *
     * @return the discrete Fourier transform (DFT) of this matrix.
     */
    DenseDComplexMatrix1D getFft() {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        final List<double> elems;
        if (_isNoView == true) {
            elems = _elements;
        } else {
            elems = this.copy().elements() as List<double>;
        }
        DenseDComplexMatrix1D c = new DenseDComplexMatrix1D(_size);
        final List<double> elementsC = (c).elements();
        System.arraycopy(elems, 0, elementsC, 0, _size);
        if (_fft == null) {
            _fft = new DoubleFFT_1D(_size);
        }
        _fft.realForwardFull(elementsC);
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
        return c;
    }

    /**
     * Returns new complex matrix which is the inverse of the discrete Fourier
     * (IDFT) transform of this matrix.
     *
     * @return the inverse of the discrete Fourier transform (IDFT) of this
     *         matrix.
     */
    DenseDComplexMatrix1D getIfft(bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        final List<double> elems;
        if (_isNoView == true) {
            elems = _elements;
        } else {
            elems = this.copy().elements() as List<double>;
        }
        DenseDComplexMatrix1D c = new DenseDComplexMatrix1D(_size);
        final List<double> elementsC = (c).elements();
        System.arraycopy(elems, 0, elementsC, 0, _size);
        if (_fft == null) {
            _fft = new DoubleFFT_1D(_size);
        }
        _fft.realInverseFull(elementsC, scale);
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
        return c;
    }

    void getNonZeros(final IntArrayList indexList, final DoubleArrayList valueList) {
        bool fillIndexList = indexList != null;
        bool fillValueList = valueList != null;
        if (fillIndexList)
            indexList.clear();
        if (fillValueList)
            valueList.clear();
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

    void getPositiveValues(final IntArrayList indexList, final DoubleArrayList valueList) {
        bool fillIndexList = indexList != null;
        bool fillValueList = valueList != null;
        if (fillIndexList)
            indexList.clear();
        if (fillValueList)
            valueList.clear();
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

    void getNegativeValues(final IntArrayList indexList, final DoubleArrayList valueList) {
        bool fillIndexList = indexList != null;
        bool fillValueList = valueList != null;
        if (fillIndexList)
            indexList.clear();
        if (fillValueList)
            valueList.clear();
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

    List<double> getMaxLocation() {
        int location = 0;
        double maxValue = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
            nthreads = Math.min(nthreads, _size);
            List<Future> futures = new List<Future>(nthreads);
            List<List<double>> results = new List<List<double>>(nthreads);//[2];
            int k = _size / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstIdx = j * k;
                final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int idx = _zero + firstIdx * _stride;
                        double maxValue = _elements[idx];
                        int location = (idx - _zero) / _stride;
                        for (int i = firstIdx + 1; i < lastIdx; i++) {
                            idx += _stride;
                            if (maxValue < _elements[idx]) {
                                maxValue = _elements[idx];
                                location = (idx - _zero) / _stride;
                            }
                        }
                        return [ maxValue, location ];
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
        } else {
            maxValue = _elements[_zero];
            location = 0;
            int idx = _zero;
            for (int i = 1; i < _size; i++) {
                idx += _stride;
                if (maxValue < _elements[idx]) {
                    maxValue = _elements[idx];
                    location = (idx - _zero) / _stride;
                }
            }
        }
        return [ maxValue, location ];
    }

    List<double> getMinLocation() {
        int location = 0;
        double minValue = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
            nthreads = Math.min(nthreads, _size);
            List<Future> futures = new List<Future>(nthreads);
            List<List<double>> results = new List<List<double>>(nthreads);//[2];
            int k = _size / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstIdx = j * k;
                final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int idx = _zero + firstIdx * _stride;
                        double minValue = _elements[idx];
                        int location = (idx - _zero) / _stride;
                        for (int i = firstIdx + 1; i < lastIdx; i++) {
                            idx += _stride;
                            if (minValue > _elements[idx]) {
                                minValue = _elements[idx];
                                location = (idx - _zero) / _stride;
                            }
                        }
                        return [ minValue, location ];
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
        } else {
            minValue = _elements[_zero];
            location = 0;
            int idx = _zero;
            for (int i = 1; i < _size; i++) {
                idx += _stride;
                if (minValue > _elements[idx]) {
                    minValue = _elements[idx];
                    location = (idx - _zero) / _stride;
                }
            }
        }
        return [ minValue, location ];
    }

    double getQuick(int index) {
        return _elements[_zero + index * _stride];
    }

    /**
     * Computes the inverse of the discrete cosine transform (DCT-III) of this
     * matrix.
     *
     * @param scale
     *            if true then scaling is performed
     */
    void idct(bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        if (_dct == null) {
            _dct = new DoubleDCT_1D(_size);
        }
        if (_isNoView) {
            _dct.inverse(_elements, scale);
        } else {
            DoubleMatrix1D copy = this.copy();
            _dct.inverse(copy.elements() as List<double>, scale);
            this.assign(copy.elements() as List<double>);
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the inverse of the discrete Hartley transform (IDHT) of this
     * matrix.
     *
     * @param scale
     *            if true then scaling is performed
     */
    void idht(bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        if (_dht == null) {
            _dht = new DoubleDHT_1D(_size);
        }
        if (_isNoView) {
            _dht.inverse(_elements, scale);
        } else {
            DoubleMatrix1D copy = this.copy();
            _dht.inverse(copy.elements() as List<double>, scale);
            this.assign(copy.elements() as List<double>);
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the inverse of discrete sine transform (DST-III) of this matrix.
     *
     * @param scale
     *            if true then scaling is performed
     */
    void idst(bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        if (_dst == null) {
            _dst = new DoubleDST_1D(_size);
        }
        if (_isNoView) {
            _dst.inverse(_elements, scale);
        } else {
            DoubleMatrix1D copy = this.copy();
            _dst.inverse(copy.elements() as List<double>, scale);
            this.assign(copy.elements() as List<double>);
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    /**
     * Computes the inverse of the discrete Fourier transform (DFT) of this
     * matrix. The physical layout of the input data has to be as follows:
     *
     * <pre>
     * this[2*k] = Re[k], 0&lt;=k&lt;size/2
     * this[2*k+1] = Im[k], 0&lt;k&lt;size/2
     * this[1] = Re[size/2]
     * </pre>
     *
     * This method computes only half of the elements of the real transform. The
     * other half satisfies the symmetry condition. If you want the full real
     * inverse transform, use <code>getIfft</code>.
     *
     */
    void ifft(bool scale) {
        int oldNthreads = ConcurrencyUtils.getNumberOfThreads();
        ConcurrencyUtils.setNumberOfThreads(ConcurrencyUtils.nextPow2(oldNthreads));
        if (_fft == null) {
            _fft = new DoubleFFT_1D(_size);
        }
        if (_isNoView) {
            _fft.realInverse(_elements, scale);
        } else {
            DoubleMatrix1D copy = this.copy();
            _fft.realInverse(copy.elements() as List<double>, scale);
            this.assign(copy.elements() as List<double>);
        }
        ConcurrencyUtils.setNumberOfThreads(oldNthreads);
    }

    DoubleMatrix1D like(int size) {
        return new DenseDoubleMatrix1D(size);
    }

    DoubleMatrix2D like2D(int rows, int columns) {
        return new DenseDoubleMatrix2D(rows, columns);
    }

    DoubleMatrix2D reshape(final int rows, final int columns) {
        if (rows * columns != _size) {
            throw new IllegalArgumentException("rows*columns != size");
        }
        DoubleMatrix2D M = new DenseDoubleMatrix2D(rows, columns);
        final List<double> elementsOther = M.elements() as List<double>;
        final int zeroOther = M.index(0, 0) as int;
        final int rowStrideOther = M.rowStride();
        final int columnStrideOther = M.columnStride();
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
            nthreads = Math.min(nthreads, _size);
            List<Future> futures = new List<Future>(nthreads);
            int k = columns / nthreads;
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
        } else {
            int idxOther;
            int idx = _zero;
            for (int c = 0; c < columns; c++) {
                idxOther = zeroOther + c * columnStrideOther;
                for (int r = 0; r < rows; r++) {
                    elementsOther[idxOther] = _elements[idx];
                    idxOther += rowStrideOther;
                    idx += _stride;
                }
            }
        }
        return M;
    }

    DoubleMatrix3D reshape3D(final int slices, final int rows, final int columns) {
        if (slices * rows * columns != _size) {
            throw new IllegalArgumentException("slices*rows*columns != size");
        }
        DoubleMatrix3D M = new DenseDoubleMatrix3D(slices, rows, columns);
        final List<double> elementsOther = M.elements() as List<double>;
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
    }

    void setQuick(int index, double value) {
        _elements[_zero + index * _stride] = value;
    }

    void swap(final DoubleMatrix1D other) {
        // overriden for performance only
        if (!(other is DenseDoubleMatrix1D)) {
            super.swap(other);
        }
        DenseDoubleMatrix1D y = other as DenseDoubleMatrix1D;
        if (y == this)
            return;
        checkSize(y);
        final List<double> elementsOther = y._elements;
        if (_elements == null || elementsOther == null)
            throw new InternalError();
        final int zeroOther = other.index(0) as int;
        final int strideOther = other.stride();
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
            nthreads = Math.min(nthreads, _size);
            List<Future> futures = new List<Future>(nthreads);
            int k = _size / nthreads;
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
        } else {
            int idx = _zero;
            int idxOther = zeroOther;
            for (int k = 0; k < _size; k++) {
                double tmp = _elements[idx];
                _elements[idx] = elementsOther[idxOther];
                elementsOther[idxOther] = tmp;
                idx += _stride;
                idxOther += strideOther;
            }
        }
    }

    void toArray(List<double> values) {
        if (values.length < _size)
            throw new IllegalArgumentException("values too small");
        if (this._isNoView)
            System.arraycopy(this._elements, 0, values, 0, this._elements.length);
        else
            super.toArray(values);
    }

    double zDotProduct(DoubleMatrix1D y, int from, int length) {
        if (!(y is DenseDoubleMatrix1D)) {
            return super.zDotProduct(y, from, length);
        }
        DenseDoubleMatrix1D yy = y as DenseDoubleMatrix1D;

        int tail = from + length;
        if (from < 0 || length < 0)
            return 0;
        if (_size < tail)
            tail = _size;
        if (y.size() < tail)
            tail = y.size() as int;
        final List<double> elementsOther = yy._elements;
        int zeroThis = index(from) as int;
        int zeroOther = yy.index(from) as int;
        int strideOther = yy._stride;
        if (_elements == null || elementsOther == null)
            throw new InternalError();
        double sum = 0;
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (length >= ConcurrencyUtils.getThreadsBeginN_1D())) {
            final int zeroThisF = zeroThis;
            final int zeroOtherF = zeroOther;
            final int strideOtherF = strideOther;
            nthreads = Math.min(nthreads, length);
            List<Future> futures = new List<Future>(nthreads);
            List<double> results = new List<double>(nthreads);
            int k = length / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstIdx = j * k;
                final int lastIdx = (j == nthreads - 1) ? length : firstIdx + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        int idx = zeroThisF + firstIdx * _stride;
                        int idxOther = zeroOtherF + firstIdx * strideOtherF;
                        idx -= _stride;
                        idxOther -= strideOtherF;
                        double sum = 0;
                        int min = lastIdx - firstIdx;
                        for (int k = min / 4; --k >= 0;) {
                            sum += _elements[idx += _stride] * elementsOther[idxOther += strideOtherF]
                                    + _elements[idx += _stride] * elementsOther[idxOther += strideOtherF]
                                    + _elements[idx += _stride] * elementsOther[idxOther += strideOtherF]
                                    + _elements[idx += _stride] * elementsOther[idxOther += strideOtherF];
                        }
                        for (int k = min % 4; --k >= 0;) {
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
        } else {
            zeroThis -= _stride;
            zeroOther -= strideOther;
            int min = tail - from;
            for (int k = min / 4; --k >= 0;) {
                sum += _elements[zeroThis += _stride] * elementsOther[zeroOther += strideOther]
                        + _elements[zeroThis += _stride] * elementsOther[zeroOther += strideOther]
                        + _elements[zeroThis += _stride] * elementsOther[zeroOther += strideOther]
                        + _elements[zeroThis += _stride] * elementsOther[zeroOther += strideOther];
            }
            for (int k = min % 4; --k >= 0;) {
                sum += _elements[zeroThis += _stride] * elementsOther[zeroOther += strideOther];
            }
        }
        return sum;
    }

    double zSum() {
        double sum = 0;
        final List<double> elems = this._elements;
        if (elems == null)
            throw new InternalError();
        int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (_size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
            nthreads = Math.min(nthreads, _size);
            List<Future> futures = new List<Future>(nthreads);
            List<double> results = new List<double>(nthreads);
            int k = _size / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstIdx = j * k;
                final int lastIdx = (j == nthreads - 1) ? _size : firstIdx + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        double sum = 0;
                        int idx = _zero + firstIdx * _stride;
                        for (int i = firstIdx; i < lastIdx; i++) {
                            sum += elems[idx];
                            idx += _stride;
                        }
                        return Double.valueOf(sum);
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
        } else {
            int idx = _zero;
            for (int k = 0; k < _size; k++) {
                sum += elems[idx];
                idx += _stride;
            }
        }
        return sum;
    }

    int _cardinality(int maxCardinality) {
        int cardinality = 0;
        int index = _zero;
        List<double> elems = this._elements;
        int i = _size;
        while (--i >= 0 && cardinality < maxCardinality) {
            if (elems[index] != 0)
                cardinality++;
            index += _stride;
        }
        return cardinality;
    }

    bool _haveSharedCellsRaw(DoubleMatrix1D other) {
        if (other is SelectedDenseDoubleMatrix1D) {
            SelectedDenseDoubleMatrix1D otherMatrix = other as SelectedDenseDoubleMatrix1D;
            return this._elements == otherMatrix._elements;
        } else if (other is DenseDoubleMatrix1D) {
            DenseDoubleMatrix1D otherMatrix = other as DenseDoubleMatrix1D;
            return this._elements == otherMatrix._elements;
        }
        return false;
    }

    int index(int rank) {
        return _zero + rank * _stride;
    }

    DoubleMatrix1D _viewSelectionLike(List<int> offsets) {
        return new SelectedDenseDoubleMatrix1D(this._elements, offsets);
    }
}
