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
 * Dense 2-d matrix holding <tt>complex</tt> elements. <b>Implementation:</b>
 * <p>
 * Internally holds one single contigous one-dimensional array, addressed in row
 * major. Complex data is represented by 2 double values in sequence, i.e.
 * elements[idx] constitute the real part and elements[idx+1] constitute the
 * imaginary part, where idx = index(0,0) + row * rowStride + column *
 * columnStride. Note that this implementation is not synchronized.
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
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
class DenseDComplexMatrix2D extends DComplexMatrix2D {

  /**
   * The elements of this matrix. elements are stored in row major. Complex
   * data is represented by 2 double values in sequence, i.e. elements[idx]
   * constitute the real part and elements[idx+1] constitute the imaginary
   * part, where idx = index(0,0) + row * rowStride + column * columnStride.
   */
  Float64List _elements;

  /**
   * Constructs a matrix with a copy of the given values. <tt>values</tt> is
   * required to have the form
   * <tt>re = values[row][2*column]; im = values[row][2*column+1]</tt> and
   * have exactly the same number of rows and columns as the receiver. Due to
   * the fact that complex data is represented by 2 double values in sequence:
   * the real and imaginary parts, the new matrix will be of the size
   * values.length by values[0].length / 2.
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
  factory DenseDComplexMatrix2D.fromValues(List<Float64List> values) {
    return new DenseDComplexMatrix2D(values.length, values.length == 0 ? 0 : values[0].length / 2)
      ..setAll2D(values);
  }

  /**
   * Constructs a complex matrix with the same size as <tt>realPart</tt>
   * matrix and fills the real part of this matrix with elements of
   * <tt>realPart</tt>.
   *
   * @param realPart
   *            a real matrix whose elements become a real part of this matrix
   * @throws ArgumentError
   *             if
   *             <tt>rows<0 || columns<0 || (double)columns*rows > Integer.MAX_VALUE</tt>
   *             .
   */
  factory DenseDComplexMatrix2D.fromRealPart(DoubleMatrix2D realPart) {
    return new DenseDComplexMatrix2D(realPart.rows, realPart.columns)..setReal(realPart);
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
  /*DenseDComplexMatrix2D(int rows, int columns) {
        _setUp(rows, columns, 0, 0, 2 * columns, 2);
        this._elements = new Float64List(rows * 2 * columns);
    }*/

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
   * @param isNoView
   *            if false then the view is constructed
   *
   * @throws ArgumentError
   *             if
   *             <tt>rows<0 || columns<0 || (double)columns*rows > Integer.MAX_VALUE</tt>
   *             or flip's are illegal.
   */
  DenseDComplexMatrix2D(int rows, int columns, [Float64List elements = null, int rowZero = 0, int columnZero = 0, int rowStride = null, int columnStride = 2, bool isNoView = true]) {
    if (rowStride == null) {
      rowStride = 2 * columns;
    }
    if (elements == null) {
      elements = new Float64List(rows * 2 * columns);
    }
    _setUp(rows, columns, rowZero, columnZero, rowStride, columnStride);
    this._elements = elements;
    this._isNoView = isNoView;
  }

  Float64List reduce(final cfunc.DComplexDComplexDComplexFunction aggr, final cfunc.DComplexDComplexFunction f) {
    Float64List b = new Float64List(2);
    if (length == 0) {
      b[0] = double.NAN;
      b[1] = double.NAN;
      return b;
    }
    final int zero = index(0, 0);
    Float64List a = null;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = zero + firstRow * _rowStride;
          Float64List a = f(_elements[idx], _elements[idx + 1]);
          int d = 1;
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = d; c < _columns; c++) {
              idx = zero + r * _rowStride + c * _columnStride;
              a = aggr(a, f(_elements[idx], _elements[idx + 1]));
            }
            d = 0;
          }
          return a;
        });
      }
      a = ConcurrencyUtils.waitForCompletion(futures, aggr);
    } else {*/
    a = f(new Float64List.fromList([_elements[zero], _elements[zero + 1]]));
    int d = 1; // first cell already done
    int idx;
    for (int r = 0; r < _rows; r++) {
      for (int c = d; c < _columns; c++) {
        idx = zero + r * _rowStride + c * _columnStride;
        a = aggr(a, f(new Float64List.fromList([_elements[idx], _elements[idx + 1]])));
      }
      d = 0;
    }
    //}
    return a;
  }

  Float64List reduceMatrix(final DComplexMatrix2D other, final cfunc.DComplexDComplexDComplexFunction aggr, final cfunc.DComplexDComplexDComplexFunction f) {
    if (!(other is DenseDComplexMatrix2D)) {
      return super.reduceMatrix(other, aggr, f);
    }
    checkShape(other);
    Float64List b = new Float64List(2);
    if (length == 0) {
      b[0] = double.NAN;
      b[1] = double.NAN;
      return b;
    }
    final int zero = index(0, 0);
    final int zeroOther = other.index(0, 0);
    final int rowStrideOther = other.rowStride;
    final int columnStrideOther = other.columnStride;
    final Float64List elemsOther = other.elements() as Float64List;
    Float64List a = null;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
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
          Float64List a = f([_elements[idx], _elements[idx + 1]], [elemsOther[idxOther], elemsOther[idxOther + 1]]);
          int d = 1;
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = d; c < _columns; c++) {
              idx = zero + r * _rowStride + c * _columnStride;
              idxOther = zeroOther + r * rowStrideOther + c * columnStrideOther;
              a = aggr(a, f([_elements[idx], _elements[idx + 1]], [elemsOther[idxOther], elemsOther[idxOther + 1]]));
            }
            d = 0;
          }
          return a;
        });
      }
      a = ConcurrencyUtils.waitForCompletion(futures, aggr);
    } else {*/
    int idx;
    int idxOther;
    a = f(new Float64List.fromList([_elements[zero], _elements[zero + 1]]), new Float64List.fromList([elemsOther[zeroOther], elemsOther[zeroOther + 1]]));
    int d = 1; // first cell already done
    for (int r = 0; r < _rows; r++) {
      for (int c = d; c < _columns; c++) {
        idx = zero + r * _rowStride + c * _columnStride;
        idxOther = zeroOther + r * rowStrideOther + c * columnStrideOther;
        a = aggr(a, f(new Float64List.fromList([_elements[idx], _elements[idx + 1]]), new Float64List.fromList([elemsOther[idxOther], elemsOther[idxOther + 1]])));
      }
      d = 0;
    }
    //}
    return a;
  }

  DComplexMatrix2D forEach(final cfunc.DComplexDComplexFunction function) {
    final int zero = index(0, 0);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      if (function is DComplexMult) {
        Float64List multiplicator = (function as DComplexMult).multiplicator;
        if (multiplicator[0] == 1 && multiplicator[1] == 0) return this;
        if (multiplicator[0] == 0 && multiplicator[1] == 0) return assignValue(0, 0);
      }
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = zero + firstRow * _rowStride;
          Float64List tmp = new Float64List(2);
          if (function is DComplexMult) {
            Float64List multiplicator = (function as DComplexMult).multiplicator;
            // x[i] = mult*x[i]
            for (int r = firstRow; r < lastRow; r++) {
              for (int i = idx,
                  c = 0; c < _columns; c++) {
                tmp[0] = _elements[i];
                tmp[1] = _elements[i + 1];
                _elements[i] = tmp[0] * multiplicator[0] - tmp[1] * multiplicator[1];
                _elements[i + 1] = tmp[1] * multiplicator[0] + tmp[0] * multiplicator[1];
                i += _columnStride;
              }
              idx += _rowStride;
            }
          } else {
            for (int r = firstRow; r < lastRow; r++) {
              for (int i = idx,
                  c = 0; c < _columns; c++) {
                tmp = function(_elements[i], _elements[i + 1]);
                _elements[i] = tmp[0];
                _elements[i + 1] = tmp[1];
                i += _columnStride;
              }
              idx += _rowStride;
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idx = zero;
    Float64List tmp = new Float64List(2);
    if (function is cfunc.DComplexMult) {
      Float64List multiplicator = (function as cfunc.DComplexMult).multiplicator;
      // x[i] = mult*x[i]
      for (int r = 0; r < _rows; r++) {
        for (int i = idx,
            c = 0; c < _columns; c++) {
          tmp[0] = _elements[i];
          tmp[1] = _elements[i + 1];
          _elements[i] = tmp[0] * multiplicator[0] - tmp[1] * multiplicator[1];
          _elements[i + 1] = tmp[1] * multiplicator[0] + tmp[0] * multiplicator[1];
          i += _columnStride;
        }
        idx += _rowStride;
      }
    } else {
      for (int r = 0; r < _rows; r++) {
        for (int i = idx,
            c = 0; c < _columns; c++) {
          tmp = function(new Float64List.fromList([_elements[i], _elements[i + 1]]));
          _elements[i] = tmp[0];
          _elements[i + 1] = tmp[1];
          i += _columnStride;
        }
        idx += _rowStride;
      }
    }
    //}
    return this;
  }

  DComplexMatrix2D forEachWhere(final cfunc.DComplexProcedure cond, final cfunc.DComplexDComplexFunction function) {
    final int zero = index(0, 0);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          Float64List elem = new Float64List(2);
          int idx = zero + firstRow * _rowStride;
          for (int r = firstRow; r < lastRow; r++) {
            for (int i = idx,
                c = 0; c < _columns; c++) {
              elem[0] = _elements[i];
              elem[1] = _elements[i + 1];
              if (cond(elem) == true) {
                elem = function(elem);
                _elements[i] = elem[0];
                _elements[i + 1] = elem[1];
              }
              i += _columnStride;
            }
            idx += _rowStride;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    Float64List elem = new Float64List(2);
    int idx = zero;
    for (int r = 0; r < _rows; r++) {
      for (int i = idx,
          c = 0; c < _columns; c++) {
        elem[0] = _elements[i];
        elem[1] = _elements[i + 1];
        if (cond(elem) == true) {
          elem = function(elem);
          _elements[i] = elem[0];
          _elements[i + 1] = elem[1];
        }
        i += _columnStride;
      }
      idx += _rowStride;
    }
    //}
    return this;
  }

  DComplexMatrix2D fillWhere(final cfunc.DComplexProcedure cond, final Float64List value) {
    final int zero = index(0, 0);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = zero + firstRow * _rowStride;
          Float64List elem = new Float64List(2);
          for (int r = firstRow; r < lastRow; r++) {
            for (int i = idx,
                c = 0; c < _columns; c++) {
              elem[0] = _elements[i];
              elem[1] = _elements[i + 1];
              if (cond(elem) == true) {
                _elements[i] = value[0];
                _elements[i + 1] = value[1];
              }
              i += _columnStride;
            }
            idx += _rowStride;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    Float64List elem = new Float64List(2);
    int idx = zero;
    for (int r = 0; r < _rows; r++) {
      for (int i = idx,
          c = 0; c < _columns; c++) {
        elem[0] = _elements[i];
        elem[1] = _elements[i + 1];
        if (cond(elem) == true) {
          _elements[i] = value[0];
          _elements[i + 1] = value[1];
        }
        i += _columnStride;
      }
      idx += _rowStride;
    }
    //}
    return this;
  }

  DComplexMatrix2D forEachReal(final cfunc.DComplexRealFunction function) {
    final int zero = index(0, 0);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = zero + firstRow * _rowStride;
          Float64List tmp = new Float64List(2);
          if (function == cfunc.abs) {
            for (int r = firstRow; r < lastRow; r++) {
              for (int i = idx,
                  c = 0; c < _columns; c++) {
                tmp[0] = _elements[i];
                tmp[1] = _elements[i + 1];
                double absX = Math.abs(_elements[i]);
                double absY = Math.abs(_elements[i + 1]);
                if (absX == 0 && absY == 0) {
                  _elements[i] = 0;
                } else if (absX >= absY) {
                  double d = tmp[1] / tmp[0];
                  _elements[i] = absX * Math.sqrt(1 + d * d);
                } else {
                  double d = tmp[0] / tmp[1];
                  _elements[i] = absY * Math.sqrt(1 + d * d);
                }
                _elements[i + 1] = 0;
                i += _columnStride;
              }
              idx += _rowStride;
            }
          } else {
            for (int r = firstRow; r < lastRow; r++) {
              for (int i = idx,
                  c = 0; c < _columns; c++) {
                tmp[0] = _elements[i];
                tmp[1] = _elements[i + 1];
                tmp[0] = function(tmp);
                _elements[i] = tmp[0];
                _elements[i + 1] = 0;
                i += _columnStride;
              }
              idx += _rowStride;
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idx = zero;
    Float64List tmp = new Float64List(2);
    if (function == cfunc.abs) {
      for (int r = 0; r < _rows; r++) {
        for (int i = idx,
            c = 0; c < _columns; c++) {
          tmp[0] = _elements[i];
          tmp[1] = _elements[i + 1];
          double absX = tmp[0].abs();
          double absY = tmp[1].abs();
          if (absX == 0 && absY == 0) {
            _elements[i] = 0.0;
          } else if (absX >= absY) {
            double d = tmp[1] / tmp[0];
            _elements[i] = absX * Math.sqrt(1 + d * d);
          } else {
            double d = tmp[0] / tmp[1];
            _elements[i] = absY * Math.sqrt(1 + d * d);
          }
          _elements[i + 1] = 0.0;
          i += _columnStride;
        }
        idx += _rowStride;
      }
    } else {
      for (int r = 0; r < _rows; r++) {
        for (int i = idx,
            c = 0; c < _columns; c++) {
          tmp[0] = _elements[i];
          tmp[1] = _elements[i + 1];
          tmp[0] = function(tmp);
          _elements[i] = tmp[0];
          _elements[i + 1] = 0.0;
          i += _columnStride;
        }
        idx += _rowStride;
      }
    }
    //}
    return this;
  }

  DComplexMatrix2D copyFrom(final DComplexMatrix2D source) {
    // overriden for performance only
    if (!(source is DenseDComplexMatrix2D)) {
      super.copyFrom(source);
      return this;
    }
    DenseDComplexMatrix2D other = source as DenseDComplexMatrix2D;
    if (other == this) return this; // nothing to do
    checkShape(other);
    if (this._isNoView && other._isNoView) { // quickest
      //System.arraycopy(other._elements, 0, this._elements, 0, this._elements.length);
      this._elements.setAll(0, other._elements);
      return this;
    }
    if (_haveSharedCells(other)) {
      DComplexMatrix2D c = other.copy();
      if (!(c is DenseDComplexMatrix2D)) { // should not happen
        super.copyFrom(other);
        return this;
      }
      other = c as DenseDComplexMatrix2D;
    }

    final Float64List elemsOther = other._elements;
    if (_elements == null || elemsOther == null) {
      throw new Error();
    }
    final int columnStrideOther = other._columnStride;
    final int rowStrideOther = other._rowStride;
    final int zeroOther = other.index(0, 0);
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
          int idxOther = zeroOther + firstRow * rowStrideOther;
          for (int r = firstRow; r < lastRow; r++) {
            for (int i = idx,
                j = idxOther,
                c = 0; c < _columns; c++) {
              _elements[i] = elemsOther[j];
              _elements[i + 1] = elemsOther[j + 1];
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
        _elements[i] = elemsOther[j];
        _elements[i + 1] = elemsOther[j + 1];
        i += _columnStride;
        j += columnStrideOther;
      }
      idx += _rowStride;
      idxOther += rowStrideOther;
    }
    //}
    return this;
  }

  DComplexMatrix2D forEachMatrix(final DComplexMatrix2D y, final cfunc.DComplexDComplexDComplexFunction function) {
    // overriden for performance only
    if (!(y is DenseDComplexMatrix2D)) {
      super.forEachMatrix(y, function);
      return this;
    }
    checkShape(y);
    final Float64List elemsOther = (y as DenseDComplexMatrix2D)._elements;
    if (_elements == null || elemsOther == null) {
      throw new Error();
    }
    final int columnStrideOther = y.columnStride;
    final int rowStrideOther = y.rowStride;
    final int zeroOther = y.index(0, 0);
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
          int idxOther = zeroOther + firstRow * rowStrideOther;
          Float64List tmp1 = new Float64List(2);
          Float64List tmp2 = new Float64List(2);
          if (function == cfunc.mult) {
            for (int r = firstRow; r < lastRow; r++) {
              for (int i = idx,
                  j = idxOther,
                  c = 0; c < _columns; c++) {
                tmp1[0] = _elements[i];
                tmp1[1] = _elements[i + 1];
                tmp2[0] = elemsOther[j];
                tmp2[1] = elemsOther[j + 1];
                _elements[i] = tmp1[0] * tmp2[0] - tmp1[1] * tmp2[1];
                _elements[i + 1] = tmp1[1] * tmp2[0] + tmp1[0] * tmp2[1];
                i += _columnStride;
                j += columnStrideOther;
              }
              idx += _rowStride;
              idxOther += rowStrideOther;
            }
          } else if (function == cfunc.multConjFirst) {
            for (int r = firstRow; r < lastRow; r++) {
              for (int i = idx,
                  j = idxOther,
                  c = 0; c < _columns; c++) {
                tmp1[0] = _elements[i];
                tmp1[1] = _elements[i + 1];
                tmp2[0] = elemsOther[j];
                tmp2[1] = elemsOther[j + 1];
                _elements[i] = tmp1[0] * tmp2[0] + tmp1[1] * tmp2[1];
                _elements[i + 1] = -tmp1[1] * tmp2[0] + tmp1[0] * tmp2[1];
                i += _columnStride;
                j += columnStrideOther;
              }
              idx += _rowStride;
              idxOther += rowStrideOther;
            }

          } else if (function == cfunc.multConjSecond) {
            for (int r = firstRow; r < lastRow; r++) {
              for (int i = idx,
                  j = idxOther,
                  c = 0; c < _columns; c++) {
                tmp1[0] = _elements[i];
                tmp1[1] = _elements[i + 1];
                tmp2[0] = elemsOther[j];
                tmp2[1] = elemsOther[j + 1];
                _elements[i] = tmp1[0] * tmp2[0] + tmp1[1] * tmp2[1];
                _elements[i + 1] = tmp1[1] * tmp2[0] - tmp1[0] * tmp2[1];
                i += _columnStride;
                j += columnStrideOther;
              }
              idx += _rowStride;
              idxOther += rowStrideOther;
            }
          } else {
            for (int r = firstRow; r < lastRow; r++) {
              for (int i = idx,
                  j = idxOther,
                  c = 0; c < _columns; c++) {
                tmp1[0] = _elements[i];
                tmp1[1] = _elements[i + 1];
                tmp2[0] = elemsOther[j];
                tmp2[1] = elemsOther[j + 1];
                tmp1 = function(tmp1, tmp2);
                _elements[i] = tmp1[0];
                _elements[i + 1] = tmp1[1];
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
    Float64List tmp1 = new Float64List(2);
    Float64List tmp2 = new Float64List(2);
    int idx = zero;
    int idxOther = zeroOther;
    if (function == cfunc.mult) {
      for (int r = 0; r < _rows; r++) {
        for (int i = idx,
            j = idxOther,
            c = 0; c < _columns; c++) {
          tmp1[0] = _elements[i];
          tmp1[1] = _elements[i + 1];
          tmp2[0] = elemsOther[j];
          tmp2[1] = elemsOther[j + 1];
          _elements[i] = tmp1[0] * tmp2[0] - tmp1[1] * tmp2[1];
          _elements[i + 1] = tmp1[1] * tmp2[0] + tmp1[0] * tmp2[1];
          i += _columnStride;
          j += columnStrideOther;
        }
        idx += _rowStride;
        idxOther += rowStrideOther;
      }
    } else if (function == cfunc.multConjFirst) {
      for (int r = 0; r < _rows; r++) {
        for (int i = idx,
            j = idxOther,
            c = 0; c < _columns; c++) {
          tmp1[0] = _elements[i];
          tmp1[1] = _elements[i + 1];
          tmp2[0] = elemsOther[j];
          tmp2[1] = elemsOther[j + 1];
          _elements[i] = tmp1[0] * tmp2[0] + tmp1[1] * tmp2[1];
          _elements[i + 1] = -tmp1[1] * tmp2[0] + tmp1[0] * tmp2[1];
          i += _columnStride;
          j += columnStrideOther;
        }
        idx += _rowStride;
        idxOther += rowStrideOther;
      }

    } else if (function == cfunc.multConjSecond) {
      for (int r = 0; r < _rows; r++) {
        for (int i = idx,
            j = idxOther,
            c = 0; c < _columns; c++) {
          tmp1[0] = _elements[i];
          tmp1[1] = _elements[i + 1];
          tmp2[0] = elemsOther[j];
          tmp2[1] = elemsOther[j + 1];
          _elements[i] = tmp1[0] * tmp2[0] + tmp1[1] * tmp2[1];
          _elements[i + 1] = tmp1[1] * tmp2[0] - tmp1[0] * tmp2[1];
          i += _columnStride;
          j += columnStrideOther;
        }
        idx += _rowStride;
        idxOther += rowStrideOther;
      }
    } else {
      for (int r = 0; r < _rows; r++) {
        for (int i = idx,
            j = idxOther,
            c = 0; c < _columns; c++) {
          tmp1[0] = _elements[i];
          tmp1[1] = _elements[i + 1];
          tmp2[0] = elemsOther[j];
          tmp2[1] = elemsOther[j + 1];
          tmp1 = function(tmp1, tmp2);
          _elements[i] = tmp1[0];
          _elements[i + 1] = tmp1[1];
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

  DComplexMatrix2D fill(final double re, final double im) {
    final int zero = index(0, 0);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
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
            for (int i = idx,
                c = 0; c < _columns; c++) {
              _elements[i] = re;
              _elements[i + 1] = im;
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
        _elements[i] = re;
        _elements[i + 1] = im;
        i += _columnStride;
      }
      idx += _rowStride;
    }
    //}
    return this;
  }

  DComplexMatrix2D setAll(final Float64List values) {
    if (values.length != _rows * 2 * _columns) {
      throw new ArgumentError("Must have same length: length=${values.length} rows()*2*columns()=${rows * 2 * columns}");
    }
    if (this._isNoView) {
      //System.arraycopy(values, 0, this._elements, 0, values.length);
      this._elements.setAll(0, values);
    } else {
      final int zero = index(0, 0);
      /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
      if ((nthreads > 1) && (_rows * _columns >= ConcurrencyUtils.getThreadsBeginN_2D())) {
        nthreads = Math.min(nthreads, _rows);
        List<Future> futures = new List<Future>(nthreads);
        int k = _rows / nthreads;
        for (int j = 0; j < nthreads; j++) {
          final int firstRow = j * k;
          final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
          futures[j] = ConcurrencyUtils.submit(() {
            int idxOther = firstRow * _columns * 2;
            int idx = zero + firstRow * _rowStride;
            for (int r = firstRow; r < lastRow; r++) {
              for (int i = idx,
                  c = 0; c < _columns; c++) {
                _elements[i] = values[idxOther++];
                _elements[i + 1] = values[idxOther++];
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
          _elements[i + 1] = values[idxOther++];
          i += _columnStride;
        }
        idx += _rowStride;
      }
      //}
    }
    return this;
  }

  DComplexMatrix2D setAll2D(final List<Float64List> values) {
    if (values.length != _rows) {
      throw new ArgumentError("Must have same number of rows: rows=${values.length} rows()=${rows}");
    }
    //int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if (this._isNoView) {
      /*if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
        nthreads = Math.min(nthreads, _rows);
        List<Future> futures = new List<Future>(nthreads);
        int k = _rows / nthreads;
        for (int j = 0; j < nthreads; j++) {
          final int firstRow = j * k;
          final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
          futures[j] = ConcurrencyUtils.submit(() {
            int idx = 2 * _columns;
            int i = firstRow * _rowStride;
            for (int r = firstRow; r < lastRow; r++) {
              Float64List currentRow = values[r];
              if (currentRow.length != idx) throw new ArgumentError("Must have same number of columns in every row: columns=" + currentRow.length + "2*columns()=" + idx);
              System.arraycopy(currentRow, 0, _elements, i, idx);
              i += idx;
            }
          });
        }
        ConcurrencyUtils.waitForCompletion(futures);
      } else {*/
      int idx = 2 * _columns;
      int i = 0;
      for (int r = 0; r < _rows; r++) {
        Float64List currentRow = values[r];
        if (currentRow.length != idx) {
          throw new ArgumentError("Must have same number of columns in every row: columns=${currentRow.length} 2*columns()=$idx");
        }
        //System.arraycopy(currentRow, 0, this._elements, i, idx);
        this._elements.setAll(i, currentRow);
        i += idx;
      }
      //}
    } else {
      final int zero = index(0, 0);
      /*if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
        nthreads = Math.min(nthreads, _rows);
        List<Future> futures = new List<Future>(nthreads);
        int k = _rows / nthreads;
        for (int j = 0; j < nthreads; j++) {
          final int firstRow = j * k;
          final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
          futures[j] = ConcurrencyUtils.submit(() {
            int idx = zero + firstRow * _rowStride;
            for (int r = firstRow; r < lastRow; r++) {
              Float64List currentRow = values[r];
              if (currentRow.length != 2 * _columns) throw new ArgumentError("Must have same number of columns in every row: columns=" + currentRow.length + "2*columns()=" + 2 * columns());
              for (int i = idx,
                  c = 0; c < _columns; c++) {
                _elements[i] = currentRow[2 * c];
                _elements[i + 1] = currentRow[2 * c + 1];
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
        if (currentRow.length != 2 * _columns) {
          throw new ArgumentError("Must have same number of columns in every row: columns=${currentRow.length} 2*columns()=${2 * columns}");
        }
        for (int i = idx,
            c = 0; c < _columns; c++) {
          _elements[i] = currentRow[2 * c];
          _elements[i + 1] = currentRow[2 * c + 1];
          i += _columnStride;
        }
        idx += _rowStride;
      }
      //}
    }
    return this;
  }

  DComplexMatrix2D setImaginary(final DoubleMatrix2D other) {
    checkShape(other);
    final int columnStrideOther = other.columnStride;
    final int rowStrideOther = other.rowStride;
    final int zeroOther = other.index(0, 0);
    final int zero = index(0, 0);
    final Float64List elemsOther = (other as DenseDoubleMatrix2D).elements();
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
          int idxOther = zeroOther + firstRow * rowStrideOther;
          for (int r = firstRow; r < lastRow; r++) {
            for (int i = idx,
                j = idxOther,
                c = 0; c < _columns; c++) {
              _elements[i + 1] = elemsOther[j];
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
        _elements[i + 1] = elemsOther[j];
        i += _columnStride;
        j += columnStrideOther;
      }
      idx += _rowStride;
      idxOther += rowStrideOther;
    }
    //}
    return this;
  }

  DComplexMatrix2D setReal(final DoubleMatrix2D other) {
    checkShape(other);
    final int columnStrideOther = other.columnStride;
    final int rowStrideOther = other.rowStride;
    final int zeroOther = other.index(0, 0);
    final int zero = index(0, 0);
    final Float64List elemsOther = (other as DenseDoubleMatrix2D).elements();
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
          int idxOther = zeroOther + firstRow * rowStrideOther;
          for (int r = firstRow; r < lastRow; r++) {
            for (int i = idx,
                j = idxOther,
                c = 0; c < _columns; c++) {
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
    } else {*/
    int idx = zero;
    int idxOther = zeroOther;
    for (int r = 0; r < _rows; r++) {
      for (int i = idx,
          j = idxOther,
          c = 0; c < _columns; c++) {
        _elements[i] = elemsOther[j];
        i += _columnStride;
        j += columnStrideOther;
      }
      idx += _rowStride;
      idxOther += rowStrideOther;
    }
    //}
    return this;
  }

  int get cardinality {
    int cardinality = 0;
    final int zero = index(0, 0);
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
              if ((_elements[i] != 0.0) || (_elements[i + 1] != 0.0)) cardinality++;
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
        if ((_elements[i] != 0.0) || (_elements[i + 1] != 0.0)) {
          cardinality++;
        }
        i += _columnStride;
      }
      idx += _rowStride;
    }
    //}
    return cardinality;
  }

  DComplexMatrix2D forEachNonZero(final cfunc.IntIntDComplexFunction function) {
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
          Float64List value = new Float64List(2);
          for (int r = firstRow; r < lastRow; r++) {
            for (int i = idx,
                c = 0; c < _columns; c++) {
              value[0] = _elements[i];
              value[1] = _elements[i + 1];
              if (value[0] != 0 || value[1] != 0) {
                Float64List v = function(r, c, value);
                _elements[i] = v[0];
                _elements[i + 1] = v[1];
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
    Float64List value = new Float64List(2);
    for (int r = 0; r < _rows; r++) {
      for (int i = idx,
          c = 0; c < _columns; c++) {
        value[0] = _elements[i];
        value[1] = _elements[i + 1];
        if (value[0] != 0 || value[1] != 0) {
          Float64List v = function(r, c, value);
          _elements[i] = v[0];
          _elements[i + 1] = v[1];
        }
        i += _columnStride;
      }
      idx += _rowStride;
    }
    //}
    return this;
  }

  DComplexMatrix2D conjugateTranspose() {
    DComplexMatrix2D transpose = this.dice().copy();
    final Float64List elemsOther = (transpose as DenseDComplexMatrix2D)._elements;
    final int zeroOther = transpose.index(0, 0);
    final int columnStrideOther = transpose.columnStride;
    final int rowStrideOther = transpose.rowStride;
    final int columnsOther = transpose.columns;
    final int rowsOther = transpose.rows;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, rowsOther);
      List<Future> futures = new List<Future>(nthreads);
      int k = rowsOther ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? rowsOther : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idxOther = zeroOther + firstRow * rowStrideOther;
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = 0; c < columnsOther; c++) {
              elemsOther[idxOther + 1] = -elemsOther[idxOther + 1];
              idxOther += columnStrideOther;
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idxOther = zeroOther;
    for (int r = 0; r < rowsOther; r++) {
      for (int c = 0; c < columnsOther; c++) {
        elemsOther[idxOther + 1] = -elemsOther[idxOther + 1];
        idxOther += columnStrideOther;
      }
    }
    //}
    return transpose;
  }

  Float64List elements() {
    return _elements;
  }

  DoubleMatrix2D imaginary() {
    final DenseDoubleMatrix2D Im = new DenseDoubleMatrix2D(_rows, _columns);
    final Float64List elemsOther = Im.elements();
    final int columnStrideOther = Im.columnStride;
    final int rowStrideOther = Im.rowStride;
    final int zeroOther = Im.index(0, 0);
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
          int idxOther = zeroOther + firstRow * rowStrideOther;
          for (int r = firstRow; r < lastRow; r++) {
            for (int i = idx,
                j = idxOther,
                c = 0; c < _columns; c++) {
              elemsOther[j] = _elements[i + 1];
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
        elemsOther[j] = _elements[i + 1];
        i += _columnStride;
        j += columnStrideOther;
      }
      idx += _rowStride;
      idxOther += rowStrideOther;
    }
    //}
    return Im;
  }

  void nonZeros(final List<int> rowList, final List<int> columnList, final List<Float64List> valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    int idx = index(0, 0);
    for (int r = 0; r < _rows; r++) {
      for (int i = idx,
          c = 0; c < _columns; c++) {
        Float64List value = new Float64List(2);
        value[0] = _elements[i];
        value[1] = _elements[i + 1];
        if (value[0] != 0 || value[1] != 0) {
          rowList.add(r);
          columnList.add(c);
          valueList.add(value);
        }
        i += _columnStride;
      }
      idx += _rowStride;
    }

  }

  Float64List get(int row, int column) {
    int idx = _rowZero + row * _rowStride + _columnZero + column * _columnStride;
    return new Float64List.fromList([_elements[idx], _elements[idx + 1]]);
  }

  DoubleMatrix2D real() {
    final DenseDoubleMatrix2D R = new DenseDoubleMatrix2D(_rows, _columns);
    final Float64List elemsOther = R.elements();
    final int columnStrideOther = R.columnStride;
    final int rowStrideOther = R.rowStride;
    final int zeroOther = R.index(0, 0);
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
          int idxOther = zeroOther + firstRow * rowStrideOther;
          for (int r = firstRow; r < lastRow; r++) {
            for (int i = idx,
                j = idxOther,
                c = 0; c < _columns; c++) {
              elemsOther[j] = _elements[i];
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
        elemsOther[j] = _elements[i];
        i += _columnStride;
        j += columnStrideOther;
      }
      idx += _rowStride;
      idxOther += rowStrideOther;
    }
    //}
    return R;
  }

  DComplexMatrix2D like2D(int rows, int columns) {
    return new DenseDComplexMatrix2D(rows, columns);
  }

  DComplexMatrix1D like1D(int size) {
    return new DenseDComplexMatrix1D(size);
  }

  void setParts(int row, int column, double re, double im) {
    int idx = _rowZero + row * _rowStride + _columnZero + column * _columnStride;
    _elements[idx] = re;
    _elements[idx + 1] = im;
  }

  void set(int row, int column, Float64List value) {
    int idx = _rowZero + row * _rowStride + _columnZero + column * _columnStride;
    _elements[idx] = value[0];
    _elements[idx + 1] = value[1];
  }

  List<Float64List> toList() {
    final List<Float64List> values = new List<Float64List>.generate(_rows,
        (_) => new Float64List(2 * _columns));
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
              values[r][2 * c] = _elements[i];
              values[r][2 * c + 1] = _elements[i + 1];
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
        values[r][2 * c] = _elements[i];
        values[r][2 * c + 1] = _elements[i + 1];
        i += _columnStride;
      }
      idx += _rowStride;
    }
    //}
    return values;
  }

  DComplexMatrix1D vectorize() {
    final DComplexMatrix1D v = new DenseDComplexMatrix1D(this.length);
    final int zero = index(0, 0);
    final int zeroOther = v.index(0);
    final int strideOther = v.stride();
    final Float64List elemsOther = v.elements() as Float64List;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_rows * _columns >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _columns);
      List<Future> futures = new List<Future>(nthreads);
      int k = _columns ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstColumn = j * k;
        final int lastColumn = (j == nthreads - 1) ? _columns : firstColumn + k;
        final int firstIdx = j * k * _rows;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = 0;
          int idxOther = zeroOther + firstIdx * strideOther;
          for (int c = firstColumn; c < lastColumn; c++) {
            idx = zero + c * _columnStride;
            for (int r = 0; r < _rows; r++) {
              elemsOther[idxOther] = _elements[idx];
              elemsOther[idxOther + 1] = _elements[idx + 1];
              idx += _rowStride;
              idxOther += strideOther;
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idx = 0;
    int idxOther = zeroOther;
    for (int c = 0; c < _columns; c++) {
      idx = zero + c * _columnStride;
      for (int r = 0; r < _rows; r++) {
        elemsOther[idxOther] = _elements[idx];
        elemsOther[idxOther + 1] = _elements[idx + 1];
        idx += _rowStride;
        idxOther += strideOther;
      }
    }
    //}
    return v;
  }

  DComplexMatrix1D mult(final DComplexMatrix1D y, DComplexMatrix1D z, [Float64List alpha = null, Float64List beta = null, bool transposeA = false]) {
    if (alpha == null) {
      alpha = new Float64List.fromList([1.0, 0.0]);
    }
    if (beta == null) {
      beta = (z == null ? new Float64List.fromList([1.0, 0.0]) : new Float64List.fromList([0.0, 0.0]));
    }
    if (transposeA) {
      return conjugateTranspose().mult(y, z, alpha, beta, false);
    }
    DComplexMatrix1D zz;
    if (z == null) {
      zz = new DenseDComplexMatrix1D(this._rows);
    } else {
      zz = z;
    }
    if (_columns != y.length || _rows > zz.length) {
      throw new ArgumentError("Incompatible args: " + toStringShort() + ", " + y.toStringShort() + ", " + zz.toStringShort());
    }
    final Float64List elemsY = y.elements() as Float64List;
    final Float64List elemsZ = zz.elements() as Float64List;
    if (_elements == null || elemsY == null || elemsZ == null) {
      throw new Error();
    }
    final int strideY = y.stride();
    final int strideZ = zz.stride();
    final int zero = index(0, 0);
    final int zeroY = y.index(0);
    final int zeroZ = zz.index(0);
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
          double reS;
          double imS;
          double reA;
          double imA;
          double reY;
          double imY;
          double reZ;
          double imZ;
          for (int r = firstRow; r < lastRow; r++) {
            reS = 0;
            imS = 0;
            int idx = idxZero;
            int idxY = zeroY;
            for (int c = 0; c < _columns; c++) {
              reA = _elements[idx];
              imA = _elements[idx + 1];
              reY = elemsY[idxY];
              imY = elemsY[idxY + 1];
              reS += reA * reY - imA * imY;
              imS += imA * reY + reA * imY;
              idx += _columnStride;
              idxY += strideY;
            }
            reZ = elemsZ[idxZeroZ];
            imZ = elemsZ[idxZeroZ + 1];
            elemsZ[idxZeroZ] = reS * alpha[0] - imS * alpha[1] + reZ * beta[0] - imZ * beta[1];
            elemsZ[idxZeroZ + 1] = imS * alpha[0] + reS * alpha[1] + imZ * beta[0] + reZ * beta[1];
            idxZero += _rowStride;
            idxZeroZ += strideZ;
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    int idxZero = zero;
    int idxZeroZ = zeroZ;
    double reS;
    double imS;
    double reA;
    double imA;
    double reY;
    double imY;
    double reZ;
    double imZ;

    for (int r = 0; r < _rows; r++) {
      reS = 0.0;
      imS = 0.0;
      int idx = idxZero;
      int idxY = zeroY;
      for (int c = 0; c < _columns; c++) {
        reA = _elements[idx];
        imA = _elements[idx + 1];
        reY = elemsY[idxY];
        imY = elemsY[idxY + 1];
        reS += reA * reY - imA * imY;
        imS += imA * reY + reA * imY;
        idx += _columnStride;
        idxY += strideY;
      }
      reZ = elemsZ[idxZeroZ];
      imZ = elemsZ[idxZeroZ + 1];
      elemsZ[idxZeroZ] = reS * alpha[0] - imS * alpha[1] + reZ * beta[0] - imZ * beta[1];
      elemsZ[idxZeroZ + 1] = imS * alpha[0] + reS * alpha[1] + imZ * beta[0] + reZ * beta[1];
      idxZero += _rowStride;
      idxZeroZ += strideZ;
    }
    //}
    return zz;
  }

  DComplexMatrix2D multiply(final DComplexMatrix2D B, DComplexMatrix2D C, [Float64List alpha = null, Float64List beta = null, final bool transposeA = false, final bool transposeB = false]) {
    if (alpha == null) {
      alpha = new Float64List.fromList([1.0, 0.0]);
    }
    if (beta == null) {
      beta = (C == null ? new Float64List.fromList([1.0, 0.0]) : new Float64List.fromList([0.0, 0.0]));
    }
    final int rowsA = _rows;
    final int columnsA = _columns;
    final int rowsB = B.rows;
    final int columnsB = B.columns;
    final int rowsC = transposeA ? columnsA : rowsA;
    final int columnsC = transposeB ? rowsB : columnsB;

    if (C == null) {
      C = new DenseDComplexMatrix2D(rowsC, columnsC);
    }

    if (transposeA) {
      return conjugateTranspose().multiply(B, C, alpha, beta, false, transposeB);
    }
    if (transposeB) {
      return this.multiply(B.conjugateTranspose(), C, alpha, beta, transposeA, false);
    }
    if (B.rows != columnsA) {
      throw new ArgumentError("Matrix2D inner dimensions must agree:" + toStringShort() + ", " + B.toStringShort());
    }
    if (C.rows != rowsA || C.columns != columnsB) {
      throw new ArgumentError("Incompatibe result matrix: " + toStringShort() + ", " + B.toStringShort() + ", " + C.toStringShort());
    }
    if (this == C || B == C) {
      throw new ArgumentError("Matrices must not be identical");
    }
    int flops = 2 * rowsA * columnsA * columnsB;
    int noOfTasks = 1;//Math.min(flops / 30000, ConcurrencyUtils.getNumberOfThreads()); // each
    /* thread should process at least 30000 flops */
    bool splitB = (columnsB >= noOfTasks);
    int width = splitB ? columnsB : rowsA;
    noOfTasks = Math.min(width, noOfTasks);

    if (noOfTasks < 2) {
      return this._multiplySeq(B, C, alpha, beta, transposeA, transposeB);
    }
    // set up concurrent tasks
    int span = width ~/ noOfTasks;
    final List<Future> subTasks = new List<Future>(noOfTasks);
    for (int i = 0; i < noOfTasks; i++) {
      final int offset = i * span;
      if (i == noOfTasks - 1) span = width - span * i; // last span may be a bit larger
      DComplexMatrix2D AA, BB, CC;
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
        (AA as DenseDComplexMatrix2D)._zMultSeq(BB, CC, alpha, beta, transposeA, transposeB);
      });*/
      (AA as DenseDComplexMatrix2D)._multiplySeq(BB, CC, alpha, beta, transposeA, transposeB);
    }
    //ConcurrencyUtils.waitForCompletion(subTasks);

    return C;
  }

  DComplexMatrix2D _multiplySeq(DComplexMatrix2D B, DComplexMatrix2D C, Float64List alpha, Float64List beta, bool transposeA, bool transposeB) {
    if (transposeA) {
      return conjugateTranspose().multiply(B, C, alpha, beta, false, transposeB);
    }
    if (transposeB) {
      return this.multiply(B.conjugateTranspose(), C, alpha, beta, transposeA, false);
    }
    int m = _rows;
    int n = _columns;
    int p = B.columns;
    if (C == null) {
      C = new DenseDComplexMatrix2D(m, p);
    }
    if (!(C is DenseDComplexMatrix2D)) {
      return super.multiply(B, C, alpha, beta, transposeA, transposeB);
    }
    if (B.rows != n) {
      throw new ArgumentError("Matrix2D inner dimensions must agree:" + toStringShort() + ", " + B.toStringShort());
    }
    if (C.rows != m || C.columns != p) {
      throw new ArgumentError("Incompatibel result matrix: " + toStringShort() + ", " + B.toStringShort() + ", " + C.toStringShort());
    }
    if (this == C || B == C) {
      throw new ArgumentError("Matrices must not be identical");
    }

    DenseDComplexMatrix2D BB = B as DenseDComplexMatrix2D;
    DenseDComplexMatrix2D CC = C as DenseDComplexMatrix2D;
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
    int m_optimal = (BLOCK_SIZE - n) ~/ (n + 1);
    if (m_optimal <= 0) m_optimal = 1;
    int blocks = m ~/ m_optimal;
    int rr = 0;
    if (m % m_optimal != 0) blocks++;
    double reS;
    double imS;
    double reA;
    double imA;
    double reB;
    double imB;
    double reC;
    double imC;
    for ( ; --blocks >= 0; ) {
      int jB = BB.index(0, 0);
      int indexA = index(rr, 0);
      int jC = CC.index(rr, 0);
      rr += m_optimal;
      if (blocks == 0) m_optimal += m - rr;

      for (int j = p; --j >= 0; ) {
        int iA = indexA;
        int iC = jC;
        for (int i = m_optimal; --i >= 0; ) {
          int kA = iA;
          int kB = jB;
          reS = 0.0;
          imS = 0.0;
          // loop unrolled
          kA -= cA;
          kB -= rB;
          for (int k = n % 4; --k >= 0; ) {
            kA += cA;
            kB += rB;
            reA = AElems[kA];
            imA = AElems[kA + 1];
            reB = BElems[kB];
            imB = BElems[kB + 1];
            reS += reA * reB - imA * imB;
            imS += imA * reB + reA * imB;
          }
          for (int k = n ~/ 4; --k >= 0; ) {
            kA += cA;
            kB += rB;
            reA = AElems[kA];
            imA = AElems[kA + 1];
            reB = BElems[kB];
            imB = BElems[kB + 1];
            reS += reA * reB - imA * imB;
            imS += imA * reB + reA * imB;
            kA += cA;
            kB += rB;
            reA = AElems[kA];
            imA = AElems[kA + 1];
            reB = BElems[kB];
            imB = BElems[kB + 1];
            reS += reA * reB - imA * imB;
            imS += imA * reB + reA * imB;
            kA += cA;
            kB += rB;
            reA = AElems[kA];
            imA = AElems[kA + 1];
            reB = BElems[kB];
            imB = BElems[kB + 1];
            reS += reA * reB - imA * imB;
            imS += imA * reB + reA * imB;
            kA += cA;
            kB += rB;
            reA = AElems[kA];
            imA = AElems[kA + 1];
            reB = BElems[kB];
            imB = BElems[kB + 1];
            reS += reA * reB - imA * imB;
            imS += imA * reB + reA * imB;
          }
          reC = CElems[iC];
          imC = CElems[iC + 1];
          CElems[iC] = alpha[0] * reS - alpha[1] * imS + beta[0] * reC - beta[1] * imC;
          CElems[iC + 1] = alpha[1] * reS + alpha[0] * imS + beta[1] * reC + beta[0] * imC;
          iA += rA;
          iC += rC;
        }
        jB += cB;
        jC += cC;
      }
    }
    return C;
  }

  Float64List sum() {
    Float64List sum = new Float64List(2);
    final int zero = this.index(0, 0);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (_rows * _columns >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          Float64List sum = new Float64List(2);
          int idx = zero + firstRow * _rowStride;
          for (int r = firstRow; r < lastRow; r++) {
            for (int i = idx,
                c = 0; c < _columns; c++) {
              sum[0] += _elements[i];
              sum[1] += _elements[i + 1];
              i += _columnStride;
            }
            idx += _rowStride;
          }
          return sum;
        });
      }
      try {
        Float64List tmp;
        for (int j = 0; j < nthreads; j++) {
          tmp = futures[j].get() as Float64List;
          sum[0] = sum[0] + tmp[0];
          sum[1] = sum[1] + tmp[1];
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
        sum[0] += _elements[i];
        sum[1] += _elements[i + 1];
        i += _columnStride;
      }
      idx += _rowStride;
    }
    //}
    return sum;
  }

  bool _haveSharedCellsRaw(DComplexMatrix2D other) {
    if (other is SelectedDenseDComplexMatrix2D) {
      return this._elements == other._elements;
    } else if (other is DenseDComplexMatrix2D) {
      return this._elements == other._elements;
    }
    return false;
  }

  int index(int row, int column) {
    return _rowZero + row * _rowStride + _columnZero + column * _columnStride;
  }

  DComplexMatrix1D _like1D(int size, int zero, int stride) {
    return new DenseDComplexMatrix1D(size, this._elements, zero, stride, false);
  }

  DComplexMatrix2D _viewSelectionLike(Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedDenseDComplexMatrix2D.offsets(this._elements, rowOffsets, columnOffsets, 0);
  }

  Object clone() {
    return new DenseDComplexMatrix2D(_rows, _columns, _elements, _rowZero, _columnZero, _rowStride, _columnStride, _isNoView);
  }
}

/**
 * Selection view on dense 2-d matrices holding <tt>complex</tt> elements.
 * <b>Implementation:</b>
 * <p>
 * Objects of this class are typically constructed via <tt>viewIndexes</tt>
 * methods on some source matrix. The interface introduced in abstract super
 * classes defines everything a user can do. From a user point of view there is
 * nothing special about this class; it presents the same functionality with the
 * same signatures and semantics as its abstract superclass(es) while
 * introducing no additional functionality. Thus, this class need not be visible
 * to users.
 * <p>
 * This class uses no delegation. Its instances point directly to the data. Cell
 * addressing overhead is 1 additional int addition and 2 additional array index
 * accesses per get/set.
 * <p>
 * Note that this implementation is not synchronized.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
class SelectedDenseDComplexMatrix2D extends DComplexMatrix2D {

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
  factory SelectedDenseDComplexMatrix2D.offsets(Float64List elements, Int32List rowOffsets, Int32List columnOffsets, int offset) {
    return new SelectedDenseDComplexMatrix2D(rowOffsets.length, columnOffsets.length, elements, 0, 0, 1, 1, rowOffsets, columnOffsets, offset);
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
  SelectedDenseDComplexMatrix2D(int rows, int columns, Float64List elements, int rowZero, int columnZero, int rowStride, int columnStride, Int32List rowOffsets, Int32List columnOffsets, int offset) {
    // be sure parameters are valid, we do not check...
    _setUp(rows, columns, rowZero, columnZero, rowStride, columnStride);

    this._elements = elements;
    this._rowOffsets = rowOffsets;
    this._columnOffsets = columnOffsets;
    this._offset = offset;

    this._isNoView = false;
  }

  int _columnOffset(int absRank) {
    return _columnOffsets[absRank];
  }

  int _rowOffset(int absRank) {
    return _rowOffsets[absRank];
  }

  Float64List get(int row, int column) {
    int idxr = _rowZero + row * _rowStride;
    int idxc = _columnZero + column * _columnStride;
    return new Float64List.fromList([_elements[_offset + _rowOffsets[idxr] + _columnOffsets[idxc]], _elements[_offset + _rowOffsets[idxr] + _columnOffsets[idxc] + 1]]);
  }

  Float64List elements() {
    throw new UnsupportedError("This method is not supported.");

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
  bool _haveSharedCellsRaw(DComplexMatrix2D other) {
    if (other is SelectedDenseDComplexMatrix2D) {
      return this._elements == other._elements;
    } else if (other is DenseDComplexMatrix2D) {
      return this._elements == other._elements;
    }
    return false;
  }

  int index(int row, int column) {
    return this._offset + _rowOffsets[_rowZero + row * _rowStride] + _columnOffsets[_columnZero + column * _columnStride];
  }

  DComplexMatrix2D like2D(int rows, int columns) {
    return new DenseDComplexMatrix2D(rows, columns);
  }

  DComplexMatrix1D like1D(int size) {
    return new DenseDComplexMatrix1D(size);
  }

  DComplexMatrix1D _like1D(int size, int zero, int stride) {
    throw new Error(); // this method is never called since
    // viewRow() and viewColumn are overridden
    // properly.
  }

  void set(int row, int column, Float64List value) {
    int idxr = _rowZero + row * _rowStride;
    int idxc = _columnZero + column * _columnStride;
    _elements[_offset + _rowOffsets[idxr] + _columnOffsets[idxc]] = value[0];
    _elements[_offset + _rowOffsets[idxr] + _columnOffsets[idxc] + 1] = value[1];
  }

  DComplexMatrix1D vectorize() {
    throw new UnsupportedError("This method is not supported.");
  }

  void setParts(int row, int column, double re, double im) {
    int idxr = _rowZero + row * _rowStride;
    int idxc = _columnZero + column * _columnStride;
    _elements[_offset + _rowOffsets[idxr] + _columnOffsets[idxc]] = re;
    _elements[_offset + _rowOffsets[idxr] + _columnOffsets[idxc] + 1] = im;
  }

//  void _setUp(int rows, int columns, [int rowZero = 0, int columnZero = 0, int rowStride = null, int columnStride = 1]) {
//    super._setUp(rows, columns);
//    this._rowStride = 1;
//    this._columnStride = 2;
//    this._offset = 0;
//  }

  AbstractMatrix2D _vDice() {
    super._vDice();
    // swap
    Int32List tmp = _rowOffsets;
    _rowOffsets = _columnOffsets;
    _columnOffsets = tmp;

    this._isNoView = false;
    return this;
  }

  DComplexMatrix1D column(int column) {
    _checkColumn(column);
    int viewSize = this._rows;
    int viewZero = this._rowZero;
    int viewStride = this._rowStride;
    Int32List viewOffsets = this._rowOffsets;
    int viewOffset = this._offset + _columnOffset(_columnRank(column));
    return new SelectedDenseDComplexMatrix1D(viewSize, this._elements, viewZero, viewStride, viewOffsets, viewOffset);
  }

  DComplexMatrix1D row(int row) {
    _checkRow(row);
    int viewSize = this._columns;
    int viewZero = _columnZero;
    int viewStride = this._columnStride;
    Int32List viewOffsets = this._columnOffsets;
    int viewOffset = this._offset + _rowOffset(_rowRank(row));
    return new SelectedDenseDComplexMatrix1D(viewSize, this._elements, viewZero, viewStride, viewOffsets, viewOffset);
  }

  DComplexMatrix2D _viewSelectionLike(Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedDenseDComplexMatrix2D.offsets(this._elements, rowOffsets, columnOffsets, this._offset);
  }

  DoubleMatrix2D real() {
    final DenseDoubleMatrix2D R = new DenseDoubleMatrix2D(_rows, _columns);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
        if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
            nthreads = Math.min(nthreads, _rows);
            List<Future> futures = new List<Future>(nthreads);
            int k = _rows / nthreads;
            for (int j = 0; j < nthreads; j++) {
                final int firstRow = j * k;
                final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
                futures[j] = ConcurrencyUtils.submit(() {
                        Float64List tmp;
                        for (int r = firstRow; r < lastRow; r++) {
                            for (int c = 0; c < _columns; c++) {
                                tmp = getQuick(r, c);
                                R.setQuick(r, c, tmp[0]);
                            }
                        }
                });
            }
            ConcurrencyUtils.waitForCompletion(futures);
        } else {*/
    for (int r = 0; r < _rows; r++) {
      for (int c = 0; c < _columns; c++) {
        final tmp = get(r, c);
        R.set(r, c, tmp[0]);
      }
    }
    //}
    return R;
  }

  DoubleMatrix2D imaginary() {
    final DenseDoubleMatrix2D Im = new DenseDoubleMatrix2D(_rows, _columns);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
        nthreads = Math.min(nthreads, _rows);
        List<Future> futures = new List<Future>(nthreads);
        int k = _rows / nthreads;
        for (int j = 0; j < nthreads; j++) {
            final int firstRow = j * k;
            final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
            futures[j] = ConcurrencyUtils.submit(() {
                    Float64List tmp;
                    for (int r = firstRow; r < lastRow; r++) {
                        for (int c = 0; c < _columns; c++) {
                            tmp = getQuick(r, c);
                            Im.setQuick(r, c, tmp[1]);
                        }
                    }
            });
        }
        ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    for (int r = 0; r < _rows; r++) {
      for (int c = 0; c < _columns; c++) {
        final tmp = get(r, c);
        Im.set(r, c, tmp[1]);
      }
    }
    //}
    return Im;
  }

  Object clone() {
    return new SelectedDenseDComplexMatrix2D(_rows, _columns, _elements, _rowZero, _columnZero, _rowStride, _columnStride, _rowOffsets, _columnOffsets, _offset);
  }
}
