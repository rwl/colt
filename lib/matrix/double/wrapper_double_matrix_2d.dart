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
 * 2-d matrix holding <tt>double</tt> elements; either a view wrapping another
 * matrix or a matrix whose views are wrappers.
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 04/14/2000
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
class WrapperDoubleMatrix2D extends DoubleMatrix2D {

  /*
     * The elements of the matrix.
     */
  DoubleMatrix2D _content;

  WrapperDoubleMatrix2D(DoubleMatrix2D newContent) {
    if (newContent != null) {
      try {
        _setUp(newContent.rows(), newContent.columns());
      } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
        if ("matrix too large" != exc.message) throw exc;
      }
    }
    this._content = newContent;
  }

  DoubleMatrix2D assignValues(final Float64List values) {
    if (_content is DiagonalDoubleMatrix2D) {
      int dlength = (_content as DiagonalDoubleMatrix2D)._dlength;
      final Float64List elems = (_content as DiagonalDoubleMatrix2D)._elements;
      if (values.length != dlength) {
        throw new ArgumentError("Must have same length: length=${values.length} dlength=$dlength");
      }
      /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
      if ((nthreads > 1) && (dlength >= ConcurrencyUtils.getThreadsBeginN_2D())) {
        nthreads = Math.min(nthreads, dlength);
        List<Future> futures = new List<Future>(nthreads);
        int k = dlength ~/ nthreads;
        for (int j = 0; j < nthreads; j++) {
          final int firstIdx = j * k;
          final int lastIdx = (j == nthreads - 1) ? dlength : firstIdx + k;
          futures[j] = ConcurrencyUtils.submit(() {
            for (int i = firstIdx; i < lastIdx; i++) {
              elems[i] = values[i];
            }
          });
        }
        ConcurrencyUtils.waitForCompletion(futures);
      } else {*/
        for (int i = 0; i < dlength; i++) {
          elems[i] = values[i];
        }
      //}
      return this;
    } else {
      return super.assignValues(values);
    }
  }

  DoubleMatrix2D assignFunc(final DoubleMatrix2D y, final func.DoubleDoubleFunction function) {
    checkShape(y);
    if (y is WrapperDoubleMatrix2D) {
      final rowList = new List<int>();
      final columnList = new List<int>();
      final valueList = new List<double>();
      y.getNonZeros(rowList, columnList, valueList);
      assignFuncIndex(y, function,
          new Int32List.fromList(rowList),
          new Int32List.fromList(columnList));
    } else {
      super.assignFunc(y, function);
    }
    return this;
  }

  Object elements() {
    return _content.elements();
  }

  double getQuick(int row, int column) {
    return _content.getQuick(row, column);
  }

  bool equalsValue(double value) {
    if (_content is DiagonalDoubleMatrix2D) {
      double epsilon = DoubleProperty.DEFAULT.tolerance();
      Float64List elements = _content.elements() as Float64List;
      for (int r = 0; r < elements.length; r++) {
        double x = elements[r];
        double diff = (value - x).abs();
        if ((diff != diff) && ((value != value && x != x) || value == x)) {
          diff = 0.0;
        }
        if (!(diff <= epsilon)) {
          return false;
        }
      }
      return true;
    } else {
      return super.equalsValue(value);
    }
  }

  bool equals(Object obj) {
    if (_content is DiagonalDoubleMatrix2D && obj is DiagonalDoubleMatrix2D) {
      double epsilon = DoubleProperty.DEFAULT.tolerance();
      if (this == obj) return true;
      if (!(this != null && obj != null)) return false;
      DiagonalDoubleMatrix2D A = _content as DiagonalDoubleMatrix2D;
      DiagonalDoubleMatrix2D B = obj;
      if (A.columns() != B.columns() || A.rows() != B.rows() || A.diagonalIndex() != B.diagonalIndex() || A.diagonalLength() != B.diagonalLength()) return false;
      Float64List AElements = A.elements();
      Float64List BElements = B.elements();
      for (int r = 0; r < AElements.length; r++) {
        double x = AElements[r];
        double value = BElements[r];
        double diff = (value - x).abs();
        if ((diff != diff) && ((value != value && x != x) || value == x)) {
          diff = 0.0;
        }
        if (!(diff <= epsilon)) {
          return false;
        }
      }
      return true;
    } else {
      return super.equals(obj);
    }
  }

  DoubleMatrix2D like2D(int rows, int columns) {
    return _content.like2D(rows, columns);
  }

  DoubleMatrix1D like1D(int size) {
    return _content.like1D(size);
  }

  void setQuick(int row, int column, double value) {
    _content.setQuick(row, column, value);
  }

  DoubleMatrix1D vectorize() {
    final DenseDoubleMatrix1D v = new DenseDoubleMatrix1D(size());
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _columns);
      List<Future> futures = new List<Future>(nthreads);
      int k = _columns ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstCol = j * k;
        final int lastCol = (j == nthreads - 1) ? _columns : firstCol + k;
        final int firstidx = j * k * _rows;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = firstidx;
          for (int c = firstCol; c < lastCol; c++) {
            for (int r = 0; r < _rows; r++) {
              v.setQuick(idx++, getQuick(r, c));
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
      int idx = 0;
      for (int c = 0; c < _columns; c++) {
        for (int r = 0; r < _rows; r++) {
          v.setQuick(idx++, getQuick(r, c));
        }
      }
    //}
    return v;
  }

  DoubleMatrix1D viewColumn(int column) {
    return viewDice().viewRow(column);
  }

  DoubleMatrix2D viewColumnFlip() {
    if (_columns == 0) return this;
    WrapperDoubleMatrix2D view = new ViewColumnFlipWrapperDoubleMatrix2D(this);
    view._isNoView = false;

    return view;
  }

  DoubleMatrix2D viewDice() {
    WrapperDoubleMatrix2D view = new ViewDiceWrapperDoubleMatrix2D(this);
    view._rows = _columns;
    view._columns = _rows;
    view._isNoView = false;

    return view;
  }

  DoubleMatrix2D viewPart(final int row, final int column, int height, int width) {
    _checkBox(row, column, height, width);
    WrapperDoubleMatrix2D view = new ViewPartWrapperDoubleMatrix2D(this, row, column);
    view._rows = height;
    view._columns = width;
    view._isNoView = false;

    return view;
  }

  DoubleMatrix1D viewRow(int row) {
    _checkRow(row);
    return new DelegateDoubleMatrix1D(this, row);
  }

  DoubleMatrix2D viewRowFlip() {
    if (_rows == 0) return this;
    WrapperDoubleMatrix2D view = new ViewRowFlipWrapperDoubleMatrix2D(this);
    view._isNoView = false;
    return view;
  }

  DoubleMatrix2D viewSelection(Int32List rowIndexes, Int32List columnIndexes) {
    // check for "all"
    if (rowIndexes == null) {
      rowIndexes = new Int32List(_rows);
      for (int i = _rows; --i >= 0; ) rowIndexes[i] = i;
    }
    if (columnIndexes == null) {
      columnIndexes = new Int32List(_columns);
      for (int i = _columns; --i >= 0; ) columnIndexes[i] = i;
    }

    _checkRowIndexes(rowIndexes);
    _checkColumnIndexes(columnIndexes);
    final Int32List rix = rowIndexes;
    final Int32List cix = columnIndexes;

    WrapperDoubleMatrix2D view = new ViewSelectionWrapperDoubleMatrix2D(this, cix, rix);
    view._rows = rowIndexes.length;
    view._columns = columnIndexes.length;
    view._isNoView = false;

    return view;
  }

  DoubleMatrix2D viewStrides(final int _rowStride, final int _columnStride) {
    if (_rowStride <= 0 || _columnStride <= 0) {
      throw new RangeError("illegal stride");
    }
    WrapperDoubleMatrix2D view = new ViewStridesWrapperDoubleMatrix2D(this, _rowStride, _columnStride);
    if (_rows != 0) {
      view._rows = (_rows - 1) ~/ _rowStride + 1;
    }
    if (_columns != 0) {
      view._columns = (_columns - 1) ~/ _columnStride + 1;
    }
    view._isNoView = false;

    return view;
  }

  DoubleMatrix2D _getContent() {
    return _content;
  }

  DoubleMatrix1D _like1D(int size, int offset, int stride) {
    throw new Error(); // should never get called
  }

  DoubleMatrix2D _viewSelectionLike(Int32List rowOffsets, Int32List columnOffsets) {
    throw new Error(); // should never be called
  }

  Object clone() {
    return new WrapperDoubleMatrix2D(_content);
  }
}

class ViewColumnFlipWrapperDoubleMatrix2D extends WrapperDoubleMatrix2D {

  ViewColumnFlipWrapperDoubleMatrix2D(DoubleMatrix2D newContent) : super(newContent);

  double getQuick(int row, int column) {
      return _content.getQuick(row, _columns - 1 - column);
  }

  void setQuick(int row, int column, double value) {
      _content.setQuick(row, _columns - 1 - column, value);
  }

  double get(int row, int column) {
      return _content.get(row, _columns - 1 - column);
  }

  void set(int row, int column, double value) {
      _content.set(row, _columns - 1 - column, value);
  }

  Object clone() {
    return new ViewColumnFlipWrapperDoubleMatrix2D(_content);
  }
}

class ViewDiceWrapperDoubleMatrix2D extends WrapperDoubleMatrix2D {

  ViewDiceWrapperDoubleMatrix2D(DoubleMatrix2D newContent) : super(newContent);

  double getQuick(int row, int column) {
      return _content.getQuick(column, row);
  }

  void setQuick(int row, int column, double value) {
      _content.setQuick(column, row, value);
  }

  double get(int row, int column) {
      return _content.get(column, row);
  }

  void set(int row, int column, double value) {
      _content.set(column, row, value);
  }

  Object clone() {
    return new ViewDiceWrapperDoubleMatrix2D(_content);
  }
}

class ViewPartWrapperDoubleMatrix2D extends WrapperDoubleMatrix2D {

  final int _row;
  final int _column;

  ViewPartWrapperDoubleMatrix2D(DoubleMatrix2D newContent, this._row, this._column) : super(newContent);

  double getQuick(int i, int j) {
      return _content.getQuick(_row + i, _column + j);
  }

  void setQuick(int i, int j, double value) {
      _content.setQuick(_row + i, _column + j, value);
  }

  double get(int i, int j) {
      return _content.get(_row + i, _column + j);
  }

  void set(int i, int j, double value) {
      _content.set(_row + i, _column + j, value);
  }

  Object clone() {
    return new ViewPartWrapperDoubleMatrix2D(_content, _row, _column);
  }
}

class ViewRowFlipWrapperDoubleMatrix2D extends WrapperDoubleMatrix2D {

  ViewRowFlipWrapperDoubleMatrix2D(DoubleMatrix2D newContent) : super(newContent);

  double getQuick(int row, int column) {
      return _content.getQuick(_rows - 1 - row, column);
  }

  void setQuick(int row, int column, double value) {
      _content.setQuick(_rows - 1 - row, column, value);
  }

  double get(int row, int column) {
      return _content.get(_rows - 1 - row, column);
  }

  void set(int row, int column, double value) {
      _content.set(_rows - 1 - row, column, value);
  }

  Object clone() {
    return new ViewRowFlipWrapperDoubleMatrix2D(_content);
  }
}

class ViewSelectionWrapperDoubleMatrix2D extends WrapperDoubleMatrix2D {

  final Int32List _cix;
  final Int32List _rix;

  ViewSelectionWrapperDoubleMatrix2D(DoubleMatrix2D newContent, this._cix, this._rix) : super(newContent);

  double getQuick(int i, int j) {
      return _content.getQuick(_rix[i], _cix[j]);
  }

  void setQuick(int i, int j, double value) {
      _content.setQuick(_rix[i], _cix[j], value);
  }

  double get(int i, int j) {
      return _content.get(_rix[i], _cix[j]);
  }

  void set(int i, int j, double value) {
      _content.set(_rix[i], _cix[j], value);
  }

  Object clone() {
    return new ViewSelectionWrapperDoubleMatrix2D(_content, _cix, _rix);
  }
}

class ViewStridesWrapperDoubleMatrix2D extends WrapperDoubleMatrix2D {

  final int _rowStride;
  final int _columnStride;

  ViewStridesWrapperDoubleMatrix2D(DoubleMatrix2D newContent, this._rowStride, this._columnStride) : super(newContent);

  double getQuick(int row, int column) {
      return _content.getQuick(_rowStride * row, _columnStride * column);
  }

  void setQuick(int row, int column, double value) {
      _content.setQuick(_rowStride * row, _columnStride * column, value);
  }

  double get(int row, int column) {
      return _content.get(_rowStride * row, _columnStride * column);
  }

  void set(int row, int column, double value) {
      _content.set(_rowStride * row, _columnStride * column, value);
  }

  Object clone() {
    return new ViewStridesWrapperDoubleMatrix2D(_content, _rowStride, _columnStride);
  }
}
