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
 * 2-d matrix holding <tt>complex</tt> elements; either a view wrapping another
 * matrix or a matrix whose views are wrappers.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
class WrapperComplexMatrix extends AbstractComplexMatrix {

  /*
     * The elements of the matrix.
     */
  AbstractComplexMatrix _content;

  WrapperComplexMatrix(AbstractComplexMatrix newContent) {
    if (newContent != null) {
      _setUp(newContent.rows, newContent.columns);
    }
    this._content = newContent;
  }

  void setAll(final Float64List values) {
    if (_content is DiagonalComplexMatrix) {
      int dlength = (_content as DiagonalComplexMatrix)._dlength;
      final Float64List elems = (_content as DiagonalComplexMatrix)._elements;
      if (values.length != 2 * dlength) {
        throw new ArgumentError("Must have same length: length=${values.length} 2 * dlength=${2 * dlength}");
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
              elems[2 * i] = values[2 * i];
              elems[2 * i + 1] = values[2 * i + 1];
            }
          });
        }
        ConcurrencyUtils.waitForCompletion(futures);
      } else {*/
      for (int i = 0; i < dlength; i++) {
        elems[2 * i] = values[2 * i];
        elems[2 * i + 1] = values[2 * i + 1];
      }
      //}
      return;
    } else {
      super.setAll(values);
      return;
    }
  }

  Object elements() {
    return _content.elements();
  }

  /*bool equalsValue(Float64List value) {
    if (_content is DiagonalComplexMatrix) {
      double epsilon = ComplexProperty.DEFAULT.tolerance();
      Float64List elements = _content.elements() as Float64List;
      int dlength = (_content as DiagonalComplexMatrix)._dlength;
      Float64List x = new Float64List(2);
      Float64List diff = new Float64List(2);
      for (int i = 0; i < dlength; i++) {
        x[0] = elements[2 * i];
        x[1] = elements[2 * i + 1];
        diff[0] = (value[0] - x[0]).abs();
        diff[1] = (value[1] - x[1]).abs();
        if (((diff[0] != diff[0]) || (diff[1] != diff[1])) && ((((value[0] != value[0]) || (value[1] != value[1])) && ((x[0] != x[0]) || (x[1] != x[1])))) || (Complex.isEqual(value, x, epsilon))) {
          diff[0] = 0.0;
          diff[1] = 0.0;
        }
        if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
          return false;
        }
      }
      return true;
    } else {
      return super ==(value);
    }
  }*/

  bool operator ==(var obj) {
    if (obj is Float64List) {
      final value = obj;
      if (_content is DiagonalComplexMatrix) {
        double epsilon = ComplexProperty.DEFAULT.tolerance();
        Float64List elements = _content.elements() as Float64List;
        int dlength = (_content as DiagonalComplexMatrix)._dlength;
        Float64List x = new Float64List(2);
        Float64List diff = new Float64List(2);
        for (int i = 0; i < dlength; i++) {
          x[0] = elements[2 * i];
          x[1] = elements[2 * i + 1];
          diff[0] = (value[0] - x[0]).abs();
          diff[1] = (value[1] - x[1]).abs();
          if (((diff[0] != diff[0]) || (diff[1] != diff[1])) && ((((value[0] != value[0]) || (value[1] != value[1])) && ((x[0] != x[0]) || (x[1] != x[1])))) || (Complex.isEqual(value, x, epsilon))) {
            diff[0] = 0.0;
            diff[1] = 0.0;
          }
          if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
            return false;
          }
        }
        return true;
      } else {
        return super == (value);
      }
    }
    if (_content is DiagonalComplexMatrix && obj is DiagonalComplexMatrix) {
      DiagonalComplexMatrix other = obj;
      int dlength = (_content as DiagonalComplexMatrix)._dlength;
      double epsilon = ComplexProperty.DEFAULT.tolerance();
      if (identical(this, obj)) {
        return true;
      }
      if (!(this != null && obj != null)) {
        return false;
      }
      DiagonalComplexMatrix A = _content as DiagonalComplexMatrix;
      DiagonalComplexMatrix B = obj;
      if (A.columns != B.columns || A.rows != B.rows || A.diagonalIndex != B.diagonalIndex || A.diagonalLength != B.diagonalLength) {
        return false;
      }
      Float64List otherElements = other._elements;
      Float64List elements = (_content as DiagonalComplexMatrix)._elements;
      Float64List x = new Float64List(2);
      Float64List value = new Float64List(2);
      Float64List diff = new Float64List(2);
      for (int i = 0; i < dlength; i++) {
        x[0] = elements[2 * i];
        x[1] = elements[2 * i + 1];
        value[0] = otherElements[2 * i];
        value[1] = otherElements[2 * i + 1];
        diff[0] = (value[0] - x[0]).abs();
        diff[1] = (value[1] - x[1]).abs();
        if (((diff[0] != diff[0]) || (diff[1] != diff[1])) && ((((value[0] != value[0]) || (value[1] != value[1])) && ((x[0] != x[0]) || (x[1] != x[1])))) || (Complex.isEqual(value, x, epsilon))) {
          diff[0] = 0.0;
          diff[1] = 0.0;
        }
        if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
          return false;
        }
      }
      return true;
    } else {
      return super == (obj);
    }
  }

  Float64List get(int row, int column) {
    return _content.get(row, column);
  }

  AbstractComplexMatrix like2D(int rows, int columns) {
    return _content.like2D(rows, columns);
  }

  AbstractComplexVector like1D(int size) {
    return _content.like1D(size);
  }

  void set(int row, int column, Float64List value) {
    _content.set(row, column, value);
  }

  void setParts(int row, int column, double re, double im) {
    _content.setParts(row, column, re, im);
  }

  AbstractComplexVector vectorize() {
    final ComplexVector v = new ComplexVector(this.length);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _columns);
      List<Future> futures = new List<Future>(nthreads);
      int k = _columns ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstColumn = j * k;
        final int lastColumn = (j == nthreads - 1) ? _columns : firstColumn + k;
        futures[j] = ConcurrencyUtils.submit(() {
          int idx = firstColumn * _rows;
          for (int c = firstColumn; c < lastColumn; c++) {
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
        v.set(idx++, get(r, c));
      }
    }
    //}
    return v;
  }

  AbstractComplexVector column(int column) {
    return dice().row(column);
  }

  AbstractComplexMatrix columnFlip() {
    if (_columns == 0) return this;
    WrapperComplexMatrix view = new ColumnFlipWrapperComplexMatrix(this);
    view._isNoView = false;
    return view;
  }

  AbstractComplexMatrix dice() {
    WrapperComplexMatrix view = new DiceWrapperComplexMatrix(this);
    view._rows = _columns;
    view._columns = _rows;
    view._isNoView = false;

    return view;
  }

  AbstractComplexMatrix part(final int row, final int column, int height, int width) {
    _checkBox(row, column, height, width);
    WrapperComplexMatrix view = new PartWrapperComplexMatrix(this, row, column);
    view._rows = height;
    view._columns = width;
    view._isNoView = false;

    return view;
  }

  AbstractComplexVector row(int row) {
    _checkRow(row);
    return new DelegateComplexVector(this, row);
  }

  AbstractComplexMatrix rowFlip() {
    if (_rows == 0) return this;
    WrapperComplexMatrix view = new RowWrapperComplexMatrix(this);
    view._isNoView = false;

    return view;
  }

  AbstractComplexMatrix select(Int32List rowIndexes, Int32List columnIndexes) {
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

    WrapperComplexMatrix view = new SelectionWrapperComplexMatrix(this, cix, rix);
    view._rows = rowIndexes.length;
    view._columns = columnIndexes.length;
    view._isNoView = false;

    return view;
  }

  AbstractComplexMatrix strides(final int rowStride, final int columnStride) {
    if (rowStride <= 0 || columnStride <= 0) {
      throw new RangeError("illegal stride");
    }
    WrapperComplexMatrix view = new StridesWrapperComplexMatrix(this, rowStride, columnStride);
    if (_rows != 0) {
      view._rows = (_rows - 1) ~/ rowStride + 1;
    }
    if (_columns != 0) {
      view._columns = (_columns - 1) ~/ columnStride + 1;
    }
    view._isNoView = false;

    return view;
  }

  AbstractComplexMatrix _getContent() {
    return _content;
  }

  AbstractComplexVector _like1D(int size, int offset, int stride) {
    throw new Error(); // should never get called
  }

  AbstractComplexMatrix _viewSelectionLike(Int32List rowOffsets, Int32List columnOffsets) {
    throw new Error(); // should never be called
  }

  AbstractDoubleMatrix imaginary() {
    final LargeDoubleMatrix Im = new LargeDoubleMatrix(_rows, _columns);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, _rows);
      List<Future> futures = new List<Future>(nthreads);
      int k = _rows ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? _rows : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = 0; c < _columns; c++) {
              Im.setQuick(r, c, getQuick(r, c)[1]);
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    for (int r = 0; r < _rows; r++) {
      for (int c = 0; c < _columns; c++) {
        Im.set(r, c, get(r, c)[1]);
      }
    }
    //}
    return Im;
  }

  AbstractDoubleMatrix real() {
    final LargeDoubleMatrix Re = new LargeDoubleMatrix(_rows, _columns);
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
              Re.setQuick(r, c, getQuick(r, c)[0]);
            }
          }
        });
      }
      ConcurrencyUtils.waitForCompletion(futures);
    } else {*/
    for (int r = 0; r < _rows; r++) {
      for (int c = 0; c < _columns; c++) {
        Re.set(r, c, get(r, c)[0]);
      }
    }
    //}
    return Re;
  }

  Object clone() {
    return new WrapperComplexMatrix(_content);
  }
}

class ColumnFlipWrapperComplexMatrix extends WrapperComplexMatrix {

  ColumnFlipWrapperComplexMatrix(AbstractComplexMatrix newContent) : super(newContent);

  Float64List get(int row, int column) {
    return _content.get(row, _columns - 1 - column);
  }

  void set(int row, int column, Float64List value) {
    _content.set(row, _columns - 1 - column, value);
  }

  void setParts(int row, int column, double re, double im) {
    _content.setParts(row, _columns - 1 - column, re, im);
  }

  Float64List at(int row, int column) {
    return _content.at(row, _columns - 1 - column);
  }

  void put(int row, int column, Float64List value) {
    _content.put(row, _columns - 1 - column, value);
  }

  void putParts(int row, int column, double re, double im) {
    _content.putParts(row, _columns - 1 - column, re, im);
  }

  Object clone() {
    return new ColumnFlipWrapperComplexMatrix(_content);
  }
}

class DiceWrapperComplexMatrix extends WrapperComplexMatrix {

  DiceWrapperComplexMatrix(AbstractComplexMatrix newContent) : super(newContent);

  Float64List get(int row, int column) {
    return _content.get(column, row);
  }

  void set(int row, int column, Float64List value) {
    _content.set(column, row, value);
  }

  void setParts(int row, int column, double re, double im) {
    _content.setParts(column, row, re, im);
  }

  Float64List at(int row, int column) {
    return _content.at(column, row);
  }

  void put(int row, int column, Float64List value) {
    _content.put(column, row, value);
  }

  void putParts(int row, int column, double re, double im) {
    _content.putParts(column, row, re, im);
  }

  Object clone() {
    return new DiceWrapperComplexMatrix(_content);
  }
}

class PartWrapperComplexMatrix extends WrapperComplexMatrix {
  final int _row;
  final int _column;

  PartWrapperComplexMatrix(AbstractComplexMatrix newContent, int row, int column)
      : super(newContent),
        _row = row,
        _column = column;

  Float64List get(int i, int j) {
    return _content.get(_row + i, _column + j);
  }

  void set(int i, int j, Float64List value) {
    _content.set(_row + i, _column + j, value);
  }

  void setParts(int i, int j, double re, double im) {
    _content.setParts(_row + i, _column + j, re, im);
  }

  Float64List at(int i, int j) {
    return _content.at(_row + i, _column + j);
  }

  void put(int i, int j, Float64List value) {
    _content.put(_row + i, _column + j, value);
  }

  void putParts(int i, int j, double re, double im) {
    _content.putParts(_row + i, _column + j, re, im);
  }

  Object clone() {
    return new PartWrapperComplexMatrix(_content, _row, _column);
  }
}

class RowWrapperComplexMatrix extends WrapperComplexMatrix {

  RowWrapperComplexMatrix(AbstractComplexMatrix newContent) : super(newContent);

  Float64List get(int row, int column) {
    return _content.get(_rows - 1 - row, column);
  }

  void set(int row, int column, Float64List value) {
    _content.set(_rows - 1 - row, column, value);
  }

  void setParts(int row, int column, double re, double im) {
    _content.setParts(_rows - 1 - row, column, re, im);
  }

  Float64List at(int row, int column) {
    return _content.at(_rows - 1 - row, column);
  }

  void put(int row, int column, Float64List value) {
    _content.put(_rows - 1 - row, column, value);
  }

  void putParts(int row, int column, double re, double im) {
    _content.putParts(_rows - 1 - row, column, re, im);
  }

  Object clone() {
    return new RowWrapperComplexMatrix(_content);
  }
}

class SelectionWrapperComplexMatrix extends WrapperComplexMatrix {
  final Int32List cix;
  final Int32List rix;

  SelectionWrapperComplexMatrix(AbstractComplexMatrix newContent, Int32List cix, Int32List rix)
      : super(newContent),
        cix = cix,
        rix = rix;

  Float64List get(int i, int j) {
    return _content.get(rix[i], cix[j]);
  }

  void set(int i, int j, Float64List value) {
    _content.set(rix[i], cix[j], value);
  }

  void setParts(int i, int j, double re, double im) {
    _content.setParts(rix[i], cix[j], re, im);
  }

  Float64List at(int i, int j) {
    return _content.at(rix[i], cix[j]);
  }

  void put(int i, int j, Float64List value) {
    _content.put(rix[i], cix[j], value);
  }

  void putParts(int i, int j, double re, double im) {
    _content.putParts(rix[i], cix[j], re, im);
  }

  Object clone() {
    return new SelectionWrapperComplexMatrix(_content, cix, rix);
  }
}

class StridesWrapperComplexMatrix extends WrapperComplexMatrix {

  StridesWrapperComplexMatrix(AbstractComplexMatrix newContent, int rowStride, int columnStride) : super(newContent) {
    _rowStride = rowStride;
    _columnStride = columnStride;
  }

  Float64List get(int row, int column) {
    return _content.get(_rowStride * row, _columnStride * column);
  }

  void set(int row, int column, Float64List value) {
    _content.set(_rowStride * row, _columnStride * column, value);
  }

  void setParts(int row, int column, double re, double im) {
    _content.setParts(_rowStride * row, _columnStride * column, re, im);
  }

  Float64List at(int row, int column) {
    return _content.at(_rowStride * row, _columnStride * column);
  }

  void put(int row, int column, Float64List value) {
    _content.put(_rowStride * row, _columnStride * column, value);
  }

  void putParts(int row, int column, double re, double im) {
    _content.putParts(_rowStride * row, _columnStride * column, re, im);
  }

  Object clone() {
    return new StridesWrapperComplexMatrix(_content, _rowStride, _columnStride);
  }
}
