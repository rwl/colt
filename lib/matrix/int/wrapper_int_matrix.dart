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
 * 2-d matrix holding <tt>int</tt> elements; either a view wrapping another
 * matrix or a matrix whose views are wrappers.
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 04/14/2000
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
class WrapperIntMatrix extends AbstractIntMatrix {

  /*
     * The elements of the matrix.
     */
  AbstractIntMatrix _content;

  WrapperIntMatrix(AbstractIntMatrix newContent) {
    if (newContent != null) try {
      _setUp(newContent.rows, newContent.columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if ("matrix too large" != exc.message) {
        throw exc;
      }
    }
    this._content = newContent;
  }

  void forEachWith(final AbstractIntMatrix y, final ifunc.IntIntFunction function) {
    checkShape(y);
    if (y is WrapperIntMatrix) {
      List<int> rowList = new List<int>();
      List<int> columnList = new List<int>();
      List<int> valueList = new List<int>();
      y.nonZeros(rowList, columnList, valueList);
      forEachWithRange(y, function, rowList, columnList);
    } else {
      super.forEachWith(y, function);
    }
    return;
  }

  void setAll(final Int32List values) {
    if (_content is DiagonalIntMatrix) {
      int dlength = (_content as DiagonalIntMatrix)._dlength;
      final Int32List elems = (_content as DiagonalIntMatrix)._elements;
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
      return;
    } else {
      super.setAll(values);
      return;
    }
  }

  Object elements() {
    return _content.elements();
  }

  int get(int row, int column) {
    return _content.get(row, column);
  }

  bool equalsValue(int value) {
    if (_content is DiagonalIntMatrix) {
      Int32List elements = _content.elements() as Int32List;
      for (int r = 0; r < elements.length; r++) {
        int x = elements[r];
        int diff = value - x;
        if (diff != 0) {
          return false;
        }
      }
      return true;
    } else {
      return super.equalsValue(value);
    }
  }

  bool equals(Object obj) {
    if (_content is DiagonalIntMatrix && obj is DiagonalIntMatrix) {
      if (this == obj) {
        return true;
      }
      if (!(this != null && obj != null)) {
        return false;
      }
      DiagonalIntMatrix A = _content as DiagonalIntMatrix;
      DiagonalIntMatrix B = obj;
      if (A.columns != B.columns || A.rows != B.rows || A.diagonalIndex != B.diagonalIndex || A.diagonalLength != B.diagonalLength) {
        return false;
      }
      Int32List AElements = A.elements();
      Int32List BElements = B.elements();
      for (int r = 0; r < AElements.length; r++) {
        int x = AElements[r];
        int value = BElements[r];
        int diff = value - x;
        if (diff != 0) {
          return false;
        }
      }
      return true;
    } else {
      return super.equals(obj);
    }
  }

  AbstractIntMatrix like2D(int rows, int columns) {
    return _content.like2D(rows, columns);
  }

  AbstractIntVector like1D(int size) {
    return _content.like1D(size);
  }

  void set(int row, int column, int value) {
    _content.set(row, column, value);
  }

  AbstractIntVector vectorize() {
    final IntVector v = new IntVector(length);
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (length >= ConcurrencyUtils.getThreadsBeginN_2D())) {
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
              v.set(idx++, get(r, c));
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

  AbstractIntVector column(int column) {
    return dice().row(column);
  }

  AbstractIntMatrix columnFlip() {
    if (_columns == 0) return this;
    WrapperIntMatrix view = new FlipWrapperIntMatrix2D(this);
    view._isNoView = false;

    return view;
  }

  AbstractIntMatrix dice() {
    WrapperIntMatrix view = new DiceWrapperIntMatrix2D(this);
    view._rows = _columns;
    view._columns = _rows;
    view._isNoView = false;

    return view;
  }

  AbstractIntMatrix part(final int row, final int column, int height, int width) {
    _checkBox(row, column, height, width);
    WrapperIntMatrix view = new PartWrapperIntMatrix2D(this, column, row);
    view._rows = height;
    view._columns = width;
    view._isNoView = false;

    return view;
  }

  AbstractIntVector row(int row) {
    _checkRow(row);
    return new DelegateIntVector(this, row);
  }

  AbstractIntMatrix rowFlip() {
    if (_rows == 0) return this;
    WrapperIntMatrix view = new RowFlipWrapperIntMatrix2D(this);
    view._isNoView = false;
    return view;
  }

  AbstractIntMatrix select(Int32List rowIndexes, Int32List columnIndexes) {
    // check for "all"
    if (rowIndexes == null) {
      rowIndexes = new Int32List(_rows);
      for (int i = _rows; --i >= 0; ) {
        rowIndexes[i] = i;
      }
    }
    if (columnIndexes == null) {
      columnIndexes = new Int32List(_columns);
      for (int i = _columns; --i >= 0; ) {
        columnIndexes[i] = i;
      }
    }

    _checkRowIndexes(rowIndexes);
    _checkColumnIndexes(columnIndexes);
    final Int32List rix = rowIndexes;
    final Int32List cix = columnIndexes;

    WrapperIntMatrix view = new SelectWrapperIntMatrix2D(this, cix, rix);
    view._rows = rowIndexes.length;
    view._columns = columnIndexes.length;
    view._isNoView = false;

    return view;
  }

  AbstractIntMatrix strides(final int _rowStride, final int _columnStride) {
    if (_rowStride <= 0 || _columnStride <= 0) {
      throw new RangeError("illegal stride");
    }
    WrapperIntMatrix view = new StridesWrapperIntMatrix2D(this);
    if (_rows != 0) {
      view._rows = (_rows - 1) ~/ _rowStride + 1;
    }
    if (_columns != 0) {
      view._columns = (_columns - 1) ~/ _columnStride + 1;
    }
    view._isNoView = false;

    return view;
  }

  AbstractIntMatrix _getContent() {
    return _content;
  }

  AbstractIntVector _like1D(int size, int offset, int stride) {
    throw new Error(); // should never get called
  }

  AbstractIntMatrix _viewSelectionLike(Int32List rowOffsets, Int32List columnOffsets) {
    throw new Error(); // should never be called
  }
}

class StridesWrapperIntMatrix2D extends WrapperIntMatrix {

  StridesWrapperIntMatrix2D(AbstractIntMatrix newContent) : super(newContent);

  int get(int row, int column) {
    return _content.get(_rowStride * row, _columnStride * column);
  }

  void set(int row, int column, int value) {
    _content.set(_rowStride * row, _columnStride * column, value);
  }

  int at(int row, int column) {
    return _content.at(_rowStride * row, _columnStride * column);
  }

  void put(int row, int column, int value) {
    _content.put(_rowStride * row, _columnStride * column, value);
  }
}

class SelectWrapperIntMatrix2D extends WrapperIntMatrix {
  final Int32List cix;
  final Int32List rix;

  SelectWrapperIntMatrix2D(AbstractIntMatrix newContent, Int32List cix, Int32List rix) : super(newContent),
    cix = cix,
    rix = rix;

  int get(int i, int j) {
    return _content.get(rix[i], cix[j]);
  }

  void set(int i, int j, int value) {
    _content.set(rix[i], cix[j], value);
  }

  int at(int i, int j) {
    return _content.at(rix[i], cix[j]);
  }

  void put(int i, int j, int value) {
    _content.put(rix[i], cix[j], value);
  }
}

class RowFlipWrapperIntMatrix2D extends WrapperIntMatrix {

  RowFlipWrapperIntMatrix2D(AbstractIntMatrix newContent) : super(newContent);

  int get(int row, int column) {
    return _content.get(_rows - 1 - row, column);
  }

  void set(int row, int column, int value) {
    _content.set(_rows - 1 - row, column, value);
  }

  int at(int row, int column) {
    return _content.at(_rows - 1 - row, column);
  }

  void put(int row, int column, int value) {
    _content.put(_rows - 1 - row, column, value);
  }
}

class PartWrapperIntMatrix2D extends WrapperIntMatrix {
  final int __column;
  final int __row;

  PartWrapperIntMatrix2D(AbstractIntMatrix newContent, int column, int row) : super(newContent),
    __column = column,
    __row = row;


  int get(int i, int j) {
    return _content.get(__row + i, __column + j);
  }

  void set(int i, int j, int value) {
    _content.set(__row + i, __column + j, value);
  }

  int at(int i, int j) {
    return _content.at(__row + i, __column + j);
  }

  void put(int i, int j, int value) {
    _content.put(__row + i, __column + j, value);
  }
}

class DiceWrapperIntMatrix2D extends WrapperIntMatrix {

  DiceWrapperIntMatrix2D(AbstractIntMatrix newContent) : super(newContent);

  int get(int row, int column) {
    return _content.get(column, row);
  }

  void set(int row, int column, int value) {
    _content.set(column, row, value);
  }

  int at(int row, int column) {
    return _content.at(column, row);
  }

  void put(int row, int column, int value) {
    _content.put(column, row, value);
  }
}

class FlipWrapperIntMatrix2D extends WrapperIntMatrix {

  FlipWrapperIntMatrix2D(AbstractIntMatrix newContent) : super(newContent);

  int get(int row, int column) {
    return _content.get(row, _columns - 1 - column);
  }

  void set(int row, int column, int value) {
    _content.set(row, _columns - 1 - column, value);
  }

  int at(int row, int column) {
    return _content.at(row, _columns - 1 - column);
  }

  void put(int row, int column, int value) {
    _content.put(row, _columns - 1 - column, value);
  }
}
