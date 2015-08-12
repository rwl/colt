// Copyright (C) 1999 CERN - European Organization for Nuclear Research.
//
// Permission to use, copy, modify, distribute and sell this software and
// its documentation for any purpose is hereby granted without fee, provided
// that the above copyright notice appear in all copies and that both that
// copyright notice and this permission notice appear in supporting
// documentation.
//
// CERN makes no representations about the suitability of this software for
// any purpose. It is provided "as is" without expressed or implied warranty.
part of cern.colt.matrix.double;

/**
 * 2-d matrix holding <tt>double</tt> elements; either a view wrapping another
 * matrix or a matrix whose views are wrappers.
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 04/14/2000
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
class WrapperDoubleMatrix extends DoubleMatrix {
  DoubleMatrix _content;

  WrapperDoubleMatrix._(int rows, int columns) : super._(rows, columns);

  WrapperDoubleMatrix._wrap(DoubleMatrix newContent)
      : super._(newContent.rows, newContent.columns) {
    _content = newContent;
  }

  void setValues(Float64List values) {
    if (_content is DiagonalDoubleMatrix) {
      int dlength = (_content as DiagonalDoubleMatrix)._dlength;
      Float64List elems = (_content as DiagonalDoubleMatrix)._elements;
      if (values.length != dlength) {
        throw new ArgumentError(
            "Must have same length: length=${values.length} dlength=$dlength");
      }

      for (int i = 0; i < dlength; i++) {
        elems[i] = values[i];
      }
      return;
    } else {
      super.setValues(values);
    }
  }

  /*void assign(
      DoubleMatrix y, func.DoubleDoubleFunction fn) {
    checkShape(this, y);
    if (y is WrapperDoubleMatrix) {
      var rowList = new List<int>();
      var columnList = new List<int>();
      var valueList = new List<double>();
      y.nonzero(rowList, columnList, valueList);
      forEachWithNonZero(y, fn, new Int32List.fromList(rowList),
          new Int32List.fromList(columnList));
    } else {
      super.assign(y, fn);
    }
  }*/

  Object get elements => _content.elements;

  double get(int row, int column) => _content.get(row, column);

  bool equals(DoubleMatrix obj) {
    if (_content is DiagonalDoubleMatrix && obj is DiagonalDoubleMatrix) {
      double epsilon = EPSILON;
      if (identical(this, obj)) {
        return true;
      }
      if (!(this != null && obj != null)) {
        return false;
      }
      var A = _content as DiagonalDoubleMatrix;
      var B = obj;
      if (A.columns != B.columns ||
          A.rows != B.rows ||
          A.diagonalIndex != B.diagonalIndex ||
          A.diagonalLength != B.diagonalLength) {
        return false;
      }
      Float64List AElements = A.elements;
      Float64List BElements = B.elements;
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
      return super == (obj);
    }
  }

  DoubleMatrix like2D(int rows, int columns) {
    return _content.like2D(rows, columns);
  }

  DoubleVector like1D(int size) => _content.like1D(size);

  void set(int row, int column, double value) {
    _content.set(row, column, value);
  }

  DoubleVector column(int column) => dice().row(column);

  DoubleMatrix columnFlip() {
    if (columns == 0) {
      return this;
    }
    var view = new ColumnFlipWrapperDoubleMatrix(this);
    setIsNoView(view, false);

    return view;
  }

  DoubleMatrix dice() {
    WrapperDoubleMatrix view = new DiceWrapperDoubleMatrix(this);
    setRows(view, columns);
    setColumns(view, rows);
    setIsNoView(view, false);
    return view;
  }

  DoubleMatrix part(int row, int column, int height, int width) {
    checkBox(this, row, column, height, width);
    var view = new PartWrapperDoubleMatrix(this, row, column);
    setRows(view, height);
    setColumns(view, width);
    setIsNoView(view, false);
    return view;
  }

  DoubleVector row(int row) {
    checkRow(this, row);
    return new DelegateDoubleVector(this, row);
  }

  DoubleMatrix rowFlip() {
    if (rows == 0) {
      return this;
    }
    var view = new RowFlipWrapperDoubleMatrix(this);
    setIsNoView(view, false);
    return view;
  }

  DoubleMatrix select(Int32List rowIndexes, Int32List columnIndexes) {
    // check for "all"
    if (rowIndexes == null) {
      rowIndexes = new Int32List(rows);
      for (int i = rows; --i >= 0;) {
        rowIndexes[i] = i;
      }
    }
    if (columnIndexes == null) {
      columnIndexes = new Int32List(columns);
      for (int i = columns; --i >= 0;) {
        columnIndexes[i] = i;
      }
    }

    checkRowIndexes(this, rowIndexes);
    checkColumnIndexes(this, columnIndexes);
    var rix = new Int32List.fromList(rowIndexes);
    var cix = new Int32List.fromList(columnIndexes);

    WrapperDoubleMatrix view = new SelectionWrapperDoubleMatrix(this, cix, rix);
    setRows(view, rowIndexes.length);
    setColumns(view, columnIndexes.length);
    setIsNoView(view, false);

    return view;
  }

  DoubleMatrix strides(int _rowStride, int _columnStride) {
    if (_rowStride <= 0 || _columnStride <= 0) {
      throw new RangeError("illegal stride");
    }
    var view = new StridesWrapperDoubleMatrix(this, _rowStride, _columnStride);
    if (rows != 0) {
      setRows(view, (rows - 1) ~/ _rowStride + 1);
    }
    if (columns != 0) {
      setColumns(view, (columns - 1) ~/ _columnStride + 1);
    }
    setIsNoView(view, false);

    return view;
  }

  DoubleMatrix _getContent() => _content;

  DoubleVector _like1D(int size, int offset, int stride) {
    throw new Error(); // should never get called
  }

  DoubleMatrix _viewSelectionLike(
      Int32List rowOffsets, Int32List columnOffsets) {
    throw new Error(); // should never be called
  }

  Object clone() => new WrapperDoubleMatrix._wrap(_content);
}

class ColumnFlipWrapperDoubleMatrix extends WrapperDoubleMatrix {
  ColumnFlipWrapperDoubleMatrix(DoubleMatrix newContent)
      : super._wrap(newContent);

  double get(int row, int column) => _content.get(row, columns - 1 - column);

  void set(int row, int column, double value) {
    _content.set(row, columns - 1 - column, value);
  }

  double at(int row, int column) {
    return _content.at(row, columns - 1 - column);
  }

  void put(int row, int column, double value) {
    _content.put(row, columns - 1 - column, value);
  }

  Object clone() => new ColumnFlipWrapperDoubleMatrix(_content);
}

class DiceWrapperDoubleMatrix extends WrapperDoubleMatrix {
  DiceWrapperDoubleMatrix(DoubleMatrix newContent) : super._wrap(newContent);

  double get(int row, int column) => _content.get(column, row);

  void set(int row, int column, double value) {
    _content.set(column, row, value);
  }

  double at(int row, int column) => _content.at(column, row);

  void put(int row, int column, double value) {
    _content.put(column, row, value);
  }

  Object clone() => new DiceWrapperDoubleMatrix(_content);
}

class PartWrapperDoubleMatrix extends WrapperDoubleMatrix {
  final int _row;
  final int _column;

  PartWrapperDoubleMatrix(DoubleMatrix newContent, this._row, this._column)
      : super._wrap(newContent);

  double get(int i, int j) => _content.get(_row + i, _column + j);

  void set(int i, int j, double value) {
    _content.set(_row + i, _column + j, value);
  }

  double at(int i, int j) => _content.at(_row + i, _column + j);

  void put(int i, int j, double value) {
    _content.put(_row + i, _column + j, value);
  }

  Object clone() => new PartWrapperDoubleMatrix(_content, _row, _column);
}

class RowFlipWrapperDoubleMatrix extends WrapperDoubleMatrix {
  RowFlipWrapperDoubleMatrix(DoubleMatrix newContent) : super._wrap(newContent);

  double get(int row, int column) => _content.get(rows - 1 - row, column);

  void set(int row, int column, double value) {
    _content.set(rows - 1 - row, column, value);
  }

  double at(int row, int column) => _content.at(rows - 1 - row, column);

  void put(int row, int column, double value) {
    _content.put(rows - 1 - row, column, value);
  }

  Object clone() => new RowFlipWrapperDoubleMatrix(_content);
}

class SelectionWrapperDoubleMatrix extends WrapperDoubleMatrix {
  final Int32List _cix;
  final Int32List _rix;

  SelectionWrapperDoubleMatrix(DoubleMatrix newContent, this._cix, this._rix)
      : super._wrap(newContent);

  double get(int i, int j) => _content.get(_rix[i], _cix[j]);

  void set(int i, int j, double value) {
    _content.set(_rix[i], _cix[j], value);
  }

  double at(int i, int j) => _content.at(_rix[i], _cix[j]);

  void put(int i, int j, double value) {
    _content.put(_rix[i], _cix[j], value);
  }

  Object clone() => new SelectionWrapperDoubleMatrix(_content, _cix, _rix);
}

class StridesWrapperDoubleMatrix extends WrapperDoubleMatrix {
  final int _rowStride;
  final int _columnStride;

  StridesWrapperDoubleMatrix(
      DoubleMatrix newContent, this._rowStride, this._columnStride)
      : super._wrap(newContent);

  double get(int row, int column) {
    return _content.get(_rowStride * row, _columnStride * column);
  }

  void set(int row, int column, double value) {
    _content.set(_rowStride * row, _columnStride * column, value);
  }

  double at(int row, int column) {
    return _content.at(_rowStride * row, _columnStride * column);
  }

  void put(int row, int column, double value) {
    _content.put(_rowStride * row, _columnStride * column, value);
  }

  Object clone() {
    return new StridesWrapperDoubleMatrix(_content, _rowStride, _columnStride);
  }
}
