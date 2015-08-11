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
part of cern.colt.matrix.int;

/// 2-d matrix holding [int] elements; either a view wrapping another
/// matrix or a matrix whose views are wrappers.
class WrapperIntMatrix extends IntMatrix {
  WrapperIntMatrix._(int rows, int columns) : super._(rows, columns);

  IntMatrix _content;

  WrapperIntMatrix._wrap(IntMatrix newContent)
      : super._(newContent.rows, newContent.columns) {
    _content = newContent;
  }

  /*void assign(IntMatrix y, ifunc.IntIntFunction function) {
    checkShape(this, y);
    if (y is WrapperIntMatrix) {
      var rowList = new List();
      var columnList = new List();
      var valueList = new List();
      y.nonzero(rowList, columnList, valueList);
      forEachWithRange(y, function, rowList, columnList);
    } else {
      super.assign(y, function);
    }
    return;
  }*/

  void setAll(Int32List values) {
    if (_content is DiagonalIntMatrix) {
      int dlength = (_content as DiagonalIntMatrix)._dlength;
      Int32List elems = (_content as DiagonalIntMatrix)._elements;
      if (values.length != dlength) {
        throw new ArgumentError(
            "Must have same length: length=${values.length} dlength=$dlength");
      }
      for (int i = 0; i < dlength; i++) {
        elems[i] = values[i];
      }
      return;
    } else {
      super.setAll(values);
      return;
    }
  }

  Object get elements => _content.elements;

  int get(int row, int column) => _content.get(row, column);

  bool equals(IntMatrix obj) {
    if (_content is DiagonalIntMatrix && obj is DiagonalIntMatrix) {
      if (this == obj) {
        return true;
      }
      if (!(this != null && obj != null)) {
        return false;
      }
      DiagonalIntMatrix A = _content as DiagonalIntMatrix;
      DiagonalIntMatrix B = obj;
      if (A.columns != B.columns ||
          A.rows != B.rows ||
          A.diagonalIndex != B.diagonalIndex ||
          A.diagonalLength != B.diagonalLength) {
        return false;
      }
      Int32List AElements = A.elements;
      Int32List BElements = B.elements;
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

  IntMatrix like2D(int rows, int columns) {
    return _content.like2D(rows, columns);
  }

  IntVector like1D(int size) => _content.like1D(size);

  void set(int row, int column, int value) {
    _content.set(row, column, value);
  }

  IntVector column(int column) => dice().row(column);

  IntMatrix columnFlip() {
    if (columns == 0) {
      return this;
    }
    var view = new FlipWrapperIntMatrix2D(this);
    setIsNoView(view, false);

    return view;
  }

  IntMatrix dice() {
    var view = new DiceWrapperIntMatrix2D(this);
    setRows(view, columns);
    setColumns(view, rows);
    setIsNoView(view, false);

    return view;
  }

  IntMatrix part(int row, int column, int height, int width) {
    checkBox(this, row, column, height, width);
    var view = new PartWrapperIntMatrix2D(this, column, row);
    setRows(view, height);
    setColumns(view, width);
    setIsNoView(view, false);

    return view;
  }

  IntVector row(int row) {
    checkRow(this, row);
    return new DelegateIntVector(this, row);
  }

  IntMatrix rowFlip() {
    if (rows == 0) {
      return this;
    }
    var view = new RowFlipWrapperIntMatrix2D(this);
    setIsNoView(view, false);
    return view;
  }

  IntMatrix select(Int32List rowIndexes, Int32List columnIndexes) {
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

    var view = new SelectWrapperIntMatrix2D(this, cix, rix);
    setRows(view, rowIndexes.length);
    setColumns(view, columnIndexes.length);
    setIsNoView(view, false);

    return view;
  }

  IntMatrix strides(int _rowStride, int _columnStride) {
    if (_rowStride <= 0 || _columnStride <= 0) {
      throw new RangeError("illegal stride");
    }
    var view = new StridesWrapperIntMatrix2D(this, _rowStride, _columnStride);
    if (rows != 0) {
      setRows(view, (rows - 1) ~/ _rowStride + 1);
    }
    if (columns != 0) {
      setColumns(view, (columns - 1) ~/ _columnStride + 1);
    }
    setIsNoView(view, false);

    return view;
  }

  IntMatrix _getContent() => _content;

  IntVector _like1D(int size, int offset, int stride) {
    throw new Error(); // should never get called
  }

  IntMatrix _viewSelectionLike(Int32List rowOffsets, Int32List columnOffsets) {
    throw new Error(); // should never be called
  }

  Object clone() => new WrapperIntMatrix._wrap(_content);
}

class StridesWrapperIntMatrix2D extends WrapperIntMatrix {
  final int _rowStride;
  final int _columnStride;

  StridesWrapperIntMatrix2D(
      IntMatrix newContent, this._rowStride, this._columnStride)
      : super._wrap(newContent);

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

  Object clone() =>
      new StridesWrapperIntMatrix2D(_content, _rowStride, _columnStride);
}

class SelectWrapperIntMatrix2D extends WrapperIntMatrix {
  final Int32List cix;
  final Int32List rix;

  SelectWrapperIntMatrix2D(IntMatrix newContent, Int32List cix, Int32List rix)
      : super._wrap(newContent),
        cix = cix,
        rix = rix;

  int get(int i, int j) => _content.get(rix[i], cix[j]);

  void set(int i, int j, int value) => _content.set(rix[i], cix[j], value);

  int at(int i, int j) => _content.at(rix[i], cix[j]);

  void put(int i, int j, int value) => _content.put(rix[i], cix[j], value);

  Object clone() => new SelectWrapperIntMatrix2D(_content, cix, rix);
}

class RowFlipWrapperIntMatrix2D extends WrapperIntMatrix {
  RowFlipWrapperIntMatrix2D(IntMatrix newContent) : super._wrap(newContent);

  int get(int row, int column) => _content.get(rows - 1 - row, column);

  void set(int row, int column, int value) {
    _content.set(rows - 1 - row, column, value);
  }

  int at(int row, int column) => _content.at(rows - 1 - row, column);

  void put(int row, int column, int value) {
    _content.put(rows - 1 - row, column, value);
  }

  Object clone() => new RowFlipWrapperIntMatrix2D(_content);
}

class PartWrapperIntMatrix2D extends WrapperIntMatrix {
  final int __column;
  final int __row;

  PartWrapperIntMatrix2D(IntMatrix newContent, int column, int row)
      : super._wrap(newContent),
        __column = column,
        __row = row;

  int get(int i, int j) => _content.get(__row + i, __column + j);

  void set(int i, int j, int value) {
    _content.set(__row + i, __column + j, value);
  }

  int at(int i, int j) => _content.at(__row + i, __column + j);

  void put(int i, int j, int value) {
    _content.put(__row + i, __column + j, value);
  }

  Object clone() => new PartWrapperIntMatrix2D(_content, __column, __row);
}

class DiceWrapperIntMatrix2D extends WrapperIntMatrix {
  DiceWrapperIntMatrix2D(IntMatrix newContent) : super._wrap(newContent);

  int get(int row, int column) => _content.get(column, row);

  void set(int row, int column, int value) {
    _content.set(column, row, value);
  }

  int at(int row, int column) => _content.at(column, row);

  void put(int row, int column, int value) {
    _content.put(column, row, value);
  }

  Object clone() => new DiceWrapperIntMatrix2D(_content);
}

class FlipWrapperIntMatrix2D extends WrapperIntMatrix {
  FlipWrapperIntMatrix2D(IntMatrix newContent) : super._wrap(newContent);

  int get(int row, int column) => _content.get(row, columns - 1 - column);

  void set(int row, int column, int value) {
    _content.set(row, columns - 1 - column, value);
  }

  int at(int row, int column) => _content.at(row, columns - 1 - column);

  void put(int row, int column, int value) {
    _content.put(row, columns - 1 - column, value);
  }

  Object clone() => new FlipWrapperIntMatrix2D(_content);
}
