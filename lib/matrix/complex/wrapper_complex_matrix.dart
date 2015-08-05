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
part of cern.colt.matrix.complex;

/// 2-d matrix holding `complex` elements; either a view wrapping another
/// matrix or a matrix whose views are wrappers.
class WrapperComplexMatrix extends AbstractComplexMatrix {
  AbstractComplexMatrix _content;

  WrapperComplexMatrix._(int rows, int columns) : super(rows, columns);

  WrapperComplexMatrix._wrap(AbstractComplexMatrix newContent)
      : super(newContent.rows, newContent.columns) {
    _content = newContent;
  }

  void setAll(final Float64List values) {
    if (_content is DiagonalComplexMatrix) {
      int dlength = (_content as DiagonalComplexMatrix)._dlength;
      Float64List elems = (_content as DiagonalComplexMatrix)._elements;
      if (values.length != 2 * dlength) {
        throw new ArgumentError(
            "Must have same length: length=${values.length} 2 * dlength=${2 * dlength}");
      }
      for (int i = 0; i < dlength; i++) {
        elems[2 * i] = values[2 * i];
        elems[2 * i + 1] = values[2 * i + 1];
      }
      return;
    } else {
      super.setAll(values);
      return;
    }
  }

  Object get elements => _content.elements;

  bool equals(AbstractComplexMatrix obj) {
    if (_content is DiagonalComplexMatrix && obj is DiagonalComplexMatrix) {
      DiagonalComplexMatrix other = obj;
      int dlength = (_content as DiagonalComplexMatrix)._dlength;
      double epsilon = EPSILON;
      if (identical(this, obj)) {
        return true;
      }
      if (!(this != null && obj != null)) {
        return false;
      }
      DiagonalComplexMatrix A = _content as DiagonalComplexMatrix;
      DiagonalComplexMatrix B = obj;
      if (A.columns != B.columns ||
          A.rows != B.rows ||
          A.diagonalIndex != B.diagonalIndex ||
          A.diagonalLength != B.diagonalLength) {
        return false;
      }
      Float64List otherElements = other._elements;
      Float64List elements = (_content as DiagonalComplexMatrix)._elements;
      var x = new Float64List(2);
      var value = new Float64List(2);
      var diff = new Float64List(2);
      for (int i = 0; i < dlength; i++) {
        x[0] = elements[2 * i];
        x[1] = elements[2 * i + 1];
        value[0] = otherElements[2 * i];
        value[1] = otherElements[2 * i + 1];
        diff[0] = (value[0] - x[0]).abs();
        diff[1] = (value[1] - x[1]).abs();
        if (((diff[0] != diff[0]) || (diff[1] != diff[1])) &&
                ((((value[0] != value[0]) || (value[1] != value[1])) &&
                    ((x[0] != x[0]) || (x[1] != x[1])))) ||
            (Complex.isEqual(value, x, epsilon))) {
          diff[0] = 0.0;
          diff[1] = 0.0;
        }
        if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
          return false;
        }
      }
      return true;
    } else {
      return super.equals(obj);
    }
  }

  Float64List get(int row, int column) => _content.get(row, column);

  AbstractComplexMatrix like2D(int rows, int columns) {
    return _content.like2D(rows, columns);
  }

  AbstractComplexVector like1D(int size) => _content.like1D(size);

  void set(int row, int column, Float64List value) {
    _content.set(row, column, value);
  }

  void setParts(int row, int column, double re, double im) {
    _content.setParts(row, column, re, im);
  }

  AbstractComplexVector column(int column) => dice().row(column);

  AbstractComplexMatrix columnFlip() {
    if (columns == 0) {
      return this;
    }
    var view = new ColumnFlipWrapperComplexMatrix(this);
    setIsNoView(view, false);
    return view;
  }

  AbstractComplexMatrix dice() {
    var view = new DiceWrapperComplexMatrix(this);
    setRows(view, columns);
    setColumns(view, rows);
    setIsNoView(view, false);
    return view;
  }

  AbstractComplexMatrix part(
      final int row, final int column, int height, int width) {
    checkBox(this, row, column, height, width);
    var view = new PartWrapperComplexMatrix(this, row, column);
    setRows(view, height);
    setColumns(view, width);
    setIsNoView(view, false);
    return view;
  }

  AbstractComplexVector row(int row) {
    checkRow(this, row);
    return new DelegateComplexVector(this, row);
  }

  AbstractComplexMatrix rowFlip() {
    if (rows == 0) {
      return this;
    }
    WrapperComplexMatrix view = new RowWrapperComplexMatrix(this);
    setIsNoView(view, false);
    return view;
  }

  AbstractComplexMatrix select(Int32List rowIndexes, Int32List columnIndexes) {
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
    final rix = new Int32List.fromList(rowIndexes);
    final cix = new Int32List.fromList(columnIndexes);

    var view = new SelectionWrapperComplexMatrix(this, cix, rix);
    setRows(view, rowIndexes.length);
    setColumns(view, columnIndexes.length);
    setIsNoView(view, false);
    return view;
  }

  AbstractComplexMatrix strides(final int rowStride, final int columnStride) {
    if (rowStride <= 0 || columnStride <= 0) {
      throw new RangeError("illegal stride");
    }
    var view = new StridesWrapperComplexMatrix(this, rowStride, columnStride);
    if (rows != 0) {
      setRows(view, (rows - 1) ~/ rowStride + 1);
    }
    if (columns != 0) {
      setColumns(view, (columns - 1) ~/ columnStride + 1);
    }
    setIsNoView(view, false);
    return view;
  }

  AbstractComplexMatrix _getContent() => _content;

  AbstractComplexVector _like1D(int size, int offset, int stride) {
    throw new Error(); // should never get called
  }

  AbstractComplexMatrix _viewSelectionLike(
      Int32List rowOffsets, Int32List columnOffsets) {
    throw new Error(); // should never be called
  }

  AbstractDoubleMatrix imaginary() {
    var Im = new LargeDoubleMatrix(rows, columns);
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        Im.set(r, c, get(r, c)[1]);
      }
    }
    return Im;
  }

  AbstractDoubleMatrix real() {
    var Re = new LargeDoubleMatrix(rows, columns);
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        Re.set(r, c, get(r, c)[0]);
      }
    }
    return Re;
  }

  Object clone() => new WrapperComplexMatrix._wrap(_content);
}

class ColumnFlipWrapperComplexMatrix extends WrapperComplexMatrix {
  ColumnFlipWrapperComplexMatrix(AbstractComplexMatrix newContent)
      : super._wrap(newContent);

  Float64List get(int row, int column) {
    return _content.get(row, columns - 1 - column);
  }

  void set(int row, int column, Float64List value) {
    _content.set(row, columns - 1 - column, value);
  }

  void setParts(int row, int column, double re, double im) {
    _content.setParts(row, columns - 1 - column, re, im);
  }

  Float64List at(int row, int column) {
    return _content.at(row, columns - 1 - column);
  }

  void put(int row, int column, Float64List value) {
    _content.put(row, columns - 1 - column, value);
  }

  void putParts(int row, int column, double re, double im) {
    _content.putParts(row, columns - 1 - column, re, im);
  }

  Object clone() {
    return new ColumnFlipWrapperComplexMatrix(_content);
  }
}

class DiceWrapperComplexMatrix extends WrapperComplexMatrix {
  DiceWrapperComplexMatrix(AbstractComplexMatrix newContent)
      : super._wrap(newContent);

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

  PartWrapperComplexMatrix(
      AbstractComplexMatrix newContent, int row, int column)
      : super._wrap(newContent),
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
  RowWrapperComplexMatrix(AbstractComplexMatrix newContent)
      : super._wrap(newContent);

  Float64List get(int row, int column) {
    return _content.get(rows - 1 - row, column);
  }

  void set(int row, int column, Float64List value) {
    _content.set(rows - 1 - row, column, value);
  }

  void setParts(int row, int column, double re, double im) {
    _content.setParts(rows - 1 - row, column, re, im);
  }

  Float64List at(int row, int column) {
    return _content.at(rows - 1 - row, column);
  }

  void put(int row, int column, Float64List value) {
    _content.put(rows - 1 - row, column, value);
  }

  void putParts(int row, int column, double re, double im) {
    _content.putParts(rows - 1 - row, column, re, im);
  }

  Object clone() {
    return new RowWrapperComplexMatrix(_content);
  }
}

class SelectionWrapperComplexMatrix extends WrapperComplexMatrix {
  final Int32List cix;
  final Int32List rix;

  SelectionWrapperComplexMatrix(
      AbstractComplexMatrix newContent, Int32List cix, Int32List rix)
      : super._wrap(newContent),
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
  final int _rowStride;
  final int _columnStride;

  StridesWrapperComplexMatrix(
      AbstractComplexMatrix newContent, int rowStride, int columnStride)
      : super._wrap(newContent),
        _rowStride = rowStride,
        _columnStride = columnStride;

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
