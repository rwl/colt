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

/// 2-d matrix holding [Complex] elements; either a view wrapping another
/// matrix or a matrix whose views are wrappers.
class WrapperComplexMatrix extends ComplexMatrix {
  ComplexMatrix _content;

  WrapperComplexMatrix._(int rows, int columns) : super._(rows, columns);

  WrapperComplexMatrix._wrap(ComplexMatrix newContent)
      : super._(newContent.rows, newContent.columns) {
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

  dynamic get elements => _content.elements;

  bool equals(ComplexMatrix obj) {
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
      for (int i = 0; i < dlength; i++) {
        var x = new Complex(elements[2 * i], elements[2 * i + 1]);
        var value = new Complex(otherElements[2 * i], otherElements[2 * i + 1]);
        var diff = new Complex((value.real - x.real).abs(), (value.imaginary - x.imaginary).abs());
        if (((diff.real != diff.real) || (diff.imaginary != diff.imaginary)) &&
                ((((value.real != value.real) || (value.imaginary != value.imaginary)) &&
                    ((x.real != x.real) || (x.imaginary != x.imaginary)))) ||
            (isEqual(value, x, epsilon))) {
          diff = Complex.ZERO;
        }
        if ((diff.real > epsilon) || (diff.imaginary > epsilon)) {
          return false;
        }
      }
      return true;
    } else {
      return super.equals(obj);
    }
  }

  Complex get(int row, int column) => _content.get(row, column);

  ComplexMatrix like2D(int rows, int columns) {
    return _content.like2D(rows, columns);
  }

  ComplexVector like1D(int size) => _content.like1D(size);

  void set(int row, int column, Complex value) {
    _content.set(row, column, value);
  }

  void setParts(int row, int column, double re, double im) {
    _content.setParts(row, column, re, im);
  }

  ComplexVector column(int column) => dice().row(column);

  ComplexMatrix columnFlip() {
    if (columns == 0) {
      return this;
    }
    var view = new ColumnFlipWrapperComplexMatrix(this);
    setIsNoView(view, false);
    return view;
  }

  ComplexMatrix dice() {
    var view = new DiceWrapperComplexMatrix(this);
    setRows(view, columns);
    setColumns(view, rows);
    setIsNoView(view, false);
    return view;
  }

  ComplexMatrix part(
      final int row, final int column, int height, int width) {
    checkBox(this, row, column, height, width);
    var view = new PartWrapperComplexMatrix(this, row, column);
    setRows(view, height);
    setColumns(view, width);
    setIsNoView(view, false);
    return view;
  }

  ComplexVector row(int row) {
    checkRow(this, row);
    return new DelegateComplexVector(this, row);
  }

  ComplexMatrix rowFlip() {
    if (rows == 0) {
      return this;
    }
    WrapperComplexMatrix view = new RowWrapperComplexMatrix(this);
    setIsNoView(view, false);
    return view;
  }

  ComplexMatrix select(Int32List rowIndexes, Int32List columnIndexes) {
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

  ComplexMatrix strides(final int rowStride, final int columnStride) {
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

  ComplexMatrix _getContent() => _content;

  ComplexVector _like1D(int size, int offset, int stride) {
    throw new Error(); // should never get called
  }

  ComplexMatrix _viewSelectionLike(
      Int32List rowOffsets, Int32List columnOffsets) {
    throw new Error(); // should never be called
  }

  DoubleMatrix imaginary() {
    var Im = new LargeDoubleMatrix(rows, columns);
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        Im.set(r, c, get(r, c).imaginary);
      }
    }
    return Im;
  }

  DoubleMatrix real() {
    var Re = new LargeDoubleMatrix(rows, columns);
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        Re.set(r, c, get(r, c).real);
      }
    }
    return Re;
  }

  Object clone() => new WrapperComplexMatrix._wrap(_content);
}

class ColumnFlipWrapperComplexMatrix extends WrapperComplexMatrix {
  ColumnFlipWrapperComplexMatrix(ComplexMatrix newContent)
      : super._wrap(newContent);

  Complex get(int row, int column) {
    return _content.get(row, columns - 1 - column);
  }

  void set(int row, int column, Complex value) {
    _content.set(row, columns - 1 - column, value);
  }

  void setParts(int row, int column, double re, double im) {
    _content.setParts(row, columns - 1 - column, re, im);
  }

  Complex at(int row, int column) {
    return _content.at(row, columns - 1 - column);
  }

  void put(int row, int column, Complex value) {
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
  DiceWrapperComplexMatrix(ComplexMatrix newContent)
      : super._wrap(newContent);

  Complex get(int row, int column) {
    return _content.get(column, row);
  }

  void set(int row, int column, Complex value) {
    _content.set(column, row, value);
  }

  void setParts(int row, int column, double re, double im) {
    _content.setParts(column, row, re, im);
  }

  Complex at(int row, int column) {
    return _content.at(column, row);
  }

  void put(int row, int column, Complex value) {
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
      ComplexMatrix newContent, int row, int column)
      : super._wrap(newContent),
        _row = row,
        _column = column;

  Complex get(int i, int j) {
    return _content.get(_row + i, _column + j);
  }

  void set(int i, int j, Complex value) {
    _content.set(_row + i, _column + j, value);
  }

  void setParts(int i, int j, double re, double im) {
    _content.setParts(_row + i, _column + j, re, im);
  }

  Complex at(int i, int j) {
    return _content.at(_row + i, _column + j);
  }

  void put(int i, int j, Complex value) {
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
  RowWrapperComplexMatrix(ComplexMatrix newContent)
      : super._wrap(newContent);

  Complex get(int row, int column) {
    return _content.get(rows - 1 - row, column);
  }

  void set(int row, int column, Complex value) {
    _content.set(rows - 1 - row, column, value);
  }

  void setParts(int row, int column, double re, double im) {
    _content.setParts(rows - 1 - row, column, re, im);
  }

  Complex at(int row, int column) {
    return _content.at(rows - 1 - row, column);
  }

  void put(int row, int column, Complex value) {
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
      ComplexMatrix newContent, Int32List cix, Int32List rix)
      : super._wrap(newContent),
        cix = cix,
        rix = rix;

  Complex get(int i, int j) {
    return _content.get(rix[i], cix[j]);
  }

  void set(int i, int j, Complex value) {
    _content.set(rix[i], cix[j], value);
  }

  void setParts(int i, int j, double re, double im) {
    _content.setParts(rix[i], cix[j], re, im);
  }

  Complex at(int i, int j) {
    return _content.at(rix[i], cix[j]);
  }

  void put(int i, int j, Complex value) {
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
      ComplexMatrix newContent, int rowStride, int columnStride)
      : super._wrap(newContent),
        _rowStride = rowStride,
        _columnStride = columnStride;

  Complex get(int row, int column) {
    return _content.get(_rowStride * row, _columnStride * column);
  }

  void set(int row, int column, Complex value) {
    _content.set(_rowStride * row, _columnStride * column, value);
  }

  void setParts(int row, int column, double re, double im) {
    _content.setParts(_rowStride * row, _columnStride * column, re, im);
  }

  Complex at(int row, int column) {
    return _content.at(_rowStride * row, _columnStride * column);
  }

  void put(int row, int column, Complex value) {
    _content.put(_rowStride * row, _columnStride * column, value);
  }

  void putParts(int row, int column, double re, double im) {
    _content.putParts(_rowStride * row, _columnStride * column, re, im);
  }

  Object clone() {
    return new StridesWrapperComplexMatrix(_content, _rowStride, _columnStride);
  }
}
