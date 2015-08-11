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

/// Dense 2-d matrix holding [double] elements.
///
/// Internally holds one two-dimensional array, `elements[rows][columns]`.
class LargeDoubleMatrix extends WrapperDoubleMatrix {
  List<Float64List> _elements;

  factory LargeDoubleMatrix(int rows, int columns) {
    var elements = new List.generate(rows, (_) => new Float64List(columns));
    return new LargeDoubleMatrix._internal(rows, columns, elements);
  }

  LargeDoubleMatrix._internal(int rows, int columns, List<Float64List> elements)
      : super._(rows, columns) {
    _elements = elements;
    _content = this;
  }

  double get(int row, int column) => _elements[row][column];

  void set(int row, int column, double value) {
    _elements[row][column] = value;
  }

  List<Float64List> get elements => _elements;

  DoubleMatrix _getContent() => this;

  DoubleMatrix like2D(int rows, int columns) {
    return new LargeDoubleMatrix(rows, columns);
  }

  DoubleVector like1D(int size) => new DenseDoubleVector(size);

  Object clone() {
    return new LargeDoubleMatrix._internal(rows, columns, _elements);
  }
}
