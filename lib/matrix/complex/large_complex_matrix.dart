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

/// Dense 2-d matrix holding `complex` elements.
///
/// Internally holds one two-dimensional array, `elements[rows][2*columns]`.
/// Complex data is represented by 2 double values in sequence.
class LargeComplexMatrix extends WrapperComplexMatrix {
  List<Float64List> _elements;

  LargeComplexMatrix(int rows, int columns) : super._(rows, columns) {
    _elements = new List.generate(rows, (_) => new Float64List(2 * columns));
    _content = this;
  }

  Float64List get(int row, int column) {
    return new Float64List.fromList(
        [_elements[row][2 * column], _elements[row][2 * column + 1]]);
  }

  void set(int row, int column, Float64List value) {
    _elements[row][2 * column] = value[0];
    _elements[row][2 * column + 1] = value[1];
  }

  void setParts(int row, int column, double re, double im) {
    _elements[row][2 * column] = re;
    _elements[row][2 * column + 1] = im;
  }

  Object get elements => _elements;

  AbstractComplexMatrix _getContent() => this;

  AbstractComplexMatrix like2D(int rows, int columns) {
    return new LargeComplexMatrix(rows, columns);
  }

  AbstractComplexVector like1D(int size) => new ComplexVector(size);
}
