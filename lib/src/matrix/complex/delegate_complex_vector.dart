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

/// 1-d matrix holding [Complex] elements; either a view wrapping another
/// 2-d matrix and therefore delegating calls to it.
class DelegateComplexVector extends ComplexVector {
  ComplexMatrix _content;

  /// The row this view is bound to.
  int _row;

  DelegateComplexVector(ComplexMatrix newContent, int row)
      : super._(newContent.columns) {
    if (row < 0 || row >= newContent.rows) {
      throw new ArgumentError();
    }
    _row = row;
    _content = newContent;
  }

  Complex get(int index) => _content.get(_row, index);

  ComplexVector like1D(int size) => _content.like1D(size);

  ComplexMatrix like2D(int rows, int columns) {
    return _content.like2D(rows, columns);
  }

  void set(int index, Complex value) {
    _content.set(_row, index, value);
  }

  void setParts(int index, double re, double im) {
    _content.setParts(_row, index, re, im);
  }

  Object get elements => _content.elements;

  // This method is not supported.
  ComplexMatrix reshape(int rows, int columns) {
    throw new ArgumentError("This method is not supported.");
  }

  // This method is not supported.
  ComplexVector _viewSelectionLike(Int32List offsets) {
    throw new ArgumentError("This method is not supported.");
  }

  DoubleVector imaginary() => _content.row(_row).imaginary();

  DoubleVector real() => _content.row(_row).real();

  Object clone() => new DelegateComplexVector(_content, _row);
}
