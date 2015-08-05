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

/// 1-d matrix holding `complex` elements; either a view wrapping another
/// 2-d matrix and therefore delegating calls to it.
class DelegateComplexVector extends AbstractComplexVector {
  AbstractComplexMatrix _content;

  /// The row this view is bound to.
  int _row;

  DelegateComplexVector(AbstractComplexMatrix newContent, int row)
      : super(newContent.columns) {
    if (row < 0 || row >= newContent.rows) {
      throw new ArgumentError();
    }
    _row = row;
    _content = newContent;
  }

  Float64List get(int index) => _content.get(_row, index);

  AbstractComplexVector like1D(int size) => _content.like1D(size);

  AbstractComplexMatrix like2D(int rows, int columns) {
    return _content.like2D(rows, columns);
  }

  void set(int index, Float64List value) {
    _content.set(_row, index, value);
  }

  void setParts(int index, double re, double im) {
    _content.setParts(_row, index, re, im);
  }

  Object get elements => _content.elements;

  // This method is not supported.
  AbstractComplexMatrix reshape(int rows, int columns) {
    throw new ArgumentError("This method is not supported.");
  }

  // This method is not supported.
  AbstractComplexVector _viewSelectionLike(Int32List offsets) {
    throw new ArgumentError("This method is not supported.");
  }

  AbstractDoubleVector imaginary() => _content.row(_row).imaginary();

  AbstractDoubleVector real() => _content.row(_row).real();

  Object clone() => new DelegateComplexVector(_content, _row);
}
