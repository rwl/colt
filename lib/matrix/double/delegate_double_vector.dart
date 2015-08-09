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

/// 1-d matrix holding [double] elements; a view wrapping another 2-d
/// matrix and therefore delegating calls to it.
class DelegateDoubleVector extends DoubleVector {
  DoubleMatrix _content;

  int _row;

  DelegateDoubleVector(DoubleMatrix newContent, int row)
      : super._(newContent.columns) {
    if (row < 0 || row >= newContent.rows) {
      throw new ArgumentError();
    }
    _row = row;
    _content = newContent;
  }

  double get(int index) => _content.get(_row, index);

  DoubleVector like1D(int size) => _content.like1D(size);

  DoubleMatrix like2D(int rows, int columns) {
    return _content.like2D(rows, columns);
  }

  void set(int index, double value) => _content.set(_row, index, value);

  dynamic get elements => _content.elements;

  DoubleMatrix reshape(int rows, int columns) {
    throw new ArgumentError("This method is not supported.");
  }

  DoubleVector _viewSelectionLike(Int32List offsets) {
    throw new Error(); // should never get called
  }

  Object clone() => new DelegateDoubleVector(_content, _row);
}
