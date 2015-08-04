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

/// Dense 2-d matrix holding [int] elements.
///
/// Internally holds one two-dimensional array, `elements[rows][columns]`.
class LargeIntMatrix extends WrapperIntMatrix {
  List<Int32List> _elements;

  LargeIntMatrix(int rows, int columns) : super._(rows, columns) {
    _elements = new List.generate(rows, (_) => new Int32List(columns));
    _content = this;
  }

  int get(int row, int column) => _elements[row][column];

  void set(int row, int column, int value) {
    _elements[row][column] = value;
  }

  Object get elements => _elements;

  AbstractIntMatrix _getContent() => this;

  AbstractIntMatrix like2D(int rows, int columns) {
    return new LargeIntMatrix(rows, columns);
  }

  AbstractIntVector like1D(int size) => new IntVector(size);
}
