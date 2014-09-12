/*
Copyright (C) 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
is hereby granted without fee, provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear in supporting documentation.
CERN makes no representations about the suitability of this software for any purpose.
It is provided "as is" without expressed or implied warranty.
*/
part of cern.colt.matrix;

/**
 * Dense 2-d matrix holding <tt>int</tt> elements. First see the <a
 * href="package-summary.html">package summary</a> and javadoc <a
 * href="package-tree.html">tree view</a> to get the broad picture.
 * <p>
 * <b>Implementation:</b>
 * <p>
 * This data structure allows to store more than 2^31 elements. Internally holds
 * one two-dimensional array, elements[rows][columns]. Note that this
 * implementation is not synchronized.
 * <p>
 * <b>Time complexity:</b>
 * <p>
 * <tt>O(1)</tt> (i.e. constant time) for the basic operations <tt>get</tt>,
 * <tt>getQuick</tt>, <tt>set</tt>, <tt>setQuick</tt> and <tt>size</tt>.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
class LargeIntMatrix extends WrapperIntMatrix {

  List<List<int>> _elements;

  LargeIntMatrix(int rows, int columns) : super(null) {
    try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if ("matrix too large" != exc.message) {
        throw exc;
      }
    }
    _elements = new List<List<int>>.generate(rows,
        (_) => new Int32List(columns));
    _content = this;
  }

  int get(int row, int column) {
    return _elements[row][column];
  }

  void set(int row, int column, int value) {
    _elements[row][column] = value;
  }

  Object get elements => _elements;

  AbstractIntMatrix _getContent() {
    return this;
  }

  AbstractIntMatrix like(int rows, int columns) {
    return new LargeIntMatrix(rows, columns);
  }

  AbstractIntVector like1D(int size) {
    return new IntVector(size);
  }
}
