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

/**
 * Dense 2-d matrix holding <tt>complex</tt> elements.<br>
 * <b>Implementation:</b>
 * <p>
 * This data structure allows to store more than 2^31 elements. Internally holds
 * one two-dimensional array, elements[rows][2*columns]. Complex data is
 * represented by 2 double values in sequence, i.e. elements[row][2*column]
 * constitute the real part and elements[row][2*column+1] constitute the
 * imaginary part. Note that this implementation is not synchronized.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 *
 */
class LargeComplexMatrix extends WrapperComplexMatrix {

  List<Float64List> _elements;

  LargeComplexMatrix(int rows, int columns) : super(null) {
    try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if ("matrix too large" != exc.message) {
        throw exc;
      }
    }
    _elements = new List<Float64List>.generate(rows, (_) => new Float64List(2 * columns));
    _content = this;
  }

  Float64List get(int row, int column) {
    return new Float64List.fromList([_elements[row][2 * column], _elements[row][2 * column + 1]]);
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

  AbstractComplexMatrix _getContent() {
    return this;
  }

  AbstractComplexMatrix like2D(int rows, int columns) {
    return new LargeComplexMatrix(rows, columns);
  }

  AbstractComplexVector like1D(int size) {
    return new ComplexVector(size);
  }
}
