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
part of cern.colt.matrix.format;

/// Flexible, human readable matrix print formatting; By default decimal
/// point aligned.
///
/// A column can be broader than specified by the parameter [minColumnWidth]
/// (because a cell may not fit into that width) but a column is never smaller
/// than [minColumnWidth]. Normally one does not need to specify
/// [minColumnWidth] (default is `1`). This parameter is only interesting
/// when wanting to print two distinct matrices such that both matrices have
/// the same column width, for example, to make it easier to see which column
/// of matrix A corresponds to which column of matrix B.
///
/// Analyzes the entire matrix before producing output. Each cell is converted
/// to a [String] as indicated by the given format string. If `null` is
/// passed as format string, [double.toString()] is used instead, yielding
/// full precision.
///
/// Next, leading and trailing whitespaces are removed. For each column the
/// maximum number of characters before and after the decimal point is
/// determined. (No problem if decimal points are missing). Each cell is then
/// padded with leading and trailing blanks, as necessary to achieve decimal
/// point aligned, left justified formatting.
class DoubleFormatter extends AbstractFormatter {
  DoubleFormatter([String format = null /*"%G"*/]) {
    _format = format;
    _alignment = AbstractFormatter.DECIMAL;
  }

  /// Converts a given cell to a String; no alignment considered.
  String _formDoubleVector(
      AbstractDoubleVector matrix, int index, Former formatter) {
    return formatter.formDouble(matrix.get(index));
  }

  /// Converts a given cell to a String; no alignment considered.
  String _form(AbstractVector matrix, int index, Former formatter) {
    return this._formDoubleVector(
        matrix as AbstractDoubleVector, index, formatter);
  }

  /// Returns a string representations of all cells; no alignment considered.
  List<List<String>> format(AbstractDoubleMatrix matrix) {
    List<List<String>> strings =
        new List<List<String>>(matrix.rows); //[matrix.columns()];
    for (int row = matrix.rows; --row >= 0;) {
      strings[row] = _formatRow(matrix.row(row));
    }
    return strings;
  }

  /// Returns a string representations of all cells; no alignment considered.
  List<List<String>> _formatMatrix(AbstractMatrix matrix) {
    return format(matrix as AbstractDoubleMatrix);
  }

  /// Returns the index of the decimal point.
  int _indexOfDecimalPoint(String s) {
    int i = s.lastIndexOf('.');
    if (i < 0) i = s.lastIndexOf('e');
    if (i < 0) i = s.lastIndexOf('E');
    if (i < 0) i = s.length;
    return i;
  }

  /// Returns the number of characters before the decimal point.
  int _lead(String s) {
    if (_alignment == AbstractFormatter.DECIMAL) return _indexOfDecimalPoint(s);
    return super._lead(s);
  }

  /// Returns a string representation of the given matrix.
  String toStringDoubleVector(AbstractDoubleVector matrix) {
    AbstractDoubleMatrix easy = matrix.like2D(1, matrix.size);
    easy.row(0).copyFrom(matrix);
    return toString2D(easy);
  }

  /// Returns a string representation of the given matrix.
  String toStringDoubleMatrix(AbstractDoubleMatrix matrix) {
    return super.toStringMatrix(matrix);
  }

  /// Returns a string representation of the given matrix.
  String toString2D(AbstractMatrix matrix) {
    return toStringDoubleMatrix(matrix as AbstractDoubleMatrix);
  }

  /*Object clone() {
    return new DoubleFormatter()
      ..alignment = _alignment
      ..format = _format
      ..minColumnWidth = _minColumnWidth
      ..columnSeparator = _columnSeparator
      ..rowSeparator = _rowSeparator
      ..sliceSeparator = _sliceSeparator
      ..printShape = _printShape;
  }*/
}
