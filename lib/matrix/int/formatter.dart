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

/// Flexible, well human readable matrix print formatting; By default decimal
/// point aligned.
///
/// A column can be broader than specified by the parameter [minColumnWidth]
/// (because a cell may not fit into that width) but a column is never smaller
/// than [minColumnWidth]. Normally one does not need to specify
/// [minColumnWidth] (default is `1`). This parameter is only interesting when
/// wanting to print two distinct matrices such that both matrices have the
/// same column width, for example, to make it easier to see which column of
/// matrix A corresponds to which column of matrix B.
///
/// Analyzes the entire matrix before producing output. Each cell is converted
/// to a [String] as indicated by the given format string. If `null` is passed
/// as format string, [int.toString] is used instead.
class IntFormatter extends AbstractFormatter {
  IntFormatter([String format = null]) {
    if (format == null) {
      format = "%d";
    }
    _format = format;
    _alignment = AbstractFormatter.DECIMAL;
  }

  /// Converts a given cell to a String; no alignment considered.
  String _form(AbstractIntVector matrix, int index, Former formatter) {
    return formatter.formInt(matrix[index]);
  }

  /// Converts a given cell to a String; no alignment considered.
  /*String _formVector(AbstractVector matrix, int index, Former formatter) {
    return this._form(matrix as AbstractIntVector, index, formatter);
  }*/

  /// Returns a string representations of all cells; no alignment considered.
  List<List<String>> formatMatrix(AbstractIntMatrix matrix) {
    final strings = new List<List<String>>.generate(
        matrix.rows, (_) => new List<String>(matrix.columns));
    for (int row = matrix.rows; --row >= 0;) {
      strings[row] = _formatRow(matrix.row(row));
    }
    return strings;
  }

  /// Returns a string representations of all cells; no alignment considered.
  List<List<String>> _formatMatrix(AbstractMatrix matrix) {
    return formatMatrix(matrix as AbstractIntMatrix);
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
    if (_alignment == AbstractFormatter.DECIMAL) {
      return _indexOfDecimalPoint(s);
    }
    return super._lead(s);
  }

  /// Returns a string representation of the given matrix.
  String toStringVector(AbstractIntVector matrix) {
    AbstractIntMatrix easy = matrix.like2D(1, matrix.size);
    easy.row(0).copyFrom(matrix);
    return toStringMatrix(easy);
  }

  /// Returns a string representation of the given matrix.
  String toStringMatrix(AbstractIntMatrix matrix) {
    return super.toStringMatrix(matrix);
  }

  /// Returns a string representation of the given matrix.
  /*String _toStringMatrix(AbstractMatrix matrix) {
    return toStringMatrix(matrix as AbstractIntMatrix);
  }*/

  Object clone() {
    return new IntFormatter()
      .._alignment = _alignment
      .._format = _format
      .._minColumnWidth = _minColumnWidth
      .._columnSeparator = _columnSeparator
      .._rowSeparator = _rowSeparator
      .._sliceSeparator = _sliceSeparator
      .._printShape = _printShape;
  }
}
