/*
Copyright (C) 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
is hereby granted without fee, provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear in supporting documentation.
CERN makes no representations about the suitability of this software for any purpose.
It is provided "as is" without expressed or implied warranty.
 */
library cern.colt.matrix.format;

import 'dart:math' as Math;
import 'package:intl/intl.dart';
import 'matrix.dart';
import 'double/matrix.dart';
import 'int/matrix.dart';
import '../math.dart' show MIN_INT;

part 'double/formatter.dart';
part 'int/formatter.dart';

/// Formats a double or complex (double[]) into a string (like sprintf in C).
abstract class Former {
  /// Formats a double into a string (like sprintf in C).
  String formDouble(double value);

  /// Formats an int into a string (like sprintf in C).
  String formInt(int value);

  /// Formats a complex (double[]) into a string (like sprintf in C).
  String formComplex(List<double> value);
}

class DefaultFormer implements Former {
  NumberFormat format;

  DefaultFormer(String fmt) {
    format = new NumberFormat(fmt);
  }

  String formDouble(double value) => format.format(value);

  String formInt(int value) => format.format(value);

  String formComplex(List<double> value) {
    if (value[0] == 0 && value[1] == 0) {
      return "0";
    }
    if (value[1] == 0) {
      return format.format(value[0]);
    }
    if (value[1] < 0) {
      return format.format(value[0]) + " - " + format.format(-value[1]);
    }
    return format.format(value[0]) + " + " + format.format(value[1]);
  }
}

/// Abstract base class for flexible, human readable matrix print formatting.
/// Value type independent. A single cell is formatted via a format string.
/// Columns can be aligned left, centered, right and by decimal point.
abstract class AbstractFormatter {

  /// The alignment string aligning the cells of a column to the left.
  static final String LEFT = "left";

  /// The alignment string aligning the cells of a column to its center.
  static final String CENTER = "center";

  /// The alignment string aligning the cells of a column to the right.
  static final String RIGHT = "right";

  /// The alignment string aligning the cells of a column to the decimal point.
  static final String DECIMAL = "decimal";

  /// The default minimum number of characters a column may have; currently
  /// `1`.
  static final int DEFAULT_MIN_COLUMN_WIDTH = 1;

  /// The default string separating any two columns from another; currently
  /// `" "`.
  static final String DEFAULT_COLUMN_SEPARATOR = " ";

  /// The default string separating any two rows from another; currently
  /// `"\n"`.
  static final String DEFAULT_ROW_SEPARATOR = "\n";

  /// The default string separating any two slices from another; currently
  /// `"\n\n"`.
  static final String DEFAULT_SLICE_SEPARATOR = "\n\n";

  String _alignment = LEFT;

  String _format = null;//"%G";

  int _minColumnWidth = DEFAULT_MIN_COLUMN_WIDTH;

  String _columnSeparator = DEFAULT_COLUMN_SEPARATOR;

  String _rowSeparator = DEFAULT_ROW_SEPARATOR;

  String _sliceSeparator = DEFAULT_SLICE_SEPARATOR;

  /// Tells whether String representations are to be preceded with summary of
  /// the shape; currently `true`.
  bool _printShape = true;

  static List<String> _blanksCache = []; // for efficient String manipulations

  Function _factory = (String fmt) => new DefaultFormer(fmt);

  /*static {
    _setupBlanksCache();
  }*/

  /// Modifies the strings in a column of the string matrix to be aligned
  /// (left,centered,right,decimal).
  void _align(List<List<String>> strings) {
    int rows = strings.length;
    int columns = 0;
    if (rows > 0) {
      columns = strings[0].length;
    }

    var maxColWidth = new List<int>(columns);
    var maxColLead = null;
    bool isDecimal = _alignment == DECIMAL;
    if (isDecimal) {
      maxColLead = new List<int>(columns);
    }
    // var maxColTrail = new List<int>(columns);

    // for each column, determine alignment parameters
    for (int column = 0; column < columns; column++) {
      int maxWidth = _minColumnWidth;
      int maxLead = MIN_INT;
      // int maxTrail = Integer.MIN_VALUE;
      for (int row = 0; row < rows; row++) {
        String s = strings[row][column];
        maxWidth = Math.max(maxWidth, s.length);
        if (isDecimal) maxLead = Math.max(maxLead, _lead(s));
        // maxTrail = Math.max(maxTrail, trail(s));
      }
      maxColWidth[column] = maxWidth;
      if (isDecimal) maxColLead[column] = maxLead;
      // maxColTrail[column] = maxTrail;
    }

    // format each row according to alignment parameters
    // var total = new StringBuffer();
    for (int row = 0; row < rows; row++) {
      _alignRow(strings[row], maxColWidth, maxColLead);
    }

  }

  /// Converts a row into a string.
  /*int _alignmentCode(String alignment) {
    // {-1,0,1,2} = {left,centered,right,decimal point}
    if (alignment == LEFT) {
      return -1; }
    else if (alignment == CENTER) {
      return 0;
    } else if (alignment == RIGHT) {
      return 1;
    } else if (alignment == DECIMAL) {
      return 2;
    } else {
      throw new ArgumentError("unknown alignment: " + alignment);
    }
  }*/

  /// Modifies the strings the string matrix to be aligned
  /// (left,centered,right,decimal).
  void _alignRow(List<String> row, List<int> maxColWidth, List<int> maxColLead) {
    var s = new StringBuffer();

    int columns = row.length;
    for (int column = 0; column < columns; column++) {
      s.clear();
      String c = row[column];
      if (_alignment == RIGHT) {
        s.write(_blanks(maxColWidth[column] - s.length));
        s.write(c);
      } else if (_alignment == DECIMAL) {
        s.write(_blanks(maxColLead[column] - _lead(c)));
        s.write(c);
        s.write(_blanks(maxColWidth[column] - s.length));
      } else if (_alignment == CENTER) {
        s.write(_blanks((maxColWidth[column] - c.length) ~/ 2));
        s.write(c);
        s.write(_blanks(maxColWidth[column] - s.length));

      } else if (_alignment == LEFT) {
        s.write(c);
        s.write(_blanks(maxColWidth[column] - s.length));
      } else throw new Error();

      row[column] = s.toString();
    }
  }

  /// Returns a String with `length` blanks.
  String _blanks(int length) {
    if (length < 0) {
      length = 0;
    }
    if (length < _blanksCache.length) {
      return _blanksCache[length];
    }

    var buf = new StringBuffer(length);
    for (int k = 0; k < length; k++) {
      buf.write(' ');
    }
    return buf.toString();
  }

  /// Converts a given cell to a String; no alignment considered.
  String _form(AbstractVector matrix, int index, Former formatter);

  /// Returns a string representations of all cells; no alignment considered.
  List<List<String>> _formatMatrix(AbstractMatrix matrix);

  /// Returns a string representations of all cells; no alignment considered.
  List<String> _formatRow(AbstractVector vector) {
    Former formatter = null;
    formatter = _factory(_format);
    int s = vector.size;
    var strings = new List<String>(s);
    for (int i = 0; i < s; i++) {
      strings[i] = _form(vector, i, formatter);
    }
    return strings;
  }

  /// Returns the number of characters or the number of characters before the
  /// decimal point.
  int _lead(String s) {
    return s.length;
  }

  /// Returns a String with the given character repeated `length` times.
  /*String _repeat(String character, int length) {
    if (character == ' ') return _blanks(length);
    if (length < 0) length = 0;
    StringBuffer buf = new StringBuffer(length);
    for (int k = 0; k < length; k++) {
      buf.write(character);
    }
    return buf.toString();
  }*/

  /// Sets the column alignment (left,center,right,decimal).
  void set alignment(String alignment) {
    _alignment = alignment;
  }

  /// Sets the string separating any two columns from another.
  void set columnSeparator(String columnSeparator) {
    _columnSeparator = columnSeparator;
  }

  /// Sets the way a single cell value is to be formatted.
  void set format(String format) {
    _format = format;
  }

  /// Sets the minimum number of characters a column may have.
  void set minColumnWidth(int minColumnWidth) {
    if (minColumnWidth < 0) {
      throw new ArgumentError();
    }
    _minColumnWidth = minColumnWidth;
  }

  /// Specifies whether a string representation of a matrix is to be preceded
  /// with a summary of its shape.
  void set printShape(bool printShape) {
    _printShape = printShape;
  }

  /// Sets the string separating any two rows from another.
  void set rowSeparator(String rowSeparator) {
    _rowSeparator = rowSeparator;
  }

  /// Sets the string separating any two slices from another.
  void set sliceSeparator(String sliceSeparator) {
    _sliceSeparator = sliceSeparator;
  }

  /// Cache for faster string processing.
  /*static void _setupBlanksCache() {
    // Pre-fabricate 40 static strings with 0,1,2,..,39 blanks, for usage
    // within method blanks(length).
    // Now, we don't need to construct and fill them on demand, and garbage
    // collect them again.
    // All 40 strings share the identical char[] array, only with different
    // offset and length --> somewhat smaller static memory footprint
    int size = 40;
    _blanksCache = new List<String>(size);
    var buf = new StringBuffer(size);
    for (int i = size; --i >= 0; ) {
      buf.write(' ');
    }
    String str = buf.toString();
    for (int i = size; --i >= 0; ) {
      _blanksCache[i] = str.substring(0, i);
      // System.out.println(i+"-"+blanksCache[i]+"-");
    }
  }*/

  /// Returns a short string representation describing the shape of the matrix.
  static String shapeVector(AbstractVector matrix) {
    // return "Vector of size="+matrix.size();
    // return matrix.size()+" element matrix";
    // return "matrix("+matrix.size()+")";
    return "${matrix.size} vector";
  }

  /// Returns a short string representation describing the shape of the matrix.
  static String shapeMatrix(AbstractMatrix matrix) {
    return "${matrix.rows} x ${matrix.columns} matrix";
  }

  /// Returns a single string representation of the given string matrix.
  String _toString(List<List<String>> strings) {
    int rows = strings.length;
    int columns = strings.length <= 0 ? 0 : strings[0].length;

    var total = new StringBuffer();
    var s = new StringBuffer();
    for (int row = 0; row < rows; row++) {
      s.clear();
      for (int column = 0; column < columns; column++) {
        s.write(strings[row][column]);
        if (column < columns - 1) {
          s.write(_columnSeparator);
        }
      }
      total.write(s);
      if (row < rows - 1) {
        total.write(_rowSeparator);
      }
    }

    return total.toString();
  }

  /// Returns a string representation of the given matrix.
  String toStringMatrix(AbstractMatrix matrix) {
    List<List<String>> strings = _formatMatrix(matrix);
    _align(strings);
    String total = _toString(strings);
    if (_printShape) {
      "${shapeMatrix(matrix)}\n$total";
    }
    return total.toString();
  }
}
