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
 * Abstract base class for flexible, well human readable matrix print
 * formatting. Value type independent. A single cell is formatted via a format
 * string. Columns can be aligned left, centered, right and by decimal point.
 * <p>
 * A column can be broader than specified by the parameter
 * <tt>minColumnWidth</tt> (because a cell may not fit into that width) but a
 * column is never smaller than <tt>minColumnWidth</tt>. Normally one does not
 * need to specify <tt>minColumnWidth</tt>. Cells in a row are separated by a
 * separator string, similar separators can be set for rows and slices. For more
 * info, see the concrete subclasses.
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 09/24/99
 */
abstract class AbstractFormatter {//extends cern.colt.PersistentObject {

    /**
     * The alignment string aligning the cells of a column to the left.
     */
    static final String LEFT = "left";

    /**
     * The alignment string aligning the cells of a column to its center.
     */
    static final String CENTER = "center";

    /**
     * The alignment string aligning the cells of a column to the right.
     */
    static final String RIGHT = "right";

    /**
     * The alignment string aligning the cells of a column to the decimal point.
     */
    static final String DECIMAL = "decimal";

    /**
     * The default minimum number of characters a column may have; currently
     * <tt>1</tt>.
     */
    static final int DEFAULT_MIN_COLUMN_WIDTH = 1;

    /**
     * The default string separating any two columns from another; currently
     * <tt>" "</tt>.
     */
    static final String DEFAULT_COLUMN_SEPARATOR = " ";

    /**
     * The default string separating any two rows from another; currently
     * <tt>"\n"</tt>.
     */
    static final String DEFAULT_ROW_SEPARATOR = "\n";

    /**
     * The default string separating any two slices from another; currently
     * <tt>"\n\n"</tt>.
     */
    static final String DEFAULT_SLICE_SEPARATOR = "\n\n";

    /**
     * The default format string for formatting a single cell value; currently
     * <tt>"%G"</tt>.
     */
    String _alignment = LEFT;

    /**
     * The default format string for formatting a single cell value; currently
     * <tt>"%G"</tt>.
     */
    String _format = "%G";

    /**
     * The default minimum number of characters a column may have; currently
     * <tt>1</tt>.
     */
    int _minColumnWidth = DEFAULT_MIN_COLUMN_WIDTH;

    /**
     * The default string separating any two columns from another; currently
     * <tt>" "</tt>.
     */
    String _columnSeparator = DEFAULT_COLUMN_SEPARATOR;

    /**
     * The default string separating any two rows from another; currently
     * <tt>"\n"</tt>.
     */
    String _rowSeparator = DEFAULT_ROW_SEPARATOR;

    /**
     * The default string separating any two slices from another; currently
     * <tt>"\n\n"</tt>.
     */
    String _sliceSeparator = DEFAULT_SLICE_SEPARATOR;

    /**
     * Tells whether String representations are to be preceded with summary of
     * the shape; currently <tt>true</tt>.
     */
    bool _printShape = true;

    static List<String> _blanksCache; // for efficient String manipulations

    static final FormerFactory _factory = new FormerFactory();

    /*static {
        _setupBlanksCache();
    }*/

    /**
     * Makes this class non instantiable, but still let's others inherit from
     * it.
     */
    AbstractFormatter() {
    }

    /**
     * Modifies the strings in a column of the string matrix to be aligned
     * (left,centered,right,decimal).
     */
    void _align(List<List<String>> strings) {
        int rows = strings.length;
        int columns = 0;
        if (rows > 0)
            columns = strings[0].length;

        List<int> maxColWidth = new List<int>(columns);
        List<int> maxColLead = null;
        bool isDecimal = _alignment == DECIMAL;
        if (isDecimal)
            maxColLead = new List<int>(columns);
        // List<int> maxColTrail = new int[columns];

        // for each column, determine alignment parameters
        for (int column = 0; column < columns; column++) {
            int maxWidth = _minColumnWidth;
            int maxLead = -576460752303423487;//Integer.MIN_VALUE;
            // int maxTrail = Integer.MIN_VALUE;
            for (int row = 0; row < rows; row++) {
                String s = strings[row][column];
                maxWidth = math.max(maxWidth, s.length);
                if (isDecimal)
                    maxLead = math.max(maxLead, _lead(s));
                // maxTrail = Math.max(maxTrail, trail(s));
            }
            maxColWidth[column] = maxWidth;
            if (isDecimal)
                maxColLead[column] = maxLead;
            // maxColTrail[column] = maxTrail;
        }

        // format each row according to alignment parameters
        // StringBuffer total = new StringBuffer();
        for (int row = 0; row < rows; row++) {
            _alignRow(strings[row], maxColWidth, maxColLead);
        }

    }

    /**
     * Converts a row into a string.
     */
    int _alignmentCode(String alignment) {
        // {-1,0,1,2} = {left,centered,right,decimal point}
        if (alignment == LEFT)
            return -1;
        else if (alignment == CENTER)
            return 0;
        else if (alignment == RIGHT)
            return 1;
        else if (alignment == DECIMAL)
            return 2;
        else
            throw new ArgumentError("unknown alignment: " + alignment);
    }

    /**
     * Modifies the strings the string matrix to be aligned
     * (left,centered,right,decimal).
     */
    void _alignRow(List<String> row, List<int> maxColWidth, List<int> maxColLead) {
        StringBuffer s = new StringBuffer();

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
            } else
                throw new Error();

            row[column] = s.toString();
        }
    }

    /**
     * Returns a String with <tt>length</tt> blanks.
     */
    String _blanks(int length) {
        if (length < 0)
            length = 0;
        if (length < _blanksCache.length)
            return _blanksCache[length];

        StringBuffer buf = new StringBuffer(length);
        for (int k = 0; k < length; k++) {
            buf.write(' ');
        }
        return buf.toString();
    }

    /**
     * Converts a given cell to a String; no alignment considered.
     */
    String _form(AbstractMatrix1D matrix, int index, Former formatter);

    /**
     * Returns a string representations of all cells; no alignment considered.
     */
    List<List<String>> _format2D(AbstractMatrix2D matrix);

    /**
     * Returns a string representations of all cells; no alignment considered.
     */
    List<String> _formatRow(AbstractMatrix1D vector) {
        Former formatter = null;
        formatter = _factory.create(_format);
        int s = vector.size();
        List<String> strings = new List<String>(s);
        for (int i = 0; i < s; i++) {
            strings[i] = _form(vector, i, formatter);
        }
        return strings;
    }

    /**
     * Returns the number of characters or the number of characters before the
     * decimal point.
     */
    int _lead(String s) {
        return s.length;
    }

    /**
     * Returns a String with the given character repeated <tt>length</tt> times.
     */
    String _repeat(String character, int length) {
        if (character == ' ')
            return _blanks(length);
        if (length < 0)
            length = 0;
        StringBuffer buf = new StringBuffer(length);
        for (int k = 0; k < length; k++) {
            buf.write(character);
        }
        return buf.toString();
    }

    /**
     * Sets the column alignment (left,center,right,decimal).
     *
     * @param alignment
     *            the new alignment to be used; must be one of
     *            <tt>{LEFT,CENTER,RIGHT,DECIMAL}</tt>.
     */
    void setAlignment(String alignment) {
        this._alignment = alignment;
    }

    /**
     * Sets the string separating any two columns from another.
     *
     * @param columnSeparator
     *            the new columnSeparator to be used.
     */
    void setColumnSeparator(String columnSeparator) {
        this._columnSeparator = columnSeparator;
    }

    /**
     * Sets the way a <i>single</i> cell value is to be formatted.
     *
     * @param format
     *            the new format to be used.
     */
    void setFormat(String format) {
        this._format = format;
    }

    /**
     * Sets the minimum number of characters a column may have.
     *
     * @param minColumnWidth
     *            the new minColumnWidth to be used.
     */
    void setMinColumnWidth(int minColumnWidth) {
        if (minColumnWidth < 0)
            throw new ArgumentError();
        this._minColumnWidth = minColumnWidth;
    }

    /**
     * Specifies whether a string representation of a matrix is to be preceded
     * with a summary of its shape.
     *
     * @param printShape
     *            <tt>true</tt> shape summary is printed, otherwise not printed.
     */
    void setPrintShape(bool printShape) {
        this._printShape = printShape;
    }

    /**
     * Sets the string separating any two rows from another.
     *
     * @param rowSeparator
     *            the new rowSeparator to be used.
     */
    void setRowSeparator(String rowSeparator) {
        this._rowSeparator = rowSeparator;
    }

    /**
     * Sets the string separating any two slices from another.
     *
     * @param sliceSeparator
     *            the new sliceSeparator to be used.
     */
    void setSliceSeparator(String sliceSeparator) {
        this._sliceSeparator = sliceSeparator;
    }

    /**
     * Cache for faster string processing.
     */
    static void _setupBlanksCache() {
        // Pre-fabricate 40 static strings with 0,1,2,..,39 blanks, for usage
        // within method blanks(length).
        // Now, we don't need to construct and fill them on demand, and garbage
        // collect them again.
        // All 40 strings share the identical char[] array, only with different
        // offset and length --> somewhat smaller static memory footprint
        int size = 40;
        _blanksCache = new List<String>(size);
        StringBuffer buf = new StringBuffer(size);
        for (int i = size; --i >= 0;)
            buf.write(' ');
        String str = buf.toString();
        for (int i = size; --i >= 0;) {
            _blanksCache[i] = str.substring(0, i);
            // System.out.println(i+"-"+blanksCache[i]+"-");
        }
    }

    /**
     * Returns a short string representation describing the shape of the matrix.
     */
    static String shape(AbstractMatrix1D matrix) {
        // return "Matrix1D of size="+matrix.size();
        // return matrix.size()+" element matrix";
        // return "matrix("+matrix.size()+")";
        return "${matrix.size()} matrix";
    }

    /**
     * Returns a short string representation describing the shape of the matrix.
     */
    static String shape2D(AbstractMatrix2D matrix) {
        return "${matrix.rows()} x ${matrix.columns()} matrix";
    }

    /**
     * Returns a short string representation describing the shape of the matrix.
     */
    /*static String shape3D(AbstractMatrix3D matrix) {
        return matrix.slices() + " x " + matrix.rows() + " x " + matrix.columns() + " matrix";
    }*/

    /**
     * Returns a single string representation of the given string matrix.
     *
     * @param strings
     *            the matrix to be converted to a single string.
     */
    String _toString(List<List<String>> strings) {
        int rows = strings.length;
        int columns = strings.length <= 0 ? 0 : strings[0].length;

        StringBuffer total = new StringBuffer();
        StringBuffer s = new StringBuffer();
        for (int row = 0; row < rows; row++) {
            s.clear();
            for (int column = 0; column < columns; column++) {
                s.write(strings[row][column]);
                if (column < columns - 1)
                    s.write(_columnSeparator);
            }
            total.write(s);
            if (row < rows - 1)
                total.write(_rowSeparator);
        }

        return total.toString();
    }

    /**
     * Returns a string representation of the given matrix.
     *
     * @param matrix
     *            the matrix to convert.
     */
    String _toString2D(AbstractMatrix2D matrix) {
        List<List<String>> strings = this._format2D(matrix);
        _align(strings);
        String total = _toString(strings);
        if (_printShape)
            "${shape2D(matrix)}\n$total";
        return total.toString();
    }
}
