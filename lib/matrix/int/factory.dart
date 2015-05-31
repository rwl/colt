/*
Copyright (C) 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
is hereby granted without fee, provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear in supporting documentation.
CERN makes no representations about the suitability of this software for any purpose.
It is provided "as is" without expressed or implied warranty.
 */
library cern.colt.matrix.int.factory;

import 'dart:math' as Math;
import 'dart:typed_data';
import 'package:colt/function/int.dart' as ifunc;
import '../matrix.dart';

/**
 * Factory for convenient construction of 1-d matrices holding <tt>int</tt>
 * cells. Use idioms like <tt>IntFactory1D.dense.make(1000)</tt> to construct
 * dense matrices, <tt>IntFactory1D.sparse.make(1000)</tt> to construct sparse
 * matrices.
 *
 * If the factory is used frequently it might be useful to streamline the
 * notation. For example by aliasing:
 * <table>
 * <td class="PRE">
 *
 * <pre>
 *  IntFactory1D F = IntFactory1D.dense;
 *  F.make(1000);
 *  F.descending(10);
 *  F.random(3);
 *  ...
 * </pre>
 *
 * </td>
 * </table>
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 09/24/99
 */
//class IntFactory1D {
//
//  /**
//   * A factory producing dense matrices.
//   */
//  static final IntFactory1D dense = new IntFactory1D._internal();
//
//  /**
//   * A factory producing sparse matrices.
//   */
//  static final IntFactory1D sparse = new IntFactory1D._internal();
//
//  /**
//   * Makes this class non instantiable, but still let's others inherit from
//   * it.
//   */
//  IntFactory1D._internal();

/**
 * C = A||B; Constructs a new matrix which is the concatenation of two other
 * matrices. Example: <tt>0 1</tt> append <tt>3 4</tt> --> <tt>0 1 3 4</tt>.
 */
AbstractIntVector append(AbstractIntVector A, AbstractIntVector B, AbstractIntVector make(int)) {
  // concatenate
  AbstractIntVector matrix = make(A.length + B.length);
  matrix.part(0, A.length).copyFrom(A);
  matrix.part(A.length, B.length).copyFrom(B);
  return matrix;
}

/**
 * Constructs a matrix with cells having ascending values. For debugging
 * purposes. Example: <tt>0 1 2</tt>
 */
AbstractIntVector ascending(int size, AbstractIntVector make(int)) {
  return descending(size, make)..forEach(ifunc.chain2(ifunc.neg, ifunc.subtract(size)));
}

/**
 * Constructs a matrix with cells having descending values. For debugging
 * purposes. Example: <tt>2 1 0</tt>
 */
AbstractIntVector descending(int size, AbstractIntVector make(int)) {
  AbstractIntVector matrix = make(size);
  int v = 0;
  for (int i = size; --i >= 0; ) {
    matrix.set(i, v++);
  }
  return matrix;
}

/**
 * Constructs a matrix with the given cell values. The values are copied. So
 * subsequent changes in <tt>values</tt> are not reflected in the matrix,
 * and vice-versa.
 *
 * @param values
 *            The values to be filled into the new matrix.
 */
//  AbstractIntVector withValues(Int32List values) {
//    if (this == sparse) {
//      return new SparseIntVector.fromList(values);
//    } else {
//      return new IntVector.fromList(values);
//    }
//  }

/**
 * Constructs a matrix which is the concatenation of all given parts. Cells
 * are copied.
 */
AbstractIntVector concat(List<AbstractIntVector> parts, AbstractIntVector make(int)) {
  if (parts.length == 0) {
    return make(0);
  }

  int size = 0;
  for (int i = 0; i < parts.length; i++) {
    size += parts[i].length;
  }

  AbstractIntVector vector = make(size);
  size = 0;
  for (int i = 0; i < parts.length; i++) {
    vector.part(size, parts[i].length).copyFrom(parts[i]);
    size += parts[i].length;
  }

  return vector;
}

/**
 * Constructs a matrix with the given shape, each cell initialized with
 * zero.
 */
//  AbstractIntVector make(int size) {
//    if (this == sparse) {
//      return new SparseIntVector(size);
//    }
//    return new IntVector(size);
//  }

/**
 * Constructs a matrix with the given shape, each cell initialized with the
 * given value.
 */
AbstractIntVector fill(int size, int initialValue, AbstractIntVector make(int)) {
  return make(size)..fill(initialValue);
}

/**
 * Constructs a matrix from the values of the given list. The values are
 * copied. So subsequent changes in <tt>values</tt> are not reflected in the
 * matrix, and vice-versa.
 *
 * @param values
 *            The values to be filled into the new matrix.
 * @return a new matrix.
 */
AbstractIntVector fromList(List<int> values, AbstractIntVector make(int)) {
  int size = values.length;
  AbstractIntVector vector = make(size);
  for (int i = size; --i >= 0; ) {
    vector[i] = values[i];
  }
  return vector;
}

/**
 * Constructs a matrix with uniformly distributed values in <tt>(0,1)</tt>
 * (exclusive).
 */
AbstractIntVector random(int size, AbstractIntVector make(int)) {
  return make(size)..forEach(ifunc.random());
}

/**
 * C = A||A||..||A; Constructs a new matrix which is concatenated
 * <tt>repeat</tt> times. Example:
 *
 * <pre>
 *   0 1
 *   repeat(3) --&gt;
 *   0 1 0 1 0 1
 *
 * </pre>
 */
AbstractIntVector repeat(AbstractIntVector A, int repeat, AbstractIntVector make(int)) {
  int size = A.length;
  AbstractIntVector matrix = make(repeat * size);
  for (int i = repeat; --i >= 0; ) {
    matrix.part(size * i, size).copyFrom(A);
  }
  return matrix;
}

/**
 * Constructs a randomly sampled matrix with the given shape. Randomly picks
 * exactly <tt>Math.round(size*nonZeroFraction)</tt> cells and initializes
 * them to <tt>value</tt>, all the rest will be initialized to zero. Note
 * that this is not the same as setting each cell with probability
 * <tt>nonZeroFraction</tt> to <tt>value</tt>.
 *
 * @throws ArgumentError
 *             if <tt>nonZeroFraction < 0 || nonZeroFraction > 1</tt>.
 * @see cern.jet.random.tdouble.sampling.DoubleRandomSamplingAssistant
 */
/*IntVector sample(int size, int value, int nonZeroFraction) {
  double epsilon = 1e-09;
  if (nonZeroFraction < 0 - epsilon || nonZeroFraction > 1 + epsilon) {
    throw new ArgumentError();
  }
  if (nonZeroFraction < 0) {
    nonZeroFraction = 0;
  }
  if (nonZeroFraction > 1) {
    nonZeroFraction = 1;
  }

  IntVector matrix = make(size);

  int n = (size * nonZeroFraction).round();
  if (n == 0) return matrix;

  final sampler = new DoubleRandomSamplingAssistant(n, size, new DoubleMersenneTwister());
  for (int i = size; --i >= 0; ) {
    if (sampler.sampleNextElement()) {
      matrix.put(i, value);
    }
  }

  return matrix;
}*/

/**
 * Constructs a list from the given matrix. The values are copied. So
 * subsequent changes in <tt>values</tt> are not reflected in the list, and
 * vice-versa.
 *
 * @param values
 *            The values to be filled into the new list.
 * @return a new list.
 */
Int32List toList(AbstractIntVector values) {
  int size = values.length;
  final list = new Int32List(size);
  //list.setSize(size);
  for (int i = size; --i >= 0; ) {
    list[i] = values[i];
  }
  return list;
}
//}

/**
 * Factory for convenient construction of 2-d matrices holding <tt>int</tt>
 * cells. Also provides convenient methods to compose (concatenate) and
 * decompose (split) matrices from/to constituent blocks. </p>
 * <p>
 * &nbsp;
 * </p>
 * <table border="0" cellspacing="0">
 * <tr align="left" valign="top">
 * <td><i>Construction</i></td>
 * <td>Use idioms like <tt>IntFactory2D.dense.make(4,4)</tt> to construct dense
 * matrices, <tt>IntFactory2D.sparse.make(4,4)</tt> to construct sparse
 * matrices.</td>
 * </tr>
 * <tr align="left" valign="top">
 * <td><i> Construction with initial values </i></td>
 * <td>Use other <tt>make</tt> methods to construct matrices with given initial
 * values.</td>
 * </tr>
 * <tr align="left" valign="top">
 * <td><i> Appending rows and columns </i></td>
 * <td>Use methods {@link #appendColumns(IntMatrix,IntMatrix) appendColumns}, {@link #appendColumns(IntMatrix,IntMatrix) appendRows} and
 * {@link #repeat(IntMatrix,int,int) repeat} to append rows and columns.</td>
 * </tr>
 * <tr align="left" valign="top">
 * <td><i> General block matrices </i></td>
 * <td>Use methods {@link #compose(IntMatrix[][]) compose} and
 * {@link #decompose(IntMatrix[][],IntMatrix) decompose} to work with
 * general block matrices.</td>
 * </tr>
 * <tr align="left" valign="top">
 * <td><i> Diagonal matrices </i></td>
 * <td>Use methods {@link #diagonal(IntVector) diagonal(vector)},
 * {@link #diagonal2D(IntMatrix) diagonal(matrix)} and {@link #identity(int)
 * identity} to work with diagonal matrices.</td>
 * </tr>
 * <tr align="left" valign="top">
 * <td><i> Diagonal block matrices </i></td>
 * <td>Use method {@link #composeDiagonal(IntMatrix,IntMatrix,IntMatrix)
 * composeDiagonal} to work with diagonal block matrices.</td>
 * </tr>
 * <tr align="left" valign="top">
 * <td><i>Random</i></td>
 * <td>Use methods {@link #random(int,int) random} and
 * {@link #sample(int,int,int,int) sample} to construct random matrices.</td>
 * </tr>
 * </table>
 * <p>
 * &nbsp;
 * </p>
 * <p>
 * If the factory is used frequently it might be useful to streamline the
 * notation. For example by aliasing:
 * </p>
 * <table>
 * <td class="PRE">
 *
 * <pre>
 *  IntFactory2D F = IntFactory2D.dense;
 *  F.make(4,4);
 *  F.descending(10,20);
 *  F.random(4,4);
 *  ...
 * </pre>
 *
 * </td>
 * </table>
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 09/24/99
 */
//class IntFactory2D {
//
//  /**
//   * A factory producing dense matrices.
//   */
//  static final IntFactory2D dense = new IntFactory2D._internal();
//
//  /**
//   * A factory producing sparse hash matrices.
//   */
//  static final IntFactory2D sparse = new IntFactory2D._internal();
//
//  /**
//   * A factory producing sparse row compressed matrices.
//   */
//  static final IntFactory2D rowCompressed = new IntFactory2D._internal();
//
//  /*
//   * A factory producing sparse row compressed modified matrices.
//   */
//  // static final IntFactory2D rowCompressedModified = new
//  // IntFactory2D();
//
//  IntFactory2D._internal();

/**
 * C = A||B; Constructs a new matrix which is the column-wise concatenation
 * of two other matrices.
 *
 * <pre>
 *   0 1 2
 *   3 4 5
 *   appendColumns
 *   6 7
 *   8 9
 *   --&gt;
 *   0 1 2 6 7
 *   3 4 5 8 9
 *
 * </pre>
 */
AbstractIntMatrix appendColumns(AbstractIntMatrix A, AbstractIntMatrix B, AbstractIntMatrix make(int, int)) {
  // force both to have maximal shared number of rows.
  if (B.rows > A.rows) {
    B = B.part(0, 0, A.rows, B.columns);
  } else if (B.rows < A.rows) {
    A = A.part(0, 0, B.rows, A.columns);
  }

  // concatenate
  int ac = A.columns;
  int bc = B.columns;
  int r = A.rows;
  AbstractIntMatrix matrix = make(r, ac + bc);
  matrix.part(0, 0, r, ac).copyFrom(A);
  matrix.part(0, ac, r, bc).copyFrom(B);
  return matrix;
}

AbstractIntMatrix appendColumn(AbstractIntMatrix A, AbstractIntVector b, AbstractIntMatrix make(int, int)) {
  // force both to have maximal shared number of rows.
  if (b.length > A.rows) {
    b = b.part(0, A.rows);
  } else if (b.length < A.rows) {
    A = A.part(0, 0, b.length, A.columns);
  }

  // concatenate
  int ac = A.columns;
  int bc = 1;
  int r = A.rows;
  AbstractIntMatrix matrix = make(r, ac + bc);
  matrix.part(0, 0, r, ac).copyFrom(A);
  matrix.column(ac).copyFrom(b);
  return matrix;
}

/**
 * C = A||B; Constructs a new matrix which is the row-wise concatenation of
 * two other matrices.
 *
 * <pre>
 *   0 1
 *   2 3
 *   4 5
 *   appendRows
 *   6 7
 *   8 9
 *   --&gt;
 *   0 1
 *   2 3
 *   4 5
 *   6 7
 *   8 9
 *
 * </pre>
 */
AbstractIntMatrix appendRows(AbstractIntMatrix A, AbstractIntMatrix B, AbstractIntMatrix make(int, int)) {
  // force both to have maximal shared number of columns.
  if (B.columns > A.columns) B = B.part(0, 0, B.rows, A.columns); else if (B.columns < A.columns) A = A.part(0, 0, A.rows, B.columns);

  // concatenate
  int ar = A.rows;
  int br = B.rows;
  int c = A.columns;
  AbstractIntMatrix matrix = make(ar + br, c);
  matrix.part(0, 0, ar, c).copyFrom(A);
  matrix.part(ar, 0, br, c).copyFrom(B);
  return matrix;
}

AbstractIntMatrix appendRow(AbstractIntMatrix A, AbstractIntVector b, AbstractIntMatrix make(int, int)) {
  // force both to have maximal shared number of columns.
  if (b.length > A.columns) b = b.part(0, A.columns); else if (b.length < A.columns) A = A.part(0, 0, A.rows, b.length);

  // concatenate
  int ar = A.rows;
  int br = 1;
  int c = A.columns;
  AbstractIntMatrix matrix = make(ar + br, c);
  matrix.part(0, 0, ar, c).copyFrom(A);
  matrix.row(ar).copyFrom(b);
  return matrix;
}

/**
 * Constructs a matrix with cells having ascending values. For debugging
 * purposes. Example:
 *
 * <pre>
 *   0 1 2
 *   3 4 5
 *
 * </pre>
 */
AbstractIntMatrix ascendingMatrix(int rows, int columns, AbstractIntMatrix make(int, int)) {
  return descendingMatrix(rows, columns, make)..forEach(ifunc.chain2(ifunc.neg, ifunc.subtract(columns * rows)));
}

/**
 * Checks whether the given array is rectangular, that is, whether all rows
 * have the same number of columns.
 *
 * @throws ArgumentError
 *             if the array is not rectangular.
 */
void _checkRectangularShape(List<Int32List> array) {
  int columns = -1;
  for (int row = array.length; --row >= 0; ) {
    if (array[row] != null) {
      if (columns == -1) columns = array[row].length;
      if (array[row].length != columns) throw new ArgumentError("All rows of array must have same number of columns.");
    }
  }
}

/**
 * Checks whether the given array is rectangular, that is, whether all rows
 * have the same number of columns.
 *
 * @throws ArgumentError
 *             if the array is not rectangular.
 */
void _checkRectShape(List<List<AbstractIntMatrix>> array) {
  int columns = -1;
  for (int row = array.length; --row >= 0; ) {
    if (array[row] != null) {
      if (columns == -1) {
        columns = array[row].length;
      }
      if (array[row].length != columns) {
        throw new ArgumentError("All rows of array must have same number of columns.");
      }
    }
  }
}

AbstractIntMatrix reshape(AbstractIntVector a, int rows, int columns, AbstractIntMatrix make(int, int)) {
  if (a.length != rows * columns) {
    throw new ArgumentError("a.length != rows*columns");
  }
  AbstractIntMatrix A = make(rows, columns);
  for (int c = 0; c < columns; c++) {
    A.column(c).copyFrom(a.part(c * rows, rows));
  }
  return A;
}

/**
 * Constructs a block matrix made from the given parts. The inverse to
 * method {@link #decompose(IntMatrix[][], IntMatrix)}.
 * <p>
 * All matrices of a given column within <tt>parts</tt> must have the same
 * number of columns. All matrices of a given row within <tt>parts</tt> must
 * have the same number of rows. Otherwise an
 * <tt>ArgumentError</tt> is thrown. Note that <tt>null</tt>s
 * within <tt>parts[row,col]</tt> are an exception to this rule: they are
 * ignored. Cells are copied. Example:
 * <table border="1" cellspacing="0">
 * <tr align="left" valign="top">
 * <td><tt>Code</tt></td>
 * <td><tt>Result</tt></td>
 * </tr>
 * <tr align="left" valign="top">
 * <td>
 *
 * <pre>
 * IntMatrix[][] parts1 = { { null, make(2, 2, 1), null }, { make(4, 4, 2), null, make(4, 3, 3) },
 *         { null, make(2, 2, 4), null } };
 * System.out.println(compose(parts1));
 * </pre>
 *
 * </td>
 * <td><tt>8&nbsp;x&nbsp;9&nbsp;matrix<br>
   0&nbsp;0&nbsp;0&nbsp;0&nbsp;1&nbsp;1&nbsp;0&nbsp;0&nbsp;0<br>
   0&nbsp;0&nbsp;0&nbsp;0&nbsp;1&nbsp;1&nbsp;0&nbsp;0&nbsp;0<br>
   2&nbsp;2&nbsp;2&nbsp;2&nbsp;0&nbsp;0&nbsp;3&nbsp;3&nbsp;3<br>
   2&nbsp;2&nbsp;2&nbsp;2&nbsp;0&nbsp;0&nbsp;3&nbsp;3&nbsp;3<br>
   2&nbsp;2&nbsp;2&nbsp;2&nbsp;0&nbsp;0&nbsp;3&nbsp;3&nbsp;3<br>
   2&nbsp;2&nbsp;2&nbsp;2&nbsp;0&nbsp;0&nbsp;3&nbsp;3&nbsp;3<br>
   0&nbsp;0&nbsp;0&nbsp;0&nbsp;4&nbsp;4&nbsp;0&nbsp;0&nbsp;0<br>
   0&nbsp;0&nbsp;0&nbsp;0&nbsp;4&nbsp;4&nbsp;0&nbsp;0&nbsp;0</tt></td>
 * </tr>
 * <tr align="left" valign="top">
 * <td>
 *
 * <pre>
 * IntMatrix[][] parts3 = { { identity(3), null, }, { null, identity(3).viewColumnFlip() },
 *         { identity(3).viewRowFlip(), null } };
 * System.out.println(&quot;\n&quot; + make(parts3));
 * </pre>
 *
 * </td>
 * <td><tt>9&nbsp;x&nbsp;6&nbsp;matrix<br>
   1&nbsp;0&nbsp;0&nbsp;0&nbsp;0&nbsp;0<br>
   0&nbsp;1&nbsp;0&nbsp;0&nbsp;0&nbsp;0<br>
   0&nbsp;0&nbsp;1&nbsp;0&nbsp;0&nbsp;0<br>
   0&nbsp;0&nbsp;0&nbsp;0&nbsp;0&nbsp;1<br>
   0&nbsp;0&nbsp;0&nbsp;0&nbsp;1&nbsp;0<br>
   0&nbsp;0&nbsp;0&nbsp;1&nbsp;0&nbsp;0<br>
   0&nbsp;0&nbsp;1&nbsp;0&nbsp;0&nbsp;0<br>
   0&nbsp;1&nbsp;0&nbsp;0&nbsp;0&nbsp;0<br>
   1&nbsp;0&nbsp;0&nbsp;0&nbsp;0&nbsp;0 </tt></td>
 * </tr>
 * <tr align="left" valign="top">
 * <td>
 *
 * <pre>
 * IntMatrix A = ascending(2, 2);
 * IntMatrix B = descending(2, 2);
 * IntMatrix _ = null;
 *
 * IntMatrix[][] parts4 = { { A, _, A, _ }, { _, A, _, B } };
 * System.out.println(&quot;\n&quot; + make(parts4));
 * </pre>
 *
 * </td>
 * <td><tt>4&nbsp;x&nbsp;8&nbsp;matrix<br>
   1&nbsp;2&nbsp;0&nbsp;0&nbsp;1&nbsp;2&nbsp;0&nbsp;0<br>
   3&nbsp;4&nbsp;0&nbsp;0&nbsp;3&nbsp;4&nbsp;0&nbsp;0<br>
   0&nbsp;0&nbsp;1&nbsp;2&nbsp;0&nbsp;0&nbsp;3&nbsp;2<br>
   0&nbsp;0&nbsp;3&nbsp;4&nbsp;0&nbsp;0&nbsp;1&nbsp;0 </tt></td>
 * </tr>
 * <tr align="left" valign="top">
 * <td>
 *
 * <pre>
 * IntMatrix[][] parts2 = { { null, make(2, 2, 1), null }, { make(4, 4, 2), null, make(4, 3, 3) },
 *         { null, make(2, 3, 4), null } };
 * System.out.println(&quot;\n&quot; + Factory2D.make(parts2));
 * </pre>
 *
 * </td>
 * <td><tt>ArgumentError<br>
   A[0,1].columns != A[2,1].columns<br>
   (2 != 3)</tt></td>
 * </tr>
 * </table>
 *
 * @throws ArgumentError
 *             subject to the conditions outlined above.
 */
AbstractIntMatrix compose(List<List<AbstractIntMatrix>> parts, AbstractIntMatrix make(int, int)) {
  _checkRectShape(parts);
  int rows = parts.length;
  int columns = 0;
  if (parts.length > 0) columns = parts[0].length;
  AbstractIntMatrix empty = make(0, 0);

  if (rows == 0 || columns == 0) return empty;

  // determine maximum column width of each column
  Int32List maxWidths = new Int32List(columns);
  for (int column = columns; --column >= 0; ) {
    int maxWidth = 0;
    for (int row = rows; --row >= 0; ) {
      AbstractIntMatrix part = parts[row][column];
      if (part != null) {
        int width = part.columns;
        if (maxWidth > 0 && width > 0 && width != maxWidth) throw new ArgumentError("Different number of columns.");
        maxWidth = Math.max(maxWidth, width);
      }
    }
    maxWidths[column] = maxWidth;
  }

  // determine row height of each row
  Int32List maxHeights = new Int32List(rows);
  for (int row = rows; --row >= 0; ) {
    int maxHeight = 0;
    for (int column = columns; --column >= 0; ) {
      AbstractIntMatrix part = parts[row][column];
      if (part != null) {
        int height = part.rows;
        if (maxHeight > 0 && height > 0 && height != maxHeight) throw new ArgumentError("Different number of rows.");
        maxHeight = Math.max(maxHeight, height);
      }
    }
    maxHeights[row] = maxHeight;
  }

  // shape of result
  int resultRows = 0;
  for (int row = rows; --row >= 0; ) resultRows += maxHeights[row];
  int resultCols = 0;
  for (int column = columns; --column >= 0; ) resultCols += maxWidths[column];

  AbstractIntMatrix matrix = make(resultRows, resultCols);

  // copy
  int r = 0;
  for (int row = 0; row < rows; row++) {
    int c = 0;
    for (int column = 0; column < columns; column++) {
      AbstractIntMatrix part = parts[row][column];
      if (part != null) {
        matrix.part(r, c, part.rows, part.columns).copyFrom(part);
      }
      c += maxWidths[column];
    }
    r += maxHeights[row];
  }

  return matrix;
}

/**
 * Constructs a diagonal block matrix from the given parts (the <i>direct
 * sum</i> of two matrices). That is the concatenation
 *
 * <pre>
 *   A 0
 *   0 B
 *
 * </pre>
 *
 * (The direct sum has <tt>A.rows+B.rows</tt> rows and
 * <tt>A.columns+B.columns</tt> columns). Cells are copied.
 *
 * @return a new matrix which is the direct sum.
 */
AbstractIntMatrix composeDiagonal(AbstractIntMatrix A, AbstractIntMatrix B, AbstractIntMatrix make(int, int)) {
  int ar = A.rows;
  int ac = A.columns;
  int br = B.rows;
  int bc = B.columns;
  AbstractIntMatrix sum = make(ar + br, ac + bc);
  sum.part(0, 0, ar, ac).copyFrom(A);
  sum.part(ar, ac, br, bc).copyFrom(B);
  return sum;
}

/**
 * Constructs a diagonal block matrix from the given parts. The
 * concatenation has the form
 *
 * <pre>
 *   A 0 0
 *   0 B 0
 *   0 0 C
 *
 * </pre>
 *
 * from the given parts. Cells are copied.
 */
AbstractIntMatrix composeDiag(AbstractIntMatrix A, AbstractIntMatrix B, AbstractIntMatrix C, AbstractIntMatrix make(int, int)) {
  AbstractIntMatrix diag = make(A.rows + B.rows + C.rows, A.columns + B.columns + C.columns);
  diag.part(0, 0, A.rows, A.columns).copyFrom(A);
  diag.part(A.rows, A.columns, B.rows, B.columns).copyFrom(B);
  diag.part(A.rows + B.rows, A.columns + B.columns, C.rows, C.columns).copyFrom(C);
  return diag;
}

AbstractIntMatrix composeBidiagonal(AbstractIntMatrix A, AbstractIntMatrix B, AbstractIntMatrix make(int, int)) {
  int ar = A.rows;
  int ac = A.columns;
  int br = B.rows;
  int bc = B.columns;
  AbstractIntMatrix sum = make(ar + br - 1, ac + bc);
  sum.part(0, 0, ar, ac).copyFrom(A);
  sum.part(ar - 1, ac, br, bc).copyFrom(B);
  return sum;
}

/**
 * Splits a block matrix into its constituent blocks; Copies blocks of a
 * matrix into the given parts. The inverse to method
 * {@link #compose(IntMatrix[][])}.
 * <p>
 * All matrices of a given column within <tt>parts</tt> must have the same
 * number of columns. All matrices of a given row within <tt>parts</tt> must
 * have the same number of rows. Otherwise an
 * <tt>ArgumentError</tt> is thrown. Note that <tt>null</tt>s
 * within <tt>parts[row,col]</tt> are an exception to this rule: they are
 * ignored. Cells are copied. Example:
 * <table border="1" cellspacing="0">
 * <tr align="left" valign="top">
 * <td><tt>Code</tt></td>
 * <td><tt>matrix</tt></td>
 * <td><tt>--&gt; parts </tt></td>
 * </tr>
 * <tr align="left" valign="top">
 * <td>
 *
 * <pre>
 *   IntMatrix matrix = ... ;
 *   IntMatrix _ = null;
 *   IntMatrix A,B,C,D;
 *   A = make(2,2); B = make (4,4);
 *   C = make(4,3); D = make (2,2);
 *   IntMatrix[][] parts =
 *   {
 *      { _, A, _ },
 *      { B, _, C },
 *      { _, D, _ }
 *   };
 *   decompose(parts,matrix);
 *   System.out.println(&quot;\nA = &quot;+A);
 *   System.out.println(&quot;\nB = &quot;+B);
 *   System.out.println(&quot;\nC = &quot;+C);
 *   System.out.println(&quot;\nD = &quot;+D);
 *
 * </pre>
 *
 * </td>
 * <td><tt>8&nbsp;x&nbsp;9&nbsp;matrix<br>
   9&nbsp;9&nbsp;9&nbsp;9&nbsp;1&nbsp;1&nbsp;9&nbsp;9&nbsp;9<br>
   9&nbsp;9&nbsp;9&nbsp;9&nbsp;1&nbsp;1&nbsp;9&nbsp;9&nbsp;9<br>
   2&nbsp;2&nbsp;2&nbsp;2&nbsp;9&nbsp;9&nbsp;3&nbsp;3&nbsp;3<br>
   2&nbsp;2&nbsp;2&nbsp;2&nbsp;9&nbsp;9&nbsp;3&nbsp;3&nbsp;3<br>
   2&nbsp;2&nbsp;2&nbsp;2&nbsp;9&nbsp;9&nbsp;3&nbsp;3&nbsp;3<br>
   2&nbsp;2&nbsp;2&nbsp;2&nbsp;9&nbsp;9&nbsp;3&nbsp;3&nbsp;3<br>
   9&nbsp;9&nbsp;9&nbsp;9&nbsp;4&nbsp;4&nbsp;9&nbsp;9&nbsp;9<br>
   9&nbsp;9&nbsp;9&nbsp;9&nbsp;4&nbsp;4&nbsp;9&nbsp;9&nbsp;9</tt></td>
 * <td>
 * <p>
 * <tt>A = 2&nbsp;x&nbsp;2&nbsp;matrix<br>
   1&nbsp;1<br>
   1&nbsp;1</tt>
 * </p>
 * <p>
 * <tt>B = 4&nbsp;x&nbsp;4&nbsp;matrix<br>
   2&nbsp;2&nbsp;2&nbsp;2<br>
   2&nbsp;2&nbsp;2&nbsp;2<br>
   2&nbsp;2&nbsp;2&nbsp;2<br>
   2&nbsp;2&nbsp;2&nbsp;2</tt>
 * </p>
 * <p>
 * <tt>C = 4&nbsp;x&nbsp;3&nbsp;matrix<br>
   3&nbsp;3&nbsp;3<br>
   3&nbsp;3&nbsp;3<br>
   </tt><tt>3&nbsp;3&nbsp;3<br>
   </tt><tt>3&nbsp;3&nbsp;3</tt>
 * </p>
 * <p>
 * <tt>D = 2&nbsp;x&nbsp;2&nbsp;matrix<br>
   4&nbsp;4<br>
   4&nbsp;4</tt>
 * </p>
 * </td>
 * </tr>
 * </table>
 *
 * @throws ArgumentError
 *             subject to the conditions outlined above.
 */
void decompose(List<List<AbstractIntMatrix>> parts, AbstractIntMatrix matrix) {
  _checkRectShape(parts);
  int rows = parts.length;
  int columns = 0;
  if (parts.length > 0) columns = parts[0].length;
  if (rows == 0 || columns == 0) return;

  // determine maximum column width of each column
  Int32List maxWidths = new Int32List(columns);
  for (int column = columns; --column >= 0; ) {
    int maxWidth = 0;
    for (int row = rows; --row >= 0; ) {
      AbstractIntMatrix part = parts[row][column];
      if (part != null) {
        int width = part.columns;
        if (maxWidth > 0 && width > 0 && width != maxWidth) throw new ArgumentError("Different number of columns.");
        maxWidth = Math.max(maxWidth, width);
      }
    }
    maxWidths[column] = maxWidth;
  }

  // determine row height of each row
  Int32List maxHeights = new Int32List(rows);
  for (int row = rows; --row >= 0; ) {
    int maxHeight = 0;
    for (int column = columns; --column >= 0; ) {
      AbstractIntMatrix part = parts[row][column];
      if (part != null) {
        int height = part.rows;
        if (maxHeight > 0 && height > 0 && height != maxHeight) throw new ArgumentError("Different number of rows.");
        maxHeight = Math.max(maxHeight, height);
      }
    }
    maxHeights[row] = maxHeight;
  }

  // shape of result parts
  int resultRows = 0;
  for (int row = rows; --row >= 0; ) resultRows += maxHeights[row];
  int resultCols = 0;
  for (int column = columns; --column >= 0; ) resultCols += maxWidths[column];

  if (matrix.rows < resultRows || matrix.columns < resultCols) throw new ArgumentError("Parts larger than matrix.");

  // copy
  int r = 0;
  for (int row = 0; row < rows; row++) {
    int c = 0;
    for (int column = 0; column < columns; column++) {
      AbstractIntMatrix part = parts[row][column];
      if (part != null) {
        part.copyFrom(matrix.part(r, c, part.rows, part.columns));
      }
      c += maxWidths[column];
    }
    r += maxHeights[row];
  }

}

/**
 * Demonstrates usage of this class.
 */
//void demo1() {
//  print("\n\n");
//  List<List<AbstractIntMatrix>> parts1 = [[null, fill(2, 2, 1), null], [fill(4, 4, 2), null, fill(4, 3, 3)], [null, fill(2, 2, 4), null]];
//  print("\n" + compose(parts1).toString());
//  // print("\n"+Converting.toHTML(make(parts1).toString()));
//
//  /*
//   * // illegal 2 != 3 IntMatrix[][] parts2 = { { null, make(2,2,1),
//   * null }, { make(4,4,2), null, make(4,3,3) }, { null, make(2,3,4), null } };
//   * print("\n"+make(parts2));
//   */
//
//  List<List<AbstractIntMatrix>> parts3 = [[identity(3), null,], [null, identity(3).columnFlip()], [identity(3).rowFlip(), null]];
//  print("\n" + compose(parts3).toString());
//  // print("\n"+Converting.toHTML(make(parts3).toString()));
//
//  AbstractIntMatrix A = ascending(2, 2);
//  AbstractIntMatrix B = descending(2, 2);
//  AbstractIntMatrix _ = null;
//
//  List<List<AbstractIntMatrix>> parts4 = [[A, _, A, _], [_, A, _, B]];
//  print("\n" + compose(parts4).toString());
//  // print("\n"+Converting.toHTML(make(parts4).toString()));
//
//}

/**
 * Demonstrates usage of this class.
 */
//void demo2() {
//  print("\n\n");
//  AbstractIntMatrix matrix;
//  AbstractIntMatrix A, B, C, D, E, F, G;
//  AbstractIntMatrix _ = null;
//  A = fill(2, 2, 1);
//  B = fill(4, 4, 2);
//  C = fill(4, 3, 3);
//  D = fill(2, 2, 4);
//  List<List<AbstractIntMatrix>> parts1 = [[_, A, _], [B, _, C], [_, D, _]];
//  matrix = compose(parts1);
//  print("\n" + matrix.toString());
//
//  A.fill(9);
//  B.fill(9);
//  C.fill(9);
//  D.fill(9);
//  decompose(parts1, matrix);
//  print(A);
//  print(B);
//  print(C);
//  print(D);
//  // print("\n"+Converting.toHTML(make(parts1).toString()));
//
//  /*
//   * // illegal 2 != 3 IntMatrix[][] parts2 = { { null, make(2,2,1),
//   * null }, { make(4,4,2), null, make(4,3,3) }, { null, make(2,3,4), null } };
//   * print("\n"+Factory2D.make(parts2));
//   */
//
//  /*
//   * IntMatrix[][] parts3 = { { identity(3), null, }, { null,
//   * identity(3).viewColumnFlip() }, { identity(3).viewRowFlip(), null } };
//   * System.out.println("\n"+make(parts3));
//   * //System.out.println("\n"+cern.colt.matrixpattern.Converting.toHTML(make(parts3).toString()));
//   *
//   * IntMatrix A = ascending(2,2); IntMatrix B =
//   * descending(2,2); IntMatrix _ = null;
//   *
//   * IntMatrix[][] parts4 = { { A, _, A, _ }, { _, A, _, B } };
//   * System.out.println("\n"+make(parts4));
//   * //System.out.println("\n"+cern.colt.matrixpattern.Converting.toHTML(make(parts4).toString()));
//   */
//}

/**
 * Constructs a matrix with cells having descending values. For debugging
 * purposes. Example:
 *
 * <pre>
 *   5 4 3
 *   2 1 0
 *
 * </pre>
 */
AbstractIntMatrix descendingMatrix(int rows, int columns, AbstractIntMatrix make(int, int)) {
  AbstractIntMatrix matrix = make(rows, columns);
  int v = 0;
  for (int row = rows; --row >= 0; ) {
    for (int column = columns; --column >= 0; ) {
      matrix.set(row, column, v++);
    }
  }
  return matrix;
}

/**
 * Constructs a new diagonal matrix whose diagonal elements are the elements
 * of <tt>vector</tt>. Cells values are copied. The new matrix is not a
 * view. Example:
 *
 * <pre>
 *   5 4 3 --&gt;
 *   5 0 0
 *   0 4 0
 *   0 0 3
 *
 * </pre>
 *
 * @return a new matrix.
 */
AbstractIntMatrix diagonal(AbstractIntVector vector, AbstractIntMatrix make(int, int)) {
  int size = vector.length;
  AbstractIntMatrix diag = make(size, size);
  for (int i = size; --i >= 0; ) {
    diag.set(i, i, vector.get(i));
  }
  return diag;
}

/**
 * Constructs a new diagonal matrix whose diagonal elements are the elements
 * of <tt>vector</tt>. Cells values are copied. The new matrix is not a
 * view. Example:
 *
 * <pre>
 *   5 4 3 --&gt;
 *   5 0 0
 *   0 4 0
 *   0 0 3
 *
 * </pre>
 *
 * @return a new matrix.
 */
AbstractIntMatrix diagonalList(Int32List vector, AbstractIntMatrix make(int, int)) {
  int size = vector.length;
  AbstractIntMatrix diag = make(size, size);
  for (int i = 0; i < size; i++) {
    diag.set(i, i, vector[i]);
  }
  return diag;
}

/**
 * Constructs a new vector consisting of the diagonal elements of <tt>A</tt>
 * . Cells values are copied. The new vector is not a view. Example:
 *
 * <pre>
 *   5 0 0 9
 *   0 4 0 9
 *   0 0 3 9
 *   --&gt; 5 4 3
 *
 * </pre>
 *
 * @param A
 *            the matrix, need not be square.
 * @return a new vector.
 */
AbstractIntVector diagonal2D(AbstractIntMatrix A, AbstractIntMatrix make(int, int)) {
  int min = Math.min(A.rows, A.columns);
  AbstractIntVector diag = _make1D(min, make);
  for (int i = min; --i >= 0; ) {
    diag.set(i, A.get(i, i));
  }
  return diag;
}

/**
 * Constructs an identity matrix (having ones on the diagonal and zeros
 * elsewhere).
 */
AbstractIntMatrix identity(int rowsAndColumns, AbstractIntMatrix make(int, int)) {
  AbstractIntMatrix matrix = make(rowsAndColumns, rowsAndColumns);
  for (int i = rowsAndColumns; --i >= 0; ) {
    matrix.set(i, i, 1);
  }
  return matrix;
}

/**
 * Constructs a matrix with the given cell values. <tt>values</tt> is
 * required to have the form <tt>values[row][column]</tt> and have exactly
 * the same number of columns in every row.
 * <p>
 * The values are copied. So subsequent changes in <tt>values</tt> are not
 * reflected in the matrix, and vice-versa.
 *
 * @param values
 *            The values to be filled into the new matrix.
 * @throws ArgumentError
 *             if
 *             <tt>for any 1 &lt;= row &lt; values.length: values[row].length != values[row-1].length</tt>
 *             .
 */
//  AbstractIntMatrix withValues(List<Int32List> values) {
//    if (this == sparse) {
//      return new SparseIntMatrix.fromList(values);
//    } else {
//      return new IntMatrix.fromList(values);
//    }
//  }

/**
 * Construct a matrix from a one-dimensional column-major packed array, ala
 * Fortran. Has the form
 * <tt>matrix.get(row,column) == values[row + column*rows]</tt>. The values
 * are copied.
 *
 * @param values
 *            One-dimensional array of doubles, packed by columns (ala
 *            Fortran).
 * @param rows
 *            the number of rows.
 * @exception ArgumentError
 *                <tt>values.length</tt> must be a multiple of <tt>rows</tt>
 *                .
 */
AbstractIntMatrix columnMajor(Int32List values, int rows, AbstractIntMatrix make(int, int)) {
  int columns = (rows != 0 ? values.length / rows : 0);
  if (rows * columns != values.length) {
    throw new ArgumentError("Array length must be a multiple of m.");
  }

  AbstractIntMatrix matrix = make(rows, columns);
  for (int row = 0; row < rows; row++) {
    for (int column = 0; column < columns; column++) {
      matrix.set(row, column, values[row + column * rows]);
    }
  }
  return matrix;
}

/**
 * Constructs a matrix with the given shape, each cell initialized with
 * zero.
 */
//  AbstractIntMatrix make(int rows, int columns) {
//    if (this == sparse) {
//      return new SparseIntMatrix(rows, columns);
//    }
//    if (this == rowCompressed) {
//      return new SparseRCIntMatrix(rows, columns); // if (this==rowCompressedModified) return new
//    }
//    // RCMIntMatrix(rows,columns);
//    else return new IntMatrix(rows, columns);
//  }

/**
 * Constructs a matrix with the given shape, each cell initialized with the
 * given value.
 */
AbstractIntMatrix fillMatrix(int rows, int columns, int initialValue, AbstractIntMatrix make(int, int)) {
  if (initialValue == 0) return make(rows, columns);
  return make(rows, columns)..fill(initialValue);
}

/**
 * Constructs a 1d matrix of the right dynamic type.
 */
AbstractIntVector _make1D(int size, AbstractIntMatrix make(int, int)) {
  return make(0, 0).like1D(size);
}

/**
 * Constructs a matrix with uniformly distributed values in <tt>(0,1)</tt>
 * (exclusive).
 */
AbstractIntMatrix randomMatrix(int rows, int columns, AbstractIntMatrix make(int, int)) {
  return make(rows, columns)..forEach(ifunc.random());
}

/**
 * C = A||A||..||A; Constructs a new matrix which is duplicated both along
 * the row and column dimension. Example:
 *
 * <pre>
 *   0 1
 *   2 3
 *   repeat(2,3) --&gt;
 *   0 1 0 1 0 1
 *   2 3 2 3 2 3
 *   0 1 0 1 0 1
 *   2 3 2 3 2 3
 *
 * </pre>
 */
AbstractIntMatrix repeatMatrix(AbstractIntMatrix A, int rowRepeat, int columnRepeat, AbstractIntMatrix make(int, int)) {
  int r = A.rows;
  int c = A.columns;
  AbstractIntMatrix matrix = make(r * rowRepeat, c * columnRepeat);
  for (int i = rowRepeat; --i >= 0; ) {
    for (int j = columnRepeat; --j >= 0; ) {
      matrix.part(r * i, c * j, r, c).copyFrom(A);
    }
  }
  return matrix;
}

/**
 * Constructs a randomly sampled matrix with the given shape. Randomly picks
 * exactly <tt>Math.round(rows*columns*nonZeroFraction)</tt> cells and
 * initializes them to <tt>value</tt>, all the rest will be initialized to
 * zero. Note that this is not the same as setting each cell with
 * probability <tt>nonZeroFraction</tt> to <tt>value</tt>. Note: The random
 * seed is a constant.
 *
 * @throws ArgumentError
 *             if <tt>nonZeroFraction < 0 || nonZeroFraction > 1</tt>.
 * @see cern.jet.random.tdouble.sampling.DoubleRandomSamplingAssistant
 */
/*IntMatrix sample(int rows, int columns, int value, int nonZeroFraction) {
  IntMatrix matrix = make(rows, columns);
  sample(matrix, value, nonZeroFraction);
  return matrix;
}*/

/**
 * Modifies the given matrix to be a randomly sampled matrix. Randomly picks
 * exactly <tt>Math.round(rows*columns*nonZeroFraction)</tt> cells and
 * initializes them to <tt>value</tt>, all the rest will be initialized to
 * zero. Note that this is not the same as setting each cell with
 * probability <tt>nonZeroFraction</tt> to <tt>value</tt>. Note: The random
 * seed is a constant.
 *
 * @throws ArgumentError
 *             if <tt>nonZeroFraction < 0 || nonZeroFraction > 1</tt>.
 * @see cern.jet.random.tdouble.sampling.DoubleRandomSamplingAssistant
 */
/*IntMatrix sample(IntMatrix matrix, int value, int nonZeroFraction) {
  int rows = matrix.rows;
  int columns = matrix.columns;
  double epsilon = 1e-09;
  if (nonZeroFraction < 0 - epsilon || nonZeroFraction > 1 + epsilon) throw new ArgumentError();
  if (nonZeroFraction < 0) nonZeroFraction = 0;
  if (nonZeroFraction > 1) nonZeroFraction = 1;

  matrix.fill(0);

  int size = rows * columns;
  int n = Math.round(size * nonZeroFraction);
  if (n == 0) return matrix;

  final sampler = new DoubleRandomSamplingAssistant(n, size, new DoubleMersenneTwister());
  for (int i = 0; i < size; i++) {
    if (sampler.sampleNextElement()) {
      int row = (i / columns);
      int column = (i % columns);
      matrix.put(row, column, value);
    }
  }

  return matrix;
}*/
//}
