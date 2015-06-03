library cern.colt.matrix.double.factory;

import 'dart:math' as Math;
import 'dart:typed_data';
import 'package:colt/function/double.dart' as func;
import '../matrix.dart';

/**
 * Factory for convenient construction of 1-d matrices holding <tt>double</tt>
 * cells. Use idioms like <tt>DoubleFactory1D.dense.make(1000)</tt> to construct
 * dense matrices, <tt>DoubleFactory1D.sparse.make(1000)</tt> to construct
 * sparse matrices.
 *
 * If the factory is used frequently it might be useful to streamline the
 * notation. For example by aliasing:
 * <table>
 * <td class="PRE">
 *
 * <pre>
 *  DoubleFactory1D F = DoubleFactory1D.dense;
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
//class DoubleFactory1D {
//
//  /**
//   * A factory producing dense matrices.
//   */
//  static final DoubleFactory1D dense = new DoubleFactory1D._internal();
//
//  /**
//   * A factory producing sparse matrices.
//   */
//  static final DoubleFactory1D sparse = new DoubleFactory1D._internal();
//
//  DoubleFactory1D._internal();

/**
 * C = A||B; Constructs a new matrix which is the concatenation of two other
 * matrices. Example: <tt>0 1</tt> append <tt>3 4</tt> --> <tt>0 1 3 4</tt>.
 */
AbstractDoubleVector append(AbstractDoubleVector A, AbstractDoubleVector B, AbstractDoubleVector make(int)) {
  // concatenate
  AbstractDoubleVector matrix = make(A.length + B.length);
  matrix.part(0, A.length).copyFrom(A);
  matrix.part(A.length, B.length).copyFrom(B);
  return matrix;
}

/**
 * Constructs a matrix with cells having ascending values. For debugging
 * purposes. Example: <tt>0 1 2</tt>
 */
AbstractDoubleVector ascending(int size, AbstractDoubleVector make(int)) {
  return descending(size, make)..forEach(func.chainGH(func.neg, func.subtract(size)));
}

/**
 * Constructs a matrix with cells having descending values. For debugging
 * purposes. Example: <tt>2 1 0</tt>
 */
AbstractDoubleVector descending(int size, AbstractDoubleVector make(int)) {
  AbstractDoubleVector matrix = make(size);
  double v = 0.0;
  for (int i = size; --i >= 0; ) {
    matrix.set(i, v);
    v += 1;
  }
  return matrix;
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
  /*AbstractDoubleVector makeList(Float64List values) {
    int size = values.length;
    AbstractDoubleVector vector = make(size);
    for (int i = size; --i >= 0; ) {
      vector.set(i, values[i]);
    }
    return vector;
  }*/

  /**
   * Constructs a matrix with the given cell values. The values are copied. So
   * subsequent changes in <tt>values</tt> are not reflected in the matrix,
   * and vice-versa.
   *
   * @param values
   *            The values to be filled into the new matrix.
   */
//  AbstractDoubleVector fromList(Float64List values) {
//    if (this == sparse) {
//      return new SparseDoubleVector.fromList(values);
//    } else {
//      return new DoubleVector.fromList(values);
//    }
//  }

/**
 * Constructs a matrix which is the concatenation of all given parts. Cells
 * are copied.
 */
AbstractDoubleVector concat(List<AbstractDoubleVector> parts, AbstractDoubleVector make(int)) {
  if (parts.length == 0) return make(0);

  int size = 0;
  for (int i = 0; i < parts.length; i++) size += parts[i].length;

  AbstractDoubleVector vector = make(size);
  size = 0;
  for (int i = 0; i < parts.length; i++) {
    vector.part(size, parts[i].length).copyFrom(parts[i]);
    size += parts[i].length;
  }

  return vector;
}

//  /**
//   * Constructs a matrix with the given shape, each cell initialized with
//   * zero.
//   */
//  AbstractDoubleVector make(int size) {
//    if (this == sparse) {
//      return new SparseDoubleVector(size);
//    }
//    return new DoubleVector(size);
//  }

/**
 * Constructs a matrix with the given shape, each cell initialized with the
 * given value.
 */
AbstractDoubleVector fill(int size, double initialValue, AbstractDoubleVector make(int)) {
  return make(size)..fill(initialValue);
}

/**
 * Constructs a matrix with uniformly distributed values in <tt>(0,1)</tt>
 * (exclusive).
 */
AbstractDoubleVector random(int size, AbstractDoubleVector make(int)) {
  return make(size)..forEach(func.random());
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
AbstractDoubleVector repeat(AbstractDoubleVector A, int repeat, AbstractDoubleVector make(int)) {
  int size = A.length;
  AbstractDoubleVector matrix = make(repeat * size);
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
   * @throws IllegalArgumentException
   *             if <tt>nonZeroFraction < 0 || nonZeroFraction > 1</tt>.
   * @see cern.jet.random.tdouble.sampling.DoubleRandomSampler
   */
  /*DoubleVector sample(int size, double value, double nonZeroFraction) {
    double epsilon = 1e-09;
    if (nonZeroFraction < 0 - epsilon || nonZeroFraction > 1 + epsilon) {
      throw new ArgumentError();
    }
    if (nonZeroFraction < 0) nonZeroFraction = 0.0;
    if (nonZeroFraction > 1) nonZeroFraction = 1.0;

    DoubleVector matrix = make(size);

    int n = Math.round(size * nonZeroFraction);
    if (n == 0) return matrix;

    DoubleRandomSamplingAssistant sampler = new DoubleRandomSamplingAssistant(n, size, new DoubleMersenneTwister());
    for (int i = size; --i >= 0; ) {
      if (sampler.sampleNextElement()) {
        matrix.setQuick(i, value);
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
  /*DoubleArrayList toList(DoubleVector values) {
    int size = values.size();
    DoubleArrayList list = new DoubleArrayList(size);
    list.setSize(size);
    for (int i = size; --i >= 0; ) list._setQuick(i, values.get(i));
    return list;
  }*/
//}


/**
 * Factory for convenient construction of 2-d matrices holding <tt>double</tt>
 * cells. Also provides convenient methods to compose (concatenate) and
 * decompose (split) matrices from/to constituent blocks. </p>
 * <p>
 * &nbsp;
 * </p>
 * <table border="0" cellspacing="0">
 * <tr align="left" valign="top">
 * <td><i>Construction</i></td>
 * <td>Use idioms like <tt>DoubleFactory2D.dense.make(4,4)</tt> to construct
 * dense matrices, <tt>DoubleFactory2D.sparse.make(4,4)</tt> to construct sparse
 * matrices.</td>
 * </tr>
 * <tr align="left" valign="top">
 * <td><i> Construction with initial values </i></td>
 * <td>Use other <tt>make</tt> methods to construct matrices with given initial
 * values.</td>
 * </tr>
 * <tr align="left" valign="top">
 * <td><i> Appending rows and columns </i></td>
 * <td>Use methods {@link #appendColumns(DoubleMatrix,DoubleMatrix)
 * appendColumns}, {@link #appendColumns(DoubleMatrix,DoubleMatrix)
 * appendRows} and {@link #repeat(DoubleMatrix,int,int) repeat} to append rows
 * and columns.</td>
 * </tr>
 * <tr align="left" valign="top">
 * <td><i> General block matrices </i></td>
 * <td>Use methods {@link #compose(DoubleMatrix[][]) compose} and
 * {@link #decompose(DoubleMatrix[][],DoubleMatrix) decompose} to work with
 * general block matrices.</td>
 * </tr>
 * <tr align="left" valign="top">
 * <td><i> Diagonal matrices </i></td>
 * <td>Use methods {@link #diagonal(DoubleVector) diagonal(vector)},
 * {@link #diagonal(DoubleMatrix) diagonal(matrix)} and {@link #identity(int)
 * identity} to work with diagonal matrices.</td>
 * </tr>
 * <tr align="left" valign="top">
 * <td><i> Diagonal block matrices </i></td>
 * <td>Use method
 * {@link #composeDiagonal(DoubleMatrix,DoubleMatrix,DoubleMatrix)
 * composeDiagonal} to work with diagonal block matrices.</td>
 * </tr>
 * <tr align="left" valign="top">
 * <td><i>Random</i></td>
 * <td>Use methods {@link #random(int,int) random} and
 * {@link #sample(int,int,double,double) sample} to construct random matrices.</td>
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
 *  DoubleFactory2D F = DoubleFactory2D.dense;
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
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
//class DoubleFactory2D {
//
//  /**
//   * A factory producing dense matrices.
//   */
//  static final DoubleFactory2D dense = new DoubleFactory2D._internal();
//
//  /**
//   * A factory producing sparse hash matrices.
//   */
//  static final DoubleFactory2D sparse = new DoubleFactory2D._internal();

/**
 * Checks whether the given array is rectangular, that is, whether all rows
 * have the same number of columns.
 *
 * @throws IllegalArgumentException
 *             if the array is not rectangular.
 */
void _checkRectangularShape(List<Float64List> array) {
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

/**
 * Checks whether the given array is rectangular, that is, whether all rows
 * have the same number of columns.
 *
 * @throws IllegalArgumentException
 *             if the array is not rectangular.
 */
void _checkRectangularShapeParts(List<List<AbstractDoubleMatrix>> array) {
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

//  DoubleFactory2D._internal();

/**
 * C = A||b; Constructs a new matrix which is the column-wise concatenation
 * of two other matrices.
 *
 * <pre>
 *   0 1 2
 *   3 4 5
 *   appendColumn
 *   6
 *   8
 *   --&gt;
 *   0 1 2 6
 *   3 4 5 8
 *
 * </pre>
 */
AbstractDoubleMatrix appendColumn(AbstractDoubleMatrix A, AbstractDoubleVector b, AbstractDoubleMatrix make(int r, int c)) {
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
  AbstractDoubleMatrix matrix = make(r, ac + bc);
  matrix.part(0, 0, r, ac).copyFrom(A);
  matrix.column(ac).copyFrom(b);
  return matrix;
}

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
AbstractDoubleMatrix appendColumns(AbstractDoubleMatrix A, AbstractDoubleMatrix B, AbstractDoubleMatrix make(int r, int c)) {
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
  AbstractDoubleMatrix matrix = make(r, ac + bc);
  matrix.part(0, 0, r, ac).copyFrom(A);
  matrix.part(0, ac, r, bc).copyFrom(B);
  return matrix;
}

/**
 * C = A||b; Constructs a new matrix which is the row-wise concatenation of
 * two other matrices.
 *
 * <pre>
 *   0 1
 *   2 3
 *   4 5
 *   appendRow
 *   6 7
 *   --&gt;
 *   0 1
 *   2 3
 *   4 5
 *   6 7
 *
 * </pre>
 */
AbstractDoubleMatrix appendRow(AbstractDoubleMatrix A, AbstractDoubleVector b, AbstractDoubleMatrix make(int r, int c)) {
  // force both to have maximal shared number of columns.
  if (b.length > A.columns) {
    b = b.part(0, A.columns);
  } else if (b.length < A.columns) {
    A = A.part(0, 0, A.rows, b.length);
  }

  // concatenate
  int ar = A.rows;
  int br = 1;
  int c = A.columns;
  AbstractDoubleMatrix matrix = make(ar + br, c);
  matrix.part(0, 0, ar, c).copyFrom(A);
  matrix.row(ar).copyFrom(b);
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
AbstractDoubleMatrix appendRows(AbstractDoubleMatrix A, AbstractDoubleMatrix B, AbstractDoubleMatrix make(int r, int c)) {
  // force both to have maximal shared number of columns.
  if (B.columns > A.columns) {
    B = B.part(0, 0, B.rows, A.columns);
  } else if (B.columns < A.columns) {
    A = A.part(0, 0, A.rows, B.columns);
  }

  // concatenate
  int ar = A.rows;
  int br = B.rows;
  int c = A.columns;
  AbstractDoubleMatrix matrix = make(ar + br, c);
  matrix.part(0, 0, ar, c).copyFrom(A);
  matrix.part(ar, 0, br, c).copyFrom(B);
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
AbstractDoubleMatrix ascendingMatrix(int rows, int columns, AbstractDoubleMatrix make(int r, int c)) {
  return descendingMatrix(rows, columns, make)..forEach(func.chainGH(func.neg, func.subtract(columns * rows)));
}

/**
 * Constructs a block matrix made from the given parts. The inverse to
 * method {@link #decompose(DoubleMatrix[][], DoubleMatrix)}.
 * <p>
 * All matrices of a given column within <tt>parts</tt> must have the same
 * number of columns. All matrices of a given row within <tt>parts</tt> must
 * have the same number of rows. Otherwise an
 * <tt>IllegalArgumentException</tt> is thrown. Note that <tt>null</tt>s
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
 * DoubleMatrix[][] parts1 = { { null, make(2, 2, 1), null }, { make(4, 4, 2), null, make(4, 3, 3) },
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
 * DoubleMatrix[][] parts3 = { { identity(3), null, }, { null, identity(3).viewColumnFlip() },
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
 * DoubleMatrix A = ascending(2, 2);
 * DoubleMatrix B = descending(2, 2);
 * DoubleMatrix _ = null;
 *
 * DoubleMatrix[][] parts4 = { { A, _, A, _ }, { _, A, _, B } };
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
 * DoubleMatrix[][] parts2 = { { null, make(2, 2, 1), null }, { make(4, 4, 2), null, make(4, 3, 3) },
 *         { null, make(2, 3, 4), null } };
 * System.out.println(&quot;\n&quot; + Factory2D.make(parts2));
 * </pre>
 *
 * </td>
 * <td><tt>IllegalArgumentException<br>
   A[0,1].columns != A[2,1].columns<br>
   (2 != 3)</tt></td>
 * </tr>
 * </table>
 *
 * @throws IllegalArgumentException
 *             subject to the conditions outlined above.
 */
AbstractDoubleMatrix compose(List<List<AbstractDoubleMatrix>> parts, AbstractDoubleMatrix make(int r, int c)) {
  _checkRectangularShapeParts(parts);
  int rows = parts.length;
  int columns = 0;
  if (parts.length > 0) columns = parts[0].length;
  AbstractDoubleMatrix empty = make(0, 0);

  if (rows == 0 || columns == 0) return empty;

  // determine maximum column width of each column
  Int32List maxWidths = new Int32List(columns);
  for (int column = columns; --column >= 0; ) {
    int maxWidth = 0;
    for (int row = rows; --row >= 0; ) {
      AbstractDoubleMatrix part = parts[row][column];
      if (part != null) {
        int width = part.columns;
        if (maxWidth > 0 && width > 0 && width != maxWidth) {
          throw new ArgumentError("Different number of columns.");
        }
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
      AbstractDoubleMatrix part = parts[row][column];
      if (part != null) {
        int height = part.rows;
        if (maxHeight > 0 && height > 0 && height != maxHeight) {
          throw new ArgumentError("Different number of rows.");
        }
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

  AbstractDoubleMatrix matrix = make(resultRows, resultCols);

  // copy
  int r = 0;
  for (int row = 0; row < rows; row++) {
    int c = 0;
    for (int column = 0; column < columns; column++) {
      AbstractDoubleMatrix part = parts[row][column];
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
 * Constructs a bidiagonal block matrix from the given parts. The
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
 *
 * @param A
 *            bidiagonal matrix
 * @param B
 *            bidiagonal matrix
 * @return bidiagonal matrix
 */
AbstractDoubleMatrix composeBidiagonal(AbstractDoubleMatrix A, AbstractDoubleMatrix B, AbstractDoubleMatrix make(int r, int c)) {
  int ar = A.rows;
  int ac = A.columns;
  int br = B.rows;
  int bc = B.columns;
  AbstractDoubleMatrix sum = make(ar + br - 1, ac + bc);
  sum.part(0, 0, ar, ac).copyFrom(A);
  sum.part(ar - 1, ac, br, bc).copyFrom(B);
  return sum;
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
 * (The direct sum has <tt>A.rows()+B.rows()</tt> rows and
 * <tt>A.columns()+B.columns()</tt> columns). Cells are copied.
 *
 * @return a new matrix which is the direct sum.
 */
AbstractDoubleMatrix composeDiagonal(AbstractDoubleMatrix A, AbstractDoubleMatrix B, AbstractDoubleMatrix make(int r, int c)) {
  int ar = A.rows;
  int ac = A.columns;
  int br = B.rows;
  int bc = B.columns;
  AbstractDoubleMatrix sum = make(ar + br, ac + bc);
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
AbstractDoubleMatrix composeDiagonal3(AbstractDoubleMatrix A, AbstractDoubleMatrix B, AbstractDoubleMatrix C, AbstractDoubleMatrix make(int r, int c)) {
  AbstractDoubleMatrix diag = make(A.rows + B.rows + C.rows, A.columns + B.columns + C.columns);
  diag.part(0, 0, A.rows, A.columns).copyFrom(A);
  diag.part(A.rows, A.columns, B.rows, B.columns).copyFrom(B);
  diag.part(A.rows + B.rows, A.columns + B.columns, C.rows, C.columns).copyFrom(C);
  return diag;
}

/**
 * Splits a block matrix into its constituent blocks; Copies blocks of a
 * matrix into the given parts. The inverse to method
 * {@link #compose(DoubleMatrix[][])}.
 * <p>
 * All matrices of a given column within <tt>parts</tt> must have the same
 * number of columns. All matrices of a given row within <tt>parts</tt> must
 * have the same number of rows. Otherwise an
 * <tt>IllegalArgumentException</tt> is thrown. Note that <tt>null</tt>s
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
 *   DoubleMatrix matrix = ... ;
 *   DoubleMatrix _ = null;
 *   DoubleMatrix A,B,C,D;
 *   A = make(2,2); B = make (4,4);
 *   C = make(4,3); D = make (2,2);
 *   DoubleMatrix[][] parts =
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
 * @throws IllegalArgumentException
 *             subject to the conditions outlined above.
 */
void decompose(List<List<AbstractDoubleMatrix>> parts, AbstractDoubleMatrix matrix) {
  _checkRectangularShapeParts(parts);
  int rows = parts.length;
  int columns = 0;
  if (parts.length > 0) columns = parts[0].length;
  if (rows == 0 || columns == 0) return;

  // determine maximum column width of each column
  Int32List maxWidths = new Int32List(columns);
  for (int column = columns; --column >= 0; ) {
    int maxWidth = 0;
    for (int row = rows; --row >= 0; ) {
      AbstractDoubleMatrix part = parts[row][column];
      if (part != null) {
        int width = part.columns;
        if (maxWidth > 0 && width > 0 && width != maxWidth) {
          throw new ArgumentError("Different number of columns.");
        }
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
      AbstractDoubleMatrix part = parts[row][column];
      if (part != null) {
        int height = part.rows;
        if (maxHeight > 0 && height > 0 && height != maxHeight) {
          throw new ArgumentError("Different number of rows.");
        }
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

  if (matrix.rows < resultRows || matrix.columns < resultCols) {
    throw new ArgumentError("Parts larger than matrix.");
  }

  // copy
  int r = 0;
  for (int row = 0; row < rows; row++) {
    int c = 0;
    for (int column = 0; column < columns; column++) {
      AbstractDoubleMatrix part = parts[row][column];
      if (part != null) {
        part.copyFrom(matrix.part(r, c, part.rows, part.columns));
      }
      c += maxWidths[column];
    }
    r += maxHeights[row];
  }

}

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
AbstractDoubleMatrix descendingMatrix(int rows, int columns, AbstractDoubleMatrix make(int r, int c)) {
  AbstractDoubleMatrix matrix = make(rows, columns);
  double v = 0.0;
  for (int row = rows; --row >= 0; ) {
    for (int column = columns; --column >= 0; ) {
      matrix.set(row, column, v);
      v += 1;
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
AbstractDoubleMatrix diagonal(Float64List vector, AbstractDoubleMatrix make(int r, int c)) {
  int size = vector.length;
  AbstractDoubleMatrix diag = make(size, size);
  for (int i = 0; i < size; i++) {
    diag.set(i, i, vector[i]);
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
AbstractDoubleMatrix diagonal1D(AbstractDoubleVector vector, AbstractDoubleMatrix make(int r, int c)) {
  int size = vector.length;
  AbstractDoubleMatrix diag = make(size, size);
  for (int i = size; --i >= 0; ) {
    diag.set(i, i, vector.get(i));
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
AbstractDoubleVector diagonal2D(AbstractDoubleMatrix A, AbstractDoubleMatrix make(int r, int c)) {
  int min = Math.min(A.rows, A.columns);
  AbstractDoubleVector diag = _make1D(min, make);
  for (int i = min; --i >= 0; ) {
    diag.set(i, A.get(i, i));
  }
  return diag;
}

/**
 * Constructs an identity matrix (having ones on the diagonal and zeros
 * elsewhere).
 */
AbstractDoubleMatrix identity(int rowsAndColumns, AbstractDoubleMatrix make(int r, int c)) {
  AbstractDoubleMatrix matrix = make(rowsAndColumns, rowsAndColumns);
  for (int i = rowsAndColumns; --i >= 0; ) {
    matrix.set(i, i, 1.0);
  }
  return matrix;
}

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
 * @exception IllegalArgumentException
 *                <tt>values.length</tt> must be a multiple of <tt>rows</tt>
 *                .
 */
AbstractDoubleMatrix makeColumn(Float64List values, int rows, AbstractDoubleMatrix make(int r, int c)) {
  int columns = (rows != 0 ? values.length / rows : 0);
  if (rows * columns != values.length) {
    throw new ArgumentError("Array length must be a multiple of m.");
  }

  AbstractDoubleMatrix matrix = make(rows, columns);
  for (int row = 0; row < rows; row++) {
    for (int column = 0; column < columns; column++) {
      matrix.set(row, column, values[row + column * rows]);
    }
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
   * @throws IllegalArgumentException
   *             if
   *             <tt>for any 1 &lt;= row &lt; values.length: values[row].length != values[row-1].length</tt>
   *             .
   */
//  AbstractDoubleMatrix fromList(List<Float64List> values) {
//    if (this == sparse) {
//      return new SparseDoubleMatrix.fromList(values);
//    } else {
//      return new DoubleMatrix.fromList(values);
//    }
//  }

  /**
   * Constructs a matrix with the given shape, each cell initialized with
   * zero.
   */
//  AbstractDoubleMatrix make(int rows, int columns) {
//    if (this == sparse) {
//      return new SparseDoubleMatrix(rows, columns);
//    } else {
//      return new DoubleMatrix(rows, columns);
//    }
//  }

/**
 * Constructs a matrix with the given shape, each cell initialized with the
 * given value.
 */
AbstractDoubleMatrix fillMatrix(int rows, int columns, double initialValue, AbstractDoubleMatrix make(int r, int c)) {
  if (initialValue == 0) return make(rows, columns);
  return make(rows, columns)..fill(initialValue);
}

/**
 * Constructs a matrix with uniformly distributed values in <tt>(0,1)</tt>
 * (exclusive).
 */
AbstractDoubleMatrix randomMatrix(int rows, int columns, AbstractDoubleMatrix make(int r, int c)) {
  return make(rows, columns)..forEach(func.random());
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
AbstractDoubleMatrix repeatMatrix(AbstractDoubleMatrix A, int rowRepeat, int columnRepeat, AbstractDoubleMatrix make(int r, int c)) {
  int r = A.rows;
  int c = A.columns;
  AbstractDoubleMatrix matrix = make(r * rowRepeat, c * columnRepeat);
  for (int i = rowRepeat; --i >= 0; ) {
    for (int j = columnRepeat; --j >= 0; ) {
      matrix.part(r * i, c * j, r, c).copyFrom(A);
    }
  }
  return matrix;
}

/**
 * Modifies the given matrix to be a randomly sampled matrix. Randomly picks
 * exactly <tt>Math.round(rows*columns*nonZeroFraction)</tt> cells and
 * initializes them to <tt>value</tt>, all the rest will be initialized to
 * zero. Note that this is not the same as setting each cell with
 * probability <tt>nonZeroFraction</tt> to <tt>value</tt>. Note: The random
 * seed is a constant.
 *
 * @throws IllegalArgumentException
 *             if <tt>nonZeroFraction < 0 || nonZeroFraction > 1</tt>.
 * @see cern.jet.random.tdouble.sampling.DoubleRandomSampler
 */
/*DoubleMatrix sample(DoubleMatrix matrix, double value, double nonZeroFraction) {
  int rows = matrix.rows();
  int columns = matrix.columns();
  double epsilon = 1e-09;
  if (nonZeroFraction < 0 - epsilon || nonZeroFraction > 1 + epsilon) throw new IllegalArgumentException();
  if (nonZeroFraction < 0) nonZeroFraction = 0;
  if (nonZeroFraction > 1) nonZeroFraction = 1;

  matrix.assign(0);

  int size = rows * columns;
  int n = Math.round(size * nonZeroFraction);
  if (n == 0) return matrix;

  DoubleRandomSamplingAssistant sampler = new DoubleRandomSamplingAssistant(n, size, new DoubleMersenneTwister());
  for (int i = 0; i < size; i++) {
    if (sampler.sampleNextElement()) {
      int row = (i / columns);
      int column = (i % columns);
      matrix.setQuick(row, column, value);
    }
  }

  return matrix;
}*/

/**
 * Constructs a randomly sampled matrix with the given shape. Randomly picks
 * exactly <tt>Math.round(rows*columns*nonZeroFraction)</tt> cells and
 * initializes them to <tt>value</tt>, all the rest will be initialized to
 * zero. Note that this is not the same as setting each cell with
 * probability <tt>nonZeroFraction</tt> to <tt>value</tt>. Note: The random
 * seed is a constant.
 *
 * @throws IllegalArgumentException
 *             if <tt>nonZeroFraction < 0 || nonZeroFraction > 1</tt>.
 * @see cern.jet.random.tdouble.sampling.DoubleRandomSampler
 */
/*DoubleMatrix sample(int rows, int columns, double value, double nonZeroFraction) {
  DoubleMatrix matrix = make(rows, columns);
  sample(matrix, value, nonZeroFraction);
  return matrix;
}*/

/**
 * Constructs a 1d matrix of the right dynamic type.
 */
AbstractDoubleVector _make1D(int size, AbstractDoubleMatrix make(int r, int c)) {
  return make(0, 0).like1D(size);
}
//}
