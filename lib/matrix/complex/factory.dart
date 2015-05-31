/*
Copyright (C) 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
is hereby granted without fee, provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear in supporting documentation.
CERN makes no representations about the suitability of this software for any purpose.
It is provided "as is" without expressed or implied warranty.
 */
library cern.colt.matrix.complex.factory;

import 'dart:math' as Math;
import 'dart:typed_data';
import 'package:colt/function/complex.dart' as cfunc;
import '../matrix.dart';

/**
 * Factory for convenient construction of 1-d matrices holding <tt>complex</tt>
 * cells. Use idioms like <tt>ComplexFactory1D.dense.make(1000)</tt> to
 * construct dense matrices, <tt>ComplexFactory1D.sparse.make(1000)</tt> to
 * construct sparse matrices.
 *
 * If the factory is used frequently it might be useful to streamline the
 * notation. For example by aliasing:
 * <table>
 * <td class="PRE">
 *
 * <pre>
 *  ComplexFactory1D F = ComplexFactory1D.dense;
 *  F.make(1000);
 *  F.random(3);
 *  ...
 * </pre>
 *
 * </td>
 * </table>
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
//class DComplexFactory1D {
//
//  /**
//   * A factory producing dense matrices.
//   */
//  static final DComplexFactory1D dense = new DComplexFactory1D._internal();
//
//  /**
//   * A factory producing sparse matrices.
//   */
//  static final DComplexFactory1D sparse = new DComplexFactory1D._internal();
//
//  /**
//   * Makes this class non instantiable, but still let's others inherit from
//   * it.
//   */
//  DComplexFactory1D._internal();

/**
 * C = A||B; Constructs a new matrix which is the concatenation of two other
 * matrices. Example: <tt>0 1</tt> append <tt>3 4</tt> --> <tt>0 1 3 4</tt>.
 */
AbstractComplexVector append(AbstractComplexVector A, AbstractComplexVector B, AbstractComplexVector make(int)) {
  // concatenate
  AbstractComplexVector matrix = make((A.length + B.length));
  matrix.part(0, A.length).copyFrom(A);
  matrix.part(A.length, B.length).copyFrom(B);
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
//  AbstractComplexVector withValues(List<double> values) {
//    if (this == sparse) {
//      return new SparseComplexVector.fromList(values);
//    } else {
//      return new ComplexVector.fromList(values);
//    }
//  }

/**
 * Constructs a matrix which is the concatenation of all given parts. Cells
 * are copied.
 */
AbstractComplexVector concat(List<AbstractComplexVector> parts, AbstractComplexVector make(int)) {
  if (parts.length == 0) {
    return make(0);
  }

  int size = 0;
  for (int i = 0; i < parts.length; i++) {
    size += parts[i].length;
  }

  AbstractComplexVector vector = make(size);
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
//  AbstractComplexVector make(int size) {
//    if (this == sparse) {
//      return new SparseComplexVector(size);
//    } else {
//      return new ComplexVector(size);
//    }
//  }

/**
 * Constructs a matrix with the given shape, each cell initialized with the
 * given value.
 */
AbstractComplexVector fill(int size, List<double> initialValue, AbstractComplexVector make(int)) {
  return make(size)..fill(initialValue[0], initialValue[1]);
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
AbstractComplexVector fromList(List<List<double>> values, AbstractComplexVector make(int)) {
  int size = values.length;
  AbstractComplexVector vector = make(size);
  for (int i = 0; i < size; i++) {
    vector.set(i, values[i]);
  }
  return vector;
}

/**
 * Constructs a matrix with uniformly distributed values in <tt>(0,1)</tt>
 * (exclusive).
 */
AbstractComplexVector random(int size, AbstractComplexVector make(int)) {
  return make(size)..forEach(cfunc.random);
}

/**
 * C = A||A||..||A; Constructs a new matrix which is concatenated
 * <tt>repeat</tt> times.
 */
AbstractComplexVector repeat(AbstractComplexVector A, int repeat, AbstractComplexVector make(int)) {
  int size = A.length;
  AbstractComplexVector matrix = make(repeat * size);
  for (int i = 0; i < repeat; i++) {
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
/*AbstractComplexVector sample(int size, Float64List value, double nonZeroFraction) {
  double epsilon = 1e-09;
  if (nonZeroFraction < 0 - epsilon || nonZeroFraction > 1 + epsilon) {
    throw new ArgumentError();
  }
  if (nonZeroFraction < 0) {
    nonZeroFraction = 0.0;
  }
  if (nonZeroFraction > 1) {
    nonZeroFraction = 1.0;
  }

  AbstractComplexVector matrix = make(size);

  int n = (size * nonZeroFraction).round();
  if (n == 0) {
    return matrix;
  }

  final sampler = new DoubleRandomSamplingAssistant(n, size, new DoubleMersenneTwister());
  for (int i = 0; i < size; i++) {
    if (sampler.sampleNextElement()) {
      matrix.set(i, value);
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
List<List<double>> toList(AbstractComplexVector values) {
  int size = values.length;
  List<Float64List> list = new List<Float64List>(size);
  for (int i = 0; i < size; i++) {
    list[i] = values.get(i);
  }
  return list;
}
//}

/**
 * Factory for convenient construction of 2-d matrices holding <tt>complex</tt>
 * cells. Also provides convenient methods to compose (concatenate) and
 * decompose (split) matrices from/to constituent blocks. </p>
 * <p>
 * &nbsp;
 * </p>
 * <table border="0" cellspacing="0">
 * <tr align="left" valign="top">
 * <td><i>Construction</i></td>
 * <td>Use idioms like <tt>ComplexFactory2D.dense.make(4,4)</tt> to construct
 * dense matrices, <tt>ComplexFactory2D.sparse.make(4,4)</tt> to construct
 * sparse matrices.</td>
 * </tr>
 * <tr align="left" valign="top">
 * <td><i> Construction with initial values </i></td>
 * <td>Use other <tt>make</tt> methods to construct matrices with given initial
 * values.</td>
 * </tr>
 * <tr align="left" valign="top">
 * <td><i> Appending rows and columns </i></td>
 * <td>Use methods {@link #appendColumns(AbstractComplexMatrix,AbstractComplexMatrix)
 * appendColumns}, {@link #appendColumns(AbstractComplexMatrix,AbstractComplexMatrix)
 * appendRows} and {@link #repeat(AbstractComplexMatrix,int,int) repeat} to append
 * rows and columns.</td>
 * </tr>
 * <tr align="left" valign="top">
 * <td><i> General block matrices </i></td>
 * <td>Use methods {@link #compose(List<List<AbstractComplexMatrix>>) compose} and
 * {@link #decompose(List<List<AbstractComplexMatrix>>,AbstractComplexMatrix) decompose} to work
 * with general block matrices.</td>
 * </tr>
 * <tr align="left" valign="top">
 * <td><i> Diagonal matrices </i></td>
 * <td>Use methods {@link #diagonal(AbstractComplexVector) diagonal(vector)},
 * {@link #diagonal(AbstractComplexMatrix) diagonal(matrix)} and
 * {@link #identityidentity} to work with diagonal matrices.</td>
 * </tr>
 * <tr align="left" valign="top">
 * <td><i> Diagonal block matrices </i></td>
 * <td>Use method
 * {@link #composeDiagonal(AbstractComplexMatrix,AbstractComplexMatrix,AbstractComplexMatrix)
 * composeDiagonal} to work with diagonal block matrices.</td>
 * </tr>
 * <tr align="left" valign="top">
 * <td><i>Random</i></td>
 * <td>Use methods {@link #random(int,int) random} and
 * {@link #sample(int,int,Float64List,double) sample} to construct random matrices.
 * </td>
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
 *  ComplexFactory2D F = ComplexFactory2D.dense;
 *  F.make(4,4);
 *  F.random(4,4);
 *  ...
 * </pre>
 *
 * </td>
 * </table>
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
//class DComplexFactory2D {
//
//  /**
//   * A factory producing dense matrices.
//   */
//  static final DComplexFactory2D dense = new DComplexFactory2D._internal();
//
//  /**
//   * A factory producing sparse hash matrices.
//   */
//  static final DComplexFactory2D sparse = new DComplexFactory2D._internal();
//
//  DComplexFactory2D._internal();

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
AbstractComplexMatrix appendColumns(AbstractComplexMatrix A, AbstractComplexMatrix B, AbstractComplexMatrix make(int, int)) {
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
  AbstractComplexMatrix matrix = make(r, ac + bc);
  matrix.part(0, 0, r, ac).copyFrom(A);
  matrix.part(0, ac, r, bc).copyFrom(B);
  return matrix;
}

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
AbstractComplexMatrix appendColumn(AbstractComplexMatrix A, AbstractComplexVector b, AbstractComplexMatrix make(int, int)) {
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
  AbstractComplexMatrix matrix = make(r, ac + bc);
  matrix.part(0, 0, r, ac).copyFrom(A);
  matrix.column(ac).copyFrom(b);
  return matrix;
}

/**
 * C = A||B; Constructs a new matrix which is the row-wise concatenation of
 * two other matrices.
 *
 */
AbstractComplexMatrix appendRows(AbstractComplexMatrix A, AbstractComplexMatrix B, AbstractComplexMatrix make(int, int)) {
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
  AbstractComplexMatrix matrix = make(ar + br, c);
  matrix.part(0, 0, ar, c).copyFrom(A);
  matrix.part(ar, 0, br, c).copyFrom(B);
  return matrix;
}

/**
 * C = A||b; Constructs a new matrix which is the row-wise concatenation of
 * two other matrices.
 *
 */
AbstractComplexMatrix appendRow(AbstractComplexMatrix A, AbstractComplexVector b, AbstractComplexMatrix make(int, int)) {
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
  AbstractComplexMatrix matrix = make(ar + br, c);
  matrix.part(0, 0, ar, c).copyFrom(A);
  matrix.row(ar).copyFrom(b);
  return matrix;
}

/**
 * Checks whether the given array is rectangular, that is, whether all rows
 * have the same number of columns.
 *
 * @throws ArgumentError
 *             if the array is not rectangular.
 */
void _checkRectangularShape(List<List<double>> array) {
  int columns = -1;
  for (int r = 0; r < array.length; r++) {
    if (array[r] != null) {
      if (columns == -1) columns = array[r].length;
      if (array[r].length != columns) {
        throw new ArgumentError("All rows of array must have same number of columns.");
      }
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
void _checkRectangularShape2D(List<List<AbstractComplexMatrix>> array) {
  int columns = -1;
  for (int r = 0; r < array.length; r++) {
    if (array[r] != null) {
      if (columns == -1) columns = array[r].length;
      if (array[r].length != columns) throw new ArgumentError("All rows of array must have same number of columns.");
    }
  }
}

/**
 * Constructs a block matrix made from the given parts. The inverse to
 * method {@link #decompose(List<List<AbstractComplexMatrix>>, AbstractComplexMatrix)}.
 * <p>
 * All matrices of a given column within <tt>parts</tt> must have the same
 * number of columns. All matrices of a given row within <tt>parts</tt> must
 * have the same number of rows. Otherwise an
 * <tt>ArgumentError</tt> is thrown. Note that <tt>null</tt>s
 * within <tt>parts[row,col]</tt> are an exception to this rule: they are
 * ignored. Cells are copied.
 *
 * @throws ArgumentError
 *             subject to the conditions outlined above.
 */
AbstractComplexMatrix compose(List<List<AbstractComplexMatrix>> parts, AbstractComplexMatrix make(int, int)) {
  _checkRectangularShape2D(parts);
  int rows = parts.length;
  int columns = 0;
  if (parts.length > 0) columns = parts[0].length;
  AbstractComplexMatrix empty = make(0, 0);

  if (rows == 0 || columns == 0) return empty;

  // determine maximum column width of each column
  Int32List maxWidths = new Int32List(columns);
  for (int c = 0; c < columns; c++) {
    int maxWidth = 0;
    for (int r = 0; r < rows; r++) {
      AbstractComplexMatrix part = parts[r][c];
      if (part != null) {
        int width = part.columns;
        if (maxWidth > 0 && width > 0 && width != maxWidth) throw new ArgumentError("Different number of columns.");
        maxWidth = Math.max(maxWidth, width);
      }
    }
    maxWidths[c] = maxWidth;
  }

  // determine row height of each row
  Int32List maxHeights = new Int32List(rows);
  for (int r = 0; r < rows; r++) {
    int maxHeight = 0;
    for (int c = 0; c < columns; c++) {
      AbstractComplexMatrix part = parts[r][c];
      if (part != null) {
        int height = part.rows;
        if (maxHeight > 0 && height > 0 && height != maxHeight) throw new ArgumentError("Different number of rows.");
        maxHeight = Math.max(maxHeight, height);
      }
    }
    maxHeights[r] = maxHeight;
  }

  // shape of result
  int resultRows = 0;
  for (int r = 0; r < rows; r++) resultRows += maxHeights[r];
  int resultCols = 0;
  for (int c = 0; c < columns; c++) resultCols += maxWidths[c];

  AbstractComplexMatrix matrix = make(resultRows, resultCols);

  // copy
  int idxr = 0;
  for (int r = 0; r < rows; r++) {
    int idxc = 0;
    for (int c = 0; c < columns; c++) {
      AbstractComplexMatrix part = parts[r][c];
      if (part != null) {
        matrix.part(idxr, idxc, part.rows, part.columns).copyFrom(part);
      }
      idxc += maxWidths[c];
    }
    idxr += maxHeights[r];
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
AbstractComplexMatrix composeDiagonal(AbstractComplexMatrix A, AbstractComplexMatrix B, AbstractComplexMatrix make(int, int)) {
  int ar = A.rows;
  int ac = A.columns;
  int br = B.rows;
  int bc = B.columns;
  AbstractComplexMatrix sum = make(ar + br, ac + bc);
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
AbstractComplexMatrix composeDiagonal3(AbstractComplexMatrix A, AbstractComplexMatrix B, AbstractComplexMatrix C, AbstractComplexMatrix make(int, int)) {
  AbstractComplexMatrix diag = make(A.rows + B.rows + C.rows, A.columns + B.columns + C.columns);
  diag.part(0, 0, A.rows, A.columns).copyFrom(A);
  diag.part(A.rows, A.columns, B.rows, B.columns).copyFrom(B);
  diag.part(A.rows + B.rows, A.columns + B.columns, C.rows, C.columns).copyFrom(C);
  return diag;
}

/**
 * Constructs a bidiagonal block matrix from the given parts.
 *
 * from the given parts. Cells are copied.
 */
AbstractComplexMatrix composeBidiagonal(AbstractComplexMatrix A, AbstractComplexMatrix B, AbstractComplexMatrix make(int, int)) {
  int ar = A.rows;
  int ac = A.columns;
  int br = B.rows;
  int bc = B.columns;
  AbstractComplexMatrix sum = make(ar + br - 1, ac + bc);
  sum.part(0, 0, ar, ac).copyFrom(A);
  sum.part(ar - 1, ac, br, bc).copyFrom(B);
  return sum;
}

/**
 * Splits a block matrix into its constituent blocks; Copies blocks of a
 * matrix into the given parts. The inverse to method
 * {@link #compose(List<List<AbstractComplexMatrix>>)}.
 * <p>
 * All matrices of a given column within <tt>parts</tt> must have the same
 * number of columns. All matrices of a given row within <tt>parts</tt> must
 * have the same number of rows. Otherwise an
 * <tt>ArgumentError</tt> is thrown. Note that <tt>null</tt>s
 * within <tt>parts[row,col]</tt> are an exception to this rule: they are
 * ignored. Cells are copied.
 *
 * @throws ArgumentError
 *             subject to the conditions outlined above.
 */
void decompose(List<List<AbstractComplexMatrix>> parts, AbstractComplexMatrix matrix) {
  _checkRectangularShape2D(parts);
  int rows = parts.length;
  int columns = 0;
  if (parts.length > 0) columns = parts[0].length;
  if (rows == 0 || columns == 0) return;

  // determine maximum column width of each column
  Int32List maxWidths = new Int32List(columns);
  for (int c = 0; c < columns; c++) {
    int maxWidth = 0;
    for (int r = 0; r < rows; r++) {
      AbstractComplexMatrix part = parts[r][c];
      if (part != null) {
        int width = part.columns;
        if (maxWidth > 0 && width > 0 && width != maxWidth) throw new ArgumentError("Different number of columns.");
        maxWidth = Math.max(maxWidth, width);
      }
    }
    maxWidths[c] = maxWidth;
  }

  // determine row height of each row
  Int32List maxHeights = new Int32List(rows);
  for (int r = 0; r < rows; r++) {
    int maxHeight = 0;
    for (int c = 0; c < columns; c++) {
      AbstractComplexMatrix part = parts[r][c];
      if (part != null) {
        int height = part.rows;
        if (maxHeight > 0 && height > 0 && height != maxHeight) throw new ArgumentError("Different number of rows.");
        maxHeight = Math.max(maxHeight, height);
      }
    }
    maxHeights[r] = maxHeight;
  }

  // shape of result parts
  int resultRows = 0;
  for (int r = 0; r < rows; r++) resultRows += maxHeights[r];
  int resultCols = 0;
  for (int c = 0; c < columns; c++) resultCols += maxWidths[c];

  if (matrix.rows < resultRows || matrix.columns < resultCols) throw new ArgumentError("Parts larger than matrix.");

  // copy
  int idxr = 0;
  for (int r = 0; r < rows; r++) {
    int idxc = 0;
    for (int c = 0; c < columns; c++) {
      AbstractComplexMatrix part = parts[r][c];
      if (part != null) {
        part.copyFrom(matrix.part(idxr, idxc, part.rows, part.columns));
      }
      idxc += maxWidths[c];
    }
    idxr += maxHeights[r];
  }

}

/**
 * Demonstrates usage of this class.
 */
//void demo1() {
//  print("\n\n");
//  List<List<AbstractComplexMatrix>> parts1 = [[null, fill(2, 2, [1.0, 2.0]), null],
//                                              [fill(4, 4, [3.0, 4.0]), null, fill(4, 3, [5.0, 6.0])],
//                                              [null, fill(2, 2, [7.0, 8.0]), null]];
//  print("\n${compose(parts1)}");
//}

/**
 * Demonstrates usage of this class.
 */
//void demo2() {
//  print("\n\n");
//  AbstractComplexMatrix matrix;
//  AbstractComplexMatrix A, B, C, D;
//  AbstractComplexMatrix _ = null;
//  A = fill(2, 2, new Float64List.fromList([1, 2]));
//  B = fill(4, 4, new Float64List.fromList([3, 4]));
//  C = fill(4, 3, new Float64List.fromList([5, 6]));
//  D = fill(2, 2, new Float64List.fromList([7, 8]));
//  List<List<AbstractComplexMatrix>> parts1 = [[_, A, _], [B, _, C], [_, D, _]];
//  matrix = compose(parts1);
//  print("\n${matrix}");
//
//  final nine = new Float64List.fromList([9.0, 9.0]);
//  A.setAll(nine);
//  B.setAll(nine);
//  C.setAll(nine);
//  D.setAll(nine);
//  decompose(parts1, matrix);
//  print(A);
//  print(B);
//  print(C);
//  print(D);
//}

/**
 * Constructs a new diagonal matrix whose diagonal elements are the elements
 * of <tt>vector</tt>. Cells values are copied. The new matrix is not a
 * view.
 *
 * @return a new matrix.
 */
AbstractComplexMatrix diag(AbstractComplexVector vector, AbstractComplexMatrix make(int, int)) {
  int size = vector.length;
  AbstractComplexMatrix diag = make(size, size);
  for (int i = 0; i < size; i++) {
    diag.set(i, i, vector.get(i));
  }
  return diag;
}

/**
 * Constructs a new vector consisting of the diagonal elements of <tt>A</tt>
 * . Cells values are copied. The new vector is not a view.
 *
 * @param A
 *            the matrix, need not be square.
 * @return a new vector.
 */
AbstractComplexVector diagonal(AbstractComplexMatrix A, AbstractComplexMatrix make(int, int)) {
  int min = Math.min(A.rows, A.columns);
  AbstractComplexVector diag = _make1D(min, make);
  for (int i = 0; i < min; i++) {
    diag.set(i, A.get(i, i));
  }
  return diag;
}

/**
 * Constructs an identity matrix (having ones on the diagonal and zeros
 * elsewhere).
 */
AbstractComplexMatrix identity(int rowsAndColumns, AbstractComplexMatrix make(int, int)) {
  AbstractComplexMatrix matrix = make(rowsAndColumns, rowsAndColumns);
  Float64List one = new Float64List.fromList([1.0, 0.0]);
  for (int i = rowsAndColumns; --i >= 0; ) {
    matrix.set(i, i, one);
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
//  AbstractComplexMatrix fromList(List<List<double>> values) {
//    if (this == sparse) {
//      return new SparseComplexMatrix.fromList(values);
//    } else {
//      return new ComplexMatrix.fromList(values);
//    }
//  }

/**
 * Constructs a matrix with the given shape, each cell initialized with
 * zero.
 */
//  AbstractComplexMatrix make(int rows, int columns) {
//    if (this == sparse) {
//      return new SparseComplexMatrix(rows, columns);
//    } else {
//      return new ComplexMatrix(rows, columns);
//    }
//  }

/**
 * Constructs a matrix with the given shape, each cell initialized with the
 * given value.
 */
AbstractComplexMatrix fillMatrix(int rows, int columns, List<double> initialValue, AbstractComplexMatrix make(int, int)) {
  if (initialValue[0] == 0 && initialValue[1] == 0) {
    return make(rows, columns);
  }
  return make(rows, columns)..setAll(initialValue);
}

/**
 * Constructs a 1d matrix of the right dynamic type.
 */
AbstractComplexVector _make1D(int size, AbstractComplexMatrix make(int, int)) {
  return make(0, 0).like1D(size);
}

/**
 * Constructs a matrix with uniformly distributed values in <tt>(0,1)</tt>
 * (exclusive).
 */
AbstractComplexMatrix randomMatrix(int rows, int columns, AbstractComplexMatrix make(int, int)) {
  return make(rows, columns)..forEach(cfunc.random);
}

/**
 * C = A||A||..||A; Constructs a new matrix which is duplicated both along
 * the row and column dimension.
 */
AbstractComplexMatrix repeatMatrix(AbstractComplexMatrix A, int rowRepeat, int columnRepeat, AbstractComplexMatrix make(int, int)) {
  int r = A.rows;
  int c = A.columns;
  AbstractComplexMatrix matrix = make(r * rowRepeat, c * columnRepeat);
  for (int i = 0; i < rowRepeat; i++) {
    for (int j = 0; j < columnRepeat; j++) {
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
 * @see cern.jet.random.tdouble.sampling.DoubleRandomSampler
 */
/*AbstractComplexMatrix sample(int rows, int columns, Float64List value, double nonZeroFraction) {
  AbstractComplexMatrix matrix = make(rows, columns);
  sampleMatrix(matrix, value, nonZeroFraction);
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
 * @see cern.jet.random.tdouble.sampling.DoubleRandomSampler
 */
/*AbstractComplexMatrix sampleMatrix(AbstractComplexMatrix matrix, Float64List value, double nonZeroFraction) {
  int rows = matrix.rows;
  int columns = matrix.columns;
  double epsilon = 1e-09;
  if (nonZeroFraction < 0 - epsilon || nonZeroFraction > 1 + epsilon) {
    throw new ArgumentError();
  }
  if (nonZeroFraction < 0) {
    nonZeroFraction = 0.0;
  }
  if (nonZeroFraction > 1) {
    nonZeroFraction = 1.0;
  }

  matrix.setAll(0, 0);

  int size = rows * columns;
  int n = (size * nonZeroFraction).round();
  if (n == 0) {
    return matrix;
  }

  DoubleRandomSamplingAssistant sampler = new DoubleRandomSamplingAssistant(n, size, new DoubleMersenneTwister());
  for (int i = 0; i < size; i++) {
    if (sampler.sampleNextElement()) {
      int row = (i ~/ columns);
      int column = (i % columns);
      matrix.set(row, column, value);
    }
  }

  return matrix;
}*/
//}
