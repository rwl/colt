/*
Copyright (C) 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
is hereby granted without fee, provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear in supporting documentation.
CERN makes no representations about the suitability of this software for any purpose.
It is provided "as is" without expressed or implied warranty.
 */
part of cern.colt.matrix.tdouble.algo;

/*import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.Future;

import cern.colt.matrix.AbstractFormatter;
import cern.colt.matrix.tdouble.DoubleVector;
import cern.colt.matrix.tdouble.DoubleMatrix;
import cern.colt.matrix.tdouble.DoubleMatrix3D;
import cern.colt.matrix.tdouble.impl.DenseColumnDoubleMatrix;
import cern.colt.matrix.tdouble.impl.DenseDoubleVector;
import cern.colt.matrix.tdouble.impl.DenseDoubleMatrix;
import cern.colt.matrix.tdouble.impl.SparseCCDoubleMatrix;
import cern.colt.matrix.tdouble.impl.SparseDoubleVector;
import cern.colt.matrix.tdouble.impl.SparseRCDoubleMatrix;
import cern.jet.math.tdouble.DoubleFunctions;
import edu.emory.mathcs.util.ConcurrencyUtils;*/

/**
 * Tests matrices for linear algebraic properties (equality, tridiagonality,
 * symmetry, singularity, etc).
 * <p>
 * Except where explicitly indicated, all methods involving equality tests (
 * <tt>==</tt>) allow for numerical instability, to a degree specified upon
 * instance construction and returned by method {@link #tolerance()}. The public
 * static final variable <tt>DEFAULT</tt> represents a default Property object
 * with a tolerance of <tt>1.0E-9</tt>. The static final variable
 * <tt>ZERO</tt> represents a Property object with a tolerance of <tt>0.0</tt>.
 * The static final variable <tt>TWELVE</tt> represents a Property object
 * with a tolerance of <tt>1.0E-12</tt>. As long as you are happy with these
 * tolerances, there is no need to construct Property objects. Simply use idioms
 * like <tt>Property.DEFAULT.equals(A,B)</tt>,
 * <tt>Property.ZERO.equals(A,B)</tt>, <tt>Property.TWELVE.equals(A,B)</tt>.
 * <p>
 * To work with a different tolerance (e.g. <tt>1.0E-15</tt> or <tt>1.0E-5</tt>)
 * use the constructor and/or method {@link #setTolerance(double)}. Note that
 * the static final Property objects are immutable: Is is not possible to
 * alter their tolerance. Any attempt to do so will throw an Exception.
 * <p>
 * Note that this implementation is not synchronized.
 * <p>
 * Example: <tt>equals(DoubleMatrix A, DoubleMatrix B)</tt> is defined as
 * follows
 * <table>
 * <td class="PRE">
 *
 * <pre>
 *  { some other tests not related to tolerance go here }
 *  double epsilon = tolerance();
 *  for (int row=rows; --row &gt;= 0;) {
 *     for (int column=columns; --column &gt;= 0;) {
 *        //if (!(A.getQuick(row,column) == B.getQuick(row,column))) return false;
 *        if (Math.abs(A.getQuick(row,column) - B.getQuick(row,column)) &gt; epsilon) return false;
 *     }
 *  }
 *  return true;
 * </pre>
 *
 * </td>
 * </table>
 * Here are some example properties
 * <table border="1" cellspacing="0">
 * <tr align="left" valign="top">
 * <td valign="middle" align="left"><tt>matrix</tt></td>
 * <td> <tt>4&nbsp;x&nbsp;4&nbsp;<br>
 0&nbsp;0&nbsp;0&nbsp;0<br>
 0&nbsp;0&nbsp;0&nbsp;0<br>
 0&nbsp;0&nbsp;0&nbsp;0<br>
 0&nbsp;0&nbsp;0&nbsp;0 </tt></td>
 * <td><tt>4&nbsp;x&nbsp;4<br>
 1&nbsp;0&nbsp;0&nbsp;0<br>
 0&nbsp;0&nbsp;0&nbsp;0<br>
 0&nbsp;0&nbsp;0&nbsp;0<br>
 0&nbsp;0&nbsp;0&nbsp;1 </tt></td>
 * <td><tt>4&nbsp;x&nbsp;4<br>
 1&nbsp;1&nbsp;0&nbsp;0<br>
 1&nbsp;1&nbsp;1&nbsp;0<br>
 0&nbsp;1&nbsp;1&nbsp;1<br>
 0&nbsp;0&nbsp;1&nbsp;1 </tt></td>
 * <td><tt> 4&nbsp;x&nbsp;4<br>
 0&nbsp;1&nbsp;1&nbsp;1<br>
 0&nbsp;1&nbsp;1&nbsp;1<br>
 0&nbsp;0&nbsp;0&nbsp;1<br>
 0&nbsp;0&nbsp;0&nbsp;1 </tt></td>
 * <td><tt> 4&nbsp;x&nbsp;4<br>
 0&nbsp;0&nbsp;0&nbsp;0<br>
 1&nbsp;1&nbsp;0&nbsp;0<br>
 1&nbsp;1&nbsp;0&nbsp;0<br>
 1&nbsp;1&nbsp;1&nbsp;1 </tt></td>
 * <td><tt>4&nbsp;x&nbsp;4<br>
 1&nbsp;1&nbsp;0&nbsp;0<br>
 0&nbsp;1&nbsp;1&nbsp;0<br>
 0&nbsp;1&nbsp;0&nbsp;1<br>
 1&nbsp;0&nbsp;1&nbsp;1 </tt><tt> </tt></td>
 * <td><tt>4&nbsp;x&nbsp;4<br>
 1&nbsp;1&nbsp;1&nbsp;0<br>
 0&nbsp;1&nbsp;0&nbsp;0<br>
 1&nbsp;1&nbsp;0&nbsp;1<br>
 0&nbsp;0&nbsp;1&nbsp;1 </tt></td>
 * </tr>
 * <tr align="center" valign="middle">
 * <td><tt>upperBandwidth</tt></td>
 * <td><div align="center"><tt>0</tt></div></td>
 * <td><div align="center"><tt>0</tt></div></td>
 * <td><div align="center"><tt>1</tt></div></td>
 * <td><tt>3</tt></td>
 * <td align="center" valign="middle"><tt>0</tt></td>
 * <td align="center" valign="middle"><div align="center"><tt>1</tt></div></td>
 * <td align="center" valign="middle"><div align="center"><tt>2</tt></div></td>
 * </tr>
 * <tr align="center" valign="middle">
 * <td><tt>lowerBandwidth</tt></td>
 * <td><div align="center"><tt>0</tt></div></td>
 * <td><div align="center"><tt>0</tt></div></td>
 * <td><div align="center"><tt>1</tt></div></td>
 * <td><tt>0</tt></td>
 * <td align="center" valign="middle"><tt>3</tt></td>
 * <td align="center" valign="middle"><div align="center"><tt>3</tt></div></td>
 * <td align="center" valign="middle"><div align="center"><tt>2</tt></div></td>
 * </tr>
 * <tr align="center" valign="middle">
 * <td><tt>semiBandwidth</tt></td>
 * <td><div align="center"><tt>1</tt></div></td>
 * <td><div align="center"><tt>1</tt></div></td>
 * <td><div align="center"><tt>2</tt></div></td>
 * <td><tt>4</tt></td>
 * <td align="center" valign="middle"><tt>4</tt></td>
 * <td align="center" valign="middle"><div align="center"><tt>4</tt></div></td>
 * <td align="center" valign="middle"><div align="center"><tt>3</tt></div></td>
 * </tr>
 * <tr align="center" valign="middle">
 * <td><tt>description</tt></td>
 * <td><div align="center"><tt>zero</tt></div></td>
 * <td><div align="center"><tt>diagonal</tt></div></td>
 * <td><div align="center"><tt>tridiagonal</tt></div></td>
 * <td><tt>upper triangular</tt></td>
 * <td align="center" valign="middle"><tt>lower triangular</tt></td>
 * <td align="center" valign="middle"><div align="center"><tt>unstructured</tt>
 * </div></td>
 * <td align="center" valign="middle"><div align="center"><tt>unstructured</tt>
 * </div></td>
 * </tr>
 * </table>
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.1, 28/May/2000 (fixed strange bugs involving NaN, -inf, inf)
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
class DoubleProperty {//extends cern.colt.PersistentObject {

  /**
   * The default Property object; currently has <tt>tolerance()==1.0E-9</tt>.
   */
  static final DoubleProperty DEFAULT = new DoubleProperty(1.0E-9);

  /**
   * A Property object with <tt>tolerance()==0.0</tt>.
   */
  static final DoubleProperty ZERO = new DoubleProperty(0.0);

  /**
   * A Property object with <tt>tolerance()==1.0E-12</tt>.
   */
  static final DoubleProperty TWELVE = new DoubleProperty(1.0E-12);

  double _tolerance;

  /**
   * Not instantiable by no-arg constructor.
   */
  //    private DoubleProperty() {
  //        this(1.0E-9); // just to be on the safe side
  //    }

  /**
   * Constructs an instance with a tolerance of
   * <tt>Math.abs(newTolerance)</tt>.
   */
  DoubleProperty([double newTolerance = 1.0E-9]) {
    _tolerance = newTolerance.abs();
  }

  /**
   * Returns a String with <tt>length</tt> blanks.
   */
  static String _blanks(int length) {
    if (length < 0) length = 0;
    StringBuffer buf = new StringBuffer(length);
    for (int k = 0; k < length; k++) {
      buf.write(' ');
    }
    return buf.toString();
  }

  /**
   * Checks whether the given matrix <tt>A</tt> is <i>rectangular</i>.
   *
   * @throws ArgumentError
   *             if <tt>A.rows() < A.columns()</tt>.
   */
  void checkRectangular(DoubleMatrix A) {
    if (A.rows < A.columns) {
      throw new ArgumentError("Matrix must be rectangular: " + AbstractFormatter.shape2D(A));
    }
  }

  /**
   * Checks whether the given matrix <tt>A</tt> is <i>square</i>.
   *
   * @throws ArgumentError
   *             if <tt>A.rows() != A.columns()</tt>.
   */
  void checkSquare(DoubleMatrix A) {
    if (A.rows != A.columns) {
      throw new ArgumentError("Matrix must be square: " + AbstractFormatter.shape2D(A));
    }
  }

  /*void checkDense2D(DoubleMatrix A) {
    if (!(A is DenseDoubleMatrix) && !(A is DenseColumnDoubleMatrix)) {
      throw new ArgumentError("Matrix must be dense");
    }
  }*/

  void checkDense(DoubleVector A) {
    if (!(A is DenseDoubleVector)) {
      throw new ArgumentError("Matrix must be dense");
    }
  }

  void checkSparse(DoubleVector A) {
    if (!(A is SparseDoubleVector)) {
      throw new ArgumentError("Matrix must be sparse");
    }
  }

  void checkSparse2D(DoubleMatrix A) {
    //        if (!(A is SparseDoubleMatrix) && !(A is RCDoubleMatrix) && !(A is RCMDoubleMatrix)
    //                && !(A is CCDoubleMatrix) && !(A is CCMDoubleMatrix))
    if (!(A is SparseCCDoubleMatrix) && !(A is SparseRCDoubleMatrix)) {
      throw new ArgumentError("Matrix must be sparse");
    }
  }

  /**
   * Returns the matrix's fraction of non-zero cells;
   * <tt>A.cardinality() / A.size()</tt>.
   */
  double density(DoubleMatrix A) {
    return A.cardinality / (A.length as double);
  }

  /**
   * Returns whether all cells of the given matrix <tt>A</tt> are equal to the
   * given value. The result is <tt>true</tt> if and only if
   * <tt>A != null</tt> and <tt>! (Math.abs(value - A[i]) > tolerance())</tt>
   * holds for all coordinates.
   *
   * @param A
   *            the first matrix to compare.
   * @param value
   *            the value to compare against.
   * @return <tt>true</tt> if the matrix is equal to the value; <tt>false</tt>
   *         otherwise.
   */
  bool equals(final DoubleVector A, final double value) {
    if (A == null) {
      return false;
    }
    int size = A.length;
    final double epsilon = tolerance();
    bool result = false;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, size);
      List<Future> futures = new List<Future>(nthreads);
      List<bool> results = new List<bool>(nthreads);
      int k = size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int i = firstIdx; i < lastIdx; i++) {
            double x = A.getQuick(i);
            double diff = (value - x).abs();
            if ((diff != diff) && ((value != value && x != x) || value == x)) diff = 0.0;
            if (!(diff <= epsilon)) {
              return false;
            }
          }
          return true;
        });
      }
      //try {
      for (int j = 0; j < nthreads; j++) {
        results[j] = futures[j].get() as bool;
      }
      result = results[0].booleanValue();
      for (int j = 1; j < nthreads; j++) {
        result = result && results[j].booleanValue();
      }
      /*} on ExecutionException catch (ex) {
                ex.printStackTrace();
            } on InterruptedException catch (e) {
                e.printStackTrace();
            }*/
      return result;
    } else {*/
      for (int i = 0; i < size; i++) {
        double x = A.get(i);
        double diff = (value - x).abs();
        if ((diff != diff) && ((value != value && x != x) || value == x)) {
          diff = 0.0;
        }
        if (!(diff <= epsilon)) {
          return false;
        }
      }
      return true;
    //}
  }

  /**
   * Returns whether both given matrices <tt>A</tt> and <tt>B</tt> are equal.
   * The result is <tt>true</tt> if <tt>A==B</tt>. Otherwise, the result is
   * <tt>true</tt> if and only if both arguments are <tt>!= null</tt>, have
   * the same size and <tt>! (Math.abs(A[i] - B[i]) > tolerance())</tt> holds
   * for all indexes.
   *
   * @param A
   *            the first matrix to compare.
   * @param B
   *            the second matrix to compare.
   * @return <tt>true</tt> if both matrices are equal; <tt>false</tt>
   *         otherwise.
   */
  bool equalsVector(final DoubleVector A, final DoubleVector B) {
    if (identical(A, B)) {
      return true;
    }
    if (!(A != null && B != null)) {
      return false;
    }
    int size = A.length;
    if (size != B.length) {
      return false;
    }

    final double epsilon = tolerance();
    bool result = false;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, size);
      List<Future> futures = new List<Future>(nthreads);
      List<bool> results = new List<bool>(nthreads);
      int k = size ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int i = firstIdx; i < lastIdx; i++) {
            double x = A.getQuick(i);
            double value = B.getQuick(i);
            double diff = (value - x).abs();
            if ((diff != diff) && ((value != value && x != x) || value == x)) {
              diff = 0.0;
            }
            if (!(diff <= epsilon)) {
              return false;
            }
          }
          return true;
        });
      }
      //try {
      for (int j = 0; j < nthreads; j++) {
        results[j] = futures[j].get() as bool;
      }
      result = results[0].booleanValue();
      for (int j = 1; j < nthreads; j++) {
        result = result && results[j].booleanValue();
      }
      } on ExecutionException catch (ex) {
                ex.printStackTrace();
            } on InterruptedException catch (e) {
                e.printStackTrace();
            }
      return result;
    } else {*/
      for (int i = 0; i < size; i++) {
        double x = A.get(i);
        double value = B.get(i);
        double diff = (value - x).abs();
        if ((diff != diff) && ((value != value && x != x) || value == x)) {
          diff = 0.0;
        }
        if (!(diff <= epsilon)) {
          return false;
        }
      }
      return true;
    //}
  }

  /**
   * Returns whether all cells of the given matrix <tt>A</tt> are equal to the
   * given value. The result is <tt>true</tt> if and only if
   * <tt>A != null</tt> and
   * <tt>! (Math.abs(value - A[row,col]) > tolerance())</tt> holds for all
   * coordinates.
   *
   * @param A
   *            the first matrix to compare.
   * @param value
   *            the value to compare against.
   * @return <tt>true</tt> if the matrix is equal to the value; <tt>false</tt>
   *         otherwise.
   */
  bool equalsMatrixValue(final DoubleMatrix A, final double value) {
    if (A == null) {
      return false;
    }
    final int rows = A.rows;
    final int columns = A.columns;
    bool result = false;
    final double epsilon = tolerance();
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (A.size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, A.rows());
      List<Future> futures = new List<Future>(nthreads);
      List<bool> results = new List<bool>(nthreads);
      int k = A.rows() ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? A.rows() : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = 0; c < columns; c++) {
              double x = A.getQuick(r, c);
              double diff = (value - x).abs();
              if ((diff != diff) && ((value != value && x != x) || value == x)) diff = 0;
              if (!(diff <= epsilon)) {
                return false;
              }
            }
          }
          return true;
        });
      }
      //try {
      for (int j = 0; j < nthreads; j++) {
        results[j] = futures[j].get() as bool;
      }
      result = results[0].booleanValue();
      for (int j = 1; j < nthreads; j++) {
        result = result && results[j].booleanValue();
      }
      } on ExecutionException catch (ex) {
                ex.printStackTrace();
            } on InterruptedException catch (e) {
                e.printStackTrace();
            }
      return result;
    } else {*/
      for (int r = 0; r < rows; r++) {
        for (int c = 0; c < columns; c++) {
          double x = A.get(r, c);
          double diff = (value - x).abs();
          if ((diff != diff) && ((value != value && x != x) || value == x)) {
            diff = 0.0;
          }
          if (!(diff <= epsilon)) {
            return false;
          }
        }
      }
      return true;
    //}
  }

  /**
   * Returns whether both given matrices <tt>A</tt> and <tt>B</tt> are equal.
   * The result is <tt>true</tt> if <tt>A==B</tt>. Otherwise, the result is
   * <tt>true</tt> if and only if both arguments are <tt>!= null</tt>, have
   * the same number of columns and rows and
   * <tt>! (Math.abs(A[row,col] - B[row,col]) > tolerance())</tt> holds for
   * all coordinates.
   *
   * @param A
   *            the first matrix to compare.
   * @param B
   *            the second matrix to compare.
   * @return <tt>true</tt> if both matrices are equal; <tt>false</tt>
   *         otherwise.
   */
  bool equalsMatrix(final DoubleMatrix A, final DoubleMatrix B) {
    if (identical(A, B)) {
      return true;
    }
    if (!(A != null && B != null)) {
      return false;
    }
    final int rows = A.rows;
    final int columns = A.columns;
    if (columns != B.columns || rows != B.rows) {
      return false;
    }
    bool result = false;
    final double epsilon = tolerance();
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (A.size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, A.rows());
      List<Future> futures = new List<Future>(nthreads);
      List<bool> results = new List<bool>(nthreads);
      int k = A.rows() ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? A.rows() : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = 0; c < columns; c++) {
              double x = A.getQuick(r, c);
              double value = B.getQuick(r, c);
              double diff = (value - x).abs();
              if ((diff != diff) && ((value != value && x != x) || value == x)) diff = 0;
              if (!(diff <= epsilon)) {
                return false;
              }
            }
          }
          return true;
        });
      }
      //try {
      for (int j = 0; j < nthreads; j++) {
        results[j] = futures[j].get() as bool;
      }
      result = results[0].booleanValue();
      for (int j = 1; j < nthreads; j++) {
        result = result && results[j].booleanValue();
      }
      } on ExecutionException catch (ex) {
                ex.printStackTrace();
            } on InterruptedException catch (e) {
                e.printStackTrace();
            }
      return result;
    } else {*/
      for (int r = 0; r < rows; r++) {
        for (int c = 0; c < columns; c++) {
          double x = A.get(r, c);
          double value = B.get(r, c);
          double diff = (value - x).abs();
          if ((diff != diff) && ((value != value && x != x) || value == x)) diff = 0.0;
          if (!(diff <= epsilon)) {
            return false;
          }
        }
      }
      return true;
    //}
  }

  /**
   * Returns whether all cells of the given matrix <tt>A</tt> are equal to the
   * given value. The result is <tt>true</tt> if and only if
   * <tt>A != null</tt> and
   * <tt>! (Math.abs(value - A[slice,row,col]) > tolerance())</tt> holds for
   * all coordinates.
   *
   * @param A
   *            the first matrix to compare.
   * @param value
   *            the value to compare against.
   * @return <tt>true</tt> if the matrix is equal to the value; <tt>false</tt>
   *         otherwise.
   */
  /*bool equalsMatrix3DValue(final DoubleMatrix3D A, final double value) {
    if (A == null) return false;
    final int slices = A.slices();
    final int rows = A.rows();
    final int columns = A.columns();
    bool result = false;
    final double epsilon = tolerance();
    int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (A.size() >= ConcurrencyUtils.getThreadsBeginN_3D())) {
      nthreads = Math.min(nthreads, slices);
      List<Future> futures = new List<Future>(nthreads);
      List<bool> results = new List<bool>(nthreads);
      int k = slices ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstSlice = j * k;
        final int lastSlice = (j == nthreads - 1) ? slices : firstSlice + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int s = firstSlice; s < lastSlice; s++) {
            for (int r = 0; r < rows; r++) {
              for (int c = 0; c < columns; c++) {
                double x = A.getQuick(s, r, c);
                double diff = (value - x).abs();
                if ((diff != diff) && ((value != value && x != x) || value == x)) diff = 0.0;
                if (!(diff <= epsilon)) {
                  return false;
                }
              }
            }
          }
          return true;
        });
      }
      //try {
      for (int j = 0; j < nthreads; j++) {
        results[j] = futures[j].get() as bool;
      }
      result = results[0].booleanValue();
      for (int j = 1; j < nthreads; j++) {
        result = result && results[j].booleanValue();
      }
      } on ExecutionException catch (ex) {
                ex.printStackTrace();
            } on InterruptedException catch (e) {
                e.printStackTrace();
            }
      return result;
    } else {
      for (int s = 0; s < slices; s++) {
        for (int r = 0; r < rows; r++) {
          for (int c = 0; c < columns; c++) {
            double x = A.getQuick(s, r, c);
            double diff = (value - x).abs();
            if ((diff != diff) && ((value != value && x != x) || value == x)) diff = 0.0;
            if (!(diff <= epsilon)) {
              return false;
            }
          }
        }
      }
      return true;
    }
  }*/

  /**
   * Returns whether both given matrices <tt>A</tt> and <tt>B</tt> are equal.
   * The result is <tt>true</tt> if <tt>A==B</tt>. Otherwise, the result is
   * <tt>true</tt> if and only if both arguments are <tt>!= null</tt>, have
   * the same number of columns, rows and slices, and
   * <tt>! (Math.abs(A[slice,row,col] - B[slice,row,col]) > tolerance())</tt>
   * holds for all coordinates.
   *
   * @param A
   *            the first matrix to compare.
   * @param B
   *            the second matrix to compare.
   * @return <tt>true</tt> if both matrices are equal; <tt>false</tt>
   *         otherwise.
   */
  /*bool equalsMatrix3D(final DoubleMatrix3D A, final DoubleMatrix3D B) {
    if (A == B) return true;
    if (!(A != null && B != null)) return false;
    final int slices = A.slices();
    final int rows = A.rows();
    final int columns = A.columns();
    if (columns != B.columns() || rows != B.rows() || slices != B.slices()) return false;
    bool result = false;
    final double epsilon = tolerance();
    int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (A.size() >= ConcurrencyUtils.getThreadsBeginN_3D())) {
      nthreads = Math.min(nthreads, slices);
      List<Future> futures = new List<Future>(nthreads);
      List<bool> results = new List<bool>(nthreads);
      int k = slices ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstSlice = j * k;
        final int lastSlice = (j == nthreads - 1) ? slices : firstSlice + k;
        futures[j] = ConcurrencyUtils.submit(() {
          for (int s = firstSlice; s < lastSlice; s++) {
            for (int r = 0; r < rows; r++) {
              for (int c = 0; c < columns; c++) {
                double x = A.getQuick(s, r, c);
                double value = B.getQuick(s, r, c);
                double diff = (value - x).abs();
                if ((diff != diff) && ((value != value && x != x) || value == x)) diff = 0.0;
                if (!(diff <= epsilon)) {
                  return false;
                }
              }
            }
          }
          return true;
        });
      }
      //try {
      for (int j = 0; j < nthreads; j++) {
        results[j] = futures[j].get() as bool;
      }
      result = results[0].booleanValue();
      for (int j = 1; j < nthreads; j++) {
        result = result && results[j].booleanValue();
      }
      } on ExecutionException catch (ex) {
                ex.printStackTrace();
            } on InterruptedException catch (e) {
                e.printStackTrace();
            }
      return result;
    } else {
      for (int s = 0; s < slices; s++) {
        for (int r = 0; r < rows; r++) {
          for (int c = 0; c < columns; c++) {
            double x = A.getQuick(s, r, c);
            double value = B.getQuick(s, r, c);
            double diff = (value - x).abs();
            if ((diff != diff) && ((value != value && x != x) || value == x)) diff = 0.0;
            if (!(diff <= epsilon)) {
              return false;
            }
          }
        }
      }
      return true;
    }
  }*/

  /**
   * Modifies the given matrix square matrix <tt>A</tt> such that it is
   * diagonally dominant by row and column, hence non-singular, hence
   * invertible. For testing purposes only.
   *
   * @param A
   *            the square matrix to modify.
   * @throws ArgumentError
   *             if <tt>!isSquare(A)</tt>.
   */
  void generateNonSingular(DoubleMatrix A) {
    checkSquare(A);
    //cern.jet.math.tdouble.DoubleFunctions F = cern.jet.math.tdouble.DoubleFunctions.functions;
    int min = Math.min(A.rows, A.columns);
    for (int i = min; --i >= 0; ) {
      A.set(i, i, 0.0);
    }
    for (int i = min; --i >= 0; ) {
      double rowSum = A.row(i).reduce(plus, abs);
      double colSum = A.column(i).reduce(plus, abs);
      A.set(i, i, Math.max(rowSum, colSum) + i + 1);
    }
  }

  static String _get(List list, int index) {
    return list[index] as String;
  }

  /**
   * A matrix <tt>A</tt> is <i>diagonal</i> if <tt>A[i,j] == 0</tt> whenever
   * <tt>i != j</tt>. Matrix may but need not be square.
   */
  bool isDiagonal(DoubleMatrix A) {
    double epsilon = tolerance();
    int rows = A.rows;
    int columns = A.columns;
    for (int row = rows; --row >= 0; ) {
      for (int column = columns; --column >= 0; ) {
        if (row != column && !((A.get(row, column)).abs() <= epsilon)) return false;
      }
    }
    return true;
  }

  /**
   * A matrix <tt>A</tt> is <i>diagonally dominant by column</i> if the
   * absolute value of each diagonal element is larger than the sum of the
   * absolute values of the off-diagonal elements in the corresponding column.
   *
   * <tt>returns true if for all i: abs(A[i,i]) &gt; Sum(abs(A[j,i])); j != i.</tt>
   * Matrix may but need not be square.
   * <p>
   * Note: Ignores tolerance.
   */
  bool isDiagonallyDominantByColumn(DoubleMatrix A) {
    //cern.jet.math.tdouble.DoubleFunctions F = cern.jet.math.tdouble.DoubleFunctions.functions;
    int min = Math.min(A.rows, A.columns);
    for (int i = min; --i >= 0; ) {
      double diag = (A.get(i, i)).abs();
      diag += diag;
      if (diag <= A.column(i).reduce(plus, abs)) return false;
    }
    return true;
  }

  /**
   * A matrix <tt>A</tt> is <i>diagonally dominant by row</i> if the absolute
   * value of each diagonal element is larger than the sum of the absolute
   * values of the off-diagonal elements in the corresponding row.
   * <tt>returns true if for all i: abs(A[i,i]) &gt; Sum(abs(A[i,j])); j != i.</tt>
   * Matrix may but need not be square.
   * <p>
   * Note: Ignores tolerance.
   */
  bool isDiagonallyDominantByRow(DoubleMatrix A) {
    //cern.jet.math.tdouble.DoubleFunctions F = cern.jet.math.tdouble.DoubleFunctions.functions;
    int min = Math.min(A.rows, A.columns);
    for (int i = min; --i >= 0; ) {
      double diag = (A.get(i, i)).abs();
      diag += diag;
      if (diag <= A.row(i).reduce(plus, abs)) return false;
    }
    return true;
  }

  /**
   * A matrix <tt>A</tt> is an <i>identity</i> matrix if <tt>A[i,i] == 1</tt>
   * and all other cells are zero. Matrix may but need not be square.
   */
  bool isIdentity(DoubleMatrix A) {
    double epsilon = tolerance();
    int rows = A.rows;
    int columns = A.columns;
    for (int row = rows; --row >= 0; ) {
      for (int column = columns; --column >= 0; ) {
        double v = A.get(row, column);
        if (row == column) {
          if (!((1 - v).abs() < epsilon)) return false;
        } else if (!(v.abs() <= epsilon)) return false;
      }
    }
    return true;
  }

  /**
   * A matrix <tt>A</tt> is <i>lower bidiagonal</i> if <tt>A[i,j]==0</tt>
   * unless <tt>i==j || i==j+1</tt>. Matrix may but need not be square.
   */
  bool isLowerBidiagonal(DoubleMatrix A) {
    double epsilon = tolerance();
    int rows = A.rows;
    int columns = A.columns;
    for (int row = rows; --row >= 0; ) {
      for (int column = columns; --column >= 0; ) {
        if (!(row == column || row == column + 1)) {
          if (!((A.get(row, column)).abs() <= epsilon)) return false;
        }
      }
    }
    return true;
  }

  /**
   * A matrix <tt>A</tt> is <i>lower triangular</i> if <tt>A[i,j]==0</tt>
   * whenever <tt>i &lt; j</tt>. Matrix may but need not be square.
   */
  bool isLowerTriangular(DoubleMatrix A) {
    double epsilon = tolerance();
    int rows = A.rows;
    int columns = A.columns;
    for (int column = columns; --column >= 0; ) {
      for (int row = Math.min(column, rows); --row >= 0; ) {
        if (!((A.get(row, column)).abs() <= epsilon)) return false;
      }
    }
    return true;
  }

  /**
   * A matrix <tt>A</tt> is <i>non-negative</i> if <tt>A[i,j] &gt;= 0</tt>
   * holds for all cells.
   * <p>
   * Note: Ignores tolerance.
   */
  bool isNonNegative(DoubleMatrix A) {
    int rows = A.rows;
    int columns = A.columns;
    for (int row = rows; --row >= 0; ) {
      for (int column = columns; --column >= 0; ) {
        if (!(A.get(row, column) >= 0)) return false;
      }
    }
    return true;
  }

  /**
   * A square matrix <tt>A</tt> is <i>orthogonal</i> if
   * <tt>A*transpose(A) = I</tt>.
   *
   * @throws ArgumentError
   *             if <tt>!isSquare(A)</tt>.
   */
  bool isOrthogonal(DoubleMatrix A) {
    checkSquare(A);
    return equalsMatrix(A.multiply(A, null, 1.0, 0.0, false, true), DoubleFactory2D.dense.identity(A.rows));
  }

  /**
   * A matrix <tt>A</tt> is <i>positive</i> if <tt>A[i,j] &gt; 0</tt> holds
   * for all cells.
   * <p>
   * Note: Ignores tolerance.
   */
  bool isPositive(DoubleMatrix A) {
    int rows = A.rows;
    int columns = A.columns;
    for (int row = rows; --row >= 0; ) {
      for (int column = columns; --column >= 0; ) {
        if (!(A.get(row, column) > 0)) return false;
      }
    }
    return true;
  }

  /**
   * A matrix <tt>A</tt> is <i>singular</i> if it has no inverse, that is, iff
   * <tt>det(A)==0</tt>.
   */
  /*bool isSingular(DoubleMatrix A) {
    return !(DenseDoubleAlgebra.DEFAULT.det(A).abs() >= tolerance());
  }*/

  /**
   * A square matrix <tt>A</tt> is <i>skew-symmetric</i> if
   * <tt>A = -transpose(A)</tt>, that is <tt>A[i,j] == -A[j,i]</tt>.
   *
   * @throws ArgumentError
   *             if <tt>!isSquare(A)</tt>.
   */
  bool isSkewSymmetric(DoubleMatrix A) {
    checkSquare(A);
    double epsilon = tolerance();
    int rows = A.rows;
    for (int row = rows; --row >= 0; ) {
      for (int column = rows; --column >= 0; ) {
        if (!((A.get(row, column) + A.get(column, row)).abs() <= epsilon)) return false;
      }
    }
    return true;
  }

  /**
   * A matrix <tt>A</tt> is <i>square</i> if it has the same number of rows
   * and columns.
   */
  bool isSquare(DoubleMatrix A) {
    return A.rows == A.columns;
  }

  /**
   * A matrix <tt>A</tt> is <i>strictly lower triangular</i> if
   * <tt>A[i,j]==0</tt> whenever <tt>i &lt;= j</tt>. Matrix may but need not
   * be square.
   */
  bool isStrictlyLowerTriangular(DoubleMatrix A) {
    double epsilon = tolerance();
    int rows = A.rows;
    int columns = A.columns;
    for (int column = columns; --column >= 0; ) {
      for (int row = Math.min(rows, column + 1); --row >= 0; ) {
        if (!(A.get(row, column).abs() <= epsilon)) return false;
      }
    }
    return true;
  }

  /**
   * A matrix <tt>A</tt> is <i>strictly triangular</i> if it is triangular and
   * its diagonal elements all equal 0. Matrix may but need not be square.
   */
  bool isStrictlyTriangular(DoubleMatrix A) {
    if (!isTriangular(A)) return false;

    double epsilon = tolerance();
    for (int i = Math.min(A.rows, A.columns); --i >= 0; ) {
      if (!(A.get(i, i).abs() <= epsilon)) return false;
    }
    return true;
  }

  /**
   * A matrix <tt>A</tt> is <i>strictly upper triangular</i> if
   * <tt>A[i,j]==0</tt> whenever <tt>i &gt;= j</tt>. Matrix may but need not
   * be square.
   */
  bool isStrictlyUpperTriangular(DoubleMatrix A) {
    double epsilon = tolerance();
    int rows = A.rows;
    int columns = A.columns;
    for (int column = columns; --column >= 0; ) {
      for (int row = rows; --row >= column; ) {
        if (!(A.get(row, column).abs() <= epsilon)) return false;
      }
    }
    return true;
  }

  /**
   * A matrix <tt>A</tt> is <i>symmetric</i> if <tt>A = tranpose(A)</tt>, that
   * is <tt>A[i,j] == A[j,i]</tt>.
   *
   * @throws ArgumentError
   *             if <tt>!isSquare(A)</tt>.
   */
  bool isSymmetric(DoubleMatrix A) {
    checkSquare(A);
    return equalsMatrix(A, A.dice());
  }

  /**
   * A matrix <tt>A</tt> is <i>triangular</i> iff it is either upper or lower
   * triangular. Matrix may but need not be square.
   */
  bool isTriangular(DoubleMatrix A) {
    return isLowerTriangular(A) || isUpperTriangular(A);
  }

  /**
   * A matrix <tt>A</tt> is <i>tridiagonal</i> if <tt>A[i,j]==0</tt> whenever
   * <tt>Math.abs(i-j) > 1</tt>. Matrix may but need not be square.
   */
  bool isTridiagonal(DoubleMatrix A) {
    double epsilon = tolerance();
    int rows = A.rows;
    int columns = A.columns;
    for (int row = rows; --row >= 0; ) {
      for (int column = columns; --column >= 0; ) {
        if ((row - column).abs() > 1) {
          if (!(A.get(row, column).abs() <= epsilon)) return false;
        }
      }
    }
    return true;
  }

  /**
   * A matrix <tt>A</tt> is <i>unit triangular</i> if it is triangular and its
   * diagonal elements all equal 1. Matrix may but need not be square.
   */
  bool isUnitTriangular(DoubleMatrix A) {
    if (!isTriangular(A)) return false;

    double epsilon = tolerance();
    for (int i = Math.min(A.rows, A.columns); --i >= 0; ) {
      if (!((1 - A.get(i, i)).abs() <= epsilon)) return false;
    }
    return true;
  }

  /**
   * A matrix <tt>A</tt> is <i>upper bidiagonal</i> if <tt>A[i,j]==0</tt>
   * unless <tt>i==j || i==j-1</tt>. Matrix may but need not be square.
   */
  bool isUpperBidiagonal(DoubleMatrix A) {
    double epsilon = tolerance();
    int rows = A.rows;
    int columns = A.columns;
    for (int row = rows; --row >= 0; ) {
      for (int column = columns; --column >= 0; ) {
        if (!(row == column || row == column - 1)) {
          if (!(A.get(row, column).abs() <= epsilon)) return false;
        }
      }
    }
    return true;
  }

  /**
   * A matrix <tt>A</tt> is <i>upper triangular</i> if <tt>A[i,j]==0</tt>
   * whenever <tt>i &gt; j</tt>. Matrix may but need not be square.
   */
  bool isUpperTriangular(DoubleMatrix A) {
    double epsilon = tolerance();
    int rows = A.rows;
    int columns = A.columns;
    for (int column = columns; --column >= 0; ) {
      for (int row = rows; --row > column; ) {
        if (!(A.get(row, column).abs() <= epsilon)) return false;
      }
    }
    return true;
  }

  /**
   * A matrix <tt>A</tt> is <i>zero</i> if all its cells are zero.
   */
  bool isZero(DoubleMatrix A) {
    return equalsMatrixValue(A, 0.0);
  }

  /**
   * The <i>lower bandwidth</i> of a square matrix <tt>A</tt> is the maximum
   * <tt>i-j</tt> for which <tt>A[i,j]</tt> is nonzero and <tt>i &gt; j</tt>.
   * A <i>banded</i> matrix has a "band" about the diagonal. Diagonal,
   * tridiagonal and triangular matrices are special cases.
   *
   * @param A
   *            the square matrix to analyze.
   * @return the lower bandwith.
   * @throws ArgumentError
   *             if <tt>!isSquare(A)</tt>.
   * @see #semiBandwidth(DoubleMatrix)
   * @see #upperBandwidth(DoubleMatrix)
   */
  int lowerBandwidth(DoubleMatrix A) {
    checkSquare(A);
    double epsilon = tolerance();
    int rows = A.rows;

    for (int k = rows; --k >= 0; ) {
      for (int i = rows - k; --i >= 0; ) {
        int j = i + k;
        if (!(A.get(j, i).abs() <= epsilon)) return k;
      }
    }
    return 0;
  }

  /**
   * Returns the <i>semi-bandwidth</i> of the given square matrix <tt>A</tt>.
   * A <i>banded</i> matrix has a "band" about the diagonal. It is a matrix
   * with all cells equal to zero, with the possible exception of the cells
   * along the diagonal line, the <tt>k</tt> diagonal lines above the
   * diagonal, and the <tt>k</tt> diagonal lines below the diagonal. The
   * <i>semi-bandwith l</i> is the number <tt>k+1</tt>. The <i>bandwidth p</i>
   * is the number <tt>2*k + 1</tt>. For example, a tridiagonal matrix
   * corresponds to <tt>k=1, l=2, p=3</tt>, a diagonal or zero matrix
   * corresponds to <tt>k=0, l=1, p=1</tt>,
   * <p>
   * The <i>upper bandwidth</i> is the maximum <tt>j-i</tt> for which
   * <tt>A[i,j]</tt> is nonzero and <tt>j &gt; i</tt>. The <i>lower
   * bandwidth</i> is the maximum <tt>i-j</tt> for which <tt>A[i,j]</tt> is
   * nonzero and <tt>i &gt; j</tt>. Diagonal, tridiagonal and triangular
   * matrices are special cases.
   * <p>
   * Examples:
   * <table border="1" cellspacing="0">
   * <tr align="left" valign="top">
   * <td valign="middle" align="left"><tt>matrix</tt></td>
   * <td> <tt>4&nbsp;x&nbsp;4&nbsp;<br>
     0&nbsp;0&nbsp;0&nbsp;0<br>
     0&nbsp;0&nbsp;0&nbsp;0<br>
     0&nbsp;0&nbsp;0&nbsp;0<br>
     0&nbsp;0&nbsp;0&nbsp;0 </tt></td>
   * <td><tt>4&nbsp;x&nbsp;4<br>
     1&nbsp;0&nbsp;0&nbsp;0<br>
     0&nbsp;0&nbsp;0&nbsp;0<br>
     0&nbsp;0&nbsp;0&nbsp;0<br>
     0&nbsp;0&nbsp;0&nbsp;1 </tt></td>
   * <td><tt>4&nbsp;x&nbsp;4<br>
     1&nbsp;1&nbsp;0&nbsp;0<br>
     1&nbsp;1&nbsp;1&nbsp;0<br>
     0&nbsp;1&nbsp;1&nbsp;1<br>
     0&nbsp;0&nbsp;1&nbsp;1 </tt></td>
   * <td><tt> 4&nbsp;x&nbsp;4<br>
     0&nbsp;1&nbsp;1&nbsp;1<br>
     0&nbsp;1&nbsp;1&nbsp;1<br>
     0&nbsp;0&nbsp;0&nbsp;1<br>
     0&nbsp;0&nbsp;0&nbsp;1 </tt></td>
   * <td><tt> 4&nbsp;x&nbsp;4<br>
     0&nbsp;0&nbsp;0&nbsp;0<br>
     1&nbsp;1&nbsp;0&nbsp;0<br>
     1&nbsp;1&nbsp;0&nbsp;0<br>
     1&nbsp;1&nbsp;1&nbsp;1 </tt></td>
   * <td><tt>4&nbsp;x&nbsp;4<br>
     1&nbsp;1&nbsp;0&nbsp;0<br>
     0&nbsp;1&nbsp;1&nbsp;0<br>
     0&nbsp;1&nbsp;0&nbsp;1<br>
     1&nbsp;0&nbsp;1&nbsp;1 </tt><tt> </tt></td>
   * <td><tt>4&nbsp;x&nbsp;4<br>
     1&nbsp;1&nbsp;1&nbsp;0<br>
     0&nbsp;1&nbsp;0&nbsp;0<br>
     1&nbsp;1&nbsp;0&nbsp;1<br>
     0&nbsp;0&nbsp;1&nbsp;1 </tt></td>
   * </tr>
   * <tr align="center" valign="middle">
   * <td><tt>upperBandwidth</tt></td>
   * <td><div align="center"><tt>0</tt></div></td>
   * <td><div align="center"><tt>0</tt></div></td>
   * <td><div align="center"><tt>1</tt></div></td>
   * <td><tt>3</tt></td>
   * <td align="center" valign="middle"><tt>0</tt></td>
   * <td align="center" valign="middle"><div align="center"><tt>1</tt></div></td>
   * <td align="center" valign="middle"><div align="center"><tt>2</tt></div></td>
   * </tr>
   * <tr align="center" valign="middle">
   * <td><tt>lowerBandwidth</tt></td>
   * <td><div align="center"><tt>0</tt></div></td>
   * <td><div align="center"><tt>0</tt></div></td>
   * <td><div align="center"><tt>1</tt></div></td>
   * <td><tt>0</tt></td>
   * <td align="center" valign="middle"><tt>3</tt></td>
   * <td align="center" valign="middle"><div align="center"><tt>3</tt></div></td>
   * <td align="center" valign="middle"><div align="center"><tt>2</tt></div></td>
   * </tr>
   * <tr align="center" valign="middle">
   * <td><tt>semiBandwidth</tt></td>
   * <td><div align="center"><tt>1</tt></div></td>
   * <td><div align="center"><tt>1</tt></div></td>
   * <td><div align="center"><tt>2</tt></div></td>
   * <td><tt>4</tt></td>
   * <td align="center" valign="middle"><tt>4</tt></td>
   * <td align="center" valign="middle"><div align="center"><tt>4</tt></div></td>
   * <td align="center" valign="middle"><div align="center"><tt>3</tt></div></td>
   * </tr>
   * <tr align="center" valign="middle">
   * <td><tt>description</tt></td>
   * <td><div align="center"><tt>zero</tt></div></td>
   * <td><div align="center"><tt>diagonal</tt></div></td>
   * <td><div align="center"><tt>tridiagonal</tt></div></td>
   * <td><tt>upper triangular</tt></td>
   * <td align="center" valign="middle"><tt>lower triangular</tt></td>
   * <td align="center" valign="middle"><div align="center">
   * <tt>unstructured</tt></div></td>
   * <td align="center" valign="middle"><div align="center">
   * <tt>unstructured</tt></div></td>
   * </tr>
   * </table>
   *
   * @param A
   *            the square matrix to analyze.
   * @return the semi-bandwith <tt>l</tt>.
   * @throws ArgumentError
   *             if <tt>!isSquare(A)</tt>.
   * @see #lowerBandwidth(DoubleMatrix)
   * @see #upperBandwidth(DoubleMatrix)
   */
  int semiBandwidth(DoubleMatrix A) {
    checkSquare(A);
    double epsilon = tolerance();
    int rows = A.rows;

    for (int k = rows; --k >= 0; ) {
      for (int i = rows - k; --i >= 0; ) {
        int j = i + k;
        if (!(A.get(j, i).abs() <= epsilon)) return k + 1;
        if (!(A.get(i, j).abs() <= epsilon)) return k + 1;
      }
    }
    return 1;
  }

  /**
   * Sets the tolerance to <tt>Math.abs(newTolerance)</tt>.
   *
   * @throws UnsupportedOperationException
   *             if <tt>this==DEFAULT || this==ZERO || this==TWELVE</tt>.
   */
  void setTolerance(double newTolerance) {
    if (this == DEFAULT || this == ZERO || this == TWELVE) {
      throw new ArgumentError("Attempted to modify immutable object.");
    }
    _tolerance = newTolerance.abs();
  }

  /**
   * Returns the current tolerance.
   */
  double tolerance() {
    return _tolerance;
  }

  /**
   * Returns summary information about the given matrix <tt>A</tt>. That is a
   * String with (propertyName, propertyValue) pairs. Useful for debugging or
   * to quickly get the rough picture of a matrix. For example,
   *
   * <pre>
   * 	 density                      : 0.9
   * 	 isDiagonal                   : false
   * 	 isDiagonallyDominantByRow    : false
   * 	 isDiagonallyDominantByColumn : false
   * 	 isIdentity                   : false
   * 	 isLowerBidiagonal            : false
   * 	 isLowerTriangular            : false
   * 	 isNonNegative                : true
   * 	 isOrthogonal                 : Illegal operation or error: Matrix must be square.
   * 	 isPositive                   : true
   * 	 isSingular                   : Illegal operation or error: Matrix must be square.
   * 	 isSkewSymmetric              : Illegal operation or error: Matrix must be square.
   * 	 isSquare                     : false
   * 	 isStrictlyLowerTriangular    : false
   * 	 isStrictlyTriangular         : false
   * 	 isStrictlyUpperTriangular    : false
   * 	 isSymmetric                  : Illegal operation or error: Matrix must be square.
   * 	 isTriangular                 : false
   * 	 isTridiagonal                : false
   * 	 isUnitTriangular             : false
   * 	 isUpperBidiagonal            : false
   * 	 isUpperTriangular            : false
   * 	 isZero                       : false
   * 	 lowerBandwidth               : Illegal operation or error: Matrix must be square.
   * 	 semiBandwidth                : Illegal operation or error: Matrix must be square.
   * 	 upperBandwidth               : Illegal operation or error: Matrix must be square.
   *
   * </pre>
   */
  String toString2D(DoubleMatrix A) {
    final List names = new List();
    final List values = new List();
    String unknown = "Illegal operation or error: ";

    // determine properties
    names.add("density");
    try {
      values.add(density(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    // determine properties
    names.add("isDiagonal");
    try {
      values.add(isDiagonal(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    // determine properties
    names.add("isDiagonallyDominantByRow");
    try {
      values.add(isDiagonallyDominantByRow(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    // determine properties
    names.add("isDiagonallyDominantByColumn");
    try {
      values.add(isDiagonallyDominantByColumn(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    names.add("isIdentity");
    try {
      values.add(isIdentity(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    names.add("isLowerBidiagonal");
    try {
      values.add(isLowerBidiagonal(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    names.add("isLowerTriangular");
    try {
      values.add(isLowerTriangular(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    names.add("isNonNegative");
    try {
      values.add(isNonNegative(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    names.add("isOrthogonal");
    try {
      values.add(isOrthogonal(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    names.add("isPositive");
    try {
      values.add(isPositive(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    /*names.add("isSingular");
    try {
      values.add(isSingular(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }*/

    names.add("isSkewSymmetric");
    try {
      values.add(isSkewSymmetric(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    names.add("isSquare");
    try {
      values.add(isSquare(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    names.add("isStrictlyLowerTriangular");
    try {
      values.add(isStrictlyLowerTriangular(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    names.add("isStrictlyTriangular");
    try {
      values.add(isStrictlyTriangular(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    names.add("isStrictlyUpperTriangular");
    try {
      values.add(isStrictlyUpperTriangular(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    names.add("isSymmetric");
    try {
      values.add(isSymmetric(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    names.add("isTriangular");
    try {
      values.add(isTriangular(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    names.add("isTridiagonal");
    try {
      values.add(isTridiagonal(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    names.add("isUnitTriangular");
    try {
      values.add(isUnitTriangular(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    names.add("isUpperBidiagonal");
    try {
      values.add(isUpperBidiagonal(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    names.add("isUpperTriangular");
    try {
      values.add(isUpperTriangular(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    names.add("isZero");
    try {
      values.add(isZero(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    names.add("lowerBandwidth");
    try {
      values.add(lowerBandwidth(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    names.add("semiBandwidth");
    try {
      values.add(semiBandwidth(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    names.add("upperBandwidth");
    try {
      values.add(upperBandwidth(A).toString());
    } on ArgumentError catch (exc) {
      values.add(unknown + exc.message);
    }

    // sort ascending by property name
    comp(int a, int b) {
      return _get(names, a).compareTo(_get(names, b));
    }
    void swapper(int a, int b) {
      Object tmp;
      tmp = names[a];
      names[a] = names[b];
      names[b] = tmp;
      tmp = values[a];
      values[a] = values[b];
      values[b] = tmp;
    }
//    GenericSorting.quickSort(0, names.length, comp, swapper);

    // determine padding for nice formatting
    int maxLength = 0;
    for (int i = 0; i < names.length; i++) {
      int length = (names[i] as String).length;
      maxLength = Math.max(length, maxLength);
    }

    // finally, format properties
    StringBuffer buf = new StringBuffer();
    for (int i = 0; i < names.length; i++) {
      String name = (names[i] as String);
      buf.write(name);
      buf.write(_blanks(maxLength - name.length));
      buf.write(" : ");
      buf.write(values[i]);
      if (i < names.length - 1) buf.write('\n');
    }

    return buf.toString();
  }

  /**
   * The <i>upper bandwidth</i> of a square matrix <tt>A</tt> is the maximum
   * <tt>j-i</tt> for which <tt>A[i,j]</tt> is nonzero and <tt>j &gt; i</tt>.
   * A <i>banded</i> matrix has a "band" about the diagonal. Diagonal,
   * tridiagonal and triangular matrices are special cases.
   *
   * @param A
   *            the square matrix to analyze.
   * @return the upper bandwith.
   * @throws ArgumentError
   *             if <tt>!isSquare(A)</tt>.
   * @see #semiBandwidth(DoubleMatrix)
   * @see #lowerBandwidth(DoubleMatrix)
   */
  int upperBandwidth(DoubleMatrix A) {
    checkSquare(A);
    double epsilon = tolerance();
    int rows = A.rows;

    for (int k = rows; --k >= 0; ) {
      for (int i = rows - k; --i >= 0; ) {
        int j = i + k;
        if (!(A.get(i, j).abs() <= epsilon)) return k;
      }
    }
    return 0;
  }
}
