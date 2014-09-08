/*
Copyright (C) 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
is hereby granted without fee, provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear in supporting documentation.
CERN makes no representations about the suitability of this software for any purpose.
It is provided "as is" without expressed or implied warranty.
 */
part of cern.colt.matrix.complex.algo;

/**
 * Tests matrices for equality.
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
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.1, 28/May/2000 (fixed strange bugs involving NaN, -inf, inf)
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
class ComplexProperty {

  /**
   * The default Property object; currently has <tt>tolerance()==1.0E-9</tt>.
   */
  static final ComplexProperty DEFAULT = new ComplexProperty(1.0E-9);

  /**
   * A Property object with <tt>tolerance()==0.0</tt>.
   */
  static final ComplexProperty ZERO = new ComplexProperty(0.0);

  /**
   * A Property object with <tt>tolerance()==1.0E-12</tt>.
   */
  static final ComplexProperty TWELVE = new ComplexProperty(1.0E-12);

  double _tolerance;

  /**
   * Constructs an instance with a tolerance of
   * <tt>Math.abs(newTolerance)</tt>.
   */
  ComplexProperty(double newTolerance) {
    _tolerance = newTolerance.abs();
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

  void checkDense(AbstractComplexVector A) {
    if (!(A is ComplexVector)) {
      throw new ArgumentError("Matrix must be dense");
    }
  }

  /**
   * Checks whether the given matrix <tt>A</tt> is <i>square</i>.
   *
   * @throws ArgumentError
   *             if <tt>A.rows() != A.columns()</tt>.
   */
  void checkSquare(AbstractComplexMatrix A) {
    if (A.rows != A.columns) {
      throw new ArgumentError("Matrix must be square: " + AbstractFormatter.shape2D(A));
    }
  }

  void checkSparse(AbstractComplexMatrix A) {
    if (!(A is SparseCCComplexMatrix) && !(A is SparseRCComplexMatrix)) {
      throw new ArgumentError("Matrix must be sparse");
    }
  }

  /**
   * Returns whether all cells of the given matrix <tt>A</tt> are equal to the
   * given value.
   *
   * @param A
   *            the first matrix to compare.
   * @param value
   *            the value to compare against.
   * @return <tt>true</tt> if the matrix is equal to the value; <tt>false</tt>
   *         otherwise.
   */
  bool equalsValue1D(final AbstractComplexVector A, final Float64List value) {
    if (A == null) {
      return false;
    }
    final double epsilon = tolerance();
    bool result = false;
    int size = A.length;
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
      nthreads = Math.min(nthreads, size);
      List<Future> futures = new List<Future>(nthreads);
      List<Boolean> results = new List<Boolean>(nthreads);
      int k = size / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          Float64List diff = new Float64List(2);
          for (int i = firstIdx; i < lastIdx; i++) {
            Float64List x = A.getQuick(i);
            diff[0] = Math.abs(value[0] - x[0]);
            diff[1] = Math.abs(value[1] - x[1]);
            if (((diff[0] != diff[0]) || (diff[1] != diff[1])) && ((((value[0] != value[0]) || (value[1] != value[1])) && ((x[0] != x[0]) || (x[1] != x[1])))) || (Complex.isEqual(value, x, epsilon))) {
              diff[0] = 0;
              diff[1] = 0;
            }
            if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
              return false;
            }
          }
          return true;
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          results[j] = futures[j].get();
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
      Float64List diff = new Float64List(2);
      for (int i = 0; i < A.length; i++) {
        Float64List x = A.get(i);
        diff[0] = (value[0] - x[0]).abs();
        diff[1] = (value[1] - x[1]).abs();
        if (((diff[0] != diff[0]) || (diff[1] != diff[1])) && ((((value[0] != value[0]) || (value[1] != value[1])) && ((x[0] != x[0]) || (x[1] != x[1])))) || (Complex.isEqual(value, x, epsilon))) {
          diff[0] = 0.0;
          diff[1] = 0.0;
        }
        if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
          return false;
        }
      }
      return true;
    //}
  }

  /**
   * Returns whether both given matrices <tt>A</tt> and <tt>B</tt> are equal.
   *
   * @param A
   *            the first matrix to compare.
   * @param B
   *            the second matrix to compare.
   * @return <tt>true</tt> if both matrices are equal; <tt>false</tt>
   *         otherwise.
   */
  bool equalsVector(final AbstractComplexVector A, final AbstractComplexVector B) {
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
      List<Boolean> results = new List<Boolean>(nthreads);
      int k = size / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstIdx = j * k;
        final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
        futures[j] = ConcurrencyUtils.submit(() {
          Float64List diff = new Float64List(2);
          for (int i = firstIdx; i < lastIdx; i++) {
            Float64List x = A.getQuick(i);
            Float64List value = B.getQuick(i);
            diff[0] = Math.abs(value[0] - x[0]);
            diff[1] = Math.abs(value[1] - x[1]);
            if (((diff[0] != diff[0]) || (diff[1] != diff[1])) && ((((value[0] != value[0]) || (value[1] != value[1])) && ((x[0] != x[0]) || (x[1] != x[1])))) || (Complex.isEqual(value, x, epsilon))) {
              diff[0] = 0;
              diff[1] = 0;
            }
            if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
              return false;
            }
          }
          return true;
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          results[j] = futures[j].get();
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
      Float64List diff = new Float64List(2);
      for (int i = 0; i < size; i++) {
        Float64List x = A.get(i);
        Float64List value = B.get(i);
        diff[0] = (value[0] - x[0]).abs();
        diff[1] = (value[1] - x[1]).abs();
        if (((diff[0] != diff[0]) || (diff[1] != diff[1])) && ((((value[0] != value[0]) || (value[1] != value[1])) && ((x[0] != x[0]) || (x[1] != x[1])))) || (Complex.isEqual(value, x, epsilon))) {
          diff[0] = 0.0;
          diff[1] = 0.0;
        }
        if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
          return false;
        }
      }
      return true;
    //}
  }

  /**
   * Returns whether all cells of the given matrix <tt>A</tt> are equal to the
   * given value.
   *
   * @param A
   *            the first matrix to compare.
   * @param value
   *            the value to compare against.
   * @return <tt>true</tt> if the matrix is equal to the value; <tt>false</tt>
   *         otherwise.
   */
  bool equalsValue2D(final AbstractComplexMatrix A, final Float64List value) {
    if (A == null) {
      return false;
    }
    int rows = A.rows;
    int columns = A.columns;
    bool result = false;
    final double epsilon = tolerance();
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (A.size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, A.rows());
      List<Future> futures = new List<Future>(nthreads);
      List<Boolean> results = new List<Boolean>(nthreads);
      int k = A.rows() / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? A.rows() : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          Float64List diff = new Float64List(2);
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = 0; c < A.columns(); c++) {
              Float64List x = A.getQuick(r, c);
              diff[0] = Math.abs(value[0] - x[0]);
              diff[1] = Math.abs(value[1] - x[1]);
              if (((diff[0] != diff[0]) || (diff[1] != diff[1])) && ((((value[0] != value[0]) || (value[1] != value[1])) && ((x[0] != x[0]) || (x[1] != x[1])))) || (Complex.isEqual(value, x, epsilon))) {
                diff[0] = 0;
                diff[1] = 0;
              }
              if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
                return false;
              }
            }
          }
          return true;
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          results[j] = futures[j].get();
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
      Float64List diff = new Float64List(2);
      for (int r = 0; r < rows; r++) {
        for (int c = 0; c < columns; c++) {
          Float64List x = A.get(r, c);
          diff[0] = (value[0] - x[0]).abs();
          diff[1] = (value[1] - x[1]).abs();
          if (((diff[0] != diff[0]) || (diff[1] != diff[1])) && ((((value[0] != value[0]) || (value[1] != value[1])) && ((x[0] != x[0]) || (x[1] != x[1])))) || (Complex.isEqual(value, x, epsilon))) {
            diff[0] = 0.0;
            diff[1] = 0.0;
          }
          if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
            return false;
          }
        }
      }
      return true;
    //}
  }

  /**
   * Returns whether both given matrices <tt>A</tt> and <tt>B</tt> are equal.
   *
   * @param A
   *            the first matrix to compare.
   * @param B
   *            the second matrix to compare.
   * @return <tt>true</tt> if both matrices are equal; <tt>false</tt>
   *         otherwise.
   */
  bool equalsMatrix(final AbstractComplexMatrix A, final AbstractComplexMatrix B) {
    if (identical(A, B)) {
      return true;
    }
    if (!(A != null && B != null)) {
      return false;
    }
    int rows = A.rows;
    int columns = A.columns;
    if (columns != B.columns || rows != B.rows) {
      return false;
    }
    bool result = false;
    final double epsilon = tolerance();
    /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (A.size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
      nthreads = Math.min(nthreads, A.rows());
      List<Future> futures = new List<Future>(nthreads);
      List<Boolean> results = new List<Boolean>(nthreads);
      int k = A.rows() / nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstRow = j * k;
        final int lastRow = (j == nthreads - 1) ? A.rows() : firstRow + k;
        futures[j] = ConcurrencyUtils.submit(() {
          Float64List diff = new Float64List(2);
          for (int r = firstRow; r < lastRow; r++) {
            for (int c = 0; c < A.columns(); c++) {
              Float64List x = A.getQuick(r, c);
              Float64List value = B.getQuick(r, c);
              diff[0] = Math.abs(value[0] - x[0]);
              diff[1] = Math.abs(value[1] - x[1]);
              if (((diff[0] != diff[0]) || (diff[1] != diff[1])) && ((((value[0] != value[0]) || (value[1] != value[1])) && ((x[0] != x[0]) || (x[1] != x[1])))) || (Complex.isEqual(value, x, epsilon))) {
                diff[0] = 0;
                diff[1] = 0;
              }
              if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
                return false;
              }
            }
          }
          return true;
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          results[j] = futures[j].get();
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
      Float64List diff = new Float64List(2);
      for (int r = 0; r < rows; r++) {
        for (int c = 0; c < columns; c++) {
          Float64List x = A.get(r, c);
          Float64List value = B.get(r, c);
          diff[0] = (value[0] - x[0]).abs();
          diff[1] = (value[1] - x[1]).abs();
          if (((diff[0] != diff[0]) || (diff[1] != diff[1])) && ((((value[0] != value[0]) || (value[1] != value[1])) && ((x[0] != x[0]) || (x[1] != x[1])))) || (Complex.isEqual(value, x, epsilon))) {
            diff[0] = 0.0;
            diff[1] = 0.0;
          }
          if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
            return false;
          }
        }
      }
      return true;
    //}
  }

  /**
   * Returns whether all cells of the given matrix <tt>A</tt> are equal to the
   * given value.
   *
   * @param A
   *            the first matrix to compare.
   * @param value
   *            the value to compare against.
   * @return <tt>true</tt> if the matrix is equal to the value; <tt>false</tt>
   *         otherwise.
   */
  /*bool equalsValue3D(final ComplexMatrix3D A, final Float64List value) {
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
      List<Boolean> results = new List<Boolean>(nthreads);
      int k = slices ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int firstSlice = j * k;
        final int lastSlice = (j == nthreads - 1) ? slices : firstSlice + k;
        futures[j] = ConcurrencyUtils.submit(() {
          Float64List diff = new Float64List(2);
          for (int s = firstSlice; s < lastSlice; s++) {
            for (int r = 0; r < rows; r++) {
              for (int c = 0; c < columns; c++) {
                Float64List x = A.getQuick(s, r, c);
                diff[0] = Math.abs(value[0] - x[0]);
                diff[1] = Math.abs(value[1] - x[1]);
                if (((diff[0] != diff[0]) || (diff[1] != diff[1])) && ((((value[0] != value[0]) || (value[1] != value[1])) && ((x[0] != x[0]) || (x[1] != x[1])))) || (Complex.isEqual(value, x, epsilon))) {
                  diff[0] = 0;
                  diff[1] = 0;
                }
                if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
                  return false;
                }
              }
            }
          }
          return true;
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          results[j] = futures[j].get();
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
      Float64List diff = new Float64List(2);
      for (int s = 0; s < slices; s++) {
        for (int r = 0; r < rows; r++) {
          for (int c = 0; c < columns; c++) {
            Float64List x = A.getQuick(s, r, c);
            diff[0] = Math.abs(value[0] - x[0]);
            diff[1] = Math.abs(value[1] - x[1]);
            if (((diff[0] != diff[0]) || (diff[1] != diff[1])) && ((((value[0] != value[0]) || (value[1] != value[1])) && ((x[0] != x[0]) || (x[1] != x[1])))) || (Complex.isEqual(value, x, epsilon))) {
              diff[0] = 0;
              diff[1] = 0;
            }
            if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
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
   *
   * @param A
   *            the first matrix to compare.
   * @param B
   *            the second matrix to compare.
   * @return <tt>true</tt> if both matrices are equal; <tt>false</tt>
   *         otherwise.
   */
  /*bool equalsMatrix3D(final ComplexMatrix3D A, final ComplexMatrix3D B) {
    if (A == B) return true;
    if (!(A != null && B != null)) return false;
    bool result = false;
    final int slices = A.slices();
    final int rows = A.rows();
    final int columns = A.columns();
    if (columns != B.columns() || rows != B.rows() || slices != B.slices()) return false;
    final double epsilon = tolerance();
    int nthreads = ConcurrencyUtils.getNumberOfThreads();
    if ((nthreads > 1) && (A.size() >= ConcurrencyUtils.getThreadsBeginN_3D())) {
      nthreads = Math.min(nthreads, slices);
      List<Future> futures = new List<Future>(nthreads);
      List<Boolean> results = new List<Boolean>(nthreads);
      int k = slices ~/ nthreads;
      for (int j = 0; j < nthreads; j++) {
        final int startslice = j * k;
        final int stopslice;
        if (j == nthreads - 1) {
          stopslice = slices;
        } else {
          stopslice = startslice + k;
        }
        futures[j] = ConcurrencyUtils.submit(() {
          Float64List diff = new Float64List(2);
          for (int s = startslice; s < stopslice; s++) {
            for (int r = 0; r < rows; r++) {
              for (int c = 0; c < columns; c++) {
                Float64List x = A.getQuick(s, r, c);
                Float64List value = B.getQuick(s, r, c);
                diff[0] = Math.abs(value[0] - x[0]);
                diff[1] = Math.abs(value[1] - x[1]);
                if (((diff[0] != diff[0]) || (diff[1] != diff[1])) && ((((value[0] != value[0]) || (value[1] != value[1])) && ((x[0] != x[0]) || (x[1] != x[1])))) || (Complex.isEqual(value, x, epsilon))) {
                  diff[0] = 0;
                  diff[1] = 0;
                }
                if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
                  return false;
                }
              }
            }
          }
          return true;
        });
      }
      try {
        for (int j = 0; j < nthreads; j++) {
          results[j] = futures[j].get();
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
      Float64List diff = new Float64List(2);
      for (int s = 0; s < slices; s++) {
        for (int r = 0; r < rows; r++) {
          for (int c = 0; c < columns; c++) {
            Float64List x = A.getQuick(s, r, c);
            Float64List value = B.getQuick(s, r, c);
            diff[0] = Math.abs(value[0] - x[0]);
            diff[1] = Math.abs(value[1] - x[1]);
            if (((diff[0] != diff[0]) || (diff[1] != diff[1])) && ((((value[0] != value[0]) || (value[1] != value[1])) && ((x[0] != x[0]) || (x[1] != x[1])))) || (Complex.isEqual(value, x, epsilon))) {
              diff[0] = 0;
              diff[1] = 0;
            }
            if ((diff[0] > epsilon) || (diff[1] > epsilon)) {
              return false;
            }
          }
        }
      }
      return true;
    }
  }*/

}
