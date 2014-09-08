part of cern.colt.matrix.tint.algo;

/**
 * Checks whether the given matrix <tt>A</tt> is <i>rectangular</i>.
 *
 * @throws ArgumentError
 *             if <tt>A.rows < A.columns</tt>.
 */
void checkRectangular(AbstractIntMatrix A) {
  if (A.rows < A.columns) {
    throw new ArgumentError("Matrix must be rectangular: " + AbstractFormatter.shape2D(A));
  }
}

/**
 * Checks whether the given matrix <tt>A</tt> is <i>square</i>.
 *
 * @throws ArgumentError
 *             if <tt>A.rows != A.columns</tt>.
 */
void checkSquare(AbstractIntMatrix A) {
  if (A.rows != A.columns) {
    throw new ArgumentError("Matrix must be square: " + AbstractFormatter.shape2D(A));
  }
}

/**
 * Returns the matrix's fraction of non-zero cells;
 * <tt>A.cardinality() / A.size()</tt>.
 */
int density(AbstractIntMatrix A) {
  return A.cardinality() ~/ A.length;
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
bool equalsVectorValue(final AbstractIntVector A, final int value) {
  if (A == null) {
    return false;
  }
  int size = A.length;
  bool result = false;
  /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
  if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
    nthreads = Math.min(nthreads, size);
    List<Future> futures = new List<Future>(nthreads);
    List<bool> results = new List<bool>(nthreads);
    int k = size / nthreads;
    for (int j = 0; j < nthreads; j++) {
      final int firstIdx = j * k;
      final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
      futures[j] = ConcurrencyUtils.submit(() {
        for (int i = firstIdx; i < lastIdx; i++) {
          if (!(A.get(i) == value)) return false;
        }
        return true;
      });
    }
    try {
      for (int j = 0; j < nthreads; j++) {
        results[j] = futures[j].get();
      }
      result = results[0];
      for (int j = 1; j < nthreads; j++) {
        result = result && results[j];
      }
    } on ExecutionException catch (ex) {
      ex.printStackTrace();
    } on InterruptedException catch (e) {
      e.printStackTrace();
    }
    return result;
  } else {*/
  for (int i = 0; i < size; i++) {
    if (!(A.get(i) == value)) {
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
bool equalsVector(final AbstractIntVector A, final AbstractIntVector B) {
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

  bool result = false;
  /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
  if ((nthreads > 1) && (size >= ConcurrencyUtils.getThreadsBeginN_1D())) {
    nthreads = Math.min(nthreads, size);
    List<Future> futures = new List<Future>(nthreads);
    List<bool> results = new List<bool>(nthreads);
    int k = size / nthreads;
    for (int j = 0; j < nthreads; j++) {
      final int firstIdx = j * k;
      final int lastIdx = (j == nthreads - 1) ? size : firstIdx + k;
      futures[j] = ConcurrencyUtils.submit(() {
        for (int i = firstIdx; i < lastIdx; i++) {
          if (!(A.get(i) == B.get(i))) return false;
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
  for (int i = 0; i < size; i++) {
    if (!(A.get(i) == B.get(i))) return false;
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
bool equalsMatrixValue(final AbstractIntMatrix A, final int value) {
  if (A == null) return false;
  final int rows = A.rows;
  final int columns = A.columns;
  bool result = false;
  /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
  if ((nthreads > 1) && (A.length >= ConcurrencyUtils.getThreadsBeginN_2D())) {
    nthreads = Math.min(nthreads, A.rows);
    List<Future> futures = new List<Future>(nthreads);
    List<bool> results = new List<bool>(nthreads);
    int k = A.rows / nthreads;
    for (int j = 0; j < nthreads; j++) {
      final int firstRow = j * k;
      final int lastRow = (j == nthreads - 1) ? A.rows : firstRow + k;
      futures[j] = ConcurrencyUtils.submit(() {
        for (int r = firstRow; r < lastRow; r++) {
          for (int c = 0; c < columns; c++) {
            if (!(A.get(r, c) == value)) return false;
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
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < columns; c++) {
      if (!(A.get(r, c) == value)) {
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
bool equalsMatrix(final AbstractIntMatrix A, final AbstractIntMatrix B) {
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
  /*int nthreads = ConcurrencyUtils.getNumberOfThreads();
  if ((nthreads > 1) && (A.size() >= ConcurrencyUtils.getThreadsBeginN_2D())) {
    nthreads = Math.min(nthreads, A.rows);
    List<Future> futures = new List<Future>(nthreads);
    List<bool> results = new List<bool>(nthreads);
    int k = A.rows / nthreads;
    for (int j = 0; j < nthreads; j++) {
      final int firstRow = j * k;
      final int lastRow = (j == nthreads - 1) ? A.rows : firstRow + k;
      futures[j] = ConcurrencyUtils.submit(() {
        for (int r = firstRow; r < lastRow; r++) {
          for (int c = 0; c < columns; c++) {
            if (!(A.get(r, c) == B.get(r, c))) return false;
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
  for (int r = 0; r < rows; r++) {
    for (int c = 0; c < columns; c++) {
      if (!(A.get(r, c) == B.get(r, c))) {
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
/*bool equalsCubeValue(final IntMatrix3D A, final int value) {
  if (A == null) return false;
  final int slices = A.slices();
  final int rows = A.rows;
  final int columns = A.columns;
  bool result = false;
  int nthreads = ConcurrencyUtils.getNumberOfThreads();
  if ((nthreads > 1) && (A.size() >= ConcurrencyUtils.getThreadsBeginN_3D())) {
    nthreads = Math.min(nthreads, slices);
    List<Future> futures = new List<Future>(nthreads);
    List<bool> results = new List<bool>(nthreads);
    int k = slices / nthreads;
    for (int j = 0; j < nthreads; j++) {
      final int startslice = j * k;
      final int stopslice;
      if (j == nthreads - 1) {
        stopslice = slices;
      } else {
        stopslice = startslice + k;
      }
      futures[j] = ConcurrencyUtils.submit(() {
        for (int s = startslice; s < stopslice; s++) {
          for (int r = 0; r < rows; r++) {
            for (int c = 0; c < columns; c++) {
              if (!(A.getQuick(s, r, c) == value)) return false;
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
    for (int s = 0; s < slices; s++) {
      for (int r = 0; r < rows; r++) {
        for (int c = 0; c < columns; c++) {
          if (!(A.getQuick(s, r, c) == value)) return false;
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
/*bool equalsCube(final IntMatrix3D A, final IntMatrix3D B) {
  if (A == B) return true;
  if (!(A != null && B != null)) return false;
  final int slices = A.slices();
  final int rows = A.rows;
  final int columns = A.columns;
  if (columns != B.columns || rows != B.rows || slices != B.slices()) return false;
  bool result = false;
  int nthreads = ConcurrencyUtils.getNumberOfThreads();
  if ((nthreads > 1) && (A.size() >= ConcurrencyUtils.getThreadsBeginN_3D())) {
    nthreads = Math.min(nthreads, slices);
    List<Future> futures = new List<Future>(nthreads);
    List<bool> results = new List<bool>(nthreads);
    int k = slices / nthreads;
    for (int j = 0; j < nthreads; j++) {
      final int startslice = j * k;
      final int stopslice;
      if (j == nthreads - 1) {
        stopslice = slices;
      } else {
        stopslice = startslice + k;
      }
      futures[j] = ConcurrencyUtils.submit(() {
        for (int s = startslice; s < stopslice; s++) {
          for (int r = 0; r < rows; r++) {
            for (int c = 0; c < columns; c++) {
              if (!(A.getQuick(s, r, c) == B.getQuick(s, r, c))) return false;
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
    for (int s = 0; s < slices; s++) {
      for (int r = 0; r < rows; r++) {
        for (int c = 0; c < columns; c++) {
          if (!(A.getQuick(s, r, c) == B.getQuick(s, r, c))) return false;
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
void generateNonSingular(AbstractIntMatrix A) {
  checkSquare(A);
  int min = math.min(A.rows, A.columns);
  for (int i = min; --i >= 0; ) {
    A.set(i, i, 0);
  }
  for (int i = min; --i >= 0; ) {
    int rowSum = A.row(i).reduce(ifunc.plus, ifunc.abs);
    int colSum = A.column(i).reduce(ifunc.plus, ifunc.abs);
    A.set(i, i, math.max(rowSum, colSum) + i + 1);
  }
}

String _get(List list, int index) {
  return list[index];
}

/**
 * A matrix <tt>A</tt> is <i>diagonal</i> if <tt>A[i,j] == 0</tt> whenever
 * <tt>i != j</tt>. Matrix may but need not be square.
 */
bool isDiagonal(AbstractIntMatrix A) {
  int rows = A.rows;
  int columns = A.columns;
  for (int row = rows; --row >= 0; ) {
    for (int column = columns; --column >= 0; ) {
      if (row != column && A.get(row, column) != 0) return false;
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
bool isDiagonallyDominantByColumn(AbstractIntMatrix A) {
  int min = math.min(A.rows, A.columns);
  for (int i = min; --i >= 0; ) {
    int diag = A.get(i, i).abs();
    diag += diag;
    if (diag <= A.column(i).reduce(ifunc.plus, ifunc.abs)) return false;
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
bool isDiagonallyDominantByRow(AbstractIntMatrix A) {
  int min = math.min(A.rows, A.columns);
  for (int i = min; --i >= 0; ) {
    int diag = A.get(i, i).abs();
    diag += diag;
    if (diag <= A.row(i).reduce(ifunc.plus, ifunc.abs)) return false;
  }
  return true;
}

/**
 * A matrix <tt>A</tt> is an <i>identity</i> matrix if <tt>A[i,i] == 1</tt>
 * and all other cells are zero. Matrix may but need not be square.
 */
bool isIdentity(AbstractIntMatrix A) {
  int rows = A.rows;
  int columns = A.columns;
  for (int row = rows; --row >= 0; ) {
    for (int column = columns; --column >= 0; ) {
      int v = A.get(row, column);
      if (row == column) {
        if (v != 1) return false;
      } else if (v != 0) return false;
    }
  }
  return true;
}

/**
 * A matrix <tt>A</tt> is <i>lower bidiagonal</i> if <tt>A[i,j]==0</tt>
 * unless <tt>i==j || i==j+1</tt>. Matrix may but need not be square.
 */
bool isLowerBidiagonal(AbstractIntMatrix A) {
  int rows = A.rows;
  int columns = A.columns;
  for (int row = rows; --row >= 0; ) {
    for (int column = columns; --column >= 0; ) {
      if (!(row == column || row == column + 1)) {
        if (A.get(row, column) != 0) return false;
      }
    }
  }
  return true;
}

/**
 * A matrix <tt>A</tt> is <i>lower triangular</i> if <tt>A[i,j]==0</tt>
 * whenever <tt>i &lt; j</tt>. Matrix may but need not be square.
 */
bool isLowerTriangular(AbstractIntMatrix A) {
  int rows = A.rows;
  int columns = A.columns;
  for (int column = columns; --column >= 0; ) {
    for (int row = math.min(column, rows); --row >= 0; ) {
      if (A.get(row, column) != 0) return false;
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
bool isNonNegative(AbstractIntMatrix A) {
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
bool isOrthogonal(AbstractIntMatrix A) {
  checkSquare(A);
  return equalsMatrix(A.multiply(A, null, 1, 0, false, true),
      IntFactory2D.dense.identity(A.rows));
}

/**
 * A matrix <tt>A</tt> is <i>positive</i> if <tt>A[i,j] &gt; 0</tt> holds
 * for all cells.
 * <p>
 * Note: Ignores tolerance.
 */
bool isPositive(AbstractIntMatrix A) {
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
//    bool isSingular(IntMatrix A) {
//        return !(Math.abs(IntAlgebra.DEFAULT.det(A)) >= tolerance());
//    }
/**
 * A square matrix <tt>A</tt> is <i>skew-symmetric</i> if
 * <tt>A = -transpose(A)</tt>, that is <tt>A[i,j] == -A[j,i]</tt>.
 *
 * @throws ArgumentError
 *             if <tt>!isSquare(A)</tt>.
 */
bool isSkewSymmetric(AbstractIntMatrix A) {
  checkSquare(A);
  int rows = A.rows;
  int columns = A.columns;
  for (int row = rows; --row >= 0; ) {
    for (int column = rows; --column >= 0; ) {
      if (A.get(row, column) != -A.get(column, row)) return false;
    }
  }
  return true;
}

/**
 * A matrix <tt>A</tt> is <i>square</i> if it has the same number of rows
 * and columns.
 */
bool isSquare(AbstractIntMatrix A) {
  return A.rows == A.columns;
}

/**
 * A matrix <tt>A</tt> is <i>strictly lower triangular</i> if
 * <tt>A[i,j]==0</tt> whenever <tt>i &lt;= j</tt>. Matrix may but need not
 * be square.
 */
bool isStrictlyLowerTriangular(AbstractIntMatrix A) {
  int rows = A.rows;
  int columns = A.columns;
  for (int column = columns; --column >= 0; ) {
    for (int row = math.min(rows, column + 1); --row >= 0; ) {
      if (A.get(row, column) != 0) return false;
    }
  }
  return true;
}

/**
 * A matrix <tt>A</tt> is <i>strictly triangular</i> if it is triangular and
 * its diagonal elements all equal 0. Matrix may but need not be square.
 */
bool isStrictlyTriangular(AbstractIntMatrix A) {
  if (!isTriangular(A)) return false;

  for (int i = math.min(A.rows, A.columns); --i >= 0; ) {
    if (A.get(i, i) != 0) return false;
  }
  return true;
}

/**
 * A matrix <tt>A</tt> is <i>strictly upper triangular</i> if
 * <tt>A[i,j]==0</tt> whenever <tt>i &gt;= j</tt>. Matrix may but need not
 * be square.
 */
bool isStrictlyUpperTriangular(AbstractIntMatrix A) {
  int rows = A.rows;
  int columns = A.columns;
  for (int column = columns; --column >= 0; ) {
    for (int row = rows; --row >= column; ) {
      if (A.get(row, column) != 0) return false;
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
bool isSymmetric(AbstractIntMatrix A) {
  checkSquare(A);
  return equalsMatrix(A, A.dice());
}

/**
 * A matrix <tt>A</tt> is <i>triangular</i> iff it is either upper or lower
 * triangular. Matrix may but need not be square.
 */
bool isTriangular(AbstractIntMatrix A) {
  return isLowerTriangular(A) || isUpperTriangular(A);
}

/**
 * A matrix <tt>A</tt> is <i>tridiagonal</i> if <tt>A[i,j]==0</tt> whenever
 * <tt>Math.abs(i-j) > 1</tt>. Matrix may but need not be square.
 */
bool isTridiagonal(AbstractIntMatrix A) {
  int rows = A.rows;
  int columns = A.columns;
  for (int row = rows; --row >= 0; ) {
    for (int column = columns; --column >= 0; ) {
      if ((row - column).abs() > 1) {
        if (A.get(row, column) != 0) return false;
      }
    }
  }
  return true;
}

/**
 * A matrix <tt>A</tt> is <i>unit triangular</i> if it is triangular and its
 * diagonal elements all equal 1. Matrix may but need not be square.
 */
bool isUnitTriangular(AbstractIntMatrix A) {
  if (!isTriangular(A)) return false;

  for (int i = math.min(A.rows, A.columns); --i >= 0; ) {
    if (A.get(i, i) != 1) return false;
  }
  return true;
}

/**
 * A matrix <tt>A</tt> is <i>upper bidiagonal</i> if <tt>A[i,j]==0</tt>
 * unless <tt>i==j || i==j-1</tt>. Matrix may but need not be square.
 */
bool isUpperBidiagonal(AbstractIntMatrix A) {
  int rows = A.rows;
  int columns = A.columns;
  for (int row = rows; --row >= 0; ) {
    for (int column = columns; --column >= 0; ) {
      if (!(row == column || row == column - 1)) {
        if (A.get(row, column) != 0) return false;
      }
    }
  }
  return true;
}

/**
 * A matrix <tt>A</tt> is <i>upper triangular</i> if <tt>A[i,j]==0</tt>
 * whenever <tt>i &gt; j</tt>. Matrix may but need not be square.
 */
bool isUpperTriangular(AbstractIntMatrix A) {
  int rows = A.rows;
  int columns = A.columns;
  for (int column = columns; --column >= 0; ) {
    for (int row = rows; --row > column; ) {
      if (A.get(row, column) != 0) return false;
    }
  }
  return true;
}

/**
 * A matrix <tt>A</tt> is <i>zero</i> if all its cells are zero.
 */
bool isZero(AbstractIntMatrix A) {
  return equalsMatrixValue(A, 0);
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
 * @see #semiBandwidth(IntMatrix)
 * @see #upperBandwidth(IntMatrix)
 */
int lowerBandwidth(AbstractIntMatrix A) {
  checkSquare(A);
  int rows = A.rows;

  for (int k = rows; --k >= 0; ) {
    for (int i = rows - k; --i >= 0; ) {
      int j = i + k;
      if (A.get(j, i) != 0) return k;
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
 * @see #lowerBandwidth(IntMatrix)
 * @see #upperBandwidth(IntMatrix)
 */
int semiBandwidth(AbstractIntMatrix A) {
  checkSquare(A);
  int rows = A.rows;

  for (int k = rows; --k >= 0; ) {
    for (int i = rows - k; --i >= 0; ) {
      int j = i + k;
      if (A.get(j, i) != 0) return k + 1;
      if (A.get(i, j) != 0) return k + 1;
    }
  }
  return 1;
}

/**
 * Returns summary information about the given matrix <tt>A</tt>. That is a
 * String with (propertyName, propertyValue) pairs. Useful for debugging or
 * to quickly get the rough picture of a matrix. For example,
 *
 * <pre>
 *   density                      : 0.9
 *   isDiagonal                   : false
 *   isDiagonallyDominantByRow    : false
 *   isDiagonallyDominantByColumn : false
 *   isIdentity                   : false
 *   isLowerBidiagonal            : false
 *   isLowerTriangular            : false
 *   isNonNegative                : true
 *   isOrthogonal                 : Illegal operation or error: Matrix must be square.
 *   isPositive                   : true
 *   isSingular                   : Illegal operation or error: Matrix must be square.
 *   isSkewSymmetric              : Illegal operation or error: Matrix must be square.
 *   isSquare                     : false
 *   isStrictlyLowerTriangular    : false
 *   isStrictlyTriangular         : false
 *   isStrictlyUpperTriangular    : false
 *   isSymmetric                  : Illegal operation or error: Matrix must be square.
 *   isTriangular                 : false
 *   isTridiagonal                : false
 *   isUnitTriangular             : false
 *   isUpperBidiagonal            : false
 *   isUpperTriangular            : false
 *   isZero                       : false
 *   lowerBandwidth               : Illegal operation or error: Matrix must be square.
 *   semiBandwidth                : Illegal operation or error: Matrix must be square.
 *   upperBandwidth               : Illegal operation or error: Matrix must be square.
 *
 * </pre>
 */
String toString(AbstractIntMatrix A) {
  final names = new List();
  final values = new List();
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

  //        names.add("isSingular");
  //        try {
  //            values.add(isSingular(A)));
  //        } on ArgumentError catch (exc) {
  //            values.add(unknown + exc.message);
  //        }

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

  // TODO sort ascending by property name
  /*int comp(int a, int b) {
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
  GenericSorting.quickSort(0, names.length, comp, swapper);*/

  // determine padding for nice formatting
  int maxLength = 0;
  for (int i = 0; i < names.length; i++) {
    int length = (names[i]).length;
    maxLength = math.max(length, maxLength);
  }

  // finally, format properties
  StringBuffer buf = new StringBuffer();
  for (int i = 0; i < names.length; i++) {
    String name = names[i];
    buf.write(name);
    buf.write(_blanks(maxLength - name.length));
    buf.write(" : ");
    buf.write(values[i]);
    if (i < names.length - 1) {
      buf.write('\n');
    }
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
 * @see #semiBandwidth(IntMatrix)
 * @see #lowerBandwidth(IntMatrix)
 */
int upperBandwidth(AbstractIntMatrix A) {
  checkSquare(A);
  int rows = A.rows;

  for (int k = rows; --k >= 0; ) {
    for (int i = rows - k; --i >= 0; ) {
      int j = i + k;
      if (A.get(i, j) != 0) return k;
    }
  }
  return 0;
}

/** Returns a String with <tt>length</tt> blanks. */
String _blanks(int length) {
  if (length < 0) length = 0;
  StringBuffer buf = new StringBuffer(length);
  for (int k = 0; k < length; k++) {
    buf.write(' ');
  }
  return buf.toString();
}
