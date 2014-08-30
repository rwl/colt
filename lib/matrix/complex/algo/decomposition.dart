library cern.colt.matrix.complex.algo.decomposition;

import 'dart:typed_data';
import 'package:csparse/complex/csparse.dart';
import 'package:klu/complex.dart' as klu;
import 'package:btf/btf.dart' as btf;

import '../../../function/complex.dart' show multiply;
import '../../matrix.dart';
import 'algo.dart';

/**
 * For a square matrix <tt>A</tt>, the LU decomposition is an unit lower
 * triangular matrix <tt>L</tt>, an upper triangular matrix <tt>U</tt>, and a
 * permutation vector <tt>piv</tt> so that <tt>A(piv,:) = L*U</tt>
 * <P>
 * The LU decomposition with pivoting always exists, even if the matrix is
 * singular. The primary use of the LU decomposition is in the solution of
 * square systems of simultaneous linear equations. This will fail if
 * <tt>isNonsingular()</tt> returns false.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
class SparseDComplexLUDecomposition {
  DZcss S;
  DZcsn N;
  DComplexMatrix2D L;
  DComplexMatrix2D U;
  bool rcMatrix = false;
  bool isNonSingular = true;
  /**
   * Row and column dimension (square matrix).
   */
  int n;

  /**
   * Constructs and returns a new LU Decomposition object; The decomposed
   * matrices can be retrieved via instance methods of the returned
   * decomposition object.
   *
   * @param A
   *            Square matrix
   * @param order
   *            ordering option (0 to 3); 0: natural ordering, 1: amd(A+A'),
   *            2: amd(S'*S), 3: amd(A'*A)
   * @param checkIfSingular
   *            if true, then the singularity test (based on
   *            Dulmage-Mendelsohn decomposition) is performed.
   * @throws ArgumentError
   *             if <tt>A</tt> is not square or is not sparse.
   * @throws ArgumentError
   *             if <tt>order</tt> is not in [0,3]
   */
  SparseDComplexLUDecomposition(DComplexMatrix2D A, int order, bool checkIfSingular) {
    DComplexProperty.DEFAULT.checkSquare(A);
    DComplexProperty.DEFAULT.checkSparse(A);

    if (order < 0 || order > 3) {
      throw new ArgumentError("order must be a number between 0 and 3");
    }
    DZcs dcs;
    if (A is SparseRCDComplexMatrix2D) {
      rcMatrix = true;
      dcs = (A).getColumnCompressed().elements();
    } else {
      dcs = A.elements() as DZcs;
    }
    n = A.rows();

    S = DZcs_sqr.cs_sqr(order, dcs, false);
    if (S == null) {
      throw new ArgumentError("Exception occured in cs_sqr()");
    }
    N = DZcs_lu.cs_lu(dcs, S, 1);
    if (N == null) {
      throw new ArgumentError("Exception occured in cs_lu()");
    }
    if (checkIfSingular) {
      DZcsd D = DZcs_dmperm.cs_dmperm(dcs, 1);
      /* check if matrix is singular */
      if (D != null && D.rr[3] < n) {
        isNonSingular = false;
      }
    }
  }

  /**
   * Returns the determinant, <tt>det(A)</tt>.
   *
   */
  Float64List det() {
    if (!isNonsingular()) return new Float64List.fromList([0, 0]); // avoid rounding errors
    int pivsign = 1;
    for (int i = 0; i < n; i++) {
      if (N.pinv[i] != i) {
        pivsign = -pivsign;
      }
    }
    if (U == null) {
      U = new SparseCCDComplexMatrix2D(N.U);
      if (rcMatrix) {
        U = (U as SparseCCDComplexMatrix2D).getRowCompressed();
      }
    }
    Float64List det = new Float64List.fromList([pivsign, 0]);
    for (int j = 0; j < n; j++) {
      det = multiply(det)(U.getQuick(j, j));
    }
    return det;
  }

  /**
   * Returns the lower triangular factor, <tt>L</tt>.
   *
   * @return <tt>L</tt>
   */
  DComplexMatrix2D getL() {
    if (L == null) {
      L = new SparseCCDComplexMatrix2D(N.L);
      if (rcMatrix) {
        L = (L as SparseCCDComplexMatrix2D).getRowCompressed();
      }
    }
    return L.copy();
  }

  /**
   * Returns a copy of the pivot permutation vector.
   *
   * @return piv
   */
  Int32List getPivot() {
    if (N.pinv == null) {
      return null;
    }
    Int32List pinv = new Int32List(N.pinv.length);
    //System.arraycopy(N.pinv, 0, pinv, 0, pinv.length);
    pinv.setAll(0, N.pinv);
    return pinv;
  }

  /**
   * Returns the upper triangular factor, <tt>U</tt>.
   *
   * @return <tt>U</tt>
   */
  DComplexMatrix2D getU() {
    if (U == null) {
      U = new SparseCCDComplexMatrix2D(N.U);
      if (rcMatrix) {
        U = (U as SparseCCDComplexMatrix2D).getRowCompressed();
      }
    }
    return U.copy();
  }

  /**
   * Returns a copy of the symbolic LU analysis object
   *
   * @return symbolic LU analysis
   */
  DZcss getSymbolicAnalysis() {
    DZcss S2 = new DZcss();
    S2.cp = S.cp != null ? S.cp.clone() : null;
    S2.leftmost = S.leftmost != null ? S.leftmost.clone() : null;
    S2.lnz = S.lnz;
    S2.m2 = S.m2;
    S2.parent = S.parent != null ? S.parent.clone() : null;
    S2.pinv = S.pinv != null ? S.pinv.clone() : null;
    S2.q = S.q != null ? S.q.clone() : null;
    S2.unz = S.unz;
    return S2;
  }

  /**
   * Returns whether the matrix is nonsingular (has an inverse).
   *
   * @return true if <tt>U</tt>, and hence <tt>A</tt>, is nonsingular; false
   *         otherwise.
   */
  bool isNonsingular() {
    return isNonSingular;
  }

  /**
   * Solves <tt>A*x = b</tt>(in-place). Upon return <tt>b</tt> is overridden
   * with the result <tt>x</tt>.
   *
   * @param b
   *            A vector with of size A.rows();
   * @exception ArgumentError
   *                if <tt>b.size() != A.rows()</tt> or if A is singular.
   */
  void solve(DComplexMatrix1D b) {
    if (b.size() != n) {
      throw new ArgumentError("b.size() != A.rows()");
    }
    if (!isNonsingular()) {
      throw new ArgumentError("A is singular");
    }
    DComplexProperty.DEFAULT.checkDense(b);
    DZcsa y = new DZcsa(n);
    DZcsa x;
    if (b.isView()) {
      x = new DZcsa(b.copy().elements() as Float64List);
    } else {
      x = new DZcsa(b.elements() as Float64List);
    }
    cs_ipvec(N.pinv, x, y, n);
    /* y = b(p) */
    cs_lsolve(N.L, y);
    /* y = L\y */
    cs_usolve(N.U, y);
    /* y = U\y */
    cs_ipvec(S.q, y, x, n);
    /* b(q) = x */

    if (b.isView()) {
      b.assignValues(x.x);
    }
  }
}
