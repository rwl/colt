library cern.colt.matrix.double.algo.decomposition;

import 'dart:typed_data';
//import 'package:csparse/double/csparse.dart';
//import 'package:klu/double.dart' as klu;
//import 'package:btf/btf.dart' as btf;

import '../../matrix.dart';
import 'algo.dart';

typedef void solve(AbstractDoubleMatrix A, AbstractDoubleVector b);

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
abstract class SparseDoubleLUDecomposition {

  /**
   * Returns the determinant, <tt>det(A)</tt>.
   *
   */
  double det();

  /**
   * Returns the lower triangular factor, <tt>L</tt>.
   *
   * @return <tt>L</tt>
   */
  AbstractDoubleMatrix getL();

  /**
   * Returns a copy of the pivot permutation vector.
   *
   * @return piv
   */
  Int32List getPivot();

  /**
   * Returns the upper triangular factor, <tt>U</tt>.
   *
   * @return <tt>U</tt>
   */
  AbstractDoubleMatrix getU();

  /**
   * Returns a copy of the symbolic LU analysis object
   *
   * @return symbolic LU analysis
   */
  Object getSymbolicAnalysis();

  /**
   * Returns whether the matrix is nonsingular (has an inverse).
   *
   * @return true if <tt>U</tt>, and hence <tt>A</tt>, is nonsingular; false
   *         otherwise.
   */
  bool isNonsingular();

  /**
   * Solves <tt>A*x = b</tt>(in-place). Upon return <tt>b</tt> is overridden
   * with the result <tt>x</tt>.
   *
   * @param b
   *            A vector with of size A.rows();
   * @exception ArgumentError
   *                if <tt>b.size() != A.rows()</tt> or if A is singular.
   */
  void solve(AbstractDoubleVector b);

}

/**
 * LU decomposition implemented using CSparseJ.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
class CSparseDoubleLUDecomposition implements SparseDoubleLUDecomposition {
  Symbolic _S;
  Numeric _N;
  AbstractDoubleMatrix _L;
  AbstractDoubleMatrix _U;
  bool _rcMatrix = false;
  bool _isNonSingular = true;
  /**
     * Row and column dimension (square matrix).
     */
  int _n;

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
  CSparseDoubleLUDecomposition(AbstractDoubleMatrix A, int order, bool checkIfSingular) {
    DoubleProperty.DEFAULT.checkSquare(A);
    DoubleProperty.DEFAULT.checkSparse2D(A);

    if (order < 0 || order > 3) {
      throw new ArgumentError("order must be a number between 0 and 3");
    }
    Matrix dcs;
    if (A is SparseRCDoubleMatrix) {
      _rcMatrix = true;
      dcs = A.getColumnCompressed().elements;
    } else {
      dcs = A.elements as Matrix;
    }
    _n = A.rows;

    _S = sqr(order, dcs, false);
    if (_S == null) {
      throw new ArgumentError("Exception occured in cs_sqr()");
    }
    _N = lu(dcs, _S, 1.0);
    if (_N == null) {
      throw new ArgumentError("Exception occured in cs_lu()");
    }
    if (checkIfSingular) {
      Decomposition D = dmperm(dcs, 1);
      /* check if matrix is singular */
      if (D != null && D.rr[3] < _n) {
        _isNonSingular = false;
      }
    }
  }

  double det() {
    if (!isNonsingular()) {
      return 0.0; // avoid rounding errors
    }
    int pivsign = 1;
    for (int i = 0; i < _n; i++) {
      if (_N.pinv[i] != i) {
        pivsign = -pivsign;
      }
    }
    if (_U == null) {
      _U = new SparseCCDoubleMatrix(_N.U);
      if (_rcMatrix) {
        _U = (_U as SparseCCDoubleMatrix).rowCompressed();
      }
    }
    double det = pivsign.toDouble();
    for (int j = 0; j < _n; j++) {
      det *= _U.get(j, j);
    }
    return det;
  }

  AbstractDoubleMatrix getL() {
    if (_L == null) {
      _L = new SparseCCDoubleMatrix(_N.L);
      if (_rcMatrix) {
        _L = (_L as SparseCCDoubleMatrix).rowCompressed();
      }
    }
    return _L.copy();
  }

  Int32List getPivot() {
    if (_N.pinv == null) {
      return null;
    }
    Int32List pinv = new Int32List(_N.pinv.length);
    //System.arraycopy(_N.pinv, 0, pinv, 0, pinv.length);
    pinv.setAll(0, _N.pinv);
    return pinv;
  }

  AbstractDoubleMatrix getU() {
    if (_U == null) {
      _U = new SparseCCDoubleMatrix(_N.U);
      if (_rcMatrix) {
        _U = (_U as SparseCCDoubleMatrix).rowCompressed();
      }
    }
    return _U.copy();
  }

  Symbolic getSymbolicAnalysis() {
    Symbolic S2 = new Symbolic();
    S2.cp = _S.cp != null ? _S.cp.clone() : null;
    S2.leftmost = _S.leftmost != null ? _S.leftmost.clone() : null;
    S2.lnz = _S.lnz;
    S2.m2 = _S.m2;
    S2.parent = _S.parent != null ? _S.parent.clone() : null;
    S2.pinv = _S.pinv != null ? _S.pinv.clone() : null;
    S2.q = _S.q != null ? _S.q.clone() : null;
    S2.unz = _S.unz;
    return S2;
  }

  bool isNonsingular() {
    return _isNonSingular;
  }

  void solve(AbstractDoubleVector b) {
    if (b.length != _n) {
      throw new ArgumentError("b.size() != A.rows()");
    }
    if (!isNonsingular()) {
      throw new ArgumentError("A is singular");
    }
    DoubleProperty.DEFAULT.checkDense(b);
    Float64List y = new Float64List(_n);
    Float64List x;
    if (b.isView) {
      x = b.copy().elements as Float64List;
    } else {
      x = b.elements as Float64List;
    }
    ipvec(_N.pinv, x, y, _n);
    /* y = b(p) */
    lsolve(_N.L, y);
    /* y = L\y */
    usolve(_N.U, y);
    /* y = U\y */
    ipvec(_S.q, y, x, _n);
    /* b(q) = x */

    if (b.isView) {
      b.setValues(x);
    }
  }
}




/**
 * LU decomposition implemented using JKLU.
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @author Richard Lincoln (r.w.lincoln@gmail.com)
 */
class SparseDoubleKLUDecomposition implements SparseDoubleLUDecomposition {
  klu.KLU_symbolic _S;
  klu.KLU_numeric _N;
  klu.KLU_common _Common;
  AbstractDoubleMatrix _L;
  AbstractDoubleMatrix _U;
  bool _rcMatrix = false;
  bool _isNonSingular = true;
  /**
     * Row and column dimension (square matrix).
     */
  int _n;

  /**
     * Constructs and returns a new LU Decomposition object; The decomposed
     * matrices can be retrieved via instance methods of the returned
     * decomposition object.
     *
     * @param A
     *            Square matrix
     * @param order
     *            ordering option (0 to 1); 0: AMD, 1: COLAMD
     * @param checkIfSingular
     *            if true, then the singularity test (based on
     *            BTFJ) is performed.
     * @param preOrder use BTF pre-ordering, or not
     * @throws ArgumentError
     *             if <tt>A</tt> is not square or is not sparse.
     * @throws ArgumentError
     *             if <tt>order</tt> is not in [0,1]
     */
  SparseDoubleKLUDecomposition(AbstractDoubleMatrix A, int order, bool checkIfSingular, [bool preOrder = true]) {
    DoubleProperty.DEFAULT.checkSquare(A);
    DoubleProperty.DEFAULT.checkSparse2D(A);

    if (order < 0 || order > 3) {
      throw new ArgumentError("order must be a number between 0 and 3");
    }

    _Common = new klu.KLU_common();
    klu.defaults(_Common);
    _Common.ordering = order;
    _Common.btf = preOrder ? 1 : 0;

    Matrix dcs;
    if (A is SparseRCDoubleMatrix) {
      _rcMatrix = true;
      dcs = A.getColumnCompressed().elements;
    } else {
      dcs = A.elements as Matrix;
    }
    _n = A.rows;
    Int32List Ap = dcs.p,
        Ai = dcs.i;
    Float64List Ax = dcs.x;

    _S = klu.analyze(_n, Ap, Ai, _Common);
    if (_S == null) {
      throw new ArgumentError("Exception occured in klu_analyze()");
    }
    _N = klu.factor(Ap, Ai, Ax, _S, _Common);
    if (_N == null) {
      throw new ArgumentError("Exception occured in klu_factor()");
    }
    if (checkIfSingular) {
      /* check if matrix is singular */
      int sprank = btf.maxtrans(_n, _n, Ap, Ai, _Common.maxwork, new Float64List(1), new Int32List(_n));
      if (sprank < _n) {
        _isNonSingular = false;
      }
    }
  }

  double det() {
    if (!isNonsingular()) {
      return 0.0; // avoid rounding errors
    }
    int pivsign = 1;
    for (int i = 0; i < _n; i++) {
      if (_N.Pinv[i] != i) {
        pivsign = -pivsign;
      }
    }
    if (_U == null) {
      _U = getU();
    }
    double det = pivsign.toDouble();
    for (int j = 0; j < _n; j++) {
      det *= _U.get(j, j);
    }
    return det;
  }

  AbstractDoubleMatrix getL() {
    if (_L == null) {
      Int32List Lp = new Int32List(_N.lnz + 1);
      Int32List Li = new Int32List(_N.lnz);
      Float64List Lx = new Float64List(_N.lnz);
      klu.extract(_N, _S, Lp, Li, Lx, null, null, null, null, null, null, null, null, null, null, _Common);
      _L = new SparseCCDoubleMatrix.SparseCCDoubleMatrix(_n, _n, Li, Lp, Lx);
      if (_rcMatrix) {
        _L = (_L as SparseCCDoubleMatrix).rowCompressed();
      }
    }
    return _L.copy();
  }

  Int32List getPivot() {
    if (_N.Pinv == null) return null;
    Int32List pinv = new Int32List(_N.Pinv.length);
    //System.arraycopy(_N.Pinv, 0, pinv, 0, pinv.length);
    pinv.setAll(0, _N.Pinv);
    return pinv;
  }

  AbstractDoubleMatrix getU() {
    if (_U == null) {
      Int32List Up = new Int32List(_N.unz + 1);
      Int32List Ui = new Int32List(_N.unz);
      Float64List Ux = new Float64List(_N.unz);
      klu.extract(_N, _S, null, null, null, Up, Ui, Ux, null, null, null, null, null, null, null, _Common);
      _U = new SparseCCDoubleMatrix.SparseCCDoubleMatrix(_n, _n, Ui, Up, Ux);
      if (_rcMatrix) {
        _U = (_U as SparseCCDoubleMatrix).rowCompressed();
      }
    }
    return _U.copy();
  }

  klu.KLU_symbolic getSymbolicAnalysis() {
    klu.KLU_symbolic S2 = new klu.KLU_symbolic();
    S2.symmetry = _S.symmetry;
    S2.est_flops = _S.est_flops;
    S2.lnz = _S.lnz;
    S2.unz = _S.unz;
    S2.Lnz = _S.Lnz.clone();
    S2.n = _S.n;
    S2.nz = _S.nz;
    S2.nzoff = _S.nzoff;
    S2.nblocks = _S.nblocks;
    S2.maxblock = _S.maxblock;
    S2.ordering = _S.ordering;
    S2.do_btf = _S.do_btf;
    S2.P = _S.P.clone();
    S2.Q = _S.Q.clone();
    S2.R = _S.R.clone();
    S2.structural_rank = _S.structural_rank;
    return S2;
  }

  bool isNonsingular() {
    return _isNonSingular;
  }

  void solve(AbstractDoubleVector b) {
    if (b.length != _n) {
      throw new ArgumentError("b.size() != A.rows()");
    }
    if (!isNonsingular()) {
      throw new ArgumentError("A is singular");
    }
    DoubleProperty.DEFAULT.checkDense(b);
    Float64List x;
    if (b.isView) {
      x = b.copy().elements as Float64List;
    } else {
      x = b.elements as Float64List;
    }
    klu.solve(_S, _N, _n, 1, x, 0, _Common);

    if (b.isView) {
      b.setValues(x);
    }
  }
}
