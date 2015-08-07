// Copyright (C) 1999 CERN - European Organization for Nuclear Research.
//
// Permission to use, copy, modify, distribute and sell this software and
// its documentation for any purpose is hereby granted without fee, provided
// that the above copyright notice appear in all copies and that both that
// copyright notice and this permission notice appear in supporting
// documentation.
//
// CERN makes no representations about the suitability of this software for
// any purpose. It is provided "as is" without expressed or implied warranty.
part of cern.colt.matrix.complex;

/// Dense 2-d matrix holding `complex` elements.
///
/// Internally holds one single contigous one-dimensional array, addressed in
/// row major order. Complex data is represented by 2 double values in sequence,
/// i.e. `elements[idx]` constitute the real part and `elements[idx+1]`
/// constitute the imaginary part, where
/// `idx = index(0,0) + row * rowStride + column * columnStride`.
class ComplexMatrix extends AbstractComplexMatrix {
  Float64List _elements;

  /// Constructs a complex matrix with the same size as [realPart] matrix and
  /// fills the real part of this matrix with elements of [realPart].
  factory ComplexMatrix.fromReal(AbstractDoubleMatrix realPart) {
    return new ComplexMatrix(realPart.rows, realPart.columns)
      ..setReal(realPart);
  }

  /// Constructs a matrix with a given number of rows and columns. All entries
  /// are initially `0`.
  factory ComplexMatrix(int rows, int columns) {
    final elements = new Float64List(rows * 2 * columns);
    return new ComplexMatrix._internal(
        rows, columns, elements, 0, 0, 2 * columns, 2, true);
  }

  ComplexMatrix._internal(int rows, int columns, Float64List elements,
      int rowZero, int columnZero, int rowStride, int columnStride,
      bool isNoView)
      : super(rows, columns, rowZero, columnZero, rowStride, columnStride,
          isNoView) {
    _elements = elements;
  }

  static ComplexMatrix create(int rows, int columns) {
    return new ComplexMatrix(rows, columns);
  }

  Float64List aggregate(final cfunc.ComplexComplexComplexFunction aggr,
      final cfunc.ComplexComplexFunction fn) {
    if (size == 0) {
      var b = new Float64List(2);
      b[0] = double.NAN;
      b[1] = double.NAN;
      return b;
    }
    int zero = index(0, 0);
    Float64List a =
        fn(new Float64List.fromList([_elements[zero], _elements[zero + 1]]));
    int d = 1; // first cell already done
    int idx;
    for (int r = 0; r < rows; r++) {
      for (int c = d; c < columns; c++) {
        idx = zero + r * rowStride + c * columnStride;
        a = aggr(a,
            fn(new Float64List.fromList([_elements[idx], _elements[idx + 1]])));
      }
      d = 0;
    }
    return a;
  }

  void apply(final cfunc.ComplexComplexFunction fn) {
    final int zero = index(0, 0);
    int idx = zero;
    var tmp = new Float64List(2);
    if (fn is cfunc.ComplexMult) {
      Float64List multiplicator = fn.multiplicator;
      // x[i] = mult*x[i]
      for (int r = 0; r < rows; r++) {
        for (int i = idx, c = 0; c < columns; c++) {
          tmp[0] = _elements[i];
          tmp[1] = _elements[i + 1];
          _elements[i] = tmp[0] * multiplicator[0] - tmp[1] * multiplicator[1];
          _elements[i + 1] =
              tmp[1] * multiplicator[0] + tmp[0] * multiplicator[1];
          i += columnStride;
        }
        idx += rowStride;
      }
    } else {
      for (int r = 0; r < rows; r++) {
        for (int i = idx, c = 0; c < columns; c++) {
          tmp = fn(new Float64List.fromList([_elements[i], _elements[i + 1]]));
          _elements[i] = tmp[0];
          _elements[i + 1] = tmp[1];
          i += columnStride;
        }
        idx += rowStride;
      }
    }
  }

  void applyReal(final cfunc.ComplexRealFunction fn) {
    final int zero = index(0, 0);
    int idx = zero;
    var tmp = new Float64List(2);
    if (fn == cfunc.abs) {
      for (int r = 0; r < rows; r++) {
        for (int i = idx, c = 0; c < columns; c++) {
          tmp[0] = _elements[i];
          tmp[1] = _elements[i + 1];
          double absX = tmp[0].abs();
          double absY = tmp[1].abs();
          if (absX == 0 && absY == 0) {
            _elements[i] = 0.0;
          } else if (absX >= absY) {
            double d = tmp[1] / tmp[0];
            _elements[i] = absX * Math.sqrt(1 + d * d);
          } else {
            double d = tmp[0] / tmp[1];
            _elements[i] = absY * Math.sqrt(1 + d * d);
          }
          _elements[i + 1] = 0.0;
          i += columnStride;
        }
        idx += rowStride;
      }
    } else {
      for (int r = 0; r < rows; r++) {
        for (int i = idx, c = 0; c < columns; c++) {
          tmp[0] = _elements[i];
          tmp[1] = _elements[i + 1];
          tmp[0] = fn(tmp);
          _elements[i] = tmp[0];
          _elements[i + 1] = 0.0;
          i += columnStride;
        }
        idx += rowStride;
      }
    }
  }

  void copyFrom(final AbstractComplexMatrix source) {
    // overriden for performance only
    if (source is! ComplexMatrix) {
      super.copyFrom(source);
      return;
    }
    ComplexMatrix other = source as ComplexMatrix;
    if (other == this) {
      return; // nothing to do
    }
    checkShape(this, other);
    if (!isView && !other.isView) {
      // quickest
      _elements.setAll(0, other._elements);
      return;
    }
    if (_haveSharedCells(other)) {
      AbstractComplexMatrix c = other.copy();
      if (c is! ComplexMatrix) {
        // should not happen
        super.copyFrom(other);
        return;
      }
      other = c as ComplexMatrix;
    }

    Float64List elemsOther = other._elements;
    if (_elements == null || elemsOther == null) {
      throw new Error();
    }
    int columnStrideOther = other.columnStride;
    int rowStrideOther = other.rowStride;
    int zeroOther = other.index(0, 0);
    int zero = index(0, 0);
    int idx = zero;
    int idxOther = zeroOther;
    for (int r = 0; r < rows; r++) {
      for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
        _elements[i] = elemsOther[j];
        _elements[i + 1] = elemsOther[j + 1];
        i += columnStride;
        j += columnStrideOther;
      }
      idx += rowStride;
      idxOther += rowStrideOther;
    }
  }

  void assign(final AbstractComplexMatrix y,
      final cfunc.ComplexComplexComplexFunction fn) {
    // overriden for performance only
    if (y is! ComplexMatrix) {
      super.assign(y, fn);
      return;
    }
    checkShape(this, y);
    Float64List elemsOther = y.elements;
    if (_elements == null || elemsOther == null) {
      throw new Error();
    }
    int columnStrideOther = y.columnStride;
    int rowStrideOther = y.rowStride;
    int zeroOther = y.index(0, 0);
    int zero = index(0, 0);
    Float64List tmp1 = new Float64List(2);
    Float64List tmp2 = new Float64List(2);
    int idx = zero;
    int idxOther = zeroOther;
    if (fn == cfunc.mult) {
      for (int r = 0; r < rows; r++) {
        for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
          tmp1[0] = _elements[i];
          tmp1[1] = _elements[i + 1];
          tmp2[0] = elemsOther[j];
          tmp2[1] = elemsOther[j + 1];
          _elements[i] = tmp1[0] * tmp2[0] - tmp1[1] * tmp2[1];
          _elements[i + 1] = tmp1[1] * tmp2[0] + tmp1[0] * tmp2[1];
          i += columnStride;
          j += columnStrideOther;
        }
        idx += rowStride;
        idxOther += rowStrideOther;
      }
    } else if (fn == cfunc.multConjFirst) {
      for (int r = 0; r < rows; r++) {
        for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
          tmp1[0] = _elements[i];
          tmp1[1] = _elements[i + 1];
          tmp2[0] = elemsOther[j];
          tmp2[1] = elemsOther[j + 1];
          _elements[i] = tmp1[0] * tmp2[0] + tmp1[1] * tmp2[1];
          _elements[i + 1] = -tmp1[1] * tmp2[0] + tmp1[0] * tmp2[1];
          i += columnStride;
          j += columnStrideOther;
        }
        idx += rowStride;
        idxOther += rowStrideOther;
      }
    } else if (fn == cfunc.multConjSecond) {
      for (int r = 0; r < rows; r++) {
        for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
          tmp1[0] = _elements[i];
          tmp1[1] = _elements[i + 1];
          tmp2[0] = elemsOther[j];
          tmp2[1] = elemsOther[j + 1];
          _elements[i] = tmp1[0] * tmp2[0] + tmp1[1] * tmp2[1];
          _elements[i + 1] = tmp1[1] * tmp2[0] - tmp1[0] * tmp2[1];
          i += columnStride;
          j += columnStrideOther;
        }
        idx += rowStride;
        idxOther += rowStrideOther;
      }
    } else {
      for (int r = 0; r < rows; r++) {
        for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
          tmp1[0] = _elements[i];
          tmp1[1] = _elements[i + 1];
          tmp2[0] = elemsOther[j];
          tmp2[1] = elemsOther[j + 1];
          tmp1 = fn(tmp1, tmp2);
          _elements[i] = tmp1[0];
          _elements[i + 1] = tmp1[1];
          i += columnStride;
          j += columnStrideOther;
        }
        idx += rowStride;
        idxOther += rowStrideOther;
      }
    }
  }

  void fill(final double re, final double im) {
    final int zero = index(0, 0);
    int idx = zero;
    for (int r = 0; r < rows; r++) {
      for (int i = idx, c = 0; c < columns; c++) {
        _elements[i] = re;
        _elements[i + 1] = im;
        i += columnStride;
      }
      idx += rowStride;
    }
  }

  void setAll(final Float64List values) {
    if (values.length != rows * 2 * columns) {
      throw new ArgumentError(
          "Must have same length: length=${values.length} rows()*2*columns()=${rows * 2 * columns}");
    }
    if (!isView) {
      _elements.setAll(0, values);
    } else {
      final int zero = index(0, 0);
      int idxOther = 0;
      int idx = zero;
      for (int r = 0; r < rows; r++) {
        for (int i = idx, c = 0; c < columns; c++) {
          _elements[i] = values[idxOther++];
          _elements[i + 1] = values[idxOther++];
          i += columnStride;
        }
        idx += rowStride;
      }
    }
  }

  void setImaginary(final AbstractDoubleMatrix other) {
    checkShape(this, other);
    int columnStrideOther = other.columnStride;
    int rowStrideOther = other.rowStride;
    int zeroOther = other.index(0, 0);
    int zero = index(0, 0);
    Float64List elemsOther = (other as DoubleMatrix).elements;
    int idx = zero;
    int idxOther = zeroOther;
    for (int r = 0; r < rows; r++) {
      for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
        _elements[i + 1] = elemsOther[j];
        i += columnStride;
        j += columnStrideOther;
      }
      idx += rowStride;
      idxOther += rowStrideOther;
    }
  }

  void setReal(final AbstractDoubleMatrix other) {
    checkShape(this, other);
    int columnStrideOther = other.columnStride;
    int rowStrideOther = other.rowStride;
    int zeroOther = other.index(0, 0);
    int zero = index(0, 0);
    Float64List elemsOther = (other as DoubleMatrix).elements;
    int idx = zero;
    int idxOther = zeroOther;
    for (int r = 0; r < rows; r++) {
      for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
        _elements[i] = elemsOther[j];
        i += columnStride;
        j += columnStrideOther;
      }
      idx += rowStride;
      idxOther += rowStrideOther;
    }
  }

  int get cardinality {
    int cardinality = 0;
    final int zero = index(0, 0);
    int idx = zero;
    for (int r = 0; r < rows; r++) {
      for (int i = idx, c = 0; c < columns; c++) {
        if ((_elements[i] != 0.0) || (_elements[i + 1] != 0.0)) {
          cardinality++;
        }
        i += columnStride;
      }
      idx += rowStride;
    }
    return cardinality;
  }

  void forEachNonZero(final cfunc.IntIntComplexFunction function) {
    final int zero = index(0, 0);
    int idx = zero;
    var value = new Float64List(2);
    for (int r = 0; r < rows; r++) {
      for (int i = idx, c = 0; c < columns; c++) {
        value[0] = _elements[i];
        value[1] = _elements[i + 1];
        if (value[0] != 0 || value[1] != 0) {
          Float64List v = function(r, c, value);
          _elements[i] = v[0];
          _elements[i + 1] = v[1];
        }
        i += columnStride;
      }
      idx += rowStride;
    }
  }

  AbstractComplexMatrix conjugateTranspose() {
    AbstractComplexMatrix transpose = dice().copy();
    Float64List elemsOther = (transpose as ComplexMatrix)._elements;
    int zeroOther = transpose.index(0, 0);
    int columnStrideOther = transpose.columnStride;
    //int rowStrideOther = transpose.rowStride;
    int columnsOther = transpose.columns;
    int rowsOther = transpose.rows;
    int idxOther = zeroOther;
    for (int r = 0; r < rowsOther; r++) {
      for (int c = 0; c < columnsOther; c++) {
        elemsOther[idxOther + 1] = -elemsOther[idxOther + 1];
        idxOther += columnStrideOther;
      }
    }
    return transpose;
  }

  Float64List get elements => _elements;

  AbstractDoubleMatrix imaginary() {
    var Im = new DoubleMatrix(rows, columns);
    Float64List elemsOther = Im.elements;
    int columnStrideOther = Im.columnStride;
    int rowStrideOther = Im.rowStride;
    int zeroOther = Im.index(0, 0);
    int zero = index(0, 0);
    int idx = zero;
    int idxOther = zeroOther;
    for (int r = 0; r < rows; r++) {
      for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
        elemsOther[j] = _elements[i + 1];
        i += columnStride;
        j += columnStrideOther;
      }
      idx += rowStride;
      idxOther += rowStrideOther;
    }
    return Im;
  }

  void nonzero(final List<int> rowList, final List<int> columnList,
      final List<List<double>> valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    int idx = index(0, 0);
    var value = new Float64List(2);
    for (int r = 0; r < rows; r++) {
      for (int i = idx, c = 0; c < columns; c++) {
        value[0] = _elements[i];
        value[1] = _elements[i + 1];
        if (value[0] != 0 || value[1] != 0) {
          rowList.add(r);
          columnList.add(c);
          valueList.add(value);
        }
        i += columnStride;
      }
      idx += rowStride;
    }
  }

  Float64List get(int row, int column) {
    int idx = rowZero + row * rowStride + columnZero + column * columnStride;
    return new Float64List.fromList([_elements[idx], _elements[idx + 1]]);
  }

  AbstractDoubleMatrix real() {
    var R = new DoubleMatrix(rows, columns);
    Float64List elemsOther = R.elements;
    int columnStrideOther = R.columnStride;
    int rowStrideOther = R.rowStride;
    int zeroOther = R.index(0, 0);
    int zero = index(0, 0);
    int idx = zero;
    int idxOther = zeroOther;
    for (int r = 0; r < rows; r++) {
      for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
        elemsOther[j] = _elements[i];
        i += columnStride;
        j += columnStrideOther;
      }
      idx += rowStride;
      idxOther += rowStrideOther;
    }
    return R;
  }

  AbstractComplexMatrix like2D(int rows, int columns) {
    return new ComplexMatrix(rows, columns);
  }

  AbstractComplexVector like1D(int size) => new ComplexVector(size);

  void setParts(int row, int column, double re, double im) {
    int idx = rowZero + row * rowStride + columnZero + column * columnStride;
    _elements[idx] = re;
    _elements[idx + 1] = im;
  }

  void set(int row, int column, Float64List value) {
    int idx = rowZero + row * rowStride + columnZero + column * columnStride;
    _elements[idx] = value[0];
    _elements[idx + 1] = value[1];
  }

  AbstractComplexVector mult(final AbstractComplexVector y,
      [AbstractComplexVector z = null, Float64List alpha = null,
      Float64List beta = null, bool transposeA = false]) {
    if (alpha == null) {
      alpha = new Float64List.fromList([1.0, 0.0]);
    }
    if (beta == null) {
      beta = (z == null
          ? new Float64List.fromList([1.0, 0.0])
          : new Float64List.fromList([0.0, 0.0]));
    }
    if (transposeA) {
      return conjugateTranspose().mult(y, z, alpha, beta, false);
    }
    AbstractComplexVector zz;
    if (z == null) {
      zz = new ComplexVector(rows);
    } else {
      zz = z;
    }
    if (columns != y.size || rows > zz.size) {
      throw new ArgumentError("Incompatible args: " +
          toStringShort() +
          ", " +
          y.toStringShort() +
          ", " +
          zz.toStringShort());
    }
    Float64List elemsY = y.elements as Float64List;
    Float64List elemsZ = zz.elements as Float64List;
    if (_elements == null || elemsY == null || elemsZ == null) {
      throw new Error();
    }
    int strideY = y.stride;
    int strideZ = zz.stride;
    int zero = index(0, 0);
    int zeroY = y.index(0);
    int zeroZ = zz.index(0);
    int idxZero = zero;
    int idxZeroZ = zeroZ;

    for (int r = 0; r < rows; r++) {
      var reS = 0.0;
      var imS = 0.0;
      int idx = idxZero;
      int idxY = zeroY;
      for (int c = 0; c < columns; c++) {
        var reA = _elements[idx];
        var imA = _elements[idx + 1];
        var reY = elemsY[idxY];
        var imY = elemsY[idxY + 1];
        reS += reA * reY - imA * imY;
        imS += imA * reY + reA * imY;
        idx += columnStride;
        idxY += strideY;
      }
      var reZ = elemsZ[idxZeroZ];
      var imZ = elemsZ[idxZeroZ + 1];
      elemsZ[idxZeroZ] =
          reS * alpha[0] - imS * alpha[1] + reZ * beta[0] - imZ * beta[1];
      elemsZ[idxZeroZ + 1] =
          imS * alpha[0] + reS * alpha[1] + imZ * beta[0] + reZ * beta[1];
      idxZero += rowStride;
      idxZeroZ += strideZ;
    }
    return zz;
  }

  AbstractComplexMatrix multiply(AbstractComplexMatrix B,
      [AbstractComplexMatrix C = null, Float64List alpha = null,
      Float64List beta = null, final bool transposeA = false,
      final bool transposeB = false]) {
    if (alpha == null) {
      alpha = new Float64List.fromList([1.0, 0.0]);
    }
    if (beta == null) {
      beta = (C == null
          ? new Float64List.fromList([1.0, 0.0])
          : new Float64List.fromList([0.0, 0.0]));
    }
    if (transposeA) {
      return conjugateTranspose().multiply(
          B, C, alpha, beta, false, transposeB);
    }
    if (transposeB) {
      return this.multiply(
          B.conjugateTranspose(), C, alpha, beta, transposeA, false);
    }
    int m = rows;
    int n = columns;
    int p = B.columns;
    if (C == null) {
      C = new ComplexMatrix(m, p);
    }
    if (C is! ComplexMatrix) {
      return super.multiply(B, C, alpha, beta, transposeA, transposeB);
    }
    if (B.rows != n) {
      throw new ArgumentError("Matrix inner dimensions must agree:" +
          toStringShort() +
          ", " +
          B.toStringShort());
    }
    if (C.rows != m || C.columns != p) {
      throw new ArgumentError("Incompatibel result matrix: " +
          toStringShort() +
          ", " +
          B.toStringShort() +
          ", " +
          C.toStringShort());
    }
    if (this == C || B == C) {
      throw new ArgumentError("Matrices must not be identical");
    }

    ComplexMatrix BB = B as ComplexMatrix;
    ComplexMatrix CC = C as ComplexMatrix;
    Float64List AElems = _elements;
    Float64List BElems = BB._elements;
    Float64List CElems = CC._elements;
    if (AElems == null || BElems == null || CElems == null) {
      throw new Error();
    }

    int cA = columnStride;
    int cB = BB.columnStride;
    int cC = CC.columnStride;

    int rA = rowStride;
    int rB = BB.rowStride;
    int rC = CC.rowStride;

    // A is blocked to hide memory latency xxxxxxx B xxxxxxx xxxxxxx A xxx
    // xxxxxxx C xxx xxxxxxx --- ------- xxx xxxxxxx xxx xxxxxxx --- -------
    // xxx xxxxxxx
    const BLOCK_SIZE = 30000; // * 8 == Level 2 cache in bytes
    int m_optimal = (BLOCK_SIZE - n) ~/ (n + 1);
    if (m_optimal <= 0) {
      m_optimal = 1;
    }
    int blocks = m ~/ m_optimal;
    int rr = 0;
    if (m % m_optimal != 0) {
      blocks++;
    }
    for (; --blocks >= 0;) {
      int jB = BB.index(0, 0);
      int indexA = index(rr, 0);
      int jC = CC.index(rr, 0);
      rr += m_optimal;
      if (blocks == 0) {
        m_optimal += m - rr;
      }

      for (int j = p; --j >= 0;) {
        int iA = indexA;
        int iC = jC;
        for (int i = m_optimal; --i >= 0;) {
          int kA = iA;
          int kB = jB;
          var reS = 0.0;
          var imS = 0.0;
          // loop unrolled
          kA -= cA;
          kB -= rB;
          for (int k = n % 4; --k >= 0;) {
            kA += cA;
            kB += rB;
            var reA = AElems[kA];
            var imA = AElems[kA + 1];
            var reB = BElems[kB];
            var imB = BElems[kB + 1];
            reS += reA * reB - imA * imB;
            imS += imA * reB + reA * imB;
          }
          for (int k = n ~/ 4; --k >= 0;) {
            kA += cA;
            kB += rB;
            var reA = AElems[kA];
            var imA = AElems[kA + 1];
            var reB = BElems[kB];
            var imB = BElems[kB + 1];
            reS += reA * reB - imA * imB;
            imS += imA * reB + reA * imB;
            kA += cA;
            kB += rB;
            reA = AElems[kA];
            imA = AElems[kA + 1];
            reB = BElems[kB];
            imB = BElems[kB + 1];
            reS += reA * reB - imA * imB;
            imS += imA * reB + reA * imB;
            kA += cA;
            kB += rB;
            reA = AElems[kA];
            imA = AElems[kA + 1];
            reB = BElems[kB];
            imB = BElems[kB + 1];
            reS += reA * reB - imA * imB;
            imS += imA * reB + reA * imB;
            kA += cA;
            kB += rB;
            reA = AElems[kA];
            imA = AElems[kA + 1];
            reB = BElems[kB];
            imB = BElems[kB + 1];
            reS += reA * reB - imA * imB;
            imS += imA * reB + reA * imB;
          }
          var reC = CElems[iC];
          var imC = CElems[iC + 1];
          CElems[iC] =
              alpha[0] * reS - alpha[1] * imS + beta[0] * reC - beta[1] * imC;
          CElems[iC + 1] =
              alpha[1] * reS + alpha[0] * imS + beta[1] * reC + beta[0] * imC;
          iA += rA;
          iC += rC;
        }
        jB += cB;
        jC += cC;
      }
    }
    return C;
  }

  Float64List sum() {
    var sum = new Float64List(2);
    int zero = index(0, 0);
    int idx = zero;
    for (int r = 0; r < rows; r++) {
      for (int i = idx, c = 0; c < columns; c++) {
        sum[0] += _elements[i];
        sum[1] += _elements[i + 1];
        i += columnStride;
      }
      idx += rowStride;
    }
    return sum;
  }

  bool _haveSharedCellsRaw(AbstractComplexMatrix other) {
    if (other is SelectedDenseComplexMatrix) {
      return _elements == other._elements;
    } else if (other is ComplexMatrix) {
      return _elements == other._elements;
    }
    return false;
  }

  int index(int row, int column) {
    return rowZero + row * rowStride + columnZero + column * columnStride;
  }

  AbstractComplexVector _like1D(int size, int zero, int stride) {
    return new ComplexVector._internal(
        size, this._elements, zero, stride, false);
  }

  AbstractComplexMatrix _viewSelectionLike(
      Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedDenseComplexMatrix(
        this._elements, rowOffsets, columnOffsets, 0);
  }

  Object clone() {
    return new ComplexMatrix._internal(rows, columns, _elements, rowZero,
        columnZero, rowStride, columnStride, !isView);
  }
}

/// Selection view on dense 2-d matrices holding `complex` elements.
///
/// Objects of this class are typically constructed via `select` methods on
/// some source matrix. The interface introduced in abstract super classes
/// defines everything a user can do. From a user point of view there is
/// nothing special about this class; it presents the same functionality with
/// the same signatures and semantics as its abstract superclass(es) while
/// introducing no additional functionality. Thus, this class need not be
/// visible to users.
class SelectedDenseComplexMatrix extends AbstractComplexMatrix {
  Float64List _elements;

  /// The offsets of the visible cells of this matrix.
  Int32List _rowOffsets, _columnOffsets;

  int _offset;

  factory SelectedDenseComplexMatrix(Float64List elements, Int32List rowOffsets,
      Int32List columnOffsets, int offset) {
    return new SelectedDenseComplexMatrix._internal(rowOffsets.length,
        columnOffsets.length, elements, 0, 0, 1, 1, rowOffsets, columnOffsets,
        offset);
  }

  SelectedDenseComplexMatrix._internal(int rows, int columns,
      Float64List elements, int rowZero, int columnZero, int rowStride,
      int columnStride, Int32List rowOffsets, Int32List columnOffsets,
      int offset)
      : super(
          rows, columns, rowZero, columnZero, rowStride, columnStride, false) {
    _elements = elements;
    _rowOffsets = rowOffsets;
    _columnOffsets = columnOffsets;
    _offset = offset;
  }

  int columnOffset(int absRank) => _columnOffsets[absRank];

  int rowOffset(int absRank) => _rowOffsets[absRank];

  Float64List get(int row, int column) {
    int idxr = rowZero + row * rowStride;
    int idxc = columnZero + column * columnStride;
    return new Float64List.fromList([
      _elements[_offset + _rowOffsets[idxr] + _columnOffsets[idxc]],
      _elements[_offset + _rowOffsets[idxr] + _columnOffsets[idxc] + 1]
    ]);
  }

  // This method is not supported.
  Object get elements {
    throw new UnsupportedError("This method is not supported.");
  }

  bool _haveSharedCellsRaw(AbstractComplexMatrix other) {
    if (other is SelectedDenseComplexMatrix) {
      return _elements == other._elements;
    } else if (other is ComplexMatrix) {
      return _elements == other._elements;
    }
    return false;
  }

  int index(int row, int column) {
    return _offset +
        _rowOffsets[rowZero + row * rowStride] +
        _columnOffsets[columnZero + column * columnStride];
  }

  AbstractComplexMatrix like2D(int rows, int columns) {
    return new ComplexMatrix(rows, columns);
  }

  AbstractComplexVector like1D(int size) => new ComplexVector(size);

  AbstractComplexVector _like1D(int size, int zero, int stride) {
    // never called since row() and column() are overridden
    throw new Error();
  }

  void set(int row, int column, Float64List value) {
    int idxr = rowZero + row * rowStride;
    int idxc = columnZero + column * columnStride;
    _elements[_offset + _rowOffsets[idxr] + _columnOffsets[idxc]] = value[0];
    _elements[_offset + _rowOffsets[idxr] + _columnOffsets[idxc] + 1] =
        value[1];
  }

  // This method is not supported.
  AbstractComplexVector vectorize() {
    throw new UnsupportedError("This method is not supported.");
  }

  void setParts(int row, int column, double re, double im) {
    int idxr = rowZero + row * rowStride;
    int idxc = columnZero + column * columnStride;
    _elements[_offset + _rowOffsets[idxr] + _columnOffsets[idxc]] = re;
    _elements[_offset + _rowOffsets[idxr] + _columnOffsets[idxc] + 1] = im;
  }

  AbstractComplexMatrix dice() {
    var v = _view();
    vDice(v);
    Int32List tmp = _rowOffsets;
    _rowOffsets = _columnOffsets;
    _columnOffsets = tmp;
    setIsNoView(v, false);
    return v;
  }

  AbstractComplexVector column(int column) {
    checkColumn(this, column);
    int viewSize = rows;
    int viewZero = rowZero;
    int viewStride = rowStride;
    Int32List viewOffsets = _rowOffsets;
    int viewOffset = _offset + columnOffset(columnRank(column));
    return new SelectedDenseComplexVector._internal(viewSize, this._elements,
        viewZero, viewStride, viewOffsets, viewOffset);
  }

  AbstractComplexVector row(int row) {
    checkRow(this, row);
    int viewSize = columns;
    int viewZero = columnZero;
    int viewStride = columnStride;
    Int32List viewOffsets = _columnOffsets;
    int viewOffset = _offset + rowOffset(rowRank(row));
    return new SelectedDenseComplexVector._internal(
        viewSize, _elements, viewZero, viewStride, viewOffsets, viewOffset);
  }

  AbstractComplexMatrix _viewSelectionLike(
      Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedDenseComplexMatrix(
        this._elements, rowOffsets, columnOffsets, _offset);
  }

  AbstractDoubleMatrix real() {
    var R = new DoubleMatrix(rows, columns);
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        final tmp = get(r, c);
        R.set(r, c, tmp[0]);
      }
    }
    return R;
  }

  AbstractDoubleMatrix imaginary() {
    var Im = new DoubleMatrix(rows, columns);
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        final tmp = get(r, c);
        Im.set(r, c, tmp[1]);
      }
    }
    return Im;
  }

  Object clone() {
    return new SelectedDenseComplexMatrix._internal(rows, columns, _elements,
        rowZero, columnZero, rowStride, columnStride, _rowOffsets,
        _columnOffsets, _offset);
  }
}
