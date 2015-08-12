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

/// Dense 2-d matrix holding [Complex] elements.
///
/// Internally holds one single contigous one-dimensional array, addressed in
/// row major order. Complex data is represented by 2 double values in sequence,
/// i.e. `elements[idx]` constitute the real part and `elements[idx+1]`
/// constitute the imaginary part, where
/// `idx = index(0,0) + row * rowStride + column * columnStride`.
class DenseComplexMatrix extends ComplexMatrix {
  Float64List _elements;

  /// Constructs a complex matrix with the same size as [realPart] matrix and
  /// fills the real part of this matrix with elements of [realPart].
  factory DenseComplexMatrix.fromReal(DoubleMatrix realPart) {
    return new DenseComplexMatrix(realPart.rows, realPart.columns)
      ..setReal(realPart);
  }

  /// Constructs a matrix with a given number of rows and columns. All entries
  /// are initially `0`.
  factory DenseComplexMatrix(int rows, int columns) {
    var elements = new Float64List(rows * 2 * columns);
    return new DenseComplexMatrix._internal(
        rows, columns, elements, 0, 0, 2 * columns, 2, true);
  }

  DenseComplexMatrix._internal(int rows, int columns, Float64List elements,
      int rowZero, int columnZero, int rowStride, int columnStride,
      bool isNoView)
      : super._(rows, columns, rowZero, columnZero, rowStride, columnStride,
          isNoView) {
    _elements = elements;
  }

  static DenseComplexMatrix create(int rows, int columns) {
    return new DenseComplexMatrix(rows, columns);
  }

  Complex aggregate(cfunc.ComplexComplexComplexFunction aggr,
      cfunc.ComplexComplexFunction fn) {
    if (size == 0) {
      return Complex.NAN;
    }
    int zero = index(0, 0);
    var a = fn(new Complex(_elements[zero], _elements[zero + 1]));
    int d = 1; // first cell already done
    int idx;
    for (int r = 0; r < rows; r++) {
      for (int c = d; c < columns; c++) {
        idx = zero + r * rowStride + c * columnStride;
        a = aggr(a, fn(new Complex(_elements[idx], _elements[idx + 1])));
      }
      d = 0;
    }
    return a;
  }

  void apply(cfunc.ComplexComplexFunction fn) {
    var zero = index(0, 0);
    int idx = zero;
    if (fn is cfunc.ComplexMult) {
      var multiplicator = fn.multiplicator;
      // x[i] = mult*x[i]
      for (int r = 0; r < rows; r++) {
        for (int i = idx, c = 0; c < columns; c++) {
          var tmp = new Complex(_elements[i], _elements[i + 1]);
          _elements[i] = tmp.real * multiplicator.real -
              tmp.imaginary * multiplicator.imaginary;
          _elements[i + 1] = tmp.imaginary * multiplicator.real +
              tmp.real * multiplicator.imaginary;
          i += columnStride;
        }
        idx += rowStride;
      }
    } else {
      for (int r = 0; r < rows; r++) {
        for (int i = idx, c = 0; c < columns; c++) {
          var tmp = fn(new Complex(_elements[i], _elements[i + 1]));
          _elements[i] = tmp.real;
          _elements[i + 1] = tmp.imaginary;
          i += columnStride;
        }
        idx += rowStride;
      }
    }
  }

  void applyReal(cfunc.ComplexRealFunction fn) {
    var zero = index(0, 0);
    int idx = zero;
    if (fn == cfunc.abs) {
      for (int r = 0; r < rows; r++) {
        for (int i = idx, c = 0; c < columns; c++) {
          var tmp = new Complex(_elements[i], _elements[i + 1]);
          double absX = tmp.real.abs();
          double absY = tmp.imaginary.abs();
          if (absX == 0 && absY == 0) {
            _elements[i] = 0.0;
          } else if (absX >= absY) {
            double d = tmp.imaginary / tmp.real;
            _elements[i] = absX * Math.sqrt(1 + d * d);
          } else {
            double d = tmp.real / tmp.imaginary;
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
          var tmp = new Complex(_elements[i], _elements[i + 1]);
          var re = fn(tmp);
          _elements[i] = re;
          _elements[i + 1] = 0.0;
          i += columnStride;
        }
        idx += rowStride;
      }
    }
  }

  void copyFrom(ComplexMatrix source) {
    // overriden for performance only
    if (source is! DenseComplexMatrix) {
      super.copyFrom(source);
      return;
    }
    DenseComplexMatrix other = source as DenseComplexMatrix;
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
      ComplexMatrix c = other.copy();
      if (c is! DenseComplexMatrix) {
        // should not happen
        super.copyFrom(other);
        return;
      }
      other = c as DenseComplexMatrix;
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

  void assign(ComplexMatrix y, cfunc.ComplexComplexComplexFunction fn) {
    // overriden for performance only
    if (y is! DenseComplexMatrix) {
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
    int idx = zero;
    int idxOther = zeroOther;
    if (fn == cfunc.mult) {
      for (int r = 0; r < rows; r++) {
        for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
          var tmp1 = new Complex(_elements[i], _elements[i + 1]);
          var tmp2 = new Complex(elemsOther[j], elemsOther[j + 1]);
          _elements[i] =
              tmp1.real * tmp2.real - tmp1.imaginary * tmp2.imaginary;
          _elements[i + 1] =
              tmp1.imaginary * tmp2.real + tmp1.real * tmp2.imaginary;
          i += columnStride;
          j += columnStrideOther;
        }
        idx += rowStride;
        idxOther += rowStrideOther;
      }
    } else if (fn == cfunc.multConjFirst) {
      for (int r = 0; r < rows; r++) {
        for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
          var tmp1 = new Complex(_elements[i], _elements[i + 1]);
          var tmp2 = new Complex(elemsOther[j], elemsOther[j + 1]);
          _elements[i] =
              tmp1.real * tmp2.real + tmp1.imaginary * tmp2.imaginary;
          _elements[i + 1] =
              -tmp1.imaginary * tmp2.real + tmp1.real * tmp2.imaginary;
          i += columnStride;
          j += columnStrideOther;
        }
        idx += rowStride;
        idxOther += rowStrideOther;
      }
    } else if (fn == cfunc.multConjSecond) {
      for (int r = 0; r < rows; r++) {
        for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
          var tmp1 = new Complex(_elements[i], _elements[i + 1]);
          var tmp2 = new Complex(elemsOther[j], elemsOther[j + 1]);
          _elements[i] =
              tmp1.real * tmp2.real + tmp1.imaginary * tmp2.imaginary;
          _elements[i + 1] =
              tmp1.imaginary * tmp2.real - tmp1.real * tmp2.imaginary;
          i += columnStride;
          j += columnStrideOther;
        }
        idx += rowStride;
        idxOther += rowStrideOther;
      }
    } else {
      for (int r = 0; r < rows; r++) {
        for (int i = idx, j = idxOther, c = 0; c < columns; c++) {
          var tmp1 = new Complex(_elements[i], _elements[i + 1]);
          var tmp2 = new Complex(elemsOther[j], elemsOther[j + 1]);
          tmp1 = fn(tmp1, tmp2);
          _elements[i] = tmp1.real;
          _elements[i + 1] = tmp1.imaginary;
          i += columnStride;
          j += columnStrideOther;
        }
        idx += rowStride;
        idxOther += rowStrideOther;
      }
    }
  }

  void fill(double re, double im) {
    var zero = index(0, 0);
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

  void setValues(Float64List values) {
    if (values.length != rows * 2 * columns) {
      throw new ArgumentError(
          "Must have same length: length=${values.length} rows()*2*columns()=${rows * 2 * columns}");
    }
    if (!isView) {
      _elements.setAll(0, values);
    } else {
      var zero = index(0, 0);
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

  void setImaginary(DoubleMatrix other) {
    checkShape(this, other);
    int columnStrideOther = other.columnStride;
    int rowStrideOther = other.rowStride;
    int zeroOther = other.index(0, 0);
    int zero = index(0, 0);
    Float64List elemsOther = (other as DenseDoubleMatrix).elements;
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

  void setReal(DoubleMatrix other) {
    checkShape(this, other);
    int columnStrideOther = other.columnStride;
    int rowStrideOther = other.rowStride;
    int zeroOther = other.index(0, 0);
    int zero = index(0, 0);
    Float64List elemsOther = (other as DenseDoubleMatrix).elements;
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
    var zero = index(0, 0);
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

  void forEachNonZero(cfunc.IntIntComplexFunction fn) {
    var zero = index(0, 0);
    int idx = zero;
    for (int r = 0; r < rows; r++) {
      for (int i = idx, c = 0; c < columns; c++) {
        var value = new Complex(_elements[i], _elements[i + 1]);
        if (value.real != 0 || value.imaginary != 0) {
          var v = fn(r, c, value);
          _elements[i] = v.real;
          _elements[i + 1] = v.imaginary;
        }
        i += columnStride;
      }
      idx += rowStride;
    }
  }

  ComplexMatrix conjugateTranspose() {
    ComplexMatrix transpose = dice().copy();
    Float64List elemsOther = (transpose as DenseComplexMatrix)._elements;
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

  dynamic get elements => _elements;

  DoubleMatrix imaginary() {
    var Im = new DenseDoubleMatrix(rows, columns);
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

  void nonzero(
      List<int> rowList, List<int> columnList, List<Complex> valueList) {
    rowList.clear();
    columnList.clear();
    valueList.clear();
    int idx = index(0, 0);
    for (int r = 0; r < rows; r++) {
      for (int i = idx, c = 0; c < columns; c++) {
        var value = new Complex(_elements[i], _elements[i + 1]);
        if (value.real != 0 || value.imaginary != 0) {
          rowList.add(r);
          columnList.add(c);
          valueList.add(value);
        }
        i += columnStride;
      }
      idx += rowStride;
    }
  }

  Complex get(int row, int column) {
    int idx = rowZero + row * rowStride + columnZero + column * columnStride;
    return new Complex(_elements[idx], _elements[idx + 1]);
  }

  DoubleMatrix real() {
    var R = new DenseDoubleMatrix(rows, columns);
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

  ComplexMatrix like2D(int rows, int columns) {
    return new DenseComplexMatrix(rows, columns);
  }

  ComplexVector like1D(int size) => new DenseComplexVector(size);

  void setParts(int row, int column, double re, double im) {
    int idx = rowZero + row * rowStride + columnZero + column * columnStride;
    _elements[idx] = re;
    _elements[idx + 1] = im;
  }

  void set(int row, int column, Complex value) {
    int idx = rowZero + row * rowStride + columnZero + column * columnStride;
    _elements[idx] = value.real;
    _elements[idx + 1] = value.imaginary;
  }

  ComplexVector mult(ComplexVector y, [ComplexVector z = null,
      Complex alpha = null, Complex beta = null, bool transposeA = false]) {
    if (alpha == null) {
      alpha = Complex.ONE;
    }
    if (beta == null) {
      beta = (z == null ? Complex.ONE : Complex.ZERO);
    }
    if (transposeA) {
      return conjugateTranspose().mult(y, z, alpha, beta, false);
    }
    ComplexVector zz;
    if (z == null) {
      zz = new DenseComplexVector(rows);
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
      elemsZ[idxZeroZ] = reS * alpha.real -
          imS * alpha.imaginary +
          reZ * beta.real -
          imZ * beta.imaginary;
      elemsZ[idxZeroZ + 1] = imS * alpha.real +
          reS * alpha.imaginary +
          imZ * beta.real +
          reZ * beta.imaginary;
      idxZero += rowStride;
      idxZeroZ += strideZ;
    }
    return zz;
  }

  ComplexMatrix multiply(ComplexMatrix B, [ComplexMatrix C = null,
      Complex alpha = null, Complex beta = null, bool transposeA = false,
      bool transposeB = false]) {
    if (alpha == null) {
      alpha = Complex.ONE;
    }
    if (beta == null) {
      beta = (C == null ? Complex.ONE : Complex.ZERO);
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
      C = new DenseComplexMatrix(m, p);
    }
    if (C is! DenseComplexMatrix) {
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

    DenseComplexMatrix BB = B as DenseComplexMatrix;
    DenseComplexMatrix CC = C as DenseComplexMatrix;
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
          CElems[iC] = alpha.real * reS -
              alpha.imaginary * imS +
              beta.real * reC -
              beta.imaginary * imC;
          CElems[iC + 1] = alpha.imaginary * reS +
              alpha.real * imS +
              beta.imaginary * reC +
              beta.real * imC;
          iA += rA;
          iC += rC;
        }
        jB += cB;
        jC += cC;
      }
    }
    return C;
  }

  Complex sum() {
    var sum = Complex.ZERO;
    int zero = index(0, 0);
    int idx = zero;
    for (int r = 0; r < rows; r++) {
      for (int i = idx, c = 0; c < columns; c++) {
        sum += new Complex(_elements[i], _elements[i + 1]);
        i += columnStride;
      }
      idx += rowStride;
    }
    return sum;
  }

  bool _haveSharedCellsRaw(ComplexMatrix other) {
    if (other is SelectedDenseComplexMatrix) {
      return _elements == other._elements;
    } else if (other is DenseComplexMatrix) {
      return _elements == other._elements;
    }
    return false;
  }

  int index(int row, int column) {
    return rowZero + row * rowStride + columnZero + column * columnStride;
  }

  ComplexVector _like1D(int size, int zero, int stride) {
    return new DenseComplexVector._internal(
        size, this._elements, zero, stride, false);
  }

  ComplexMatrix _viewSelectionLike(
      Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedDenseComplexMatrix(
        this._elements, rowOffsets, columnOffsets, 0);
  }

  Object clone() {
    return new DenseComplexMatrix._internal(rows, columns, _elements, rowZero,
        columnZero, rowStride, columnStride, !isView);
  }
}

/// Selection view on dense 2-d matrices holding [Complex] elements.
///
/// Objects of this class are typically constructed via `select` methods on
/// some source matrix. The interface introduced in abstract super classes
/// defines everything a user can do. From a user point of view there is
/// nothing special about this class; it presents the same functionality with
/// the same signatures and semantics as its abstract superclass(es) while
/// introducing no additional functionality. Thus, this class need not be
/// visible to users.
class SelectedDenseComplexMatrix extends ComplexMatrix {
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
      : super._(
          rows, columns, rowZero, columnZero, rowStride, columnStride, false) {
    _elements = elements;
    _rowOffsets = rowOffsets;
    _columnOffsets = columnOffsets;
    _offset = offset;
  }

  int columnOffset(int absRank) => _columnOffsets[absRank];

  int rowOffset(int absRank) => _rowOffsets[absRank];

  Complex get(int row, int column) {
    int idxr = rowZero + row * rowStride;
    int idxc = columnZero + column * columnStride;
    return new Complex(
        _elements[_offset + _rowOffsets[idxr] + _columnOffsets[idxc]],
        _elements[_offset + _rowOffsets[idxr] + _columnOffsets[idxc] + 1]);
  }

  // This method is not supported.
  Object get elements {
    throw new UnsupportedError("This method is not supported.");
  }

  bool _haveSharedCellsRaw(ComplexMatrix other) {
    if (other is SelectedDenseComplexMatrix) {
      return _elements == other._elements;
    } else if (other is DenseComplexMatrix) {
      return _elements == other._elements;
    }
    return false;
  }

  int index(int row, int column) {
    return _offset +
        _rowOffsets[rowZero + row * rowStride] +
        _columnOffsets[columnZero + column * columnStride];
  }

  ComplexMatrix like2D(int rows, int columns) {
    return new DenseComplexMatrix(rows, columns);
  }

  ComplexVector like1D(int size) => new DenseComplexVector(size);

  ComplexVector _like1D(int size, int zero, int stride) {
    // never called since row() and column() are overridden
    throw new Error();
  }

  void set(int row, int column, Complex value) {
    int idxr = rowZero + row * rowStride;
    int idxc = columnZero + column * columnStride;
    _elements[_offset + _rowOffsets[idxr] + _columnOffsets[idxc]] = value.real;
    _elements[_offset + _rowOffsets[idxr] + _columnOffsets[idxc] + 1] =
        value.imaginary;
  }

  // This method is not supported.
  ComplexVector vectorize() {
    throw new UnsupportedError("This method is not supported.");
  }

  void setParts(int row, int column, double re, double im) {
    int idxr = rowZero + row * rowStride;
    int idxc = columnZero + column * columnStride;
    _elements[_offset + _rowOffsets[idxr] + _columnOffsets[idxc]] = re;
    _elements[_offset + _rowOffsets[idxr] + _columnOffsets[idxc] + 1] = im;
  }

  ComplexMatrix dice() {
    var v = _view();
    vDice(v);
    Int32List tmp = _rowOffsets;
    _rowOffsets = _columnOffsets;
    _columnOffsets = tmp;
    setIsNoView(v, false);
    return v;
  }

  ComplexVector column(int column) {
    checkColumn(this, column);
    int viewSize = rows;
    int viewZero = rowZero;
    int viewStride = rowStride;
    Int32List viewOffsets = _rowOffsets;
    int viewOffset = _offset + columnOffset(columnRank(column));
    return new SelectedDenseComplexVector._internal(viewSize, this._elements,
        viewZero, viewStride, viewOffsets, viewOffset);
  }

  ComplexVector row(int row) {
    checkRow(this, row);
    int viewSize = columns;
    int viewZero = columnZero;
    int viewStride = columnStride;
    Int32List viewOffsets = _columnOffsets;
    int viewOffset = _offset + rowOffset(rowRank(row));
    return new SelectedDenseComplexVector._internal(
        viewSize, _elements, viewZero, viewStride, viewOffsets, viewOffset);
  }

  ComplexMatrix _viewSelectionLike(
      Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedDenseComplexMatrix(
        this._elements, rowOffsets, columnOffsets, _offset);
  }

  DoubleMatrix real() {
    var R = new DenseDoubleMatrix(rows, columns);
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        var tmp = get(r, c);
        R.set(r, c, tmp.real);
      }
    }
    return R;
  }

  DoubleMatrix imaginary() {
    var Im = new DenseDoubleMatrix(rows, columns);
    for (int r = 0; r < rows; r++) {
      for (int c = 0; c < columns; c++) {
        var tmp = get(r, c);
        Im.set(r, c, tmp.imaginary);
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
