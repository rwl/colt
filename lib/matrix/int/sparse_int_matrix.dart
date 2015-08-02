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
part of cern.colt.matrix.int;

/**
 * Sparse hashed 2-d matrix holding <tt>int</tt> elements. First see the <a
 * href="package-summary.html">package summary</a> and javadoc <a
 * href="package-tree.html">tree view</a> to get the broad picture.
 * <p>
 * <b>Implementation:</b>
 * <p>
 * Note that this implementation is not synchronized. Uses a
 * {@link cern.colt.map.tint.OpenIntIntHashMap}, which is a compact and
 * performant hashing technique.
 * <p>
 * <b>Memory requirements:</b>
 * <p>
 * Cells that
 * <ul>
 * <li>are never set to non-zero values do not use any memory.
 * <li>switch from zero to non-zero state do use memory.
 * <li>switch back from non-zero to zero state also do use memory. However,
 * their memory is automatically reclaimed from time to time. It can also
 * manually be reclaimed by calling {@link #trimToSize()}.
 * </ul>
 * <p>
 * worst case: <tt>memory [bytes] = (1/minLoadFactor) * nonZeros * 13</tt>. <br>
 * best case: <tt>memory [bytes] = (1/maxLoadFactor) * nonZeros * 13</tt>. <br>
 * Where <tt>nonZeros = cardinality()</tt> is the number of non-zero cells.
 * Thus, a 1000 x 1000 matrix with minLoadFactor=0.25 and maxLoadFactor=0.5 and
 * 1000000 non-zero cells consumes between 25 MB and 50 MB. The same 1000 x 1000
 * matrix with 1000 non-zero cells consumes between 25 and 50 KB.
 * <p>
 * <b>Time complexity:</b>
 * <p>
 * This class offers <i>expected</i> time complexity <tt>O(1)</tt> (i.e.
 * constant time) for the basic operations <tt>get</tt>, <tt>getQuick</tt>,
 * <tt>set</tt>, <tt>setQuick</tt> and <tt>size</tt> assuming the hash function
 * disperses the elements properly among the buckets. Otherwise, pathological
 * cases, although highly improbable, can occur, degrading performance to
 * <tt>O(N)</tt> in the worst case. As such this sparse class is expected to
 * have no worse time complexity than its dense counterpart
 * {@link IntMatrix}. However, constant factors are considerably larger.
 * <p>
 * Cells are internally addressed in row-major. Performance sensitive
 * applications can exploit this fact. Setting values in a loop row-by-row is
 * quicker than column-by-column, because fewer hash collisions occur. Thus
 *
 * <pre>
 * for (int row = 0; row &lt; rows; row++) {
 *     for (int column = 0; column &lt; columns; column++) {
 *         matrix.setQuick(row, column, someValue);
 *     }
 * }
 * </pre>
 *
 * is quicker than
 *
 * <pre>
 * for (int column = 0; column &lt; columns; column++) {
 *     for (int row = 0; row &lt; rows; row++) {
 *         matrix.setQuick(row, column, someValue);
 *     }
 * }
 * </pre>
 *
 * @see cern.colt.map
 * @see cern.colt.map.tint.OpenIntIntHashMap
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 09/24/99
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 */
class SparseIntMatrix extends AbstractIntMatrix {

  /*
   * The elements of the matrix.
   */
  Map<int, int> _elements;

  /**
   * Constructs a matrix with a copy of the given values. <tt>values</tt> is
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
  factory SparseIntMatrix.fromList(List<Int32List> values) {
    return new SparseIntMatrix(values.length, values.length == 0 ? 0 : values[0].length)..setAll2D(values);
  }

  /**
   * Constructs a matrix with a given number of rows and columns and default
   * memory usage. All entries are initially <tt>0</tt>.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @throws ArgumentError
   *             if
   *             <tt>rows<0 || columns<0 || (double)columns*rows > Integer.MAX_VALUE</tt>
   *             .
   */
  factory SparseIntMatrix(int rows, int columns) {
    return new SparseIntMatrix._internal(rows, columns, new Map<int, int>(), 0, 0, null, 1);
  }

  /**
   * Constructs a matrix with a given number of rows and columns using memory
   * as specified. All entries are initially <tt>0</tt>. For details related
   * to memory usage see {@link cern.colt.map.tlong.OpenLongIntHashMap}.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @param initialCapacity
   *            the initial capacity of the hash map. If not known, set
   *            <tt>initialCapacity=0</tt> or small.
   * @param minLoadFactor
   *            the minimum load factor of the hash map.
   * @param maxLoadFactor
   *            the maximum load factor of the hash map.
   * @throws ArgumentError
   *             if
   *
   *             <tt>initialCapacity < 0 || (minLoadFactor < 0.0 || minLoadFactor >= 1.0) || (maxLoadFactor <= 0.0 || maxLoadFactor >= 1.0) || (minLoadFactor >= maxLoadFactor)</tt>
   *             .
   * @throws ArgumentError
   *             if
   *             <tt>rows<0 || columns<0 || (double)columns*rows > Integer.MAX_VALUE</tt>
   *             .
   */
  /*SparseIntMatrix(int rows, int columns, int initialCapacity, double minLoadFactor, double maxLoadFactor) {
        try {
            _setUp(rows, columns);
        } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
            if (!"matrix too large".equals(exc.getMessage()))
                throw exc;
        }
        this._elements = new OpenLongIntHashMap(initialCapacity, minLoadFactor, maxLoadFactor);
    }*/

  /**
   * Constructs a matrix with a copy of the given indexes and a single value.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @param rowIndexes
   *            row indexes
   * @param columnIndexes
   *            column indexes
   * @param value
   *            numerical value
   */
  factory SparseIntMatrix.withValue(int rows, int columns, Int32List rowIndexes, Int32List columnIndexes, int value) {
    /*try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if (!"matrix too large".equals(exc.getMessage())) throw exc;
    }
    this._elements = new Map<int, int>();//rowIndexes.length);*/
    return new SparseIntMatrix(rows, columns).._insert(rowIndexes, columnIndexes, value);
  }

  /**
   * Constructs a matrix with a copy of the given indexes and values.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @param rowIndexes
   *            row indexes
   * @param columnIndexes
   *            column indexes
   * @param values
   *            numerical values
   */
  factory SparseIntMatrix.withValues(int rows, int columns, Int32List rowIndexes, Int32List columnIndexes, Int32List values) {
    /*try {
      _setUp(rows, columns);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if (!"matrix too large".equals(exc.getMessage())) throw exc;
    }
    this._elements = new OpenLongIntHashMap(rowIndexes.length);*/
    return new SparseIntMatrix(rows, columns).._insertValues(rowIndexes, columnIndexes, values);
  }

  /**
   * Constructs a matrix from MatrixVectorReader.
   *
   * @param reader
   *            matrix reader
   * @throws IOException
   */
  /*factory SparseIntMatrix.read(MatrixVectorReader reader) {
    MatrixInfo info;
    if (reader.hasInfo()) {
      info = reader.readMatrixInfo();
    } else {
      info = new MatrixInfo(true, MatrixInfo.MatrixField.Real, MatrixInfo.MatrixSymmetry.General);
    }

    if (info.isPattern()) {
      throw new UnsupportedOperationException("Pattern matrices are not supported");
    }
    if (info.isDense()) {
      throw new UnsupportedOperationException("Dense matrices are not supported");
    }
    if (info.isComplex()) {
      throw new UnsupportedOperationException("Complex matrices are not supported");
    }

    MatrixSize size = reader.readMatrixSize(info);
    try {
      _setUp(size.numRows(), size.numColumns());
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if (!"matrix too large".equals(exc.getMessage())) throw exc;
    }
    int numEntries = size.numEntries();
    Int32List columnIndexes = new Int32List(numEntries);
    Int32List rowIndexes = new Int32List(numEntries);
    Int32List values = new Int32List(numEntries);
    reader.readCoordinate(rowIndexes, columnIndexes, values);
    if (info.isSymmetric() || info.isSkewSymmetric()) {
      this._elements = new OpenLongIntHashMap(2 * rowIndexes.length);
    } else {
      this._elements = new OpenLongIntHashMap(rowIndexes.length);
    }
    _insertValues(rowIndexes, columnIndexes, values);

    if (info.isSymmetric()) {
      for (int i = 0; i < numEntries; i++) {
        if (rowIndexes[i] != columnIndexes[i]) {
          put(columnIndexes[i], rowIndexes[i], values[i]);
        }
      }
    } else if (info.isSkewSymmetric()) {
      for (int i = 0; i < numEntries; i++) {
        if (rowIndexes[i] != columnIndexes[i]) {
          put(columnIndexes[i], rowIndexes[i], -values[i]);
        }
      }
    }
  }*/

  /**
   * Constructs a view with the given parameters.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @param elements
   *            the cells.
   * @param rowZero
   *            the position of the first element.
   * @param columnZero
   *            the position of the first element.
   * @param rowStride
   *            the number of elements between two rows, i.e.
   *            <tt>index(i+1,j)-index(i,j)</tt>.
   * @param columnStride
   *            the number of elements between two columns, i.e.
   *            <tt>index(i,j+1)-index(i,j)</tt>.
   * @throws ArgumentError
   *             if
   *             <tt>rows<0 || columns<0 || (double)columns*rows > Integer.MAX_VALUE</tt>
   *             or flip's are illegal.
   */
  SparseIntMatrix._internal(int rows, int columns, Map<int, int> elements, int rowZero, int columnZero, int rowStride, int columnStride)
      : super(rows, columns, rowZero, columnZero, rowStride, columnStride) {
    try {
      _setUp(rows, columns, rowZero, columnZero, rowStride, columnStride);
    } on ArgumentError catch (exc) { // we can hold rows*columns>Integer.MAX_VALUE cells !
      if ("matrix too large" != exc.message) {
        throw exc;
      }
    }
    this._elements = elements;
    this._isNoView = false;
  }

  /*void forEach(ifunc.IntFunction function) {
    if (this._isNoView && function is ifunc.IntMult) { // x[i] = mult*x[i]
      this._elements.forEach((_, int v) => function(v));
    } else {
      super.forEach(function);
    }
    return;
  }*/

  void fill(int value) {
    // overriden for performance only
    if (this._isNoView && value == 0) {
      this._elements.clear();
    } else {
      super.fill(value);
    }
    return;
  }

  void copyFrom(AbstractIntMatrix source) {
    // overriden for performance only
    if (!(source is SparseIntMatrix)) {
      super.copyFrom(source);
      return;
    }
    SparseIntMatrix other = source as SparseIntMatrix;
    if (other == this) {
      return; // nothing to do
    }
    checkShape(other);

    if (this._isNoView && other._isNoView) { // quickest
      this._elements.addAll(other._elements);
      return;
    }
    super.copyFrom(source);
  }

  void forEachWith(final AbstractIntMatrix y, ifunc.IntIntFunction function) {
    if (!this._isNoView) {
      super.forEachWith(y, function);
      return;
    }

    checkShape(y);

    if (function is ifunc.IntPlusMultSecond) { // x[i] = x[i] + alpha*y[i]
      final int alpha = (function as ifunc.IntPlusMultSecond).multiplicator;
      if (alpha == 0) {
        return; // nothing to do
      }
      y.forEachNonZero((int i, int j, int value) {
        set(i, j, get(i, j) + alpha * value);
        return value;
      });
    } else if (function == ifunc.mult) { // x[i] = x[i] * y[i]
      this._elements.forEach((int key, int value) {
        int i = key ~/ _columns;
        int j = key % _columns;
        int r = value * y.get(i, j);
        if (r != value) {
          _elements[key] = r;
        }
        return true;
      });
    } else if (function == ifunc.div) { // x[i] = x[i] / y[i]
      this._elements.forEach((int key, int value) {
        int i = key ~/ _columns;
        int j = key % _columns;
        int r = value ~/ y.get(i, j);
        if (r != value) {
          _elements[key] = r;
        }
        return true;
      });
    } else {
      super.forEachWith(y, function);
    }
    return;
  }

  /**
   * Assigns the result of a function to each cell;
   * <tt>x[row,col] = function(x[row,col],y[row,col])</tt>, where y is given
   * in the coordinate form with single numerical value.
   *
   * @param rowIndexes
   *            row indexes of y
   * @param columnIndexes
   *            column indexes of y
   * @param value
   *            numerical value of y
   * @param function
   *            a function object taking as first argument the current cell's
   *            value of <tt>this</tt>, and as second argument the current
   *            cell's value of <tt>y</tt>,
   * @return <tt>this</tt> (for convenience only).
   */
  SparseIntMatrix assign(final Int32List rowIndexes, final Int32List columnIndexes, final int value, final ifunc.IntIntFunction function) {
    int size = rowIndexes.length;
    if (function == ifunc.plus) { // x[i] = x[i] + y[i]
      for (int i = 0; i < size; i++) {
        int row = rowIndexes[i];
        int column = columnIndexes[i];
        if (row >= _rows || column >= _columns) {
          throw new RangeError("row: $row, column: $column");
        }
        int index = _rowZero + row * _rowStride + _columnZero + column * _columnStride;
        int elem = _elements[index];
        int sum = elem + value;
        if (sum != 0) {
          _elements[index] = sum;
        } else {
          _elements.remove(index);
        }
      }
    } else {
      for (int i = 0; i < size; i++) {
        int row = rowIndexes[i];
        int column = columnIndexes[i];
        if (row >= _rows || column >= _columns) {
          throw new RangeError("row: $row, column: $column");
        }
        int index = _rowZero + row * _rowStride + _columnZero + column * _columnStride;
        int elem = _elements[index];
        int result = function(elem, value);
        if (result != 0) {
          _elements[index] = result;
        } else {
          _elements.remove(index);
        }
      }
    }
    return this;
  }

  /**
   * Assigns the result of a function to each cell;
   * <tt>x[row,col] = function(x[row,col],y[row,col])</tt>, where y is given
   * in the coordinate form.
   *
   * @param rowIndexes
   *            row indexes of y
   * @param columnIndexes
   *            column indexes of y
   * @param values
   *            numerical values of y
   * @param function
   *            a function object taking as first argument the current cell's
   *            value of <tt>this</tt>, and as second argument the current
   *            cell's value of <tt>y</tt>,
   * @return <tt>this</tt> (for convenience only).
   */
  SparseIntMatrix setValues(final Int32List rowIndexes, final Int32List columnIndexes, final Int32List values, final ifunc.IntIntFunction function) {
    int size = rowIndexes.length;
    if (function == ifunc.plus) { // x[i] = x[i] + y[i]
      for (int i = 0; i < size; i++) {
        int value = values[i];
        int row = rowIndexes[i];
        int column = columnIndexes[i];
        if (row >= _rows || column >= _columns) {
          throw new RangeError("row: $row, column: $column");
        }
        int index = _rowZero + row * _rowStride + _columnZero + column * _columnStride;
        int elem = _elements[index];
        value += elem;
        if (value != 0) {
          _elements[index] = value;
        } else {
          _elements.remove(index);
        }
      }
    } else {
      for (int i = 0; i < size; i++) {
        int value = values[i];
        int row = rowIndexes[i];
        int column = columnIndexes[i];
        if (row >= _rows || column >= _columns) {
          throw new RangeError("row: $row, column: $column");
        }
        int index = _rowZero + row * _rowStride + _columnZero + column * _columnStride;
        int elem = _elements[index];
        value = function(elem, value);
        if (value != 0) {
          _elements[index] = value;
        } else {
          _elements.remove(index);
        }
      }
    }
    return this;
  }

  int get cardinality {
    if (this._isNoView) {
      return this._elements.length;
    } else {
      return super.cardinality;
    }
  }

  /**
   * Returns a new matrix that has the same elements as this matrix, but is in
   * a column-compressed form. This method creates a new object (not a view),
   * so changes in the returned matrix are NOT reflected in this matrix.
   *
   * @param sortRowIndexes
   *            if true, then row indexes in column compressed matrix are
   *            sorted
   *
   * @return this matrix in a column-compressed form
   */
  SparseCCIntMatrix columnCompressed([bool sortRowIndexes=false]) {
    int nnz = cardinality;
    Int32List keys = _elements.keys;
    Int32List values = _elements.values;
    Int32List rowIndexes = new Int32List(nnz);
    Int32List columnIndexes = new Int32List(nnz);

    for (int k = 0; k < nnz; k++) {
      int key = keys[k];
      rowIndexes[k] = key ~/ _columns;
      columnIndexes[k] = key % _columns;
    }
    return new SparseCCIntMatrix.withValues(_rows, _columns, rowIndexes, columnIndexes, values, false, false, sortRowIndexes);
  }

  /**
   * Returns a new matrix that has the same elements as this matrix, but is in
   * a column-compressed modified form. This method creates a new object (not
   * a view), so changes in the returned matrix are NOT reflected in this
   * matrix.
   *
   * @return this matrix in a column-compressed modified form
   */
  /*SparseCCMIntMatrix columnCompressedModified() {
    SparseCCMIntMatrix A = new SparseCCMIntMatrix(_rows, _columns);
    int nnz = cardinality;
    Int32List keys = _elements.keys;
    Int32List values = _elements.values;
    for (int i = 0; i < nnz; i++) {
      int row = keys[i] ~/ _columns;
      int column = keys[i] % _columns;
      A.set(row, column, values[i]);
    }
    return A;
  }*/

  /**
   * Returns a new matrix that has the same elements as this matrix, but is in
   * a row-compressed form. This method creates a new object (not a view), so
   * changes in the returned matrix are NOT reflected in this matrix.
   *
   * @param sortColumnIndexes
   *            if true, then column indexes in row compressed matrix are
   *            sorted
   *
   * @return this matrix in a row-compressed form
   */
  SparseRCIntMatrix rowCompressed([bool sortColumnIndexes=false]) {
    int nnz = cardinality;
    Int32List keys = _elements.keys;
    Int32List values = _elements.values;
    final Int32List rowIndexes = new Int32List(nnz);
    final Int32List columnIndexes = new Int32List(nnz);
    for (int k = 0; k < nnz; k++) {
      int key = keys[k];
      rowIndexes[k] = key ~/ _columns;
      columnIndexes[k] = key % _columns;
    }
    return new SparseRCIntMatrix.withValues(_rows, _columns, rowIndexes, columnIndexes, values, false, false, sortColumnIndexes);
  }

  /**
   * Returns a new matrix that has the same elements as this matrix, but is in
   * a row-compressed modified form. This method creates a new object (not a
   * view), so changes in the returned matrix are NOT reflected in this
   * matrix.
   *
   * @return this matrix in a row-compressed modified form
   */
  /*SparseRCMIntMatrix rowCompressedModified() {
    SparseRCMIntMatrix A = new SparseRCMIntMatrix(_rows, _columns);
    int nnz = cardinality;
    Int32List keys = _elements.keys;
    Int32List values = _elements.values;
    for (int i = 0; i < nnz; i++) {
      int row = keys[i] ~/ _columns;
      int column = keys[i] % _columns;
      A.set(row, column, values[i]);
    }
    return A;
  }*/

  Object get elements => _elements;

  /*void ensureCapacity(int minCapacity) {
    this._elements.ensureCapacity(minCapacity);
  }*/

  AbstractIntMatrix forEachNonZero(final ifunc.IntIntIntFunction function) {
    if (this._isNoView) {
      this._elements.forEach((int key, int value) {
        int i = key ~/ _columns;
        int j = key % _columns;
        int r = function(i, j, value);
        if (r != value) {
          _elements[key] = r;
        }
        return true;
      });
    } else {
      super.forEachNonZero(function);
    }
    return this;
  }

  int get(int row, int column) {
    return this._elements[_rowZero + row * _rowStride + _columnZero + column * _columnStride];
  }

  int index(int row, int column) {
    return _rowZero + row * _rowStride + _columnZero + column * _columnStride;
  }

  AbstractIntMatrix like2D(int rows, int columns) {
    return new SparseIntMatrix(rows, columns);
  }

  AbstractIntVector like1D(int size) {
    return new SparseIntVector(size);
  }

  void set(int row, int column, int value) {
    int index = _rowZero + row * _rowStride + _columnZero + column * _columnStride;
    if (value == 0) {
      this._elements.remove(index);
    } else {
      this._elements[index] = value;
    }
  }

  String toString() {
    StringBuffer builder = new StringBuffer();
    builder
        ..write(_rows)
        ..write(" x ")
        ..write(_columns)
        ..write(" sparse matrix, nnz = ")
        ..write(cardinality)
        ..write('\n');
    for (int r = 0; r < _rows; r++) {
      for (int c = 0; c < _columns; c++) {
        int elem = get(r, c);
        if (elem != 0) {
          builder
              ..write('(')
              ..write(r)
              ..write(',')
              ..write(c)
              ..write(')')
              ..write('\t')
              ..write(elem)
              ..write('\n');
        }
      }
    }
    return builder.toString();
  }

  /*void trimToSize() {
    this._elements.trimToSize();
  }*/

  AbstractIntVector vectorize() {
    SparseIntVector v = new SparseIntVector(length);
    int idx = 0;
    for (int c = 0; c < _columns; c++) {
      for (int r = 0; r < _rows; r++) {
        int elem = get(r, c);
        v.set(idx++, elem);
      }
    }
    return v;
  }

  AbstractIntVector mult(AbstractIntVector y, [AbstractIntVector z = null, final int alpha = 1, int beta = null, final bool transposeA = false]) {
    if (beta == null) {
      beta = z == null ? 1 : 0;
    }
    int rowsA = _rows;
    int columnsA = _columns;
    if (transposeA) {
      rowsA = _columns;
      columnsA = _rows;
    }

    bool ignore = (z == null);
    if (z == null) {
      z = new IntVector(rowsA);
    }

    if (!(this._isNoView && y is IntVector && z is IntVector)) {
      return super.mult(y, z, alpha, beta, transposeA);
    }

    if (columnsA != y.length || rowsA > z.length) {
      throw new ArgumentError("Incompatible args: " + ((transposeA ? dice() : this).toStringShort()) + ", " + y.toStringShort() + ", " + z.toStringShort());
    }

    if (!ignore) {
      z.forEach(ifunc.multiply(beta));
    }

    IntVector zz = z as IntVector;
    final Int32List elementsZ = zz._elements;
    final int strideZ = zz.stride();
    final int zeroZ = z.index(0);

    IntVector yy = y as IntVector;
    final Int32List elementsY = yy._elements;
    final int strideY = yy.stride();
    final int zeroY = y.index(0);

    if (elementsY == null || elementsZ == null) {
      throw new Error();
    }

    this._elements.forEach((int key, int value) {
      int i = key ~/ _columns;
      int j = key % _columns;
      if (transposeA) {
        int tmp = i;
        i = j;
        j = tmp;
      }
      elementsZ[zeroZ + strideZ * i] += alpha * value * elementsY[zeroY + strideY * j];
      return true;
    });

    return z;
  }

  AbstractIntMatrix multiply(AbstractIntMatrix B, [AbstractIntMatrix C = null, final int alpha = 1, int beta = null, final bool transposeA = false, final bool transposeB = false]) {
    if (beta == null) {
      beta = C == null ? 1 : 0;
    }
    if (!(this._isNoView)) {
      return super.multiply(B, C, alpha, beta, transposeA, transposeB);
    }
    if (transposeB) B = B.dice();
    int rowsA = _rows;
    int columnsA = _columns;
    if (transposeA) {
      rowsA = _columns;
      columnsA = _rows;
    }
    int p = B.columns;
    bool ignore = (C == null);
    if (C == null) C = new IntMatrix(rowsA, p);

    if (B.rows != columnsA) {
      throw new ArgumentError("Matrix2D inner dimensions must agree:" + toStringShort() + ", " + (transposeB ? B.dice() : B).toStringShort());
    }
    if (C.rows != rowsA || C.columns != p) {
      throw new ArgumentError("Incompatibel result matrix: " + toStringShort() + ", " + (transposeB ? B.dice() : B).toStringShort() + ", " + C.toStringShort());
    }
    if (this == C || B == C) throw new ArgumentError("Matrices must not be identical");

    if (!ignore) {
      C.forEach(ifunc.multiply(beta));
    }

    // cache views
    final List<AbstractIntVector> Brows = new List<AbstractIntVector>(columnsA);
    for (int i = columnsA; --i >= 0; ) {
      Brows[i] = B.row(i);
    }
    final List<AbstractIntVector> Crows = new List<AbstractIntVector>(rowsA);
    for (int i = rowsA; --i >= 0; ) {
      Crows[i] = C.row(i);
    }

    final ifunc.IntPlusMultSecond fun = new ifunc.IntPlusMultSecond.plusMult(0);

    this._elements.forEach((int key, int value) {
      int i = key ~/ _columns;
      int j = key % _columns;
      fun.multiplicator = value * alpha;
      if (!transposeA) {
        Crows[i].forEachWith(Brows[j], fun);
      } else {
        Crows[j].forEachWith(Brows[i], fun);
      }
      return true;
    });

    return C;
  }

  void _insert(Int32List rowIndexes, Int32List columnIndexes, int value) {
    int size = rowIndexes.length;
    for (int i = 0; i < size; i++) {
      int row = rowIndexes[i];
      int column = columnIndexes[i];
      if (row >= _rows || column >= _columns) {
        throw new RangeError("row: $row, column: $column");
      }
      if (value != 0) {
        int index = _rowZero + row * _rowStride + _columnZero + column * _columnStride;
        int elem = _elements[index];
        if (elem == 0) {
          _elements[index] = value;
        } else {
          int sum = elem + value;
          if (sum == 0) {
            _elements.remove(index);
          } else {
            _elements[index] = sum;
          }
        }
      }
    }
  }

  void _insertValues(Int32List rowIndexes, Int32List columnIndexes, Int32List values) {
    int size = rowIndexes.length;
    for (int i = 0; i < size; i++) {
      int value = values[i];
      int row = rowIndexes[i];
      int column = columnIndexes[i];
      if (row >= _rows || column >= _columns) {
        throw new RangeError("row: $row, column: $column");
      }
      if (value != 0) {
        int index = _rowZero + row * _rowStride + _columnZero + column * _columnStride;
        int elem = _elements[index];
        if (elem == 0) {
          _elements[index] = value;
        } else {
          int sum = elem + value;
          if (sum == 0) {
            _elements.remove(index);
          } else {
            _elements[index] = sum;
          }
        }
      }
    }
  }

  bool _haveSharedCellsRaw(AbstractIntMatrix other) {
    if (other is SelectedSparseIntMatrix) {
      return this._elements == other._elements;
    } else if (other is SparseIntMatrix) {
      return this._elements == other._elements;
    }
    return false;
  }

  AbstractIntVector _like1D(int size, int offset, int stride) {
    return new SparseIntVector._internal(size, this._elements, offset, stride);
  }

  AbstractIntMatrix _viewSelectionLike(Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedSparseIntMatrix(this._elements, rowOffsets, columnOffsets, 0);
  }

  Object clone() {
    return new SparseIntMatrix._internal(rows, columns, _elements, _rowZero, _columnZero, _rowStride, _columnStride);
  }

}


/**
 * Selection view on sparse 2-d matrices holding <tt>int</tt> elements. First
 * see the <a href="package-summary.html">package summary</a> and javadoc <a
 * href="package-tree.html">tree view</a> to get the broad picture.
 * <p>
 * <b>Implementation:</b>
 * <p>
 * Objects of this class are typically constructed via <tt>viewIndexes</tt>
 * methods on some source matrix. The interface introduced in abstract super
 * classes defines everything a user can do. From a user point of view there is
 * nothing special about this class; it presents the same functionality with the
 * same signatures and semantics as its abstract superclass(es) while
 * introducing no additional functionality. Thus, this class need not be visible
 * to users. By the way, the same principle applies to concrete DenseXXX and
 * SparseXXX classes: they presents the same functionality with the same
 * signatures and semantics as abstract superclass(es) while introducing no
 * additional functionality. Thus, they need not be visible to users, either.
 * Factory methods could hide all these concrete types.
 * <p>
 * This class uses no delegation. Its instances point directly to the data. Cell
 * addressing overhead is 1 additional int addition and 2 additional array index
 * accesses per get/set.
 * <p>
 * Note that this implementation is not synchronized.
 * <p>
 * <b>Memory requirements:</b>
 * <p>
 * <tt>memory [bytes] = 4*(rowIndexes.length+columnIndexes.length)</tt>. Thus,
 * an index view with 1000 x 1000 indexes additionally uses 8 KB.
 * <p>
 * <b>Time complexity:</b>
 * <p>
 * Depends on the parent view holding cells.
 * <p>
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 09/24/99
 */
class SelectedSparseIntMatrix extends AbstractIntMatrix {

  /*
   * The elements of the matrix.
   */
  Map<int, int> _elements;

  /**
   * The offsets of the visible cells of this matrix.
   */
  Int32List _rowOffsets;

  Int32List _columnOffsets;

  /**
   * The offset.
   */
  int _offset;

  /**
   * Constructs a matrix view with the given parameters.
   *
   * @param elements
   *            the cells.
   * @param rowOffsets
   *            The row offsets of the cells that shall be visible.
   * @param columnOffsets
   *            The column offsets of the cells that shall be visible.
   * @param offset
   */
  factory SelectedSparseIntMatrix(Map<int, int> elements, Int32List rowOffsets, Int32List columnOffsets, int offset) {
    return new SelectedSparseIntMatrix._internal(rowOffsets.length, columnOffsets.length, elements, 0, 0, 1, 1, rowOffsets, columnOffsets, offset);
  }

  /**
   * Constructs a matrix view with the given parameters.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @param elements
   *            the cells.
   * @param rowZero
   *            the position of the first element.
   * @param columnZero
   *            the position of the first element.
   * @param rowStride
   *            the number of elements between two rows, i.e.
   *            <tt>index(i+1,j)-index(i,j)</tt>.
   * @param columnStride
   *            the number of elements between two columns, i.e.
   *            <tt>index(i,j+1)-index(i,j)</tt>.
   * @param rowOffsets
   *            The row offsets of the cells that shall be visible.
   * @param columnOffsets
   *            The column offsets of the cells that shall be visible.
   * @param offset
   */
  SelectedSparseIntMatrix._internal(int rows, int columns, Map<int, int> elements, int rowZero, int columnZero, int rowStride, int columnStride, Int32List rowOffsets, Int32List columnOffsets, int offset)
      : super(rows, columns, rowZero, columnZero, rowStride, columnStride) {
    // be sure parameters are valid, we do not check...
    _setUp(rows, columns, rowZero, columnZero, rowStride, columnStride);

    this._elements = elements;
    this._rowOffsets = rowOffsets;
    this._columnOffsets = columnOffsets;
    this._offset = offset;

    this._isNoView = false;
  }

  Object get elements {
    throw new ArgumentError("This method is not supported.");
  }

  /**
   * Returns the matrix cell value at coordinate <tt>[row,column]</tt>.
   *
   * <p>
   * Provided with invalid parameters this method may return invalid objects
   * without throwing any exception. <b>You should only use this method when
   * you are absolutely sure that the coordinate is within bounds.</b>
   * Precondition (unchecked):
   * <tt>0 &lt;= column &lt; columns() && 0 &lt;= row &lt; rows()</tt>.
   *
   * @param row
   *            the index of the row-coordinate.
   * @param column
   *            the index of the column-coordinate.
   * @return the value at the specified coordinate.
   */
  int get(int row, int column) {
    // if (debug) if (column<0 || column>=columns || row<0 || row>=rows)
    // throw new IndexOutOfBoundsException("row:"+row+", column:"+column);
    // return elements.get(index(row,column));
    // manually inlined:
    return _elements[_offset + _rowOffsets[_rowZero + row * _rowStride] + _columnOffsets[_columnZero + column * _columnStride]];
  }

  /**
   * Returns the position of the given coordinate within the (virtual or
   * non-virtual) internal 1-dimensional array.
   *
   * @param row
   *            the index of the row-coordinate.
   * @param column
   *            the index of the column-coordinate.
   */
  int index(int row, int column) {
    // return this.offset + super.index(row,column);
    // manually inlined:
    return this._offset + _rowOffsets[_rowZero + row * _rowStride] + _columnOffsets[_columnZero + column * _columnStride];
  }

  /**
   * Construct and returns a new empty matrix <i>of the same dynamic type</i>
   * as the receiver, having the specified number of rows and columns. For
   * example, if the receiver is an instance of type <tt>DenseIntMatrix2D</tt>
   * the new matrix must also be of type <tt>DenseIntMatrix2D</tt>, if the
   * receiver is an instance of type <tt>SparseIntMatrix2D</tt> the new matrix
   * must also be of type <tt>SparseIntMatrix2D</tt>, etc. In general, the new
   * matrix should have internal parametrization as similar as possible.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @return a new empty matrix of the same dynamic type.
   */
  AbstractIntMatrix like2D(int rows, int columns) {
    return new SparseIntMatrix(rows, columns);
  }

  /**
   * Construct and returns a new 1-d matrix <i>of the corresponding dynamic
   * type</i>, entirelly independent of the receiver. For example, if the
   * receiver is an instance of type <tt>DenseIntMatrix2D</tt> the new matrix
   * must be of type <tt>DenseIntMatrix1D</tt>, if the receiver is an instance
   * of type <tt>SparseIntMatrix2D</tt> the new matrix must be of type
   * <tt>SparseIntMatrix1D</tt>, etc.
   *
   * @param size
   *            the number of cells the matrix shall have.
   * @return a new matrix of the corresponding dynamic type.
   */
  AbstractIntVector like1D(int size) {
    return new SparseIntVector(size);
  }

  /**
   * Sets the matrix cell at coordinate <tt>[row,column]</tt> to the specified
   * value.
   *
   * <p>
   * Provided with invalid parameters this method may access illegal indexes
   * without throwing any exception. <b>You should only use this method when
   * you are absolutely sure that the coordinate is within bounds.</b>
   * Precondition (unchecked):
   * <tt>0 &lt;= column &lt; columns() && 0 &lt;= row &lt; rows()</tt>.
   *
   * @param row
   *            the index of the row-coordinate.
   * @param column
   *            the index of the column-coordinate.
   * @param value
   *            the value to be filled into the specified cell.
   */
  void set(int row, int column, int value) {
    // if (debug) if (column<0 || column>=columns || row<0 || row>=rows)
    // throw new IndexOutOfBoundsException("row:"+row+", column:"+column);
    // int index = index(row,column);
    // manually inlined:
    int index = _offset + _rowOffsets[_rowZero + row * _rowStride] + _columnOffsets[_columnZero + column * _columnStride];

    if (value == 0) {
      this._elements.remove(index);
    } else {
      this._elements[index] = value;
    }
  }

  /**
   * Returns a vector obtained by stacking the columns of the matrix on top of
   * one another.
   *
   * @return
   */
  AbstractIntVector vectorize() {
    SparseIntVector v = new SparseIntVector(length);
    int idx = 0;
    for (int c = 0; c < _columns; c++) {
      for (int r = 0; r < _rows; r++) {
        v.set(idx++, get(c, r));
      }
    }
    return v;
  }

  /**
   * Constructs and returns a new <i>slice view</i> representing the rows of
   * the given column. The returned view is backed by this matrix, so changes
   * in the returned view are reflected in this matrix, and vice-versa. To
   * obtain a slice view on subranges, construct a sub-ranging view (
   * <tt>viewPart(...)</tt>), then apply this method to the sub-range view.
   * <p>
   * <b>Example:</b>
   * <table border="0">
   * <tr nowrap>
   * <td valign="top">2 x 3 matrix: <br>
   * 1, 2, 3<br>
   * 4, 5, 6</td>
   * <td>viewColumn(0) ==></td>
   * <td valign="top">Matrix1D of size 2:<br>
   * 1, 4</td>
   * </tr>
   * </table>
   *
   * @param the
   *            column to fix.
   * @return a new slice view.
   * @throws IllegalArgumentException
   *             if <tt>column < 0 || column >= columns()</tt>.
   * @see #row(int)
   */
  AbstractIntVector column(int column) {
    _checkColumn(column);
    int viewSize = this._rows;
    int viewZero = this._rowZero;
    int viewStride = this._rowStride;
    Int32List viewOffsets = this._rowOffsets;
    int viewOffset = this._offset + _columnOffset(_columnRank(column));
    return new SelectedSparseIntVector._internal(viewSize, this._elements, viewZero, viewStride, viewOffsets, viewOffset);
  }

  /**
   * Constructs and returns a new <i>slice view</i> representing the columns
   * of the given row. The returned view is backed by this matrix, so changes
   * in the returned view are reflected in this matrix, and vice-versa. To
   * obtain a slice view on subranges, construct a sub-ranging view (
   * <tt>viewPart(...)</tt>), then apply this method to the sub-range view.
   * <p>
   * <b>Example:</b>
   * <table border="0">
   * <tr nowrap>
   * <td valign="top">2 x 3 matrix: <br>
   * 1, 2, 3<br>
   * 4, 5, 6</td>
   * <td>viewRow(0) ==></td>
   * <td valign="top">Matrix1D of size 3:<br>
   * 1, 2, 3</td>
   * </tr>
   * </table>
   *
   * @param the
   *            row to fix.
   * @return a new slice view.
   * @throws IndexOutOfBoundsException
   *             if <tt>row < 0 || row >= rows()</tt>.
   * @see #column(int)
   */
  AbstractIntVector row(int row) {
    _checkRow(row);
    int viewSize = this._columns;
    int viewZero = _columnZero;
    int viewStride = this._columnStride;
    Int32List viewOffsets = this._columnOffsets;
    int viewOffset = this._offset + _rowOffset(_rowRank(row));
    return new SelectedSparseIntVector._internal(viewSize, this._elements, viewZero, viewStride, viewOffsets, viewOffset);
  }

  /**
   * Returns the position of the given absolute rank within the (virtual or
   * non-virtual) internal 1-dimensional array. Default implementation.
   * Override, if necessary.
   *
   * @param rank
   *            the absolute rank of the element.
   * @return the position.
   */
  int _columnOffset(int absRank) {
    return _columnOffsets[absRank];
  }

  /**
   * Returns the position of the given absolute rank within the (virtual or
   * non-virtual) internal 1-dimensional array. Default implementation.
   * Override, if necessary.
   *
   * @param rank
   *            the absolute rank of the element.
   * @return the position.
   */
  int _rowOffset(int absRank) {
    return _rowOffsets[absRank];
  }

  /**
   * Returns <tt>true</tt> if both matrices share common cells. More formally,
   * returns <tt>true</tt> if <tt>other != null</tt> and at least one of the
   * following conditions is met
   * <ul>
   * <li>the receiver is a view of the other matrix
   * <li>the other matrix is a view of the receiver
   * <li><tt>this == other</tt>
   * </ul>
   */
  bool _haveSharedCellsRaw(AbstractIntMatrix other) {
    if (other is SelectedSparseIntMatrix) {
      return this._elements == other._elements;
    } else if (other is SparseIntMatrix) {
      return this._elements == other._elements;
    }
    return false;
  }

  /**
   * Construct and returns a new 1-d matrix <i>of the corresponding dynamic
   * type</i>, sharing the same cells. For example, if the receiver is an
   * instance of type <tt>DenseIntMatrix2D</tt> the new matrix must be of type
   * <tt>DenseIntMatrix1D</tt>, if the receiver is an instance of type
   * <tt>SparseIntMatrix2D</tt> the new matrix must be of type
   * <tt>SparseIntMatrix1D</tt>, etc.
   *
   * @param size
   *            the number of cells the matrix shall have.
   * @param zero
   *            the index of the first element.
   * @param stride
   *            the number of indexes between any two elements, i.e.
   *            <tt>index(i+1)-index(i)</tt>.
   * @return a new matrix of the corresponding dynamic type.
   */
  AbstractIntVector _like1D(int size, int zero, int stride) {
    // this method is never called since
    // row() and column() are overridden properly.
    throw new Error();
  }

  /**
   * Sets up a matrix with a given number of rows and columns.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @throws IllegalArgumentException
   *             if <tt>(double)columns*rows > Integer.MAX_VALUE</tt>.
   */
  /*void _setUp(int rows, int columns) {
    super._setUp(rows, columns);
    this._rowStride = 1;
    this._columnStride = 1;
    this._offset = 0;
  }*/

  /**
   * Self modifying version of viewDice().
   */
  void _vDice() {
    super._vDice();
    // swap
    Int32List tmp = _rowOffsets;
    _rowOffsets = _columnOffsets;
    _columnOffsets = tmp;

    // flips stay unaffected

    this._isNoView = false;
  }

  /**
   * Construct and returns a new selection view.
   *
   * @param rowOffsets
   *            the offsets of the visible elements.
   * @param columnOffsets
   *            the offsets of the visible elements.
   * @return a new view.
   */
  AbstractIntMatrix _viewSelectionLike(Int32List rowOffsets, Int32List columnOffsets) {
    return new SelectedSparseIntMatrix(this._elements, rowOffsets, columnOffsets, this._offset);
  }

  Object clone() {
    return new SelectedSparseIntMatrix._internal(rows, columns, _elements, _rowZero, _columnZero, _rowStride, _columnStride, _rowOffsets, _columnOffsets, _offset);
  }
}
