/*
Copyright (C) 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
is hereby granted without fee, provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear in supporting documentation.
CERN makes no representations about the suitability of this software for any purpose.
It is provided "as is" without expressed or implied warranty.
 */
part of cern.colt.matrix;

/**
 * Sparse hashed 1-d matrix (aka <i>vector</i>) holding <tt>int</tt> elements.
 * First see the <a href="package-summary.html">package summary</a> and javadoc
 * <a href="package-tree.html">tree view</a> to get the broad picture.
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
 * Thus, a 1000000 matrix with minLoadFactor=0.25 and maxLoadFactor=0.5 and
 * 1000000 non-zero cells consumes between 25 MB and 50 MB. The same 1000000
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
 * {@link DenseIntMatrix1D}. However, constant factors are considerably larger.
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 09/24/99
 *
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @version 1.1, 08/22/2007
 */
class SparseIntVector extends AbstractIntVector {
  /*
   * The elements of the matrix.
   */
  Map<int, int> _elements;

  /**
   * Constructs a matrix with a copy of the given values. The values are
   * copied. So subsequent changes in <tt>values</tt> are not reflected in the
   * matrix, and vice-versa.
   *
   * @param values
   *            The values to be filled into the new matrix.
   */
  factory SparseIntVector.fromList(List<int> values) {
    return new SparseIntVector(values.length)
      ..setAll(0, values);
  }

  /**
   * Constructs a matrix with a given number of cells. All entries are
   * initially <tt>0</tt>.
   *
   * @param size
   *            the number of cells the matrix shall have.
   * @throws IllegalArgumentException
   *             if <tt>size<0</tt>.
   */
  factory SparseIntVector(int size) {
    return new SparseIntVector._internal(size, new Map<int, int>(), 0, 1);
  }

  /**
   * Constructs a matrix with a given number of parameters. All entries are
   * initially <tt>0</tt>. For details related to memory usage see
   * {@link cern.colt.map.tlong.OpenLongIntHashMap}.
   *
   * @param size
   *            the number of cells the matrix shall have.
   * @param initialCapacity
   *            the initial capacity of the hash map. If not known, set
   *            <tt>initialCapacity=0</tt> or small.
   * @param minLoadFactor
   *            the minimum load factor of the hash map.
   * @param maxLoadFactor
   *            the maximum load factor of the hash map.
   * @throws IllegalArgumentException
   *             if
   *
   *             <tt>initialCapacity < 0 || (minLoadFactor < 0.0 || minLoadFactor >= 1.0) || (maxLoadFactor <= 0.0 || maxLoadFactor >= 1.0) || (minLoadFactor >= maxLoadFactor)</tt>
   *             .
   * @throws IllegalArgumentException
   *             if <tt>size<0</tt>.
   */
//  SparseIntMatrix1D(int size, int initialCapacity, double minLoadFactor, double maxLoadFactor) {
//    _setUp(size);
//    this._elements = new OpenLongIntHashMap(initialCapacity, minLoadFactor, maxLoadFactor);
//  }

  /**
   * Constructs a matrix view with a given number of parameters.
   *
   * @param size
   *            the number of cells the matrix shall have.
   * @param elements
   *            the cells.
   * @param offset
   *            the index of the first element.
   * @param stride
   *            the number of indexes between any two elements, i.e.
   *            <tt>index(i+1)-index(i)</tt>.
   * @throws IllegalArgumentException
   *             if <tt>size<0</tt>.
   */
  SparseIntVector._internal(int size, Map<int, int> elements, int offset, int stride) {
    _setUp(size, offset, stride);
    this._elements = elements;
    this._isNoView = false;
  }

  /**
   * Sets all cells to the state specified by <tt>value</tt>.
   *
   * @param value
   *            the value to be filled into the cells.
   * @return <tt>this</tt> (for convenience only).
   */
  void fill(int value) {
    if (this._isNoView && value == 0) {
      this._elements.clear();
    } else {
      super.fill(value);
    }
  }

  /**
   * Returns the number of cells having non-zero values.
   */
  int get cardinality {
    if (this._isNoView) {
      return this._elements.length;
    } else {
      return super.cardinality;
    }
  }

  /**
   * Returns the elements of this matrix.
   *
   * @return the elements
   */
  Object get elements => _elements;

  /**
   * Ensures that the receiver can hold at least the specified number of
   * non-zero cells without needing to allocate new internal memory. If
   * necessary, allocates new internal memory and increases the capacity of
   * the receiver.
   * <p>
   * This method never need be called; it is for performance tuning only.
   * Calling this method before tt>set()</tt>ing a large number of non-zero
   * values boosts performance, because the receiver will grow only once
   * instead of potentially many times and hash collisions get less probable.
   *
   * @param minCapacity
   *            the desired minimum number of non-zero cells.
   */
  /*void ensureCapacity(int minCapacity) {
    this._elements.ensureCapacity(minCapacity);
  }*/

  /**
   * Returns the matrix cell value at coordinate <tt>index</tt>.
   *
   * <p>
   * Provided with invalid parameters this method may return invalid objects
   * without throwing any exception. <b>You should only use this method when
   * you are absolutely sure that the coordinate is within bounds.</b>
   * Precondition (unchecked): <tt>index&lt;0 || index&gt;=size()</tt>.
   *
   * @param index
   *            the index of the cell.
   * @return the value of the specified cell.
   */
  int get(int index) {
    // if (debug) if (index<0 || index>=size) checkIndex(index);
    // return this.elements.get(index(index));
    // manually inlined:
    final i = _zero + index * _stride;
    if (_elements.containsKey(i)) {
      return _elements[i];
    } else {
      return 0;
    }
  }

  /**
   * Returns the position of the element with the given relative rank within
   * the (virtual or non-virtual) internal 1-dimensional array. You may want
   * to override this method for performance.
   *
   * @param rank
   *            the rank of the element.
   */
  int index(int rank) {
    // overriden for manual inlining only
    // return __offset(_rank(rank));
    return _zero + rank * _stride;
  }

  /**
   * Construct and returns a new empty matrix <i>of the same dynamic type</i>
   * as the receiver, having the specified size. For example, if the receiver
   * is an instance of type <tt>DenseIntMatrix1D</tt> the new matrix must also
   * be of type <tt>DenseIntMatrix1D</tt>, if the receiver is an instance of
   * type <tt>SparseIntMatrix1D</tt> the new matrix must also be of type
   * <tt>SparseIntMatrix1D</tt>, etc. In general, the new matrix should have
   * internal parametrization as similar as possible.
   *
   * @param size
   *            the number of cell the matrix shall have.
   * @return a new empty matrix of the same dynamic type.
   */
  AbstractIntVector like1D(int size) {
    return new SparseIntVector(size);
  }

  /**
   * Construct and returns a new 2-d matrix <i>of the corresponding dynamic
   * type</i>, entirelly independent of the receiver. For example, if the
   * receiver is an instance of type <tt>DenseIntMatrix1D</tt> the new matrix
   * must be of type <tt>DenseIntMatrix2D</tt>, if the receiver is an instance
   * of type <tt>SparseIntMatrix1D</tt> the new matrix must be of type
   * <tt>SparseIntMatrix2D</tt>, etc.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @return a new matrix of the corresponding dynamic type.
   */
  AbstractIntMatrix like2D(int rows, int columns) {
    return new SparseIntMatrix(rows, columns);
  }

  AbstractIntMatrix reshape(int rows, int columns) {
    if (rows * columns != _size) {
      throw new ArgumentError("rows*columns != size");
    }
    final M = new SparseIntMatrix(rows, columns);
    int idx = 0;
    for (int c = 0; c < columns; c++) {
      for (int r = 0; r < rows; r++) {
        int elem = get(idx++);
        if (elem != 0) {
          M.set(r, c, elem);
        }
      }
    }
    return M;
  }

  /*IntMatrix3D reshape3D(int slices, int rows, int columns) {
    if (slices * rows * columns != _size) {
      throw new IllegalArgumentException("slices*rows*columns != size");
    }
    IntMatrix3D M = new SparseIntMatrix3D(slices, rows, columns);
    int idx = 0;
    for (int s = 0; s < slices; s++) {
      for (int c = 0; c < columns; c++) {
        for (int r = 0; r < rows; r++) {
          int elem = get(idx++);
          if (elem != 0) {
            M.setQuick(s, r, c, elem);
          }
        }
      }
    }
    return M;
  }*/

  /**
   * Sets the matrix cell at coordinate <tt>index</tt> to the specified value.
   *
   * <p>
   * Provided with invalid parameters this method may access illegal indexes
   * without throwing any exception. <b>You should only use this method when
   * you are absolutely sure that the coordinate is within bounds.</b>
   * Precondition (unchecked): <tt>index&lt;0 || index&gt;=size()</tt>.
   *
   * @param index
   *            the index of the cell.
   * @param value
   *            the value to be filled into the specified cell.
   */
  void set(int index, int value) {
    // if (debug) if (index<0 || index>=size) checkIndex(index);
    // int i = index(index);
    // manually inlined:
    int i = _zero + index * _stride;
    if (value == 0) {
      this._elements.remove(i);
    } else {
      this._elements[i] = value;
    }
  }

  String toString() {
    StringBuffer builder = new StringBuffer();
    builder
        ..write("1 x ")
        ..write(_size)
        ..write(" sparse matrix, nnz = ")
        ..write(cardinality)
        ..write('\n');
    for (int i = 0; i < _size; i++) {
      int elem = get(i);
      if (elem != 0) {
        builder
            ..write('(')
            ..write(i)
            ..write(')')
            ..write('\t')
            ..write(elem)
            ..write('\n');
      }
    }
    return builder.toString();
  }

  /*void trimToSize() {
    this._elements.trimToSize();
  }*/

  /**
   * Returns <tt>true</tt> if both matrices share at least one identical cell.
   */
  bool _haveSharedCellsRaw(AbstractIntVector other) {
    if (other is SelectedSparseIntVector) {
      return this._elements == other._elements;
    } else if (other is SparseIntVector) {
      return this._elements == other._elements;
    }
    return false;
  }

  /**
   * Construct and returns a new selection view.
   *
   * @param offsets
   *            the offsets of the visible elements.
   * @return a new view.
   */
  AbstractIntVector _viewSelectionLike(List<int> offsets) {
    return new SelectedSparseIntVector(this._elements, offsets);
  }
}

/**
 * Selection view on sparse 1-d matrices holding <tt>int</tt> elements. First
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
 * to users. By the way, the same principle applies to concrete DenseXXX,
 * SparseXXX classes: they presents the same functionality with the same
 * signatures and semantics as abstract superclass(es) while introducing no
 * additional functionality. Thus, they need not be visible to users, either.
 * Factory methods could hide all these concrete types.
 * <p>
 * This class uses no delegation. Its instances point directly to the data. Cell
 * addressing overhead is 1 additional array index access per get/set.
 * <p>
 * Note that this implementation is not synchronized.
 * <p>
 * <b>Memory requirements:</b>
 * <p>
 * <tt>memory [bytes] = 4*indexes.length</tt>. Thus, an index view with 1000
 * indexes additionally uses 4 KB.
 * <p>
 * <b>Time complexity:</b>
 * <p>
 * Depends on the parent view holding cells.
 * <p>
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 09/24/99
 */
class SelectedSparseIntVector extends AbstractIntVector {

  /*
   * The elements of the matrix.
   */
  Map<int, int> _elements;

  /**
   * The offsets of visible indexes of this matrix.
   */
  List<int> _offsets;

  /**
   * The offset.
   */
  int __offset;

  /**
   * Constructs a matrix view with the given parameters.
   *
   * @param elements
   *            the cells.
   * @param indexes
   *            The indexes of the cells that shall be visible.
   */
  factory SelectedSparseIntVector(Map<int, int> elements, List<int> offsets) {
    return new SelectedSparseIntVector._internal(offsets.length, elements, 0, 1, offsets, 0);
  }

  /**
   * Constructs a matrix view with the given parameters.
   *
   * @param size
   *            the number of cells the matrix shall have.
   * @param elements
   *            the cells.
   * @param zero
   *            the index of the first element.
   * @param stride
   *            the number of indexes between any two elements, i.e.
   *            <tt>index(i+1)-index(i)</tt>.
   * @param offsets
   *            the offsets of the cells that shall be visible.
   * @param offset
   */
  SelectedSparseIntVector._internal(int size, Map<int, int> elements, int zero, int stride, List<int> offsets, int offset) {
    _setUp(size, zero, stride);

    this._elements = elements;
    this._offsets = offsets;
    this.__offset = offset;
    this._isNoView = false;
  }

  Object get elements {
    throw new ArgumentError("This method is not supported.");
  }

  /**
   * Returns the matrix cell value at coordinate <tt>index</tt>.
   *
   * <p>
   * Provided with invalid parameters this method may return invalid objects
   * without throwing any exception. <b>You should only use this method when
   * you are absolutely sure that the coordinate is within bounds.</b>
   * Precondition (unchecked): <tt>index&lt;0 || index&gt;=size()</tt>.
   *
   * @param index
   *            the index of the cell.
   * @return the value of the specified cell.
   */
  int get(int index) {
    // if (debug) if (index<0 || index>=size) checkIndex(index);
    // return elements.get(index(index));
    // manually inlined:
    return _elements[__offset + _offsets[_zero + index * _stride]];
  }

  /**
   * Returns the position of the element with the given relative rank within
   * the (virtual or non-virtual) internal 1-dimensional array. You may want
   * to override this method for performance.
   *
   * @param rank
   *            the rank of the element.
   */
  int index(int rank) {
    // return this.offset + super.index(rank);
    // manually inlined:
    return __offset + _offsets[_zero + rank * _stride];
  }

  /**
   * Construct and returns a new empty matrix <i>of the same dynamic type</i>
   * as the receiver, having the specified size. For example, if the receiver
   * is an instance of type <tt>DenseIntMatrix1D</tt> the new matrix must also
   * be of type <tt>DenseIntMatrix1D</tt>, if the receiver is an instance of
   * type <tt>SparseIntMatrix1D</tt> the new matrix must also be of type
   * <tt>SparseIntMatrix1D</tt>, etc. In general, the new matrix should have
   * internal parametrization as similar as possible.
   *
   * @param size
   *            the number of cell the matrix shall have.
   * @return a new empty matrix of the same dynamic type.
   */
  AbstractIntVector like1D(int size) {
    return new SparseIntVector(size);
  }

  /**
   * Construct and returns a new 2-d matrix <i>of the corresponding dynamic
   * type</i>, entirelly independent of the receiver. For example, if the
   * receiver is an instance of type <tt>DenseIntMatrix1D</tt> the new matrix
   * must be of type <tt>DenseIntMatrix2D</tt>, if the receiver is an instance
   * of type <tt>SparseIntMatrix1D</tt> the new matrix must be of type
   * <tt>SparseIntMatrix2D</tt>, etc.
   *
   * @param rows
   *            the number of rows the matrix shall have.
   * @param columns
   *            the number of columns the matrix shall have.
   * @return a new matrix of the corresponding dynamic type.
   */
  AbstractIntMatrix like2D(int rows, int columns) {
    return new SparseIntMatrix(rows, columns);
  }

  AbstractIntMatrix reshape(int rows, int columns) {
    throw new ArgumentError("This method is not supported.");
  }

  /*IntMatrix3D reshape3D(int slices, int rows, int columns) {
    throw new IllegalArgumentException("This method is not supported.");
  }*/

  /**
   * Sets the matrix cell at coordinate <tt>index</tt> to the specified value.
   *
   * <p>
   * Provided with invalid parameters this method may access illegal indexes
   * without throwing any exception. <b>You should only use this method when
   * you are absolutely sure that the coordinate is within bounds.</b>
   * Precondition (unchecked): <tt>index&lt;0 || index&gt;=size()</tt>.
   *
   * @param index
   *            the index of the cell.
   * @param value
   *            the value to be filled into the specified cell.
   */
  void set(int index, int value) {
    // if (debug) if (index<0 || index>=size) checkIndex(index);
    // int i = index(index);
    // manually inlined:
    int i = __offset + _offsets[_zero + index * _stride];
    if (value == 0) {
      this._elements.remove(i);
    } else {
      this._elements[i] = value;
    }
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
  int _offset(int absRank) {
    return _offsets[absRank];
  }

  /**
   * Returns <tt>true</tt> if both matrices share at least one identical cell.
   */
  bool _haveSharedCellsRaw(AbstractIntVector other) {
    if (other is SelectedSparseIntVector) {
      return this._elements == other._elements;
    } else if (other is SparseIntVector) {
      return this._elements == other._elements;
    }
    return false;
  }

  /**
   * Sets up a matrix with a given number of cells.
   *
   * @param size
   *            the number of cells the matrix shall have.
   */
  /*void _setUp(int size) {
    super._setUp(size);
    this._stride = 1;
    this.__offset = 0;
  }*/

  /**
   * Construct and returns a new selection view.
   *
   * @param offsets
   *            the offsets of the visible elements.
   * @return a new view.
   */
  AbstractIntVector _viewSelectionLike(List<int> offsets) {
    return new SelectedSparseIntVector(this._elements, offsets);
  }
}
