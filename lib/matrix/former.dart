/*
Copyright (C) 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
is hereby granted without fee, provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear in supporting documentation.
CERN makes no representations about the suitability of this software for any purpose.
It is provided "as is" without expressed or implied warranty.
 */
library cern.colt.matrix.format;

import 'dart:math' as Math;
import 'package:intl/intl.dart';
import 'matrix.dart';

part 'former_factory.dart';
part 'abstract_formatter.dart';
part 'double/algo/formatter.dart';

/**
 * Formats a double or complex (double[]) into a string (like sprintf in C).
 *
 * @author wolfgang.hoschek@cern.ch
 * @author Piotr Wendykier (piotr.wendykier@gmail.com)
 * @version 1.0, 21/07/00
 * @see java.util.Comparator
 * @see cern.colt
 * @see cern.colt.Sorting
 */
abstract class Former {
  /**
   * Formats a double into a string (like sprintf in C).
   *
   * @param value
   *            the number to format
   * @return the formatted string
   * @exception IllegalArgumentException
   *                if bad argument
   */
  String formDouble(double value);

  /**
   * Formats a float into a string (like sprintf in C).
   *
   * @param value
   *            the number to format
   * @return the formatted string
   * @exception IllegalArgumentException
   *                if bad argument
   */
  //String form(float value);

  /**
   * Formats an int into a string (like sprintf in C).
   *
   * @param value
   *            the number to format
   * @return the formatted string
   * @exception IllegalArgumentException
   *                if bad argument
   */
  String formInt(int value);

  /**
   * Formats an long into a string (like sprintf in C).
   *
   * @param value
   *            the number to format
   * @return the formatted string
   * @exception IllegalArgumentException
   *                if bad argument
   */
  //String form(long value);

  /**
   * Formats a complex (double[]) into a string (like sprintf in C).
   *
   * @param value
   *            the number to format
   * @return the formatted string
   * @exception IllegalArgumentException
   *                if bad argument
   */
  String formComplex(List<double> value);

  /**
   * Formats a complex (float[]) into a string (like sprintf in C).
   *
   * @param value
   *            the number to format
   * @return the formatted string
   * @exception IllegalArgumentException
   *                if bad argument
   */
  //String form(List<float> value);

}
