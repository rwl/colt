/*
Copyright (C) 1999 CERN - European Organization for Nuclear Research.
Permission to use, copy, modify, distribute and sell this software and its documentation for any purpose
is hereby granted without fee, provided that the above copyright notice appear in all copies and
that both that copyright notice and this permission notice appear in supporting documentation.
CERN makes no representations about the suitability of this software for any purpose.
It is provided "as is" without expressed or implied warranty.
 */
library cern.jet.math;

/**
 * Only for performance tuning of compute intensive linear algebraic
 * computations. Constructs functions that return one of
 * <ul>
 * <li><tt>a * constant</tt>
 * <li><tt>a / constant</tt>
 * </ul>
 * <tt>a</tt> is variable, <tt>constant</tt> is fixed, but for performance
 * reasons publicly accessible. Intended to be passed to
 * <tt>matrix.assign(function)</tt> methods.
 */
class DoubleMult {
  /**
   * Public read/write access to avoid frequent object construction.
   */
  double multiplicator;

  /**
   * Insert the method's description here. Creation date: (8/10/99 19:12:09)
   */
  DoubleMult(final double multiplicator) {
    this.multiplicator = multiplicator;
  }

  /**
   * Returns the result of the function evaluation.
   */
  double call(double a) {
    return a * multiplicator;
  }

  /**
   * <tt>a / constant</tt>.
   */
  factory DoubleMult.div(final double constant) {
    return new DoubleMult.mult(1 / constant);
  }

  /**
   * <tt>a * constant</tt>.
   */
  factory DoubleMult.mult(final double constant) {
    return new DoubleMult(constant);
  }
}

/**
 * Only for performance tuning of compute intensive linear algebraic
 * computations. Constructs functions that return one of
 * <ul>
 * <li><tt>a*constant + b</tt>
 * <li><tt>a*constant - b</tt>
 * <li><tt>a/constant + b</tt>
 * <li><tt>a/constant - b</tt>
 * </ul>
 * <tt>a</tt> and <tt>b</tt> are variables, <tt>constant</tt> is fixed, but for
 * performance reasons publicly accessible. Intended to be passed to
 * <tt>matrix.assign(otherMatrix,function)</tt> methods.
 */
class DoublePlusMultFirst {
  /**
   * Public read/write access to avoid frequent object construction.
   */
  double multiplicator;

  DoublePlusMultFirst(final double multiplicator) {
    this.multiplicator = multiplicator;
  }

  /**
   * Returns the result of the function evaluation.
   */
  double call(double a, double b) {
    return a * multiplicator + b;
  }

  /**
   * <tt>a - b/constant</tt>.
   */
  factory DoublePlusMultFirst.minusDiv(final double constant) {
    return new DoublePlusMultFirst(-1 / constant);
  }

  /**
   * <tt>a - b*constant</tt>.
   */
  factory DoublePlusMultFirst.minusMult(final double constant) {
    return new DoublePlusMultFirst(-constant);
  }

  /**
   * <tt>a + b/constant</tt>.
   */
  factory DoublePlusMultFirst.plusDiv(final double constant) {
    return new DoublePlusMultFirst(1 / constant);
  }

  /**
   * <tt>a + b*constant</tt>.
   */
  factory DoublePlusMultFirst.plusMult(final double constant) {
    return new DoublePlusMultFirst(constant);
  }
}

/**
 * Only for performance tuning of compute intensive linear algebraic
 * computations. Constructs functions that return one of
 * <ul>
 * <li><tt>a + b*constant</tt>
 * <li><tt>a - b*constant</tt>
 * <li><tt>a + b/constant</tt>
 * <li><tt>a - b/constant</tt>
 * </ul>
 * <tt>a</tt> and <tt>b</tt> are variables, <tt>constant</tt> is fixed, but for
 * performance reasons publicly accessible. Intended to be passed to
 * <tt>matrix.assign(otherMatrix,function)</tt> methods.
 */
class DoublePlusMultSecond {//implements cern.colt.function.tdouble.DoubleDoubleFunction {
  /**
   * Public read/write access to avoid frequent object construction.
   */
  double multiplicator;

  /**
   * Insert the method's description here. Creation date: (8/10/99 19:12:09)
   */
  DoublePlusMultSecond(final double multiplicator) {
    this.multiplicator = multiplicator;
  }

  /**
   * Returns the result of the function evaluation.
   */
  double call(double a, double b) {
    return a + b * multiplicator;
  }

  /**
   * <tt>a - b/constant</tt>.
   */
  factory DoublePlusMultSecond.minusDiv(final double constant) {
    return new DoublePlusMultSecond(-1 / constant);
  }

  /**
   * <tt>a - b*constant</tt>.
   */
  factory DoublePlusMultSecond.minusMult(final double constant) {
    return new DoublePlusMultSecond(-constant);
  }

  /**
   * <tt>a + b/constant</tt>.
   */
  factory DoublePlusMultSecond.plusDiv(final double constant) {
    return new DoublePlusMultSecond(1 / constant);
  }

  /**
   * <tt>a + b*constant</tt>.
   */
  factory DoublePlusMultSecond.plusMult(final double constant) {
    return new DoublePlusMultSecond(constant);
  }
}
