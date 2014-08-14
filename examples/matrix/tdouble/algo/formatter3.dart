import 'dart:io';
import 'package:colt/colt.dart';

/**
 * Demonstrates how to use this class.
 */
main(int size, double value) {
  Timer timer = new Timer();
  String s;
  StringBuffer buf;
  DoubleMatrix2D matrix = DoubleFactory2D.dense.make(size, size, value);

  timer.reset().start();
  buf = new StringBuffer();
  for (int i = size; --i >= 0; ) {
    for (int j = size; --j >= 0; ) {
      buf.append(matrix.getQuick(i, j));
    }
  }
  buf = null;
  timer.stop().display();

  timer.reset().start();
  Former format = new FormerFactory().create("%G");
  buf = new StringBuffer();
  for (int i = size; --i >= 0; ) {
    for (int j = size; --j >= 0; ) {
      buf.append(format.form(matrix.getQuick(i, j)));
    }
  }
  buf = null;
  timer.stop().display();

  timer.reset().start();
  s = new DoubleFormatter(null).toString(matrix);
  System.out.println(s);
  s = null;
  timer.stop().display();

  timer.reset().start();
  s = new DoubleFormatter("%G").toString(matrix);
  System.out.println(s);
  s = null;
  timer.stop().display();
}
