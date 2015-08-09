//import 'dart:io';
//import 'package:colt/colt.dart';

main() {
//  // parameters
//  final values = [// 5, 0.0, -0.0, -double.NaN, double.NaN, 0.0/0.0,
//    // double.NEGATIVE_INFINITY, double.POSITIVE_INFINITY,
//    // double.MIN_VALUE,
//    // double.MAX_VALUE
//    5, 0.0, -0.0, -double.NAN, double.NAN, 0.0 / 0.0, double.MIN_POSITIVE,
//    double.MAX_FINITE, double.NEGATIVE_INFINITY, double.INFINITY
//    // double.MIN_VALUE, double.MAX_VALUE //, double.NEGATIVE_INFINITY,
//    // double.POSITIVE_INFINITY
//  ];
//  // List<String> formats = {"%G", "%1.10G", "%f", "%1.2f", "%0.2e"};
//  List<String> formats = ["%G", "%1.19G"];
//
//  // now the processing
//  int size = formats.length;
//  DoubleVector matrix = new DenseDoubleVector(values);
//
//  List<String> strings = new List<String>(size);
//  // List<String> javaStrings = new String[size];
//
//  for (int i = 0; i < size; i++) {
//    String format = formats[i];
//    strings[i] = new DoubleFormatter(format).toString(matrix);
//    for (int j = 0; j < matrix.length; j++) {
//      stdout.writeln(matrix.get(j).toString());
//    }
//  }
//
//  stdout.writeln("original:\n" + new DoubleFormatter().toString(matrix));
//
//  for (int i = 0; i < size; i++) {
//    stdout.writeln("\nstring(" + formats[i] + "):\n" + strings[i]);
//  }
}
