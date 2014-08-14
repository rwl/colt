import 'dart:io';
import 'package:colt/colt.dart';

/**
 * Demonstrates how to use this class.
 */
main() {
    // parameters
    final values = [ [ 3, 0, -3.4, 0 ], [ 5.1, 0, 3.0123456789, 0 ], [ 16.37, 0.0, 2.5, 0 ],
            [ -16.3, 0, -3.012345678E-4, -1 ], [ 1236.3456789, 0, 7, -1.2 ] ];
    /*
     * double[][] values = [ [3, 1, ], [5.1 ,16.37, ] ];
     */
    // List<String> columnNames = [ "he", "", "he", "four" ];
    // List<String> rowNames = [ "hello", "du", null, "abcdef", "five" ];
    List<String> columnNames = [ "0.1", "0.3", "0.5", "0.7" ];
    List<String> rowNames = [ "SunJDK1.2.2 classic", "IBMJDK1.1.8", "SunJDK1.3 Hotspot", "other1", "other2" ];
    // List<String> columnNames = [ "0.1", "0.3" ];
    // List<String> rowNames = [ "SunJDK1.2.2 classic", "IBMJDK1.1.8"];

    stdout.writeln(cern.colt.matrix.tdouble.DoubleFactory2D.dense.make(values));
    stdout.writeln(new DoubleFormatter("%G")._toTitleString(DoubleFactory2D.dense
            .make(values), rowNames, columnNames, "vendor", "density", "title"));
}