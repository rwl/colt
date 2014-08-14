import 'dart:io';
import 'package:colt/colt.dart';

/**
 * Demonstrates how to use this class.
 */
main() {
    final values = [ [ 5, 10, 20, 40 ], [ 7, 8, 6, 7 ], [ 12, 10, 20, 19 ], [ 3, 1, 5, 6 ] ];
    List<String> columnNames = [ "1996", "1997", "1998", "1999" ];
    List<String> rowNames = [ "PowerBar", "Benzol", "Mercedes", "Sparcling" ];
    String rowAxisName = "CPU";
    String columnAxisName = "Year";
    String title = "CPU performance over time [nops/sec]";
    DoubleBinFunctions1D F = hep.aida.tdouble.bin.DoubleBinFunctions1D.functions;
    List<DoubleBinFunction1D> aggr = [ F.mean, F.rms, F.quantile(0.25), F.median,
            F.quantile(0.75), F.stdDev, F.min, F.max ];
    String format = "%1.2G";
    System.out.println(new DoubleFormatter(format).toTitleString(cern.colt.matrix.tdouble.DoubleFactory2D.dense
            .make(values), rowNames, columnNames, rowAxisName, columnAxisName, title, aggr));
}
