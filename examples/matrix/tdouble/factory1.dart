import 'dart:io';
import 'package:colt/colt.dart';

void main() {
  stdout.writeln("\n\n");
  List<List<DoubleMatrix>> parts1 = [[null, make(2, 2, 1), null], [make(4, 4, 2), null, make(4, 3, 3)],
                                       [null, make(2, 2, 4), null]];
  stdout.writeln("\n" + DoubleFactory2D.dense.compose(parts1).toString());
  // stdout.writeln("\n"+matrixpattern.Converting.toHTML(make(parts1).toString()));

  /*
   * // illegal 2 != 3 DoubleMatrix[][] parts2 = [ [ null, make(2,2,1),
   * null ], [ make(4,4,2), null, make(4,3,3) ], [ null, make(2,3,4), null ] ];
   * stdout.writeln("\n"+make(parts2));
   */

  List<List<DoubleMatrix>> parts3 = [[identity(3), null,], [null, identity(3).columnFlip()],
                                       [identity(3).rowFlip(), null]];
  stdout.writeln("\n" + DoubleFactory2D.compose(parts3));
  // stdout.writeln("\n"+matrixpattern.Converting.toHTML(make(parts3).toString()));

  DoubleMatrix A = DoubleFactory2D.dense.ascending(2, 2);
  DoubleMatrix B = DoubleFactory2D.dense.descending(2, 2);
  DoubleMatrix _ = null;

  List<List<DoubleMatrix>> parts4 = [[A, _, A, _], [_, A, _, B]];
  stdout.writeln("\n" + DoubleFactory2D.dense.compose(parts4).toString());
  // stdout.writeln("\n"+matrixpattern.Converting.toHTML(make(parts4).toString()));

}
