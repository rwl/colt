import 'dart:io';
import 'package:colt/colt.dart';

void main() {
  stdout.writeln("\n\n");
  DoubleMatrix2D matrix;
  DoubleMatrix2D A, B, C, D;
  DoubleMatrix2D _ = null;

  A = make(2, 2, 1);
  B = make(4, 4, 2);
  C = make(4, 3, 3);
  D = make(2, 2, 4);
  List<List<DoubleMatrix2D>> parts1 = [[_, A, _], [B, _, C], [_, D, _]];
  matrix = DoubleFactory2D.compose(parts1);
  stdout.writeln("\n" + matrix);

  A.fill(9);
  B.forEach(9);
  C.forEach(9);
  D.forEach(9);
  DoubleFactory2D.decompose(parts1, matrix);
  stdout.writeln(A);
  stdout.writeln(B);
  stdout.writeln(C);
  stdout.writeln(D);
  // stdout.writeln("\n"+cern.colt.matrixpattern.Converting.toHTML(make(parts1).toString()));

  /*
   * // illegal 2 != 3 DoubleMatrix2D[][] parts2 = [ [ null, make(2,2,1),
   * null ], [ make(4,4,2), null, make(4,3,3) ], [ null, make(2,3,4), null ] ];
   * stdout.writeln("\n"+Factory2D.make(parts2));
   */

  /*
   * DoubleMatrix2D[][] parts3 = [ [ identity(3), null, ], [ null,
   * identity(3).viewColumnFlip() ], [ identity(3).viewRowFlip(), null ] ];
   * stdout.writeln("\n"+make(parts3));
   * //stdout.writeln("\n"+cern.colt.matrixpattern.Converting.toHTML(make(parts3).toString()));
   *
   * DoubleMatrix2D A = ascending(2,2); DoubleMatrix2D B =
   * descending(2,2); DoubleMatrix2D _ = null;
   *
   * DoubleMatrix2D[][] parts4 = [ [ A, _, A, _ ], [ _, A, _, B ] ];
   * stdout.writeln("\n"+make(parts4));
   * //stdout.writeln("\n"+cern.colt.matrixpattern.Converting.toHTML(make(parts4).toString()));
   */
}
