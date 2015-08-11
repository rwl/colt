import 'dart:typed_data';
import 'package:colt/colt.dart';

main() {
  final values = [
    3.0,
    0.0,
    -3.4,
    0.0,
    5.1,
    0.0,
    3.0123456789,
    0.0,
    16.37,
    0.0,
    2.5,
    0.0,
    -16.3,
    0.0,
    -3.012345678E-4,
    -1.0,
    1236.3456789,
    0.0,
    7.0,
    -1.2
  ];
  List<String> formats = ["%G", "%1.10G", "%f", "%1.2f", "%0.2e", null];

  var matrix = new DoubleMatrix(5, 4)..setAll(new Float64List.fromList(values));

  for (var format in formats) {
    print(new DoubleFormatter(format).toStringMatrix(matrix));
  }
}
