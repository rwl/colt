import 'dart:typed_data';
import 'package:colt/colt.dart';

main() {
  var values = new Float64List.fromList([
    5.0,
    0.0,
    -0.0,
    -double.NAN,
    double.NAN,
    0.0 / 0.0,
    double.MIN_POSITIVE,
    double.MAX_FINITE,
    double.NEGATIVE_INFINITY,
    double.INFINITY
  ]);
  var formats = ["%G", "%1.10G", "%f", "%1.2f", "%0.2e"];

  var vector = new DenseDoubleVector.fromList(values);

  for (var format in formats) {
    print(new DoubleFormatter(format).toStringDoubleVector(vector));
  }
}
