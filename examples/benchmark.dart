import 'dart:typed_data';
import 'package:benchmark_harness/benchmark_harness.dart';
import 'package:complex/complex.dart';

Float64List conj(Float64List a) {
  var ret = new Float64List(2);
  ret[0] = a[0];
  ret[1] = -a[1];
  return ret;
}

class ComplexBenchmark extends BenchmarkBase {
  const ComplexBenchmark() : super("Complex");

  static void main() {
    new ComplexBenchmark().report();
  }

  void run() {
    var a = new Complex(3.0, 4.0);
    var b = a.conjugate();
    var _ = b.real;
  }
}

class TypedDataBenchmark extends BenchmarkBase {
  const TypedDataBenchmark() : super("TypedData");

  static void main() {
    new TypedDataBenchmark().report();
  }

  void run() {
    var a = new Float64List(2);
    a[0] = 3.0;
    a[1] = 4.0;
//    var a = new Float64List.fromList([3.0, 4.0]);
    var c = conj(a);
    var _ = c[0];
  }
}

main() {
  ComplexBenchmark.main();
  TypedDataBenchmark.main();
}