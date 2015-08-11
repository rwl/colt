# Colt

A [Dart][] library for vector and matrix algebra.

Includes fixed sized (non-resizable) dense and sparse matrices for `int`,
`double` and `Complex` types.

A Dart translation of [ParallelColt][] which in turn is based on the
original [Colt][] project from CERN.

### Example

```dart
import 'package:colt/colt.dart';

main() {
  var c = new ComplexVector(17);
  c.fill(3.0, 4.0);
  print(c.abs());
}
```

## np

An interface to Colt resembling [NumPy][].

### Example

```dart
import 'package:colt/np.dart' as np;

main() {
  var a = np.array([1, 2, 3]);
  var b = a * np.sin(a);
}
```

[Dart]: https://www.dartlang.org/
[ParallelColt]: https://sites.google.com/site/piotrwendykier/software/parallelcolt
[Colt]: http://dst.lbl.gov/ACSSoftware/colt/
[NumPy]: http://www.numpy.org/
