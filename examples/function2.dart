

void main(int size) {
//  //cern.jet.math.tdouble.DoubleFunctions F = cern.jet.math.tdouble.DoubleFunctions.functions;
//  print("\n\n");
//  double a = 0.0;
//  double b = 0.0;
//  double v = (Math.sin(a) + Math.pow(Math.cos(b), 2)).abs();
//  // double v = Math.sin(a) + Math.pow(Math.cos(b),2);
//  // double v = a + b;
//  print(v);
//
//  // DoubleDoubleFunction f = F.chain(F.plus,F.identity,F.identity);
//  DoubleDoubleFunction f = chain(abs, chainFGH(plus, sin, chainGH(square, cos)));
//  // DoubleDoubleFunction f =
//  // F.chain(F.plus,F.sin,F.chain(F.square,F.cos));
//  // DoubleDoubleFunction f = F.plus;
//
//  print(f(a, b));
//  double g(double x, double y) {
//    return (Math.sin(x) + Math.pow(Math.cos(y), 2)).abs();
//  };
//  print(g(a, b));
//
//  // emptyLoop
//  //colt.Timer emptyLoop = new colt.Timer().start();
//  a = 0.0;
//  b = 0.0;
//  double sum = 0.0;
//  for (int i = size; --i >= 0;) {
//    sum += a;
//    a++;
//    b++;
//  }
//  //emptyLoop.stop().display();
//  print("empty sum=$sum");
//
//  //colt.Timer timer = new colt.Timer().start();
//  a = 0.0;
//  b = 0.0;
//  sum = 0.0;
//  for (int i = size; --i >= 0;) {
//    sum += (Math.sin(a) + Math.pow(Math.cos(b), 2)).abs();
//    // sum += a + b;
//    a++;
//    b++;
//  }
//  //timer.stop().display();
//  //print("evals / sec = ${size / timer.minus(emptyLoop).seconds()}");
//  print("sum=$sum");
//
//  //timer.reset().start();
//  a = 0.0;
//  b = 0.0;
//  sum = 0.0;
//  for (int i = size; --i >= 0;) {
//    sum += f(a, b);
//    a++;
//    b++;
//  }
//  //timer.stop().display();
//  //print("evals / sec = ${size / timer.minus(emptyLoop).seconds()}");
//  print("sum=${sum}");
//
//  //timer.reset().start();
//  a = 0.0;
//  b = 0.0;
//  sum = 0.0;
//  for (int i = size; --i >= 0;) {
//    sum += g(a, b);
//    a++;
//    b++;
//  }
//  //timer.stop().display();
//  //print("evals / sec = ${size / timer.minus(emptyLoop).seconds()}");
//  print("sum=$sum");

}
