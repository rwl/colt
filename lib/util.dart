library edu.emory.mathcs.util;

import 'dart:async';
import 'function/double.dart';

typedef void Computation();

class ConcurrencyUtils {

  static int _NTHREADS = getNumberOfProcessors();

  static int _THREADS_BEGIN_N_1D = 32768;

  static int _THREADS_BEGIN_N_2D = 65536;

  static int _THREADS_BEGIN_N_3D = 65536;

  /**
   * Returns the number of available processors.
   */
  static int getNumberOfProcessors() {
    return 1;
  }

  /**
   * Returns the current number of threads.
   *
   * @return the current number of threads.
   */
  static int getNumberOfThreads() {
    return _NTHREADS;
  }

  /**
   * Returns the minimal size of 1D data for which threads are used.
   */
  static int getThreadsBeginN_1D() {
    return _THREADS_BEGIN_N_1D;
  }

  static Future<double> submitDouble(Computation f) {
    return new Future(f);
  }

  static Future submit(Computation f) {
    return new Future(f);
  }

  /**
   * Waits for all threads to complete computation and aggregates the result.
   *
   * @param futures
   *            handles to running threads
   * @param aggr
   *            an aggregation function
   * @return the result of aggregation
   */
  static double waitForCompletionAggr(List<Future<double>> futures, DoubleDoubleFunction aggr) {
    int size = futures.length;
    List<double> results = new List<double>(size);
    double a = 0.0;
    //try {
    /*for (int j = 0; j < size; j++) {
                //results[j] = futures[j].get() as double;
                futures[j].then((double value) {
                  results[j] = value;
                });
            }*/
    Future.wait(futures).then((List<double> values) {
      results = values;
    });
    a = results[0];
    for (int j = 1; j < size; j++) {
      a = aggr(a, results[j]);
    }
    /*} on ExecutionException catch (ex) {
            ex.printStackTrace();
        } on InterruptedException catch (e) {
            e.printStackTrace();
        }*/
    return a;
  }

  /**
   * Waits for all threads to complete computation.
   *
   * @param futures
   *            handles to running threads
   */
  static void waitForCompletion(List<Future> futures) {
    Future.wait(futures);
  }
}
