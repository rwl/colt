part of cern.colt.matrix;

/**
 * Factory producing implementations of {@link cern.colt.matrix.Former} via
 * method create(); Serves to isolate the interface of String formatting from
 * the actual implementation. If you want to plug in a different String
 * formatting implementation, simply replace this class with your alternative.
 *
 * @author wolfgang.hoschek@cern.ch
 * @version 1.0, 21/07/00
 */
class FormerFactory {
  /**
   * Constructs and returns a new format instance.
   *
   * @param format
   *            the format string following printf conventions. The string has
   *            a prefix, a format code and a suffix. The prefix and suffix
   *            become part of the formatted output. The format code directs
   *            the formatting of the (single) parameter to be formatted. The
   *            code has the following structure
   *            <ul>
   *            <li> a % (required)
   *            <li> a modifier (optional)
   *            <dl>
   *            <dt> +
   *            <dd> forces display of + for positive numbers
   *            <dt> 0
   *            <dd> show leading zeroes
   *            <dt> -
   *            <dd> align left in the field
   *            <dt> space
   *            <dd> prepend a space in front of positive numbers
   *            <dt> #
   *            <dd> use "alternate" format. Add 0 or 0x for octal or
   *            hexadecimal numbers. Don't suppress trailing zeroes in general
   *            floating point format.
   *            </dl>
   *            <li> an integer denoting field width (optional)
   *            <li> a period followed by an integer denoting precision
   *            (optional)
   *            <li> a format descriptor (required)
   *            <dl>
   *            <dt>f
   *            <dd> floating point number in fixed format
   *            <dt>e, E
   *            <dd> floating point number in exponential notation (scientific
   *            format). The E format results in an uppercase E for the
   *            exponent (1.14130E+003), the e format in a lowercase e.
   *            <dt>g, G
   *            <dd> floating point number in general format (fixed format for
   *            small numbers, exponential format for large numbers). Trailing
   *            zeroes are suppressed. The G format results in an uppercase E
   *            for the exponent (if any), the g format in a lowercase e.
   *            <dt>d, i
   *            <dd> integer in decimal
   *            <dt>x
   *            <dd> integer in hexadecimal
   *            <dt>o
   *            <dd> integer in octal
   *            <dt>s
   *            <dd> string
   *            <dt>c
   *            <dd> character
   *            </dl>
   *            </ul>
   * @exception IllegalArgumentException
   *                if bad format
   */
  Former create(final String format) {
    return new DefaultFormer(format);
  }
}

class DefaultFormer implements Former {
  NumberFormat format;

  DefaultFormer(String fmt) {
    format = new NumberFormat(fmt);
  }

  String formDouble(double value) {
    return format.format(value);
  }

  /*String form(float value) {
      return String.format(format, value);
  }*/

  String formInt(int value) {
    return format.format(value);
  }

  /*String form(long value) {
      return String.format(format, value);
  }*/

  String formComplex(List<double> value) {
    if (value[0] == 0 && value[1] == 0) {
      return "0";
    }
    if (value[1] == 0) {
      return format.format(value[0]);
    }
    if (value[1] < 0) {
      return format.format(value[0]) + " - " + format.format(-value[1]);
    }
    return format.format(value[0]) + " + " + format.format(value[1]);
  }

  /*String form(float[] value) {
      if (value[0] == 0 && value[1] == 0) {
          return "0";
      }
      if (value[1] == 0) {
          return String.format(format, value[0]);
      }
      if (value[1] < 0) {
          return String.format(format, value[0]) + " - " + String.format(format, -value[1]);
      }
      return String.format(format, value[0]) + " + " + String.format(format, value[1]);
  }*/
}
