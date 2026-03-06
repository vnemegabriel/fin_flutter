import 'dart:math';

/// Complex number arithmetic.
class Complex {
  final double re;
  final double im;

  const Complex(this.re, this.im);
  const Complex.real(double re) : this(re, 0.0);
  static const Complex zero = Complex(0.0, 0.0);
  static const Complex one = Complex(1.0, 0.0);
  static const Complex i = Complex(0.0, 1.0);

  Complex operator +(Complex other) => Complex(re + other.re, im + other.im);
  Complex operator -(Complex other) => Complex(re - other.re, im - other.im);
  Complex operator -() => Complex(-re, -im);

  Complex operator *(Complex other) =>
      Complex(re * other.re - im * other.im, re * other.im + im * other.re);

  Complex operator /(Complex other) {
    final denom = other.re * other.re + other.im * other.im;
    return Complex(
      (re * other.re + im * other.im) / denom,
      (im * other.re - re * other.im) / denom,
    );
  }

  Complex scale(double s) => Complex(re * s, im * s);

  double get magnitude => sqrt(re * re + im * im);
  double get magnitudeSquared => re * re + im * im;
  double get argument => atan2(im, re);

  Complex get conjugate => Complex(re, -im);

  Complex get sqrt => () {
        final r = magnitude;
        final rSqrt = math_sqrt(r);
        final halfArg = argument / 2.0;
        return Complex(rSqrt * cos(halfArg), rSqrt * sin(halfArg));
      }();

  static double math_sqrt(double x) => sqrt(x);

  @override
  String toString() {
    if (im >= 0) return '${re.toStringAsFixed(4)} + ${im.toStringAsFixed(4)}i';
    return '${re.toStringAsFixed(4)} - ${(-im).toStringAsFixed(4)}i';
  }

  @override
  bool operator ==(Object other) =>
      other is Complex && re == other.re && im == other.im;

  @override
  int get hashCode => Object.hash(re, im);
}
