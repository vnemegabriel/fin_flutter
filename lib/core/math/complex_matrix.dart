import 'dart:math';
import 'complex_number.dart';
import 'matrix.dart';

/// Dense complex matrix for flutter eigenvalue problems.
class ComplexMatrix {
  final int rows;
  final int cols;
  // Stored as flat lists for real and imaginary parts.
  final List<double> _re;
  final List<double> _im;

  ComplexMatrix(this.rows, this.cols)
      : _re = List<double>.filled(rows * cols, 0.0),
        _im = List<double>.filled(rows * cols, 0.0);

  factory ComplexMatrix.fromReal(Matrix m) {
    final cm = ComplexMatrix(m.rows, m.cols);
    for (var i = 0; i < m.rows; i++) {
      for (var j = 0; j < m.cols; j++) {
        cm.set(i, j, Complex(m.get(i, j), 0.0));
      }
    }
    return cm;
  }

  factory ComplexMatrix.identity(int n) {
    final m = ComplexMatrix(n, n);
    for (var i = 0; i < n; i++) {
      m.set(i, i, Complex.one);
    }
    return m;
  }

  Complex get(int i, int j) {
    final idx = i * cols + j;
    return Complex(_re[idx], _im[idx]);
  }

  void set(int i, int j, Complex v) {
    final idx = i * cols + j;
    _re[idx] = v.re;
    _im[idx] = v.im;
  }

  ComplexMatrix operator +(ComplexMatrix other) {
    assert(rows == other.rows && cols == other.cols);
    final result = ComplexMatrix(rows, cols);
    for (var k = 0; k < _re.length; k++) {
      result._re[k] = _re[k] + other._re[k];
      result._im[k] = _im[k] + other._im[k];
    }
    return result;
  }

  ComplexMatrix operator -(ComplexMatrix other) {
    assert(rows == other.rows && cols == other.cols);
    final result = ComplexMatrix(rows, cols);
    for (var k = 0; k < _re.length; k++) {
      result._re[k] = _re[k] - other._re[k];
      result._im[k] = _im[k] - other._im[k];
    }
    return result;
  }

  ComplexMatrix operator *(ComplexMatrix other) {
    assert(cols == other.rows);
    final result = ComplexMatrix(rows, other.cols);
    for (var i = 0; i < rows; i++) {
      for (var k = 0; k < cols; k++) {
        final aRe = _re[i * cols + k];
        final aIm = _im[i * cols + k];
        if (aRe == 0.0 && aIm == 0.0) continue;
        for (var j = 0; j < other.cols; j++) {
          final bRe = other._re[k * other.cols + j];
          final bIm = other._im[k * other.cols + j];
          result._re[i * other.cols + j] += aRe * bRe - aIm * bIm;
          result._im[i * other.cols + j] += aRe * bIm + aIm * bRe;
        }
      }
    }
    return result;
  }

  ComplexMatrix scale(Complex s) {
    final result = ComplexMatrix(rows, cols);
    for (var k = 0; k < _re.length; k++) {
      final re = _re[k];
      final im = _im[k];
      result._re[k] = re * s.re - im * s.im;
      result._im[k] = re * s.im + im * s.re;
    }
    return result;
  }

  ComplexMatrix scaleReal(double s) {
    final result = ComplexMatrix(rows, cols);
    for (var k = 0; k < _re.length; k++) {
      result._re[k] = _re[k] * s;
      result._im[k] = _im[k] * s;
    }
    return result;
  }

  ComplexMatrix transpose() {
    final result = ComplexMatrix(cols, rows);
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < cols; j++) {
        result.set(j, i, get(i, j));
      }
    }
    return result;
  }

  ComplexMatrix conjugateTranspose() {
    final result = ComplexMatrix(cols, rows);
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < cols; j++) {
        result.set(j, i, get(i, j).conjugate);
      }
    }
    return result;
  }

  /// Solve AX = b (single RHS) using Gaussian elimination with partial pivoting.
  List<Complex> solve(List<Complex> b) {
    assert(rows == cols && cols == b.length);
    final n = rows;
    final aug = ComplexMatrix(n, n);
    for (var i = 0; i < n; i++) {
      for (var j = 0; j < n; j++) {
        aug.set(i, j, get(i, j));
      }
    }
    final x = List<Complex>.from(b);

    for (var k = 0; k < n; k++) {
      // Partial pivot
      var maxMag = aug.get(k, k).magnitude;
      var maxRow = k;
      for (var i = k + 1; i < n; i++) {
        final mag = aug.get(i, k).magnitude;
        if (mag > maxMag) {
          maxMag = mag;
          maxRow = i;
        }
      }
      if (maxRow != k) {
        for (var j = 0; j < n; j++) {
          final tmp = aug.get(k, j);
          aug.set(k, j, aug.get(maxRow, j));
          aug.set(maxRow, j, tmp);
        }
        final tx = x[k];
        x[k] = x[maxRow];
        x[maxRow] = tx;
      }
      final ukk = aug.get(k, k);
      if (ukk.magnitude < 1e-14) continue;
      for (var i = k + 1; i < n; i++) {
        final factor = aug.get(i, k) / ukk;
        x[i] = x[i] - factor * x[k];
        for (var j = k + 1; j < n; j++) {
          aug.set(i, j, aug.get(i, j) - factor * aug.get(k, j));
        }
        aug.set(i, k, Complex.zero);
      }
    }

    for (var i = n - 1; i >= 0; i--) {
      for (var j = i + 1; j < n; j++) {
        x[i] = x[i] - aug.get(i, j) * x[j];
      }
      x[i] = x[i] / aug.get(i, i);
    }
    return x;
  }

  /// Multiply by column vector.
  List<Complex> multiplyVector(List<Complex> v) {
    assert(cols == v.length);
    final result = List<Complex>.filled(rows, Complex.zero);
    for (var i = 0; i < rows; i++) {
      Complex sum = Complex.zero;
      for (var j = 0; j < cols; j++) {
        sum = sum + get(i, j) * v[j];
      }
      result[i] = sum;
    }
    return result;
  }

  @override
  String toString() {
    final buf = StringBuffer('ComplexMatrix($rows x $cols)\n');
    for (var i = 0; i < rows; i++) {
      buf.write('  [');
      for (var j = 0; j < cols; j++) {
        buf.write(get(i, j).toString());
        if (j < cols - 1) buf.write(', ');
      }
      buf.writeln(']');
    }
    return buf.toString();
  }
}
