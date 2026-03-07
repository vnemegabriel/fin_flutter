import 'dart:math';
import 'dart:typed_data';

/// Dense real matrix backed by a flat Float64List (row-major).
class Matrix {
  final int rows;
  final int cols;
  final Float64List _data;

  Matrix(this.rows, this.cols) : _data = Float64List(rows * cols);

  Matrix._internal(this.rows, this.cols, this._data);

  factory Matrix.identity(int n) {
    final m = Matrix(n, n);
    for (var i = 0; i < n; i++) {
      m.set(i, i, 1.0);
    }
    return m;
  }

  factory Matrix.fromList(List<List<double>> data) {
    final rows = data.length;
    final cols = data.isEmpty ? 0 : data[0].length;
    final m = Matrix(rows, cols);
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < cols; j++) {
        m.set(i, j, data[i][j]);
      }
    }
    return m;
  }

  factory Matrix.zeros(int rows, int cols) => Matrix(rows, cols);

  factory Matrix.diagonal(List<double> diag) {
    final n = diag.length;
    final m = Matrix(n, n);
    for (var i = 0; i < n; i++) {
      m.set(i, i, diag[i]);
    }
    return m;
  }

  double get(int i, int j) => _data[i * cols + j];
  void set(int i, int j, double v) => _data[i * cols + j] = v;

  double operator [](int index) => _data[index];

  Matrix operator +(Matrix other) {
    assert(rows == other.rows && cols == other.cols);
    final result = Matrix(rows, cols);
    for (var k = 0; k < _data.length; k++) {
      result._data[k] = _data[k] + other._data[k];
    }
    return result;
  }

  Matrix operator -(Matrix other) {
    assert(rows == other.rows && cols == other.cols);
    final result = Matrix(rows, cols);
    for (var k = 0; k < _data.length; k++) {
      result._data[k] = _data[k] - other._data[k];
    }
    return result;
  }

  Matrix operator *(Matrix other) {
    assert(cols == other.rows, 'Dimension mismatch: ($rows x $cols) * (${other.rows} x ${other.cols})');
    final result = Matrix(rows, other.cols);
    for (var i = 0; i < rows; i++) {
      for (var k = 0; k < cols; k++) {
        final aik = get(i, k);
        if (aik == 0.0) continue;
        for (var j = 0; j < other.cols; j++) {
          result._data[i * other.cols + j] += aik * other.get(k, j);
        }
      }
    }
    return result;
  }

  Matrix scale(double s) {
    final result = Matrix(rows, cols);
    for (var k = 0; k < _data.length; k++) {
      result._data[k] = _data[k] * s;
    }
    return result;
  }

  Matrix transpose() {
    final result = Matrix(cols, rows);
    for (var i = 0; i < rows; i++) {
      for (var j = 0; j < cols; j++) {
        result.set(j, i, get(i, j));
      }
    }
    return result;
  }

  /// Extract submatrix from [r0,c0] to [r1,c1) (exclusive end).
  Matrix subMatrix(int r0, int c0, int r1, int c1) {
    final nr = r1 - r0;
    final nc = c1 - c0;
    final result = Matrix(nr, nc);
    for (var i = 0; i < nr; i++) {
      for (var j = 0; j < nc; j++) {
        result.set(i, j, get(r0 + i, c0 + j));
      }
    }
    return result;
  }

  /// Multiply matrix by column vector (List<double>), return List<double>.
  List<double> multiplyVector(List<double> v) {
    assert(cols == v.length);
    final result = List<double>.filled(rows, 0.0);
    for (var i = 0; i < rows; i++) {
      double sum = 0.0;
      for (var j = 0; j < cols; j++) {
        sum += get(i, j) * v[j];
      }
      result[i] = sum;
    }
    return result;
  }

  /// LU decomposition with partial pivoting (Doolittle).
  /// Returns (L, U, pivot) where pivot[i] is the original row that went to row i.
  ({Matrix L, Matrix U, List<int> pivot}) luDecompose() {
    assert(rows == cols, 'LU decomposition requires square matrix');
    final n = rows;
    final lu = Matrix(n, n);
    // Copy data
    for (var k = 0; k < _data.length; k++) {
      lu._data[k] = _data[k];
    }
    final pivot = List<int>.generate(n, (i) => i);

    for (var k = 0; k < n; k++) {
      // Find pivot
      var maxVal = lu.get(k, k).abs();
      var maxRow = k;
      for (var i = k + 1; i < n; i++) {
        final v = lu.get(i, k).abs();
        if (v > maxVal) {
          maxVal = v;
          maxRow = i;
        }
      }
      // Swap rows
      if (maxRow != k) {
        for (var j = 0; j < n; j++) {
          final tmp = lu.get(k, j);
          lu.set(k, j, lu.get(maxRow, j));
          lu.set(maxRow, j, tmp);
        }
        final tmp = pivot[k];
        pivot[k] = pivot[maxRow];
        pivot[maxRow] = tmp;
      }
      // Eliminate
      final ukk = lu.get(k, k);
      if (ukk.abs() < 1e-14) continue;
      for (var i = k + 1; i < n; i++) {
        final factor = lu.get(i, k) / ukk;
        lu.set(i, k, factor);
        for (var j = k + 1; j < n; j++) {
          lu.set(i, j, lu.get(i, j) - factor * lu.get(k, j));
        }
      }
    }

    // Extract L and U
    final L = Matrix.identity(n);
    final U = Matrix(n, n);
    for (var i = 0; i < n; i++) {
      for (var j = 0; j < n; j++) {
        if (i > j) {
          L.set(i, j, lu.get(i, j));
        } else {
          U.set(i, j, lu.get(i, j));
        }
      }
    }
    return (L: L, U: U, pivot: pivot);
  }

  /// Solve Ax = b using LU decomposition.
  List<double> solve(List<double> b) {
    assert(rows == cols && cols == b.length);
    final n = rows;
    final lu = Matrix(n, n);
    for (var k = 0; k < _data.length; k++) {
      lu._data[k] = _data[k];
    }
    final x = List<double>.from(b);

    // LU factorize with partial pivoting
    final pivot = List<int>.generate(n, (i) => i);
    for (var k = 0; k < n; k++) {
      var maxVal = lu.get(k, k).abs();
      var maxRow = k;
      for (var i = k + 1; i < n; i++) {
        final v = lu.get(i, k).abs();
        if (v > maxVal) {
          maxVal = v;
          maxRow = i;
        }
      }
      if (maxRow != k) {
        for (var j = 0; j < n; j++) {
          final tmp = lu.get(k, j);
          lu.set(k, j, lu.get(maxRow, j));
          lu.set(maxRow, j, tmp);
        }
        final tmp = pivot[k];
        pivot[k] = pivot[maxRow];
        pivot[maxRow] = tmp;
        final tx = x[k];
        x[k] = x[maxRow];
        x[maxRow] = tx;
      }
      final ukk = lu.get(k, k);
      if (ukk.abs() < 1e-14) continue;
      for (var i = k + 1; i < n; i++) {
        final factor = lu.get(i, k) / ukk;
        lu.set(i, k, factor);
        x[i] -= factor * x[k];
        for (var j = k + 1; j < n; j++) {
          lu.set(i, j, lu.get(i, j) - factor * lu.get(k, j));
        }
      }
    }

    // Back substitution
    for (var i = n - 1; i >= 0; i--) {
      for (var j = i + 1; j < n; j++) {
        x[i] -= lu.get(i, j) * x[j];
      }
      x[i] /= lu.get(i, i);
    }
    return x;
  }

  /// Solve AX = B where B is a matrix (multiple right-hand sides).
  Matrix solveMatrix(Matrix B) {
    assert(rows == cols && cols == B.rows);
    final result = Matrix(rows, B.cols);
    for (var j = 0; j < B.cols; j++) {
      final col = List<double>.generate(rows, (i) => B.get(i, j));
      final sol = solve(col);
      for (var i = 0; i < rows; i++) {
        result.set(i, j, sol[i]);
      }
    }
    return result;
  }

  /// Compute inverse using LU decomposition.
  Matrix inverse() {
    assert(rows == cols);
    return solveMatrix(Matrix.identity(rows));
  }

  /// Compute Frobenius norm.
  double norm() {
    double sum = 0.0;
    for (final v in _data) {
      sum += v * v;
    }
    return sqrt(sum);
  }

  /// Compute determinant via LU decomposition.
  double determinant() {
    assert(rows == cols);
    final n = rows;
    final lu = Matrix(n, n);
    for (var k = 0; k < _data.length; k++) {
      lu._data[k] = _data[k];
    }
    var sign = 1.0;

    for (var k = 0; k < n; k++) {
      var maxVal = lu.get(k, k).abs();
      var maxRow = k;
      for (var i = k + 1; i < n; i++) {
        final v = lu.get(i, k).abs();
        if (v > maxVal) {
          maxVal = v;
          maxRow = i;
        }
      }
      if (maxRow != k) {
        for (var j = 0; j < n; j++) {
          final tmp = lu.get(k, j);
          lu.set(k, j, lu.get(maxRow, j));
          lu.set(maxRow, j, tmp);
        }
        sign = -sign;
      }
      final ukk = lu.get(k, k);
      if (ukk.abs() < 1e-14) return 0.0;
      for (var i = k + 1; i < n; i++) {
        final factor = lu.get(i, k) / ukk;
        lu.set(i, k, factor);
        for (var j = k + 1; j < n; j++) {
          lu.set(i, j, lu.get(i, j) - factor * lu.get(k, j));
        }
      }
    }

    var det = sign;
    for (var i = 0; i < n; i++) {
      det *= lu.get(i, i);
    }
    return det;
  }

  /// Convert to List<List<double>>.
  List<List<double>> toListList() {
    return List.generate(rows, (i) => List.generate(cols, (j) => get(i, j)));
  }

  /// Column as list.
  List<double> column(int j) => List.generate(rows, (i) => get(i, j));

  /// Row as list.
  List<double> row(int i) => List.generate(cols, (j) => get(i, j));

  @override
  String toString() {
    final buf = StringBuffer('Matrix($rows x $cols):\n');
    for (var i = 0; i < rows; i++) {
      buf.write('  [');
      for (var j = 0; j < cols; j++) {
        buf.write(get(i, j).toStringAsFixed(6));
        if (j < cols - 1) buf.write(', ');
      }
      buf.writeln(']');
    }
    return buf.toString();
  }
}

/// Helper: dot product of two vectors.
double dotProduct(List<double> a, List<double> b) {
  assert(a.length == b.length);
  double sum = 0.0;
  for (var i = 0; i < a.length; i++) {
    sum += a[i] * b[i];
  }
  return sum;
}

/// Helper: vector norm.
double vectorNorm(List<double> v) {
  double sum = 0.0;
  for (final x in v) {
    sum += x * x;
  }
  return sqrt(sum);
}
