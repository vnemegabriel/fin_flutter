import 'dart:math';
import 'matrix.dart';

/// Sparse matrix in Compressed Sparse Row (CSR) format.
/// Used for large FEA global stiffness/mass matrices.
class SparseMatrix {
  final int rows;
  final int cols;
  final List<double> values;    // Non-zero values
  final List<int> colIndices;   // Column index of each value
  final List<int> rowPointers;  // rowPointers[i] = start of row i in values

  SparseMatrix({
    required this.rows,
    required this.cols,
    required this.values,
    required this.colIndices,
    required this.rowPointers,
  });

  factory SparseMatrix.fromDense(Matrix m) {
    final vals = <double>[];
    final cols = <int>[];
    final rowPtrs = <int>[0];
    for (var i = 0; i < m.rows; i++) {
      for (var j = 0; j < m.cols; j++) {
        final v = m.get(i, j);
        if (v.abs() > 1e-15) {
          vals.add(v);
          cols.add(j);
        }
      }
      rowPtrs.add(vals.length);
    }
    return SparseMatrix(
      rows: m.rows,
      cols: m.cols,
      values: vals,
      colIndices: cols,
      rowPointers: rowPtrs,
    );
  }

  /// Convert to dense Matrix.
  Matrix toDense() {
    final m = Matrix(rows, cols);
    for (var i = 0; i < rows; i++) {
      for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++) {
        m.set(i, colIndices[k], values[k]);
      }
    }
    return m;
  }

  /// Multiply by vector.
  List<double> multiplyVector(List<double> v) {
    assert(cols == v.length);
    final result = List<double>.filled(rows, 0.0);
    for (var i = 0; i < rows; i++) {
      double sum = 0.0;
      for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++) {
        sum += values[k] * v[colIndices[k]];
      }
      result[i] = sum;
    }
    return result;
  }

  /// Get element (slow - O(nnz per row)).
  double get(int i, int j) {
    for (var k = rowPointers[i]; k < rowPointers[i + 1]; k++) {
      if (colIndices[k] == j) return values[k];
    }
    return 0.0;
  }
}

/// Accumulator for building sparse matrices element-by-element.
class SparseMatrixBuilder {
  final int rows;
  final int cols;
  final Map<int, Map<int, double>> _entries = {};

  SparseMatrixBuilder({required this.rows, required this.cols});

  void add(int i, int j, double v) {
    _entries.putIfAbsent(i, () => {})[j] =
        (_entries[i]?[j] ?? 0.0) + v;
  }

  void set(int i, int j, double v) {
    _entries.putIfAbsent(i, () => {})[j] = v;
  }

  SparseMatrix build() {
    final vals = <double>[];
    final colIdx = <int>[];
    final rowPtrs = <int>[0];

    for (var i = 0; i < rows; i++) {
      final rowMap = _entries[i];
      if (rowMap != null) {
        final sortedCols = rowMap.keys.toList()..sort();
        for (final j in sortedCols) {
          vals.add(rowMap[j]!);
          colIdx.add(j);
        }
      }
      rowPtrs.add(vals.length);
    }
    return SparseMatrix(
      rows: rows,
      cols: cols,
      values: vals,
      colIndices: colIdx,
      rowPointers: rowPtrs,
    );
  }

  Matrix toDense() => build().toDense();
}
