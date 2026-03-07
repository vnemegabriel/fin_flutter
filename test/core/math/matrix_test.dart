import 'package:flutter_test/flutter_test.dart';
import 'package:fin_flutter/core/math/matrix.dart';

void main() {
  group('Matrix', () {
    test('identity matrix', () {
      final I = Matrix.identity(3);
      for (var i = 0; i < 3; i++) {
        for (var j = 0; j < 3; j++) {
          expect(I.get(i, j), i == j ? 1.0 : 0.0);
        }
      }
    });

    test('matrix multiplication', () {
      final A = Matrix.fromList([
        [1.0, 2.0],
        [3.0, 4.0],
      ]);
      final B = Matrix.fromList([
        [5.0, 6.0],
        [7.0, 8.0],
      ]);
      final C = A * B;
      expect(C.get(0, 0), closeTo(19.0, 1e-10));
      expect(C.get(0, 1), closeTo(22.0, 1e-10));
      expect(C.get(1, 0), closeTo(43.0, 1e-10));
      expect(C.get(1, 1), closeTo(50.0, 1e-10));
    });

    test('transpose', () {
      final A = Matrix.fromList([
        [1.0, 2.0, 3.0],
        [4.0, 5.0, 6.0],
      ]);
      final At = A.transpose();
      expect(At.rows, 3);
      expect(At.cols, 2);
      expect(At.get(0, 0), 1.0);
      expect(At.get(2, 1), 6.0);
    });

    test('LU solve Ax=b', () {
      final A = Matrix.fromList([
        [2.0, 1.0, 1.0],
        [4.0, 3.0, 3.0],
        [8.0, 7.0, 9.0],
      ]);
      final b = [1.0, 2.0, 4.0];
      final x = A.solve(b);

      // Verify A*x ≈ b
      final Ax = A.multiplyVector(x);
      for (var i = 0; i < 3; i++) {
        expect(Ax[i], closeTo(b[i], 1e-8));
      }
    });

    test('inverse: A * A_inv ≈ I', () {
      final A = Matrix.fromList([
        [2.0, 1.0],
        [1.0, 3.0],
      ]);
      final Ainv = A.inverse();
      final product = A * Ainv;
      for (var i = 0; i < 2; i++) {
        for (var j = 0; j < 2; j++) {
          expect(product.get(i, j), closeTo(i == j ? 1.0 : 0.0, 1e-10));
        }
      }
    });

    test('determinant of 3x3 matrix', () {
      final A = Matrix.fromList([
        [1.0, 2.0, 3.0],
        [0.0, 1.0, 4.0],
        [5.0, 6.0, 0.0],
      ]);
      // det = 1*(1*0 - 4*6) - 2*(0*0 - 4*5) + 3*(0*6 - 1*5)
      //     = 1*(-24) - 2*(-20) + 3*(-5) = -24 + 40 - 15 = 1
      expect(A.determinant(), closeTo(1.0, 1e-8));
    });

    test('identity matrix scale', () {
      final I = Matrix.identity(3);
      final S = I.scale(5.0);
      for (var i = 0; i < 3; i++) {
        expect(S.get(i, i), 5.0);
      }
    });
  });
}
