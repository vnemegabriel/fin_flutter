import 'package:flutter/material.dart';
import '../../core/models/fin_geometry.dart';
import '../theme/color_palette.dart';

/// CustomPainter widget displaying a 2D fin outline to scale.
class FinCanvas extends StatelessWidget {
  final FinGeometry geometry;
  final double width;
  final double height;

  const FinCanvas({
    super.key,
    required this.geometry,
    this.width = 300,
    this.height = 200,
  });

  @override
  Widget build(BuildContext context) {
    return SizedBox(
      width: width,
      height: height,
      child: CustomPaint(
        painter: _FinPainter(geometry),
        size: Size(width, height),
      ),
    );
  }
}

class _FinPainter extends CustomPainter {
  final FinGeometry geo;

  _FinPainter(this.geo);

  @override
  void paint(Canvas canvas, Size size) {
    const padding = 20.0;
    final drawWidth = size.width - 2 * padding;
    final drawHeight = size.height - 2 * padding;

    // Scale to fit
    final maxX = geo.trailingEdgeX(geo.span);
    final maxY = geo.span;
    final scaleX = drawWidth / maxX;
    final scaleY = drawHeight / maxY;
    final scale = scaleX < scaleY ? scaleX : scaleY;

    // Offset to center
    final offsetX = padding + (drawWidth - maxX * scale) / 2;
    final offsetY = padding + drawHeight; // Y increases downward, span goes up

    Offset toCanvas(double x, double y) =>
        Offset(offsetX + x * scale, offsetY - y * scale);

    // Fin outline points
    final rootLE = toCanvas(0, 0);
    final rootTE = toCanvas(geo.rootChord, 0);
    final tipTE = toCanvas(geo.trailingEdgeX(geo.span), geo.span);
    final tipLE = toCanvas(geo.leadingEdgeX(geo.span), geo.span);

    final path = Path()
      ..moveTo(rootLE.dx, rootLE.dy)
      ..lineTo(rootTE.dx, rootTE.dy)
      ..lineTo(tipTE.dx, tipTE.dy)
      ..lineTo(tipLE.dx, tipLE.dy)
      ..close();

    // Fill
    canvas.drawPath(path, Paint()..color = ColorPalette.primary.withOpacity(0.15));

    // Outline
    canvas.drawPath(
        path,
        Paint()
          ..color = ColorPalette.primary
          ..style = PaintingStyle.stroke
          ..strokeWidth = 2.0);

    // Root label
    final textPainter = TextPainter(
      text: TextSpan(
        text: 'root',
        style: TextStyle(color: Colors.grey[600], fontSize: 10),
      ),
      textDirection: TextDirection.ltr,
    )..layout();
    textPainter.paint(canvas, Offset(rootLE.dx + 4, rootLE.dy - 14));

    // Dimension annotations
    _drawDimension(
      canvas,
      Offset(offsetX - 12, offsetY),
      Offset(offsetX - 12, offsetY - geo.span * scale),
      '${(geo.span * 100).toStringAsFixed(1)}cm',
      true,
    );
  }

  void _drawDimension(
      Canvas canvas, Offset start, Offset end, String label, bool vertical) {
    final paint = Paint()
      ..color = Colors.grey
      ..strokeWidth = 1.0;
    canvas.drawLine(start, end, paint);

    final textPainter = TextPainter(
      text: TextSpan(
        text: label,
        style: const TextStyle(color: Colors.grey, fontSize: 9),
      ),
      textDirection: TextDirection.ltr,
    )..layout();

    final midX = (start.dx + end.dx) / 2;
    final midY = (start.dy + end.dy) / 2;
    textPainter.paint(canvas, Offset(midX - textPainter.width / 2, midY));
  }

  @override
  bool shouldRepaint(_FinPainter oldDelegate) => oldDelegate.geo != geo;
}
