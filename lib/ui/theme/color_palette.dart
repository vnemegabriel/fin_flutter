import 'package:flutter/material.dart';

/// Application color palette.
class ColorPalette {
  // Primary colors - aerospace engineering theme
  static const primary = Color(0xFF1565C0);       // Deep blue
  static const primaryVariant = Color(0xFF003C8F); // Dark blue
  static const secondary = Color(0xFFFF8F00);      // Amber accent
  static const secondaryVariant = Color(0xFFC56000);

  // Background
  static const background = Color(0xFFF5F7FA);
  static const surface = Colors.white;
  static const surfaceVariant = Color(0xFFEDF2F7);

  // Status colors
  static const safe = Color(0xFF2E7D32);
  static const warning = Color(0xFFF57F17);
  static const marginal = Color(0xFFE65100);
  static const critical = Color(0xFFC62828);

  // Chart colors (V-g diagram modes)
  static const List<Color> modeColors = [
    Color(0xFF1565C0), // Mode 1: Blue
    Color(0xFFAD1457), // Mode 2: Pink
    Color(0xFF2E7D32), // Mode 3: Green
    Color(0xFFE65100), // Mode 4: Orange
    Color(0xFF6A1B9A), // Mode 5: Purple
    Color(0xFF00838F), // Mode 6: Teal
  ];

  // Pressure map gradient
  static const pressureLow = Color(0xFF0D47A1);
  static const pressureMid = Color(0xFFFFFFFF);
  static const pressureHigh = Color(0xFFB71C1C);
}
