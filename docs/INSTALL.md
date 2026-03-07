# Installation Guide

## Prerequisites

| Dependency | Minimum Version | Notes |
|------------|----------------|-------|
| Flutter SDK | 3.10.0 | Includes Dart 3.0+ |
| Dart SDK | 3.0.0 | Bundled with Flutter |
| Git | 2.x | For cloning |

Download Flutter: https://docs.flutter.dev/get-started/install

---

## Quick Start

```bash
git clone https://github.com/vnemegabriel/fin_flutter.git
cd fin_flutter
flutter pub get
flutter run -d linux          # Desktop (Linux)
flutter run -d macos          # Desktop (macOS)
flutter run -d windows        # Desktop (Windows)
flutter run -d chrome         # Web
```

---

## Platform-Specific Setup

### Linux

Install GTK development libraries:

```bash
# Ubuntu / Debian
sudo apt update
sudo apt install libgtk-3-dev libblkid-dev liblzma-dev

# Fedora / RHEL
sudo dnf install gtk3-devel
```

Enable the Linux desktop target:

```bash
flutter config --enable-linux-desktop
flutter doctor   # verify Linux toolchain shows ✓
```

Build a release binary:

```bash
flutter build linux --release
# Output: build/linux/x64/release/bundle/fin_flutter
```

### macOS

Requires Xcode 14 or later (from the Mac App Store).

```bash
sudo xcode-select --install
xcodebuild -runFirstLaunch   # accept license

flutter config --enable-macos-desktop
flutter doctor               # verify Xcode toolchain ✓

flutter build macos --release
# Output: build/macos/Build/Products/Release/fin_flutter.app
```

### Windows

Requires Visual Studio 2022 (Community edition is free) with these workloads:
- Desktop development with C++
- Windows 10 SDK

```powershell
flutter config --enable-windows-desktop
flutter doctor               # verify Visual Studio toolchain ✓

flutter build windows --release
# Output: build\windows\x64\runner\Release\fin_flutter.exe
```

### Web

No extra system dependencies. Serve the output with any static file server.

```bash
flutter build web --release
# Output: build/web/

# Serve locally with Python
cd build/web && python3 -m http.server 8080
# Then open http://localhost:8080
```

> **Note:** The full analysis pipeline uses `dart:isolate` which is not available in web Workers
> under all browsers. If running on Web, the compute service falls back to running on the main
> thread. For production deployments, the Linux or macOS desktop build is recommended.

---

## Running Tests

No additional setup is required. The test suite uses only Flutter's built-in test framework.

```bash
# All tests
flutter test

# Specific subsystem
flutter test test/core/math/
flutter test test/core/materials/
flutter test test/core/fea/
flutter test test/core/cfd/
flutter test test/modules/optimization/
flutter test test/modules/openrocket/
flutter test test/integration/

# With coverage (requires lcov for HTML report)
flutter test --coverage
genhtml coverage/lcov.info -o coverage/html
open coverage/html/index.html
```

---

## IDE Setup

### VS Code

Install the [Flutter extension](https://marketplace.visualstudio.com/items?itemName=Dart-Code.flutter).
The workspace is automatically configured. Use `F5` to run.

### Android Studio / IntelliJ

Install the **Flutter** and **Dart** plugins via *Settings → Plugins*.
Open the project root directory; IntelliJ will detect `pubspec.yaml` automatically.

---

## Dependency Overview

Core Flutter packages used by this project (`pubspec.yaml`):

| Package | Use |
|---------|-----|
| `flutter_riverpod` | Reactive state management |
| `go_router` | Declarative navigation |
| `file_picker` | .ork file import dialog |
| `archive` | ZIP extraction for .ork parsing |
| `xml` | XML parsing for OpenRocket data |
| `fl_chart` | V-g diagram and optimization charts |
| `pdf` | Results export to PDF |
| `shared_preferences` | Project persistence |

---

## Troubleshooting

**`flutter doctor` shows missing tools**
Run `flutter doctor --verbose` and follow the links provided for each failing check.

**`flutter pub get` fails with network errors**
Set a pub proxy if behind a corporate firewall:
```bash
export PUB_HOSTED_URL=https://pub.flutter-io.cn  # China mirror
flutter pub get
```

**Linux build: `libgtk-3-dev not found`**
The GTK dev package was not installed. Re-run:
```bash
sudo apt install libgtk-3-dev
```

**Web: isolate compute not available**
The web platform has restricted isolate support. For computationally intensive analyses
(fine mesh, many modes), prefer the desktop build.

**LU decomposition `NaN` or `Inf` in results**
This usually indicates a near-singular stiffness matrix. Check:
1. Boundary conditions are applied (at least one clamped node)
2. Mesh does not contain collapsed elements (check span > 0, chords > 0)
3. Material D matrix is positive definite (non-zero ply thicknesses)
