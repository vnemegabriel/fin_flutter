# Installation Guide

fin_flutter is a C++/Python computational platform. This guide covers building the C++ core and setting up the Python orchestration layer.

---

## Prerequisites

### C++ Core

| Dependency | Minimum | Notes |
|-----------|---------|-------|
| **CMake** | 3.16 | Build system |
| **C++ Compiler** | C++17 | g++, clang, MSVC |
| **Eigen3** | 3.3 | Linear algebra library |
| **Git** | 2.x | Version control |

### Python Layer (Optional)

| Dependency | Minimum | Notes |
|-----------|---------|-------|
| **Python** | 3.10 | For orchestration |
| **pip** | 21.0 | Package manager |

---

## Quick Start

### 1. Clone Repository

```bash
git clone https://github.com/vnemegabriel/fin_flutter.git
cd fin_flutter
```

### 2. Build C++ Core

```bash
cd core/cpp
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)

# Run benchmark tests
./build/tests/run_tests
```

**Expected output:**
```
=== Results: 23 passed, 1 failed ===
Case 4 — VLM rectangular flat plate AR=5, alpha=5deg
  [FAIL] CL in [0.36, 0.42] CL=0.007746
```

(VLM magnitude undersizing under investigation; sign is correct.)

### 3. Setup Python Layer (Optional)

```bash
cd core/python
pip install -e ".[dev]"
pytest
```

---

## Platform-Specific Setup

### Linux (Ubuntu / Debian)

#### Install Dependencies

```bash
# Update package lists
sudo apt update

# C++ build tools and Eigen3
sudo apt install -y \
  build-essential \
  cmake \
  libeigen3-dev \
  git

# Python (if using orchestration layer)
sudo apt install -y \
  python3.10 \
  python3-pip \
  python3-venv
```

#### Verify Installation

```bash
cmake --version      # CMake >= 3.16
g++ --version        # GCC with C++17 support
python3 --version    # Python >= 3.10
```

#### Build & Test

```bash
cd fin_flutter/core/cpp
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
./build/tests/run_tests
```

---

### macOS

#### Install Dependencies

```bash
# Homebrew (if not installed)
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"

# C++ build tools and Eigen3
brew install cmake eigen

# Python (if using orchestration layer)
brew install python@3.10
```

#### Verify Installation

```bash
cmake --version
clang++ --version
python3 --version
```

#### Build & Test

```bash
cd fin_flutter/core/cpp
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
./build/tests/run_tests
```

---

### Windows (MSVC)

#### Install Dependencies

1. **Visual Studio Community 2022**
   - Download: https://visualstudio.microsoft.com/
   - Ensure "Desktop development with C++" workload is selected

2. **CMake**
   - Download: https://cmake.org/download/
   - Add to PATH during installation

3. **Eigen3**
   - Download: https://eigen.tuxfamily.org/index.php
   - Extract to a known location (e.g., `C:\Libraries\eigen-3.4.0`)

4. **Python** (optional)
   - Download: https://www.python.org/downloads/
   - Add to PATH during installation

#### Build & Test

```bash
cd fin_flutter\core\cpp
cmake -B build -DCMAKE_BUILD_TYPE=Release -DEIGEN3_INCLUDE_DIR=C:\Libraries\eigen-3.4.0
cmake --build build --config Release
build\tests\Release\run_tests.exe
```

(Adjust `EIGEN3_INCLUDE_DIR` path as needed.)

---

## Troubleshooting

### CMake Error: "Could not find a package configuration file provided by 'Eigen3'"

**Solution:** Install Eigen3 via package manager or specify path:

```bash
# Ubuntu/Debian
sudo apt install libeigen3-dev

# Or manually specify
cmake -DEIGEN3_INCLUDE_DIR=/path/to/eigen3 ..
```

### CMake Error: "The C compiler 'cc' is not able to compile a simple test program"

**Solution:** Install C++ compiler:

```bash
# Ubuntu/Debian
sudo apt install build-essential

# macOS
xcode-select --install

# Windows: Install Visual Studio (see above)
```

### Build Fails with "fatal error: Eigen/Dense: No such file or directory"

**Solution:** Eigen3 not found. Verify installation:

```bash
# Linux/macOS
find /usr -name "Eigen" 2>/dev/null
find /opt -name "Eigen" 2>/dev/null

# macOS with Homebrew
brew list eigen

# If installed, try:
cmake -DEIGEN3_INCLUDE_DIR=$(brew --prefix eigen)/include/eigen3 ..
```

### Tests Pass but Report "CL not in range [0.36, 0.42]"

This is Case 4 (VLM flat plate test). The sign is now correct (positive CL), but magnitude is undersized by ~50×. See [docs/test_cases.md](test_cases.md) for details. This is a known issue under investigation.

---

## Building with Different Configurations

### Release (Optimized)

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j$(nproc)
```

### Debug (With Symbols & Assertions)

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j$(nproc)
./build/tests/run_tests
```

### Verbose Build Output

```bash
cmake --build build --verbose
```

---

## Python Setup (Optional)

### Create Virtual Environment

```bash
cd core/python
python3 -m venv venv
source venv/bin/activate      # Linux/macOS
# or
venv\Scripts\activate          # Windows
```

### Install Package in Development Mode

```bash
pip install -e ".[dev]"
```

### Run Tests

```bash
pytest
```

### Run CLI Help

```bash
python -m fin_flutter.cli --help
```

(CLI not yet implemented; placeholder for future release.)

---

## Verifying Installation

### C++ Core

```bash
cd core/cpp
./build/tests/run_tests
```

Should output:
```
=== Results: 23 passed, 1 failed ===
```

### Python Environment

```bash
python3 -c "import numpy; import matplotlib; print('OK')"
```

---

## Next Steps

- **Read [docs/theory.md](theory.md)** for mathematical background
- **Explore [core/cpp/tests/test_main.cpp](../core/cpp/tests/test_main.cpp)** for usage examples
- **Check [docs/architecture.md](architecture.md)** for module structure
- **See [CLAUDE.md](../CLAUDE.md)** for development guidelines

---

## Support

For issues:
1. Check [Troubleshooting](#troubleshooting) section above
2. Verify all prerequisites are installed
3. Review [docs/test_cases.md](test_cases.md) for known limitations
4. Open an issue on GitHub with CMake output (run with `-DCMAKE_VERBOSE_MAKEFILE=ON`)
