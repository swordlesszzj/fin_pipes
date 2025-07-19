# Pype System Solver
A C++ project for solving pipe systems.

---
## Table of Contents
<!-- - [Introduction](#introduction) -->
<!-- - [Features](#features) -->
<!-- - [Installation](#installation) -->
- [Usage](#usage)
- [Project Structure](#project-structure)
- [License](#license)

---

## Usage
```bash
./build/pipe_solver Default.json
```

---

## Project Structure
```text
.
├── Default.txt
├── LICENSE
├── README.md
├── external
│   └── json.hpp
├── include
│   ├── Matrix.h
│   ├── MatrixUtils.h
│   ├── Pipe.h
│   ├── PipeSystem.h
│   └── SymmetricTridiagonalQR.h
├── src
│   ├── Matrix.cpp
│   ├── MatrixUtils.cpp
│   ├── PipeSystem.cpp
│   ├── SymmetricTridiagonalQR.cpp
│   └── main.cpp
└── build
    └── pipe_solver

```

---

## License
This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.