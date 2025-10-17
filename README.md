# Subtree Mode with Applications

Official implementation of "Subtree Mode with Applications".

## Overview

This project addresses the **Subtree Mode** problem: finding the most frequent value (the *mode*) within the leaves of any given subtree in a leaf-colored tree. This is a fundamental challenge with direct applications in string-processing, such as document retrieval, pattern mining, and sequence-to-database searches.

Our core contribution is a highly efficient algorithm that calculates the mode for **every subtree** in an input tree of $N$ nodes in optimal **$O(N)$ time**. This repository provides the implementation of this algorithm (`SCM`) along with several baselines for comparison.

Our experiments on massive real-world datasets (with trees up to 7.3 billion nodes) show that our approach is not only an order of magnitude faster than baselines but also significantly more memory-efficient.

---

## Requirements

* **OS**: A GNU/Linux system (e.g., Ubuntu, CentOS).
* **Compiler**: A modern C++17 ready compiler (e.g., GCC 7+ or Clang 6+).
* **Dependencies**: This project requires `sdsl-lite` and `libdivsufsort`. Please follow these steps to install them:

    1.  **Install `sdsl-lite`**
        ```bash
        git clone [https://github.com/simongog/sdsl-lite.git](https://github.com/simongog/sdsl-lite.git)
        cd sdsl-lite
        ./install.sh
        cd ..
        ```

    2.  **Install `libdivsufsort`**
        ```bash
        git clone [https://github.com/y-256/libdivsufsort.git](https://github.com/y-256/libdivsufsort.git)
        cd libdivsufsort
        mkdir build && cd build
        cmake ..
        make
        sudo make install
        cd ../..
        ```

---

## Compilation

All programs should be compiled with the C++17 standard and optimization flags for best performance.

```bash
# Baseline 1
g++ -std=c++17 -O3 -o baseline1 Baseline1.cpp

# Baseline 2
g++ -std=c++17 -O3 -o baseline2 Baseline2.cpp

# Baseline 3
g++ -std=c++17 -O3 -o baseline3 Baseline3.cpp -ldivsufsort64 -lsdsl

# SCM (Our Algorithm)
g++ -std=c++17 -O3 -D_USE_64 -I. SCM.cpp MultiColorESA.cpp -ldivsufsort64 -lsdsl -o scm
```

---

## Datasets

The datasets used in our experiments can be found in the following Google Drive folder:

[https://drive.google.com/drive/folders/1Hmp5RKhoZTOX-McSD_2VqPi0K4PcazgH?usp=sharing](https://drive.google.com/drive/folders/1Hmp5RKhoZTOX-McSD_2VqPi0K4PcazgH?usp=sharing)

---

## Execution

The compiled programs share a common usage pattern.

#### General Usage
```
./<program_name> <input_file> <output_file>
```

#### Example
To run the `scm` program with `input.txt` and write results to `output.txt`:
```
./scm input.txt output.txt
```