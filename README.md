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
* **Dependencies**: 
  - **For SCM_k-most**: Only `libdivsufsort` is required
  - **For baselines**: Both `sdsl-lite` and `libdivsufsort` are required
  
  Please follow these steps to install them:

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

# SCM_k-most (Our Algorithm - computes top-k most frequent colors)
# Note: When k=1, this is equivalent to the original SCM algorithm
g++ -std=c++17 -O3 -DUSE_DIVSUFSORT -o scm_k_most SCM_k-most.cpp -ldivsufsort
```

---

## Datasets

The datasets used in our experiments can be found in the following Google Drive folder:

[https://drive.google.com/drive/folders/1Hmp5RKhoZTOX-McSD_2VqPi0K4PcazgH?usp=sharing](https://drive.google.com/drive/folders/1Hmp5RKhoZTOX-McSD_2VqPi0K4PcazgH?usp=sharing)

---

## Execution

#### Baselines
The baseline programs share a common usage pattern:
```
./<program_name> <input_file> <output_file>
```

Example:
```
./baseline1 input.txt output.txt
```

#### SCM_k-most (Our Algorithm)
The `SCM_k-most` program requires an additional parameter `k` to specify the number of most frequent colors:
```
./scm_k_most <input_file> <output_file> <k>
```

Examples:
```bash
# Compute the mode (k=1, equivalent to original SCM algorithm)
./scm_k_most input.txt output.txt 1

# Compute top-3 most frequent colors
./scm_k_most input.txt output.txt 3
```

---

## Sequence Classification with Random Forest

This repository also includes a sequence classification tool using Random Forest with k-mer features.

#### Requirements

Python 3 with the following packages:
```bash
pip install numpy scikit-learn
```

#### Usage

The script `sequence_classification_rf.py` performs binary classification (SARS-CoV-2 vs. non-SARS-CoV-2) using 7-mer features:

```bash
python sequence_classification_rf.py
```