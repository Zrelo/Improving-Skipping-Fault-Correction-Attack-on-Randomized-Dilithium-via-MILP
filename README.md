# Improving Skipping Fault-Correction Attacks on Randomized Dilithium via MILP

This repository contains the implementation and experimental artifacts for the paper:

> Improving Skipping Fault-Correction Attacks on Randomized Dilithium via MILP

The code is based on the original Dilithium fault-injection framework:

https://github.com/Chair-for-Security-Engineering/dilithium-faults

The present repository contains modifications developed on 2025-10-03 and is intended to support reproduction and extension of the corresponding experimental results.

## Overview

This project extends the original Dilithium implementation with fault-injection and correction-attack experiments targeting randomized signing. The code includes the modified attack logic and the evaluation setup used to assess the impact of skipping faults and related correction strategies.

## Repository structure

```text
.
├── LICENSE
├── README.md
├── Makefile
├── correction_attacks.cpp
├── correction_attacks_changes.cpp
├── statistics_attack.cpp
├── dilithium/
│   └── ...
├── Improving Skipping Fault Correction Attack on Randomized Dilithium via MILP/
│   └── ...
└── other generated artifacts and experiment files
```

## License and attribution

This repository is released under the Creative Commons Attribution 4.0 International License (CC BY 4.0).

See `LICENSE` for the full text.

The original upstream implementation was developed by Elisabeth Krahmer and Georg Land and is attributed as follows:

- https://github.com/Chair-for-Security-Engineering/dilithium-faults
- Original license statement: Creative Commons Attribution 4.0 International License (CC BY 4.0)

Please preserve attribution and license notices when reusing or redistributing this code.

## Requirements

This project has been developed and tested in a Linux environment with the following tools and libraries:

- GCC / G++
- Make
- fplll
- NTL
- GMP
- MPFR
- Gurobi (C++ and C API)

The root `Makefile` includes hard-coded local library paths such as:

- `/home/haiouc/fplll`
- `/opt/gurobi1300/linux64`

Before building, update these paths in the `Makefile` to match your local installation.

## Build

From the repository root, run:

```bash
make
```

This builds the attack binaries used for the skipping-fault experiments.

## Main targets

The default `make` target builds the following attack executables:

- `test_2_skip`
- `test_3_skip`
- `test_5_skip`
- `test_2_skip_shuff`
- `test_3_skip_shuff`
- `test_5_skip_shuff`

To remove generated binaries and intermediate objects:

```bash
make clean
```

## Running the experiments

After compilation, run the generated binaries directly from the repository root, for example:

```bash
./test_2_skip
./test_3_skip
./test_5_skip
```

and for shuffled variants:

```bash
./test_2_skip_shuff
./test_3_skip_shuff
./test_5_skip_shuff
```

The exact runtime, outputs, and parameter settings depend on the local environment and the selected solver configuration.

## Reproducibility notes

For reproducibility, the following should be documented alongside the paper:

- compiler version
- optimization flags
- operating system
- Gurobi version
- fplll / NTL / GMP / MPFR versions
- random seed or experiment configuration
- number of trials used for statistical evaluation

## Citation

If you use this code in a paper or project, please cite the associated publication and include the repository URL and commit/version used.

## Contact

For questions regarding the implementation or reproduction of the results, please contact the repository owner or the corresponding author of the associated paper.
