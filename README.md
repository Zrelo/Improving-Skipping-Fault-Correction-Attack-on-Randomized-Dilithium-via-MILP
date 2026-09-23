# Improving Skipping Fault Correction Attacks on Randomized Dilithium via MILP

This repository contains the complete runnable implementation and experimental artifacts for the paper:

> **Improving Skipping Fault-Correction Attacks on Randomized Dilithium via MILP**

The implementation is based on the original Dilithium fault-injection framework:

https://github.com/Chair-for-Security-Engineering/dilithium-faults

Our modifications were implemented on **October 3, 2025**. The repository is organized so that the complete runnable project is kept in the directory named:

```text
Improving Skipping Fault Correction Attack on Randomized Dilithium via MILP/
```

## Repository contents

There are two different types of content in this repository:

1. **Complete runnable implementation**

   The directory

   ```text
   Improving Skipping Fault Correction Attack on Randomized Dilithium via MILP/
   ```

   contains the full source tree, build files, Dilithium implementation, attack implementation, and scripts required to compile and run the experiments.

2. **Extracted description of the modifications**

   The root-level file

   ```text
   correction_attacks_changes.cpp
   ```

   is not a separate replacement project. It records the changes made to the original `correction_attacks.cpp` implementation. In particular, it documents the Gurobi-based MILP recovery path and the related recovery-function modifications. It is provided to make the differences from the upstream implementation easier to inspect.

   The complete runnable implementation should be built from the project directory, not from the root-level extracted changes file.

## Directory structure

```text
.
├── LICENSE
├── README.md
├── correction_attacks_changes.cpp
└── Improving Skipping Fault Correction Attack on Randomized Dilithium via MILP/
    ├── README
    ├── Makefile
    ├── run.sh
    ├── correction_attacks.cpp
    ├── correction_attacks.h
    ├── statistics_attack.cpp
    └── dilithium/
        └── ...
```

Within the complete runnable project:

- `correction_attacks.cpp` contains the attack implementation used by the build;
- `statistics_attack.cpp` contains the statistical experiment driver;
- `run.sh` provides the experiment-running entry point;
- `dilithium/` contains the Dilithium implementation and its build files;
- `Makefile` builds the skipping-fault and shuffled-skipping-fault executables.

## Relation to the original implementation

This work is based on the original implementation by Elisabeth Krahmer and Georg Land:

- Repository: https://github.com/Chair-for-Security-Engineering/dilithium-faults
- Relevant original implementation: `skippingFault/correction_attacks.cpp`

The main modification is the replacement or extension of the original recovery procedure with a Gurobi-based integer linear programming (MILP/ILP) recovery path. The extracted file `correction_attacks_changes.cpp` describes these changes, including the solver configuration parameters and the modified recovery functions.

## Requirements

The complete runnable implementation has been developed and tested in a Linux environment with:

- GCC / G++
- GNU Make
- fplll
- NTL
- GMP
- MPFR
- Gurobi, including the C++ and C libraries

The project Makefile currently contains machine-specific paths. Before compiling, update them to match your local installation. For example, the following paths may need to be changed:

```text
/home/haiouc/fplll
/opt/gurobi1300/linux64
```

A valid Gurobi license may also be required to run the MILP experiments.

## Build the complete runnable implementation

Change into the directory containing the complete implementation:

```bash
cd "Improving Skipping Fault Correction Attack on Randomized Dilithium via MILP"
```

Review and update the compiler, include, and library paths in `Makefile`, then build all default targets:

```bash
make
```

To remove generated binaries and intermediate object files:

```bash
make clean
```

## Executables

The default build produces the following programs:

### Plain skipping-fault experiments

```text
test_2_skip
test_3_skip
test_5_skip
```

### Shuffled skipping-fault experiments

```text
test_2_skip_shuff
test_3_skip_shuff
test_5_skip_shuff
```

The numbers correspond to the Dilithium parameter sets used by the experiments.

## Run the experiments

After a successful build, run the programs from the complete project directory, for example:

```bash
./test_2_skip
./test_3_skip
./test_5_skip
```

For the shuffled variants:

```bash
./test_2_skip_shuff
./test_3_skip_shuff
./test_5_skip_shuff
```

The `run.sh` script in the same directory can also be used as the experiment entry point. Inspect the script before running it and adapt it to your local environment if necessary.

## Reproducibility

For results reported in the paper, record the following information together with the repository version or commit:

- operating system;
- GCC and G++ versions;
- compiler and optimization flags;
- Gurobi version;
- fplll, NTL, GMP, and MPFR versions;
- Gurobi license and thread configuration;
- values of the experiment parameters;
- random seed, if applicable;
- number of trials;
- hardware and approximate runtime.

Because MILP solver versions, hardware, and solver settings can affect runtime and search behavior, small differences in execution time may occur across systems.

## License and attribution

The original work and the modifications in this repository are released under the Creative Commons Attribution 4.0 International License (CC BY 4.0). See `LICENSE` for the applicable attribution and license information.

When reusing or redistributing this repository, please:

- cite the associated publication;
- credit the original implementation by Elisabeth Krahmer and Georg Land;
- credit the modifications made for this work;
- preserve the relevant copyright and license notices.

## Citation

If you use this implementation or the experimental artifacts, please cite the associated publication and specify the repository commit or release used for your results.

## Responsible use

This repository is provided for academic research, security evaluation, and reproducibility. Use it only in environments for which you have explicit authorization and comply with all applicable laws, regulations, and institutional policies.
