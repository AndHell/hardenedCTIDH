# Hardened CTIDH: Dummy-Free and Deterministic CTIDH
This repository contains auxiliary material for the paper *[Hardened CTIDH: Dummy-Free and Deterministic CTIDH](https://eprint.iacr.org/2025/1645/)*.

Authors:
- [Gustavo~Banegas](https://cryptme.in/) `<gustavo@cryptme.in>`
- [Andreas Hellenbrand](https://www.andhell.de/) `<andreas.hellenbrand@hs-rm.de>`
- Matheus Saldanha `<matheus.saldanha@posgrad.ufsc.br>`


This repository is a fork of ["dCTIDH: Fast & Deterministic CTIDH"](https://eprint.iacr.org/2025/107) by:
- [Fabio Campos](https://www.sopmac.org/) `<campos@sopmac.de>`
- [Andreas Hellenbrand](https://www.andhell.de/) `<andreas.hellenbrand@hs-rm.de>`
- [Michael Meyer](https://www.uni-regensburg.de/informatics-data-science/qpc/team/dr-michael-meyer/index.html) `<michael@random-oracles.org>`
- [Krijn Reijnders](https://krijnreijnders.com/) `<krijn@q1q1.nl>`

# Overview

## Building
We tested our code with GCC 12 on Debian 12.
Furthermore, the implementation makes use of the ADX (ADOX and ADCX) instructions, 
so you need an Intel Broadwell/AMD ZEN CPU or newer.

```sh
# Only necessary first time (generally)
mkdir build && cd build
cmake ..

# If you want with instrumentation for constant-time behavior testing, 
#the default value is OFF. Valgrind development files are used for this build option.
cmake -DENABLE_CT_TESTING=ON ..

# Building
make
```
This builds the executeables for 2 versions:

- 2047m4l205
- 2047m6l194

## benchmarking

### Automated Benchmarking

The project includes automated benchmark targets that make it easy to run and 
analyze benchmarks for all enabled parameter sets:

```sh
# Run benchmarks for a specific parameter set
make benchmark-ctidh-2047m4l205

# Run all benchmarks and display a summary
make benchmark

# Show just the summary of previously run benchmarks 
make benchmark-summary
```

By default, benchmarks run with 100 iterations, which will take several hours. 
You can change this by setting the `SECSIDH_BENCHMARK_RUNS` option:

```sh
# Configure with 5 benchmark runs
cmake -DSECSIDH_BENCHMARK_RUNS=5 ..

```

The benchmark results are saved to files in the build directory:
   - Raw logs: `benchmark-ctidh-<param_set>.log`
   - Analysis results: `benchmark-ctidh-<param_set>-analysis.log`

### Manual Benchmarking

You can also run benchmarks manually using the executable options:
when in `build`:
```sh
usage: 	
    ./main/ctidh-2047m4l205.main                            # for a quick test
	./main//ctidh-2047m4l205.main -bact [number of runs]    # run benchmark for the action
	./main//ctidh-2047m4l205.main -bfp [number of runs]     # run benchmark for fp arithmetic
```

Each version contains benchmarking tools for the action, as well as the finite-field arithmetic,
which can be used with `-bact`, resp. `-bfp`.

The action benchmarks can be analyzed using the `analyze_bench.py` script:
```sh
./main/ctidh-2047m4l205.main -bact 100 > bench_action.out
python3 ../analyze_bench.py < bench_action.out 
```

The analyze_bench.py script supports different output formats:
```sh
# Default grid format for terminal viewing
python3 ../analyze_bench.py < bench_action.out

# CSV format for importing into spreadsheets
python3 ../analyze_bench.py --format=csv < bench_action.out

# LaTeX format
python3 ../analyze_bench.py --format=latex < bench_action.out
```

## constant-time check
If `DENABLE_CT_TESTING=ON`, `checkct` versions of the executable are created 
for all versions, which can be validated with `valgrind`.

when in `build`:
```sh 
cmake -DENABLE_CT_TESTING=ON ..

make  # creates all versions

make checkct-2047m4l205.main # for single version
make checkct-2047m6l194.main

# run valgrind test
valgrind ./main/checkct-2047m4l205.main
valgrind ./main/checkct-2047m6l194.main
```


**Remark**: There seems to be a Valgrind issue with some combinations of GCC versions and modern CPUs due to missing AVX instructions. See the details [here](https://sourceware.org/git/?p=valgrind.git;a=blob;f=docs/internals/3_15_BUGSTATUS.txt;h=88d5466f6b799bf7b57c3ca6be0a269fb82df30f;hb=HEAD#l103).
If you encounter issues, we recommend trying again with GCC 12, as used in our setup.

## parameter searchs
We use greedy to find optimal configurations. The script explors the keyspace 
for primes with 151 to 226 ell_i and 1 to 18 batches.
We recomend to split up the search, as this will take a while 
(up to a month using 4 jobs with 48 threads each).

```sh
cd scripts/greedy/
./greedywombats.py
```


# Licenses

Code in this repository that does not indicate otherwise is placed in the public domain.
The code in this repository is based on [dCTIDH](https://github.com/PaZeZeVaAt/dCTIDH, 
which is based on [secsidh](https://github.com/kemtls-secsidh/secsidh).
Both use the same license as this work: [secsidh License (CC0)](https://github.com/kemtls-secsidh/secsidh/blob/main/code/LICENSE.md)
