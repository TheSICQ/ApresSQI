# AprèsSQI

Implementations related to the SQIsign signature scheme, using extension-fields in signing. This is code accompanying the paper *[AprèsSQI: Extra Fast Verification for SQIsign Using Extension-Field Signing](https://eprint.iacr.org/2023/1559)* by Maria Corte-Real Santos, Jonathan Komada Eriksen, Michael Meyer and Krijn Reijnders.

## ApresSQI-Sage
The implementation of the signing and an unoptimized verification procedure is available in the folder ApresSQI-Sage. To run, see README.md in that folder.

## BenchmarkVerification
The code used for benchmarking verification is available in VerificationBenchmarking. To (re)run benchmarking, go into the relevant folder and run

```bash
python benchmark.py All 100
```

Where 'All' denotes which versions to benchmark, and '100' denotes the number of runs per prime.

(tested with python 3.11.1)
