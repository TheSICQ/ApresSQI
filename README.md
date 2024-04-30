# AprèsSQI

This repository contains the code related to the SQIsign signature scheme, using extension-fields in signing. This is code accompanying the paper *[AprèsSQI: Extra Fast Verification for SQIsign Using Extension-Field Signing](https://eprint.iacr.org/2023/1559)* by Maria Corte-Real Santos, Jonathan Komada Eriksen, Michael Meyer and Krijn Reijnders.

All code in this repository was written in `Python` and [`SageMath`](https://www.sagemath.org/), and is made publicly available under the MIT license. The dependencies are:
- Python 3.11.1 (with matplotlib)
- SageMath version 10.1

## Contents

This repository contains two main directories:
- `ApresSQI-Sage`: this contains the implementation of the signing (using extension fields) and an unoptimized verification procedure. The code in this directory is written in `SageMath` and `python`.
- `BenchmarkVerification`: this contains the code used for benchmarking verification. The code in this directory is written in `python`.

We now describe the files in each of these directories and how to run their contents. 

### ApresSQI-Sage

The directory `ApresSQI-Sage` provides a proof of concept implementation of the SQIsign signature algorithm, using field extensions when signing, allowing for lightning fast verification. 


The files contained in the main directory are as follows:
- `ApresSQI.py`: contains a class for the SQIsign digital signature scheme, using multiple field extensions when signing, as well as classes for the SQIsign public and private keys. 
- `benchmark.py`: this file can be run to time signing and (re)obtain the benchmarks given in Table 1 of the accompanying paper. Simply run 
```bash
sage benchmark.py
```
- `ec.py`: contains algorithms involving elliptic curves and their arithmetic
- `example.py`: this file contains an example that can be run as per the instructions below. 
- `id2iso.py`: contains functions used to perform ideal to isogeny translations via the Deuring correspondence. 
- `klpt.py`: contains functions to run the [KLPT algorithm](https://www.cambridge.org/core/journals/lms-journal-of-computation-and-mathematics/article/on-the-quaternion-def-xmlpi-1def-mathsfbi-1boldsymbol-mathsf-1let-le-leqslant-let-leq-leqslant-let-ge-geqslant-let-geq-geqslant-def-pr-mathit-prdef-fr-mathit-frdef-rey-mathit-reell-isogeny-path-problem/3E24D9DC2C6D30684EDFC647AD254F30).
- `param_loader.py`: the main function in this file loads the parameters needed for a prime $p$. These parameters should be precomputed and found in the sub-directory `Precomputed`, otherwise the function will output an error. 
- `quaternion.py`: contains functions needed to perform quaternion arithmetic. 

The directory `ApresSQI-Sage` also contains two sub-directories:
- `Precomputed`: contains precomputed parameters for certain primes $p$ that are loaded by `param_loader.py`. The necessary parameters for a particular prime $p$ can be precomputed using the file `precompute.py`.
- `Privkeys`: contains private keys for certain primes $p$ that are loaded by SQIsign signature call in `ApresSQI.py`

We give the precomputed parameters in `Precomputed` and private keys in `Privkeys` for the following primes:
- `toy`: $2$<sup>$33$</sup>$3$<sup>$10$</sup> - 1
- `NIST`: $2$<sup>$75$</sup>$3$<sup>$36$</sup>$23$<sup>$2$</sup>$59$<sup>$2$</sup>$101$<sup>$2$</sup>$109$<sup>$2$</sup>$197$<sup>$2$</sup>$491$<sup>$2$</sup>$743$<sup>$2$</sup>$1913$<sup>$2$</sup> - 1 (denoted as $p$<sub>$1973$</sub> in the accompanying paper)
- `7-block`: $2$<sup>$145$</sup>$3$<sup>$9$</sup>$59$<sup>$3$</sup>$311$<sup>$3$</sup>$317$<sup>$3$</sup>$503$<sup>$3$</sup> - 1 (denoted as $p$<sub>$7$</sub> in the accompanying paper)
- `4-block`: $2$<sup>$242$</sup>$* 3 * 67 - 1$ (denoted as $p$<sub>$4$</sub> in the accompanying paper)

If other primes are found in these sub-directories, they have been requested and added for other applications (beyond AprèsSQI).


#### Example usage

```python
from ApresSQI import SQIsign

# Signer
msg = "Kani came in like a wrecking ball"
signer = SQIsign()
pk = signer.KeyGen()
signer.sk.export('priv.key')

sigma = signer.Sign(msg)

# Verifier
verifier = SQIsign()
verified = verifier.verify(msg, sigma, pk)
print('Signature was ' + ('CORRECT' if verified else 'WRONG'))
```

Or, if you just want to test a ready-made example, you can run
```bash
sage example.py param seeds compressed
```
which allows for signing with any of the four parameters by setting `param` = `toy`, `NIST`, `7-block`, `4-block`, signing both with and without seeds by setting `seeds` = `True` or `False` (respectively), as well as using compressed or uncompressed signatures by setting `compressed` = `True` or `False` (respectively). 

By default, running
```bash
sage example.py
```
will use the `toy` parameter, signing with seeds and using compressed signatures. 

#### Acknowledgements
The code in this directory is based on the [documentation](http://sqisign.org/spec/sqisign-20230601.pdf) from the NIST submission of SQIsign. Further, it uses many ideas introduced in [Deuring for the People](https://eprint.iacr.org/2023/106), and some code found here is also taken from the [code](https://github.com/friends-of-quaternions/deuring), made by Jonathan Komada Eriksen, Lorenz Panny, Jana Sotáková and Mattia Veroni. Some other code (specifically, lattice enumeration) is also taken from the [code](https://github.com/LearningToSQI/SQISign-SageMath) accompanying the project [Learning to SQI](https://learningtosqi.github.io/), made by Maria Corte-Real Santos, Jonathan Komada Eriksen, Michael Meyer and Giacomo Pope. Finally, the implementation of sqrt_velu in xonly_velusqrt.py is from the SageMath source code, adapted to work with our custom xPoint class (doing x-only arithmetic).

### Benchmark Verification

The directory `BenchmarkVerification` allows us to (re)run benchmarking found in Section 7 of the paper. It contains two subdirectories:

- `Apres`: benchmarks ApresSQI variants of SQIsign
- `NISTandLWXZ`: benchmarks NIST and LWXZ variants of SQIsign

#### AprèsSQI directory

The files in the `Apres` sub-directory are as follows:
- `all_primes.py`: contains a list of primes used to benchmark.
- `ApresSQI.py`: contains the AprèsSQI class for SQIsign verification.
- `benchmark.py`: this is the main file used to run benchmarks for the 4 variants of AprèsSQI verification. We detail how to use this to run benchmarks below.
- `ec_misc.py`: contains miscellaneous functions relating to elliptic curves and their arithmetic
- `fp_arith.py`: contains the functions needed to perform operations in **F**<sub>p</sub> and **F**<sub>p<sup>2</sup></sub>. This contains fast square-root computation taken from the [SuperSolver code](https://github.com/microsoft/SuperSolver) following Scott's [*Tricks of the Trade*](https://eprint.iacr.org/2020/1497).
- `isogenychains.py`: contains functions needed to compute isogenies from the torsion basis.
- `mont_xonly.py`: contains functions needed to do x-only Montgomery scalar multiplication.
- `precomp.py`: contains functions needed in to perform the pre-computation.
- `strategies_sike.py`: contains code taken from [SIBC](https://github.com/JJChiDguez/sibc), implementing the algorithm for computing optimal strategies from the SIKE specification. 
- `testing.py`: this file can be used for testing purposes
- `tors_basis.py`: contains functions to compute the torsion basis.

To benchmark ApresSQI variants of SQIsign go into the sub-directory `Apres` and run

```bash
python benchmark.py version N
```

where 
- `version`  = `APRESwswp` (variant with seeds and pushing), `APRESwsnp` (variant with pushing and no seeds), `APRESnswp` (variant with no seeds and with pushing), `APRESnsnp` (variant with no seeds and no pushing), `uncompressed`, or `All` (runs all variants available). 
See Appendix C.1 of the accompanying paper for more details on the AprèsSQI variants.
- `N` denotes the number of runs per prime.

**NOTE:** Throughout, by `python` we mean `python3` if your terminal requires this (e.g., some MacOS versions may require this.) Furthermore, if you want to output plots, uncomment the line `import matplotlib.pyplot as plt`.

If you want to run the benchmarks for a specific prime, add the prime to the array `SpecificPrimes` on line 343, and then run 
```bash
python benchmark.py version N True
```
where `version` and `N` as above. 

#### NIST and LWXZ directory

The files contained in this directory have identical descriptions barring the following differences:
- `benchmark.py`: this is the main file used to run benchmarks for the NIST and LWXZ variants of SQIsign verification. We detail how to use this to run benchmarks below.
- `fp_arith.py`: the same as above except does not use accelerated square-root computations.
- `SQIsign.py`: this replaces the `ApresSQI.py` above, and contains the class for SQIsign verification used in NIST (and LWXZ).


To benchmark NIST and LWXZ variants of SQIsign go into the sub-directory `NISTandLWXZ` and run

```bash
python benchmark.py version N
```

where 
- `version`  = `NIST` (SQIsign as per the [NIST submission](https://sqisign.org/)), `LWXZ` (SQIsign as per [Lin, Wang, Xu and Zhao](https://eprint.iacr.org/2023/753)), or `All` (runs all variants available).
- `N` denotes the number of runs per prime.

If you want to run the benchmarks for a specific prime, add the prime to the array `SpecificPrimes` on line 162, and then run 
```bash
python benchmark.py version N True
```
where `version` and `N` as above. 

#### Example usage

For example, to benchmark the NIST version of SQIsign with 100 runs per prime, navigate to the `NISTandLWXZ` subdirectory and run:

```bash
python benchmark.py NIST 100
```

## Licenses

Code in this repository that does not indicate otherwise is published under the MIT license.

For the third party code see their licenses:
- [Deuring for the People](https://github.com/friends-of-quaternions/deuring): [https://github.com/friends-of-quaternions/deuring/blob/main/LICENSE.txt](https://github.com/friends-of-quaternions/deuring/blob/main/LICENSE.txt)
- [LearningToSQI](https://github.com/LearningToSQI/SQISign-SageMath): [https://github.com/LearningToSQI/SQISign-SageMath/blob/main/LICENCE](https://github.com/LearningToSQI/SQISign-SageMath/blob/main/LICENCE)
- [SageMath](https://github.com/sagemath/sage): [https://github.com/sagemath/sage?tab=License-1-ov-file](https://github.com/sagemath/sage?tab=License-1-ov-file)
- [sibc](https://github.com/JJChiDguez/sibc): [https://github.com/JJChiDguez/sibc/blob/master/LICENSE](https://github.com/JJChiDguez/sibc/blob/master/LICENSE)
- [SuperSolver](https://github.com/microsoft/SuperSolver): [https://github.com/microsoft/SuperSolver/blob/main/LICENSE](https://github.com/microsoft/SuperSolver/blob/main/LICENSE)
