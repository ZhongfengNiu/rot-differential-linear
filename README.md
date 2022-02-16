## Practical (Rotational) Differential-Linear Distinguishers of ARX Ciphers with Arbitrary Output Linear Masks

This repository contains the source code of the paper *Practical (Rotational) Differential-Linear Distinguishers of ARX Ciphers with Arbitrary Output Linear Masks*. 

- The SageMath Notebook file [Notations.ipynb](https://github.com/rdlattack/rot-differential-linear/blob/main/Notations.ipynb) helps the reader to familiarize the notations employed in the paper. 

- The SageMath Notebook file [R-DL-Modulo-Addition.ipynb](https://github.com/rdlattack/rot-differential-linear/blob/main/R-DL-Modulo-Addition.ipynb) is used to compute the exact (rotational) differential-linear correlations of modulo additions. We note that the results are exact and in rational numbers. In the paper, we some times convert these results to floating point numbers to save space.
- The python scripts used to compute the correlations of (rotational) differential-linear correlations of ARX ciphers like ``Alzette``, ``Chacha``, ``SipHahs``, and ``Speck`` can be found in the directory ``Compute-RDL``. 
- Codes used for verification are placed in the  directory ``Crypoto_verify``.
