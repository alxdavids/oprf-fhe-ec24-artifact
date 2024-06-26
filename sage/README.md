# Sage implementation

Instructions on how to use sage implementation.

## Functionality

- [oprf.py](./oprf.py): Core code for running oblivious pseudorandom function protocol
- [tfhe.py](./tfhe.py): Reimplementation of TFHE with modified test polynomial, for allowing modular reduction by non-power-of-two moduli
- [cpbs.py](./cpbs.py): Circuit-private bootstrapping functionality
- [compression.py](./compression.py): Public-key compression functionality
- [gadget.py](./gadget.py): Gadget matrix utility functions

## Example run of OPRF

```bash
sage # enter sage environment
┌────────────────────────────────────────────────────────────────────┐
│ SageMath version 10.2, Release Date: 2023-12-03                    │
│ Using Python 3.11.8. Type "help()" for help.                       │
└────────────────────────────────────────────────────────────────────┘
sage: set_random_seed(1337)
sage: from tfhe import LWE
sage: from oprf import OPRF
sage: oprf = OPRF(LWE(4, 3*7681, "binary", 3.0, p=2))
sage: oprf([0]*8) # Standard non-blinded evaluation of PRF
(1, 0, 2, 2, 1, 1, 0, 2, 2, 1, 2, 2) # Evaluation results
sage: c = oprf.blind_eval([0]*8) # Run blind evaluation of PRF over encrypted ciphertexts
sage: vector([oprf.msbs.lwe_o.decrypt(c_) for c_ in c]) # Decrypt evaluated ciphertexts
(1, 0, 2, 2, 1, 1, 0, 2, 2, 1, 2, 2) # Should be equivalent to non-blinded evaluation results
```


## Running tests

``` bash
PYTHONPATH=`pwd` sage -tp *.py
```

## Sage version

We tested our code with SageMath 10.2, installed in a Conda environment on a Debian GNU/Linux trixie system. 
