# Software Artifact: "Crypto Dark Matter on the Torus: Oblivious PRFs from shallow PRFs and FHE"

**Description**: Software implementations of the oblivious pseudorandom function protocol

- Sage (`sage/`): Near fully-featured implementation of protocol, demonstrating functionality and correctness of approach.
- Rust (`rust/`): Benchmarking implementation for estimating performance of core oblivious pseudorandom function protocol.

**Paper** (full version): <https://eprint.iacr.org/2023/232>

## Authors

- [Martin R. Albrecht](https://malb.io)
- [Alex Davidson](https://alxdavids.xyz)
- [Amit Deo](https://scholar.google.com/citations?user=TPREbisAAAAJ&hl=en)
- [Daniel Gardham](https://www.surrey.ac.uk/people/daniel-gardham)

## Instructions for use

- [Sage implementation](sage/README.md)
- [Rust benchmarking implementation](rust/README.md)

## Missing functionality

The following summarises missing functionality from our implementations.

**Sage code**:
- Zero-knowledge proofs (i.e. not verifiable)

**Rust code**:
- Zero-knowledge proofs (i.e. not verifiable)
- Circuit-private bootstrapping (since non-power-of-two `q` is not supported in tfhe-rs v0.3.1)
- Depth-one correctness (since modified test polynomials and `p != 3` are not supported in tfhe-rs v0.3.1)
- Public-key compression
