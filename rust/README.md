# Rust implementation <!-- omit in toc -->

Instructions on how to run Rust benchmarking software, using [tfhe-rs](https://github.com/zama-ai/tfhe-rs/releases/tag/0.3.1) (version 0.3.1).

- [Environment setup](#environment-setup)
- [Building, testing, documentation](#building-testing-documentation)
- [Benchmarking](#benchmarking)
- [Example run](#example-run)

## Environment setup

First, [install](https://www.rust-lang.org/tools/install) and/or update your version of rust to the latest version:
```bash
rustup update
```

## Building, testing, documentation

The following commands can be used for building and running the code.
Note that the `--release` flag is used to speed up compiled code.

If using an ARM processor, please change the tfhe feature "x86_64-unix" to "aarch64-unix" in [Cargo.toml](https://github.com/alxdavids/oprf-fhe-ec24-artifact/blob/main/rust/Cargo.toml).

Build code binary:
```bash
cargo build --release
```

Run all tests:
```bash
cargo test --release
```

Open documentation in web browser:
```bash
cargo doc --open --no-deps
```

## Benchmarking

Our benchmarks include performance measurements for client query generation, server homomorphic PRF evaluation, and client finalisation of the output.

To benchmark the implementation (**WARNING: this could take a long time, and require large amounts of CPU and memory**), you can use the following command:
```bash
make bench
```

The default number of threads used is 64. To alter the number of threads (e.g. to 16) that are used in the benchmark, you can use:
```bash
make THREADS=16 bench
```

The benchmarking tool will generate a text file, where a standard output for a single piece of functionality takes the form below.
```
FHE OPRF benchmarks/Client: generate encrypted request (λ: "100")
                        time:   [32.467 ms 32.755 ms 33.081 ms]
```
The numbers above can be interpreted as minimum, average, maximum. Each benchmark will run for `λ: "100"` and `λ: "128"` bits of security.

To run the benchmarks in the same way that we used for the submission (i.e. comparing 64 threads vs single-threaded execution), use the command:
```bash
make bench-paper
```

## Example run

Debug run (including output of various runtime parameters):
```bash
<OPTIONAL_ENV_VARIABLE>=<VALUE> cargo run --release
```

The `cargo run` command can be run with the following environment variables:

- `METADATA`: changes the metadata string that is used (default: "some_metadata");
- `CLIENT_ONLY`: boolean to trigger only client-side functionality;
- `SERVER_ONLY`: boolean to trigger only server-side functionality.

The `cargo run` command also saves various cryptographic material to the
`data/` folder to speed up subsequent run-throughs. Delete this folder
to trigger sampling of a new key.
