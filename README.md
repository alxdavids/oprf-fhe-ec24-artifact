# FHE OPRF implementation

Rust implementation of FHE-based Oblivious Pseudorandom Function
protocol functionality, using the
[tfhe-rs](https://github.com/zama-ai/tfhe-rs)
FHE crate.

## Quickstart

First, update your version of rust to the latest version:
```
rustup update
```

The following commands can be used for building and running the code.
Note that the `--release` flag is used to speed up compiled code.

Build:
```
cargo build --release
```

Run tests:
```
cargo test --release
```

Open documentation in web browser:
```
cargo doc --open --no-deps
```

Benchmarking:
```
make bench-paper
```

## Example run

Debug run (including output of various runtime parameters):
```
<OPTIONAL_ENV_VARIABLE>=<VALUE> cargo run --release
```

The `cargo run` command can be run with the following environment variables:

- `METADATA`: changes the metadata string that is used (default: "some_metadata");
- `CLIENT_ONLY`: boolean to trigger only client-side functionality;
- `SERVER_ONLY`: boolean to trigger only server-side functionality.

The `cargo run` command also saves various cryptographic material to the
`data/` folder to speed up subsequent run-throughs. Delete this folder
to trigger sampling of a new key.

## Parallelism

By default, this implementation will attempt to parallelize vector-matrix
multiplications, by processing each matrix column in a separate thread.
If you would like to change the number of threads that should be used
you prefix your `cargo` command with `RAYON_NUM_THREADS=x` where `x` is
the number of permitted threads.

