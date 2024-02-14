use criterion::{criterion_group, criterion_main, Criterion};

use oprf_fhe::objects::{sample_bin_vec, BinMatrix};
use oprf_fhe::prf::{
  client_finalise, client_setup, server_blind_eval, server_blind_eval_partial,
  ClientRequest, PublicParams,
};

struct ParamSet {
  p: tfhe::shortint::parameters::ClassicPBSParameters,
  n: String,
}

fn criterion_benchmark(c: &mut Criterion) {
  let mut group = c.benchmark_group("FHE OPRF benchmarks");
  group.sample_size(10);

  // Secure parameters
  let input_dim = 128;
  let output_dim = 128;
  let sk_rows = 2*input_dim + 1;
  let sk_cols = 2*input_dim;
  let bench_params = [ParamSet {p: oprf_fhe::params::PARAM_100, n: String::from("100")}, ParamSet {p: oprf_fhe::params::PARAM_128, n: String::from("128")}];

  for params in bench_params {
    // Generate public parameters and client keys at the start
    let pp =
      PublicParams::new(input_dim, output_dim, sk_rows, sk_cols, &params.p, false);
    let keys = client_setup(params.p);
    let t: Vec<u8> = "public_metadata".as_bytes().to_vec();

    group.bench_function(format!("Client: generate encrypted request (λ: {:?})", params.n), |b| {
      let x = sample_bin_vec(pp.input_dim);
      b.iter(|| ClientRequest::new(&pp, &keys.1, &keys.0, t.clone(), x.clone()));
    });

    group.bench_function(format!("Server: blind evaluation (λ: {:?})", params.n), |b| {
      let oprf_key = BinMatrix::sample(pp.sk_rows, pp.sk_cols);
      let x = sample_bin_vec(pp.input_dim);
      let req = ClientRequest::new(&pp, &keys.1, &keys.0, vec![], x);
      b.iter(|| server_blind_eval(&pp, &oprf_key, &req));
    });

    group.bench_function(format!("Server: blind evaluation (with metadata) (λ: {:?})", params.n), |b| {
      let oprf_key = BinMatrix::sample(pp.sk_rows, pp.sk_cols);
      let x = sample_bin_vec(pp.input_dim);
      let req = ClientRequest::new(&pp, &keys.1, &keys.0, t.clone(), x);
      b.iter(|| server_blind_eval_partial(&pp, &oprf_key, &req));
    });

    group.bench_function(format!("Client: finalise (λ: {:?})", params.n), |b| {
      let oprf_key = BinMatrix::sample(pp.sk_rows, pp.sk_cols);
      let x = sample_bin_vec(pp.input_dim);
      let req = ClientRequest::new(&pp, &keys.1, &keys.0, t.clone(), x.clone());
      let rep = server_blind_eval_partial(&pp, &oprf_key, &req);
      b.iter(|| client_finalise(&pp, &keys.0, t.clone(), x.clone(), rep.clone()));
    });

    group.bench_function(format!("Client: generate keys (λ: {:?})", params.n), |b| {
      b.iter(|| client_setup(params.p))
    });
  }
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);
