//! This binary runs an end-to-end test, and provides information on a
//! number of performance criteria, including: size of bootstrapping
//! key, number of bootstrappings, ciphertext size, and peak RAM usage.
#![allow(non_snake_case)]
use std::fs;
use std::path::Path;

use tfhe::shortint::Ciphertext;
use oprf_fhe::objects::{sample_bin_vec, BinMatrix};
use oprf_fhe::prf::{
  client_finalise, client_setup, server_blind_eval_partial, ClientRequest,
  FHEPublicKey, FHESecretKey, PublicParams,
};

fn main() {
  let input_dim = 128;
  let output_dim = 128;
  let sk_rows = input_dim * 2 + 1;
  let sk_cols = input_dim * 2;
  let params = oprf_fhe::params::PARAM_100;
  let bounded = std::env::var("BOUNDED").is_ok();
  let pp =
    PublicParams::new(input_dim, output_dim, sk_rows, sk_cols, &params, bounded);
  let data_path = Path::new("./data");
  let sk_path = Path::new("./data/private");
  let pk_path = Path::new("./data/public");
  if !data_path.exists() {
    fs::create_dir(data_path).unwrap();
  }
  let server_only = std::env::var("SERVER_ONLY").is_ok();
  let check_output = std::env::var("CHECK").is_ok();

  let (sk, pk, serialized_sk, serialized_pk) = if !sk_path.exists() {
    println!("Generating client key pair...");
    let client_keys = client_setup(params);
    let serialized_sk = bincode::serialize(&client_keys.0).unwrap();
    fs::write(sk_path, &serialized_sk).unwrap();
    let serialized_pk = bincode::serialize(&client_keys.1).unwrap();
    fs::write(pk_path, &serialized_pk).unwrap();
    println!(
      "Wrote client key pair to {}",
      sk_path.parent().unwrap().to_str().unwrap()
    );
    (
      Some(client_keys.0),
      client_keys.1,
      Some(serialized_sk),
      serialized_pk,
    )
  } else {
    println!("Reading client key pair from file...");
    let pk_bytes = fs::read(pk_path).unwrap();
    let pk: FHEPublicKey = bincode::deserialize(&pk_bytes).unwrap();
    let (sk_wrap, sk_bytes_wrap) = if !server_only {
      let sk_bytes = fs::read(sk_path).unwrap();
      let sk: FHESecretKey = bincode::deserialize(&sk_bytes).unwrap();
      (Some(sk), Some(sk_bytes))
    } else {
      (None, None)
    };
    (sk_wrap, pk, sk_bytes_wrap, pk_bytes)
  };

  // generate server OPRF key
  println!("Generating server OPRF key...");
  let oprf_key = BinMatrix::sample(pp.sk_rows, sk_cols);
  // generate client input
  let client_input = sample_bin_vec(input_dim);
  println!("Client input: {:?}", client_input.clone());
  // some metadata tag
  let t: Result<String, String> =
    std::env::var("METADATA").or(Ok("some_metadata".to_string()));
  let md_str = t.unwrap();
  println!("Metadata used: {}", md_str);
  let md_bytes = md_str.as_bytes().to_vec();

  // oblivious evaluation
  let req_path = Path::new("./data/ciphertexts");
  let (req, serialized_req) = if !server_only {
    println!("Creating client request...");
    let req = ClientRequest::new(
      &pp,
      &pk,
      sk.as_ref().unwrap(),
      md_bytes.clone(),
      client_input.clone(),
    );
    let serialized_req = bincode::serialize(&req.ct).unwrap();
    fs::write(req_path, &serialized_req).unwrap();
    (req, serialized_req)
  } else {
    println!("Reading client request from file...");
    let serialized_req = fs::read(req_path).unwrap();
    let ct: Vec<Ciphertext> = bincode::deserialize(&serialized_req).unwrap();
    let req = ClientRequest {
      ct,
      pk,
      t: md_bytes.clone(),
    };
    (req, serialized_req)
  };

  let w_enc = if std::env::var("CLIENT_ONLY").is_ok() {
    // this is a hack to give the client something to finalize
    println!("Ignoring server computation...");
    req.ct
  } else {
    println!("Running server blind evaluation...");
    server_blind_eval_partial(&pp, &oprf_key, &req)
  };
  let serialized_resp = bincode::serialize(&w_enc).unwrap();

  if !server_only {
    println!("Running client finalise...");
    let z_enc =
      client_finalise(&pp, sk.as_ref().unwrap(), md_bytes.clone(), client_input.clone(), w_enc);
    println!("Client PRF output: {:?}", z_enc);

    // check plain evaluation
    if check_output {
      println!("Checking final output...");
      let w_clear = oprf_fhe::prf::server_eval_partial(&pp, &oprf_key, md_bytes.clone(), client_input.clone());
      let z_clear = oprf_fhe::prf::client_finalise_clear(&pp, md_bytes, client_input, w_clear);
      assert_eq!(z_clear, z_enc);
      println!("Success!");
    }
  }

  println!("\n#### Sizes (bytes) ####\n");
  if !server_only {
    println!(
      "\t* Total size of client secret key: {:?}",
      serialized_sk.unwrap().len()
    );
  }
  println!(
    "\t* Total size of client public key: {:?}",
    serialized_pk.len()
  );
  println!(
    "\t* Total size of client request (ciphertexts): {:?}",
    serialized_req.len()
  );
  println!(
    "\t* Total size of server response: {:?}",
    serialized_resp.len()
  );
  println!(
    "\t* Total number of bootstraps (per ciphertext): {:?}",
    pp.n_bootstraps
  );
}
