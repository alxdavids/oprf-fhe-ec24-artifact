//! The `prf` module contains all of the core functionality for the PRF
//! and OPRF that are constructed in this work.
//!
//! # Example usage
//!
//! ```no_run
//! # use oprf_fhe::prf::*;
//! # use oprf_fhe::objects::{BinMatrix, sample_bin_vec};
//! // appropriate params choice for concrete
//! use tfhe::shortint::parameters::PARAM_MESSAGE_2_CARRY_0;
//!
//! // The following dimensions are deemed to be secure based on the
//! // security of the Crypto Dark Matter PRF
//! let input_dim = 128;
//! let output_dim = 128;
//! let sk_rows = 257;
//! let sk_cols = 256;
//! let bounded = true;
//!
//! // Construct public parameters
//! let pp = PublicParams::new(
//!   input_dim,
//!   output_dim,
//!   sk_rows,
//!   sk_cols,
//!   &PARAM_MESSAGE_2_CARRY_0,
//!   bounded,
//! );
//! // generate client keyapir
//! let client_keys = client_setup(PARAM_MESSAGE_2_CARRY_0);
//!
//! // server OPRF key
//! let A = BinMatrix::sample(pp.sk_rows, sk_cols);
//! // client input
//! let x = sample_bin_vec(input_dim);
//! // public metadata
//! let t = "some_metadata".as_bytes().to_vec();
//!
//! // client request
//! let req = ClientRequest::new(
//!   &pp,
//!   &client_keys.1,
//!   &client_keys.0,
//!   t.clone(),
//!   x.clone(),
//! );
//! // server blind evaluation
//! let w_enc = server_blind_eval(&pp, &A, &req);
//! // client finalise
//! let z_enc = client_finalise(&pp, &client_keys.0, t, x, w_enc);
//! ```
use crate::objects::*;
use bitvec::prelude::*;
use tfhe::shortint::{gen_keys, Ciphertext, parameters::ClassicPBSParameters as FHEParams};
pub use tfhe::shortint::{
  ClientKey as FHESecretKey, ServerKey as FHEPublicKey,
};
use rayon::prelude::*;
use sha2::{Digest, Sha256, Sha384};

/// The `PublicParams` struct is used for storing the values and
/// parameters that are necessary for instantiating the overall
/// environment.
pub struct PublicParams {
  pub input_dim: usize,
  pub output_dim: usize,
  pub sk_cols: usize,
  pub sk_rows: usize,
  pub unchecked_bin_adds: usize,
  pub unchecked_ter_adds: usize,
  pub n_bootstraps: usize, // number of bootstraps per ciphertext
  pub g_inp: TernaryMatrix,
  pub g_out: TernaryMatrix,
  pub bounded: bool,
}
impl PublicParams {
  pub fn new(
    input_dim: usize,
    output_dim: usize,
    sk_rows: usize,
    sk_cols: usize,
    fhe_params: &FHEParams,
    bounded: bool,
  ) -> Self {
    let ptxt_mod: usize = fhe_params.message_modulus.0;
    // there is no algebraic reason for this restriction, but purely noise management
    let unchecked_bin_adds: usize = ptxt_mod - 1;
    let n_bootstraps_A_mul: usize = (sk_rows / 2 / unchecked_bin_adds) + 1;
    // e.g: ptxt_mod = 4, worst case input: 2, we may do one addition (input is in {0,1})
    let unchecked_ter_adds: usize = ptxt_mod - 3;
    let n_bootstraps_G_out_mul: usize = (sk_cols / unchecked_ter_adds) + 1;
    let g_inp: TernaryMatrix = TernaryMatrix::sample(input_dim, input_dim/2);
    let g_out: TernaryMatrix = TernaryMatrix::sample(sk_cols, output_dim);
    Self {
      input_dim,
      output_dim,
      sk_rows,
      sk_cols,
      unchecked_bin_adds,
      unchecked_ter_adds,
      n_bootstraps: n_bootstraps_A_mul
        + n_bootstraps_G_out_mul,
      g_inp,
      g_out,
      bounded,
    }
  }
}

/// The `client_setup` function generates the FHE keypair that is
/// required for encrypting and computing over their OPRF input
pub fn client_setup(fhe_pp: FHEParams) -> (FHESecretKey, FHEPublicKey) {
  gen_keys(fhe_pp)
}

/// The `ClientRequest` struct generates an encrypted client request
/// based on an input value `x`, public metadata `t`, and the client FHE
/// keypair.
pub struct ClientRequest {
  pub pk: FHEPublicKey,
  pub ct: Vec<Ciphertext>,
  pub t: Vec<u8>,
}
impl ClientRequest {
  pub fn new(
    pp: &PublicParams,
    pk: &FHEPublicKey,
    sk: &FHESecretKey,
    t: Vec<u8>,
    x: Vec<u8>,
  ) -> Self {
    if x.len() != pp.input_dim {
      panic!(
        "Invalid input length: {}, should be: {}",
        x.len(),
        pp.input_dim
      );
    }
    let xv = Vector::from(x);
    let mut y = plain_lut_x1(pp, &xv, &plain_G_INP_mul(&xv, &pp.g_inp));
    // Add extra one entry to vector
    y.append_one();
    let ct = y.encrypt(&sk);
    Self {
      pk: pk.clone(),
      ct,
      t: t.clone(),
    }
  }
}

/// The `server_eval_partial` function generates a partial PRF server
/// key `At`, using the public metadata `t` and a standard server PRF
/// key `A`. It then runs the standard `server_eval` algorithm using
/// `At` as the server PRF key.
pub fn server_eval_partial(
  pp: &PublicParams,
  A: &BinMatrix,
  t: Vec<u8>,
  x: Vec<u8>,
) -> Vec<u8> {
  let At = random_oracle_key(pp, A, &t);
  server_eval(pp, &At, x)
}

/// The `server_eval` runs the PRF detailed in Algorithm 1 of the paper.
/// In other words, it computes the following steps:
///
/// ```
/// // y <-- bit_decompose(G_inp * x (mod 3))
/// // r <-- y * A (mod 2)
/// // w <-- G_out * r (mod 3)  
/// ```
pub fn server_eval(pp: &PublicParams, A: &BinMatrix, x: Vec<u8>) -> Vec<u8> {
  if x.len() != pp.input_dim {
    panic!("Input length x: {}, expected: {}", x.len(), pp.input_dim);
  }
  let xv = Vector::from(x);
  let mut y = plain_lut_x1(pp, &xv, &plain_G_INP_mul(&xv, &pp.g_inp));
  y.append_one();
  let r = plain_lut_x3(pp, &plain_A_x2_mul(&y, A));
  let w = plain_lut_x5(pp, &plain_G_OUT_mul(&r, &pp.g_out));
  w.as_vec().into_iter().map(|&x| x as u8).collect()
}

/// The `client_finalise_clear` function, takes the cleartext output of
/// `server_eval`, and computes a ROM hash function evaluation over it
/// to receive `z`: the final output of the PRF.
pub fn client_finalise_clear(
  _pp: &PublicParams,
  t: Vec<u8>,
  x: Vec<u8>,
  w: Vec<u8>,
) -> Vec<u8> {
  let z = random_oracle_finalise(&t, &x, &w);
  z
}

/// The `server_blind_eval_partial` function runs the blinded version of
/// the `server_eval_partial` algorithm by operating over FHE
/// ciphertexts sent by the client.
pub fn server_blind_eval_partial(
  pp: &PublicParams,
  A: &BinMatrix,
  req: &ClientRequest,
) -> Vec<Ciphertext> {
  let ClientRequest { t, .. } = req;
  let At = random_oracle_key(pp, A, &t);
  server_blind_eval(pp, &At, req)
}

/// The `server_blind_eval` function runs the blinded version of the
/// `server_eval` algorithm by operating over FHE ciphertexts sent by
/// the client.
///
/// In other words, the client sends the encrypted value
/// `y<--bit_decompose(G_inp * x (mod 3))` and the server computes `w`
/// using FHE and the public key provided by the client.
pub fn server_blind_eval(
  pp: &PublicParams,
  A: &BinMatrix,
  req: &ClientRequest,
) -> Vec<Ciphertext> {
  let ClientRequest { pk, ct, .. } = req;
  let r = fhe_A_y_mul_mod2(pp, &pk, &ct, &A);
  let w = fhe_G_OUT_mul_mod3(pp, &pk, &r, &pp.g_out);
  w
}

/// The `client_finalise` function, takes the encrypted output of
/// `server_blind_eval`, decrypts it, and then computes a ROM hash
/// function evaluation over it to receive `z`: the final output of the
/// PRF.
pub fn client_finalise(
  pp: &PublicParams,
  sk: &FHESecretKey,
  t: Vec<u8>,
  x: Vec<u8>,
  rep: Vec<Ciphertext>,
) -> Vec<u8> {
  let w: Vec<u8> = rep.iter().map(|c| {
    sk.decrypt(&c) as u8
  }).collect();
  client_finalise_clear(pp, t, x, w)
}

fn random_oracle_key(pp: &PublicParams, A: &BinMatrix, t: &[u8]) -> BinMatrix {
  let At_cols: Vec<Vec<u64>> = A
    .cols()
    .iter()
    .map(|col| {
      let mut data = vec![];
      data.extend_from_slice(t);
      for &x in col.iter() {
        data.push(x as u8);
      }
      data.extend_from_slice("oprf_fhe_rom_key".as_bytes());
      let digest: Vec<u8> = Sha384::digest(&data).to_vec();
      let bv = BitVec::<_, Lsb0>::from_vec(digest);
      if bv.len() < pp.sk_rows {
        panic!(
          "Bit vector length ({}) is shorter than required ({}).",
          bv.len(),
          pp.sk_rows
        );
      }

      if let Some(bits) = bv.get(0..pp.sk_rows) {
        return bits.to_bitvec().iter().map(|b| *b as u64).collect();
      }
      panic!("Unable to recover bits up to length: {}", pp.sk_rows);
    })
    .collect();
  BinMatrix::new(At_cols, pp.sk_rows, pp.sk_cols)
}

fn random_oracle_finalise(t: &[u8], x: &[u8], w: &[u8]) -> Vec<u8> {
  let mut ro_2 = Sha256::new();
  ro_2.update(t);
  ro_2.update(x);
  ro_2.update(w);
  ro_2.update("oprf_fhe_rom_finalise".as_bytes());
  ro_2.finalize().to_vec()
}

fn plain_G_INP_mul(x0: &Vector, g_inp: &TernaryMatrix) -> Vector {
  x0 * g_inp
}

fn plain_lut_x1(pp: &PublicParams, x0: &Vector, x1: &Vector) -> Vector {
  // we treat G_inp as a input_dim x 3*input_dim/2 matrix, where the first input_dim x input_dim is
  // an identity matrix, and the remaining input_dim x input_dim/2 matrix is a random
  // ternary matrix.
  // 
  // The decomposition is then applied only on the result of multiplying
  // the client input with this second matrix portion, to expand the
  // client input to the required.
  let mut x2 = Vector::from(vec![0u64; 2*pp.input_dim]);
  // first input_dim bits are simply the client input
  for i in 0..pp.input_dim {
    x2.set(i, x0.get(i));
  }
  // second input_dim bits are expanded from the subsequent randomly
  // distributed input_dim/2 bits
  for i in 0..pp.g_inp.width() {
    let x2idx = (2*i) + pp.input_dim;
    match modulo(x1.get(i), 3) {
      0 => {
        x2.set(x2idx, 0);
        x2.set(x2idx + 1, 0);
      }
      1 => {
        x2.set(x2idx, 1);
        x2.set(x2idx + 1, 0);
      }
      2 => {
        x2.set(x2idx, 0);
        x2.set(x2idx + 1, 1);
      }
      _ => panic!("Bad value received: {}", modulo(x1.get(i), 3)),
    }
  }
  x2
}

fn plain_A_x2_mul(x2: &Vector, A: &BinMatrix) -> Vector {
  x2 * A
}

fn plain_lut_x3(pp: &PublicParams, x3: &Vector) -> Vector {
  let mut x4 = Vector::from(vec![0u64; pp.sk_cols]);
  for j in 0..pp.sk_cols {
    x4.set(j, modulo(x3.get(j), 2) as u64);
  }
  x4
}

fn plain_G_OUT_mul(x4: &Vector, g_out: &TernaryMatrix) -> Vector {
  x4 * g_out
}

fn plain_lut_x5(pp: &PublicParams, x5: &Vector) -> Vector {
  Vector::from(
    (0..pp.output_dim)
      .into_iter()
      .map(|i| modulo(x5.get(i), 3) as u64)
      .collect::<Vec<u64>>(),
  )
}

fn fhe_A_y_mul_mod2(
  pp: &PublicParams,
  pk: &FHEPublicKey,
  cx_in: &[Ciphertext],
  A: &BinMatrix,
) -> Vec<Ciphertext> {
  let mod_2 = pk.generate_lookup_table(|x| x % 2);
  let cx_out = (0..A.cols().len())
    .into_par_iter()
    .map(|i| {
      let col = &A.cols()[i];
      let mut entry = pk.create_trivial(0);
      let mut count = 0;
      for j in 0..pp.sk_rows {
        if col[j] == 0 {
          continue;
        } 
        pk.unchecked_add_assign(&mut entry, &cx_in[j]);
        count += 1;
        if count >= pp.unchecked_bin_adds && pp.bounded {
          pk.apply_lookup_table_assign(&mut entry, &mod_2);
          count = 0;
        }
      }
      pk.apply_lookup_table_assign(&mut entry, &mod_2);
      entry
    })
    .collect();
  cx_out
}

fn fhe_G_OUT_mul_mod3(
  pp: &PublicParams,
  pk: &FHEPublicKey,
  cx_in: &[Ciphertext],
  g_out: &TernaryMatrix,
) -> Vec<Ciphertext> {
  let mod_3 = pk.generate_lookup_table(|x| x % 3);
  let cx_out = (0..g_out.cols().len())
    .into_par_iter()
    .map(|i| {
      let col = &g_out.cols()[i];
      let mut entry = pk.create_trivial(0);
      let mut count = 0;
      for j in 0..pp.sk_cols {
        if col[j] == 0 {
          continue;
        }
        pk.unchecked_add_assign(&mut entry, &cx_in[j]);
        count += 1;
        // we need: 2 + count < ptxt_mod
        if count >= pp.unchecked_ter_adds && pp.bounded {
          pk.apply_lookup_table_assign(&mut entry, &mod_3);
          count = 0;
        }
        if col[j] == 2 {
          pk.unchecked_add_assign(&mut entry, &cx_in[j]);
          count += 1;
          if count >= pp.unchecked_ter_adds && pp.bounded {
            pk.apply_lookup_table_assign(&mut entry, &mod_3);
            count = 0;
          }
        }
      }
      pk.apply_lookup_table_assign(&mut entry, &mod_3);
      entry
    })
    .collect();
  cx_out
}

#[cfg(test)]
mod test {
  use super::*;
  use tfhe::shortint::parameters::PARAM_MESSAGE_2_CARRY_0;

  #[test]
  fn test_end_to_end() {
    let input_dim = 16;
    let output_dim = 16;
    let sk_rows = 2*input_dim + 1;
    let sk_cols = 2*input_dim;
    let params = PARAM_MESSAGE_2_CARRY_0;
    let pp = PublicParams::new(
      input_dim,
      output_dim,
      sk_rows,
      sk_cols,
      &params,
      true,
    );
    let client_keys = client_setup(params);

    for _ in 0..4 {
      let A = BinMatrix::sample(pp.sk_rows, sk_cols);
      let x = sample_bin_vec(input_dim);
      // no metadata
      let t = vec![];

      // clear evaluation
      let w_clear = server_eval(&pp, &A, x.clone());
      // oblivious evaluation
      let req = ClientRequest::new(
        &pp,
        &client_keys.1,
        &client_keys.0,
        t.clone(),
        x.clone(),
      );
      let w_enc = server_blind_eval(&pp, &A, &req);

      let z_clear = client_finalise_clear(&pp, t.clone(), x.clone(), w_clear);
      let z_enc = client_finalise(&pp, &client_keys.0, t, x, w_enc);

      assert_eq!(z_clear, z_enc);
    }
  }

  #[test]
  fn test_end_to_end_partial() {
    let input_dim = 16;
    let output_dim = 16;
    let sk_rows = 2*input_dim + 1;
    let sk_cols = 2*input_dim;
    let params = PARAM_MESSAGE_2_CARRY_0;
    let pp = PublicParams::new(
      input_dim,
      output_dim,
      sk_rows,
      sk_cols,
      &params,
      true,
    );
    let client_keys = client_setup(params);

    for _ in 0..4 {
      let A = BinMatrix::sample(pp.sk_rows, sk_cols);
      let x = sample_bin_vec(input_dim);
      let t = "some_metadata".as_bytes().to_vec();

      // clear evaluation
      let w_clear = server_eval_partial(&pp, &A, t.clone(), x.clone());
      // oblivious evaluation
      let req = ClientRequest::new(
        &pp,
        &client_keys.1,
        &client_keys.0,
        t.clone(),
        x.clone(),
      );
      let w_enc = server_blind_eval_partial(&pp, &A, &req);

      let z_clear = client_finalise_clear(&pp, t.clone(), x.clone(), w_clear);
      let z_enc = client_finalise(&pp, &client_keys.0, t, x, w_enc);

      assert_eq!(z_clear, z_enc);
    }
  }
}
