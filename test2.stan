// generated with brms 2.19.0
functions {
  // Emission likelihood given that the true state is zero
  real emission_0_fp(array[] real zl) {
    // zl gives the likelihood of the observation given a true zero.
    real out = sum(log(zl)); // the likelihood when the true history is all zeros
    return out;
  }
  // emission likelihood given that state equals one
  real emission_1_fp(row_vector det, array[] real zl, array[] real ol) {
    // det gives logit-detection probabilities
    // zl gives the likelihood of the observation given a true zero
    // ol gives the likelihood of the observation given a true one
    
    int n = size(det); // number of reps
    
    real out = 0;
    
    for (i in 1 : n) {
      real ll_true_zero = log(zl[i]) + bernoulli_logit_lpmf(0 | det[i]);
      real ll_true_one = log(ol[i]) + bernoulli_logit_lpmf(1 | det[i]);
      out += log_sum_exp(ll_true_zero, ll_true_one);
    }
    
    return out;
  }
  real occupancy_single_fp_lpdf(vector fp, // the fp probs
                                vector mu,
                                // lin pred for detection
                                vector occ,
                                // lin pred for occupancy. Elements after vint1[1] irrelevant.
                                array[] int vint1,
                                // # units (n_unit). Elements after 1 irrelevant.
                                array[] int vint2,
                                // # sampling events per unit (n_rep). Elements after vint1[1] irrelevant.
                                array[] int vint3,
                                // Indicator for > 0 certain detections (Q). Elements after vint1[1] irrelevant.
                                array[] int vint4,
                                // number of blocks. elements after 1 irrelevant
                                array[] int vint5,
                                // integer identifying the block 
                                // indices for jth repeated sampling event to each unit (elements after vint1[1] irrelevant):
                                array[] int vint6, array[] int vint7,
                                array[] int vint8, array[] int vint9,
                                array[] real vreal1,
                                // expected fraction of obs 1s that are true. Elements after vint4[1] irrelevant
                                array[] real vreal2, // P
                                array[] real vreal3 // QQ
                                ) {
    // Create array of the rep indices that correspond to each unit.
    array[vint1[1], 4] int index_array;
    index_array[ : , 1] = vint6[1 : vint1[1]];
    index_array[ : , 2] = vint7[1 : vint1[1]];
    index_array[ : , 3] = vint8[1 : vint1[1]];
    index_array[ : , 4] = vint9[1 : vint1[1]];
    
    // find expected number of real ones
    vector[vint4[1]] eno;
    vector[vint4[1]] total_size;
    for (i in 1 : vint4[1]) {
      eno[i] = 0;
      total_size[i] = 0;
    }
    
    for (i in 1 : vint1[1]) {
      array[vint2[i]] int indices = index_array[i, 1 : vint2[i]];
      eno[vint5[i]] += sum(inv_logit(mu[indices])) * inv_logit(occ[i]);
      total_size[vint5[i]] += vint2[i];
    }
    
    // expected number of real zeros
    vector[vint4[1]] enz = total_size - eno;
    
    // expected number of false ones
    vector[vint4[1]] excess_ones = (eno ./ to_vector(vreal1[1 : vint4[1]]))
                                   - eno;
    
    // likelihoods given true zero for 0 and 1
    vector[vint4[1]] vr1 = excess_ones ./ enz; // population prop of true 0s that are misclassified
    vector[vint4[1]] vr0 = 1 - vr1; // population prop of true 0s that are not misclassified
    
    int N = size(mu);
    array[N] real vr;
    for (i in 1 : N) {
      if (fp[i] == 0) {
        vr[i] = vr0[vint5[i]];
      } else {
        vr[i] = vr1[vint5[i]] * vreal3[i];
      }
    }
    
    // Likelihood given a true one for 0 and 1
    array[N] real vor;
    for (i in 1 : N) {
      if (fp[i] == 0) {
        vor[i] = 0;
      } else {
        vor[i] = vreal2[i];
      }
    }
    
    // Initialize and compute log-likelihood
    real lp = 0;
    for (i in 1 : vint1[1]) {
      array[vint2[i]] int indices = index_array[i, 1 : vint2[i]];
      lp += log_sum_exp(bernoulli_logit_lpmf(1 | occ[i])
                        + emission_1_fp(to_row_vector(mu[indices]),
                                        vr[indices], vor[indices]),
                        bernoulli_logit_lpmf(0 | occ[i])
                        + emission_0_fp(vr[indices]));
    }
    return lp;
  }
}
data {
  int<lower=1> N; // total number of observations
  vector[N] Y; // response variable
  // data for custom real vectors
  array[N] real vreal1;
  // data for custom real vectors
  array[N] real vreal2;
  // data for custom real vectors
  array[N] real vreal3;
  // data for custom integer vectors
  array[N] int vint1;
  // data for custom integer vectors
  array[N] int vint2;
  // data for custom integer vectors
  array[N] int vint3;
  // data for custom integer vectors
  array[N] int vint4;
  // data for custom integer vectors
  array[N] int vint5;
  // data for custom integer vectors
  array[N] int vint6;
  // data for custom integer vectors
  array[N] int vint7;
  // data for custom integer vectors
  array[N] int vint8;
  // data for custom integer vectors
  array[N] int vint9;
  int<lower=1> K; // number of population-level effects
  matrix[N, K] X; // population-level design matrix
  int<lower=1> K_occ; // number of population-level effects
  matrix[N, K_occ] X_occ; // population-level design matrix
  int prior_only; // should the likelihood be ignored?
}
transformed data {
  
}
parameters {
  vector[K] b; // population-level effects
  vector[K_occ] b_occ; // population-level effects
}
transformed parameters {
  real lprior = 0; // prior contributions to the log posterior
}
model {
  // likelihood including constants
  if (!prior_only) {
    // initialize linear predictor term
    vector[N] mu = rep_vector(0.0, N);
    // initialize linear predictor term
    vector[N] occ = rep_vector(0.0, N);
    mu += X * b;
    occ += X_occ * b_occ;
    target += occupancy_single_fp_lpdf(Y | mu, occ, vint1, vint2, vint3, vint4, vint5, vint6, vint7, vint8, vint9, vreal1, vreal2, vreal3);
  }
  // priors including constants
  target += lprior;
}
generated quantities {
  
}