// Actual compositions known
// Diagonal covariance matrix
// Uses the QR approach for beta described in section 1.7 of the Stan Users
// Guide version 2.19
functions {
}
data {
    // Number of taxa 
    int<lower=0> K;
    // Number of samples 
    int<lower=0> N;
    // Count data
    int<lower=0> observed[N, K];
    // Actual relative abundances
    vector<lower=0>[K] actual[N];
    // Scale for prior on spread of mean log efficiencies
    real<lower=0> scale_sigma_beta;
    // Scale for prior on lognormal distribution spread param
    real<lower=0> scale_Sigma;
}
transformed data {
    vector[K] log_actual[N];
    matrix[K, K] M = diag_matrix(rep_vector(1,K));
    matrix[K, K-1] M_qr;
    for (i in 1:N) {
        log_actual[i] = log(actual[i]);
    }
    for (i in 1:K-1) M[K,i] = -1;
    M[K,K] = 0;
    M_qr = qr_Q(M)[ , 1:(K-1)];
}
parameters {
    // Realized log efficiencies in each sample
    vector[K] log_B[N];
    // Mean log efficiencies of taxa 1:K-1 relative to taxon K
    vector[K-1] beta_raw;
    // Std. dev of realized log efficiencies
    vector<lower=0>[K] sigma_log_B;
    // Spread of beta_i's
    real<lower = 0> sigma_beta; 
}
transformed parameters {
    // Mean log efficiencies of taxa 1:K relative to all species
    vector[K] beta =  M_qr * beta_raw;
    // Covariance matrix for log efficiencies
    cov_matrix[K] Sigma;
    Sigma = diag_matrix(sigma_log_B .* sigma_log_B);
}
model {
    // Prior on log efficiencies; see https://mc-stan.org/docs/2_19/stan-users-guide/parameterizing-centered-vectors.html
    beta ~ normal(0, sigma_beta * inv(sqrt(1 - inv(K))));
    sigma_beta ~ exponential(scale_sigma_beta);
    // Prior on efficiency variance
    sigma_log_B ~ exponential(scale_Sigma);
    // Observed counts after random multiplicative error and multinomial
    // sampling
    for (i in 1:N) {
        log_B[i] ~ multi_normal(beta, Sigma);
        observed[i] ~ multinomial( 
            softmax( log_actual[i] + log_B[i] )
        );
    } 
}
generated quantities {
    vector[K] mean_clr_B;
    mean_clr_B = beta - mean(beta);
}
