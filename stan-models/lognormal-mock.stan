// Actual compositions known
// Diagonal covariance matrix
// TODO: compare types (matrix, vector, row_vector to *-cov.stan
// Uses the "Soft Centering" approach for beta described in section 1.7 of the
// Stan Users Guide version 2.19
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
    matrix<lower=0>[N, K] actual;
    // Scale for prior on spread of mean log efficiencies
    real<lower=0> scale_sigma_beta;
    // Scale for prior on lognormal distribution spread param
    real<lower=0> scale_Sigma;
}
transformed data {
    // TODO: consider turning actual and log_actual into vector arrays; may
    // improve speed
    matrix[N, K] log_actual;
    for (i in 1:N) {
        log_actual[i] = log(actual[i]);
    }
}
parameters {
    // Realized log efficiencies in each sample
    matrix[N,K] log_B;
    // Mean log efficiencies of taxa 1:K relative to all species
    row_vector[K] beta;
    // Std. dev of realized log efficiencies
    vector<lower=0>[K] sigma_log_B;
    // Spread of beta_i's
    real<lower = 0> sigma_beta; 
}
transformed parameters {
    // Covariance matrix for log efficiencies
    cov_matrix[K] Sigma;
    Sigma = diag_matrix(sigma_log_B .* sigma_log_B);  
}
model {
    // Prior on biases
    // TODO: consider allowing fatter tails
    beta ~ normal(0, sigma_beta);
    sigma_beta ~ exponential(scale_sigma_beta);
    // Prior on efficiency variance
    sigma_log_B ~ exponential(scale_Sigma);
    // Observed counts after random multiplicative error and multinomial
    // sampling
    for (i in 1:N) {
        log_B[i] ~ multi_normal(beta, Sigma);
        observed[i] ~ multinomial( 
            softmax( log_actual[i]' + log_B[i]' )
        );
    } 
}
generated quantities {
    row_vector[K] mean_clr_B;
    mean_clr_B = beta - mean(beta);
}
