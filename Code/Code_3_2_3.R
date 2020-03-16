
# Gibbs sampler for LDA ---------------------------------------------------


# 利用パッケージ
library(tidyverse)


# パラメータ設定 -----------------------------------------------------------------


# トピック数
K <- 4


# ハイパーパラメータ
alpha_k <- vector("double", length = K)
beta_v <- vector("double", length = V)


# 単語分布の初期値
phi_kv <- matrix(sample(seq(0, 1, by = 0.01), size = K * V, replace = TRUE), nrow = K, ncol = V)
phi_kv <- phi_kv / apply(phi_kv, 1, sum)

# トピック分布の初期値
theta_dk <- matrix(sample(seq(0, 1, by = 0.01), size = M * K, replace = TRUE), nrow = M, ncol = K)
theta_dk <- theta_dk / apply(theta_dk, 1, sum)

# 潜在トピック集合
z_dv <- matrix(0, nrow = M, ncol = V)


# サンプリング回数
S <- 1000

p <- vector("double", length = K)
for(s in 1:S) {
  
  for(d in 1:M) { # 各文書
    
    for(v in 1:V) { # 各語彙
      
      # サンプリング確率
      p <- phi_kv[, v] * theta_dk[d,] / sum(phi_kv[, v] * theta_dk[d, ])
    }
    
    # トピック分布の更新
    C_theta <- lgamma(sum(n_dk[d, ] + alpha_k)) - sum(lgamma(n_dk[d, ] + alpha_k))
    theta_dk[d, ] <- exp(C_theta + (n_dk + alpha_k - 1) * theta_dk[d, ])
  }
  
  for(k in 1:K) { # 各トピック
    
    # 単語分布の更新
    C_phi <- lgamma(sum(n_kv[k, ] + beta_v)) - sum(lgamma(n_kv[k, ] + beta_v))
    phi_kv[k, ] <- exp(C_phi + (n_kv[k, ] + beta_v - 1) * phi_kv[k, ])
  }
}

