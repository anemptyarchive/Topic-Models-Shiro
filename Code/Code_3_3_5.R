
# Variational Bayes for LDA (1) -------------------------------------------

# 利用パッケージ
library(tidyverse)


# パラメータの設定 ----------------------------------------------------------------


# 事前分布のパラメータ
alpha_k <- rep(2, K)
beta_v  <- rep(2, V)


# トピック分布の初期値
theta_dk <- seq(0, 1, by = 0.01) %>% 
            sample(size = M * K, replace = TRUE) %>% 
            matrix(nrow = M, ncol = K)
# 正規化
theta_dk <- theta_dk / apply(theta_dk, 1, sum)


# 単語分布の初期値
phi_kv <- seq(0, 1, by = 0.01) %>% 
          sample(size = K * V, replace = TRUE) %>% 
          matrix(nrow = K, ncol = V)
# 正規化
phi_kv <- phi_kv / apply(phi_kv, 1, sum)


# イタレーション数
Iter <- 10


# 変分ベイズ -------------------------------------------------------------------


for(I in 1:Iter) { # 試行回数
  
  for(d in 1:M) { # 各文書
    
    for(v in 1:V) { # 各語彙
      if(n_dv[d, v] > 0) {
        for(n in 1:n_dv[d, v]) { # 各単語
          
          # 潜在トピック集合の事後分布:式(3.99)
          term1 <- digamma(eta_kv[, v]) - digamma(apply(eta_kv, 1, sum))
          term2 <- digamma(eta_dk[d, ]) - digamma(apply(eta_dk, 2, sum))
          q_z <- exp(term1) * exp(term2)
        }
      }
    }
    
    # 期待値
    n_dk
    n_kv
    
    
    # 事後分布のパラメータ:式(3.89)
    eta_dk[d, ] <- n_dk[d, ] + alpha_k

    # トピック分布の事後分布を計算:式(3.90)
    theta_dk <- theta_dk[d, ] ^ (eta_dk[d, ] - 1)
    
    for(k in 1:K) { # 各トピック
      
      # 事後分布のパラメータ:式(3.95)
      eta_kv[k, ] <- n_kv[k, ] + beta_v
      
      # 単語分布の事後分布を計算:式(3.96)
      phi_kv[k, ] <- phi_kv[k, ] ^ (eta_kv[k, ] - 1)
    }
  }
  
  # パラメータを正規化
  theta_dk <- theta_dk / apply(theta_dk, 1, sum)
  phi_kv   <- phi_kv / apply(phi_kv, 1, sum)
}






