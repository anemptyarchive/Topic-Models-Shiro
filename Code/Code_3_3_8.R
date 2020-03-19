
# Collapsed Variational Bayes Method for LDA -------------------------------------------

# 利用パッケージ
library(tidyverse)


# パラメータの設定 ----------------------------------------------------------------


# トピック数
K <- 5


# カウント
z_dvk.di <- z_dvk
z_dvk.di[d, v, ] <- 0

# E[n_{d,k}^{\d,i}]
E_n_dk.di <- apply(z_dvk.di, c(1, 3), sum)

# V[n_{d,k}^{\d,i}]
V_n_dk.di <- apply(z_dvk.di * (1 - z_dvk.di), c(1, 3), sum)

tmp_z_dvk.di <- array(0, c(M, V, K))
for(k in 1:K) {
  tmp_z_dvk.di[, , k] <- z_dvk.di[, , k] * n_dv
}
# E[n_{k,v}^{\d,i}]
E_n_kv.di <- apply(tmp_z_dvk.di, c(3, 2), sum)

# V[n_{k,v}^{\d,i}]
V_n_kv.di <- apply(tmp_z_dvk.di * (1 - tmp_z_dvk.di), c(3, 2), sum)


# イタレーション数
Iter <- 10


# 周辺化変分ベイズ ----------------------------------------------------------------


for(I in 1:Iter) {
  
  for(d in 1:M) { # 各文書
    
    for(v in 1:V) { # 各語彙
      if(n_dv[d, v] > 0) {
        for(n in 1:n_dv[d, v]) { # 各単語
          
          # 潜在トピック集合の分布を更新:式(3.130)
          term1 <- (E_n_kv.di[, v] + beta_v[v]) / apply(t(E_n_kv.di) + beta_v, 2, sum) * (E_n_dk.di[d, ] + alpha_k)
          term2.1 <- V_n_kv.di[, v] / (2 * (E_n_kv.di[, v] + beta_v[v])^2)
          term2.2 <- V_n_dk.di[d, ] / (2 * (E_n_dk.di[d, ] + alpha_k)^2)
          term3 <- apply(V_n_kv.di, 1, sum) / (2 * apply(t(E_n_kv.di) + beta_v, 2, sum)^2)
          z_dvk[d, v, ] <- term1 * exp(term2.1 - term2.2) * exp(term3)
          
          for(k in 1:K) { # 各トピック
            
            # 
            
          }
        }
      }
    }
  }
}



