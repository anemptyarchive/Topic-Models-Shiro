
# Particle filter for LDA -------------------------------------------------

# 利用パッケージ
library(tidyverse)



# パラメータの設定 ----------------------------------------------------------------

# サンプリング数
S <- 10

# トピック数
K <- 5


# 潜在トピック集合の分布
z_di_k <- array(0, dim = c(M, V, max(n_dv), K))
for(d in 1:M) {
  for(v in 1:V) {
    if(n_dv[d, v] > 0) {
      for(n in 1:n_dv[d, v]) {
        # ランダムに値を生成
        tmp_z <- sample(seq(0, 1, by = 0.01), size = K, replace = TRUE)
        # 正規化
        z_di_k[d, v, n, ] <- tmp_z / sum(tmp_z)
      }
    }
  }
}

# パーティクルフィルタ --------------------------------------------------------------


for(d in 1:M) { ## (各文書)
  
  for(v in 1:V) { ## (各語彙)
    if(n_dv[d, v] > 0) {
      for(n in 1:n_dv[d, v]) { ## (各単語)
        
        for(s in 1:S) { ## (サンプリング)
          
          # 潜在トピック集合の分布を計算:式(3.172)
          term1 <- (n_kv.di[, v] + beta_v[v]) / apply(t(n_kv.di) + beta_v, 2, sum)
          term2 <- (n_dk.di[d, ] + alpha_k) / apply(t(n_dk.di) + alpha_k, 2, sum)
          z_di_k[d, v, n, ] <- term1 * term2
          
          # :式(3.176)
          term1 <- (n_kv.di[, v] + beta_v[v]) / apply(t(n_kv.di) + beta_v, 2, sum)
          term2 <- (n_dk.di[d, ] + alpha_k) / apply(t(n_dk.di) + alpha_k, 2, sum)
          p_wz <- term1 * term2
          
          # 重みを計算:式(3.175)
          w_z_di.s[d, v, n, ] <- w_z_di.s[d, v, n, ] * sum(p_wz)
          
        } ## (/サンプリング)
        
      } ## (/各単語)
    }
  } ## (/各語彙)
} ## (/各文書)

