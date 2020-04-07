
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
        tmp_q_z <- sample(seq(0, 1, by = 0.01), size = K, replace = TRUE)
        # 正規化
        z_di_k[d, v, n, ] <- tmp_q_z / sum(tmp_q_z)
      }
    }
  }
}

# 潜在トピック集合の受け皿
new_z_di.s <- array(0, dim = c(M, V, max(n_dv), s))

# 重みの受け皿
w_z_di.s <- array(0, dim = c(M, V, max(n_dv), S))



# パーティクルフィルタ --------------------------------------------------------------


for(d in 1:M) { ## (各文書)
  
  for(v in 1:V) { ## (各語彙)
    if(n_dv[d, v] > 0) {
      for(n in 1:n_dv[d, v]) { ## (各単語)
        
        for(s in 1:S) { ## (サンプリング)
          
          # 潜在トピック集合の分布を計算:式(3.172)
          term1 <- (n_kv.di[, v] + beta_v[v]) / apply(t(n_kv.di) + beta_v, 2, sum)
          term2 <- (n_dk.di[d, ] + alpha_k) / apply(t(n_dk.di) + alpha_k, 2, sum)
          tmp_q_z <- term1 * term2
          z_di_k[d, v, n, ] <- tmp_q_z / sum(tmp_q_z)
          
          # サンプリング
          new_z_di.s[d, v, n, s] <- sample(1:K, size = 1, prob = z_di_k[d, v, n, ])
          
          # カウント
          tmp_n_k <- rep(0, K)
          tmp_n_k[new_z_di.s[d, v, n, s]] <- 1
          n_dk.di.s[d, , s] <- tmp_n_k
          n_kv.di.s[, v, s] <- tmp_n_k
          
          ## 重みを更新
          # 式(3.176)の計算
          term1 <- (n_kv.di.s[, v, s] + beta_v[v]) / apply(t(n_kv.di.s[, , s]) + beta_v, 2, sum)
          term2 <- (n_dk.di.s[d, , s] + alpha_k) / apply(t(n_dk.di.s[, , s]) + alpha_k, 2, sum)
          p_wz.s <- term1 * term2
          
          # 重みを計算:式(3.175)
          w_z_di.s[d, v, n, s] <- w_z_di.s[d, v, n, s] * sum(p_wz.s)
          
        } ## (/サンプリング)
        
        # 重みを正規化
        w_z_di.s[d, v, n, ] <- w_z_di.s[d, v, n, ] / sum(w_z_di.s[d, v, n, ])
        
      } ## (/各単語)
    }
  } ## (/各語彙)
} ## (/各文書)



# 推定結果の確認 -----------------------------------------------------------------


