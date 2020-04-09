
# Particle filter for LDA -------------------------------------------------

# 利用パッケージ
library(tidyverse)



# パラメータの設定 ----------------------------------------------------------------

# サンプリング数
S <- 10

# リサンプリング数
R <- 5

# しきい値
threshold <- 0.1


# トピック数
K <- 5

# 事前分布のパラメータ
alpha_k <- rep(2, K)
beta_v <- rep(2, V)

# 潜在トピック集合の事後分布
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
z_di_s <- array(0, dim = c(M, V, max(n_dv), S))

# 重みの受け皿
w_z_di_s <- array(0, dim = c(M, V, max(n_dv), S))

# 割り当てられたトピックに関する単語数の初期値
n_dk_s <- array(0, dim = c(M, K, S))
n_kv_s <- array(0, dim = c(K, V, S))


# パーティクルフィルタ --------------------------------------------------------------

d <- v <- n <- k <- s <- 1
for(d in 1:M) { ## (各文書)
  
  for(v in 1:V) { ## (各語彙)
    if(n_dv[d, v] > 0) {
      for(n in 1:n_dv[d, v]) { ## (各単語)
        
        for(s in 1:S) { ## (サンプリング)
          
          # 潜在トピック集合の分布を計算:式(3.172)
          term1 <- (n_kv_s[, v, s] + beta_v[v]) / apply(t(n_kv_s[, , s]) + beta_v, 2, sum)
          term2 <- (n_dk_s[d, , s] + alpha_k) / sum(n_dk_s[d, , s] + alpha_k)
          tmp_q_z <- term1 * term2
          z_di_k[d, v, n, ] <- tmp_q_z / sum(tmp_q_z)
          
          # 潜在トピックをサンプリング
          z_di_s[d, v, n, s] <- sample(1:K, size = 1, prob = z_di_k[d, v, n, ])
          
          # カウント
          tmp_n_k <- rep(0, K)
          tmp_n_k[z_di_s[d, v, n, s]] <- 1
          n_dk_s[d, , s] <- n_dk_s[d, , s] + tmp_n_k
          n_kv_s[, v, s] <- n_kv_s[, v, s] + tmp_n_k
          
          ## 重みを計算:式(3.175)
          # 式(3.176)の計算
          term1 <- (n_kv_s[, v, s] + beta_v[v]) / apply(t(n_kv_s[, , s]) + beta_v, 2, sum)
          term2 <- (n_dk_s[d, , s] + alpha_k) / sum(n_dk_s[d, , s] + alpha_k)
          p_wz_s <- term1 * term2
          
          # 重みを計算:式(3.175)
          w_z_di_s[d, v, n, s] <- w_z_di_s[d, v, n, s] * sum(p_wz_s)
          
        } ## (/サンプリング)
        
        # 重みを正規化
        w_z_di_s[d, v, n, ] <- w_z_di_s[d, v, n, ] / sum(w_z_di_s[d, v, n, ])
        # NaN対策(仮)
        if(is.nan(w_z_di_s[d, v, n, 1])) {
          w_z_di_s[d, v, n, ] <- rep(1 / K, K)
        }
        
        # ESSを計算
        ESS <- 1 / apply(w_z_di_s ^ 2, c(1, 2), sum)
        
        if(ESS < threshold) {
          
          for(r in 1:R) {
            
            l <- sample(1:M, size = 1)
            m_v <- sample(1:V, size = 1, prob = n_dv[l, ])
            m_n <- sample(1:n_dv[l, m_v], size = 1)
            
            for(s in 1:S) {
              
              # カウント
              n_dk.lm_s <- rep(0, K)
              n_kv.lm_s <- matrix(0, nrow = K, ncol = V)
              for(k in 1:K) {
                n_dk.lm_s[k] <- sum(z_di_s[1:l, 1:m_v, 1:m_n, s] == k)
                n_kv.lm_s[k, ] <- apply(z_di_s[1:l, , , s] == k, 2, sum)
              }
              tmp_n_k <- rep(0, K)
              tmp_n_k[z_di_s[l, m_v, m_n, s]] <- 1
              n_dk.lm_s <- n_dk.lm_s - tmp_n_k
              n_kv.lm_s <- n_kv.lm_s - tmp_n_k
              
              # リサンプリング確率を計算:式(3.179)
              term1 <- (n_kv.lm_s + beta_v[m_v]) / apply(t(n_kv.lm_s[, , s]) + beta_v, 2, sum)
              term2 <- (n_dk.lm_s + alpha_k) / sum(n_dk.lm_s + alpha_k)
              tmp_q_z <- term1 * term2
              
              # 潜在トピックをサンプリング
              z_di_s[l, m_v, m_n, s] <- sample(1:K, size = 1, prob = tmp_q_z)
            
            } ## (/)
          } ## (/)
          
          w_z_di_s <- array(1 / S, dim = c(M, V, max(n_dv), S))
        }
        
      } ## (/各単語)
    }
  } ## (/各語彙)
} ## (/各文書)



# 推定結果の確認 -----------------------------------------------------------------


