
# Collapsed Gibbs sampler for LDA -----------------------------------------

# 利用パッケージ
library(tidyverse)


# パラメータの設定 -----------------------------------------------------------------


# トピック数
K <- 4


# ハイパーパラメータ
alpha_k <- rep(2, K)
beta_v  <- rep(2, V)


# 潜在トピック集合の初期値
z_di <- array(0, dim = c(M, V, max(n_dv)))

# 各文書において各トピックが割り当てられた単語数
n_dk <- matrix(0, nrow = M, ncol = K)

# 全文書において各トピックが割り当てられた単語数
n_kv <- matrix(0, nrow = K, ncol = V)


# サンプリング回数
S <- 2


# ギブスサンプリング ----------------------------------------------------------------------


for(s in 1:S) { # サンプリング回数
  
  # カウントを初期化
  new_n_kv <- matrix(0, nrow = K, ncol = V)
  new_n_dk <- matrix(0, nrow = M, ncol = K)
  
  for(d in 1:M) { # 各文書
    
    for(v in 1:V) { # 各語彙
      if(n_dv[d, v] > 0) {
        for(n in 1:n_dv[d, v]) { # 各単語
          
          # トピック番号を取り出す
          k <- z_di[d, v, n]
          
          # d,i要素のみのカウント
          count_kv <- matrix(0, nrow = K, ncol = V)
          count_dk <- matrix(0, nrow = M, ncol = K)
          count_kv[k, v] <- 1
          count_dk[d, k] <- 1
          
          # d,i要素を除いたカウント
          n_kv.di <- n_kv - count_kv
          n_dk.di <- n_dk - count_dk
          
          
          # サンプリング確率を計算
          term1 <- (n_kv.di[, v] + beta_v[v]) / apply(t(n_kv.di) + beta_v, 2, sum)
          term2 <- (n_dk.di[d, ] + alpha_k) / sum(n_dk.di[d, ] + alpha_k)
          p_z <- term1 * term2
          
          
          # 潜在トピックを割り当て
          res_z <- rmultinom(n = 1, size = 1, prob = p_z)
          z_di[d, v, n] <- which(res_z == 1)
        }
      }
    }
  }
  
  # カウントを更新
  for(k in 1:K) {
    n_kv[k, ] <- apply(z_di == k, 2, sum)
    n_dk[, k] <- apply(z_di == k, 1, sum)
  }
}

sum(apply(n_dk, 1, sum) == n_d) == M
sum(apply(n_kv, 2, sum) == apply(n_dv, 2, sum)) == V
sum(n_dk) == sum(n_kv)


