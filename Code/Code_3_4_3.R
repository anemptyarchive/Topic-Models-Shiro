
# Stochastic Variational Bayes for LDA ------------------------------------

# 利用パッケージ
library(tidyverse)



# パラメータの設定 ----------------------------------------------------------------


# サンプリング数
S <- 10

# ステップサイズ
nu <- 1

# トピック数
K <- 5

# 事前分布のパラメータ
alpha_k <- rep(2, K)
beta_v <- rep(2, V)

# 事後分布のパラメータの初期値
xi_kv <- seq(1, 5) %>% 
  sample(size = K * V, replace = TRUE) %>% 
  matrix(nrow = K, ncol = V)

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

# 文書dにおいてトピックkが割り当てられた単語数の期待値
E_n_dk <- apply(z_di_k, c(1, 4), sum)

# 文書dの語彙vにおいてトピックkが割り当てられた単語数の期待値
E_n_dkv <- apply(z_di_k, c(1, 4, 2), sum)

# 動作確認
sum(E_n_dk) == sum(n_dv)
sum(E_n_dkv) == sum(n_dv)


# 確率的変分ベイズ推定 --------------------------------------------------------------

for(d in 1:M) { ## (各文書)
  for(s in 1:S) {
    
    for(II in 1:InIter) {
      
      for(i in 1:n_d) {
        
        # 潜在トピック集合の近似事後分布を計算:式(3.99)
        term1 <- exp(digamma(xi_kv)) / exp(digamma(apply(xi_kv, 1, sum)))
        term2 <- exp(digamma(xi_dk)) / exp(digamma(apply(xi_dk, 1, sum)))
        z_di_k[d, v, i, ] <- term1 * term2
        
      }
      
      # 文書dにおいてトピックkが割り当てられた単語数の期待値を計算
      E_n_dk[d, ] <- apply(z_di_k[d, , , ], 3, sum)
      
      # トピック分布の近似事後分布のパラメータを計算:式(3.89)
      xi_dk[d, ] <- E_n_dk[d, ] + alpha_k
      
    }
    
    # 文書dの語彙vにおいてトピックkが割り当てられた単語数の期待値を計算
    E_n_dkv <- apply(z_di_k, c(1, 4, 2), sum)
    
    for(k in 1:K) { ## (各トピック)
      
      # 文書番号dをサンプリング
      d <- sample(1:M, size = 1)
      
      # 単語分布の近似事後分布のパラメータを計算:式(3.159)
      xi_kv[k, ] <- xi_kv[k, ] + nu * (M * E_n_dkv[d, k, ] + beta_v - xi_kv[k, ])
      
    } ## (/各トピック)
  }
} ## (/各文書)
