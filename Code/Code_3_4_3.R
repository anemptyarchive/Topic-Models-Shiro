
# Stochastic Variational Bayes for LDA ------------------------------------

# 利用パッケージ
library(tidyverse)



# パラメータの設定 ----------------------------------------------------------------

# サンプリング数
S <- 10

# 試行回数
InIter <- 3

# ステップサイズ
nu <- 1

# トピック数
K <- 5

# 事前分布のパラメータ
alpha_k <- rep(2, K)
beta_v <- rep(2, V)

# 事後分布のパラメータの初期値
xi_dk <- seq(1, 5) %>% 
  sample(size = M * K, replace = TRUE) %>% 
  matrix(nrow = M, ncol = K)
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

# 推移の確認用データフレームを作成
trace_xi_theta <- array(0, dim = c(M, K, S + 1))
trace_xi_theta[, , 1] <- xi_dk
trace_xi_phi <- array(0, dim = c(K, V, S + 1))
trace_xi_phi[, , 1] <- xi_kv

for(d in 1:M) { ## (各文書)
  
  for(s in 1:S) {
    
    for(II in 1:InIter) {
      
      for(v in 1:V) {
        if(n_dv[d, v] > 0) {
          for(n in 1:n_dv[d, v]) {
            
            # 潜在トピック集合の近似事後分布を計算:式(3.99)
            term1 <- exp(digamma(xi_kv[, v])) / exp(digamma(apply(xi_kv, 1, sum)))
            term2 <- exp(digamma(xi_dk[d, ])) / exp(digamma(apply(xi_dk, 2, sum)))
            z_di_k[d, v, n, ] <- term1 * term2
            
          }
        }
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
  
  # 推移の確認用データフレームを作成
  trace_xi_theta[, , s + 1] <- xi_dk
  trace_xi_phi[, , s + 1] <- xi_kv
  
} ## (/各文書)



# 推定結果の確認 -----------------------------------------------------------------


# 
library(gganimate)

trace_xi_theta_WideDF <- cbind(
  as_tibble(trace_xi_theta[, , 1]), 
  doc = as.factor(1:M), 
  S = 0
)
trace_xi_phi_WideDF <- cbind(
  as_tibble(trace_xi_phi[, , 1]), 
  topic = as.factor(1:K), 
  S = 0
)
for(s in 1:S) {
  tmp_theta_df <- cbind(
    as_tibble(trace_xi_theta[, , s + 1]), 
    doc = as.factor(1:M), 
    S = s
  )
  trace_xi_theta_WideDF <- rbind(trace_xi_theta_WideDF, tmp_theta_df)
  
  tmp_trace_xi_phi_df <- cbind(
    as_tibble(trace_xi_phi[, , s + 1]), 
    topic = as.factor(1:K), 
    S = s
  )
  trace_xi_phi_WideDF <- rbind(trace_xi_phi_WideDF, tmp_trace_xi_phi_df)
}

# long型に変換
trace_xi_theta_LongDF <- pivot_longer(
  trace_xi_theta_WideDF, 
  cols = -c(doc, S), 
  names_to = "topic", 
  names_prefix = "V", 
  names_ptypes = list(topic = factor()), 
  values_to = "value"
)
trace_xi_phi_LongDF <- pivot_longer(
  trace_xi_phi_WideDF, 
  cols = -c(topic, S), 
  names_to = "word", 
  names_prefix = "V", 
  names_ptypes = list(word = factor()), 
  values_to = "value"
)

# 作図
graph_xi_theta <- ggplot(trace_xi_theta_LongDF, aes(topic, value, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~ doc, labeller = label_both) + 
  transition_manual(S) + 
  labs(title = "Stochastic Variational Bayes for LDA")

animate(graph_xi_theta, nframes = S + 1, fps = 10)

graph_xi_phi <- ggplot(trace_xi_phi_LongDF, aes(word, value, fill = word)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~ topic, labeller = label_both) + 
  theme(legend.position = "none") + 
  scale_x_discrete(breaks = seq(1, V, by = 10)) + 
  transition_manual(S) + 
  labs(title = "Stochastic Variational Bayes for LDA")

animate(graph_xi_phi, nframes = S + 1, fps = 10)
