
# Stochastic Variational Bayes for LDA ------------------------------------

# 利用パッケージ
library(tidyverse)

## 簡易文書データ
# 文書数
M <- 10
# 語彙数
V <- 20
# 文書ごとの各語彙数
n_dv <- matrix(sample(1:5, M * V, replace = TRUE), M, V)
# 文書ごとの単語数
n_d <- apply(n_dv, 1, sum)


# パラメータの設定 ----------------------------------------------------------------

# サンプリング数を指定
S <- 1000

# 試行回数を指定
InIter <- 25

# ステップサイズを指定
nu <- 1

# トピック数を指定
K <- 5

# 事前分布のパラメータを指定
alpha_k <- rep(2, K)
beta_v <- rep(2, V)

# 事後分布のパラメータの初期値
xi_dk.theta <- xi_dk.theta <- rep(n_d / K, K) + rep(alpha_k, each = M) %>% # E_n_dk + alpha_dk
  matrix(nrow = M, ncol = K)
xi_kv_s.phi <- xi_kv.phi <- seq(1, 5) %>% # とり得る値の範囲を指定
  sample(size = K * V, replace = TRUE) %>% 
  matrix(nrow = K, ncol = V)

# 潜在トピック集合の分布
z_dv_k <- array(0, dim = c(M, V, K))
for(d in 1:M) {
  for(v in 1:V) {
    if(n_dv[d, v] > 0) {
      
      # ランダムに値を生成
      tmp_q_z <- sample(seq(0, 1, by = 0.01), size = K, replace = TRUE)
      
      # 正規化
      z_dv_k[d, v, ] <- tmp_q_z / sum(tmp_q_z)
    }
  }
}

# カウントの期待値用受け皿
tmp_n <- array(0, dim = c(M, V, K))

# 文書dの語彙vにおいてトピックkが割り当てられた単語数の期待値の受け皿
E_n_dkv <- array(0, dim = c(M, K, V))


# 確率的変分ベイズ推定 --------------------------------------------------------------

# 推移の確認用
trace_xi_theta <- array(0, dim = c(M, K, S + 1))
trace_xi_phi <- array(0, dim = c(K, V, S + 1))
# 初期値を代入
trace_xi_theta[, , 1] <- xi_dk.theta
trace_xi_phi[, , 1] <- xi_kv_s.phi

for(s in 1:S) { ## (サンプリング)
  
  # 動作確認
  start_time <- Sys.time()
  
  # 文書をサンプリング
  d <- sample(1:M, size = 1)
  
  for(II in 1:InIter) { ## (イタレーション)
    
    for(v in 1:V) { ## (各語彙)
      if(n_dv[d, v] > 0) {
        
        # 潜在トピック集合の近似事後分布を計算:式(3.99)
        term1 <- exp(digamma(xi_kv_s.phi[, v])) / exp(digamma(apply(xi_kv_s.phi, 1, sum)))
        term2 <- exp(digamma(xi_dk.theta[d, ])) / exp(digamma(apply(xi_dk.theta, 2, sum)))
        tmp_q_z <- term1 * term2
        
        # 正規化
        z_dv_k[d, v, ] <- tmp_q_z / sum(tmp_q_z)
        
      }
    } ## (/各語彙)
    
    ## カウントの期待値
    tmp_n[d, , ] <- z_dv_k[d, , ] * n_dv[d, ]
    
    # 文書dにおいてトピックkが割り当てられた単語数の期待値を計算
    E_n_dk <- apply(tmp_n[d, , ], 2, sum) ## (ベクトルで扱うことに注意)
    
    # トピック分布の近似事後分布のパラメータを計算:式(3.89)
    xi_dk.theta[d, ] <- E_n_dk + alpha_k
    
  } ## (/イタレーション)
  
  for(k in 1:K) { ## (各トピック)
    
    # 文書dの語彙vにおいてトピックkが割り当てられた単語数の期待値を計算
    E_n_dkv[d, k, ] <- tmp_n[d, , k]
    
    # 単語分布の近似事後分布のパラメータを計算:式(3.159)
    xi_kv_s.phi[k, ] <- xi_kv_s.phi[k, ] + nu * (M * E_n_dkv[d, k, ] + beta_v - xi_kv_s.phi[k, ])
    
  } ## (/各トピック)
  
  # 
  xi_kv.phi <- xi_kv.phi + xi_kv_s.phi
  
  # 推移の確認用配列に代入
  trace_xi_theta[, , s + 1] <- xi_dk.theta
  trace_xi_phi[, , s + 1] <- xi_kv_s.phi
  
  # 動作確認
  print(paste0("s=", s, "(", round(s / S * 100, 1), "%)", "...", round(Sys.time() - start_time, 2)))

} ## (/サンプリング)

# 
xi_kv.phi <- xi_kv.phi / S


# 推定結果の確認 -----------------------------------------------------------------

## トピック分布(期待値)
# thetaの期待値を計算:式(2.10)
theta_dk <- xi_dk.theta / apply(xi_dk.theta, 1, sum)

# 作図用のデータフレームを作成
theta_WideDF <- cbind(
  as.data.frame(theta_dk), 
  doc = as.factor(1:M) # 文書番号
)

# データフレームをlong型に変換
theta_LongDF <- pivot_longer(
  theta_WideDF, 
  cols = -doc, 
  names_to = "topic", 
  names_prefix = "V", 
  names_ptypes = list(topic = factor()), 
  values_to = "prob"
)

# 作図
ggplot(theta_LongDF, aes(x = topic, y = prob, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  facet_wrap(~ doc, labeller = label_both) + # グラフの分割
  labs(title = "Stochastic Variational Bayes for LDA", 
       subtitle = expression(xi^theta)) # ラベル


## 単語分布(期待値)
# phiの期待値を計算:式(2.10)
phi_kv <- xi_kv.phi / apply(xi_kv.phi, 1, sum)

# 作図用のデータフレームを作成
phi_WideDF <- cbind(
  as.data.frame(phi_kv), 
  topic = as.factor(1:K) # トピック番号
)

# データフレームをlong型に変換
phi_LongDF <- pivot_longer(
  phi_WideDF, 
  cols = -topic, 
  names_to = "word", 
  names_prefix = "V", 
  names_ptypes = list(word = factor()), 
  values_to = "prob"
)

# 作図
ggplot(phi_LongDF, aes(x = word, y = prob, fill = word, color = word)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  facet_wrap(~ topic, labeller = label_both) + # グラフの分割
  theme(legend.position = "none") + # 凡例
  scale_x_discrete(breaks = seq(1, V, by = 10)) + # x軸目盛(連続値)
  labs(title = "Stochastic Variational Bayes for LDA", 
       subtitle = expression(xi^phi)) # ラベル


# 推移の確認 -------------------------------------------------------------------

# 利用パッケージ
library(gganimate)


## 作図用データフレームを作成
# 配列をデータフレームに変換
trace_xi_theta_WideDF <- data.frame()
trace_xi_phi_WideDF <- data.frame()
for(s in 1:(S + 1)) {
  # データフレームを作成
  tmp_theta_df <- cbind(
    as_tibble(trace_xi_theta[, , s]), 
    doc = as.factor(1:M), # 文書番号
    S = s - 1 # 試行回数
  )
  tmp_trace_xi_phi_df <- cbind(
    as_tibble(trace_xi_phi[, , s]), 
    topic = as.factor(1:K), # トピック番号
    S = s - 1 # 試行回数
  )
  # 結合
  trace_xi_theta_WideDF <- rbind(trace_xi_theta_WideDF, tmp_theta_df)
  trace_xi_phi_WideDF <- rbind(trace_xi_phi_WideDF, tmp_trace_xi_phi_df)
}

# データフレームをlong型に変換
trace_xi_theta_LongDF <- pivot_longer(
  trace_xi_theta_WideDF, 
  cols = -c(doc, S), # 変換しない列
  names_to = "topic", # 現列名を格納する列の名前
  names_prefix = "V", # 現列名から取り除く文字列
  names_ptypes = list(topic = factor()), # 現列名を格納する際のデータ型
  values_to = "value" # 現セルを格納する列の名前
)
trace_xi_phi_LongDF <- pivot_longer(
  trace_xi_phi_WideDF, 
  cols = -c(topic, S), # 変換しない列
  names_to = "word", # 現列名を格納する列の名前
  names_prefix = "V", # 現列名から取り除く文字列
  names_ptypes = list(word = factor()), # 現列名を格納する際のデータ型
  values_to = "value" # 現セルを格納する列の名前
)


### 棒グラフ(gif)で可視化
## トピック分布のパラメータ
# 作図
graph_xi_theta <- ggplot(trace_xi_theta_LongDF, aes(topic, value, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  facet_wrap(~ doc, labeller = label_both) + # グラフの分割
  transition_manual(S) + # フレーム
  labs(title = "Stochastic Variational Bayes for LDA", 
       subtitle = "S={current_frame}") # ラベル

# 描画
animate(graph_xi_theta, nframes = S + 1, fps = 20)


## 単語分布のパラメータ
# 作図
graph_xi_phi <- ggplot(trace_xi_phi_LongDF, aes(word, value, fill = word, color = word)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  facet_wrap(~ topic, labeller = label_both) + # グラフの分割
  theme(legend.position = "none") + # 凡例
  scale_x_discrete(breaks = seq(1, V, by = 10)) + # x軸目盛(連続値)
  transition_manual(S) + # フレーム
  labs(title = "Stochastic Variational Bayes for LDA", 
       subtitle = "S={current_frame}") # ラベル

# 描画
animate(graph_xi_phi, nframes = S + 1, fps = 20)


### 折れ線グラフで可視化
## トピック分布のパラメータ
# 文書番号を指定
doc_num <- 1

# 作図
trace_xi_theta_LongDF %>% 
  filter(doc == doc_num) %>% 
  ggplot(aes(x = S, y = value, color = topic)) + 
    geom_line() + # 折れ線グラフ
    labs(title = "Stochastic Variational Bayes for LDA", 
         subtitle = paste0("d=", doc_num)) # ラベル


## 単語分布のパラメータ
# トピック番号を指定
topic_num <- 1

# 作図
trace_xi_phi_LongDF %>% 
  filter(topic == topic_num) %>% 
  ggplot(aes(x = S, y = value, color = word)) + 
    geom_line(alpha = 0.5) + # 折れ線グラフ
    theme(legend.position = "none") + # 凡例
    labs(title = "Stochastic Variational Bayes for LDA", 
         subtitle = paste0("k=", topic_num)) # ラベル


