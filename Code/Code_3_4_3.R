
# Stochastic Variational Bayes for LDA ------------------------------------

# 利用パッケージ
library(tidyverse)

## 簡易文書データ
# 文書数
M <- 10
# 語彙数
V <- 20
# 文書ごとの各語彙数
n_dv <- matrix(sample(0:10, M * V, replace = TRUE), M, V)
# 文書ごとの単語数
n_d <- apply(n_dv, 1, sum)


# パラメータの設定 ----------------------------------------------------------------

# サンプリング数を指定
S <- 1000

# イタレーション数を指定
InIter <- 25

# ステップサイズ(学習率)を指定
nu <- 0.01

# トピック数を指定
K <- 5

# 事前分布のパラメータを指定
alpha_k <- rep(2, K)
beta_v  <- rep(2, V)

# 事後分布のパラメータを指定
xi_theta_dk <- rep(n_d / K, K) + rep(alpha_k, each = M) %>% # E_n_dk + alpha_dk
  matrix(nrow = M, ncol = K)
xi_phi_kv <- seq(1, 5) %>% # とり得る値の範囲を指定
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


# 確率的変分ベイズ推定 --------------------------------------------------------------

# 文書dの語彙vにおいてトピックkが割り当てられた単語数の期待値の受け皿
E_n_dkv <- matrix(0, nrow = K, ncol = V) ## (マトリクスで扱うことに注意)

# 推移の確認用
trace_xi_theta <- array(0, dim = c(M, K, S + 1))
trace_xi_phi   <- array(0, dim = c(K, V, S + 1))
# 初期値を代入
trace_xi_theta[, , 1] <- xi_theta_dk
trace_xi_phi[, , 1]   <- xi_phi_kv

for(s in 1:S) { ## (サンプリング)
  
  # 動作確認用に開始時間を記録
  start_time <- Sys.time()
  
  # 文書(番号)をサンプリング
  d <- sample(1:M, size = 1)
  
  for(II in 1:InIter) { ## (内側ループ)
    
    for(v in 1:V) { ## (各語彙)
      if(n_dv[d, v] > 0) {
        
        # 潜在トピック集合の近似事後分布を計算:式(3.99)
        term1 <- digamma(xi_phi_kv[, v]) - digamma(apply(xi_phi_kv, 1, sum))
        term2 <- digamma(xi_theta_dk[d, ]) - digamma(sum(xi_theta_dk[d, ]))
        tmp_q_z <- exp(term1 + term2)
        z_dv_k[d, v, ] <- tmp_q_z / sum(tmp_q_z) # 正規化
        
      }
    } ## (/各語彙)
    
    # 文書dにおいてトピックkが割り当てられた単語数の期待値を計算
    tmp_n_vk <- z_dv_k[d, , ] * n_dv[d, ]
    E_n_dk <- apply(z_dv_k[d, , ] * n_dv[d, ], 2, sum) ## (K次元ベクトルで扱うことに注意)
    
    # トピック分布の近似事後分布のパラメータを更新:式(3.89)
    xi_theta_dk[d, ] <- E_n_dk + alpha_k
    
  } ## (/内側ループ)
  
  for(k in 1:K) { ## (各トピック)
    
    # 文書dの語彙vにおいてトピックkが割り当てられた単語数の期待値を計算
    E_n_dkv[k, ] <- tmp_n_vk[, k]
    
    # 単語分布の近似事後分布のパラメータを更新:式(3.159)
    xi_phi_kv[k, ] <- xi_phi_kv[k, ] + nu * (M * E_n_dkv[k, ] + beta_v - xi_phi_kv[k, ])
    
  } ## (/各トピック)
  
  # 推移の確認用に値を保存
  trace_xi_theta[, , s + 1] <- xi_theta_dk
  trace_xi_phi[, , s + 1]   <- xi_phi_kv
  
  # 動作確認
  print(paste0("s=", s, "(", round(s / S * 100, 1), "%)", "...", round(Sys.time() - start_time, 3)))
}

# 処理の検証
tmp_E_n <- array(0, dim = c(M, K, V))
for(k in 1:K) {
  tmp_E_n[, k, ] <- z_dv_k[, , k] * n_dv
}
sum(near(apply(tmp_E_n, 1, sum), n_d)) == M
sum(near(apply(tmp_E_n, 3, sum), apply(n_dv, 2, sum))) == V


# 推定結果の確認 -----------------------------------------------------------------

## トピック分布(期待値)
# thetaの期待値を計算:式(2.10)
theta_dk <- xi_theta_dk / apply(xi_theta_dk, 1, sum)

# 作図用のデータフレームを作成
theta_df_wide <- cbind(
  as.data.frame(theta_dk), 
  doc = as.factor(1:M) # 文書番号
)

# データフレームをlong型に変換
theta_df_long <- pivot_longer(
  theta_df_wide, 
  cols = -doc, # 変換せずにそのまま残す現列名
  names_to = "topic", # 現列名を格納する新しい列の名前
  names_prefix = "V", # 現列名から取り除く文字列
  names_ptypes = list(topic = factor()), # 現列名を要素とする際の型
  values_to = "prob" # 現要素を格納する新しい列の名前
)

# 作図
ggplot(theta_df_long, aes(x = topic, y = prob, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  facet_wrap(~ doc, labeller = label_both) + # グラフの分割
  labs(title = "Stochastic Variational Bayes for LDA", 
       subtitle = expression(xi^theta)) # ラベル


## 単語分布(期待値)
# phiの期待値を計算:式(2.10)
phi_kv <- xi_phi_kv / apply(xi_phi_kv, 1, sum)

# 作図用のデータフレームを作成
phi_df_wide <- cbind(
  as.data.frame(phi_kv), 
  topic = as.factor(1:K)
)

# データフレームをlong型に変換
phi_df_long <- pivot_longer(
  phi_df_wide, 
  cols = -topic, # 変換せずにそのまま残す現列名
  names_to = "word", # 現列名を格納する新しい列の名前
  names_prefix = "V", # 現列名から取り除く文字列
  names_ptypes = list(word = factor()), # 現列名を要素とする際の型
  values_to = "prob" # 現要素を格納する新しい列の名前
)

# 作図
ggplot(phi_df_long, aes(x = word, y = prob, fill = word, color = word)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  facet_wrap(~ topic, labeller = label_both) + # グラフの分割
  theme(legend.position = "none") + # 凡例
  scale_x_discrete(breaks = seq(1, V, by = 10)) + # x軸目盛(連続値)
  labs(title = "Stochastic Variational Bayes for LDA", 
       subtitle = expression(xi^phi)) # ラベル


# 更新値の推移を確認:折れ線グラフ --------------------------------------------------------

## トピック分布のパラメータ
# 作図用データフレームを作成
trace_xi_theta_df_wide <- data.frame()
for(s in 1:(S + 1)) {
  # データフレームに変換
  tmp_theta_df <- cbind(
    as_tibble(trace_xi_theta[, , s]), 
    doc = as.factor(1:M), 
    sampling = s - 1
  )
  # 結合
  trace_xi_theta_df_wide <- rbind(trace_xi_theta_df_wide, tmp_theta_df)
}

# データフレームをlong型に変換
trace_xi_theta_df_long <- pivot_longer(
  trace_xi_theta_df_wide, 
  cols = -c(doc, sampling), # 変換しない列
  names_to = "topic", # 現列名を格納する列の名前
  names_prefix = "V", # 現列名から取り除く文字列
  names_ptypes = list(topic = factor()), # 現列名を格納する際のデータ型
  values_to = "value" # 現セルを格納する列の名前
)

# 文書番号を指定
doc_num <- 10

# 作図
trace_xi_theta_df_long %>% 
  filter(doc == doc_num) %>% 
  ggplot(aes(x = sampling, y = value, color = topic)) + 
    geom_line() + # 折れ線グラフ
    labs(title = "Stochastic Variational Bayes for LDA", 
         subtitle = paste0("d=", doc_num)) # ラベル


## 単語分布のパラメータ
# 作図用データフレームを作成
trace_xi_phi_df_wide <- data.frame()
for(s in 1:(S + 1)) {
  # データフレームに変換
  tmp_trace_xi_phi_df <- cbind(
    as_tibble(trace_xi_phi[, , s]), 
    topic = as.factor(1:K), 
    sampling = s - 1
  )
  # 結合
  trace_xi_phi_df_wide <- rbind(trace_xi_phi_df_wide, tmp_trace_xi_phi_df)
}

# データフレームをlong型に変換
trace_xi_phi_df_long <- pivot_longer(
  trace_xi_phi_df_wide, 
  cols = -c(topic, sampling), # 変換しない列
  names_to = "word", # 現列名を格納する列の名前
  names_prefix = "V", # 現列名から取り除く文字列
  names_ptypes = list(word = factor()), # 現列名を格納する際のデータ型
  values_to = "value" # 現セルを格納する列の名前
)

# トピック番号を指定
topic_num <- 1

# 作図
trace_xi_phi_df_long %>% 
  filter(topic == topic_num) %>% 
  ggplot(aes(x = sampling, y = value, color = word)) + 
    geom_line(alpha = 0.5) + # 折れ線グラフ
    theme(legend.position = "none") + # 凡例
    labs(title = "Stochastic Variational Bayes for LDA", 
         subtitle = paste0("k=", topic_num)) # ラベル


# 更新値の推移の確認:gif -------------------------------------------------------------------

# 利用パッケージ
library(gganimate)


## トピック分布のパラメータ
# 作図
graph_xi_theta <- ggplot(trace_xi_theta_df_long, aes(x = topic, y = value, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  facet_wrap(~ doc, labeller = label_both) + # グラフの分割
  transition_manual(sampling) + # フレーム
  labs(title = "Stochastic Variational Bayes for LDA", 
       subtitle = "s={current_frame}") # ラベル

# gif画像を作成
animate(graph_xi_theta, nframes = S + 1, fps = 20)


## 単語分布のパラメータ
# 作図
graph_xi_phi <- ggplot(trace_xi_phi_df_long, aes(x = word, y = value, fill = word, color = word)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  facet_wrap(~ topic, labeller = label_both) + # グラフの分割
  theme(legend.position = "none") + # 凡例
  scale_x_discrete(breaks = seq(0, V, by = 10)) + # x軸目盛(連続値)
  transition_manual(sampling) + # フレーム
  labs(title = "Stochastic Variational Bayes for LDA", 
       subtitle = "s={current_frame}") # ラベル

# gif画像を作成
animate(graph_xi_phi, nframes = S + 1, fps = 20)


