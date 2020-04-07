
# Stochastic Variational Bayes for LDA ------------------------------------

# 利用パッケージ
library(tidyverse)


# パラメータの設定 ----------------------------------------------------------------

# サンプリング数
S <- 100

# 試行回数
InIter <- 20

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

# 処理の確認
#sum(E_n_dk) == sum(n_dv)
#sum(E_n_dkv) == sum(n_dv)


# 確率的変分ベイズ推定 --------------------------------------------------------------

## 推移の確認用配列を作成
# 受け皿
trace_xi_theta <- array(0, dim = c(M, K, S + 1))
trace_xi_phi <- array(0, dim = c(K, V, S + 1))
# 初期値を代入
trace_xi_theta[, , 1] <- xi_dk
trace_xi_phi[, , 1] <- xi_kv

for(d in 1:M) { ## (各文書)
  
  # 動作確認
  start_time <- Sys.time()
  
  for(s in 1:S) { ## (サンプリング)
    
    for(II in 1:InIter) { ## (イタレーション)
      
      for(v in 1:V) { ## (各語彙)
        if(n_dv[d, v] > 0) {
          for(n in 1:n_dv[d, v]) { ## (各単語)
            
            # 潜在トピック集合の近似事後分布を計算:式(3.99)
            term1 <- exp(digamma(xi_kv[, v])) / exp(digamma(apply(xi_kv, 1, sum)))
            term2 <- exp(digamma(xi_dk[d, ])) / exp(digamma(apply(xi_dk, 2, sum)))
            tmp_q_z <- term1 * term2
            z_di_k[d, v, n, ] <- tmp_q_z / sum(tmp_q_z)
            
          } ## (/各単語)
        }
      } ## (/各語彙)
      
      # 文書dにおいてトピックkが割り当てられた単語数の期待値を計算
      E_n_dk[d, ] <- apply(z_di_k[d, , , ], 3, sum)
      
      # トピック分布の近似事後分布のパラメータを計算:式(3.89)
      xi_dk[d, ] <- E_n_dk[d, ] + alpha_k
      
    } ## (/イタレーション)
    
    # 文書dの語彙vにおいてトピックkが割り当てられた単語数の期待値を計算
    E_n_dkv <- apply(z_di_k, c(1, 4, 2), sum)
    
    for(k in 1:K) { ## (各トピック)
      
      # 文書番号dをサンプリング
      #d2 <- sample(1:M, size = 1)
      
      # 単語分布の近似事後分布のパラメータを計算:式(3.159)
      xi_kv[k, ] <- xi_kv[k, ] + nu * (M * E_n_dkv[d, k, ] + beta_v - xi_kv[k, ])
      
    } ## (/各トピック)
    
    # 推移の確認用配列に代入
    trace_xi_theta[d, , s + 1] <- xi_dk[d, ]
    trace_xi_phi[, , s + 1] <- trace_xi_phi[, , s + 1] + xi_kv
    
  } ## (/サンプリング)
  
  # 動作確認
  print(paste0("d=", d, ":", start_time - Sys.time()))
} ## (/各文書)



# 推定結果の確認 -----------------------------------------------------------------

## トピック分布(期待値)
# thetaの期待値を計算:式(2.10)
theta_dk <- xi_dk / apply(xi_dk, 1, sum)

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
phi_kv <- xi_kv / apply(xi_kv, 1, sum)

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
ggplot(phi_LongDF, aes(x = word, y = prob, fill = word)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  facet_wrap(~ topic, labeller = label_both) + # グラフの分割
  theme(legend.position = "none") + # 凡例
  scale_x_discrete(breaks = seq(1, V, by = 10)) + # x軸目盛(連続値)
  labs(title = "Stochastic Variational Bayes for LDA", 
       subtitle = expression(xi^phi)) # ラベル





# 推移の確認 -------------------------------------------------------------------


# 利用パッケージ
library(gganimate)

## 配列をデータフレームに変換
# 受け皿として初期値のみのデータフレームを作成
trace_xi_theta_WideDF <- cbind(
  as_tibble(trace_xi_theta[, , 1]), 
  doc = as.factor(1:M), # 文書番号
  S = 0 # 初期値
)
trace_xi_phi_WideDF <- cbind(
  as_tibble(trace_xi_phi[, , 1]), 
  topic = as.factor(1:K), # トピック番号
  S = 0 # 初期値
)

# 同様にデータフレームを行方向に結合
for(s in 1:S) {
  # データフレームを作成
  tmp_theta_df <- cbind(
    as_tibble(trace_xi_theta[, , s + 1]), 
    doc = as.factor(1:M), # 文書番号
    S = s # 試行回数
  )
  tmp_trace_xi_phi_df <- cbind(
    as_tibble(trace_xi_phi[, , s + 1]), 
    topic = as.factor(1:K), # トピック番号
    S = s # 試行回数
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


## トピック分布のパラメータの推移
# 作図
graph_xi_theta <- ggplot(trace_xi_theta_LongDF, aes(topic, value, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  facet_wrap(~ doc, labeller = label_both) + # グラフの分割
  transition_manual(S) + # フレーム
  labs(title = "Stochastic Variational Bayes for LDA", 
       subtitle = "S={current_frame}") # ラベル

# 描画
animate(graph_xi_theta, nframes = S + 1, fps = 10)


## 折れ線グラフで可視化
# 文書番号を指定
doc_num <- 10

# 作図
trace_xi_theta_LongDF %>% 
  filter(doc == doc_num) %>% 
  ggplot(aes(x = S, y = value, color = topic)) + 
    geom_line() + # 折れ線グラフ
    labs(title = "Stochastic Variational Bayes for LDA", 
         subtitle = paste0("d=", doc_num)) # ラベル


## 単語分布のパラメータの推移
# 作図
graph_xi_phi <- ggplot(trace_xi_phi_LongDF, aes(word, value, fill = word)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  facet_wrap(~ topic, labeller = label_both) + # グラフの分割
  theme(legend.position = "none") + # 凡例
  scale_x_discrete(breaks = seq(1, V, by = 10)) + # x軸目盛(連続値)
  transition_manual(S) + # フレーム
  labs(title = "Stochastic Variational Bayes for LDA", 
       subtitle = "S={current_frame}") # ラベル

# 描画
animate(graph_xi_phi, nframes = S + 1, fps = 10)

## 折れ線グラフで可視化
# トピック番号を指定
topic_num <- 3

# 作図
trace_xi_phi_LongDF %>% 
  filter(topic == topic_num) %>% 
  ggplot(aes(x = S, y = value, color = word)) + 
  geom_line(alpha = 0.5) + # 折れ線グラフ
  theme(legend.position = "none") + # 凡例
  labs(title = "Stochastic Variational Bayes for LDA", 
       subtitle = paste0("k=", topic_num)) # ラベル

