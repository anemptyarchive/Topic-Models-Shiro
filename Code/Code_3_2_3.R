
# Gibbs sampler for LDA ---------------------------------------------------


# 利用パッケージ
library(tidyverse)

## 動作検証用
# 文書数
M <- 10
# 語彙数
V <- 20
# 文書ごとの各語彙数
n_dv <- matrix(sample(1:3, M * V, replace = TRUE), M, V)


# パラメータの設定 -----------------------------------------------------------------


# サンプリング回数(試行回数)
S <- 1000

# トピック数
K <- 5

# 事前分布のパラメータ:(1以上の値)
alpha_k <- rep(2, K)
beta_v  <- rep(2, V)

# トピック分布の初期値
theta_dk <- seq(0, 1, by = 0.01) %>% 
            sample(size = M * K, replace = TRUE) %>% 
            matrix(nrow = M, ncol = K)
# 正規化
theta_dk <- theta_dk / apply(theta_dk, 1, sum)

# 単語分布の初期値
phi_kv <- seq(0, 1, by = 0.01) %>% 
          sample(size = K * V, replace = TRUE) %>% 
          matrix(nrow = K, ncol = V)
# 正規化
phi_kv <- phi_kv / apply(phi_kv, 1, sum)

# 潜在トピック集合の初期値
z_di <- array(0, dim = c(M, V, max(n_dv)))

# 各文書において各トピックが割り当てられた単語数
n_dk <- new_n_dk <- matrix(0, nrow = M, ncol = K)

# 全文書において各トピックが割り当てられた単語数
n_kv <- new_n_kv <- matrix(0, nrow = K, ncol = V)


# ギブスサンプリング ----------------------------------------------------------------------

# 推移の確認用データフレームを作成
trace_theta <- cbind(
  as.data.frame(theta_dk), 
  doc = as.factor(1:M),  # 文書番号
  S = 0 # 試行回数
)
trace_phi <- cbind(
  as.data.frame(phi_kv), 
  topic = as.factor(1:K),  # トピック番号
  S = 0 # サンプリング回数
)

# 推定
for(s in 1:S) { ## (イタレーション)
  
  for(d in 1:M) { ## (各文書)
    
    for(v in 1:V) { ## (各語彙)
      if(n_dv[d, v] > 0) {
        for(n in 1:n_dv[d, v]) { ## (各単語)
          
          # サンプリング確率を計算：式(3.29)
          q_z <- phi_kv[, v] * theta_dk[d, ] / sum(phi_kv[, v] * theta_dk[d, ])
          q_z[is.na(q_z)] <- 1 / K # (アンダーフロー対策(仮))
        
          # 潜在トピックを割り当て
          res_z <- rmultinom(n = 1, size = 1, prob = q_z)
          z_di[d, v, n] <- which(res_z == 1)
          
        } ## (/各単語)
      }
    } ## (/各語彙)
    
    # 事後分布のパラメータを更新：式(3.36)の指数部分
    new_alpha_k <- n_dk[d, ] + alpha_k
    
    # theta_{d,k}をサンプリング
    theta_dk[d, ] <- MCMCpack::rdirichlet(1, alpha = new_alpha_k - 1) %>% 
      as.vector()
    
  } ## (/各文書)
  
  for(k in 1:K) { ## (各トピック)
    
    # 事後分布のパラメータを更新：式(3.37)の指数部分
    new_beta_v <- n_kv[k, ] + beta_v
    
    # phi_{k,v}をサンプリング
    phi_kv[k, ] <- MCMCpack::rdirichlet(1, alpha = new_beta_v - 1) %>% 
      as.vector()
    
    
    # カウントを更新
    n_dk[, k] <- apply(z_di == k, 1, sum)
    n_kv[k, ] <- apply(z_di == k, 2, sum)
    
  } ## (/各トピック)
  
  # パラメータを正規化
  phi_kv   <- phi_kv / apply(phi_kv, 1, sum)
  theta_dk <- theta_dk / apply(theta_dk, 1, sum)
  
  
  # 推移の確認用データフレームを作成
  tmp_trace_theta <- cbind(
    as.data.frame(theta_dk), 
    doc = as.factor(1:M),  # 文書番号
    S = s # 試行回数
  )
  tmp_trace_phi <- cbind(
    as.data.frame(phi_kv), 
    topic = as.factor(1:K),  # トピック番号
    S = s # 試行回数
  )
  
  # データフレームを結合
  trace_theta <- rbind(trace_theta, tmp_trace_theta)
  trace_phi <- rbind(trace_phi, tmp_trace_phi)
}

# 処理の検証
sum(apply(n_dk, 1, sum) == apply(n_dv, 1, sum)) == M
sum(apply(n_kv, 2, sum) == apply(n_dv, 2, sum)) == V
apply(phi_kv, 1, sum)
apply(theta_dk, 1, sum)


# 推定結果の確認 ----------------------------------------------------------------------


## トピック分布
# 作図用のデータフレームを作成
theta_WideDF <- cbind(
  as.data.frame(theta_dk), 
  doc = as.factor(1:M) # 文書番号
)

# データフレームをlong型に変換
theta_LongDF <- pivot_longer(
  theta_WideDF, 
  cols = -doc,         # 変換せずにそのまま残す現列名
  names_to = "topic",  # 現列名を格納する新しい列の名前
  names_prefix = "V",  # 現列名から取り除く文字列
  names_ptypes = list(topic = factor()),  # 現列名を要素とする際の型
  values_to = "prob"   # 現要素を格納する新しい列の名前
)

# 作図
ggplot(theta_LongDF, aes(x = topic, y = prob, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  facet_wrap( ~ doc, labeller = label_both) + # グラフの分割
  labs(title = "Gibbs sampler for LDA", 
       subtitle = expression(Theta)) # ラベル


## 単語分布
# 作図用のデータフレームを作成
phi_WideDF <- cbind(
  as.data.frame(phi_kv), 
  topic = as.factor(1:K) # トピック番号
)

# データフレームをlong型に変換
phi_LongDF <- pivot_longer(
  phi_WideDF, 
  cols = -topic,       # 変換せずにそのまま残す現列名
  names_to = "word",   # 現列名を格納する新しい列の名前
  names_prefix = "V",  # 現列名から取り除く文字列
  names_ptypes = list(word = factor()),  # 現列名を要素とする際の型
  values_to = "prob"   # 現要素を格納する新しい列の名前
)

# 作図
ggplot(phi_LongDF, aes(x = word, y = prob, fill = word)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  facet_wrap( ~ topic, labeller = label_both) + # グラフの分割
  scale_x_discrete(breaks = seq(1, V, by = 10)) + # x軸目盛
  theme(legend.position = "none") + # 凡例
  labs(title = "Gibbs sampler for LDA", 
       subtitle = expression(Phi)) # ラベル


# 推移の確認用gif ---------------------------------------------------------------------

# 利用パッケージ
library(gganimate)


## トピック分布
# データフレームをlong型に変換
trace_theta_LongDF <- pivot_longer(
  trace_theta, 
  cols = -c(doc, S),   # 変換せずにそのまま残す現列名
  names_to = "topic",  # 現列名を格納する新しい列の名前
  names_prefix = "V",  # 現列名から取り除く文字列
  names_ptypes = list(topic = factor()),  # 現列名を要素とする際の型
  values_to = "prob"   # 現要素を格納する新しい列の名前
)

# 作図
graph_theta <- ggplot(trace_theta_LongDF, aes(x = topic, y = prob, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") +  # 棒グラフ
  facet_wrap( ~ doc, labeller = label_both) +        # グラフの分割
  transition_manual(S) + 
  labs(title = "Gibbs sampler for LDA", 
       subtitle = "S={current_frame}") # ラベル

# 描画
animate(graph_theta, nframes = S + 1, fps = 10)


## 単語分布
# 作図用のデータフレームを作成
trace_phi_LongDF <- pivot_longer(
  trace_phi, 
  cols = -c(topic, S), # 変換せずにそのまま残す現列名
  names_to = "word",   # 現列名を格納する新しい列の名前
  names_prefix = "V",  # 現列名から取り除く文字列
  names_ptypes = list(word = factor()),  # 現列名を要素とする際の型
  values_to = "prob"   # 現要素を格納する新しい列の名前
)

# 作図
graph_phi <- ggplot(trace_phi_LongDF, aes(x = word, y = prob, fill = word)) + 
  geom_bar(stat = "identity", position = "dodge") +  # 棒グラフ
  facet_wrap( ~ topic, labeller = label_both) +      # グラフの分割
  scale_x_discrete(breaks = seq(1, V, by = 10)) +    # x軸目盛
  theme(legend.position = "none") +                  # 凡例
  transition_manual(S) + 
  labs(title = "Gibbs sampler for LDA", 
       subtitle = "S={current_frame}") # ラベル

# 描画
animate(graph_phi, nframes = S + 1, fps = 10)


# try ---------------------------------------------------------------------

df <- tibble()
N <- 3000
for(n in 1:N) {
  tmp_df <- tibble(
    x = sample(1:10, 10, replace = TRUE), 
    y = sample(1:10, 10, replace = TRUE), 
    N = as.numeric(n)
  )
  df <- rbind(df, tmp_df)
}
df$N <- as.numeric(df$N)
graph <- ggplot(df, aes(x, y)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  transition_manual(N) + 
  labs(title = "N={current_frame}")
animate(graph, duration = N / 10)

## n_dk
# 作図用のデータフレームを作成
n_dk_WideDF <- cbind(
  as.data.frame(n_dk), 
  doc = as.factor(1:M)
)

# データフレームをlong型に変換
n_dk_LongDF <- pivot_longer(
  n_dk_WideDF, 
  cols = -doc,        # 変換せずにそのまま残す現列名
  names_to = "topic",    # 現列名を格納する新しい列の名前
  names_prefix = "V",  # 現列名から取り除く文字列
  names_ptypes = list(word = factor()),  # 現列名を要素とする際の型
  values_to = "count"    # 現要素を格納する新しい列の名前
)

# 作図
ggplot(n_dk_LongDF, aes(topic, count, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  facet_wrap(~ doc, labeller = label_both) + # 画面の分割
  labs(title = "Gibbs sampler for LDA", 
       subtitle = expression(n[dk])) # ラベル


## n_kv
# 作図用のデータフレームを作成
n_kv_WideDF <- cbind(
  as.data.frame(n_kv), 
  topic = as.factor(1:K)
)

# データフレームをlong型に変換
n_kv_LongDF <- pivot_longer(
  n_kv_WideDF, 
  cols = -topic,        # 変換せずにそのまま残す現列名
  names_to = "word",    # 現列名を格納する新しい列の名前
  names_prefix = "V",  # 現列名から取り除く文字列
  names_ptypes = list(word = factor()),  # 現列名を要素とする際の型
  values_to = "count"    # 現要素を格納する新しい列の名前
)

# 作図
ggplot(n_kv_LongDF, aes(word, count, fill = word)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  facet_wrap(~ topic, labeller = label_both) + # 画面の分割
  scale_x_discrete(breaks = seq(1, V, by = 10)) + # x軸目盛
  theme(legend.position = "none") + # 凡例
  labs(title = "Gibbs sampler for LDA", 
       subtitle = expression(n[kv])) # ラベル


