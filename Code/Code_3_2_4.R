
# Collapsed Gibbs Sampler for LDA -----------------------------------------

# 利用パッケージ
library(tidyverse)

## 動作検証用
# 文書数
M <- 10
# 語彙数
V <- 20
# 文書ごとの各語彙数
n_dv <- matrix(sample(0:10, M * V, replace = TRUE), M, V)


# パラメータの設定 -----------------------------------------------------------------

# サンプリング回数(試行回数)を指定
S <- 500

# トピック数
K <- 5

# 事前分布のパラメータ:(0より大きい値)
alpha_k <- rep(2, K)
beta_v  <- rep(2, V)

# 潜在トピック集合の初期値
z_dvn <- array(0, dim = c(M, V, max(n_dv)))

# 各文書において各トピックが割り当てられた単語数の初期値
n_dk <- matrix(0, nrow = M, ncol = K)

# 全文書において各トピックが割り当てられた語彙ごとの単語数の初期値
n_kv <- matrix(0, nrow = K, ncol = V)


# ギブスサンプリング ----------------------------------------------------------------------

# 推移の確認用
trace_alpha <- matrix(0, nrow = K, ncol = S + 1)
trace_beta <- matrix(0, nrow = V, ncol = S + 1)
# 初期値を代入
trace_alpha[, 1] <- alpha_k
trace_beta[, 1] <- beta_v

for(s in 1:S) { ## (イタレーション)
  
  # 動作確認用
  star_time <- Sys.time()
  
  for(d in 1:M) { ## (各文書)
    
    for(v in 1:V) { ## (各語彙)
      if(n_dv[d, v] > 0) { ## (出現回数が0のときは行わない)
        for(n in 1:n_dv[d, v]) { ## (各単語)
          
          # 更新した値を移す
          n_dk.di <- n_dk
          n_kv.di <- n_kv
          
          if(s > 1) { ## (初回以外)
            
            # d,i要素のトピック番号を取り出す
            k <- z_dvn[d, v, n]
            
            # d,i要素に関する値を除く
            n_dk.di[d, k] <- n_dk.di[d, k] - 1
            n_kv.di[k, v] <- n_kv.di[k, v] - 1
            
          }
          
          # 潜在トピックのサンプリング確率を計算:式(3.38)
          tmp_kv <- log(t(t(n_kv.di) + beta_v))
          term_beta <- exp(tmp_kv[, v] - apply(tmp_kv, 1, max)) / apply(exp(tmp_kv - apply(tmp_kv, 1, max)), 1, sum) # (アンダーフロー対策)
          tmp_k <- log(n_dk.di[d, ] + alpha_k)
          term_alpha <- exp(tmp_k - max(tmp_k)) / sum(exp(tmp_k - max(tmp_k))) # (アンダーフロー対策)
          tmp_q_z_k <- term_beta * term_alpha
          q_z_k <- tmp_q_z_k / sum(tmp_q_z_k) # 正規化
          
          # 潜在トピックを割り当てる
          res_z <- rmultinom(n = 1, size = 1, prob = q_z_k)
          z_dvn[d, v, n] <- which(res_z == 1)
          
        } ## (/各単語)
      }
    } ## (/各語彙)
  } ## (/各文書)
  
  # 割り当てられたトピックに関する単語数を更新
  for(k in 1:K) {
    n_dk[, k] <- apply(z_dvn == k, 1, sum)
    n_kv[k, ] <- apply(z_dvn == k, 2, sum)
  }
  
  # トピック分布のパラメータを更新:式(3.191)
  alpha_numer <- apply(digamma(t(n_dk) + alpha_k) - digamma(alpha_k), 1, sum)
  alpha_denom <- sum(digamma(apply(t(n_dk) + alpha_k, 2, sum)) - digamma(sum(alpha_k)))
  alpha_k <- alpha_k * alpha_numer / alpha_denom
  
  # 単語分布のパラメータを更新:式(3.192)
  beta_numer <- apply(digamma(t(n_kv) + beta_v) - digamma(beta_v), 1, sum)
  beta_denom <- sum(digamma(apply(t(n_kv) + beta_v, 2, sum)) - digamma(sum(beta_v)))
  beta_v <- beta_v * beta_numer / beta_denom
  
  # 更新値の推移の確認用に値を保存
  trace_alpha[, s + 1] <- alpha_k
  trace_beta[, s + 1]  <- beta_v
  
  # 動作確認
  print(paste0(s, "th Sample...", round(Sys.time() - star_time, 3)))
}

# 処理の検証
sum(apply(n_dk, 1, sum) == apply(n_dv, 1, sum)) == M
sum(apply(n_kv, 2, sum) == apply(n_dv, 2, sum)) == V
sum(n_dk) == sum(n_kv)


# 推定結果の確認 -----------------------------------------------------------------

## トピック分布(期待値)
# thetaの期待値を計算:式(2.10)
theta_k <- alpha_k / sum(alpha_k)

# 作図用のデータフレームを作成
theta_df <- tibble(
  prob = theta_k, 
  topic = as.factor(1:K) # トピック番号
)

# 作図
ggplot(theta_df, aes(x = topic, y = prob, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  labs(title = "Collapsed Gibbs Sampler for LDA", 
       subtitle = expression(theta)) # ラベル


## 単語分布(期待値)
# phiの期待値を計算:式(2.10)
phi_v <- beta_v / sum(beta_v)

# 作図用のデータフレームを作成
phi_df <- tibble(
  prob = phi_v, 
  word = as.factor(1:V) # 語彙番号
)

# 作図
ggplot(phi_df, aes(x = word, y = prob, fill = word, color = word)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  scale_x_discrete(breaks = seq(1, V, by = 10)) + # x軸目盛
  theme(legend.position = "none") + # 凡例
  labs(title = "Collapsed Gibbs Sampler for LDA", 
       subtitle = expression(phi)) # ラベル


# 更新値の推移の確認 ---------------------------------------------------------------------

### 折れ線グラフ
## トピック分布のパラメータ
# 作図用のデータフレームを作成
trace_alpha_df_wide <- cbind(
  as.data.frame(trace_alpha), 
  topic = as.factor(1:K)
)

# データフレームをlong型に変換
trace_alpha_df_long <- pivot_longer(
  trace_alpha_df_wide, 
  cols = -topic, # 変換せずにそのまま残す現列名
  names_to = "sample", # 現列名を格納する新しい列の名前
  names_prefix = "V", # 現列名から取り除く文字列
  names_ptypes = list(sample = numeric()), # 現列名を要素とする際の型
  values_to = "value" # 現要素を格納する新しい列の名前
)

# 作図
ggplot(trace_alpha_df_long, aes(x = sample, y = value, color = topic)) + 
  geom_line() + # 折れ線グラフ
  labs(title = "Collapsed Gibbs Sampler for LDA") # ラベル


## 単語分布のパラメータ
# 作図用のデータフレームを作成
trace_beta_df_wide <- cbind(
  as.data.frame(trace_beta), 
  word = as.factor(1:V)
)

# データフレームをlong型に変換
trace_beta_df_long <- pivot_longer(
  trace_beta_df_wide, 
  cols = -word, # 変換せずにそのまま残す現列名
  names_to = "sample", # 現列名を格納する新しい列の名前
  names_prefix = "V", # 現列名から取り除く文字列
  names_ptypes = list(sample = numeric()), # 現列名を要素とする際の型
  values_to = "value" # 現要素を格納する新しい列の名前
)

# 作図
ggplot(trace_beta_df_long, aes(x = sample, y = value, color = word)) + 
  geom_line(alpha = 0.5) + # 折れ線グラフ
  theme(legend.position = "none") + # 凡例
  labs(title = "Collapsed Gibbs Sampler for LDA", 
       subtitle = expression(beta)) # ラベル


# 利用パッケージ
library(gganimate)


## トピック分布(期待値)
# 作図用のデータフレームを作成
trace_theta_df_wide <- cbind(
  as.data.frame(t(trace_alpha) / apply(trace_alpha, 2, sum)), 
  sampling = 0:S
)

# データフレームをlong型に変換
trace_theta_df_long <- pivot_longer(
  trace_theta_df_wide, 
  cols = -sampling, # 変換せずにそのまま残す現列名
  names_to = "topic", # 現列名を格納する新しい列の名前
  names_prefix = "V", # 現列名から取り除く文字列
  names_ptypes = list(sampling = numeric()), # 現列名を要素とする際の型
  values_to = "prob" # 現要素を格納する新しい列の名前
)

# 作図
graph_theta <- ggplot(trace_theta_df_long, aes(topic, prob, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") +  # 棒グラフ
  transition_manual(sampling) + # フレーム
  labs(title = "Gibbs sampler for LDA", 
       subtitle = "s={current_frame}") # ラベル

# gif画像を作成
animate(graph_theta, nframes = S + 1, fps = 10)


## 単語分布(期待値)
# 作図用のデータフレームを作成
trace_phi_df_wide <- cbind(
  as.data.frame(t(trace_beta) / apply(trace_beta, 2, sum)), 
  sampling = 0:S
)

# データフレームをlong型に変換
trace_phi_df_long <- pivot_longer(
  trace_phi_df_wide, 
  cols = -sampling, # 変換せずにそのまま残す現列名
  names_to = "word", # 現列名を格納する新しい列の名前
  names_prefix = "V", # 現列名から取り除く文字列
  names_ptypes = list(sampling = numeric()), # 現列名を要素とする際の型
  values_to = "prob" # 現要素を格納する新しい列の名前
)

# 作図
graph_phi <- ggplot(trace_phi_df_long, aes(word, prob, fill = word, color = word)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  scale_x_discrete(breaks = seq(1, V, by = 10)) + # x軸目盛
  theme(legend.position = "none") + # 凡例
  transition_manual(sampling) + # フレーム
  labs(title = "Gibbs sampler for LDA", 
       subtitle = "s={current_frame}") # ラベル

# gif画像を作成
animate(graph_phi, nframes = S + 1, fps = 10)


