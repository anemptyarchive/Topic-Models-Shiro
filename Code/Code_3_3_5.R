
# Variational Bayes for LDA (1) -------------------------------------------

# 利用パッケージ
library(tidyverse)


# パラメータの設定 ----------------------------------------------------------------


# トピック数
K <- 5


# 事前分布のパラメータ
alpha_k <- rep(2, K)
beta_v  <- rep(2, V)


# トピック分布の初期値
theta_dk <- seq(0, 1, by = 0.01) %>% 
            sample(size = M * K, replace = TRUE) %>% 
            matrix(nrow = M, ncol = K)
# 正規化
theta_dk <- theta_dk / apply(theta_dk, 1, sum)
theta_dk <- matrix(1 / K, nrow = M, ncol = K)

# 単語分布の初期値
phi_kv <- seq(0, 1, by = 0.01) %>% 
          sample(size = K * V, replace = TRUE) %>% 
          matrix(nrow = K, ncol = V)
# 正規化
phi_kv <- phi_kv / apply(phi_kv, 1, sum)
phi_kv <- matrix(1 / V, nrow = K, ncol = V)


# トピック集合
z_dvk <- array(1 / K, dim = c(M, V, K))

# 期待値
tmp_z <- array(0, dim = c(M, V, K))
for(k in 1:K) {
  tmp_z[, , k] <- z_dvk[, , k] * n_dv
}
n_dk <- apply(z_dvk, c(1, 3), sum)
n_kv <- apply(tmp_z, c(3, 2), sum)
sum(n_dv)
sum(n_dk)
sum(n_kv)


# 事後分布パラメータの初期値(=事前分布のパラメータ)
eta_dk <- t(t(n_dk) + alpha_k)
eta_kv <- t(t(n_kv) + beta_v)


# イタレーション数
Iter <- 10


# 変分ベイズ -------------------------------------------------------------------

# 推移の確認用
trace_theta <- cbind(
  as.data.frame(theta_dk), 
  doc = as.factor(1:M),  # 文書番号
  Iter = 0 # サンプリング回数
)
trace_phi <- cbind(
  as.data.frame(phi_kv), 
  topic = as.factor(1:K),  # トピック番号
  Iter = 0 # サンプリング回数
)

for(I in 1:Iter) { # 試行回数
  
  for(d in 1:M) { # 各文書
    
    for(v in 1:V) { # 各語彙
      if(n_dv[d, v] > 0) {
        for(n in 1:n_dv[d, v]) { # 各単語
          
          # 潜在トピック集合の事後分布:式(3.99)
          term1 <- digamma(eta_kv[, v]) - digamma(apply(eta_kv, 1, sum))
          term2 <- digamma(eta_dk[d, ]) - digamma(apply(eta_dk, 2, sum))
          z_dvk[d, v, ] <- exp(term1) * exp(term2)
        }
      }
    }
    
    # 期待値
    for(k in 1:K) {
      tmp_z[, , k] <- z_dvk[, , k] * n_dv
    }
    n_dk <- apply(z_dvk, c(1, 3), sum)
    n_kv <- apply(tmp_z, c(3, 2), sum)
    
    
    # 事後分布のパラメータ:式(3.89)
    eta_dk[d, ] <- n_dk[d, ] + alpha_k

    # トピック分布の事後分布を計算:式(3.90)
    theta_dk[d, ] <- theta_dk[d, ] ^ (eta_dk[d, ] - 1)
    
    for(k in 1:K) { # 各トピック
      
      # 事後分布のパラメータ:式(3.95)
      eta_kv[k, ] <- n_kv[k, ] + beta_v
      
      # 単語分布の事後分布を計算:式(3.96)
      phi_kv[k, ] <- phi_kv[k, ] ^ (eta_kv[k, ] - 1)
    }
  }
  
  # パラメータを正規化
  theta_dk <- theta_dk / apply(theta_dk, 1, sum)
  phi_kv   <- phi_kv / apply(phi_kv, 1, sum)
  
  
  # 推移の確認用
  tmp_trace_theta <- cbind(
    as.data.frame(theta_dk), 
    doc = as.factor(1:M),  # 文書番号
    Iter = I # サンプリング回数
  )
  tmp_trace_phi <- cbind(
    as.data.frame(phi_kv), 
    topic = as.factor(1:K),  # トピック番号
    Iter = I # サンプリング回数
  )
  trace_theta <- rbind(trace_theta, tmp_trace_theta)
  trace_phi <- rbind(trace_phi, tmp_trace_phi)
}



# 推定結果の確認 ----------------------------------------------------------------------


## トピック分布
# 作図用のデータフレームを作成
theta_WideDF <- cbind(
  as.data.frame(theta_dk), 
  doc = as.factor(1:M)
)

# データフレームをlong型に変換
theta_LongDF <- pivot_longer(
  theta_WideDF, 
  cols = -doc,          # 変換せずにそのまま残す現列名
  names_to = "topic",   # 現列名を格納する新しい列の名前
  names_prefix = "V",  # 現列名から取り除く文字列
  names_ptypes = list(topic = factor()),  # 現列名を要素とする際の型
  values_to = "prob"    # 現要素を格納する新しい列の名前
)

# 作図
ggplot(theta_LongDF, aes(x = topic, y = prob, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") +  # 棒グラフ
  facet_wrap( ~ doc, labeller = label_both) +        # グラフの分割
  labs(title = "Variational Bayes for LDA (1)", 
       subtitle = expression(Theta)) # ラベル


## 単語分布
# 作図用のデータフレームを作成
phi_WideDF <- cbind(
  as.data.frame(phi_kv), 
  topic = as.factor(1:K)
)

# データフレームをlong型に変換
phi_LongDF <- pivot_longer(
  phi_WideDF, 
  cols = -topic,        # 変換せずにそのまま残す現列名
  names_to = "word",    # 現列名を格納する新しい列の名前
  names_prefix = "V",  # 現列名から取り除く文字列
  names_ptypes = list(word = factor()),  # 現列名を要素とする際の型
  values_to = "prob"    # 現要素を格納する新しい列の名前
)

# 作図
ggplot(phi_LongDF, aes(x = word, y = prob, fill = word)) + 
  geom_bar(stat = "identity", position = "dodge") +  # 棒グラフ
  facet_wrap( ~ topic, labeller = label_both) +      # グラフの分割
  scale_x_discrete(breaks = seq(1, V, by = 10)) +    # x軸目盛
  theme(legend.position = "none") +                  # 凡例
  labs(title = "Variational Bayes for LDA (1)", 
       subtitle = expression(Phi)) # ラベル


# gif ---------------------------------------------------------------------


# 利用パッケージ
library(gganimate)


## トピック分布
# データフレームをlong型に変換
trace_theta_LongDF <- pivot_longer(
  trace_theta, 
  cols = -c(doc, Iter),          # 変換せずにそのまま残す現列名
  names_to = "topic",   # 現列名を格納する新しい列の名前
  names_prefix = "V",  # 現列名から取り除く文字列
  names_ptypes = list(topic = factor()),  # 現列名を要素とする際の型
  values_to = "prob"    # 現要素を格納する新しい列の名前
)

# 作図
graph_theta <- ggplot(trace_theta_LongDF, aes(x = topic, y = prob, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") +  # 棒グラフ
  facet_wrap( ~ doc, labeller = label_both) +        # グラフの分割
  transition_manual(Iter) + 
  labs(title = "Variational Bayes for LDA (1)", 
       subtitle = "Iter={current_frame}") # ラベル

# 描画
animate(graph_theta, fps = (Iter + 1) * 2)


## 単語分布
# 作図用のデータフレームを作成
trace_phi_LongDF <- pivot_longer(
  trace_phi, 
  cols = -c(topic, Iter),        # 変換せずにそのまま残す現列名
  names_to = "word",    # 現列名を格納する新しい列の名前
  names_prefix = "V",  # 現列名から取り除く文字列
  names_ptypes = list(word = factor()),  # 現列名を要素とする際の型
  values_to = "prob"    # 現要素を格納する新しい列の名前
)

# 作図
graph_phi <- ggplot(trace_phi_LongDF, aes(x = word, y = prob, fill = word)) + 
  geom_bar(stat = "identity", position = "dodge") +  # 棒グラフ
  facet_wrap( ~ topic, labeller = label_both) +      # グラフの分割
  scale_x_discrete(breaks = seq(1, V, by = 10)) +    # x軸目盛
  theme(legend.position = "none") +                  # 凡例
  transition_manual(Iter) + 
  labs(title = "Variational Bayes for LDA (1)", 
       subtitle = "Iter={current_frame}") # ラベル

# 描画
animate(graph_phi, fps = (Iter + 1) * 2)





