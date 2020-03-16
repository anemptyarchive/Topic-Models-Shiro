
# Gibbs sampler for LDA ---------------------------------------------------


# 利用パッケージ
library(tidyverse)


# パラメータ設定 -----------------------------------------------------------------


# トピック数
K <- 4


# ハイパーパラメータ
alpha_k <- rep(2, by = K)
beta_v <- rep(2, by = V)


# 単語分布の初期値
phi_kv <- seq(0, 1, by = 0.01) %>% 
          sample(size = K * V, replace = TRUE) %>% 
          matrix(nrow = K, ncol = V)
# 正規化
phi_kv <- phi_kv / apply(phi_kv, 1, sum)


# トピック分布の初期値
theta_dk <- seq(0, 1, by = 0.01) %>% 
            sample(size = M * K, replace = TRUE) %>% 
            matrix(nrow = M, ncol = K)
# 正規化
theta_dk <- theta_dk / apply(theta_dk, 1, sum)


# 潜在トピック集合の初期値
z_dn <- array(0, dim = c(M, V, max(N_dv)))

# 潜在トピックにおけるカウント
# 文書dでトピックkが割り当てられた単語数
n_dk <- matrix(0, nrow = M, ncol = K)

# 全文書でトピックkが割り当てられた単語数
n_kv <- matrix(0, nrow = K, ncol = V)


# サンプリング回数
S <- 1000


for(s in 1:S) {
  
  for(d in 1:M) { # 各文書
    
    for(v in 1:V) { # 各語彙
      if(N_dv[d, v] > 0) {
        for(i in 1:N_dv[d, v]) { # 各単語
          
          # サンプリング確率を計算
          p_z <- phi_kv[, v] * theta_dk[d, ] / sum(phi_kv[, v] * theta_dk[d, ])
          p_z[is.na(p_z)] <- 1
        
          # 潜在トピックを割り当て
          z_dn[d, v, i] <- sample(1:K, size = 1, prob = p_z)
        }
      }
    }
    
    # トピック分布を更新
    theta_dk[d, ] <- theta_dk[d, ] ^ (n_dk[d, ] + alpha_k - 1)
  }
  
  for(k in 1:K) { # 各トピック
    
    # 単語分布を更新
    phi_kv[k, ] <- phi_kv[k, ] ^ (n_kv[k, ] + beta_v - 1)
    
    
    ## パラメータを正規化
    phi_kv <- phi_kv / apply(phi_kv, 1, sum)
    theta_dk <- theta_dk / apply(theta_dk, 1, sum)
    
    
    ## カウントを更新
    n_dk[, k] <- apply(z_dn == k, 1, sum)
    n_kv[k, ] <- apply(z_dn == k, 2, sum)
  }
}

apply(phi_kv, 1, sum)
apply(theta_dk, 1, sum)



# 推定結果を確認 ----------------------------------------------------------------------


## 単語分布
# 
phi_WideDF <- cbind(
  as.data.frame(phi_kv), 
  topic = as.factor(1:K)
)

# 
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
  labs(title = "Gibbs sampler for LDA", 
       subtitle = expression(Phi)) # タイトル


## トピック分布
# データフレームを作成
theta_WideDF <- cbind(as.data.frame(theta_dk), 
                      doc = as.factor(1:M))

# データフレームをlong型に変換
theta_LongDF <- pivot_longer(
  theta_WideDF, 
  cols = -doc,          # 変換せずにそのまま残す現列名
  names_to = "topic",   # 現列名を格納する新しい列の名前
  names_prefix = "V",  # 現列名から取り除く文字列
  names_ptypes = list(topic = factor()),  # 現列名を要素とする際の型
  values_to = "prob"    # 現要素を格納する新しい列の名前
)

# 描画
ggplot(theta_LongDF, aes(x = topic, y = prob, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") +  # 棒グラフ
  facet_wrap( ~ doc, labeller = label_both) +        # グラフの分割
  labs(title = "Gibbs sampler for LDA", 
       subtitle = expression(Theta)) # タイトル
