
# Collapsed Variational Bayes Method for LDA -------------------------------------------

# 利用パッケージ
library(tidyverse)


# パラメータの設定 ----------------------------------------------------------------


# トピック数
K <- 5


# 事前分布のパラメータ
alpha_k <- rep(2, K)
beta_v  <- rep(2, V)


# 潜在トピック集合の分布
z_dvk <- array(1 / K, dim = c(M, V, K))


# イタレーション数
Iter <- 10


# 周辺化変分ベイズ ----------------------------------------------------------------


# 受け皿
new_z_dvk <- array(0, c(M, V, K))
tmp_z_dvk.di <- array(0, c(M, V, K))

for(I in 1:Iter) {
  
  for(d in 1:M) { # 各文書
    
    for(v in 1:V) { # 各語彙
      if(n_dv[d, v] > 0) {
        for(n in 1:n_dv[d, v]) { # 各単語
          
          # カウントを更新
          z_dvk.di <- z_dvk
          z_dvk.di[d, v, ] <- 0
          
          E_n_dk.di <- apply(z_dvk.di, c(1, 3), sum)
          V_n_dk.di <- apply(z_dvk.di * (1 - z_dvk.di), c(1, 3), sum)
          
          for(k in 1:K) {
            tmp_z_dvk.di[, , k] <- z_dvk.di[, , k] * n_dv
          }
          E_n_kv.di <- apply(tmp_z_dvk.di, c(3, 2), sum)
          V_n_kv.di <- apply(tmp_z_dvk.di * (1 - tmp_z_dvk.di), c(3, 2), sum)
          
          
          # 潜在トピック集合の分布を計算:式(3.130)
          term1 <- (E_n_kv.di[, v] + beta_v[v]) / apply(t(E_n_kv.di) + beta_v, 2, sum) * (E_n_dk.di[d, ] + alpha_k)
          term2.1 <- V_n_kv.di[, v] / (2 * (E_n_kv.di[, v] + beta_v[v])^2)
          term2.2 <- V_n_dk.di[d, ] / (2 * (E_n_dk.di[d, ] + alpha_k)^2)
          term3 <- apply(V_n_kv.di, 1, sum) / (2 * apply(t(E_n_kv.di) + beta_v, 2, sum)^2)
          new_z_dvk[d, v, ] <- term1 * exp(term2.1 - term2.2) * exp(term3)
        } # /各単語
      }
    } # /各語彙
  } # /各文書
  
  # 潜在トピック集合の分布を更新
  z_dvk <- new_z_dvk
}



# 推定結果の確認 -----------------------------------------------------------------


# 
z_WideDF <- NULL
for(k in 1:k) {
  tmp_z_df <- cbind(
    as.data.frame(z_dvk), 
    doc = as.factor(1:M),  # 文書番号
    topic = as.factor(k) # トピック番号
  )
  z_WideDF <- rbind(z_WideDF, tmp_z_df)
}

# データフレームをlong型に変換
z_LongDF <- pivot_longer(
  z_WideDF, 
  cols = -c(doc, topic),  # 変換せずにそのまま残す現列名
  names_to = "word",      # 現列名を格納する新しい列の名前
  names_prefix = "V",     # 現列名から取り除く文字列
  names_ptypes = list(word = factor()),  # 現列名を要素とする際の型
  values_to = "prob"      # 現要素を格納する新しい列の名前
)

# 作図
z_LongDF %>% 
  filter(doc %in% 1:5) %>% 
  filter(word %in% 1:5) %>% 
ggplot(aes(x = topic, y = prob, fill = topic)) + 
  geom_bar(stat = "identity", position = "dodge") + # 棒グラフ
  facet_wrap(doc ~ word, labeller = label_both) + # グラフの分割
#  scale_x_discrete(breaks = seq(1, V, by = 10)) + # x軸目盛
  theme(legend.position = "none") + # 凡例
  labs(title = "Collapsed Variational Bayes Method for LDA", 
       subtitle = expression(z[dn])) # ラベル



