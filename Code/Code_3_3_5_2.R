
# Variational Bayes for LDA (2) -------------------------------------------

# 利用パッケージ
library(tidyverse)


# パラメータの設定 ----------------------------------------------------------------


# トピック数
K <- 5




# 変分ベイズ -------------------------------------------------------------------


for(OI in 1:OutIter) {
  
  for(d in 1:M) { # 各文書
    
    for(II in 1:InIter) {
      
      for(v in 1:V) { # 各語彙
        if(n_dv[d, v] > 0) {
          for(n in 1:n_dv[d, v]) { # 各単語
            
            # 潜在トピック集合の分布を更新:式(3.99)
            
          }
        }
      }
      
      # トピック分布を更新:式(3.90)
      
    }
  }
  
  # 単語分布を更新:式(3.96)
  
}

