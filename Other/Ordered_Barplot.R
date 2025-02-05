df <- read.csv(text=
                 "trt,gene,freq,cols
M6,ALDH16A1,100.0000000,red
M6,Others,0.0000000,lightgrey
M12,ALDH16A1,64.6638015,red
M12,GBE1,2.0074865,#4C00FF
M12,ZNF598,1.5832525,#004CFF
M12,CHMP6,1.3503397,#00E5FF
M12,C20orf27,1.2033828,#00FF4D
M12,NEGR1,0.9676972,#4DFF00
M12,TNFAIP6,0.9122418,#E6FF00
M12,ZSCAN25,0.7375572,#FFFF00
M12,BCL2,0.6848745,#FFDE59
M12,CBL,0.6765562,#FFE0B3
M12,Others,25.2128102,lightgrey
M18,ALDH16A1,42.4503581,red
M18,ATF2,2.2360682,#4C00FF
M18,DIAPH1,1.5256507,#004CFF
M18,SESTD1,1.2053805,#00E5FF
M18,TFCP2,1.1587958,#00FF4D
M18,SCAPER,1.1180341,#4DFF00
M18,CUX1,1.0306877,#E6FF00
M18,TEX10,0.9841030,#FFFF00
M18,C6orf89,0.9666337,#FFDE59
M18,PTTG1IP,0.9258720,#FFE0B3
M18,Others,46.3984161,lightgrey")

install.packages("tidyverse")

library(tidyverse)

df$trt <- factor(df$trt,levels=unique(as.character(df$trt)))
df$gene <- factor(df$gene,levels = unique(as.character(df$gene)))
# Reorder factor
df$gene <- forcats::fct_relevel(df$gene, "Others", after = 0)
df$gene <- forcats::fct_rev(df$gene)

# named vector of fill colors
cols <- select(df, gene, cols) %>% 
  distinct() %>% 
  deframe()

p <- ggplot(df, aes(x = trt, y = freq, fill = gene)) + 
  geom_bar(stat = "identity", color = "black") + 
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size = 4))
p
