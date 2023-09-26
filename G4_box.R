###Before boxplot, merge hg38_150_BGcount_boost.csv and hg38_CR_Down_count_boost.csv to sum.csv
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
piptide_file <- read.csv("/disk/yt/make_seq/result/sum.csv", row.names = NULL)
colnames(piptide_file) <- c("Type", 'count')
piptide_file <- data.frame(Type = piptide_file$Type, 
                           count = as.numeric(piptide_file$count), 
                           stringsAsFactors = FALSE)
names(piptide_file)<-c("Type","count")
my_comparisons = list(c("BG(150)", "CR_D(150)"))
my_comparisons = list(c("BG(150)", "CR_D(150)"))
cols <- c("black")
piptide_file$Type <- factor(piptide_file$Type,  levels = c("BG(150)", "CR_D(150)"))
options(scipen = 200)
c <- ggplot(piptide_file, aes(x=Type, y=count))+
  #scale_fill_manual(values = c(1_2 = "red", 1_3 = "red"))
  geom_boxplot(color="black")+
  theme_bw()+
  theme(panel.grid = element_blank(), legend.position = c('none'), 
        axis.text.x = element_text(size = 20, face = "bold.italic"), 
        axis.text.y = element_text(size = 20, face = "bold.italic"), 
        axis.title.y = element_text(size = 20))+
  stat_summary(fun.y = mean, geom = "point", shape = 23)+
  stat_compare_means(
    comparisons = my_comparisons, 
    method = "t.test", 
    label= "p.signif",
    na.rm = T
  )+
  xlab("")
ggplot2::ggsave('D:/My_data/RNA/CRIF/g4/result/box_result.pdf', c,limitsize = FALSE)