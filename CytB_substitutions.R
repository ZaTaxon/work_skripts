
'********************************median dN/dS values*******************************************'

values <- read.table(file = "/home/olga/Documents/Rodents/CytB+COI/cytB_62_spes/cyt_b_dNdS_medians", 
                     sep='\t',h=T)


ggplot(values,aes(x=Style,y=median,fill=Style))+  geom_boxplot() + 
  ylab('dN/dS values') + theme(axis.title=element_text(size=25)) + theme_bw()+ 
  scale_x_discrete(labels = " ") + xlab('') +
  scale_fill_discrete(labels=c("Subterranean","Terrestrial"), 
                      name="Rodents lifestyle")


'********************************substitutuins_density*******************************************'

subst_data <- read.csv('/home/olga/Documents/Rodents/CytB+COI/cytB_62_spes/subs_density/TreeSAAP_density',
                       h=T, sep = '\t')
subst_data$Altogether <- ifelse(subst_data$Style == 'terrestrial',54,8)
subst_data$Proportion <- subst_data$Count / subst_data$Altogether

library(ggplot2)


ggplot(data = subst_data, aes(Position, Proportion)) + geom_point(aes(color = Style))


subst_data_table <- read.csv('/home/olga/Documents/Rodents/CytB+COI/cytB_62_spes/subs_density/TreeSAAP_common_density_table',
                       h=T, sep = '\t')
nwe_df <- subset(subst_data_table, Status != 'nonsyn' & Terrestrial !=0 | Underground != 0)

library(reshape)
subst_melted <- melt(nwe_df, id = c('Status', 'Position'))
colnames(subst_melted) <- c('Status', 'Position', 'Style', 'Density')
subst_melted$Altogether <- ifelse(subst_melted$Style == 'Terrestrial',54,8)
subst_melted$Proportion <- subst_melted$Density / subst_melted$Altogether

ggplot(data = subset(subst_melted, Status = 'syn'), aes(x=Position, y= Proportion, shape = Style, color = Style)) + 
         geom_point()

ggplot(data = subst_melted, aes(x=Position, y= Proportion, color = Status , shape = Style, jitter )) + 
  geom_point()

#просто смотрим на плотности каждого  замен у наземных и подземных
ggplot(data = subset(subst_melted, Status == 'syn'), aes(x=Position, y= Proportion, color = Style)) + 
  geom_point() + ggtitle('Synonymous substitutions')

ggplot(data = subset(subst_melted, Status == 'nonsyn'), aes(x=Position, y= Proportion, color = Style)) + 
  geom_point() + ggtitle('Nonsynonymous substitutions')


#добавляем разметку по доменам белка

ggplot(data = subset(subst_melted, Status = 'nonsyn'), aes(x=Position, y= Proportion, color = Style)) + 
  geom_rect(data=NULL,aes(xmin=33,xmax=53,ymin=-Inf,ymax=Inf),fill="lightgreen") + 
  geom_point() + geom_vline(xintercept = c(33,53,77,98,113,133,178,198,226,246,288,308,320,340,347,367))

ggplot(data = subset(subst_melted, Status == 'syn' & Proportion > 0), aes(x=Position, y= Proportion, color = Style)) +
  theme_bw() + 
  geom_rect(data=NULL,aes(xmin=33,xmax=53,ymin=-Inf,ymax=Inf),fill="gainsboro", size=0) + 
  geom_rect(data=NULL,aes(xmin=77,xmax=98,ymin=-Inf,ymax=Inf),fill="gainsboro", size=0) + 
  geom_rect(data=NULL,aes(xmin=113,xmax=133,ymin=-Inf,ymax=Inf),fill="gainsboro", size=0) + 
  geom_rect(data=NULL,aes(xmin=178,xmax=198,ymin=-Inf,ymax=Inf),fill="gainsboro", size=0) + 
  geom_rect(data=NULL,aes(xmin=226,xmax=246,ymin=-Inf,ymax=Inf),fill="gainsboro", size=0) + 
  geom_rect(data=NULL,aes(xmin=288,xmax=308,ymin=-Inf,ymax=Inf),fill="gainsboro", size=0) + 
  geom_rect(data=NULL,aes(xmin=320,xmax=340,ymin=-Inf,ymax=Inf),fill="gainsboro", size=0) + 
  geom_rect(data=NULL,aes(xmin=347,xmax=367,ymin=-Inf,ymax=Inf),fill="gainsboro", size=0) + 
  geom_point() + ggtitle('Synonymous substitutions')

ggplot(data = subset(subst_melted, Status == 'nonsyn' & Proportion > 0), aes(x=Position, y= Proportion, color = Style)) +
  theme_bw() + 
  geom_rect(data=NULL,aes(xmin=33,xmax=53,ymin=-Inf,ymax=Inf),fill="gainsboro", size=0) + 
  geom_rect(data=NULL,aes(xmin=77,xmax=98,ymin=-Inf,ymax=Inf),fill="gainsboro", size=0) + 
  geom_rect(data=NULL,aes(xmin=113,xmax=133,ymin=-Inf,ymax=Inf),fill="gainsboro", size=0) + 
  geom_rect(data=NULL,aes(xmin=178,xmax=198,ymin=-Inf,ymax=Inf),fill="gainsboro", size=0) + 
  geom_rect(data=NULL,aes(xmin=226,xmax=246,ymin=-Inf,ymax=Inf),fill="gainsboro", size=0) + 
  geom_rect(data=NULL,aes(xmin=288,xmax=308,ymin=-Inf,ymax=Inf),fill="gainsboro", size=0) + 
  geom_rect(data=NULL,aes(xmin=320,xmax=340,ymin=-Inf,ymax=Inf),fill="gainsboro", size=0) + 
  geom_rect(data=NULL,aes(xmin=347,xmax=367,ymin=-Inf,ymax=Inf),fill="gainsboro", size=0) + 
  geom_point() + ggtitle('Nonsynonymous substitutions')

#просчитываем достоверность разницы частот в каждой позиции для всех вариантов замен

nwe_df$fisher <- NA
for (i in 1:nrow(nwe_df)) {
  table_data <- matrix(c(nwe_df[i,3],nwe_df[i,4], (8-nwe_df[i,3]), (54 - nwe_df[i,4])), 2,2, byrow = T)
  f <-fisher.test(x = table_data)
  nwe_df[i,5] <- f$p.value
}
nwe_df$Singif <- ifelse( nwe_df$fisher < 0.05, 'singnif', 'N.S.')

#Далем поправку Холма на множественные сравнения

nwe_df$holm <- p.adjust(nwe_df$fisher, method = 'holm')
nwe_df$holm_sign <- ifelse( nwe_df$holm < 0.05, 'singnif', 'N.S.')
subset(nwe_df, holm_sign == 'singnif')

#считаем частоты замен по доменам в целом

nwe_df$Underground_NM <- 8 - nwe_df$Underground
nwe_df$Terrestrial_NM <- 52 - nwe_df$Terrestrial

nwe_df$Domain <- NA

for (l in 1:nrow(nwe_df)) {
  pos <-  nwe_df[l, 2]
  if (pos < 33) {nwe_df[l, 11] <- 'Memb1'
  } else if (pos >= 33 & pos < 54) {
    nwe_df[l, 11] <- 'TM1'
  } else if (pos >= 54 & pos < 77) {
    nwe_df[l, 11] <- 'Memb2'
  } else if (pos >= 77 & pos < 99) {
    nwe_df[l, 11] <- 'TM2'
  } else if (pos >= 99 & pos < 113) {
    nwe_df[l, 11] <- 'Memb3'
  } else if (pos >= 113 & pos < 134) {
    nwe_df[l, 11] <- 'TM3'
  } else if (pos >= 134 & pos < 178) {
    nwe_df[l, 11] <- 'Memb4'
  } else if (pos >= 178 & pos < 199) {
    nwe_df[l, 11] <- 'TM4'
  } else if (pos >= 199 & pos < 226) {
    nwe_df[l, 11] <- 'Memb5'
  } else if (pos >= 226 & pos < 247) {
    nwe_df[l, 11] <- 'TM5'
  } else if (pos >= 247 & pos < 288) {
    nwe_df[l, 11] <- 'Memb6'
  } else if (pos >= 288 & pos < 309) {
    nwe_df[l, 11] <- 'TM6'
  } else if (pos >= 309 & pos < 320) {
    nwe_df[l, 11] <- 'Memb7'
  } else if (pos >= 320 & pos < 341) {
    nwe_df[l, 11] <- 'TM7'
  } else if (pos >= 341 & pos < 347) {
    nwe_df[l, 11] <- 'Memb8'
  } else if (pos >= 347 & pos < 368) {
    nwe_df[l, 11] <- 'TM8'
  } else if (pos >= 368) {
    nwe_df[l, 11] <- 'Memb9'
  }
}


#суммируем по доменам

domain_table_nonsyn <- subset(nwe_df[,c(1,3,4,9,10,11)], Status == 'nonsyn')
domain_table_nonsyn <- domain_table_nonsyn[,c(2:6)]
domain_summary_nonsyn <- aggregate(. ~ Domain, domain_table_nonsyn, sum)

domain_table_syn <- subset(nwe_df[,c(1,3,4,9,10,11)], Status == 'syn')
domain_table_syn <- domain_table_syn[,c(2:6)]
domain_summary_syn <- aggregate(. ~ Domain, domain_table_syn, sum)

domain_summary_nonsyn$Status <- 'nonsyn'
domain_summary_syn$Status <- 'syn'

domain_sum <- rbind(domain_summary_nonsyn, domain_summary_syn)

#считаем достоверность разницы в заменах, делаем поправку Холма

domain_sum$fisher <- NA
for (i in 1:nrow(domain_sum)) {
  table_data <- matrix(c(domain_sum[i,2],domain_sum[i,3], domain_sum[i,4], domain_sum[i,5]), 2,2, byrow = T)
  f <-fisher.test(x = table_data)
  domain_sum[i,7] <- f$p.value
}
domain_sum$Singif <- ifelse(domain_sum$fisher < 0.05, 'singnif', 'N.S.')


domain_sum$holm <- p.adjust(domain_sum$fisher, method = 'holm')
domain_sum$holm_sign <- ifelse(domain_sum$holm < 0.05, 'singnif', 'N.S.')
subset(domain_sum, holm_sign == 'singnif')

#делаем красивый рисунок для замен

domain_sum$subterranian <- domain_sum$Underground / domain_sum$Underground_NM
domain_sum$terrestrial <- domain_sum$Terrestrial / domain_sum$Terrestrial_NM
domain_dens <- domain_sum[,c(1,6,9,10,11,12)]

library(reshape)
domain_dens_melted <- melt(domain_dens, id = c('Domain', 'Status', 'holm', 'holm_sign'))
colnames(domain_dens_melted) <- c('Domain', 'Status', 'holm', 'holm_sign', 'Lifestyle', 'Dens')

domain_dens_melted$D_N <- ifelse(domain_dens_melted$Status == 'nonsyn', 
                                 domain_dens_melted$Dens, (-domain_dens_melted$Dens))

domain_dens_melted$stars <- NA
domain_dens_melted$stars[domain_dens_melted$holm < 0.05] <- '*'
domain_dens_melted$stars[domain_dens_melted$holm < 0.01] <- '**'
domain_dens_melted$stars[domain_dens_melted$holm < 0.001] <- '***'

ggplot(domain_dens_melted,aes(x=Domain,y=D_N,fill=Lifestyle))+
  geom_bar(stat = "identity", position = "dodge") + geom_hline(yintercept=0) + 
  scale_x_discrete(limits = c("Memb1","TM1", "Memb2","TM2", "Memb3","TM3", "Memb4","TM4", 
                              "Memb5","TM5", "Memb6","TM6", "Memb7","TM7", "Memb8","TM8", "Memb9"))


ggplot(domain_dens_melted,aes(x=Domain,y=D_N,fill=Lifestyle))+
  geom_bar(stat = "identity", position = "dodge") + geom_hline(yintercept=0) + 
  scale_x_discrete(limits = c("Memb1","TM1", "Memb2","TM2", "Memb3","TM3", "Memb4","TM4", 
                              "Memb5","TM5", "Memb6","TM6", "Memb7","TM7", "Memb8","TM8", "Memb9")) + 
  geom_text(data= subset(domain_dens_melted, Lifestyle == 'under_density'),
            aes(label=stars), size = 7) + xlab('CytB domains') + 
  ylab('Nonsynonymous / Synonymous substitutions')

#####Делаем красивые рисунки отдельно графиков

library(ggplot2)

ggplot(subset(domain_dens_melted , Status == 'nonsyn'),aes(x=Domain,y=Dens,fill=Lifestyle))+
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete(limits = c("Memb1","TM1", "Memb2","TM2", "Memb3","TM3", "Memb4","TM4", 
                              "Memb5","TM5", "Memb6","TM6", "Memb7","TM7", "Memb8","TM8", "Memb9")) + 
  geom_text(data= subset(domain_dens_melted , Status == 'nonsyn'& Lifestyle == 'subterranian'),
            aes(label=stars), size = 7) + xlab('Cytochrome B protein domains') + 
  ylab('Nonsynonymous substitution density') + theme(axis.title=element_text(size=22)) + theme_bw()+
  scale_fill_discrete(labels=c("Subterranean","Terrestrial"), 
                    name="Rodents lifestyle")

ggplot(subset(domain_dens_melted , Status == 'syn'),aes(x=Domain,y=Dens,fill=Lifestyle))+
  geom_bar(stat = "identity", position = "dodge") +
  scale_x_discrete(limits = c("Memb1","TM1", "Memb2","TM2", "Memb3","TM3", "Memb4","TM4", 
                              "Memb5","TM5", "Memb6","TM6", "Memb7","TM7", "Memb8","TM8", "Memb9")) + 
  geom_text(data= subset(domain_dens_melted , Status == 'syn'& Lifestyle == 'subterranian'),
            aes(label=stars), size = 7) + xlab('Cytochrome B protein domains') + 
  ylab('Synonymous substitution density') + theme(axis.title=element_text(size=22)) + theme_bw()+
  scale_fill_discrete(labels=c("Subterranean","Terrestrial"), 
                      name="Rodents lifestyle")

#сохраняем обе таблицы с заменами

write.csv(domain_sum, '/home/olga/Documents/Rodents/CytB+COI/cytB_62_spes/subs density/domain_analysis.csv',
          row.names=FALSE)

write.csv(nwe_df[,c(1:8)], '/home/olga/Documents/Rodents/CytB+COI/cytB_62_spes/subs density/site_analysis.csv',
          row.names=FALSE)


#***************обсчет паттернов замен по каждой позиции**************************************

aa_data <-  read.table('/home/olga/Documents/Rodents/CytB+COI/cytB_62_spes/aa_pattern', h=T)
str(aa_data)
Func <- function(x){
  tryCatch(fisher.test(as.matrix(x[,-1]), workspace = 10**7, simulate.p.value=TRUE)$p.val)}

Func(aa_data[aa_data$Position == 68, -2])
#aa_data[aa_data$Position %in% 2,]  %>% select(-Lifestyle) %>% group_by(Position) %>% summarize(Func(x = .))
min(p.vals)


aa_data <-  read.table('/home/olga/Documents/Rodents/CytB+COI/cytB_62_spes/aa_pattern', h=T)
p.vals <- numeric()
for(i in c(1:67, 69:75, 77:86, 88:127, 130:135, 137:138, 140:141, 143:150, 152:168, 170:175,
           177:183, 185:196, 198:203, 205:206, 208:209, 212:219, 221:222, 226:250, 252:257, 
           259:263, 265:267, 269:273, 275:277, 279:288, 291:296, 298:335, 338:345, 347:350,
           352:354, 356:380)){
  print(aa_data[aa_data$Position == i, -2])
  p.vals <- c(p.vals, Func(x = aa_data[aa_data$Position == i, -2]))}

p.val_true <- data.frame(p.vals, c(1:67, 69:75, 77:86, 88:127, 130:135, 137:138, 140:141, 143:150, 152:168, 170:175,
                        177:183, 185:196, 198:203, 205:206, 208:209, 212:219, 221:222, 226:250, 252:257, 
                        259:263, 265:267, 269:273, 275:277, 279:288, 291:296, 298:335, 338:345, 347:350,
                        352:354, 356:380))
p.val_true$adjasted <-p.adjust(p.val_true$p.vals, method = 'fdr')
p.val_true$level <- ifelse(p.val_true$adjasted < 0.05, 'sign', 'not')

table(p.val_true$adjasted)

write.csv(p.val_true, '/home/olga/Documents/Rodents/CytB+COI/cytB_62_spes/aa_pattern_pvals.csv',
          row.names=FALSE)

#делаем таблицу с частотами использования аа

aa_under <- subset(aa_data, aa_data$Lifestyle == 'Underground')
aa_terr <- subset(aa_data, aa_data$Lifestyle == 'Terrestrial')

aa_usage_density_under <- apply(aa_under[,3:22], 1, function(x) {x/sum(x)})

aa_under_dens <- as.data.frame(aa_usage_density_under)
library(data.table)
aa_under_dens_tr <- transpose(aa_under_dens)
vec <- colnames(aa_under)
colnames(aa_under_dens_tr) <- vec[3:22]
aa_under_dens_tr$position <- aa_under$Position
aa_under_dens_tr$lfst <- 'underground'

aa_usage_density_terr <- apply(aa_terr[,3:22], 1, function(x) {x/sum(x)})

aa_under_terr <- transpose(as.data.frame(aa_usage_density_terr))
colnames(aa_under_terr) <- vec[3:22]
aa_under_terr$position <- aa_under$Position
aa_under_terr$lfst <- 'terrestrial'

aa_pattern <- rbind(aa_under_terr, aa_under_dens_tr)

colnames(p.val_true) <- c('p.vals', 'position', 'adjasted', 'level')

aa_pattern_sign <- merge(aa_pattern, p.val_true, by = 'position')

write.csv(aa_pattern_sign, '/home/olga/Documents/Rodents/CytB+COI/cytB_62_spes/aa_pattern_all.csv',
          row.names=FALSE)

# table for the article
library(reshape2)
help(melt)
str(aa_pattern_sign)
aa_pattern_sign_m <- melt(aa_pattern_sign, id = c('position', 'lfst','p.vals','adjasted','level'))
aa_pattern_sign_m <- aa_pattern_sign_m[aa_pattern_sign_m$adjasted < 0.05 & aa_pattern_sign_m$value > 0,]
aa_pattern_sign_m$value <- round(aa_pattern_sign_m$value, 5)
aa_pattern_sign_m$AA_freq <- paste(aa_pattern_sign_m$variable, aa_pattern_sign_m$value, sep=':')
aa_pattern_sign_m$id <- paste(aa_pattern_sign_m$position, aa_pattern_sign_m$AA_freq, sep=':')
aa_pattern_sign_clean <- aa_pattern_sign_m[,c(1,2,4,5,8,9)]
#install.packages('tidyr')
library(tidyr)
tmp <- spread(data = aa_pattern_sign_clean, value = 'AA_freq', key = 'lfst')
tmp2 <- tmp%>% group_by(position) %>% summarize_all(~paste(unique(na.omit(.)), collapse = ', '))
write.csv(tmp2, '/home/olga/Documents/Rodents/CytB+COI/cytB_62_spes/aa_pattern_all_art.csv',
          row.names=FALSE)
    