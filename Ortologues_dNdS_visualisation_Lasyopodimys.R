#Делаем красивые графики и сводим таблицы обсчета codeml

#*********************2 species (Lasyopodimys)*****************************************
#######################free branch ratio################################################

free_table <- read.table('/home/olga/Documents/Rodents/Under_transcropts_story/all scpecies/Lasyo/free_branch_table', h=F)
library(tidyr)

colnames(free_table) <- c('gene', 'spes', 'value')
free_table_wide <- pivot_wider(free_table, names_from = spes, values_from = value)
free_Lasyo <- free_table_wide[,c(1,5,6)]

#убираем все значения, которые больше 2. Скорее всего, это ошибка рассчета

free_Lasyo_clean <- subset(free_Lasyo, L.gregalis < 2 & L.raddei < 2)

#считаем соотношение значений, чтобы понять, есть ли разница. И визуализируем это

free_Lasyo_clean$Ratio <- free_Lasyo_clean$L.gregalis / free_Lasyo_clean$L.raddei

table(free_Lasyo_clean$Ratio)
hist(free_Lasyo_clean$Ratio)
hist(free_Lasyo_clean$Ratio, breaks = seq(0,max(free_Lasyo_clean$Ratio), 0.1), xlim = c(0,5.5), 
     main = 'w (free branch model) values ratio', xlab = 'ratio of L.gregalis/L.raddei w values')

plot(free_Lasyo_clean$L.gregalis~free_Lasyo_clean$L.raddei, main = 'w (free branch model) values ratio',
     xlab='w values L.raddei', ylab = 'w values L.gregalis')
abline(a=0, b=1, col='blue')

######################branch site model##########################################
#считываем таблицы с данными (каждый вид считался отдельно) и сводим их воедино

raddei_table <- read.table('/home/olga/Documents/Rodents/Under_transcropts_story/all scpecies/Lasyo/branch_model_table_Raddei', h=F)

colnames(raddei_table) <- c('gene', 'spes', 'value')
raddei_table_wide <- pivot_wider(raddei_table, names_from = spes, values_from = value)
raddei_cut <- raddei_table_wide[,c(1,6)]

gregalis_table <- read.table('/home/olga/Documents/Rodents/Under_transcropts_story/all scpecies/Lasyo/branch_model_table_Gregalis', h=F)
colnames(gregalis_table) <- c('gene', 'spes', 'value')
greg_table_wide <- pivot_wider(gregalis_table, names_from = spes, values_from = value)
greg_cut <- raddei_table_wide[,c(1,5)]

brach_Lasyo <- merge(greg_cut, raddei_cut, by='gene')

#Высчитываем соотношения и наносим на график

branch_Lasyo_clean <- subset(brach_Lasyo, L.gregalis < 2 & L.raddei < 2)
branch_Lasyo_clean$Ratio <- branch_Lasyo_clean$L.gregalis / branch_Lasyo_clean$L.raddei

table(branch_Lasyo_clean$Ratio)
hist(branch_Lasyo_clean$Ratio)
hist(branch_Lasyo_clean$Ratio, breaks = seq(0,max(branch_Lasyo_clean$Ratio), 0.1), xlim = c(0,6.5), 
     main = 'w (branch model) values ratio', xlab = 'ratio of L.gregalis/L.raddei w values')

plot(branch_Lasyo_clean$L.gregalis~branch_Lasyo_clean$L.raddei, main = 'w (branch model) values ratio',
     xlab='w values L.raddei', ylab = 'w values L.gregalis')
abline(a=0, b=1, col='blue')

############################################QUARTILES################################################
#рассчитываем квартили для соотношений, чтобы взять 1 и 4
quantile(free_Lasyo_clean$Ratio)
quantile(branch_Lasyo_clean$Ratio)


free_Lasyo_clean$quantile <- ifelse(free_Lasyo_clean$Ratio < 7.811407e-04 | free_Lasyo_clean$Ratio > 1.493144e+00,
                                    'yes', 'no')
free_Lasyo_quantiles <- subset(free_Lasyo_clean, free_Lasyo_clean$quantile =='yes')
free_Lasyo_quantiles$order <- ifelse(free_Lasyo_quantiles$Ratio > 1.493144e+00, 4, 1)


branch_Lasyo_clean$quantile <- ifelse(branch_Lasyo_clean$Ratio < 4.106560e-01 | branch_Lasyo_clean$Ratio > 7.812500e+02,
                                      'yes', 'no')
branch_Lasyo_quantiles <- subset(branch_Lasyo_clean, branch_Lasyo_clean$quantile =='yes')
branch_Lasyo_quantiles$order <- ifelse(branch_Lasyo_quantiles$Ratio > 7.812500e+02, 4, 1)

#сохраняем таблицы для обсчета

write.table(free_Lasyo_quantiles, '/home/olga/Documents/Rodents/Under_transcropts_story/all scpecies/Lasyo/free_model_guartiles.txt', 
            sep='\t', col.names = T, row.names = F, quote = FALSE)

write.table(branch_Lasyo_quantiles, '/home/olga/Documents/Rodents/Under_transcropts_story/all scpecies/Lasyo/branch_model_guartiles.txt', 
            sep='\t', col.names = T, row.names = F, quote = FALSE)

#отрисовываем диаграмму Виенна для пересечения значений по моделям

library(VennDiagram)
grid.newpage()
draw.pairwise.venn(area1 = 53, area2 = 45, cross.area = 14, category = c("Branch model", "Free model"), 
                   lty = "blank", fill = c("skyblue", "mediumorchid"), cat.pos = c(0, 0))
