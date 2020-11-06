#1.
library(tidyverse)
m = readRDS('data_BrainSpan.rds')

g = m$g
e = m$e
s = m$s %>% mutate(col_id = paste(structure_acronym, donor_id,gender,age,sep= "_"))
colnames(e) <- s$col_id
rownames(e) <- g$gene_symbol
e$gene_symbol <- g$gene_symbol
e1 <- e %>% gather(col_id,expression,-gene_symbol)
e1$region <- NA
e1$region <- ifelse(grepl("MFC", e1$col_id),"MFC","MD")
e1 %>% ggplot(aes(region,expression)) + geom_boxplot() + geom_point() + scale_y_continuous(trans = 'log2')

#2.
e1$gene_type <- NA
e1$gene_type <- ifelse(grepl("protein_coding",g$gene_type),"protein_coding","non_protein_coding")
e1 %>% ggplot(aes(gene_type, log2(expression+1))) + geom_boxplot() 

#3
e1$age <- s$age
e1$age <- ifelse(grepl('pcw', e1$age),as.numeric(gsub('pcw','',e1$age)),e1$age)
e1$age <- ifelse(grepl('mos',e1$age),as.numeric(gsub('mos','',e1$age)) * 30 + 365 , e1$age)
e1$age <- ifelse(grepl('yrs',e1$age),as.numeric(gsub('yrs','',e1$age)) * 365 + 365, e1$age)
e1$age <- as.numeric(e1$age)
e1 %>% filter(gene_symbol=='SCN2A') %>% arrange(age) %>% ggplot(aes(age, log2(expression+1))) + geom_point() + coord_flip()

#4

e1 %>% filter(gene_symbol=='SCN2A') %>% arrange(as.numeric(age)) %>% ggplot(aes(age, log2(expression+1))) + geom_point() + coord_flip() + facet_grid(region~.,scale='free')

#5

e2 <- tibble(age = e1$age[grepl("^SCN",e1$gene_symbol)], sodium_expression = e1$expression[grepl("^SCN",e1$gene_symbol)])
e3 <- tibble(age = e1$age[grepl("^SCNN",e1$gene_symbol)], sodium_expression = e1$expression[grepl("^SCNN", e1$gene_symbol)])
e4 <- tibble(age = e1$age[grepl("^NALCN",e1$gene_symbol)], sodium_expression = e1$expression[grepl("^NALCN", e1$gene_symbol)])
e5 <- tibble(age = e1$age[grepl("^ASIC",e1$gene_symbol)], sodium_expression = e1$expression[grepl("^ASIC", e1$gene_symbol)])
e6 <- rbind(e2,e3,e4,e5)
e6 %>% ggplot(aes(age, log2(sodium_expression+1))) + geom_point()
#Reference : https://www.genenames.org/data/genegroup/#!/group/179

#6

e7 <- tibble(age = e1$age[grepl("^KCN",e1$gene_symbol)], potassium_expression = e1$expression[grepl("^KCN",e1$gene_symbol)])
e7 %>% ggplot(aes(age,log2(potassium_expression+1))) + geom_point()
#Reference : https://www.genenames.org/data/genegroup/#!/group/183

#7

e1$sex <- NA
e1$sex <- ifelse(grepl("_M_",e1$col_id),"Male","Female")
e1 %>% filter(gene_symbol=="XIST") %>% ggplot(aes(age, log2(expression+1))) + geom_point() + facet_wrap(sex~.)

#8
table(e1$sex)
e1$expression <- as.numeric(e1$expression, na.rm =T)
Male <- e1 %>% filter(sex == "Male") %>% group_by(gene_symbol) %>% summarize(med = median(expression, na.rm = T))
Female <- e1 %>% filter(sex == "Female") %>% group_by(gene_symbol) %>% summarize(med = median(expression, na.rm = T))
diff <- abs(Male$med - Female$med)
top <- order(desc(diff))[1:3]
diff[top]
top3 <- Male$gene_symbol[top]
e1 %>% filter(gene_symbol %in% top3) %>% ggplot(aes(gene_symbol, log2(expression+1))) + geom_point() + facet_grid(.~sex)

Male %>% filter(gene_symbol == "XIST") %>% summarize(me = median(expression, na.rm = T))














