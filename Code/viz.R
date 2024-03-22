library(tibble)
library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
data_summary <- function(data, formula){
  ## formula: value ~ variables to be aggregated
  error.info <- aggregate(as.formula(formula), data, function(x) mean(is.na(x)))
  m <- aggregate(as.formula(formula), data, function(x) mean(x[!is.na(x)]))
  med <- aggregate(as.formula(formula), data, function(x) median(x[!is.na(x)]))
  se <- aggregate(as.formula(formula), data, function(x) {
    ind <- !is.na(x)
    sd(x[ind]) / sqrt(length(x[ind]))})
  sd <- aggregate(as.formula(formula), data, function(x) {
    ind <- !is.na(x)
    sd(x[ind])})
  ymin <- m[, ncol(m)] - sd[, ncol(m)]
  ymin <- ifelse(ymin > 0, ymin, 0)
  ymax <- m[, ncol(m)] + sd[, ncol(m)]
  sum <- cbind(m, median = med[, ncol(med)],SD = sd[, ncol(sd)], ymax=ymax, ymin=ymin,SE=se[, ncol(se)], ErrRate=error.info[, ncol(m)]) 
  return(sum)
}


## quantiles(93,95,97); methods(edgeR,DESeq2); datasets(1..13) 
setwd('/research/bsi/projects/staff_analysis/m216453/2022_04_08_Winsorization/Result1/')
dts <- list.files(pattern = '-permute.1.DESeq.rst.tsv$')
dts <- gsub('-permute.*','',dts)

edgeR.files <- list.files(pattern = '.edgeR.rst.tsv')
DESeq.files <- list.files(pattern = '.DESeq.rst.tsv')
edgeR_w.files <- list.files(pattern = '.edgeR_winsor.rst.tsv')
edgeR_w.files <- edgeR_w.files[grep('qt93|qt95|qt97',edgeR_w.files)]
DESeq_w.files <- list.files(pattern = '.DESeq_winsor.rst.tsv')
DESeq_w.files <- DESeq_w.files[grep('qt93|qt95|qt97',DESeq_w.files)]


compile_res <- function(files, cutoff = 0.05){
  res <- NULL
  for(file in files){
    tmp <- NULL
    if(class(structure(try({tmp = read.table(file, header = T, row.names = 1, stringsAsFactors = F,check.names = F)}))) == "try-error"){
      tmp <- 0
    }else{
      tmp <- nrow(read.table(file, header = T, row.names = 1, stringsAsFactors = F,check.names = F) %>% dplyr::filter(FDR <= cutoff))
    }
    res <- c(res, tmp)
  }
  names(res) <- files
  return(res)
}

edgeR.res <- compile_res(files = edgeR.files)
DESeq.res <- compile_res(files = DESeq.files)
edgeR_w.res <- compile_res(files = edgeR_w.files)
DESeq_w.res <- compile_res(files = DESeq_w.files)

DESeq.res0 <- as.data.frame(DESeq.res) %>% dplyr::rename(ct = DESeq.res) %>% rownames_to_column('file.name') %>% mutate(dataset = gsub('-permute.*','',file.name),tmp =  gsub('.rst.tsv','',gsub('.*-permute.','',file.name))) %>% 
  tidyr::separate(tmp, c('iter','method'),sep = "\\.") %>% 
  dplyr::select(-file.name) %>% 
  group_by(dataset, iter, method)  %>% mutate(qt = 'none',size.factor = 'none')
edgeR.res0 <- as.data.frame(edgeR.res) %>% dplyr::rename(ct = edgeR.res) %>% rownames_to_column('file.name') %>% mutate(dataset = gsub('-permute.*','',file.name),tmp =  gsub('.rst.tsv','',gsub('.*-permute.','',file.name))) %>% 
  tidyr::separate(tmp, c('iter','method'),sep = "\\.") %>% 
  dplyr::select(-file.name) %>% 
  group_by(dataset, iter, method)   %>% mutate(qt = 'none',size.factor = 'none')
DESeq_w.res0 <- as.data.frame(DESeq_w.res) %>% dplyr::rename(ct = DESeq_w.res) %>% rownames_to_column('file.name') %>% mutate(dataset = gsub('-permute.*','',file.name),tmp =  gsub('.rst.tsv','',gsub('.*-permute.','',file.name))) %>% 
  tidyr::separate(tmp, c('qt','size.factor','iter','method'),sep = "\\.") %>% 
  dplyr::select(-file.name) %>% 
  group_by(dataset, qt, size.factor, iter, method) 
edgeR_w.res0 <- as.data.frame(edgeR_w.res) %>% dplyr::rename(ct = edgeR_w.res) %>% rownames_to_column('file.name') %>% mutate(dataset = gsub('-permute.*','',file.name),tmp =  gsub('.rst.tsv','',gsub('.*-permute.','',file.name))) %>% 
  tidyr::separate(tmp, c('qt','size.factor','iter','method'),sep = "\\.") %>% 
  dplyr::select(-file.name) %>% 
  group_by(dataset, qt, size.factor, iter, method) 

res <- rbind(edgeR_w.res0,DESeq_w.res0,DESeq.res0,edgeR.res0)
res$qt <- gsub('qt93','93% quantile',res$qt)
res$qt <- gsub('qt95','95% quantile',res$qt)
res$qt <- gsub('qt97','97% quantile',res$qt)
res$method <- gsub('DESeq$','DESeq2',res$method)
res$method <- gsub('DESeq_winsor$','DESeq2',res$method)
res$method <- gsub('edgeR_winsor$','edgeR',res$method)
res$qt <- gsub('none','no-winzorization',res$qt)
save(res, file = '/research/bsi/projects/staff_analysis/m216453/2022_04_08_Winsorization/Result/res.Rdata')


############# Figure 1 ###############
load('/research/bsi/projects/staff_analysis/m216453/2022_04_08_Winsorization/Result/res.Rdata')
res$qt <- gsub('quantile','percentile',res$qt)
res <- res %>% dplyr::filter(size.factor !='TSS')
res2 <- data_summary(res, formula = 'ct ~ dataset + qt + size.factor + method') 
cols <- c( "DESeq2"= brewer.pal(8,'Reds')[7], "edgeR"= brewer.pal(8,'Greens')[c(7)])
p1 <- ggplot(res2, aes(x = method, y = ct)) + 
  geom_boxplot(aes(color = method),outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.7) + 
  scale_fill_manual(values = cols) +
  scale_color_manual(values = c(cols)) +
  facet_grid(. ~ qt, scales = 'free') + 
  theme_bw() + 
  labs(x = '', y = '# of identified DEGs \n from permuted data', fill = '') + 
  theme(axis.text.x = element_blank(),
        strip.text = element_text(size = 14),
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = 'black', size = 14),
        axis.title = element_text(color = 'black', size = 14),
        legend.title = element_text(color = 'black', size = 14),
        legend.text = element_text(color = 'black', size = 14),
        legend.position = 'none') + 
  guides(fill = FALSE)

res.ct <- res %>% dplyr::filter(size.factor !='TSS') %>% mutate(ct1 = ifelse(ct>0,1,0))
res.ct.sum <- aggregate(ct1 ~ dataset + qt + size.factor + method, data = res.ct, function(x) (sum(x)/1000)*100)
res.ct.sum <- res.ct.sum[order(res.ct.sum$ct1),]

summary((res.ct.sum %>% dplyr::filter(qt == '93% percentile' & method =='DESeq2'))$ct1)
summary((res.ct.sum %>% dplyr::filter(qt == '93% percentile' & method =='edgeR'))$ct1)
summary((res.ct.sum %>% dplyr::filter(qt == 'no-winzorization' & method =='DESeq2'))$ct1)
summary((res.ct.sum %>% dplyr::filter(qt == 'no-winzorization' & method =='edgeR'))$ct1)

p2 <- ggplot(res.ct.sum, aes(x = method, y = ct1)) + 
  geom_boxplot(aes(color = method),outlier.shape = NA) +
  geom_hline(yintercept = 5, color = 'blue', linetype='dotted') + 
  geom_jitter(width = 0.2, size = 0.7) + 
  scale_fill_manual(values = cols) +
  scale_color_manual(values = c(cols)) +
  facet_grid(. ~ qt, scales = 'free') + 
  theme_bw() + 
  labs(x = '', y = '% of permutated datasets\n with positive findings', fill = '') + 
  theme(axis.text.x = element_text(size = 14, angle = 90,hjust=0.95,vjust=0.2),
        strip.text = element_text(size = 14),
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = 'black', size = 14),
        axis.title = element_text(color = 'black', size = 14),
        legend.title = element_text(color = 'black', size = 14),
        legend.text = element_text(color = 'black', size = 14),
        legend.position = 'none') +
  ggtitle('')

p.ab <- ggarrange(p1, p2, labels = c('A','B'), nrow = 2, heights = c(0.7,1),font.label = list(size = 24, color = "black", face = "bold", family = NULL))
annotate_figure(p.ab, fig.lab = '',fig.lab.pos = "top.left",fig.lab.size = 24, fig.lab.face = 'bold') + theme(plot.margin = margin(1,0.1,0.1,0.1, "cm"))

ggsave(file = '/research/bsi/projects/staff_analysis/m216453/2022_04_08_Winsorization/plots/Figure1_n.pdf', width = 7.5, height =6)

  



## Figure 2
## extract the common findings
setwd('/research/bsi/projects/staff_analysis/m216453/2022_04_08_Winsorization/Result/original1/')
dts.w <- list.files(pattern = '.Wilcox.rst.tsv$')
ffs <- gsub('.Wilcox.rst.tsv','',dts.w)
eds <- c('-winsor.qt93.32.DESeq.rst.tsv','-winsor.qt97.32.DESeq.rst.tsv','-winsor.qt95.32.DESeq.rst.tsv',# winzorize at 0.93
         '-winsor.qt93.32.edgeR.rst.tsv','-winsor.qt97.32.edgeR.rst.tsv','-winsor.qt95.32.edgeR.rst.tsv',
         '.DESeq.rst.tsv','.edgeR.rst.tsv','.Wilcox.rst.tsv') # no winzorize
ress <- list()
for(ff in ffs){
  original <- list()
  for(ed in eds){
    cat(paste0(ff,ed),'\n')
    res.tmp <- NULL
    try({res.tmp <- read.table(paste0(ff,ed), header = T, row.names = 1, stringsAsFactors = F, check.names = F)
      if(length(grep('DESeq',ed))>0){res.tmp <- res.tmp %>% dplyr::rename(FDR = padj)}
      res.tmp <- res.tmp[,'FDR',drop=F]
      nm <- gsub('^\\.','',gsub('.rst.tsv|32\\.|\\-winsor','',ed))
      colnames(res.tmp) <- nm
      original[[nm]] <- res.tmp %>% rownames_to_column('genes')
    })
  }
  
  original <- original %>% purrr::reduce(full_join, by = "genes")
  original[is.na(original)] <- 1
  original <- original %>% column_to_rownames('genes')
  full.col <- c("qt93.DESeq","qt97.DESeq","qt95.DESeq","qt93.edgeR","qt97.edgeR","qt95.edgeR","DESeq","edgeR","Wilcox")
  if(sum(full.col %in% colnames(original)) != 9){
    miss.col <- full.col[!(full.col %in% colnames(original))]
    for(k in miss.col){
      original[,k] <- 1 
    }
  }
  ress[[ff]] <- original
}

save(ress, file = '/research/bsi/projects/staff_analysis/m216453/2022_04_08_Winsorization/Result/ress.Rdata')

load('/research/bsi/projects/staff_analysis/m216453/2022_04_08_Winsorization/Result/ress.Rdata')
comm <- nm <- xx <- ww <- NULL
for(k in names(ress)){
  df <- ress[[k]]
  for(i in c("qt93.DESeq","qt97.DESeq","qt95.DESeq","qt93.edgeR","qt97.edgeR","qt95.edgeR")){
    df1 <- df[,c(i,"Wilcox")]
    w <- sum(df1[,"Wilcox"] < 0.05)
    x <- sum(df1[,i] < 0.05)
    bt <- sum(apply(df1, 1, function(x) sum(x < 0.05)==2))
    comm <- c(comm,bt)
    # comm <- c(comm,ifelse(w!=0,(sum(apply(df1, 1, function(x) sum(x < 0.05)==2))/w) *100,100))
    nm <- c(nm, paste0(k,'-',i))
    cat('wilcox',w,' ;targte',x,'; comm:',bt, '(',paste0(k,'-',i), ')\n')
    xx <- c(xx, x)
    ww <- c(ww, w)
  }
}
names(comm) <- names(xx) <- names(ww) <- nm


comm1 <- cbind.data.frame(findings=xx, wilcox = ww, comm =comm) %>% rownames_to_column('file.name') %>% 
  mutate(file.name = gsub('TCGA\\-','TCGA\\.',file.name)) %>% 
  separate(file.name, into = c('dataset','method'), sep = '-') %>% 
  separate(method, into = c('quantile','method')) %>% 
  dplyr::mutate(quantile = gsub('qt93','93% percentile',quantile), 
                quantile = gsub('qt95','95% percentile',quantile), 
                quantile = gsub('qt97','97% percentile',quantile),
                method = gsub('DESeq','DESeq2',method), 
                ratio = ifelse(wilcox ==0 & findings ==0,1,findings/wilcox),
                comm.pct = ifelse(wilcox ==0 & comm ==0,1,comm/wilcox))
cols <- c( "DESeq2"= brewer.pal(8,'Reds')[7], "edgeR"= brewer.pal(8,'Greens')[c(7)])

p0 <- ggplot(comm1, aes(x = method, y = comm.pct)) + 
  geom_boxplot(aes(color = method),outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.7) + 
  scale_color_manual(values = c(cols)) +
  theme_bw() + 
  facet_grid(. ~ quantile) + 
  labs(x = '', y = '% of wilcox findings also found by\nwinzorization+Deseq2/edgeR', fill = '') + 
  theme(axis.text.x = element_text(size = 14, angle = 90, vjust = 0.35, hjust = 1),
        strip.text = element_text(size = 14),
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = 'black', size = 14),
        axis.title = element_text(color = 'black', size = 14),
        legend.title = element_text(color = 'black', size = 14),
        legend.text = element_text(color = 'black', size = 14),
        legend.position = 'none') + 
  guides(fill = FALSE)
p0

cols <- c( "edgeR"= brewer.pal(8,'Greens')[7], "Wilcox"= brewer.pal(8,'Blues')[c(7)])
p00 <- ggplot(comm1 %>% mutate(quantile = gsub(' percentile','', quantile)) %>% dplyr::filter(method =='edgeR'), aes(x = quantile, y = ratio)) + 
  geom_boxplot(aes(color = method),outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.7) + 
  scale_color_manual(values = c(cols)) +
  theme_bw() + 
  labs(x = 'Percentile', y = 'Ratio (# by edgeR/ # by wilcox)', fill = '') + 
  theme(axis.text.x = element_text(size = 14),
        strip.text = element_text(size = 14),
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = 'black', size = 14),
        axis.title = element_text(color = 'black', size = 14),
        legend.title = element_text(color = 'black', size = 14),
        legend.text = element_text(color = 'black', size = 14),
        legend.position = 'none') + 
  guides(fill = FALSE)
p00
library(ggpubr)
p.ab <- ggarrange(p0, p00, labels = c('A','B'), nrow = 1, widths = c(1,0.6),font.label = list(size = 24, color = "black", face = "bold", family = NULL))
annotate_figure(p.ab, fig.lab = '',fig.lab.pos = "top.left",fig.lab.size = 24, fig.lab.face = 'bold') + 
  theme(plot.margin = margin(2,0.1,0.1,0.1, "cm"))

ggsave(file = '/research/bsi/projects/staff_analysis/m216453/2022_04_08_Winsorization/plots/Figure2_n.pdf', 
       width = 9, height =5)








