# Seurat v2
for(j in 0:8){
n=table(pbmc_c@ident)
ident=j
n[names(n) %in% c(0:ident)]=0
n=n[-1]
o=n[n>0]
ident2=names(o[1])
bim <- FindMarkers(pbmc_c, ident.1 = ident,  ident.2 =ident2, only.pos = F, test.use = "bimod")
deg=bim[bim$p_val_adj < 0.01,]
deg1=deg[abs(deg$avg_logFC) >= 1,]
degs=dim(deg)
degs1=dim(deg1)

for(i in c(names(o[-1]))){
bim <- FindMarkers(pbmc_c, ident.1 = ident,  ident.2 = i, only.pos = F, test.use = "bimod")
 deg=bim[bim$p_val_adj < 0.01,]
deg1=deg[abs(deg$avg_logFC) >= 1,]
degs=cbind(degs,dim(deg))
degs1=cbind(degs1,dim(deg1))
}
print(degs1)
}

#Seurat v3
compareClusters<-function(obj){
  len<-length(levels(Idents(obj)))
  for(j in 0:len){
    n=table(Idents(obj))
    ident=j
    n[names(n) %in% c(0:ident)]=0
    n=n[-1]
    o=n[n>0]
    ident2=names(o[1])
    bim <- FindMarkers(obj, ident.1 = ident,  ident.2 =ident2, only.pos = F, test.use = "bimod")
    deg=bim[bim$p_val_adj < 0.01,]
    deg1=deg[abs(deg$avg_logFC) >= 1,]
    degs=dim(deg)
    degs1=dim(deg1)

    for(i in c(names(o[-1]))){
      bim <- FindMarkers(obj, ident.1 = ident,  ident.2 = i, only.pos = F, test.use = "bimod")
       deg=bim[bim$p_val_adj < 0.01,]
      deg1=deg[abs(deg$avg_logFC) >= 1,]
      degs=cbind(degs,dim(deg))
      degs1=cbind(degs1,dim(deg1))
    }
    print(degs1)
  }
}
