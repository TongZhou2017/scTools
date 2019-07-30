library(limma)
# function
venn4markerlist<-function(df_1, df_2){
	len_1<-length(levels(df_1$cluster))
	len_2<-length(levels(df_2$cluster))
	for(i in c(1:len_1)){
		gene_1<-subset(df_1,cluster==levels(df_1$cluster)[i])$gene;
		value_1<-rep(1,length(subset(df_1,cluster==levels(df_1$cluster)[i])$gene));
		dataL<-data.frame(gene=gene_1,value=value_1);
		for(j in c(1:len_2)){
			gene_2<-subset(df_2,cluster==levels(df_2$cluster)[j])$gene;
			value_2<-rep(1,length(subset(df_2,cluster==levels(df_2$cluster)[j])$gene));
			dataR<-data.frame(gene=gene_2,value=value_2);
			data_m<-merge(dataL,dataR,by="gene",all=T);
			row.names(data_m)<-data_m$gene;
			data_m<-data_m[,-1];
			data_m[is.na(data_m)]<-0;
			venn_counts <- vennCounts(data_m);
			filename=paste("venn_",levels(df_1$cluster)[i],"_",levels(df_2$cluster)[j],".pdf",sep="",collapse=NULL)
			pdf(filename);
			vennDiagram(venn_counts, cex=c(2,1),names = c(levels(df_1$cluster)[i],levels(df_2$cluster)[j]),circle.col = c("red","blue"),counts.col=NULL,show.include=NULL);
			dev.off()
		}
	}
}

# all in one
venn4markerlist.allinone<-function(df_1, df_2){
	len_1<-length(levels(df_1$cluster))
	len_2<-length(levels(df_2$cluster))
	filename="venn.allinone.pdf";
	pdf(filename);
	par(mfrow=c(len_1,len_2));
	for(i in c(1:len_1)){
		gene_1<-subset(df_1,cluster==levels(df_1$cluster)[i])$gene;
		value_1<-rep(1,length(subset(df_1,cluster==levels(df_1$cluster)[i])$gene));
		dataL<-data.frame(gene=gene_1,value=value_1);
		for(j in c(1:len_2)){
			gene_2<-subset(df_2,cluster==levels(df_2$cluster)[j])$gene;
			value_2<-rep(1,length(subset(df_2,cluster==levels(df_2$cluster)[j])$gene));
			dataR<-data.frame(gene=gene_2,value=value_2);
			data_m<-merge(dataL,dataR,by="gene",all=T);
			row.names(data_m)<-data_m$gene;
			data_m<-data_m[,-1];
			data_m[is.na(data_m)]<-0;
			venn_counts <- vennCounts(data_m);
			vennDiagram(venn_counts,names = c(levels(df_1$cluster)[i],levels(df_2$cluster)[j]),circle.col = c("red","blue"),counts.col=NULL,show.include=NULL,cex=c(0.5,0.5,0.5));	
		}
	}
	dev.off()
}
