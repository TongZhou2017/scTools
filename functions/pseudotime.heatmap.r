library(monocle)
diff_test_res <- differentialGeneTest(HSMM,fullModelFormulaStr = "~State")
top50<-head(diff_test_res[order(diff_test_res$qval),],n=50)$gene_short_name
marker_genes<-top50
pdf("pseudotime.heatmap.50.pdf")
plot_pseudotime_heatmap(HSMM[marker_genes,],num_clusters = 3,cores = 1,show_rownames = T)
dev.off()
