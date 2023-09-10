library(spacexr)
library(spASE)
library(Matrix)

myrctd <- readRDS('results/rctd_hippo_1.rds')
writeMM(myrctd@originalSpatialRNA@maternalCounts, 'hippo_1_slideseq_CAST.mtx')
writeMM(myrctd@originalSpatialRNA@paternalCounts, 'hippo_1_slideseq_129.mtx')
write.table(data.frame(gene_id = seq(1,nrow(myrctd@originalSpatialRNA@maternalCounts)),
                       gene_name = rownames(myrctd@originalSpatialRNA@maternalCounts)),
            file = 'hippo_1_slideseq_genes.tsv', sep='\t', row.names=F, col.names=T, quote=F)
write.table(data.frame(pixel_id = seq(1,ncol(myrctd@originalSpatialRNA@maternalCounts)),
                       pixel_name = colnames(myrctd@originalSpatialRNA@maternalCounts)),
            file = 'hippo_1_slideseq_barcodes.tsv', sep='\t', row.names=F, col.names=T, quote=F)

myrctd <- readRDS('results/rctd_hippo_2.rds')
writeMM(myrctd@originalSpatialRNA@maternalCounts, 'hippo_2_slideseq_CAST.mtx')
writeMM(myrctd@originalSpatialRNA@paternalCounts, 'hippo_2_slideseq_129.mtx')
write.table(data.frame(gene_id = seq(1,nrow(myrctd@originalSpatialRNA@maternalCounts)),
                       gene_name = rownames(myrctd@originalSpatialRNA@maternalCounts)),
            file = 'hippo_2_slideseq_genes.tsv', sep='\t', row.names=F, col.names=T, quote=F)
write.table(data.frame(pixel_id = seq(1,ncol(myrctd@originalSpatialRNA@maternalCounts)),
                       pixel_name = colnames(myrctd@originalSpatialRNA@maternalCounts)),
            file = 'hippo_2_slideseq_barcodes.tsv', sep='\t', row.names=F, col.names=T, quote=F)

myrctd <- readRDS('results/rctd_hippo_3.rds')
writeMM(myrctd@originalSpatialRNA@maternalCounts, 'hippo_3_slideseq_CAST.mtx')
writeMM(myrctd@originalSpatialRNA@paternalCounts, 'hippo_3_slideseq_129.mtx')
write.table(data.frame(gene_id = seq(1,nrow(myrctd@originalSpatialRNA@maternalCounts)),
                       gene_name = rownames(myrctd@originalSpatialRNA@maternalCounts)),
            file = 'hippo_3_slideseq_genes.tsv', sep='\t', row.names=F, col.names=T, quote=F)
write.table(data.frame(pixel_id = seq(1,ncol(myrctd@originalSpatialRNA@maternalCounts)),
                       pixel_name = colnames(myrctd@originalSpatialRNA@maternalCounts)),
            file = 'hippo_3_slideseq_barcodes.tsv', sep='\t', row.names=F, col.names=T, quote=F)

myrctd <- readRDS('results/rctd_cere_3.rds')
writeMM(myrctd@originalSpatialRNA@maternalCounts, 'cere_3_slideseq_CAST.mtx')
writeMM(myrctd@originalSpatialRNA@paternalCounts, 'cere_3_slideseq_129.mtx')
write.table(data.frame(gene_id = seq(1,nrow(myrctd@originalSpatialRNA@maternalCounts)),
                       gene_name = rownames(myrctd@originalSpatialRNA@maternalCounts)),
            file = 'cere_3_slideseq_genes.tsv', sep='\t', row.names=F, col.names=T, quote=F)
write.table(data.frame(pixel_id = seq(1,ncol(myrctd@originalSpatialRNA@maternalCounts)),
                       pixel_name = colnames(myrctd@originalSpatialRNA@maternalCounts)),
            file = 'cere_3_slideseq_barcodes.tsv', sep='\t', row.names=F, col.names=T, quote=F)

myrctd <- readRDS('results/rctd_cere_4_visium.rds')
writeMM(myrctd@originalSpatialRNA@maternalCounts, 'cere_4_visium_CAST.mtx')
writeMM(myrctd@originalSpatialRNA@paternalCounts, 'cere_4_visium_129.mtx')
write.table(data.frame(gene_id = seq(1,nrow(myrctd@originalSpatialRNA@maternalCounts)),
                       gene_name = rownames(myrctd@originalSpatialRNA@maternalCounts)),
            file = 'cere_4_visium_genes.tsv', sep='\t', row.names=F, col.names=T, quote=F)
write.table(data.frame(pixel_id = seq(1,ncol(myrctd@originalSpatialRNA@maternalCounts)),
                       pixel_name = colnames(myrctd@originalSpatialRNA@maternalCounts)),
            file = 'cere_4_visium_barcodes.tsv', sep='\t', row.names=F, col.names=T, quote=F)
