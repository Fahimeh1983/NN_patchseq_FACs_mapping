##########################################################################################
### Setting up the libraries: ############################################################
##########################################################################################
script_rep="//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_script_repository/"
.libPaths("/home/fahimehb/R/x86_64-redhat-linux-gnu-library/3.5")

library(matrixStats)
library(feather)
library(tibble)
library(dendextend)
source(file.path(script_rep,"patchseq/heatmap.R"))
source(file.path(script_rep,"patchseq/de.genes.R"))
source(file.path(script_rep,"patchseq/dendro.R"))
source(file.path(script_rep,"patchseq/patchseq.R"))
source(file.path(script_rep,"patchseq/Mapping_helper_functions.R"))
source("/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/MY_R/Utils.R")
options(stringsAsFactors = F)

#TODO
batch_date="20191001_BT014-RSC-225"

######################################################################################################
### Setting up some paths ############################################################################
######################################################################################################

### TODO: Change these if you changed them in the mapping.R
ref.dir = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/facs_seq/mouse_V1_ALM_20180520/"
res.dir = "//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20180626_collapsed40_cpm/"
query.dir = "//allen/programs/celltypes/workgroups/rnaseqanalysis/SMARTer/STAR/Mouse/patchseq/R_Object/"
work.dir = "/allen//programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/"
patchseq.dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_current/"
latest.mapping.dir = "/allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/patch_seq/star/mouse_patchseq_VISp_20191001_collapsed40_cpm/"

######################################################################################################
### Reading ref data #################################################################################
######################################################################################################

### TODO: the following two lines should be modified if the number of FACS cells has changed  
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2.3.1/REF_mapping_probability.rda")
load("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2.3.1/select.markers.rda")

ggplot(data = melt(Tree_mapping_probability), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+ theme(axis.text = element_text(size=7)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("Real label") + ylab("Mapped labels") +  
  scale_fill_gradient(low = "white", high = "red")

bp.collapse.th = 40
bp.name.add = NULL
if (!is.null(bp.collapse.th)) {
  bp.name.add = paste0(".with.bp.", bp.collapse.th)
}

###load reference data and tree
tmp.load1 = load(file=file.path(res.dir, "ref.data.rda")) # should include cl, cl.df, norm.dat. # The loaded cl is not used because it only includes cluster ids but not cluster labels 
tmp.load2 = load(file.path(file=res.dir, file=paste0("V1.dend", bp.name.add,".rda"))) # should include the pruned V1 tree
tmp.load3 = load(file.path(res.dir, file=paste0("V1.dend.list", bp.name.add,".rda"))) # should include dend.list

plot(dend)

rownames(cl.df)=cl.df$cluster_id
cltmp=cl.df[as.character(cl),"cluster_label"]
names(cltmp)=names(cl)
cl=factor(cltmp)

######################################################################################################
### Loading the query data ###########################################################################
######################################################################################################

tmp<-load(paste0(query.dir,batch_date,"_mouse_patchseq_star2.0_cpm.Rdata"))
query.dat = cpmR

# loading samp.dat object
tmp<-load(paste0(query.dir,batch_date,"_mouse_patchseq_star2.0_samp.dat.Rdata"))

keepcells = which(samp.dat$Region=="VISp" & samp.dat$Type=="patch_seq")
samp.dat = samp.dat[c(keepcells, which(samp.dat$Region=="TCx"),which(samp.dat$Region=="FCx"),which(samp.dat$Region=="MOp"),which(samp.dat$Region=="TEa")   ),]   #FCx is for Brian.  Rat samples mapped in mouse

query.dat = query.dat[,as.character(samp.dat$exp_component_name)]
colnames(query.dat)=as.character(samp.dat$patched_cell_container)

query.dat.norm = log2(as.matrix(query.dat+1))
idx=match(rownames(norm.dat), rownames(query.dat.norm))
query.dat.norm=query.dat.norm[idx,]

patchseq_anno <- read_feather(paste0(patchseq.dir, "/anno.feather"))

#Patchseq Cells of interests
locked_data = read.csv(paste0(work.dir, "mouse_met_Nov11_ccf_soma_locations.csv"), check.names=FALSE, row.names = 1 )
locked_cells_spec_id = locked_data[,1]
locked_cells_sample_id = patchseq_anno[patchseq_anno$spec_id_label %in% locked_cells_spec_id, "sample_id"]
locked_cells_sample_id <- locked_cells_sample_id$sample_id
length(locked_cells_spec_id) == length(locked_cells_sample_id)
locked_cells_sample_id <- locked_cells_sample_id[locked_cells_sample_id %in% colnames(query.dat.norm)]
dim(query.dat.norm[select.markers, locked_cells_sample_id])
#write.csv(query.dat.norm, file=paste0(work.dir, "/query_dat_norm.csv"))

##########################################################################################
######################################## Mapping data using Tree #########################
##########################################################################################

set.seed(1983)
load(paste0(latest.mapping.dir, "/mapping.memb.with.bp.40.rda"))
Patchseq_Tree_memb <- memb
load(paste0(latest.mapping.dir, "/mapping.df.with.bp.40.rda"))
Patchseq_Tree_mapping.df <- mapping.df
rm(memb)
rm(mapping.df)

Patchseq_Tree_memb <- Patchseq_Tree_memb[locked_cells_sample_id,]
Patchseq_Tree_mapping.df <- Patchseq_Tree_mapping.df[locked_cells_sample_id,]
dim(Patchseq_Tree_memb)
dim(Patchseq_Tree_mapping.df)

FACs.cells <- colnames(norm.dat)

##########################################################################################
### Read mapped NN results  ##############################################################
##########################################################################################

select.cl <- labels(dend)
Patchseq_NN_memb <- read.csv(paste0(work.dir, "/NN_imb_1000epoch_500batch_patchseq_membership.csv"), check.names=FALSE, row.names = 1)
FACs_NN_memb <- read.csv(paste0(work.dir,"/NN_imb_1000epoch_500batch_FACS_membership.csv") , check.names=FALSE, row.names = 1)

#NN_FACs_memb = NN_FACs_memb / rowSums(FACs_NN_memb)
#NN_memb = NN_memb / rowSums(Patchseq_NN_memb)

colnames(Patchseq_NN_memb) == select.cl
colnames(FACs_NN_memb) == select.cl
colnames(Patchseq_NN_memb) <- select.cl
colnames(FACs_NN_memb) <- select.cl

##########################################################################################
### Building mapping probability matrix ##################################################
##########################################################################################

set.seed(1983)
#Tree_mapping_probability = compute_mapping_probability(memb = FACS_Tree_memb, select.cells = FACs.cells, 
#                                                       select.cl = select.cl, ref.cl= cl)

#save(Tree_mapping_probability, file=file.path(paste0("Final_matrices_locked_data/REF_Tree_mapping_probability.rda")))

load("//allen/programs/celltypes/workgroups/rnaseqanalysis/shiny/Taxonomies/AIT2.3.1/REF_mapping_probability.rda")

NN_mapping_probability = compute_mapping_probability(memb = FACs_NN_memb, select.cells = FACs.cells, 
                                                     select.cl = select.cl, ref.cl = cl)

ggplot(data = melt(NN_mapping_probability[select.cl, select.cl]), aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+ theme(axis.text = element_text(size=7)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  xlab("Real label") + ylab("Mapped label") +  
  scale_fill_gradient(low = "white", high = "red")


###########################################################################################
### KLdiv #################################################################################
###########################################################################################
set.seed(1983)

Patchseq_Tree_KLdiv = compute_KLdiv(select.cl = select.cl, 
                                    select.cells = locked_cells_sample_id, 
                                    mapping_probability = Tree_mapping_probability, 
                                    memb = Patchseq_Tree_memb)

Patchseq_NN_KLdiv = compute_KLdiv(select.cl = select.cl, 
                                  select.cells = locked_cells_sample_id, 
                                  mapping_probability = NN_mapping_probability, 
                                  memb = as.matrix(Patchseq_NN_memb))


##########################################################################################
### Correlation ##########################################################################
##########################################################################################
FACs_anno <- read_feather(paste0(ref.dir, "anno.feather"))
FACs_anno <- as.data.frame(FACs_anno)
rownames(FACs_anno) <- FACs_anno$sample_id
FACs_anno <- FACs_anno[FACs.cells, ] #Ask zizhen why some cells are not in norm.dat

FACs.cl.med <- Compute_median_gene_expression(anno_file = FACs_anno, norm.dat = norm.dat , markers = select.markers)


Patchseq_FACs_cor <- Compute_correlation_mat(markers = select.markers, cells = colnames(query.dat.norm), #locked_cells_sample_id, 
                                             query.dat.norm = query.dat.norm, train.cl.med = FACs.cl.med)

#FACs_FACs_cor <- Compute_correlation_mat(markers = select.markers, cells = FACs.cells, 
#                                         query.dat.norm = norm.dat, train.cl.med = FACs.cl.med)

###########################################################################################
### Aggreagte results #####################################################################
###########################################################################################

Patchseq_Tree_memb <- Patchseq_Tree_memb[locked_cells_sample_id,select.cl]
Patchseq_NN_memb <- Patchseq_NN_memb[locked_cells_sample_id, select.cl]

Tree_3_cls <- Get_3_best_cl(Patchseq_Tree_memb, select.cl)
colnames(Tree_3_cls) <- c("Tree_first_cl", "Tree_second_cl", "Tree_third_cl")

Tree_3_bts <- Get_3_best_bt(Patchseq_Tree_memb, select.cl)
colnames(Tree_3_bts) <- c("Tree_first_bt", "Tree_second_bt", "Tree_third_bt")

NN_3_cls <- Get_3_best_cl(Patchseq_NN_memb, select.cl)
colnames(NN_3_cls) <- c("NN_first_cl", "NN_second_cl", "NN_third_cl")

NN_3_bts <- Get_3_best_bt(Patchseq_NN_memb, select.cl)
colnames(NN_3_bts) <- c("NN_first_bt", "NN_second_bt", "NN_third_bt")


Tree_3_KL <- Get_3_best_KL(memb = Patchseq_Tree_memb, ref.cl = select.cl, KLdiv = Patchseq_Tree_KLdiv)
colnames(Tree_3_KL) <-c("Tree_first_KL", "Tree_second_KL", "Tree_third_KL")

NN_3_KL <- Get_3_best_KL(memb = Patchseq_NN_memb, ref.cl = select.cl, KLdiv = Patchseq_NN_KLdiv)
colnames(NN_3_KL) <- c("NN_first_KL", "NN_second_KL", "NN_third_KL")

Tree_3_cor <- Get_3_best_cor(memb = Patchseq_Tree_memb, ref.cl = select.cl, cor = Patchseq_FACs_cor)
colnames(Tree_3_cor) <- c("Tree_first_cor", "Tree_second_cor", "Tree_third_cor")

NN_3_cor <- Get_3_best_cor(memb = Patchseq_NN_memb, ref.cl = select.cl, cor = Patchseq_FACs_cor)
colnames(NN_3_cor) <- c("NN_first_cor", "NN_second_cor", "NN_third_cor")


##########################################################################################
### Assign cell identities: ##############################################################
##########################################################################################
patchseq_anno <- as.data.frame(patchseq_anno)
rownames(patchseq_anno) <- patchseq_anno$sample_id
cells <- locked_cells_sample_id
results <- cbind.data.frame(Tree_3_cls[cells,],
                            NN_3_cls[cells,],
                            Tree_3_bts[cells,],
                            NN_3_bts[cells,],
                            Tree_3_KL[cells,],
                            NN_3_KL[cells,],
                            Tree_3_cor[cells,],
                            NN_3_cor[cells,],
                            patchseq_anno[cells, c("topLeaf_id", "topLeaf_label", "topLeaf_color")])#,

Original_cols <- colnames(results)

results <- results %>% 
  rownames_to_column("id") %>%
  mutate(Tree_not_finall_call = ifelse(Tree_first_cor > 0.5  & Tree_first_KL < 2, "Good", "PoorQ")) %>%
  mutate(Tree_call = case_when(Tree_not_finall_call == "Good" & Tree_first_bt >= 0.9 ~ "Core",
                               Tree_not_finall_call == "Good" & Tree_first_bt < 0.9 &
                                 Tree_first_bt + Tree_second_bt >= 0.7 &
                                 Tree_first_bt / Tree_second_bt >= 2 ~ "I1", 
                               Tree_not_finall_call == "Good" & Tree_first_bt < 0.9 &
                                 Tree_first_bt + Tree_second_bt >= 0.7 &
                                 Tree_first_bt / Tree_second_bt < 2 ~ "I2",
                               Tree_not_finall_call == "Good" & Tree_first_bt < 0.9 &
                                 Tree_first_bt + Tree_second_bt < 0.7 ~ "I3",
                               Tree_not_finall_call == "PoorQ" ~ "PoorQ",
                               TRUE ~ "Other")) %>%
  mutate(NN_not_finall_call = ifelse(NN_first_cor > 0.5 & NN_first_KL < 2, "Good", "PoorQ")) %>%
  mutate(NN_call = case_when(NN_not_finall_call == "Good" & NN_first_bt >= 0.85 ~ "Core",
                             NN_not_finall_call == "Good" & NN_first_bt < 0.85 &
                               NN_first_bt + NN_second_bt >= 0.6 &
                               NN_first_bt / NN_second_bt >= 2 ~ "I1", 
                             NN_not_finall_call == "Good" & NN_first_bt < 0.85 &
                               NN_first_bt + NN_second_bt >= 0.6 &
                               NN_first_bt / NN_second_bt < 2 ~ "I2",
                             NN_not_finall_call == "Good" & NN_first_bt < 0.85 &
                               NN_first_bt + NN_second_bt < 0.6 ~ "I3",
                             NN_not_finall_call == "PoorQ" ~ "PoorQ",
                             TRUE ~ "Other")) %>%
  mutate(Tree_id = case_when(Tree_call == "Core" ~ 1,
                             Tree_call == "I1" ~ 2,
                             Tree_call == "I2" ~ 3,
                             Tree_call == "I3" ~ 4,
                             Tree_call == "PoorQ" ~ 5)) %>%
  mutate(Tree_color = case_when(Tree_call == "Core" ~ "Green",
                                Tree_call == "I1" ~ "Blue",
                                Tree_call == "I2" ~ "red",
                                Tree_call == "I3" ~ "Orange",
                                Tree_call == "PoorQ" ~ "Purple")) %>%
  mutate(NN_id = case_when(NN_call == "Core" ~ 1,
                           NN_call == "I1" ~ 2,
                           NN_call == "I2" ~ 3,
                           NN_call == "I3" ~ 4,
                           NN_call == "PoorQ" ~ 5)) %>%
  mutate(NN_color = case_when(NN_call == "Core" ~ "Green",
                              NN_call == "I1" ~ "Blue",
                              NN_call == "I2" ~ "red",
                              NN_call == "I3" ~ "Orange",
                              NN_call == "PoorQ" ~ "Purple")) %>%
  column_to_rownames("id") 

results <- results[,c(Original_cols, "Tree_call", "NN_call", "Tree_color", "NN_color", "Tree_id", "NN_id")]
#sum(results$Tree_first_cl == results$NN_first_cl & results$Tree_call=="Core"  & results$NN_call=="Core" )
#sum(results$Tree_first_cl == results$NN_first_cl & (results$Tree_call=="Core" | results$Tree_call=="I1") & (results$NN_call=="Core" | results$NN_call=="I1"))
#sum(results$Tree_first_cl == results$NN_first_cl & (results$Tree_call=="Core" | results$Tree_call=="I1" | results$Tree_call=="I2") & 
#      (results$NN_call=="Core" | results$NN_call=="I1" | results$NN_call=="I2"))
#sum(results$Tree_first_cl == results$NN_first_cl & (results$Tree_call=="Core" | results$Tree_call=="I1" | results$Tree_call=="I2" | results$Tree_call=="I3") & 
#      (results$NN_call=="Core" | results$NN_call=="I1" | results$NN_call=="I2" | results$NN_call=="I3"))

#sum(results$Tree_call!="PoorQ"  & results$NN_call!="PoorQ" )


#sum(results$Tree_call!="PoorQ" & results$NN_call!="PoorQ")

results[locked_cells_sample_id,"Old_cluster"] <- Patchseq_Tree_mapping.df[locked_cells_sample_id, "cl"]
Leaf_node_cells <- rownames(Patchseq_Tree_mapping.df)[Patchseq_Tree_mapping.df$resolution.index==1]
Internal_node_cells <- rownames(Patchseq_Tree_mapping.df)[(Patchseq_Tree_mapping.df$resolution.index < 1 &Patchseq_Tree_mapping.df$resolution.index > 0.7)]
PoorQ_cells <- setdiff(rownames(Patchseq_Tree_mapping.df), c(Leaf_node_cells, Internal_node_cells))
Leaf_node_cells <- intersect(Leaf_node_cells, locked_cells_sample_id)
Internal_node_cells <- intersect(Internal_node_cells, locked_cells_sample_id)
PoorQ_cells <- intersect(PoorQ_cells, locked_cells_sample_id)

results[Leaf_node_cells, "Old_call"] <- c("Leaf_node")
results[Internal_node_cells, "Old_call"] <- c("Internal_node")
results[PoorQ_cells, "Old_call"] <- c("PoorQ")

results[Leaf_node_cells, "Old_color"] <- c("Green")
results[Internal_node_cells, "Old_color"] <- c("Blue")
results[PoorQ_cells, "Old_color"] <- c("Purple")

results[Leaf_node_cells, "Old_id"] <- c(1)
results[Internal_node_cells, "Old_id"] <- c(2)
results[PoorQ_cells, "Old_id"] <- c(3)

ggplot(results, aes(Tree_first_bt , fill = Tree_call)) + 
  geom_density(alpha = 0.3) + xlim(c(0,1)) + 
  xlab("first cluster call confidence") + ylab("Density")

ggplot(results, aes(Tree_second_bt  , fill = Tree_call)) + 
  geom_density(alpha = 0.3) + xlim(c(0,1)) +
  xlab("second cluster call confidence") + ylab("Density")

ggplot(results, aes(NN_first_bt , fill = NN_call)) + 
  geom_density(alpha = 0.3) + xlim(c(0,1)) + 
  xlab("first cluster call confidence") + ylab("Density")

ggplot(results, aes(NN_second_bt  , fill = NN_call)) + 
  geom_density(alpha = 0.3) + xlim(c(0,1)) +
  xlab("second cluster call confidence") + ylab("Density")

write.csv(results, file = "/allen/programs/celltypes/workgroups/rnaseqanalysis/Fahimehb/Manuscript_patchseq_2019/NN_Tree_anno_results_Jan2020.csv")


levels <- cl.df$cluster_label
rownames(cl.df) <-cl.df$cluster_label
NN.cl <- set_names(results$NN_first_cl, rownames(results))
NN.cl = droplevels(factor(NN.cl, levels = levels))
Tree.cl <- setNames(results$Tree_first_cl, rownames(results))
Tree.cl = droplevels(factor(Tree.cl, levels = levels))

library(scrattch.hicat)
compare_annotate(cl= NN.cl, ref.cl = Tree.cl, 
                                    ref.cl.df = cl.df, 
                                    reorder = TRUE, 
                                    rename = FALSE)

