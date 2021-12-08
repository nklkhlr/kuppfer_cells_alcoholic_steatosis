library(mosbi)
library(ggsci)
library(igraph)

# ========================= #
#### Auxiliary Functions ####
# ========================= #
get_bic_idxs_communities <- function(comm) {
  idxs <- gsub("bicluster", "", colnames(comm@adjacency_matrix))
  return(as.integer(idxs))
}

check_pie <- function(bic_pie){

  is_all_zero <- function(x)all(x==0)

  if(any(sapply(bic_pie,is_all_zero))){
    stop("In each bicluster at least one member must have one class affiliation.")
  }
}

bicluster_pie <- function(bic, types, MARGIN="column", named=T){

  attrs <- c()
  if(MARGIN=="column"){
    if(named){
      attrs <- bic@colname
    } else{
      attrs <- bic@column
    }
  } else{
    if(named){
      attrs <- bic@rowname
    } else{
      attrs <- bic@row
    }
  }
  mapped <- types[attrs]
  if(any(is.na(mapped))){
    stop("Every element in a bicluster must have an entry in the class_vector")
  }
  num_attrs <- sort(unique(types))
  out <- rep(0,length(num_attrs))
  tmp <- table(mapped)

  for(i in 1:length(num_attrs)){
    if(!is.na(tmp[num_attrs[i]]))
      out[i] <- tmp[num_attrs[i]]
  }

  return(out)
}

biclusters_pie <- function(bics, types, MARGIN="column"){
  tmp <- lapply(bics, bicluster_pie, types=types, MARGIN=MARGIN)
  check_pie(tmp)
  return(tmp)
}

get_bicluster_indices <- function(x) {
  return(
    as.numeric(
      gsub("bicluster", "", colnames(x@adjacency_matrix))
    )
  )
}


# ================ #
#### Data Setup ####
# ================ #
# pre-computed from DataPrep.R
load("lipids.RData")
# Renaming groups for readibility
group_rename <- c("Kontrolle PBS"="Co",
                  "Kontrolle Clodronat"="Co+clo",
                  "Ethanol Diät PBS"="LDC",
                  "Ethanol Diät Clodronat"="LDC+clo")
groups <- factor(group_rename[sample_meta$Group],
                 levels=group_rename)
names(groups) <- rownames(sample_meta)

clo_mask <- grepl("clo", groups)
names(clo_mask) <- names(groups)

LDC_mask <- grepl("LDC", groups)
names(LDC_mask) <- names(groups)

# Group Colouring
group_colours <- c(
  "Co"="#B5F1BB",
  "Co+clo"="#00A658",
  "LDC"="#FFBABD",
  "LDC+clo"="#DD2A1B"
)

lipid_class <- lipid_info$`Lipid Class`
names(lipid_class) <- rownames(lipid_info)

class_colours <- sapply(pal_d3("category20c", alpha=.7)(14), substr,
                        start=1, stop=7)
names(class_colours) <- unique(lipid_class)

lipid_class_colours <- sapply(pal_d3("category20")(length(unique(lipid_info$`Lipid Class`))),
                              substr, start=1, stop=7)
names(lipid_class_colours) <- unique(lipid_info$`Lipid Class`)

lipid_classes <- factor(lipid_info$`Lipid Class`,
                        levels=sort(names(lipid_class_colours)))
names(lipid_classes) <- rownames(lipid_info)

# =================================== #
#####   Full Bicluster Networks   #####
# =================================== #
# pre-computed from MoSBi.R
load("biclusters.RData")
load("lipid_bic_communities.RData")
load("non_neutral_biclusters.RData")
load("non_neutral_bic_communities.RData")

lipid_net <- bicluster_network(lipid_bics,
                               scaled,
                               n_randomizations = 100,
                               MARGIN="both",
                               metric=4,
                               n_steps=1000,
                               plot_edge_dist=TRUE)
# removing DG, TG and CE species
neutral_mask <- grepl("DG|TG|CE", rownames(scaled))
non_neutral_net <- bicluster_network(non_neutral_bics,
                                     scaled[!neutral_mask,],
                                     n_randomizations = 100,
                                     MARGIN="both",
                                     metric=4,
                                     n_steps=1000,
                                     plot_edge_dist=TRUE)

# NOTE: Louvain modularity does not find communities in the same order.
#       If you re-run our analysis the same communities should
#       appear, but will likely have a different order, hence
#       the below selection by index might select the wrong ones

# Choosing best communities for annotation
full_comms <- c(1, 2, 3, 4)
nn_comms <- c(1, 3, 4, 5)

pdf("PaperPlots/bicluster_plot.pdf",
    width=16, height=16)
par(mar=c(c(3, 1, 3, 1) + 0.5), xpd=TRUE,
    mfrow=c(2, 1))

# a) Full bicluster graph
full_layout <- plot_piechart_bicluster_network(lipid_net, lipid_bics, groups,
                                               group_colours, vertex.label=NA,
                                               vertex.size=log2(colhistogram(lipid_bics) * 40),
                                               mark.groups=lapply(lipid_louv[full_comms], get_bic_idxs_communities),
                                               mark.border=rep("black", 8),
                                               mark.col=grDevices::rainbow(8, alpha=0),
                                               mark.shape=1)

legend("topright", legend=levels(groups),
       fill=group_colours, box.lty=0, cex=2,
       title="Group")
# b) non-neutral bicluster graph
nn_layout <- plot_piechart_bicluster_network(non_neutral_net, non_neutral_bics, groups,
                                             group_colours, vertex.label=NA,
                                             vertex.size=log2(colhistogram(lipid_bics) * 40),
                                             mark.groups=lapply(non_neutral_louv[nn_comms], get_bic_idxs_communities),
                                             mark.border=rep("black", 8),
                                             mark.col=grDevices::rainbow(8, alpha=0),
                                             mark.shape=1)
dev.off()


# ======================================== #
#####   Bicluster Community Networks   #####
# ======================================== #
# Communities with neutral lipids
signatures <- lapply(
  full_comms,
  function(i) {
    tmp_ <- select_biclusters_from_bicluster_network(lipid_louv[[i]],
                                                     lipid_bics)
    bic_idxs <- get_bicluster_indices(lipid_louv[[i]])
    node_sizes <- log2(colhistogram(lipid_bics[bic_idxs])) * 4

    pdf(paste0("PaperPlots/communities/full_", i, "_network.pdf"),
        width=16, height=12)
    par(mar=c(1,0,1,0), mfrow=c(1, 2))
    plot_piechart_bicluster_network(lipid_louv[[i]], tmp_, groups,
                                    group_colours, vertex.label=NA,
                                    new_layout=FALSE,
                                    layout=full_layout[bic_idxs,],
                                    vertex.size=12)

    plot_piechart_bicluster_network(lipid_louv[[i]], tmp_, lipid_classes,
                                    lipid_class_colours, vertex.label=NA,
                                    new_layout=FALSE,
                                    MARGIN="row",
                                    layout=full_layout[bic_idxs,],
                                    vertex.size=12)
    dev.off()
  }
)
# Communities without neutral lipids
nn_comms <- c(nn_comms, 2)
nn_signatures <- lapply(
  nn_comms,
  function(i) {
    tmp_ <- select_biclusters_from_bicluster_network(non_neutral_louv[[i]],
                                                     non_neutral_bics)
    bic_idxs <- get_bicluster_indices(non_neutral_louv[[i]])
    node_sizes <- log2(colhistogram(non_neutral_bics[bic_idxs])) * 8

    pdf(paste0("PaperPlots/communities/non_neutral_", i, "_network.pdf"),
        width=16, height=12)
    par(mar=c(1,0,1,0), mfrow=c(1, 2))
    plot_piechart_bicluster_network(non_neutral_louv[[i]], tmp_, groups,
                                    group_colours, vertex.label=NA,
                                    vertex.size=12,
                                    new_layout=FALSE,
                                    layout=nn_layout[bic_idxs,])

    plot_piechart_bicluster_network(non_neutral_louv[[i]], tmp_, lipid_classes,
                                    lipid_class_colours, vertex.label=NA,
                                    new_layout=FALSE,
                                    MARGIN="row",
                                    layout=nn_layout[bic_idxs,],
                                    vertex.size=12)
    dev.off()
  }
)
