library(mosbi)
library(tidyverse)

# ========================= #
#### Auxiliary Functions ####
# ========================= #
biclust.order <- c("fabia", "isa", "biclust-qubic", "biclust-plaid",
                   "biclust-quest", "biclust-spectral")

biclust.methods <- function(data,
                            minRow=2, minCol=2,
                            unique.only=FALSE) {
  data <- as.matrix(data)

  fabi <- tryCatch(
    fabia::fabia(data, p=50) %>%
      get_biclusters(data, method="fabia") %>%
      filter_bicluster_size(minRow=2, minCol=2),
    error=function(err) {
      cat("fabia failed!\n")
      return(list())
    }
  )

  isa <- tryCatch(
    isa2::isa(data) %>%
      get_biclusters(data, method="isa") %>%
      filter_bicluster_size(minRow=2, minCol=2),
    error=function(err) {
      cat("isa failed!\n")
      return(list())
    }
  )

  qubic <- tryCatch(
    biclust::biclust(data, method=QUBIC::BCQU()) %>%
      get_biclusters(data, method="biclust-qubic") %>%
      filter_bicluster_size(minRow=2, minCol=2),
    error=function(err) {
      cat("qubic failed!\n")
      return(list())
    }
  )

  plaid <- tryCatch(
    biclust::biclust(data, method=biclust::BCPlaid()) %>%
      get_biclusters(data, method="biclust-plaid") %>%
      filter_bicluster_size(minRow=2, minCol=2),
    error=function(err) {
      cat("plaid failed!\n")
      return(list())
    }
  )

  quest <- tryCatch(
    biclust::biclust(data, method=biclust::BCQuest()) %>%
      get_biclusters(data, method="biclust-quest") %>%
      filter_bicluster_size(minRow=2, minCol=2),
    error=function(err) {
      cat("quest failed!\n")
      return(list())
    }
  )

  spect <- tryCatch(
    biclust::biclust(data, method=biclust::BCSpectral()) %>%
      get_biclusters(data, method="biclust-spectral") %>%
      filter_bicluster_size(minRow=2, minCol=2),
    error=function(err) {
      cat("spectral failed!\n")
      return(list())
    }
  )

  biclusts <- c(fabi, isa, qubic, plaid, quest, spect)

  if (unique.only) {
    return(filter_subsets(biclusts))
  } else{
    return(biclusts)
  }
}

# ================ #
#####   Main   #####
# ================ #
# computed in DataPrep.R
load("lipids.RData")

### All lipids
lipid_bics <- biclust.methods(scaled)
save(lipid_bics, file="biclusters.RData")
# community computation
# bicluster network
lipid_net <- bicluster_network(lipid_bics,
                               scaled,
                               n_randomizations = 100,
                               MARGIN="both",
                               metric=4,
                               n_steps=1000,
                               plot_edge_dist=TRUE)
# louvain communities
lipid_louv <- get_louvain_communities(lipid_net, min_size=3, bics=lipid_bics)
save(lipid_bics, file="lipid_bic_communities.RData")

### Without neutral lipids
neutral_mask <- grepl("DG|TG|CE", rownames(scaled))
non_neutral_bics <- biclust.methods(scaled[!neutral_mask,])
save(non_neutral_bics, file="non_neutral_biclusters.RData")
# community computation
# bicluster network
non_neutral_net <- bicluster_network(non_neutral_bics,
                                     scaled[!neutral_mask,],
                                     n_randomizations = 100,
                                     MARGIN="both",
                                     metric=4,
                                     n_steps=1000,
                                     plot_edge_dist=TRUE)
# louvain communities
non_neutral_louv <- get_louvain_communities(non_neutral_net, min_size=3,
                                            bics=non_neutral_bics)
save(lipid_bics, file="non_neutral_bic_communities.RData")
