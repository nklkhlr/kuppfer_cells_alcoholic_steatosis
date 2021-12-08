library(readxl)
library(tidyverse)

# ========================== #
#### Processing Functions ####
# ========================== #
glog <- function(x, l=1e-5, log.fun=log2) log.fun(x + sqrt(x^2 + l))

z_score <- function(data, sample.orient="column",
                    na.rm=TRUE) {
  switch(sample.orient,
         "column"={m <- 1},
         "row"={m <- 2},
         stop(paste("sample.orient must be 'column' or 'row' not", sample.orient))
  )
  
  
  rows <- rownames(data)
  cols <- colnames(data)
  
  scaled.data <- apply(data, m, function(x) (x-mean(x, na.rm=na.rm))/sd(x, na.rm=na.rm))
  if (sample.orient == "column") scaled.data <- t(scaled.data)
  
  rownames(scaled.data) <- rows
  colnames(scaled.data) <- cols
  
  return(scaled.data)
}

to_species_df <- function(data, is_sum_data=FALSE) {
  data_ <- as.data.frame(data)
  if (is_sum_data) {
    rownames(data_) <- paste(data$`Lipid Class`,
                             gsub("Sum ", "", data_$Species),
                             sep="_")
    return(data_[,-grep("Species|Lipid Class", colnames(data_))])
  } else {
    rownames(data_) <- data_$Species
    return(data_[,-which(colnames(data_) == "Species")])
  }
}

# ===================== #
#### Data Processing ####
# ===================== #
# for subsetting to samples used in this used (removing tg mice)
wt_sample_data <- read_excel(
  "Lipidomics MS-Daten_Alle Proben+WT only.xlsx",
  sheet="Only wt samples_paper"
)

read_excel("Sum Short 0011LI.xlsx") %>%
  filter(!is.na(.$`Lipid Class`)) -> all_data
lipid_info <- dplyr::select(
    all_data,
    matches("Species|Lipid Class| Number C-atoms|Number DBs|Number OHs")
  ) %>% filter(!grepl("Sum", Species))

# absolute species data
all_data %>%
  filter(!is.na(.$`Sum Formula`) & .$`Value Type` == "Concentration") %>%
  dplyr::select(matches("^0011LI_00|Species")) %>%
  to_species_df() %>%
  dplyr::select(wt_sample_data$`Sample Identifier`) %>%
  mutate_at(colnames(.), as.numeric) -> raw_data
# relative species
all_data %>%
  filter(!is.na(.$`Sum Formula`) & grepl("%", .$`Value Type`)) %>%
  dplyr::select(matches("^0011LI_00|Species")) %>%
  to_species_df() %>%
  dplyr::select(wt_sample_data$`Sample Identifier`) %>%
  mutate_at(colnames(.), as.numeric) -> rel_data
# absolute sum
all_data %>%
  filter(grepl("Sum", .$Species) & .$`Value Type` == "Concentration") %>%
  dplyr::select(matches("^0011LI_00|Species|Lipid Class")) %>%
  to_species_df(is_sum_data=TRUE) %>%
  dplyr::select(wt_sample_data$`Sample Identifier`) %>%
  rbind(raw_data["FC",]) %>%
  mutate_at(colnames(.), as.numeric) -> abs_sum
# relative sum
all_data %>%
  filter(grepl("Sum", .$Species) & grepl("%", .$`Value Type`)) %>%
  dplyr::select(matches("^0011LI_00|Species|Lipid Class")) %>%
  to_species_df(is_sum_data=TRUE) %>%
  dplyr::select(wt_sample_data$`Sample Identifier`) %>%
  rbind(rel_data["FC",]) %>%
  mutate_at(colnames(.), as.numeric) -> rel_sum

lipid_info <- as.data.frame(dplyr::filter(lipid_info,
                                          !duplicated(lipid_info$Species)))
rownames(lipid_info) <- lipid_info$Species

group_rename <- c("Kontrolle PBS"="Co",
                  "Kontrolle Clodronat"="Co+clo",
                  "Ethanol Diät PBS"="LDC",
                  "Ethanol Diät Clodronat"="LDC+clo")
sample_meta <- read.delim("sample_info.txt")
rownames(sample_meta) <- sample_meta$Sample.Identifier
sample_meta <- sample_meta[wt_sample_data$`Sample Identifier`,]
sample_meta$Group <- factor(group_rename[sample_meta$Group],
                            levels=c("Co", "Co+clo",
                                     "LDC", "LDC+clo"))

scaled <- z_score(raw_data)
rownames(scaled) <- rownames(raw_data)


save(lipid_info, sample_meta,
     raw_data, rel_data, abs_sum, rel_sum,
     scaled, file="lipids.RData")
