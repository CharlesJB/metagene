################################################################################
################################################################################
############################## RNA-Seq UT ######################################
################################################################################
################################################################################

### {{{ --- Test setup ---

if(FALSE) {
    library( "RUnit" )
    library( "metagene" )
    library( "data.table" )
	library( "purrr" )
}

### }}}

#########################
## Init of RNA-Seq assay 
#########################

bam_files <- 
  c(system.file("extdata/c_al4_945MLM_demo_sorted.bam", package="metagene"),
    system.file("extdata/c_al3_362PYX_demo_sorted.bam", package="metagene"),
    system.file("extdata/n_al4_310HII_demo_sorted.bam", package="metagene"),
    system.file("extdata/n_al3_588WMR_demo_sorted.bam", package="metagene"))

regions <- 
  c(system.file("extdata/ENCFF355RXX_DPM1less.bed", package="metagene"),
    system.file("extdata/ENCFF355RXX_NDUFAB1less.bed", package="metagene"),
    system.file("extdata/ENCFF355RXX_SLC25A5less.bed", package="metagene"))

mg_ori <- metagene$new(regions = regions, bam_files = bam_files, 
                       assay = 'rnaseq')
nb_nuctot <- sum(unlist(width(mg_ori$get_regions())))
nb_bam <- length(bam_files)
nb_rg <- length(regions)

mydesign <- matrix(c(1,1,1,0,0,0,0,1),ncol=2, byrow=FALSE)
mydesign <- cbind(c("c_al4_945MLM_demo_sorted.bam",
                    "c_al3_362PYX_demo_sorted.bam",
                    "n_al4_310HII_demo_sorted.bam",
                    "n_al3_588WMR_demo_sorted.bam"), mydesign)
colnames(mydesign) <- c('Samples', 'cyto', 'nucleo')
mydesign <- data.frame(mydesign)
mydesign[,2] <- as.numeric(mydesign[,2])-1
mydesign[,3] <- as.numeric(mydesign[,3])-1
#mydesign


###############################
##    produce_data_table UT
###############################

test.metagene.rna_produce_data_table_dim_checking_default_params <- function() {
  mg <- mg_ori$clone()
  mg$produce_table()
  tab <- mg$get_table()
  #nb of columns
  colnames_tab <- c("region",
                    "exon",
                    "bam",
                    "design",
                    "nuc",
                    "nuctot",
                    "exonsize",
                    "regionstartnuc",
                    "regionsize",
                    "value",
                    "strand")
  checkIdentical(dim(tab)[2], length(colnames_tab))
  #colnames
  checkIdentical(colnames(tab), colnames_tab)
  #nb of rows
  checkIdentical(nrow(tab), nb_nuctot*nb_bam)
}

test.metagene.rna_produce_data_table_dim_checking_bin100 <- function() {
  mg <- mg_ori$clone()
  mg$produce_table(bin_count = 100)
  tab <- mg$get_table()
  #nb of columns
  colnames_tab <- c("region",
                    "exon",
                    "bam",
                    "design",
                    "nuc",
                    "bin",
                    "nuctot",
                    "exonsize",
                    "regionstartnuc",
                    "regionsize",
                    "value",
                    "strand")
  checkIdentical(dim(tab)[2], length(colnames_tab))
  #colnames
  checkIdentical(colnames(tab), colnames_tab)
  #nb of rows
  checkIdentical(nrow(tab), nb_nuctot*nb_bam)
  #nb of bin
  checkIdentical(unique(tab$bin), 1:100)
  checkIdentical(dim(tab[,.N,by=c('bam','region','bin')])[1] ==  
                   nb_bam * nb_rg * 100, TRUE)
}

test.metagene.rna_produce_data_table_checking_noise_removal <- function() {
  mg <- mg_ori$clone()
  mg$produce_table()
  tab1 <- mg$get_table()
  mg$produce_table(noise_removal = "NCIS")
  tab2 <- mg$get_table()
  checkIdentical(tab1$value, tab2$value)
}

test.metagene.rna_produce_data_table_checking_normalization <- function() {
  mg <- mg_ori$clone()
  mg$produce_table()
  tab1 <- mg$get_table()
  mg$produce_table(normalization = "RPM")
  tab2 <- mg$get_table()
  checkIdentical(sum(tab1$value == tab2$value) == nb_nuctot*nb_bam, FALSE)
}

test.metagene.rna_produce_data_table_checking_flip <- function() {
  mg <- mg_ori$clone()
  mg$produce_table()
  tab1 <- mg$get_table()
  mg$produce_table(flip_regions = TRUE)
  tab2 <- mg$get_table()
  mg$produce_table()
  tab3 <- mg$get_table()
  checkIdentical(sum(tab1$nuc == tab2$nuc) == nb_nuctot*nb_bam, FALSE)
  checkIdentical(sum(tab2$nuc == tab3$nuc) == nb_nuctot*nb_bam, FALSE)
  checkIdentical(sum(tab1$nuc == tab3$nuc) == nb_nuctot*nb_bam, TRUE)
}

test.metagene.rna_produce_data_table_checking_design <- function() {
  mg <- mg_ori$clone()
  mg$produce_table()
  tab3 <- mg$get_table()
  checkIdentical(length(unique(tab3$design)), length(colnames(mg$get_design())[-1]))
  mg$produce_table(design = mydesign)
  tab1 <- mg$get_table()
  checkIdentical(length(unique(tab1$design)), length(colnames(mydesign)[-1]))
  checkIdentical(dim(tab1[design == 'nucleo',])[1] == nb_nuctot * nb_bam / 4, TRUE)
}

test.metagene.rna_produce_data_table_checking_nb_lines_by_gene <- function() {
  mg <- mg_ori$clone()
  mg$produce_table()
  tab <- mg$get_table()
  gene_length <- as.vector(rep(unlist(map(1:nb_rg, ~sum(width(mg$get_regions())[.x]))), length(bam_files)))
  gene_line <- tab[, .N, by=c('region','bam')]$N
  checkIdentical(gene_length, gene_line)
}

test.metagene.rna_produce_data_table_checking_nb_lines_by_exon <- function() {
  mg <- mg_ori$clone()
  mg$produce_table()
  tab <- mg$get_table()
  exon_length <- as.vector(rep(unlist(width(mg$get_regions())), length(bam_files)))
  exon_line <- tab[, .N, by=c('exon','region','bam')]$N
  checkIdentical(exon_length, exon_line)
}

test.metagene.rna_produce_data_table_checking_nb_lines_by_bam <- function() {
  mg <- mg_ori$clone()
  mg$produce_table()
  tab <- mg$get_table()
  bam_length <- as.vector(rep(sum(unlist(width(mg$get_regions()))), length(bam_files)))
  bam_line <- tab[, .N, by=c('bam')]$N
  checkIdentical(bam_length, bam_line)
}

###############################
##    produce_data_frame UT
###############################


test.metagene.rna_produce_data_frame_wo_design <- function() {
  mg <- mg_ori$clone()
  mg$produce_table()
  tab <- mg$get_table()
  mg$produce_data_frame()
  df <- mg$get_data_frame()
  checkIdentical(all(df[,c(-12,-13,-14)] == as.data.frame(tab)), TRUE)
}



test.metagene.rna_produce_data_frame_w_design <- function() {
  mg <- mg_ori$clone()
  mg$produce_table(design = mydesign)
  tab <- mg$get_table()
  mg$produce_data_frame()
  df <- mg$get_data_frame()
  checkIdentical(dim(df[,c(-12,-13,-14)])[1] == dim(as.data.frame(tab))[1]/2, TRUE)
}

test.metagene.rna_produce_data_frame_w_design_checking_nb_lines_by_gene <- function() {
  mg <- mg_ori$clone()
  mg$produce_table(design = mydesign)
  tab <- mg$get_table()
  mg$produce_data_frame()
  df <- mg$get_data_frame()
  gene_length <- as.vector(rep(unlist(map(1:nb_rg, ~sum(width(mg$get_regions())[.x]))), length(colnames(mydesign[-1]))))
  gene_line <- data.table(df)[, .N, by=c('region','bam')]$N
  checkIdentical(gene_length, gene_line)
}

test.metagene.rna_produce_data_frame_w_design_checking_nb_lines_by_exon <- function() {
  mg <- mg_ori$clone()
  mg$produce_table(design = mydesign)
  tab <- mg$get_table()
  mg$produce_data_frame()
  df <- mg$get_data_frame()
  exon_length <- as.vector(rep(unlist(width(mg$get_regions())), length(colnames(mydesign[-1]))))
  exon_line <- data.table(df)[, .N, by=c('exon','region','bam')]$N
  checkIdentical(exon_length, exon_line)
}

test.metagene.rna_produce_data_frame_w_design_checking_nb_lines_by_bam <- function() {
  mg <- mg_ori$clone()
  mg$produce_table(design = mydesign)
  tab <- mg$get_table()
  mg$produce_data_frame()
  df <- mg$get_data_frame()
  bam_length <- as.vector(rep(sum(unlist(width(mg$get_regions()))), length(colnames(mydesign[-1]))))
  bam_line <- data.table(df)[, .N, by=c('bam')]$N
  checkIdentical(bam_length, bam_line)
}

#ckecks for wrong arguments

test.metagene.rna_produce_data_frame_w_design_invalid_avoid_gaps <- function() {
  mg <- mg_ori$clone()
  obs <- tryCatch(mg$produce_data_frame(avoid_gaps = 'test'),
                  error = conditionMessage)
  exp <- "is.logical(avoid_gaps) is not TRUE"
  checkIdentical(obs, exp)
}

test.metagene.rna_produce_data_frame_w_design_invalid_bam_name_char <- function() {
  mg <- mg_ori$clone()
  obs <- tryCatch(mg$produce_data_frame(avoid_gaps = TRUE, 
                                        bam_name = 1234),
                  error = conditionMessage)
  exp <- "is.character(bam_name) is not TRUE"
  checkIdentical(obs, exp)
}

test.metagene.rna_produce_data_frame_w_design_bam_name_not_found <- function() {
  mg <- mg_ori$clone()
  obs <- tryCatch(mg$produce_data_frame(avoid_gaps = TRUE, 
                                        bam_name = 'test'),
                  error = conditionMessage)
  exp <- paste("bam_name argument is no one of bam_names",
               "provided to the metagene object")
  checkIdentical(obs, exp)
}

test.metagene.rna_produce_data_frame_w_design_invalid_threshold <- function() {
  mg <- mg_ori$clone()
  obs <- tryCatch(mg$produce_data_frame(avoid_gaps = TRUE, 
                                        gaps_threshold = -1),
                  error = conditionMessage)
  exp <- "gaps_threshold >= 0 is not TRUE"
  checkIdentical(obs, exp)
}

# checks for valid data table production with gaps_avoid = TRUE

test.metagene.rna_produce_data_frame_w_design_w_avoid_gaps <- function() {
  mg <- mg_ori$clone()
  mg$produce_table(design = mydesign)
  mg$produce_data_frame(avoid_gaps = TRUE)
  df_obs <- mg$get_data_frame()
  
  gaps_threshold <- 0
  bam_name <- "c_al4_945MLM_demo_sorted"
  tab <- mg$get_table()
  dfdt <- data.table(tab)
  nb_nuc_removed <- dfdt[value <= gaps_threshold 
                         & bam == bam_name, length(value),
                         by=c('exon', 'region')]
  
  nb_nuc_removed_by_gene <- dfdt[value <= gaps_threshold 
                                 & bam == bam_name, length(value),
                                 by=c('region')]
  
  checkIdentical(sum(nb_nuc_removed$V1) == sum(nb_nuc_removed_by_gene$V1), 
                 TRUE)
  
  #check nb lines identity 
  nb_design <- length(colnames(mydesign[-1]))
  nb_lines_left <- dim(dfdt)[1] - sum(nb_nuc_removed_by_gene$V1) * nb_bam
  checkIdentical(nb_lines_left/nb_bam == dim(df_obs)[1]/nb_design, TRUE)
  #idem for sum(nb_nuc_removed_by_gene) because = sum(nb_nuc_removed)
  
  exp_exon_length <- dfdt[, .N, by=c('region', 'bam', 'exon')]
  selected <- map2(nb_nuc_removed$exon,
                   nb_nuc_removed$region, 
                   ~which(exp_exon_length$exon == .x 
                          & exp_exon_length$region == .y))
  for (i in 1:length(selected)){
    exp_exon_length$N[selected[[i]]] <- exp_exon_length$N[selected[[i]]] -
      map(nb_nuc_removed$V1, ~rep(.x,each=4))[[i]]
  }
  obs_exon_length <- data.table(df_obs)[, .N, by=c('region', 'bam', 'exon')]$N
  checkIdentical(all(
    exp_exon_length$N[1:(length(exp_exon_length$N)/nb_bam*nb_design)] == 
      obs_exon_length), TRUE)
}

# checks for bin

test.metagene.rna_produce_data_frame_w_bin <- function() {
  mg <- mg_ori$clone()
  bin_count = 100
  mg$produce_table(bin_count = bin_count)
  mg$produce_data_frame()
  df <- mg$get_data_frame()
  nb_lines_exp <- bin_count * length(colnames(mg$get_design())[-1])
  nb_lines_obs <- dim(df)[1]
  checkIdentical(nb_lines_exp == nb_lines_obs, TRUE)
}

test.metagene.rna_produce_data_frame_w_bin_w_design <- function() {
  mg <- mg_ori$clone()
  bin_count = 100
  mg$produce_table(bin_count = bin_count, design = mydesign)
  mg$produce_data_frame()
  df <- mg$get_data_frame()
  nb_lines_exp <- bin_count * length(colnames(mg$get_design())[-1])
  nb_lines_obs <- dim(df)[1]
  checkIdentical(nb_lines_exp == nb_lines_obs, TRUE)
}

# checks for valid data table production with gaps_avoid = TRUE and bin

test.metagene.rna_produce_data_frame_w_bin_w_design_w_avoid_gaps <- function() {
  mg <- mg_ori$clone()
  bin_count = 100
  mg$produce_table(bin_count = bin_count, design = mydesign)
  mg$produce_data_frame(avoid_gaps = TRUE)
  df <- mg$get_data_frame()
  nb_lines_exp <- bin_count * length(colnames(mg$get_design())[-1])
  nb_lines_obs <- dim(df)[1]
  checkIdentical(nb_lines_exp == nb_lines_obs, TRUE)
}

