library(plotly)
library(shiny)
library(gapminder)
library(RColorBrewer)
library(stringr)
library(caret)
library(SummarizedExperiment)
library(tidyr)
library(dplyr)
library(DESeq2)
library(crosstalk)
library(DT)
library(umap)
library(shinyjs)
library(shinybusy)
library(shinyBS)
library(shinyWidgets)
library(RANN)
library(comprehenr)



#library(reticulate)
# library(uwot)
# library(bibliometrix)


js <- "
function(data, type, row, meta) {
  if(type === 'sort') {
    data = Math.abs(data);
  }
  return data;
}
"

collapsibleCheckboxGroupInput <- 
  function(inputId, label, i, choices = NULL, selected = NULL, width = NULL, 
           choiceNames = NULL, choiceValues = NULL){
    input <- checkboxGroupInput(inputId, label, choices = choices, 
                                selected = selected, width = width,
                                choiceNames = choiceNames, 
                                choiceValues = choiceValues)
    checkboxes <- input[[3]][[2]][[3]][[1]]
    id_btn <- paste0(inputId, "_btn")
    id_div <- paste0(inputId, "_collapsible")
    btn <- actionButton(id_btn, "More...", 
                        icon = icon("collapse-up", lib = "glyphicon"), 
                        class = "btn-primary btn-sm", 
                        `data-toggle`="collapse", 
                        `data-target` = paste0("#", id_div))
    collapsible <- div(id = id_div, class = "collapse")
    collapsible$children <- checkboxes[(i+1):length(checkboxes)]
    children <- c(checkboxes[1:i], list(btn), list(collapsible))
    input[[3]][[2]][[3]][[1]] <- children
    script <- sprintf('$(document).ready(function(){
      $("#%s_collapsible").on("hide.bs.collapse", function(){
        $("#%s_btn").html("<span class=\\\"glyphicon glyphicon-collapse-down\\\"></span> More...");
      });
      $("#%s_collapsible").on("show.bs.collapse", function(){
        $("#%s_btn").html("<span class=\\\"glyphicon glyphicon-collapse-up\\\"></span> Less...");
      });
    });', inputId, inputId, inputId, inputId)
    tagList(input, tags$script(HTML(script)))
  }


collapsibleAwesomeCheckboxGroupInput <- 
  function(inputId, label, i, choices = NULL, selected = NULL,  
           status = "primary", width = NULL){
    input <- awesomeCheckboxGroup(inputId, label, choices = choices, 
                                  selected = selected, width = width,
                                  status = status)
    checkboxes <- input[[3]][[2]][[3]][[1]]
    id_btn <- paste0(inputId, "_btn")
    id_div <- paste0(inputId, "_collapsible")
    btn <- actionBttn(id_btn, "More...", color = "primary", size = "sm", 
                      style = "minimal", icon = icon("collapse-up", lib = "glyphicon"))
    collapsible <- div(id = id_div, class = "collapse")
    collapsible$children <- checkboxes[(i+1):length(checkboxes)]
    children <- c(checkboxes[1:i], list(btn), list(collapsible))
    input[[3]][[2]][[3]][[1]] <- children
    script <- sprintf('$(document).ready(function(){
      $("#%s_btn").attr("data-target", "#%s_collapsible").attr("data-toggle", "collapse").css("margin-bottom", "11px");
      $("#%s_collapsible").on("hide.bs.collapse", function(){
        $("#%s_btn").html("<span class=\\\"glyphicon glyphicon-collapse-down\\\"></span> More...");
      });
      $("#%s_collapsible").on("show.bs.collapse", function(){
        $("#%s_btn").html("<span class=\\\"glyphicon glyphicon-collapse-up\\\"></span> Less...");
      });
    });', inputId, inputId, inputId, inputId, inputId, inputId)
    tagList(input, tags$script(HTML(script)))
  }


options(shiny.maxRequestSize = 30*1024^2)
options(warn=-1)

### FUNCTIONS ############


formatSamples <- function(df){
  # adjust sample_names format
  sample_names <- sub("X", "", colnames(df))
  sample_names <- sample_names %>% str_replace_all("\\.", "\\-") 
  return(sample_names)
}



parseInput <- function(df){
  # returns an expression dataframe 
  # rows=genes, cols=samples
  # sorted according to gene ids
  # print(df[1:3,1:5])
  
  #
  #print(dim(df))
  df <- subset(df, !(df[,1] %in% pseudogenes))
  #print(dim(df))
  #print(df[1:3,1:3])
  df <- df[order(df[,1]), ]
  #print(dim(df))
  #print(df[1:3,1:3])
  #print(dim(df))
  
  gene_ids <- df[,1] 
  expr_df <- df[1:nrow(df), 3:ncol(df)]
  row.names(expr_df) <- gene_ids
  sample_names <- formatSamples(expr_df)
  colnames(expr_df) <- sample_names
  #print(dim(expr_df))
  #print(expr_df[1:3,1:3])
  return(expr_df)
  
}



# SCALING function
scaling_method <- function(mat, method){
  if(method=="Unit variance"){
    # rows=samples, cols=genes
    scaled <- scale(mat, scale=T, center=T)
  }
  if(method=="VST"){
    # for vst transformation rows=genes, cols=samples
    scaled <- vst(mat, blind=T)
  }
  if(method=="Log"){
    # rows=samples, cols=genes
    # add pseudocounts to avoid -INF 
    scaled <- apply(mat, 2, function(x) log2(x+1))
  } 
  if(method == "MinMax"){
    pp <- preProcess(mat, method="range")
    scaled <- predict(pp, mat)
  }
  return(scaled)
}


compute_PCA <- function(x, colData, scaling_method, gene_names){
  if(scaling_method=="unitvar"){
    pca <- prcomp(x, center=T, retx=T, scale.=T)
  }
  else {
    pca <- prcomp(x, center=T, retx=T)
  }
  
  explained_variance_ratio <- summary(pca)[["importance"]]['Proportion of Variance',]
  explained_variance_ratio <- 100 * explained_variance_ratio
  
  # attach metadata
  df_pca <- as.data.frame(pca$x)[1:10]
  df_pca <- cbind(df_pca, colData)
  
  pca_rotation_named <- pca$rotation
  row.names(pca_rotation_named) <- gene_names
  
  # extract top 5 features per principal component
  topfeat <- list("pc1" = list(), 
                   "pc2" = list(),
                   "pc3" = list(), 
                   "pc4" = list(),
                   "pc5" = list(),
                   "pc6" = list(), 
                   "pc7" = list(),
                   "pc8" = list(), 
                   "pc9" = list(),
                   "pc10" = list())

  prcomps <- c(1:10)
  
  topfeat <- lapply(prcomps, function(x){
    pc <- order(abs(pca$rotation[,x]), decreasing = TRUE)[1:5]
    topfeat[x] <- pc
  })
  
  out <- list("pca_obj"=pca, "df_pca"=df_pca, "evar"=explained_variance_ratio, 
              "pca_topfeat"=topfeat, "pca_rotation"=pca$rotation, "pca_rotation_named"=pca_rotation_named)
  return(out)
}

compute_UMAP <- function(x, colData, rstate){
  umap_obj <- umap(x, method="umap-learn", random_state=rstate, n_neighbors=30)
  
  df_umap <- as.data.frame(umap_obj$layout)
  df_umap <- cbind(df_umap, colData)
  
  out <- list("umap_obj"=umap_obj, "df_umap"=df_umap)
  return(out)
}

compute_UMAP_uwot <- function(x, colData, rstate){
  set.seed <- rstate
  umap_obj <- umap(x, n_neighbors=30, metric="cosine", n_threads=2)
  
  df_umap <- as.data.frame(umap_obj)
  df_umap <- cbind(df_umap, colData)
  
  out <- list("umap_obj"=umap_obj, "df_umap"=df_umap)
  return(out)
}

# functions for colors
my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
my.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)


preProcessTestData <- function(inputfile, preProc_train, dispFunc_train, preProc_train_sc){
  testdata <- parseInput(inputfile)
  print(dim(testdata))
  testdata_df <- t(testdata)
  test_expr <- predict(preProc_train, testdata_df)
  
  test_sc <- predict(preProc_train_sc, test_expr)
  test_log2 <- apply(test_expr, 2, function(x) log2(x+1))
  pp_minmax <- preProcess(test_expr, method="range")
  test_minmax <- predict(pp_minmax, test_expr)

  
  dds_test <- DESeqDataSetFromMatrix(countData = t(test_expr),
                                     colData=as.data.frame(colnames(testdata)),
                                     design=~1) # no design
  
  dds_test <- DESeq(dds_test)
  # apply same dispersion fuction to testdata
  dispersionFunction(dds_test) <- dispFunc_train
  test_vst_s4 <- vst(dds_test, blind=F)
  test_vst <- t(assay(test_vst_s4))
  
  #rm(dds_test, testdata, testdata_df, test_expr)
  
  output=list("test_vst"=test_vst, "test_sc"=test_sc, "test_log2"=test_log2, "test_minmax"=test_minmax)  
  
  return(output)
}

projectTestData <- function(plot_type, preProc_output, scaling_method, header){
  if (scaling_method == "unitvar") {
    testdata <- preProc_output$test_sc
    if(plot_type == "PCA"){
      train_obj <- pca_results$pca_unitvar$pca_obj
    } else if(plot_type == "UMAP"){
      train_obj <- umap1$umap_obj
    }
  }
  if (scaling_method == "vst"){
    testdata <- preProc_output$test_vst
    if(plot_type == "PCA"){
      train_obj <- pca_results$pca_vst$pca_obj
    } else if(plot_type == "UMAP"){
      train_obj <- umap2$umap_obj
    }
  }
  if (scaling_method == "log") {
    testdata <- preProc_output$test_log2
    if(plot_type == "PCA"){
      train_obj <- pca_results$pca_log$pca_obj
    } else if(plot_type == "UMAP"){
      train_obj <- umap3$umap_obj
    }
  }
  if (scaling_method == "minmax") {
    testdata <- preProc_output$test_minmax
    if(plot_type == "PCA"){
      train_obj <- pca_results$pca_minmax$pca_obj
    } else if(plot_type == "UMAP"){
      train_obj <- umap4$umap_obj
    }
  }
  
  
  projection <- predict(train_obj, testdata)
  # artificial metadata frame for testdata
  meta_test <- data.frame(matrix(ncol = length(header), 
                                 nrow = nrow(projection)))
  colnames(meta_test) <- header
  meta_test[is.na(meta_test)] = "not_available"
  meta_test$Age_at_index <- rep(NA, length(row.names(meta_test)))
  meta_test$Survival_time <- rep(NA, length(row.names(meta_test)))
  # meta_test$Tumor_stage <- rep(NA, length(row.names(meta_test))) 
  meta_test$Tumor_stage %>% replace_na('not_available')
  meta_test <- mutate_if(meta_test, is.character, as.factor)
  
  # df for projected testdata with metadta
  if(plot_type == "PCA"){ testdata_projected <- cbind(projection[,1:10], meta_test) }
  else{ testdata_projected <- cbind(projection[,1:2], meta_test) }
  
  return(testdata_projected)
}

pca_rotation_table <- function(scaling_method){
  if (scaling_method == "unitvar") {return(pca_results$pca_unitvar$pca_rotation_named[,1:10])}
  if (scaling_method == "vst") {return(pca_results$pca_vst$pca_rotation_named[,1:10])}
  if (scaling_method == "log") {return(pca_results$pca_log$pca_rotation_named[,1:10])}
  if (scaling_method == "minmax") {return(pca_results$pca_minmax$pca_rotation_named[,1:10])}
}

################################## READ DATA #################################

# load("/Users/susanne/Documents/code/r_projects/anovaget_app/data/lihc_chol_liri_gtex_summarizedExperiment.RData")
# 
# # COUNT DATA
# gene_names <- rowData(sexpr)$X
# expr <- assays(sexpr)$counts
# texpr <- t(expr)
# #old_names <- row.names(texpr)
# gene_ids <- colnames(texpr)
#new_names <- to_vec(for(i in seq(length(row.names(texpr)))) paste0("id", i))
#row.names(texpr) <- new_names
#colnames(expr) <- new_names


# test if same sexpr and liver_expr
#liver_gene_names <- rowData(liver_expr)$X
#identical(gene_names,liver_gene_names)#
#liver_expr_expr <- assays(liver_expr)$counts
#all.equal(expr, liver_expr_expr, check.names=F)#
#liver_texpr <- t(liver_expr_expr)
#all.equal(texpr,liver_texpr,check.names=F)

# old_names == row.names(meta_df)

# PREPROCESSING DATA
# remove near zero variance genes
# pp_nvz <- preProcess(texpr, method = c("nzv")) # 32163 genes
# train_expr <- predict(pp_nvz, texpr)
# 
# # METADATA
# metadata <- colData(sexpr)
# meta_df <- as.data.frame(metadata)
# meta_df$Treatment_Type <- (sub('\"TACE RFA\"', 'TACE+RFA', meta_df$Treatment_Type))
# meta_df[meta_df == ''] <- 'not_available'
# meta_df[meta_df == "LIRI-JP"] <- "ICGC-LIRI-JP"
# meta_df[meta_df == "PRJNA75899"] <- "GTEX-PRJNA75899"
# meta_df <- mutate_if(meta_df, is.character, as.factor)
# row.names(meta_df) <- new_names



#  new Summarized Experiment!
#liver_expr <- SummarizedExperiment(assays=list(counts=expr),
#                              rowData=gene_ids,
#                              colData=meta_df)
#rowData(liver_expr)$X <- gene_names

#save(liver_expr, file = "/Users/susanne/Documents/code/r_projects/anovaget_app/data/lihc_chol_liri_gtex_summarizedExperiment_harmonized.RData")

load("/Users/susanne/Documents/code/r_projects/anovaget_app/data/lihc_chol_liri_gtex_summarizedExperiment_harmonized.RData")
gene_names <- rowData(liver_expr)$X
gene_matching <- rowData(liver_expr)
expr <- assays(liver_expr)$counts
texpr <- t(expr)

load(pseudogenes, file="/Users/susanne/Documents/code/r_projects/anovaget_app/data/pseudogenes.RData")

# remove near zero variance genes
# pp_nvz <- preProcess(texpr, method = c("nzv")) # 32163 genes
# train_expr <- predict(pp_nvz, texpr)

# METADATA
metadata <- colData(liver_expr)
meta_df <- as.data.frame(metadata) 
meta_df <- apply(meta_df, 2, function(x) str_replace_all(x, "Hepatocellular carcinoma", "HCC"))
meta_df <- apply(meta_df, 2, function(x) str_replace_all(x, "hepatocellular carcinoma", "HCC"))
meta_df <- apply(meta_df, 2, function(x) str_replace_all(x, "Cholangiocarcinoma", "CCA"))
meta_df <- apply(meta_df, 2, function(x) str_replace_all(x, "cholangiocarcinoma", "CCA"))
meta_df <- as.data.frame(meta_df) 
row.names(meta_df) <- row.names(metadata)
meta_df$Age_at_index <- as.numeric(meta_df$Age_at_index)
meta_df$Survival_time <- as.numeric(meta_df$Survival_time)
meta_df$Tumor_stage <- meta_df$Tumor_stage %>% replace_na('not_available')
# is.na(meta_df$Tumor_stage) <- "not_available"
meta_df <- mutate_if(meta_df, is.character, as.factor)


# DESEQ OBJECT WITHOUT DESIGN FOR FROZEN VST TRANSFORM
# running DESeq makes the dispersionFunction available for VST transformation
# count data: rows=genes,cols=samples
#dds_train <- DESeqDataSetFromMatrix(countData = t(train_expr), 
#                     colData = meta_df,
#                     design = ~ 1) # no design
#dds_train <- DESeq(dds_train)
#save(dds_train, file = "/Users/susanne/Documents/code/r_projects/anovaget_app/data/lihc_chol_liri_gtex_dds_object_new_ids.RData")

#preproc_output <- list("pp_nvz"=pp_nvz, "pp_sc"=pp_sc, "train_dispersionFunc"=train_dispersionFunc)
#save(preproc_output, file = "/Users/susanne/Documents/code/r_projects/anovaget_app/data/lihc_chol_liri_gtex_preproc.RData")
load( "/Users/susanne/Documents/code/r_projects/anovaget_app/data/lihc_chol_liri_gtex_preproc.RData")
train_expr <- predict(preproc_output$pp_nvz, texpr)

genes_to_remove <- preproc_output$pp_nvz$method$remove
gene_names_reduced <- gene_matching[!rownames(gene_matching) %in% genes_to_remove, ] 


load("/Users/susanne/Documents/code/r_projects/anovaget_app/data/lihc_chol_liri_gtex_dds_object_new_ids.RData")
train_dispersionFunc <- dispersionFunction(dds_train) 
train_vst <- vst(dds_train, blind=T) 

# pp_sc <- preProcess(train_expr, method = c("scale", "center")) # 32163 genes
train_sc <- predict(preproc_output$pp_sc, train_expr)


#####################  APPLY FUNCTIONS ##############################

# apply scaling
#unitvar <- scaling_method(train_expr, "Unit variance")
#log2_scaled <- scaling_method(train_expr, "Log")
#minmax <- scaling_method(train_expr, "MinMax")

#scaled_outputs <- list("unitvar"=unitvar,"train_vst"=train_vst,"log2_scaled"=log2_scaled,"minmax"=minmax)
#save(scaled_outputs, file = "/Users/susanne/Documents/code/r_projects/anovaget_app/data/scaled_outputs.RData")
load( "/Users/susanne/Documents/code/r_projects/anovaget_app/data/scaled_outputs.RData")

#out1 <- compute_PCA(scaled_outputs$unitvar, meta_df, "unitvar", gene_names_reduced)

## make sure to use right input data (samples =rows, genes=cols)
#out2 <- compute_PCA(t(assay(scaled_outputs$train_vst)), meta_df, "vst", gene_names_reduced)

#out3 <- compute_PCA(scaled_outputs$log2_scaled, meta_df, "log", gene_names_reduced)

#out4 <- compute_PCA(scaled_outputs$minmax, meta_df, "minmax", gene_names_reduced)

# free space
rm(liver_expr, expr, texpr, dds_train)
rm(train_expr, train_sc, train_vst)


# pca_results <- list("pca_unitvar"=out1, "pca_vst"=out2, "pca_log"=out3, "pca_minmax"=out4)
#save(pca_results, file = "/Users/susanne/Documents/code/r_projects/anovaget_app/data/pca_results_new_ids_topfeat.RData")

load("/Users/susanne/Documents/code/r_projects/anovaget_app/data/pca_results_new_ids_topfeat.RData")

#umap.uwot <- compute_UMAP_uwot(unitvar, meta_df, 42)
#umap1 <- umap.uwot
# umap_plotly_function("unitvar", "Sample_type")


#umap2 <- compute_UMAP_uwot(t(assay(train_vst)), meta_df, 42)
#umap_plotly_function("vst", "Sample_type")

umap1 <- compute_UMAP(scaled_outputs$unitvar, meta_df, 42)
umap2 <- compute_UMAP(t(assay(scaled_outputs$train_vst)), meta_df, 42)
umap3 <- compute_UMAP(scaled_outputs$log2_scaled, meta_df, 42)
umap4 <- compute_UMAP(scaled_outputs$minmax, meta_df, 42)

#umap_results <- list("umap_unitvar"=umap1, "umap_vst"=umap2, "umap_log"=umap3, "umap_minmax"=umap4)
#save(umap_results, file = "/Users/susanne/Documents/code/r_projects/anovaget_app/data/umap_results_new_ids.RData")

#rm(log2_scaled,  unitvar, minmax)

#load("/Users/susanne/Documents/code/r_projects/anovaget_app/data/umap_results_new_ids.RData")


#umap5 <- umap_results$umap_unitvar
#umap6 <- umap_results$umap_vst
#umap7 <- umap_results$umap_log
#umap8 <- umap_results$umap_minmax

# umap from pca without PC1 and PC2
# umap10 <- compute_UMAP(pca_results$pca_unitvar$pca_obj$x[,3:ncol(pca_results$pca_unitvar$pca_obj$x)], meta_df, 42)


########################## PCA COLORS ##########################################
# use annotations of liver metadata

sample_type_colors <- c(normal="#636EFA", tumor="#EF553B",not_available="#323232")
project_colors <- brewer.pal(length(unique(meta_df$Project))+1, "Set2")
treatment_type <- brewer.pal(length(unique(meta_df$Treatment_Type))+1, "Set1")
treatment_or_therapy <- brewer.pal(length(unique(meta_df$Treatment_or_Therapy))+1, "Set2") #scale_color_brewer(palette="Set1") # (length(unique(metadata$Treatment_or_Therapy)))
primary_diagnosis <- brewer.pal(length(unique(meta_df$Primary_diagnosis))+1, "Set2") # (length(unique(metadata$Primary_diagnosis)))
primary_site <- brewer.pal(length(unique(meta_df$Primary_site))+1, "Set2")# (length(unique(metadata$Primary_diagnosis)))
age_colors <- brewer.pal(n=9, name="Blues")
vital_status_colors <- brewer.pal(length(unique(meta_df$Vital_status))+1, "Set2")
sex_colors <- c(female="#EF553B", male="#636EFA", "#000000")
survival_colors <- brewer.pal(n=9, name="Blues")
# tumor_stage <- brewer.pal(length(unique(meta_df$Tumor_stage)), "Greys")
tumor_stage <- brewer.pal(length(unique(meta_df$Tumor_stage))+1, "Set2")
icd10_colors <- brewer.pal(length(unique(meta_df$Icd10))+1, "Set2")


color_list <- list("Sample_type"= sample_type_colors,
                   "Project"=project_colors, 
                   "Primary_site"=primary_site, 
                   "Treatment_or_Therapy"=treatment_or_therapy,
                   "Treatment_Type"=treatment_type,
                   "Primary_diagnosis"=primary_diagnosis,
                   "Age_at_index"=age_colors,
                   "Vital_status"=vital_status_colors,
                   "Sex"=sex_colors,
                   "Survival_time"=survival_colors,
                   "Tumor_stage"=tumor_stage,
                   "Icd10"=icd10_colors)

pca_plotly_function <- function(pcx, pcy, scaling_method, colorby){
  
  if (scaling_method == "unitvar"){
    df_out <- out1$df_pca
    explained_variance_ratio = out1$evar
  }
  if (scaling_method == "vst"){
    df_out <- out2$df_pca
    explained_variance_ratio = out2$evar
  }
  if (scaling_method == "log"){
    df_out <- out3$df_pca
    explained_variance_ratio = out3$evar
  }

  pal <- color_list[[paste(colorby)]]
  dim1 <- paste0("PC", pcx, " - " , round(explained_variance_ratio[pcx], 2), " %")
  dim2 <- paste0("PC", pcy, " - " , round(explained_variance_ratio[pcy], 2), " %")
  
  tx <- highlight_key(df_out, ~row.names(df_out), "Select a sample")
  
  fig <- plot_ly(type="scatter", mode= "markers", colors=pal) 
  fig <- fig %>% add_trace(data=tx, x = tx$data()[,pcx], y = tx$data()[,pcy], color = tx$data()[,colorby],
              text = row.names(df_out),
              hovertemplate = paste("Sample:", row.names(df_out), 
                                    "\nProject:", df_out$Project, 
                                    "\nSex:", df_out$Sex, 
                                    "\nAge:", df_out$Age_at_index,
                                    "\nPrimary diagnosis:", df_out$Primary_diagnosis,
                                    "\nSurvival time:", df_out$Survival_time,
                                    "\nTumor stage:", df_out$Tumor_stage,
                                    '<extra></extra>')) %>% highlight_key(row.names(df_out)) %>%
                                    highlight(on = "plotly_click", off = "plotly_doubleclick") 

    fig <- fig %>% layout(title= list(text = "", xanchor="right", yanchor="top", pad=list(b=50, t=20)), 
                          list(title=list(text=colorby)),
                           xaxis = list(title = list(text=dim1)),
                           yaxis = list(title = list(text=dim2)))
    return(fig)
}


plotly_plotting_function <- function(plot_type, pcx, pcy, scaling_method, colorby, row_id){
  
  if (scaling_method == "unitvar"){
    if(plot_type == "PCA"){
      df_out <- pca_results$pca_unitvar$df_pca
      explained_variance_ratio = pca_results$pca_unitvar$evar
      topfeat <- pca_results$pca_unitvar$pca_topfeat
      rotation <- pca_results$pca_unitvar$pca_rotation_named
      vector_exp <- 15000
    } else if(plot_type == "UMAP"){
      df_out <- umap1$df_umap # umap_results$umap_unitvar$
    }
  }
  if (scaling_method == "vst"){
    if(plot_type == "PCA"){
      df_out <- pca_results$pca_vst$df_pca
      explained_variance_ratio = pca_results$pca_vst$evar
      topfeat <- pca_results$pca_vst$pca_topfeat
      rotation <- pca_results$pca_vst$pca_rotation_named
      vector_exp <- 2000
    } else if(plot_type == "UMAP"){
      df_out <- umap2$df_umap #umap_results$umap_vst
    }
  }
  if (scaling_method == "log"){
    if(plot_type == "PCA"){
      df_out <- pca_results$pca_log$df_pca
      explained_variance_ratio = pca_results$pca_log$evar
      topfeat <- pca_results$pca_log$pca_topfeat
      rotation <- pca_results$pca_log$pca_rotation_named
      vector_exp <- 10000
    } else if(plot_type == "UMAP"){
      df_out <- umap3$df_umap #umap_results$umap_log
    }
  }
  if (scaling_method == "minmax"){
    if(plot_type == "PCA"){
      df_out <- pca_results$pca_minmax$df_pca
      explained_variance_ratio = pca_results$pca_minmax$evar
      topfeat <- pca_results$pca_minmax$pca_topfeat
      rotation <- pca_results$pca_minmax$pca_rotation_named
      vector_exp <- 1000
    } else if(plot_type == "UMAP"){
      df_out <- umap4$df_umap #umap_results$umap_minmax$
    }
  }
  
  pal <- color_list[[paste(colorby)]]
  tx <- highlight_key(df_out, ~row.names(df_out), "Select a sample")
  
  
  if(plot_type == "PCA"){
    
    dim1 <- paste0("PC", pcx, " - " , round(explained_variance_ratio[pcx], 2), " %")
    dim2 <- paste0("PC", pcy, " - " , round(explained_variance_ratio[pcy], 2), " %")
    
    # add arrows of top 5 featues per principal component
    # split up in 2 annotations per PC to color them differently 
    # attaches names to top features from rotation
    l.x.pcx <-  rotation[,pcx][topfeat[[pcx]]] * vector_exp
    l.y.pcx <- rotation[,pcy][topfeat[[pcx]]] * vector_exp
    
    l.x.pcy <-  rotation[,pcx][topfeat[[pcy]]] * vector_exp
    l.y.pcy <- rotation[,pcy][topfeat[[pcy]]] * vector_exp
    
    # add buttons and annotations
    pc1_annotations <- list( x = l.x.pcx, # arrow.x.end,           ## PC1 arrows are grey, PC2 arrows are blue
                             y = l.y.pcx, # ~arrow.y.end,
                             showarrow = F,
                             borderpad=4,
                             xref = 2, yref = 2,
                             axref = 2, ayref = 2,
                             text = names(l.x.pcx),
                             valign = "top",
                             align = "left",
                             xanchor = "right",
                             yanchor = "bottom",
                             ax = 2, 
                             ay = 2
    ) 
    
    pc1_arrows <- list(x = l.x.pcx, # arrow.x.end,
                       y = l.y.pcx, # ~arrow.y.end,
                       showarrow = TRUE,
                       arrowcolor="#636363",
                       text = names(l.x.pcx),
                       xref = "x", yref = "y",
                       axref = "x", ayref = "y",
                       text = "",
                       ax = 0, 
                       ay = 0)
    
    pc2_annotations <- list(
      x = l.x.pcy, # arrow.x.end,
      y = l.y.pcy, # ~arrow.y.end,
      showarrow = F,
      xref = 2, yref = 2,
      axref = 2, ayref = 2,
      text = names(l.y.pcy),
      valign = "top",
      align = "left",
      xanchor = "right",
      ax = 2, 
      ay = 2
    )
    
    pc2_arrows <- list(x = l.x.pcy, # arrow.x.end,
                       y = l.y.pcy, # ~arrow.y.end,
                       showarrow = TRUE,
                       arrowcolor="#1d33a1",
                       text = names(l.y.pcy),
                       xref = "x", yref = "y",
                       axref = "x", ayref = "y",
                       text = "",
                       ax = 0, 
                       ay = 0)
    
    ##### buttons
    # updatemenus component
    updatemenus <- list(
      list(
        active = -1,
        type = "buttons",
        direction = "right",
        xanchor = 'center',
        yanchor = "top",
        pad = list('r'= 0, 't'= 10, 'b' = 10),
        x = 0.5,
        y = 1.17,
        buttons = list(
          list(
            label = "High",
            method = "update",
            args = list(list(visible = c(FALSE, TRUE)),
                        list(title = "PCA",
                             annotations = list(c(), pc1_annotations)))),
          list(
            label = "Low",
            method = "update",
            args = list(list(visible = c(TRUE, FALSE)),
                        list(title = "PCA",
                             annotations = list(pc2_annotations, c() )))),
          list(
            label = "Both",
            method = "update",
            args = list(list(visible = c(TRUE, TRUE)),
                        list(title = "PCA",
                             annotations = list(pc1_annotations, pc2_annotations)))),
          list(
            label = "Reset",
            method = "update",
            args = list(list(visible = c(TRUE, TRUE)),
                        list(title = "PCA",
                             annotations = list(c(), c())))))
      )
    )
    
    
  #print(pc1_annotations)
  #print(pc1_arrows)
  #print(names(l.x.pcx))
  #print(names(l.y.pcy))
    
    
  }
  else if(plot_type == "UMAP"){
    dim1 <- paste0("UMAP1")
    dim2 <- paste0("UMAP2")
  }
  
  
  sorted_list <- df_out[order(df_out[,colorby]), ]
  h_marker <- ifelse(row.names(sorted_list) %in% row_id, 20, 10)
  h_border <- ifelse(row.names(sorted_list) %in% row_id, 'rgb(0,0,0)', 'rgb(255,255,255)')
  h_opa <- ifelse(row.names(sorted_list) %in% row_id, 1, 0.8)
  h_width <- ifelse(row.names(sorted_list) %in% row_id, 2, 1)
  
  fig <- plot_ly(type="scatter", mode= "markers", colors=pal) 
  fig <- fig %>% add_trace(data=tx, x = tx$data()[,pcx], y = tx$data()[,pcy], color = tx$data()[,colorby],
                           text = row.names(df_out), 
                           marker = list(size=h_marker, 
                                         line = list(color = h_border,
                                                    width = h_width),
                                        opacity = h_opa),
                           hovertemplate = paste("Sample:", row.names(df_out), 
                                                 "\nProject:", df_out$Project, 
                                                 "\nSex:", df_out$Sex, 
                                                 "\nAge:", df_out$Age_at_index,
                                                 "\nPrimary diagnosis:", df_out$Primary_diagnosis,
                                                 "\nSurvival time:", df_out$Survival_time,
                                                 "\nTumor stage:", df_out$Tumor_stage,
                                                 '<extra></extra>')) %>% highlight_key(row.names(df_out)) %>%
          highlight(on = "plotly_click", off = "plotly_doubleclick")  
  if (plot_type == "PCA"){
    # fig <- fig %>% layout(udpatemenus=updatemenus) # , annotations = list(pc1_arrows, pc2_arrows)
     fig <- fig %>% add_annotations( x = l.x.pcx, # arrow.x.end,           ## PC1 arrows are grey, PC2 arrows are blue
                                     y = l.y.pcx, # ~arrow.y.end,
                                     showarrow = F,
                                     borderpad=4,
                                    xref = 2, yref = 2,
                                     axref = 2, ayref = 2,
                                    text = names(l.x.pcx),
                                   valign = "top",
                                     align = "left",
                                     xanchor = "right", 
                                     yanchor = "bottom",
                                     ax = 2, 
                                     ay = 2) %>% 
       add_annotations( x = l.x.pcy, # arrow.x.end,
                        y = l.y.pcy, # ~arrow.y.end,
                        showarrow = F,
                        xref = 2, yref = 2,
                        axref = 2, ayref = 2,
                        text = names(l.y.pcy),
                        valign = "top",
                        align = "left",
                        xanchor = "right",
                        ax = 2, 
                        ay = 2) %>%
       add_annotations( x = l.x.pcx, # arrow.x.end,
                        y = l.y.pcx, # ~arrow.y.end,
                        showarrow = TRUE,
                        arrowcolor="#636363",
                        xref = "x", yref = "y",
                        axref = "x", ayref = "y",
                        text = "",
                        ax = 0, 
                        ay = 0)  %>%
       add_annotations( x = l.x.pcy, # arrow.x.end,
                        y = l.y.pcy, # ~arrow.y.end,
                        showarrow = TRUE,
                        arrowcolor="#1d33a1",
                        xref = "x", yref = "y",
                        axref = "x", ayref = "y",
                        text = "",
                        ax = 0, 
                        ay = 0)  
  }
  fig <- fig %>% layout(title = list(text = plot_type, xanchor="left", x=0.1), 
                        margin=list(l=20, r=20, t=30, b=30, pad=1),
                        legend = list(orientation="h", yanchor="bottom", y=-0.5,
                                      xanchor="right", x=1, entrywidth=500), # title=list(text=colorby), 
                        xaxis = list(title = list(text=dim1)),
                        yaxis = list(title = list(text=dim2)))
  return(fig)
}


TEST <- plotly_plotting_function("UMAP", 1, 2 , "unitvar", "Sample_type", "Sample_type")
plotly_plotting_function("UMAP", 1, 2 , "unitvar", "Primary_diagnosis", "Project")


plotly_plotting_function("PCA", 1, 2 , "unitvar", "Sample_type", "Sample_type")
plotly_plotting_function("PCA", 1, 2 , "vst", "Sample_type", "Sample_type")
plotly_plotting_function("PCA", 1, 2 , "log", "Sample_type", "Sample_type")
plotly_plotting_function("PCA", 1, 2 , "minmax", "Sample_type", "Sample_type")

pca_fig <- plotly_plotting_function("PCA", 1, 2 , "vst", "Sample_type", "Sample_type")

###################################################################################################
# test nearest neighbors and neighbor component analysis


# 
# # test_mat <- as.matrix(ifelse(meta_df$Sample_type %in% "normal", 1, 2))
# 
# 
samples <- row.names(scaled_outputs$unitvar)
nearest3 <- nn2(scaled_outputs$unitvar, umap.out_pp$test_sc, k=3)

nmat <- nearest3$nn.idx
dists <- apply(nearest3$nn.dist, 2, function(x) round(x, 2))

which_ones <- data.frame(nmat)
# colnames(which_ones) <- c("first", "second", "third")
which_ones <- apply(which_ones, 2, function(x) paste0("id", x))
which_ones <- cbind(which_ones, dists)
colnames(which_ones) <- c("first", "second", "third", "dist1", "dist2", "dist3")
row.names(which_ones) <- row.names(umap.out_pp$test_sc)

nearest_pca<- nn2(pca_results$pca_unitvar$df_pca[,1:10], pca_project[,1:10], k=3)
nmat <- nearest_pca$nn.idx
dists <- apply(nearest_pca$nn.dist, 2, function(x) round(x, 2))

which_ones <- data.frame(nmat)
# colnames(which_ones) <- c("first", "second", "third")
which_ones <- apply(which_ones, 2, function(x) paste0("id", x))
which_ones <- cbind(which_ones, dists)
colnames(which_ones) <- c("first", "second", "third", "dist1", "dist2", "dist3")
row.names(which_ones) <- row.names(umap.out_pp$test_sc)


# function to compute nearest neighbors and return dataframe and hoverinfo for plot
compute_NN <- function(scaling_method, test_pca, n){
  
  if (scaling_method == "unitvar") { train_pca <- pca_results$pca_unitvar$df_pca[,1:10] }
  if (scaling_method == "vst") { train_pca <- pca_results$pca_vst$df_pca[,1:10] }
  if (scaling_method == "log") { train_pca <- pca_results$pca_log$df_pca[,1:10] }
  if (scaling_method == "minmax") { train_pca <- pca_results$pca_minmax$df_pca[,1:10] }
  
  nearest_neighbors <- nn2(train_pca, test_pca, k=n)
  # data frame of rounded distances 
  nmat <- nearest_neighbors$nn.idx
  ndists <- apply(nearest_neighbors$nn.dist, 2, function(x) round(x,2))
  # label 
  neighbor_df <- data.frame(apply(nmat, 2, function(x) paste0("id", x)))
  neighbor_df <- cbind(neighbor_df, ndists)
  
  # idx of neighbor
  s <- rep(seq(1:n), 2)
  # n1 n2, n..., dist1, dist2, ...
  colnames(neighbor_df) <-  to_vec(for(i in 1:length(s)) ifelse(i<=n , paste0("n",s[i]), paste0("dist", s[i])))
  row.names(neighbor_df) <- row.names(test_pca)
  
  
  return(neighbor_df)
  
}


pca_unitvar_NN <- compute_NN("unitvar", pca_project[,1:10], 3)

seq_n <- 1:(ncol(pca_unitvar_NN)/2)
num_col <- ncol(pca_unitvar_NN)/2
for (i in seq_n){ print(paste0("\nNearest n ", pca_unitvar_NN[,i], ", dist: ", pca_unitvar_NN[,i+num_col])) }


for (i in 1:length(s)){
  if(nr == nrow(pca_unitvar_NN)){break}
  print(paste0("\nNearest Neighbor ", s[i], " :",  pca_unitvar_NN[,i], ", dist: ", pca_unitvar_NN[nr,i+s]))
  nr <- nr +1
}

###################################################################################################


test_tx <- highlight_key(pca_results$pca_unitvar$df_pca, ~row.names(out1$df_pca), "Select a sample")

umap_plotly_function <- function(scaling_method, colorby){

  if (scaling_method == "none"){
    df_umap <- umap_not_scaled$df_umap
  }
  if (scaling_method == "unitvar"){
    df_umap <- umap1$df_umap
  }
  if (scaling_method == "vst"){
    df_umap <- umap2$df_umap
  }
  if (scaling_method == "log"){
    df_umap <- umap3$df_umap
  }
  if (scaling_method == "minmax"){
    df_umap <- umap4$df_umap
  }
  
  pal <- color_list[[paste(colorby)]]
  tx <- highlight_key(df_umap, ~row.names(df_umap), "Select a sample")
  fig <- plot_ly(type="scatter", mode= "markers", colors=pal) 
  fig <- fig %>% add_trace(data=tx, x = tx$data()[,1], y = tx$data()[,2], color = tx$data()[,colorby],
                           text = row.names(df_umap),
                           hovertemplate = paste("Sample:", row.names(df_umap), 
                                                 "\nProject:", df_umap$Project, 
                                                 "\nSex:", df_umap$Sex, 
                                                 "\nAge:", df_umap$Age_at_index,
                                                 "\nPrimary diagnosis:", df_umap$Primary_diagnosis,
                                                 "\nSurvival time:", df_umap$Survival_time,
                                                 "\nTumor stage:", df_umap$Tumor_stage,
                                                 '<extra></extra>')) %>% highlight_key(row.names(df_umap)) %>%
    highlight(on = "plotly_click", off = "plotly_doubleclick") 
  
  fig <- fig %>% layout(title=list(text = "UMAP", xanchor="left", x=0.1), 
                        list(title=list(text=colorby)),
                        xaxis = list(title = list(text =paste0("UMAP 1"))),
                        yaxis = list(title = list(text =paste0("UMAP2"))))
  return(fig)
}

plotly_add_trace <- function(fig, testdata, pcx, pcy, colorby, neighbors){
  
  testdata <- cbind(testdata, neighbors)
  
  tx2 <- highlight_key(testdata, ~row.names(testdata))
  
  fig <- fig %>%
    add_trace(data=tx2, x = tx2$data()[,pcx], y = tx2$data()[,pcy], color = tx2$data()[,colorby],
              hovertemplate = paste("Sample:", row.names(testdata),
                                    "\nNearest neighbor 1:", testdata$first, " , dist:", testdata$dist1,
                                    "\nNearest neighbor 2:", testdata$second, " , dist:", testdata$dist2,
                                    "\nNearest neighbor 3:", testdata$third, " , dist:", testdata$dist3,
                                    '<extra></extra>'),
              marker=list(color="#323232")
    ) %>%
    highlight(on = "plotly_click", off = "plotly_doubleclick")
  fig
}

################ DUMMY ############################

# dummy <-  read.table("/Users/susanne/Documents/code/r_projects/anovaget_app/data/dummy_liver_featureCounts.tsv", sep="\t", header=TRUE)
dummy <-  read.table("//Users/susanne/Downloads/merged_gene_counts.txt", sep="\t", header=TRUE)
# pseudogenes <- setdiff(dummy$Geneid, row.names(expr))
# save(pseudogenes, file="/Users/susanne/Documents/code/r_projects/anovaget_app/data/pseudogenes.RData")
test <- parseInput(dummy)

umap.out_pp <- preProcessTestData(dummy, preproc_output$pp_nvz, preproc_output$train_dispersionFunc, preproc_output$pp_sc)
umap.out_project <- projectTestData("UMAP", umap.out_pp, "vst", colnames(meta_df))
umap.out.project_unitvar <- projectTestData("UMAP", umap.out_pp, "unitvar", colnames(meta_df))
umap.out_project[is.na(umap.out_project)] <- "not_available"
umap.out.project_unitvar[is.na(umap.out.project_unitvar)] <- "not_available"
umap_fig <- umap_plotly_function("unitvar", "Sample_type")


pca_project <- projectTestData("PCA", umap.out_pp, "unitvar", colnames(meta_df))
pca_project[is.na(pca_project)] <- "not_available"
# 
# plotly_add_trace(TEST, umap.out.project_unitvar, 1, 2, "Sample_type", which_ones)
test_fig <- plotly_add_trace(pca_fig, pca_project, 1, 2, "Sample_type", which_ones)
# sel <- pca_results$pca_unitvar$pca_obj$x[:,1:2]
# test_fig <- add_trace()
# 
# 
# umap_not_scaled <- compute_UMAP(train_expr, meta_df, 42)
# umap_plotly_function("none", "Sample_type")
# umap_plotly_function("none", "Project")



# Example of UI with fluidPage

ui <- fluidPage(
  useShinyjs(), 
  add_busy_spinner(spin = "fading-circle"),
  # Application title
  titlePanel("Anovaget transcriptomics"),
  br(),
  sidebarLayout(
    # Sidebar with selectable inputs
    sidebarPanel(width = 2,
      h3("Visualization options:"),
      selectInput("meta", "Color by:",
                         names(metadata)[5:length(names(metadata))],
                         selected="Sample_type"), 
      br(),
      selectInput("scaling", "Choose a method for scaling:",
                    c("Unit variance"="unitvar", "VST"="vst", "Log"="log", "MinMax"="minmax"),  
                    selected="Unit variance"),
      wellPanel(h5("PCA options:"),
      numericInput("pcx", "Principal component on x-axis:", 1, min=1, max=10, step=1),
      numericInput("pcy", "Principal component on y-axis:", 2, min=1, max=10, step=1)
      ),
      # Horizontal line ----
      tags$hr(),
      h3("Upload user data:"),
      fileInput("upload", NULL, multiple=F, width="100%",
                accept = c(".tsv"), buttonLabel="TSV"),
      numericInput("kn", "Number of nearest neighbors to compute:", 1, min=1, max=10, step=1),
      br(),
      actionButton("Add", HTML("Click to add <br/> user data in plots")),
      bsTooltip(id = "Add", title = "Please input a  raw gene-sample expression matrix", 
                placement = "bottom", trigger = "hover"),
      # Horizontal line ----
      tags$hr(),
      h3("Table options:"),
      collapsibleAwesomeCheckboxGroupInput("show_vars", "Select variables to show:\n", 3,
                         names(metadata)[3:length(names(metadata))], 
                         selected = names(metadata)[3:length(names(metadata))])
    ),
      
    # Show a plot of the generated distribution
    mainPanel(width=10,
              fluidRow(
                column(6, plotlyOutput("pcaPlot")),
                column(6, plotlyOutput("umapPlot"))
              ),
              br(),
              tags$hr(),
              tabsetPanel(type="tabs",
                          tabPanel("Data", 
                                   br(),
                                   span(textOutput(outputId  = "tab1_text"), style="font-size: 20px; font-style: bold;"),
                                   br(),
                                   wellPanel(div(style = 'overflow-x: scroll', DT::dataTableOutput("trainingMetadata"), 
                                        style = "z-index: 10; left:0; right:0; overflow-y:hidden; overflow-xy:auto"))),
                          tabPanel("PCA Loadings", 
                                   br(),
                                   span(textOutput(outputId  = "tab2_text"), style="font-size: 20px; font-style: bold;"),
                                   br(),
                                   wellPanel(div(style = 'overflow-x: scroll', DT::dataTableOutput("pcaLoadings"), 
                                                 style = "z-index: 10; left:0; right:0; overflow-y:hidden; overflow-xy:auto"))),
                          tabPanel("Neighbor Data", 
                                   br(),
                                   span(textOutput(outputId  = "tab3_text"), style="font-size: 20px; font-style: bold;"),
                                   br(),
                                   wellPanel(div(style = 'overflow-x: scroll', DT::dataTableOutput("neighbortable"), 
                                                                  style = "z-index: 10; left:0; right:0; overflow-y:hidden; overflow-xy:auto"))))
              
              
)))


  
  # Server logic
server <- function(input, output, session) {
  
  
  # reactive variables
  metaselect <- reactive(input$meta)
  scale_method <- reactive(input$scaling)
  pcx <- reactive(input$pcx)
  pcy <- reactive(input$pcy)
  kn <- reactive(input$kn)
  
  
  # TAB METADATA
  metadata_table = meta_df[1:nrow(meta_df),3:ncol(meta_df)]
  output$trainingMetadata <- DT::renderDataTable({
    DT::datatable(metadata_table[, input$show_vars, drop = FALSE], 
                  fillContainer = FALSE,
                  options = list(crollX = TRUE,
                                 autoWidth=TRUE, 
                                 columnDefs = list(
                                   list(orderSequence = c("desc", "asc"), targets = "_all"),
                                   list(className = "dt-center", targets = "_all")
                                 ),
                                 processing = FALSE,
                                 pageLength = 5,
                                 lengthMenu = list(c(5, 10, 25, 50, -1), c("5", "10", "25", "50", "All"))
                  )
      )  
      })
  
  output$tab1_text <- renderText({"Metadata annotations from TCGA, ICGC and GTEx."})
  
  
  output$pcaLoadings <- DT::renderDataTable({ DT::datatable(pca_rotation_table(scale_method()),
                                                            options = list(
                                                              "columnDefs" = list(
                                                                list(
                                                                  "targets" = "_all",
                                                                  "render"  = JS(js)
                                                                )
                                                              )
                                                            ))
                                                })
  
  output$tab2_text <- renderText({"Loadings of first 10 Principal Components."})
  
  output$tab3_text <- renderText({"Upload a dataset to compute approximate nearest neighbors to TCGA, ICGC and GTEx data. Nearest Neighbors computed on first 10 Principal Components."})
  
  id_sel <- reactiveVal()
  
  observe({
    idx <- input$trainingMetadata_rows_selected
    id_sel <- row.names(meta_df)[idx]
    id_sel_old_new <- c(id_sel(), id_sel)
    id_sel(unique(id_sel_old_new))
    
  })
  
  
  # clear the set of cars when a double-click occurs
  observeEvent(event_data("plotly_doubleclick"), {
    id_sel(NULL)
  })
  
      
    output$pcaPlot <- renderPlotly({
        
      plotly_plotting_function(
            "PCA",
            pcx(),
            pcy(),
            scaling_method=scale_method(),
            colorby=metaselect(),
            row_id=id_sel()
          )
      })
      
    output$umapPlot <- renderPlotly({
        
        umap_fig <- plotly_plotting_function(
          "UMAP",
          1,
          2,
          scaling_method=scale_method(),
          colorby=metaselect(),
          row_id=id_sel()
        )
      })

     ext_data <- reactive({
      req(input$upload)
      ext <- tools::file_ext(input$upload$name)
      validate(need(ext == "txt", "Invalid file. Please upload a tab-separated .txt file"))
      df <- read.table(input$upload$datapath, sep = "\t", header=T)
      
      out_pp <- preProcessTestData(df, preproc_output$pp_nvz, preproc_output$train_dispersionFunc, preproc_output$pp_sc)
      return(out_pp)
    })
     
     observeEvent(input$Add,{
      output$tab3_text <- renderText({"Nearest neighbors are being computed."})
      out = ext_data()
      pca_df <- projectTestData("PCA", out, scale_method(), colnames(meta_df))
      pca_df[is.na(pca_df)] <- "not_available"
      
      neighbor_df <- compute_NN(scale_method(), pca_df[,1:10], kn())
      # render neighbor dataframe
      output$neighbortable <-  DT::renderDataTable({ DT::datatable(neighbor_df) })
      num_col <- ncol(neighbor_df)/2
      
      
      umap_df <- projectTestData("UMAP", out, scale_method(), colnames(meta_df))
      umap_df[is.na(umap_df)] <- "not_available"

      tx_pca <- highlight_key(pca_df, ~row.names(pca_df)) 
      plotlyProxy("pcaPlot", session) %>%
        plotlyProxyInvoke("addTraces", x=tx_pca$data()[,pcx()], y=tx_pca$data()[,pcy()], color=tx_pca$data()[,metaselect()],
                          type="scatter", mode="markers", name=tx_pca$data()[,metaselect()][1],
                          hovertemplate = paste("Sample:", row.names(pca_df), '<br>',
                                                "Nearest Neighbor ", neighbor_df[,1], ", dist: ", neighbor_df[,1+num_col],
                                                '<extra></extra>'),
                          marker=list(color="#323232")) %>%
       highlight(on = "plotly_click", off = "plotly_doubleclick")
      
      
      tx_umap <- highlight_key(umap_df, ~row.names(umap_df)) 
      plotlyProxy("umapPlot", session) %>%
        plotlyProxyInvoke("addTraces", x=tx_umap$data()[,1], y=tx_umap$data()[,2], color=tx_umap$data()[,metaselect()],
                          type="scatter", mode="markers", name=tx_umap$data()[,metaselect()][1],
                          hovertemplate = paste("Sample:", row.names(umap_df),'<extra></extra>'),
                          marker=list(color="#323232"))  %>%
        highlight(on = "plotly_click", off = "plotly_doubleclick")
      
     })
     
     
}

    
  
  # Complete app with UI and server components
shinyApp(ui, server)

################## TEST
library(ggplot2)
library(ggfortify)
library(plotly)

#p <- autoplot(pca_results$pca_unitvar$pca_obj, data = scaled_outputs$unitvar, 
              #loadings = TRUE, loadings.colour = 'blue',
              #loadings.label = TRUE, loadings.label.size = 3)

#ggplotly(p)

  