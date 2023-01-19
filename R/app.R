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
library(umap)
library(reticulate)
print("starting script")

# shinyOptions(maxUploadSize=10000)
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
  gene_ids <- df$GeneID 
  expr_df <- df[1:nrow(df), 3:ncol(df)]
  row.names(expr_df) <- gene_ids
  sample_names <- formatSamples(expr_df)
  colnames(expr_df) <- sample_names
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
  return(scaled)
}


compute_PCA <- function(x, colData){
  pca <- prcomp(x, center=T, retx=T, scale=T)
  explained_variance_ratio <- summary(pca)[["importance"]]['Proportion of Variance',]
  explained_variance_ratio <- 100 * explained_variance_ratio
  
  # attach metadata
  df_pca <- as.data.frame(pca$x)[1:10]
  df_pca <- cbind(df_pca, colData)
  
  out <- list("pca_obj"=pca, "df_pca"=df_pca, "evar"=explained_variance_ratio, 
              "pca_center"=pca$center, "pca_scale"=pca$scale, "pca_rotation"=pca$rotation)
  return(out)
}

compute_UMAP <- function(x, colData, rstate){
  umap_obj <- umap(x, method="umap-learn", random_state=rstate, n_neighbors=30) #, metric="mahalanobis"
  
  df_umap <- as.data.frame(umap_obj$layout)
  df_umap <- cbind(df_umap, colData)
  
  out <- list("umap_obj"=umap_obj, "df_umap"=df_umap)
  return(out)
}

# functions for colors
my.max <- function(x) ifelse( !all(is.na(x)), max(x, na.rm=T), NA)
my.min <- function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA)


preProcessTestData <- function(inputfile, preProc_train, dispFunc_train, preProc_train_sc){
  testdata <- parseInput(inputfile)
  testdata_df <- t(testdata)
  test_expr <- predict(preProc_train, testdata_df)
  
  test_sc <- predict(preProc_train_sc, test_expr)
  test_log2 <- apply(test_expr, 2, function(x) log2(x+1))
  
  dds_test <- DESeqDataSetFromMatrix(countData = t(test_expr),
                                      colData=as.data.frame(colnames(testdata)),
                                      design=~1) # no design
  
  dds_test <- DESeq(dds_test)
  # apply same dispersion fuction to testdata
  dispersionFunction(dds_test) <- dispFunc_train
  test_vst_s4 <- vst(dds_test, blind=F)
  test_vst <- t(assay(test_vst_s4))
  
  #rm(dds_test, testdata, testdata_df, test_expr)
  
  output=list("test_vst"=test_vst, "test_sc"=test_sc, "test_log2"=test_log2)
  
  return(output)
}

projectTestData <- function(plot_type, preProc_output, scaling_method, header){
  if (scaling_method == "unitvar") {
    testdata <- preProc_output$test_sc
    if(plot_type == "PCA"){
      train_obj <- out1$pca_obj
    } else if(plot_type == "UMAP"){
      train_obj <- umap1$umap_obj
    }
  }
  if (scaling_method == "vst"){
    testdata <- preProc_output$test_vst
    if(plot_type == "PCA"){
      train_obj <- out2$pca_obj
    } else if(plot_type == "UMAP"){
      train_obj <- umap2$umap_obj
    }
  }
  if (scaling_method == "log") {
    testdata <- preProc_output$test_log2
    if(plot_type == "PCA"){
      train_obj <- out3$pca_obj
    } else if(plot_type == "UMAP"){
      train_obj <- umap3$umap_obj
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
  meta_test$Tumor_stage <- rep(NA, length(row.names(meta_test)))
  meta_test <- mutate_if(meta_test, is.character, as.factor)
  
  # df for projected testdata with metadata
  if(plot_type == "PCA"){ testdata_projected <- cbind(projection[,1:10], meta_test) }
  else{ testdata_projected <- cbind(projection[,1:2], meta_test) }
  
  return(testdata_projected)
}

################################## READ DATA #################################

load("/Users/susanne/Documents/repos/forks/qbic-projects/anovaget-shiny/data/lihc_chol_liri_gtex_summarizedExperiment.RData")
print("loaded summarized expr")

# COUNT DATA
gene_names <- rowData(sexpr)$X
expr <- assays(sexpr)$counts
texpr <- t(expr)

# PREPROCESSING DATA
# remove near zero variance genes
pp_nvz <- preProcess(texpr, method = c("nzv")) # 32163 genes
train_expr <- predict(pp_nvz, texpr)

# METADATA
metadata <- colData(sexpr)
meta_df <- as.data.frame(metadata)
meta_df$Treatment_Type <- (sub('\"TACE RFA\"', 'TACE+RFA', meta_df$Treatment_Type))
meta_df[meta_df == ''] <- 'not_available'
meta_df <- mutate_if(meta_df, is.character, as.factor)

# DESEQ OBJECT WITHOUT DESIGN FOR FROZEN VST TRANSFORM
# running DESeq makes the dispersionFunction available for VST transformation
# count data: rows=genes,cols=samples
#dds_train <- DESeqDataSetFromMatrix(countData = t(train_expr), 
                    #colData = meta_df,
                    #design = ~ 1) # no design
#dds_train <- DESeq(dds_train)
# save(dds_train, file = "/Users/susanne/Documents/code/r_projects/anovaget_app/data/lihc_chol_liri_gtex_dds_object.RData")


load("/Users/susanne/Documents/repos/forks/qbic-projects/anovaget-shiny/data/lihc_chol_liri_gtex_dds_object.RData")
print("loaded deseq obj")
train_dispersionFunc <- dispersionFunction(dds_train) 
train_vst <- vst(dds_train, blind=T) 

pp_sc <- preProcess(train_expr, method = c("scale", "center")) # 32163 genes
train_sc <- predict(pp_sc, train_expr)

#####################  APPLY FUNCTIONS ##############################

# apply scaling
print("scaling matrices")
unitvar <- scaling_method(train_expr, "Unit variance")
log2_scaled <- scaling_method(train_expr, "Log")

print("computing pcas")
out1 <- compute_PCA(unitvar, meta_df)
# make sure to use right input data (samples =rows, genes=cols)
out2 <- compute_PCA(t(assay(train_vst)), meta_df)
out3 <- compute_PCA(log2_scaled, meta_df)


print("computing umaps")
umap1 <- compute_UMAP(unitvar, meta_df, 42)
print("umap1 done")
umap2 <- compute_UMAP(t(assay(train_vst)), meta_df, 42)
umap3 <- compute_UMAP(log2_scaled, meta_df, 42)
print("umap done")

# free space
rm(sexpr, expr, texpr, log2_scaled,  unitvar, log2_scaled)



########################## PCA COLORS ##########################################
# use annotations of liver metadata

sample_type_colors <- c(normal="#636EFA", tumor="#EF553B",not_available="#323232")
project_colors <- brewer.pal(length(unique(meta_df$Project))+1, "Set2")
treatment_type <- brewer.pal(length(unique(meta_df$Treatment_Type))+1, "Set1")
treatment_or_therapy <- brewer.pal(length(unique(meta_df$Treatment_or_Therapy))+1, "Set1") 
primary_diagnosis <- brewer.pal(length(unique(meta_df$Primary_diagnosis))+1, "Set1") 
primary_site <- brewer.pal(length(unique(meta_df$Primary_site))+1, "Set1")
age_colors <- brewer.pal(n=9, name="Blues")
vital_status_colors <- brewer.pal(length(unique(meta_df$Vital_status))+1, "Set2")
sex_colors <- c(female="#EF553B", male="#636EFA", "#000000")
survival_colors <- brewer.pal(n=9, name="Greys")
tumor_stage <- brewer.pal(length(unique(meta_df$Tumor_stage)), "Greys")
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
print(color_list$Sample_type)

plotly_plotting_function <- function(plot_type, pcx, pcy, scaling_method, colorby){
  
  if (scaling_method == "unitvar"){
    if(plot_type == "PCA"){
      df_out <- out1$df_pca
      explained_variance_ratio = out1$evar
    } else if(plot_type == "UMAP"){
      df_out <- umap1$df_umap
    }
  }
  if (scaling_method == "vst"){
    if(plot_type == "PCA"){
      df_out <- out2$df_pca
      explained_variance_ratio = out2$evar
    } else if(plot_type == "UMAP"){
      df_out <- umap2$df_umap
    }
  }
  if (scaling_method == "log"){
    if(plot_type == "PCA"){
      df_out <- out3$df_pca
      explained_variance_ratio = out3$evar
    } else if(plot_type == "UMAP"){
      df_out <- umap3$df_umap
    }
    
  }
  
  pal <- color_list[[paste(colorby)]]
    
  tx <- highlight_key(df_out, ~row.names(df_out), "Select a sample")
  
  if(plot_type == "PCA"){
    dim1 <- paste0("PC", pcx, " - " , round(explained_variance_ratio[pcx], 2), " %")
    dim2 <- paste0("PC", pcy, " - " , round(explained_variance_ratio[pcy], 2), " %")
  }
  else if(plot_type == "UMAP"){
    dim1 <- paste0("UMAP1")
    dim2 <- paste0("UMAP2")
  }
  
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
  
  fig <- fig %>% layout(title= list(text = plot_type, xanchor="left", x=0.1), 
                        list(title=list(text=colorby)),
                        xaxis = list(title = list(text=dim1)),
                        yaxis = list(title = list(text=dim2)))
  return(fig)
}

plotly_add_trace <- function(fig, testdata, pcx, pcy, colorby){
  
  tx2 <- highlight_key(testdata, ~row.names(testdata))
  
  fig <- fig %>%
    add_trace(data=tx2, x = tx2$data()[,pcx], y = tx2$data()[,pcy], color = tx2$data()[,colorby],
              hovertemplate = paste("Sample:", row.names(testdata),'<extra></extra>'),
              marker=list(color="#323232")
    ) %>%
    highlight(on = "plotly_click", off = "plotly_doubleclick")
  fig
}
print("loaded function to plot")

#dummy <-  read.table("/Users/susanne/Documents/code/r_projects/anovaget_app/data/dummy_liver_featureCounts.tsv", sep="\t", header=TRUE)
#umap.out_pp <- preProcessTestData(dummy, pp_nvz, train_dispersionFunc, pp_sc)
#umap.out_project <- projectTestData("UMAP", umap.out_pp, "unitvar", colnames(meta_df))
#umap.out_project[is.na(umap.out_project)] <- "not_available"
#umap_fig <- umap_plotly_function("unitvar", "Sample_type")

#plotly_add_trace(umap_fig, umap.out_project, 1, 2, "Sample_type")
#umap_not_scaled <- compute_UMAP(train_expr, meta_df, 42)
#umap_plotly_function("none", "Sample_type")
#umap_plotly_function("none", "Project")



# Example of UI with fluidPage
ui <- fluidPage(
    
  # Application title
  titlePanel("Anovaget transcriptomics"),
    
  sidebarLayout(
      
    # Sidebar with a slider input
    sidebarPanel(
      br(),
      selectInput("meta", "Color by:",
                          names(metadata)[5:length(names(metadata))],
                          selected="Sample_type"), 
      br(),
      selectInput("scaling", "Choose a method for scaling:",
                    c("Unit variance"="unitvar", "VST"="vst", "Log"="log"),
                    selected="Unit variance"),
      wellPanel(h5("PCA options:"),
      numericInput("pcx", "Principal component on x-axis:", 1, min=1, max=10, step=1),
      numericInput("pcy", "Principal omponent on y-axis:", 2, min=1, max=10, step=1)
      ),
      # Horizontal line ----
      tags$hr(),
      fileInput("upload", NULL, multiple=F, width="100%",
                accept = c(".tsv"), buttonLabel="TSV"),
      actionButton("Add", "Add user data"),
      # Horizontal line ----
      tags$hr()
    ),
      
    # Show a plot of the generated distribution
    mainPanel(
      plotlyOutput("pcaPlot"),
      tags$hr(),
      plotlyOutput("umapPlot"),
      tableOutput("uploadFile")
    )
  )
)
  
  # Server logic
server <- function(input, output, session) {

      metaselect <- reactive(input$meta)
      scale_method <- reactive(input$scaling)
      pcx <- reactive(input$pcx)
      pcy <- reactive(input$pcy)
      
      output$pcaPlot <- renderPlotly({
        
        pca_fig <- plotly_plotting_function(
            "PCA",
            pcx(),
            pcy(),
            scaling_method=scale_method(),
            colorby=metaselect()
          )
      })
      
      output$umapPlot <- renderPlotly({
        
        umap_fig <- plotly_plotting_function(
          "UMAP",
          1,
          2,
          scaling_method=scale_method(),
          colorby=metaselect()
        )
      })

      ext_data <- reactive({
      req(input$upload)
      ext <- tools::file_ext(input$upload$name)
      validate(need(ext == "tsv", "Invalid file. Please upload a .tsv file"))
      df <- read.table(input$upload$datapath, sep = "\t", header=T)
      
      out_pp <- preProcessTestData(df, pp_nvz, train_dispersionFunc, pp_sc)
      return(out_pp)
    })

      observeEvent(input$Add,{
      out = ext_data()
      pca_df <- projectTestData("PCA", out, scale_method(), colnames(meta_df))
      pca_df[is.na(pca_df)] <- "not_available"
      
      umap_df <- projectTestData("UMAP", out, scale_method(), colnames(meta_df))
      umap_df[is.na(umap_df)] <- "not_available"

      tx_pca <- highlight_key(pca_df, ~row.names(pca_df)) 
      plotlyProxy("pcaPlot", session) %>%
        plotlyProxyInvoke("addTraces", x=tx_pca$data()[,pcx()], y=tx_pca$data()[,pcy()], color=tx_pca$data()[,metaselect()],
                          type="scatter", mode="markers", name=tx_pca$data()[,metaselect()][1],
                          hovertemplate = paste("Sample:", row.names(pca_df),'<extra></extra>'),
                          marker=list(color="#323232")) %>%
        highlight(on = "plotly_click", off = "plotly_doubleclick")
      
      
      tx_umap <- highlight_key(umap_df, ~row.names(umap_df)) 
      plotlyProxy("umapPlot", session) %>%
        plotlyProxyInvoke("addTraces", x=tx_umap$data()[,1], y=tx_umap$data()[,2], color=tx_umap$data()[,metaselect()],
                          type="scatter", mode="markers", name=tx_umap$data()[,metaselect()][1],
                          hovertemplate = paste("Sample:", row.names(umap_df),'<extra></extra>'),
                          marker=list(color="#323232")) %>%
        highlight(on = "plotly_click", off = "plotly_doubleclick")
      
      
      })
}

    
  
  # Complete app with UI and server components
shinyApp(ui, server)
