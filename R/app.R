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

  # outputs the pca object, a dataframe with the first 10 PCs, ...
  
  out <- list("pca_obj"=pca, "df_pca"=df_pca, "evar"=explained_variance_ratio, 
              "pca_center"=pca$center, "pca_scale"=pca$scale, "pca_rotation"=pca$rotation)
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
  test_vst <- vst(dds_test, blind=F)
  
  rm(dds_test, testdata, testdata_df, test_expr)
  
  # transformed (normalized, scaled) user testdata
  output=list("test_vst"=test_vst, "test_sc"=test_sc, "test_log2"=test_log2)
  
  return(output)
}

projectTestData <- function(preProc_output, scaling_method, header){
  if (scaling_method == "unitvar") {
    train_pca <- out1$pca_obj
    testdata <- preProc_output$test_sc
  }
  if (scaling_method == "vst"){
    train_pca <- out2$pca_obj
    testdata <- t(preProc_output$test_vst)
    }
  if (scaling_method == "log") {
    train_pca <- out3$pca_obj
    testdata <- preProc_output$test_log2
    }

  test_pca_projection <- predict(train_pca, testdata)
  # artificial metadata frame for testdata
  meta_test <- data.frame(matrix(ncol = length(header), 
                                  nrow = nrow(test_pca_projection)))
  colnames(meta_test) <- header
  
  meta_test[is.na(meta_test)] = "not_available"
  meta_test$Age_at_index <- rep(NA, length(row.names(meta_test)))
  meta_test$Survival_time <- rep(NA, length(row.names(meta_test)))
  meta_test$Tumor_stage <- rep(NA, length(row.names(meta_test)))
  
  meta_test <- mutate_if(meta_test, is.character, as.factor)
  
  # pca object with first 10 PCs for testdata
  pca_test <- cbind(test_pca_projection[,1:10], meta_test)
  
  return(pca_test)
}

################################## READ DATA #################################

# contains raw featureCounts from nf-core/rnaseq + metadata
load("data/lihc_chol_liri_gtex_summarizedExperiment.RData")

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
dds_train <- DESeqDataSetFromMatrix(countData = t(train_expr), 
                    colData = meta_df,
                    design = ~ 1) # no design
dds_train <- DESeq(dds_train)
train_dispersionFunc <- dispersionFunction(dds_train) 
train_vst <- vst(dds_train, blind=T) 

pp_sc <- preProcess(train_expr, method = c("scale", "center")) 
train_sc <- predict(pp_sc, train_expr)

#####################  APPLY FUNCTIONS ##############################

# apply scaling
unitvar <- scaling_method(train_expr, "Unit variance")
log2_scaled <- scaling_method(train_expr, "Log")

# compute PCA for all transformed data
out1 <- compute_PCA(unitvar, meta_df)
# make sure to use right input data (samples =rows, genes=cols)
out2 <- compute_PCA(t(assay(train_vst)), meta_df)
out3 <- compute_PCA(log2_scaled, meta_df)

# free space
rm(sexpr, expr, texpr, log2_scaled,  unitvar)

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

################################### PLOTTING FUNCTION ##############################################

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
  pc1 <- paste0("PC",pcx)
  pc2 <- paste0("PC",pcy)
  
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
                            xaxis = list(title = list(text =paste0("PC", pcx, " - " , round(explained_variance_ratio[pcx], 2), " %"))),
                            yaxis = list(title = list(text =paste0("PC", pcy, " - " , round(explained_variance_ratio[pcy], 2), " %"))))
    return(fig)
}


# only running in interactive mode and not in shiny! replaced by plotlyProxy
plotly_add_trace <- function(fig, pca_testdata, pcx, pcy, colorby){
  
  tx2 <- highlight_key(pca_testdata, ~row.names(pca_testdata))
  
  fig <- fig %>%
  add_trace(data=tx2, x = tx2$data()[,pcx], y = tx2$data()[,pcy], color = tx2$data()[,colorby],
            hovertemplate = paste("Sample:", row.names(pca_testdata),'<extra></extra>'),
            marker=list(color="#323232")
  ) %>%
    highlight(on = "plotly_click", off = "plotly_doubleclick")
  fig
}


############################################ UI ####################################################

# Example of UI with fluidPage
ui <- fluidPage(
    
  # Application title
  titlePanel("Anovaget transcriptomics"),
    
  sidebarLayout(
      
    # Sidebar with a dropdown menus
    sidebarPanel(
      selectInput("plot_type", "Plot type", 
                  c("PCA"="pca", "UMAP"="umap")),
      br(),
      selectInput("meta", "Color by:",
                          names(metadata)[5:length(names(metadata))],
                          selected="Sample_type"), 
      br(),
      selectInput("scaling", "Choose a method for scaling",
                    c("Unit variance"="unitvar", "VST"="vst", "Log"="log"),
                    selected="Unit variance"),
      numericInput("pcx", "Principal Component x-axis:", 1, min=1, max=10, step=1),
      numericInput("pcy", "Principal Component y-axis:", 2, min=1, max=10, step=1),
      # Horizontal line ----
      tags$hr(),
      fileInput("upload", NULL, multiple=F, width="100%",
                accept = c(".tsv"), buttonLabel="TSV"),
      actionButton("Add", "Add user data"),
      # Horizontal line ----
      tags$hr()
    ),
      
    # Show a PCA plot
    mainPanel(
      plotlyOutput("pcaPlot"),
      br(),
      # optionally add preview of user data in table format
      tableOutput("uploadFile")
    )
  )
)

########################################## SERVER FUNCTION #########################################

# Server logic
server <- function(input, output, session) {    

      metaselect <- reactive(input$meta)
      scale_method <- reactive(input$scaling)
      pcx <- reactive(input$pcx)
      pcy <- reactive(input$pcy)

      # pca plot training data
      output$pcaPlot <- renderPlotly({
        
        fig <- pca_plotly_function(
            pcx(),
            pcy(),
            scaling_method=scale_method(),
            colorby=metaselect()
          )
      })

      # user data processed upon upload
      ext_data <- reactive({
      req(input$upload)
      ext <- tools::file_ext(input$upload$name)
      validate(need(ext == "tsv", "Invalid file. Please upload a .tsv file"))
      df <- read.table(input$upload$datapath, sep = "\t", header=T)
      
      out_pp <- preProcessTestData(df, pp_nvz, train_dispersionFunc, pp_sc)
      out_project <- projectTestData(out_pp, scale_method(), colnames(meta_df))
      out_project[is.na(out_project)] <- "not_available"
      return(out_project)
    })

      # if user data available, update pca plot
      observeEvent(input$Add,{
        df = ext_data()

        tx2 <- highlight_key(df, ~row.names(df)) 
        plotlyProxy("pcaPlot", session) %>%
          plotlyProxyInvoke("addTraces", x=tx2$data()[,pcx()], y=tx2$data()[,pcy()], color=metaselect(),
                          type="scatter", mode="markers",
                          hovertemplate = paste("Sample:", row.names(df),'<extra></extra>'),
                          marker=list(color="#323232")) %>%
        highlight(on = "plotly_click", off = "plotly_doubleclick")
      })
}



########################################### RUN APP #################################################

# Complete app with UI and server components
shinyApp(ui, server)
