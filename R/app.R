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
library(DT)
library(shinyjs)
library(shinybusy)
library(shinyBS)
library(shinyWidgets)
library(RANN)
library(comprehenr)
library(here)


here()


options(shiny.maxRequestSize = 30*1024^2)
options(warn=-1)

##########################

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
    if(method == "MinMax"){
    pp <- preProcess(mat, method="range")
    scaled <- predict(pp, mat)
  }
  return(scaled)
}


compute_PCA <- function(x, colData, scaling_method){
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

preProcessTestData <- function(inputfile, preProc_train, dispFunc_train, preProc_train_sc){
  testdata <- parseInput(inputfile)
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
      train_obj <-  pca_results$pca_unitvar$pca_obj
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
  meta_test$Tumor_stage %>% replace_na('not_available')
  meta_test <- mutate_if(meta_test, is.character, as.factor)
  
  # df for projected testdata with metadata
  if(plot_type == "PCA"){ testdata_projected <- cbind(projection[,1:10], meta_test) }
  else{ testdata_projected <- cbind(projection[,1:2], meta_test) }
  
  return(testdata_projected)
}


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


################################## READ DATA #################################
#testpath <- here("data", "lihc_chol_liri_gtex_summarizedExperiment_harmonized.RData")
#print(testpath)

load(here("data", "lihc_chol_liri_gtex_summarizedExperiment_harmonized.RData"))
gene_names <- rowData(liver_expr)$X
expr <- assays(liver_expr)$counts
texpr <- t(expr)

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
meta_df <- mutate_if(meta_df, is.character, as.factor)


# # DESEQ OBJECT WITHOUT DESIGN FOR FROZEN VST TRANSFORM
# # running DESeq makes the dispersionFunction available for VST transformation
# # count data: rows=genes,cols=samples
# # dds_train <- DESeqDataSetFromMatrix(countData = t(train_expr), 
#                     #colData = meta_df,
#                     #design = ~ 1) # no design
# #dds_train <- DESeq(dds_train)


# pp_sc <- preProcess(train_expr, method = c("scale", "center")) # 32163 genes
# train_sc <- predict(pp_sc, train_expr)

#####################  APPLY FUNCTIONS ##############################
# load scaling functions
load(here("data", "lihc_chol_liri_gtex_preproc.RData"))
train_expr <- predict(preproc_output$pp_nvz, texpr)

# load dds object for VST transformation
load(here("data", "lihc_chol_liri_gtex_dds_object_new_ids.RData"))
train_dispersionFunc <- dispersionFunction(dds_train) 
train_vst <- vst(dds_train, blind=T) 
train_sc <- predict(preproc_output$pp_sc, train_expr)

# load scaled data
load(here("data", "scaled_outputs.RData"))
# scaled_outputs$unitvar, scaled_outputs$log2_scaled, scaled_outputs$minmax

load(here("data", "pca_results_new_ids.RData"))

umap1 <- compute_UMAP(scaled_outputs$unitvar, meta_df, 42)
umap2 <- compute_UMAP(t(assay(train_vst)), meta_df, 42)
umap3 <- compute_UMAP(scaled_outputs$log2_scaled, meta_df, 42)
umap4 <- compute_UMAP(scaled_outputs$minmax, meta_df, 42)
print("computed umaps")

# free space
rm(sexpr, expr, texpr, log2_scaled,  unitvar, minmax)


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
survival_colors <- brewer.pal(n=9, name="Blues")
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

plotly_plotting_function <- function(plot_type, pcx, pcy, scaling_method, colorby, row_id){
  
  if (scaling_method == "unitvar"){
    if(plot_type == "PCA"){
      df_out <- pca_results$pca_unitvar$df_pca
      explained_variance_ratio = pca_results$pca_unitvar$evar
    } else if(plot_type == "UMAP"){
      df_out <- umap1$df_umap
    }
  }
  if (scaling_method == "vst"){
    if(plot_type == "PCA"){
      df_out <- pca_results$pca_vst$df_pca
      explained_variance_ratio = pca_results$pca_vst$evar
    } else if(plot_type == "UMAP"){
      df_out <- umap2$df_umap
    }
  }
  if (scaling_method == "log"){
    if(plot_type == "PCA"){
      df_out <- pca_results$pca_log$df_pca
      explained_variance_ratio = pca_results$pca_log$evar
    } else if(plot_type == "UMAP"){
      df_out <- umap3$df_umap
    }    
  }
  if (scaling_method == "minmax"){
    if(plot_type == "PCA"){
      df_out <- pca_results$pca_minmax$df_pca
      explained_variance_ratio = pca_results$pca_minmax$evar
    } else if(plot_type == "UMAP"){
      df_out <- umap4$df_umap
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
  
  # dynamically change marker size, color, border and opacity when row is selected in datatable
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
  
  fig <- fig %>% layout(title= list(text = plot_type, xanchor="left", x=0.1), 
                        list(title=list(text=colorby)),
                        xaxis = list(title = list(text=dim1)),
                        yaxis = list(title = list(text=dim2)))
  return(fig)
}

############################################ UI ####################################################

# UI with fluidPage
ui <- fluidPage(

  useShinyjs(), 
  add_busy_spinner(spin = "fading-circle"),  
  # Application title
  titlePanel("Anovaget transcriptomics"),
  br(),  
  sidebarLayout(
      
    # Sidebar with a slider input
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
      numericInput("pcy", "Principal omponent on y-axis:", 2, min=1, max=10, step=1)
        ),
      # Horizontal line ----
      tags$hr(),
      h3("Upload user data:"),
      fileInput("upload", NULL, multiple=F, width="100%",
                accept = c(".tsv"), buttonLabel="TSV"),
      numericInput("kn", "Number of nearest neighbors to compute:", 1, min=1, max=10, step=1),
      br(),
      actionButton("Add", "Click to add user data in plots"),
      bsTooltip(id = "Add", title = "Please input a  raw gene-sample expression matrix", 
                placement = "bottom", trigger = "hover"),     
      
      # Horizontal line ----
      tags$hr(),
      h3("Table options:"),
      collapsibleAwesomeCheckboxGroupInput("show_vars", "Select variables to show:\n", 3,
                          names(metadata)[3:length(names(metadata))], 
                          selected = names(metadata)[3:length(names(metadata))])
      ),
      
    # Show plots reduced dimensionality + tables
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
                          tabPanel("Neighbor Data", 
                                    br(),
                                    span(textOutput(outputId  = "tab2_text"), style="font-size: 20px; font-style: bold;"),
                                    br(),
                                    wellPanel(div(style = 'overflow-x: scroll', DT::dataTableOutput("neighbortable"), 
                                                                  style = "z-index: 10; left:0; right:0; overflow-y:hidden; overflow-xy:auto"))))
              
              
)))

########################################## SERVER FUNCTION #########################################
  
  # Server logic
server <- function(input, output, session) {

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
  
  # text in tabsets
  output$tab1_text <- renderText({"Metadata annotations from TCGA, ICGC and GTEx."})
  output$tab2_text <- renderText({"Please upload a dataset to compute approximate nearest neighbors to TCGA, ICGC and GTEx data. Nearest Neighbors computed on first 10 Principal Components."})
  

  # value to store datatable id
  id_sel <- reactiveVal()
  
  # observe if row in datatable is selected
  observe({
    idx <- input$trainingMetadata_rows_selected
    id_sel <- row.names(meta_df)[idx]
    id_sel_old_new <- c(id_sel(), id_sel)
    id_sel(unique(id_sel_old_new))
    
  })
  
  
  # clear the set of selected datapoints when a double-click occurs
  observeEvent(event_data("plotly_doubleclick"), {
    id_sel(NULL)
  })
  
  metaselect <- reactive(input$meta)
  scale_method <- reactive(input$scaling)
  pcx <- reactive(input$pcx)
  pcy <- reactive(input$pcy)
  kn <- reactive(input$kn)
      
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
        
      plotly_plotting_function(
        "UMAP",
        1, # umap in 2 just returns first 2 dimensions
        2,
        scaling_method=scale_method(),
        colorby=metaselect(),
        row_id=id_sel()
        )
    })

    ext_data <- reactive({
      req(input$upload)
      ext <- tools::file_ext(input$upload$name)
      validate(need(ext == "tsv", "Invalid file. Please upload a .tsv file"))
      df <- read.table(input$upload$datapath, sep = "\t", header=T)
      
      out_pp <- preProcessTestData(df, preproc_output$pp_nvz, preproc_output$train_dispersionFunc, preproc_output$pp_sc)
      return(out_pp)
    })

    observeEvent(input$Add,{      
      output$tab2_text <- renderText({"Nearest neighbors are being computed."})

      out = ext_data()
      # project testdata into PCA space
      pca_df <- projectTestData("PCA", out, scale_method(), colnames(meta_df))
      pca_df[is.na(pca_df)] <- "not_available"

      # compute nearest neighbors in PCA space with first 10 PCs
      neighbor_df <- compute_NN(scale_method(), pca_df[,1:10], kn())
      # render neighbor dataframe
      output$neighbortable <-  DT::renderDataTable({ DT::datatable(neighbor_df) })
      num_col <- ncol(neighbor_df)/2
      output$tab2_text <- renderText({"Nearest Neighbors to TCGA, ICGC, GTEx data."})

      # project testdata into UMAP space      
      umap_df <- projectTestData("UMAP", out, scale_method(), colnames(meta_df))
      umap_df[is.na(umap_df)] <- "not_available"

      # update plots based on user input
      tx_pca <- highlight_key(pca_df, ~row.names(pca_df)) 
      plotlyProxy("pcaPlot", session) %>%
        plotlyProxyInvoke("addTraces", x=tx_pca$data()[,pcx()], y=tx_pca$data()[,pcy()], color=tx_pca$data()[,metaselect()],
                          type="scatter", mode="markers", name=tx_pca$data()[,metaselect()][1],
                          hovertemplate = paste("Sample:", row.names(pca_df),  '<br>',
                          "Nearest Neighbor ", neighbor_df[,1], ", dist: ", neighbor_df[,1+num_col],
                          '<extra></extra>'),
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

########################################### RUN APP #################################################
  
  # Complete app with UI and server components
shinyApp(ui, server)
