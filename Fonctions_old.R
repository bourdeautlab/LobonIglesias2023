#######################################
# Most variable genes/probes function #
#######################################

genes.selection <- function (data,
                             thres.diff,
                             thres.num,
                             probs = 0.25
) 
{
  if (missing(thres.diff) && missing(thres.num))
  {
    stop("** Stop. No method found.")
  }
  if (!missing(thres.diff) && !missing(thres.num))
  {
    stop("** Stop. Choose one of the two options - thres.diff or thres.num")
  }
  diff.quantile <- function(x, probs = 1/length(x))
  {
    vv <- quantile(x, probs = c(probs, 1 - probs), na.rm = TRUE)
    return(vv[2] - vv[1])
  }
  rangeValues <- apply(data, 1, diff.quantile, probs = probs)
  if (!missing(thres.diff))
  {
    ind <- which(rangeValues >= thres.diff)
    genesList <- rownames(data[ind, ])
  }
  else if (!missing(thres.num))
  {
    genesList <- names(sort(rangeValues, decreasing = TRUE)[1:thres.num])
  }
  return(genesList)
}

################
# PCA function #
################

pca_func <- function(matrix,
                     sampleTable,
                     highlightedVar,
                     highlightedVarName         = highlightedVar,
                     colors                     = NULL,
                     sampleLabel                = NULL,
                     sampleLabelColumnCondition = NULL,
                     sampleLabelValueCondition  = NULL,
                     fontsize                   = 40,
                     label_fontsize             = 8,
                     point_size                 = 6,
                     scale                      = FALSE 
)
{
  # Load packages
  library(factoextra)
  
  # PCA analysis
  PC        <- prcomp(x = data.frame(t(matrix)),
                      scale. = scale
  )
  eigVal    <- get_eigenvalue(PC)
  eigValPer <- round(eigVal$variance.percent, digits = 2)
  pci       <- data.frame(PC$x,sampleTable)
  
  # Plot PCA
  g <- ggplot(data = pci,
              mapping = aes(x = PC1, y = PC2)
  ) +
    geom_point(mapping = aes(fill = eval(parse(text=highlightedVar))),
               size    = point_size,
               pch     = 21
    ) +
    scale_fill_manual(name   = highlightedVarName,
                      values = colors[levels(factor(pci[,highlightedVar]))]
    ) +
    xlab(paste0("PC1 (",eigValPer[1]," %)")) +
    ylab(paste0("PC2 (",eigValPer[2]," %)")) +
    ggtitle("PCA") +
    theme_bw()+
    theme(text              = element_text(size = fontsize),
          legend.key.height = unit(1.5,"cm")
    )
  if (!is.null(sampleLabel))
  {
    if (is.null(sampleLabelColumnCondition) && is.null(sampleLabelValueCondition))
    {
      g <- g +
        geom_text_repel(mapping = aes(label = eval(parse(text=sampleLabel))),
                        size    = label_fontsize
        )
    }
    else if (!is.null(sampleLabelColumnCondition) && is.null(sampleLabelValueCondition))
    {
      stop("missing sampleLabelValueCondition, it is required since sampleLabelColumnCondition is provided", call.=TRUE)
    }
    else if (!is.null(sampleLabelColumnCondition) && !is.null(sampleLabelValueCondition))
    {
      if (!is.element(el = sampleLabelValueCondition, set = as.vector(pci[,sampleLabelColumnCondition])))
      {
        warning(paste(sampleLabelColumnCondition,"column does not contain",sampleLabelValueCondition,"value. Therefore, no sample is labeled", sep = " "), call.=TRUE)
      }
      g <- g +
        geom_text_repel(data = dplyr::filter(.data = pci,
                                             eval(parse(text=sampleLabelColumnCondition)) %in% sampleLabelValueCondition
        ),
        mapping = aes(label = eval(parse(text=sampleLabel))),
        size = label_fontsize
        )
    }
  }
  
  print(g)
}

############################
# Complex heatmap function #
############################

createTopAnnot_func <- function(sampleTable,
                                columns,
                                colorsList,
                                height               = 2,
                                width                = 2,
                                fontsize             = 25,
                                which                = "column",
                                show_legend          = TRUE,
                                show_annotation_name = TRUE
)
{
  # Load package
  library(ComplexHeatmap)
  
  # Check if vectors of column name and vectors of colors have the same length
  #if (!identical(x = length(columns), y = length(colorsList)))
  #{
   # stop("length of colors vector should be identical as the length of selected columns.",
    #     call.=TRUE
    #)
  #}
  # Create list of colors
  names(colorsList) <- columns
  heatmapAnnot <- HeatmapAnnotation(df                      = data.frame(sampleTable[columns], check.names = FALSE),
                                    col                     = colorsList,
                                    gap                     = unit(2,"mm"),
                                    show_annotation_name    = show_annotation_name,
                                    annotation_name_gp      = gpar(fontsize    = fontsize),
                                    annotation_legend_param = list(fontsize    = fontsize,
                                                                   title_gp    = gpar(fontsize = fontsize, font = 2),
                                                                   labels_gp   = gpar(fontsize = fontsize),
                                                                   grid_width  = unit(12,"mm"),
                                                                   grid_height = unit(12,"mm")
                                    ),
                                    height                  = unit(height,"cm"),
                                    width                   = unit(width,"cm"),
                                    na_col                  = "white",
                                    which                   = which,
                                    show_legend             = show_legend
  )
  return(heatmapAnnot)
}

drawHeatmap_func <- function(matrix,
                             col                     = NULL,
                             column_title            = "",
                             row_title               = "",
                             top_annot               = NULL,
                             row_annot               = NULL,
                             bottom_annot            = NULL,
                             legend_title            = "expression",
                             show_column_names       = TRUE,
                             show_row_dend           = FALSE,
                             show_row_names          = FALSE,
                             add_column_dend_axis    = FALSE,
                             columns_dist            = "pearson",
                             columns_meth            = "ward",
                             cluster_columns         = TRUE,
                             cluster_rows            = TRUE,
                             rows_dist               = "pearson",
                             rows_meth               = "ward",
                             row_title_gp            = gpar(fontsize = fontsize),
                             row_title_rot           = 90,
                             fontsize                = 25,
                             column_order               = NULL,
                             row_order               = NULL,
                             split                   = NULL,
                             annotation_legend_side  = "right",
                             column_dend_side        = "top",
                             column_subtitle         = paste0("\n Clustering (distance: ",columns_dist,", method: ",columns_meth,")"),
                             row_subtitle            = paste0("\n Clustering (distance: ",rows_dist,", method: ",rows_meth,")"),
                             heatmap_legend_side     = "left",
                             show_heatmap_legend     = TRUE,
                             annotation_legend_list  = NULL
)
{
  # Load package
  library(ComplexHeatmap)
  
  # Create heatmap object
  if (!is.null(col))
  {
    h <- Heatmap(matrix                      = matrix,
                 col                         = col,
                 name                        = legend_title,
                 cluster_columns             = cluster_columns,
                 cluster_rows                = cluster_rows,
                 clustering_distance_columns = columns_dist,
                 clustering_method_columns   = columns_meth,
                 clustering_distance_rows    = rows_dist,
                 clustering_method_rows      = rows_meth,
                 show_row_names              = show_row_names,
                 show_row_dend               = show_row_dend,
                 show_column_names           = show_column_names,
                 column_names_gp             = gpar(fontsize = fontsize),
                 row_names_gp                = gpar(fontsize = fontsize),
                 column_title                = paste0(column_title,column_subtitle),
                 row_title                   = paste0(row_title,row_subtitle),
                 top_annotation              = top_annot,
                 bottom_annotation           = bottom_annot,
                 show_heatmap_legend         = show_heatmap_legend,
                 column_dend_height          = unit(40, "mm"),
                 row_dend_width              = unit(40, "mm"),
                 column_title_gp             = gpar(fontsize = fontsize),
                 row_title_gp                = row_title_gp,
                 row_title_rot               = row_title_rot,
                 column_order                = column_order,
                 row_order                   = row_order,
                 split                       = split,
                 column_dend_side            = column_dend_side,
                 heatmap_legend_param        = list(color_bar      = "continuous",
                                                    grid_width     = unit(0.75,"cm"),
                                                    grid_height    = unit(2,"cm"),
                                                    legend_width   = unit(2,"cm"),
                                                    title          = paste0(legend_title,"\n"),
                                                    title_gp       = gpar(fontsize = fontsize),
                                                    legend_gp      = gpar(fontsize = fontsize),
                                                    title_position = "topcenter"
                 )
    )
  }
  else
  {
    h <- Heatmap(matrix                      = matrix,
                 name                        = legend_title,
                 cluster_columns             = cluster_columns,
                 cluster_rows                = cluster_rows,
                 clustering_distance_columns = columns_dist,
                 clustering_method_columns   = columns_meth,
                 clustering_distance_rows    = rows_dist,
                 clustering_method_rows      = rows_meth,
                 show_row_names              = show_row_names,
                 show_row_dend               = show_row_dend,
                 show_column_names           = show_column_names,
                 column_names_gp             = gpar(fontsize = fontsize),
                 row_names_gp                = gpar(fontsize = fontsize),
                 column_title                = paste0(column_title,column_subtitle),
                 row_title                   = paste0(row_title,row_subtitle),
                 top_annotation              = top_annot,
                 bottom_annotation           = bottom_annot,
                 show_heatmap_legend         = show_heatmap_legend,
                 column_dend_height          = unit(40, "mm"),
                 row_dend_width              = unit(40, "mm"),
                 column_title_gp             = gpar(fontsize = fontsize),
                 row_title_gp                = row_title_gp,
                 row_title_rot               = row_title_rot,
                 column_order                = column_order,
                 row_order                   = row_order,
                 split                       = split,
                 column_dend_side            = column_dend_side,
                 heatmap_legend_param        = list(color_bar      = "continuous",
                                                    grid_width     = unit(0.75,"cm"),
                                                    grid_height    = unit(2,"cm"),
                                                    legend_width   = unit(2,"cm"),
                                                    title          = paste0(legend_title,"\n"),
                                                    title_gp       = gpar(fontsize = fontsize),
                                                    legend_gp      = gpar(fontsize = fontsize),
                                                    title_position = "topcenter"
                 )
    )
  }
  if (!is.null(row_annot))
  {
    h <- h + row_annot
  }
  draw(h,
       heatmap_legend_side    = heatmap_legend_side,
       annotation_legend_side = annotation_legend_side,
       annotation_legend_list = annotation_legend_list
  )
  if (add_column_dend_axis){
    decorate_column_dend(legend_title, {grid.yaxis()})
  }
  #return(h)
}


########################################
# Plot volcano plot with Edge results  #
########################################

plot_volcano_edgeR <- function(edgeR_result,
                               title              = NULL,
                               subtitle = NULL,
                               pointSize          = 2,
                               gene_set           = "",
                               gene_set_label     = FALSE,
                               gene_set_color     = "orange3",
                               gene_set2          = "",
                               gene_set2_label    = FALSE,
                               gene_set2_color    = "magenta4",
                               gene_set3          = "",
                               gene_set3_label    = FALSE,
                               gene_set3_color    = "magenta4",
                               gene_set4          = "",
                               gene_set4_label    = FALSE,
                               gene_set4_color    = "magenta4",
                               gene_set_pointSize = 8,
                               FDR_cutoff         = 0.05,
                               log2FC_cutoff      = 1,
                               FDR_cutoff_lab     = 0,
                               log2FC_cutoff_lab  = 1000,
                               geneSup            = "",
                               ggtheme            = NULL,
                               fontsize           = 30,
                               background         = "white",
                               linetype           = 1,
                               geneSup_pointsize  = 10,
                               geneSup_fontsize   = 15,
                               colors             = c("forestgreen","blue", "grey","brown","black"),
                               ylim               = NULL
)
{
  # Load library
  library(dplyr)    # -> mutate(), filter()
  library(ggplot2)  # -> ggplot(), geom_point(), ggtitle(), ...
  library(ggrepel)  # -> geom_text_repel()
  # Add criteria column to be used to draw cut-off lines in the volcano plot
  edgeR_result <- data.frame(topTags(edgeR_result, n=Inf))
  edgeR_result$symbol <- rownames(edgeR_result)
  edgeR_result_crit <- mutate(.data  = edgeR_result,
                              crit   = ifelse(edgeR_result$FDR < FDR_cutoff &
                                                edgeR_result$logFC > log2FC_cutoff,
                                              "1",
                                              ifelse(edgeR_result$FDR < FDR_cutoff &
                                                       edgeR_result$logFC < -log2FC_cutoff,
                                                     "2",
                                                     ifelse(edgeR_result$FDR > FDR_cutoff &
                                                              (edgeR_result$logFC < log2FC_cutoff
                                                               & edgeR_result$logFC > -log2FC_cutoff) ,
                                                            "3",
                                                            ifelse(edgeR_result$FDR < FDR_cutoff &
                                                                     (edgeR_result$logFC < log2FC_cutoff
                                                                      & edgeR_result$logFC > -log2FC_cutoff),
                                                                   "4",
                                                                   "5"
                                                            )
                                                     )
                                              )
                              ))
  # Plot volcanoplot
  g <- ggplot(data = edgeR_result_crit,
              mapping = aes(x=logFC, y=-log10(FDR))
  ) +
    xlab(label=expression(paste(log[2],"(fold-change)")))+
    ylab(label=expression(paste(-log[10],"(adj. p-value)")))+
    geom_point(mapping = aes(col=crit), size = pointSize) +
    ggtitle(label = title, subtitle = subtitle) +
    scale_color_manual(values = colors) +
    geom_text_repel(data    = filter(.data = edgeR_result_crit, FDR < FDR_cutoff_lab,
                                     abs(logFC) > log2FC_cutoff_lab,
                                     !(symbol %in% geneSup)
    ),
    mapping = aes(label=symbol)#,
    #size    = 14
    ) +
    # Add gene set points
    geom_point(data    = filter(edgeR_result_crit, symbol %in% gene_set, !(symbol %in% geneSup)),
               mapping = aes(x=logFC, y=-log10(FDR)),
               col     = gene_set_color,
               size    = gene_set_pointSize
    ) +
    geom_point(data    =filter(edgeR_result_crit, symbol %in% gene_set2, !(symbol %in% geneSup)),
               mapping = aes(x=logFC, y=-log10(FDR)),
               col     = gene_set2_color,
               size    = gene_set_pointSize
    ) +
    geom_point(data    = filter(edgeR_result_crit, symbol %in% gene_set3, !(symbol %in% geneSup)),
               mapping = aes(x=logFC, y=-log10(FDR)),
               col     = gene_set3_color,
               size    = gene_set_pointSize
    ) +
    geom_point(data    = filter(edgeR_result_crit, symbol %in% gene_set4, !(symbol %in% geneSup)),
               mapping = aes(x=logFC, y=-log10(FDR)),
               col     = gene_set4_color,
               size    = gene_set_pointSize
    ) +
    # Add specific gene point and label
    geom_point(data=filter(edgeR_result_crit, symbol %in% geneSup),
               mapping = aes(x=logFC, y=-log10(FDR)),
               col="red",
               size=geneSup_pointsize
    ) +
    geom_text_repel(data=filter(edgeR_result_crit, symbol %in% geneSup),
                    mapping = aes(label=symbol),
                    size = geneSup_fontsize,
                    fontface = "bold"
    ) +
    geom_hline(yintercept = -log10(FDR_cutoff), color= "red", lty = linetype)+
    geom_vline(xintercept = log2FC_cutoff, color= "blue", lty = linetype)+
    geom_vline(xintercept = -log2FC_cutoff, color= "blue", lty = linetype)+
    theme(legend.position = "none",
          text = element_text(size = fontsize),
          plot.background  = element_rect(fill = background),
          panel.background = element_rect(fill = background),
          axis.line        = element_line(colour = "grey", linetype = 1),
          axis.text.x  = element_text(color = "black", size = 30),
          axis.title.x = element_text(color = "black", size = 40),
          axis.text.y  = element_text(color = "black", size = 30),
          axis.title.y = element_text(color = "black", size = 40),
          title = element_text(color = "black", size = 40, face = "bold")
    )
  # theme
  if (!is.null(ggtheme))
  {
    g <- g + ggtheme + theme(legend.position = "none")
  }
  # y lim
  if (!is.null(ylim))
  {
    g <- g + ylim(ylim[1],ylim[2])
  }
  # Add gene set labels
  if (gene_set_label)
  {
    g <- g + geom_text_repel(data=filter(edgeR_result_crit, symbol %in% gene_set, !(symbol %in% geneSup)), aes(label=symbol), size=15, fontface = "bold")
  }
  if (gene_set2_label)
  {
    g <- g + geom_text_repel(data=filter(edgeR_result_crit, symbol %in% gene_set2, !(symbol %in% geneSup)), aes(label=symbol), size=15, fontface = "bold")
  }
  if (gene_set3_label)
  {
    g <- g + geom_text_repel(data=filter(edgeR_result_crit, symbol %in% gene_set3, !(symbol %in% geneSup)), aes(label=symbol), size=15, fontface = "bold")
  }
  if (gene_set4_label)
  {
    g <- g + geom_text_repel(data=filter(edgeR_result_crit, symbol %in% gene_set4, !(symbol %in% geneSup)), aes(label=symbol), size=15, fontface = "bold")
  }
  # Draw volcano plot
  print(g)
}

###########################################
# Barplot function for GO analysis result #
###########################################

format_GO_labels <- function(x, prefix = "GO", wordsLength = 4,...)
{
  gsub(pattern = "((?:[^_]+_){",wordsLength,"}[^_]+)", replacement = "\\1\n", x) %>%
    gsub(pattern = prefix, replacement = "", x) %>%
    gsub(pattern = "_", replacement = " ", x)
}

plot_top_MSigDB_enriched_geneSets_func <- function(MSigDB_overlap_res,
                                                   top               = 10,
                                                   title             = NULL,
                                                   xlim              = 400,
                                                   sep               = "\t",
                                                   quote             = NULL,
                                                   selected_geneSet  = NULL,
                                                   fill              = "darkgreen",
                                                   format_labels     = FALSE,
                                                   prefixToremove   = "GO",
                                                   fontsize          = 40,
                                                   label_fontsize    = 40,
                                                   wordsLength       = 4
)
{
  # Import overlap result
  res_table <- read.table(file         = MSigDB_overlap_res,
                          header       = FALSE,
                          skip         = 10,
                          sep          = sep,
                          fill         = NA,
                          #row.names    = FALSE,
                          #col.names    = FALSE,
                          check.names  = FALSE
  )
  
  # Select pathway table result
  res_table <- res_table[,1:7]
  
  # Add colnames
  colnames(res_table) = c("pathway_ID","K","pathway_name","k","k_K","pvalue","FDR")
  
  # Ensure to have numeric data for numeric columns
  res_table$K      <- as.numeric(as.character(res_table$K))
  res_table$k_K    <- as.numeric(as.character(res_table$k_K))
  res_table$pvalue <- as.numeric(as.character(res_table$pvalue))
  res_table$FDR    <- as.numeric(as.character(res_table$FDR))
  
  # Remove non-numeric value on numeric columns (this allow us to keep only the result table)
  res_table        <- res_table[!is.na(res_table$K),]
  res_table        <- res_table[!is.na(res_table$k_K),]
  res_table        <- res_table[!is.na(res_table$pvalue),]
  res_table        <- res_table[!is.na(res_table$FDR),]
  
  # Order pathway name according to FDR 
  res_table$pathway_ID <- factor(x      = res_table$pathway_ID,
                                 levels = res_table$pathway_ID[order(res_table$FDR,
                                                                     decreasing = TRUE
                                 )
                                 ]
  )
  
  # Select top genes
  if (nrow(res_table) < top)
  {
    top <- nrow(res_table)
    res_table_top <- res_table[1:top,]
  }
  else
  {
    res_table_top <- res_table[1:top,]
  }
  
  # Give different color and face for selected gene set names
  res_table_top$y_text_color <- ifelse(res_table_top$pathway_ID %in% selected_geneSet,
                                       "black",
                                       "black"
  )
  res_table_top$y_text_face  <- ifelse(res_table_top$pathway_ID %in% selected_geneSet,
                                       "bold",
                                       "plain"
  )
  
  # Format labels
  if (format_labels)
  {
    res_table_top$pathway_ID <- format_GO_labels(res_table_top$pathway_ID, 
                                                 prefix = prefixToremove,
                                                 wordsLength = wordsLength)
  }
  
  
  # Plot barplot
  g <- ggplot(data    = res_table_top,
              mapping = aes(x = reorder(pathway_ID,-log10(FDR)),
                            y = -log10(FDR)
              )
  ) +
    geom_bar(stat = "identity", fill = fill) +
    coord_flip() +
    xlab("") +
    ggtitle(label = title) +
    ylab(expression(-log[10](FDR))) +
    ylim(0,xlim) +
    geom_text(mapping = aes(label=paste0("(",k,"/",K,")")), size = 5, hjust = -0.1) +
    theme_bw() +
    theme(title = element_text(color = "black", size = fontsize),
          panel.border = element_blank(),
          axis.text.x  = element_text(color = "black", size = fontsize),
          axis.title.x = element_text(color = "black", size = fontsize),
          axis.text.y  = element_text(color = rev(res_table_top$y_text_color),
                                      face  = rev(res_table_top$y_text_face),
                                      size  = label_fontsize
          )
    )
  return(g)    
}


##############################################
# Filtering matrix based on expression level #
##############################################

expFilter <- function(data,
                      threshold  = 3.5,
                      p          = 0.01,
                      graph      = TRUE
) 
{
  if (graph) {
    par(font.lab = 2)
    hist(data, breaks = 200, main = "Data distribution", 
         xlab = "Expression Level")
    abline(v = threshold, col = 2)
    mtext(threshold, at = threshold, side = 1, las = 1, col = 2, 
          cex = 0.6)
  }
  vecL <- apply(data, 1, function(x) {
    length(which(x > threshold))
  })
  cat("Keep genes with at least", ceiling(ncol(data) * p), 
      "sample(s) with an expression level higher than", threshold, 
      "\n")
  data <- data[which(vecL >= ceiling(ncol(data) * p)), ]
}

##########################
# Create GSEA cls format #
##########################

createGSEAclsFormat_func <- function(class,
                                     classFileName
)
{
  class      <- as.vector(class)
  class_file <- file(classFileName)
  writeLines(text = paste(as.character(length(class)),
                          as.character(length(unique(class))),
                          "1",
                          sep = "\t"
  ),
  con  = class_file
  )
  close(class_file)
  cat(c("#", unique(class), "\n"),
      sep    = "\t",
      file   = classFileName,
      append = TRUE
  )
  cat(class,
      sep    = "\t",
      file   = classFileName,
      append = TRUE
  )
}


##########################
# Create GSEA gct format #
##########################

createGSEAgctFormat_func <- function(matrix,
                                     gctFileName
)
{
  # Convert matrix to data frame
  matrix <- data.frame(matrix, check.names=FALSE)
  # Create file object
  gct_file <- file(gctFileName)
  # Write the first line of gct format (ie #1.2) in the created file
  writeLines(text = "#1.2",
             con  = gct_file
  )
  close(gct_file)
  # Add the second line of gct format (nb of genes followed by the number of samples)
  cat(c(nrow(matrix), ncol(matrix), "\n"),
      sep    = "\t",
      file   = gctFileName,
      append = TRUE
  )
  # Add NAME and Description columns to the matrix
  matrix_gct             <- matrix
  matrix_gct$NAME        <- rownames(matrix)
  matrix_gct$Description <- rownames(matrix)
  # Put NAME and Description columns at the begining
  matrix_gct <- cbind(matrix_gct[c("NAME", "Description")],
                      matrix_gct[,(3:ncol(matrix_gct)-2)]
  )
  # Export gct format matrix
  write.table(x         = matrix_gct,
              file      = gctFileName,
              quote     = FALSE,
              sep       = "\t",
              row.names = FALSE,
              append    = TRUE
  )
}


###############################
# Plot barplot of GSEA result #
###############################

barplotGSEA_func <- function(GSEAresult,
                             top              = 20,
                             #xlim             = NULL,
                             selected_geneSet = NULL#,
                             #sep              = "\t"
)
{
  stopifnot(require(readr))
  # Import GSEA result (csv format)
  GSEAresult <- read_tsv(file      = GSEAresult,
                         #sep       = sep,
                         skip      = 1,
                         col_names = c("name",
                                       "link",
                                       "details",
                                       "size",
                                       "ES",
                                       "NES",
                                       "NOMpval",
                                       "FDRqval",
                                       "FWER",
                                       "rank_at_max",
                                       "leading_edge"
                         )
  )
  
  # Convert NES to its absolute value
  GSEAresult$NES <- abs(GSEAresult$NES)
  
  # Order pathway name according to FDR 
  GSEAresult$name <- factor(x      = GSEAresult$name,
                            levels = GSEAresult$name[order(GSEAresult$NES, decreasing = FALSE)])
  
  # Select top enriched gene sets
  GSEAresult_top              <- GSEAresult[1:top,]
  
  # Give different color and face for selected gene set names
  GSEAresult_top$y_text_color <- ifelse(GSEAresult_top$name %in% selected_geneSet,"blue","black")
  GSEAresult_top$y_text_face  <- ifelse(GSEAresult_top$name %in% selected_geneSet,"bold","plain")
  
  # Plot barplot
  g <- ggplot(data = GSEAresult_top,
              mapping = aes(x = name, y = NES)
  ) +
    geom_bar(stat = "identity",
             fill = "chocolate4",
             alpha = 0.3
    ) +
    coord_flip() +
    geom_bar(data = filter(GSEAresult_top, name %in% selected_geneSet),
             stat = "identity",
             fill = "chocolate4"
    ) +
    geom_line(mapping = aes(x =name, y = -log10(FDRqval), group = 1),
              color   = "darkgreen"
    ) +
    geom_point(mapping = aes(x =name, y = -log10(FDRqval), group = 1),
               color   ="darkgreen"
    ) +
    scale_y_continuous(sec.axis = sec_axis(trans = ~.*1, name = "-log10(FDR)\n")) +
    xlab("") +
    ylab("\nNES") +
    #geom_text(mapping = aes(label=paste0("(",k,"/",K,")")), size = 7.5, hjust = -0.1) +
    theme_bw() +
    theme(text         = element_text(size=30),
          panel.border = element_blank(),
          axis.text.x  = element_text(color = "black"),
          axis.text.y  = element_text(color = rev(GSEAresult_top$y_text_color),
                                      face  = rev(GSEAresult_top$y_text_face)
          )
    )
  #if (!is.null(xlim)) {g <- g + ylim(xlim[1],xlim[2])}
  print(g)
}

######################
# Plot Venn Diagramm #
######################

plot_2_venn_diagram <- function(x,
                                y,
                                x_name ="",
                                y_name ="",
                                x_fill = "blue",
                                y_fill = "brown",
                                x_col  = "black",
                                y_col  = "black",
                                x_lwd  = 2,
                                y_lwd  = 2,
                                labPos = c(0,0),
                                ind    = TRUE,
                                dev.new = TRUE
)
{
  
  library(VennDiagram)
  
  if (dev.new)
  {
    dev.new()
  }
  
  draw.pairwise.venn(area1      = length(x),
                     area2      = length(y),
                     cross.area = length(
                       intersect(x,y)
                     ),
                     category   = c(x_name,y_name),
                     fill       = c(x_fill,y_fill), # fill color
                     lty        = rep(1,2),         # bordure des cercles
                     lwd        = c(x_lwd,y_lwd),
                     col        = c(x_col,y_col),
                     alpha      = rep(0.5,2),       # transparence des cercles
                     cat.pos    = labPos,           # label position
                     cex        = rep(3,3),
                     cat.dist   = rep(0.05,2),      # label distance
                     cat.cex    = rep(3,2),
                     cat.fontface = "bold",
                     ind          = ind
  )
}

plot_4_venn_diagram <- function(x,
                                y,
                                z,
                                u,
                                x_name ="",
                                y_name ="",
                                z_name ="",
                                u_name ="",
                                x_fill = NULL,
                                y_fill = NULL,
                                z_fill = NULL,
                                u_fill = NULL,
                                x_col  = "blue",
                                y_col  = "brown",
                                z_col  = "darkgreen",
                                u_col  = "red",
                                dev.new = TRUE
)
{
  
  library(VennDiagram)
  
  if (dev.new)
  {
    dev.new()
  }
  
  draw.quad.venn(area1         = length(x),
                 area2         = length(y),
                 area3         = length(z),
                 area4         = length(u),
                 n12           = length(intersect(x,y)),
                 n13           = length(intersect(x,z)),
                 n14           = length(intersect(x,u)),
                 n23           = length(intersect(y,z)),
                 n24           = length(intersect(y,u)),
                 n34           = length(intersect(z,u)),
                 n123          = length(Reduce(intersect, list(x,y,z))),
                 n124          = length(Reduce(intersect, list(x,y,u))),
                 n134          = length(Reduce(intersect, list(x,z,u))),
                 n234          = length(Reduce(intersect, list(y,z,u))),
                 n1234         = length(Reduce(intersect, list(x,y,z,u))),
                 category      = c(x_name,y_name,z_name,u_name),
                 fill          = c(x_fill,y_fill,z_fill,u_fill), # fill color
                 lty           = rep(1,4),         # bordure des cercles
                 lwd           = rep(2,4),
                 col           = c(x_col,y_col,z_col,u_col), # border color
                 alpha         = rep(0.5,4),       # transparence des cercles
                 cat.pos       = c(-9,-11,-1,-3),       # label position
                 cex           = rep(3,15),
                 cat.dist      = rep(0.05,4),      # label distance
                 cat.cex       = rep(3,4),
                 cat.fontface  = "bold"
  )
}


#######################################
# Compute hypergeometric test p-value #
#######################################

computeHypergeometrix_pval <- function(x,
                                       y,
                                       intersect = NULL,
                                       N)
{
  
  if (mode(x) == "character" & mode(y) == "character")
  {
    intersect <- intersect(x,y)
    
    phyper(q = length(intersect) - 1,
           m = length(x),
           n = N - length(x),
           k = length(y),
           lower.tail = FALSE
    )
  }
  
  else if (mode(x) == "numeric" & mode(y) == "numeric" & length(x) == 1 & length(y) == 1)
  {
    if (is.null(intersect)) {stop("intersect value must be provided if x and y are single numeric values", call.=FALSE)}
    
    phyper(q = intersect - 1,
           m = x,
           n = N - x,
           k = y,
           lower.tail = FALSE
    )
  }
  else
  {
    stop("x and y need to be a character verctor or a single value", call.=FALSE)
  }
  
}

#################
# UMAP function #
#################

umap_func <- function(matrix,
                      sampleTable,
                      title                      = "UMAP",
                      n_neighbors                = 15,
                      metric                     = "euclidean",
                      highlightedVar,
                      highlightedVarName         = highlightedVar,
                      colors                     = NULL,
                      seed                       = 9,
                      sampleLabel                = NULL,
                      sampleLabelColumnCondition = NULL,
                      sampleLabelValueCondition  = NULL,
                      fontsize                   = 40,
                      pointsize                  = 6,
                      stat_ellipse               = FALSE,
                      showLegend                 = TRUE,
                      x_axis_position            = "bottom"
)
{
  
  # Load packages
  require(umap)
  require(ggplot2)
  require(ggrepel)
  
  # Set seed
  set.seed(seed = seed)
  
  # Creating a projection
  umap <- umap(t(matrix),
               n_neighbors = n_neighbors,
               metric = metric
  )
  
  # Extract main component
  umap_layout <- umap$layout
  colnames(umap_layout) <- c("UMAP1","UMAP2")
  
  # Bind sample plan
  umap_layout_plus <- cbind(sampleTable,
                            umap_layout
  )
  
  # Plot UMAP
  g <- ggplot(data = umap_layout_plus,
              mapping = aes(x = UMAP1, y = UMAP2)
  ) +
    geom_point(aes(fill = eval(parse(text=highlightedVar))),
               size = pointsize,
               pch = 21
    ) +
    scale_fill_manual(name   = highlightedVarName,
                      values = colors[levels(factor(umap_layout_plus[,highlightedVar]))]
    ) +
    scale_x_continuous(position = x_axis_position) +
    xlab("UMAP1") +
    ylab("UMAP2") +
    ggtitle(title) +
    theme_bw()+
    theme(text = element_text(size = fontsize)
    )
  # Show legend or not
  if (showLegend)
  {
    g <- g + theme(legend.title = element_text(face = "bold"),
                   legend.key.height = unit(1.5,"cm")
    )
  }
  else { g <- g + theme(legend.position = "none")}
  # Sample label
  if (!is.null(sampleLabel))
  {
    if (is.null(sampleLabelColumnCondition) && is.null(sampleLabelValueCondition))
    {
      g <- g +
        geom_text_repel(mapping = aes(label = eval(parse(text=sampleLabel))),
                        size    = fontsize/5
        )
    }
    else if (!is.null(sampleLabelColumnCondition) && is.null(sampleLabelValueCondition))
    {
      stop("missing sampleLabelValueCondition, it is required since sampleLabelColumnCondition is provided", call.=TRUE)
    }
    else if (!is.null(sampleLabelColumnCondition) && !is.null(sampleLabelValueCondition))
    {
      if (!is.element(el = sampleLabelValueCondition, set = as.vector(umap_layout_plus[,sampleLabelColumnCondition])))
      {
        warning(paste(sampleLabelColumnCondition,"column does not contain",sampleLabelValueCondition,"value. Therefore, no sample is labeled", sep = " "), call.=TRUE)
      }
      g <- g +
        geom_text_repel(data = dplyr::filter(.data = umap_layout_plus,
                                             eval(parse(text=sampleLabelColumnCondition)) %in% sampleLabelValueCondition
        ),
        mapping = aes(label = eval(parse(text=sampleLabel))),
        size = fontsize/5
        )
    }
  }
  
  if (stat_ellipse)
  {
    g + stat_ellipse()
  }
  
  print(g)
}

#################
# tSNE function #
#################

tsne_func <- function(matrix,
                      sampleTable,
                      highlightedVar,
                      highlightedVarName         = highlightedVar,
                      colors                     = NULL,
                      seed                       = NULL,
                      perplexity                 = 30,
                      sampleLabel                = NULL,
                      sampleLabelColumnCondition = NULL,
                      sampleLabelValueCondition  = NULL,
                      fontsize                   = 40,
                      pointSize        = 2
)
{
  # Load packages
  library(Rtsne)
  require(ggplot2)
  require(ggrepel)
  # Set seed
  set.seed(seed = seed)
  
  # Comput t-SNE
  tsne_model <- Rtsne(X                = t(matrix),
                      dims             = 2,
                      perplexity       = perplexity,
                      theta            = 0.0,
                      pca              = TRUE,
                      check_duplicates = FALSE,
                      is_distance      = FALSE
  )
  
  # Get t-SNE matrix
  tsne_matrix <- as.data.frame(tsne_model$Y)
  
  # Bind sample plan
  tsne_matrix_plus <- cbind(sampleTable, tsne_matrix)
  
  # Plot t-SNE
  g <- ggplot(data = tsne_matrix_plus,
              mapping = aes(x = V1, y = V2)
  ) +
    geom_point(aes(fill = eval(parse(text=highlightedVar))),
               size = pointSize,
               pch = 21
    ) +
    scale_fill_manual(name   = highlightedVarName,
                      values = colors[levels(factor(tsne_matrix_plus[,highlightedVar]))]
    ) +
    xlab("tSNE1") +
    ylab("tSNE2") +
    ggtitle("tSNE") +
    theme_bw()+
    theme(text = element_text(size = fontsize),
          legend.title = element_text(face = "bold"),
          legend.key.height = unit(1.5,"cm")
    )
  if (!is.null(sampleLabel))
  {
    if (is.null(sampleLabelColumnCondition) && is.null(sampleLabelValueCondition))
    {
      g <- g +
        geom_text_repel(mapping = aes(label = eval(parse(text=sampleLabel))),
                        size    = fontsize/5
        )
    }
    else if (!is.null(sampleLabelColumnCondition) && is.null(sampleLabelValueCondition))
    {
      stop("missing sampleLabelValueCondition, it is required since sampleLabelColumnCondition is provided", call.=TRUE)
    }
    else if (!is.null(sampleLabelColumnCondition) && !is.null(sampleLabelValueCondition))
    {
      if (!is.element(el = sampleLabelValueCondition, set = as.vector(tsne_matrix_plus[,sampleLabelColumnCondition])))
      {
        warning(paste(sampleLabelColumnCondition,"column does not contain",sampleLabelValueCondition,"value. Therefore, no sample is labeled", sep = " "), call.=TRUE)
      }
      g <- g +
        geom_text_repel(data = dplyr::filter(.data = tsne_matrix_plus,
                                             eval(parse(text=sampleLabelColumnCondition)) %in% sampleLabelValueCondition
        ),
        mapping = aes(label = eval(parse(text=sampleLabel))),
        size    = fontsize/5
        )
    }
  }
  
  print(g)
  
}

###############################################
# Assign subgroups of hierarchical clustering #
###############################################

assignGroups_func <- function(matrix,
                              sampleTable,
                              distance          = "pearson",
                              method            = "ward",
                              k                 = 3,
                              group_names       = c("groupA","groupB","groupC"),
                              group_num_column  = "group_num",
                              group_name_column = "group_name",
                              silhouette        = FALSE
)
{
  
  # Hierarchical clustering
  pl <- Heatmap(matrix                      = matrix,
                clustering_distance_columns = distance,
                clustering_method_columns   = method
  )
  
  #Cut tree
  hc <- as.hclust(column_dend(pl))
  group_num <- cutree(tree = hc,
                      k    = k,
                      h    = NULL
  )
  group <- melt(group_num)
  colnames(group) <- "group_num"
  
  # Rename groups
  group$group_name <- ifelse(group$group_num == "1",group_names[1],
                             ifelse(group$group_num == "2",group_names[2],
                                    group_names[3]
                             )
  )
  
  # Add assigned group_num and group_name to sample table
  sampleTable[,group_num_column]  <- group$group_num
  sampleTable[,group_name_column] <- group$group_name
  
  if (silhouette)
  {
    # Compute silhouette
    silhouette <- silhouette_func(matrix,
                                  distance  = distance,
                                  method    = method,
                                  k         = 3
    )
    
    # Add silhouette
    sampleTable$silhouette <- silhouette
  }
  
  return(sampleTable)
}
