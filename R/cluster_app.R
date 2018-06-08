
#' Cluster_app function
#'
#' @param data 
#' @param cluster.list 
#' @param cluster.pvalue 
#' @param cluster.height 
#' @param id.col 
#' @param time.col 
#' @param x.col 
#' @param y.col 
#' @param gen.col 
#' @param cols 
#'
#' @import shiny
#' @import ggplot2
#' @import RColorBrewer
#' @import DT
#' @import ape
#' @import phangorn
#' @import ggtree
#' @import igraph
#' @import ggraph
#' @import scales
#' @import leaflet
#' @import plotly
#'
#' @return An shiny application visualizing possilbe transmission clusters
#' @export 
cluster_app<- function(data, cluster.list, cluster.pvalue, cluster.height, id.col, time.col, x.col, y.col, gen.col, cols){
        
        
        
        # Define UI for application 
        ui <- fluidPage(
                headerPanel('Visualization of possible transmission clusters in time, geospatial and genetic space'),
                fluidRow(
                        column(2, numericInput("plimit", label = h4("P-value limit"), value = 0.05), numericInput("size", label = h4("Clustersize limit of total sample"), value = 100),
                                numericInput("treelimit", label = h4("Treeheight limit (proportion of max height"), value = 1), submitButton('Update'),
                                h4("Cluster information"), dataTableOutput('table')),
                        column(3, h4("Tree heights vs p-values for all clusters"), plotOutput('plot1', height = 300), h4("Hierarchical clustering tree based on combined pairwise distances"), plotOutput('plot2', height = 300), h4("Epidemic curve"), plotOutput('plottime', height = 200)),
                        column(3, h4("Map with clusters"), leafletOutput("plot3", height= 400), h4("ML phylogenetic tree (arbitrarily rooted)"), plotOutput('plotgen', height= 400)),
                        column(3, h4("Pairwise inter-patient distance per cluster and dimension"), plotOutput('plot4', height = 400), h4("Intra-cluster correlation between dimensions"), plotlyOutput('plot5', height = 400))
                )
        )
        # Define server logic 
        server <- function(input, output) {
                
                output$table <- renderDataTable({
                        
                        data$clusters <- assign.clusters(data = data, cluster.list = cluster.list, cluster.pvalue = cluster.pvalue, 
                                cluster.height = cluster.height, p.limit = input$plimit, max.size = input$size, max.tree.height = (input$treelimit * max(cluster.height)))
                        
                        x<- as.data.frame(table(data$clusters))
                        names(x)<- c("Cluster", "Clustersize")
                        x<- datatable(x) %>% formatStyle(
                                'Cluster',
                                backgroundColor = styleEqual(c(unique(data$clusters)), c(cols[1:length(unique(data$clusters))])))
                        x
                })
                
                # time absolute distance matrix (D)
                D_time<- dist(data[, time.col], method = "manhattan", diag = TRUE, upper = TRUE)
                D_time<- as.matrix(D_time)
                
                # gen distance matrix (D)
                D_gen<- dist.gene(data[, gen.col], method = "pairwise")
                D_gen<- as.matrix(D_gen)
                
                # geo euclidean distance matrix (D)
                D_geo<- dist(data[, c(x.col, y.col)], method = "euclidean", diag = TRUE, upper = TRUE)
                D_geo<- as.matrix(D_geo)
                
                
                # pairwise dissimilaritites
                
                d_time <- diss.relative(D_time)
                d_geo <- diss.relative(D_geo)
                d_gen <- diss.relative(D_gen)
                d_combi <- d_time * d_geo * d_gen
                
                
                
                output$plot1 <- renderPlot({
                        clusters_p_size <- unlist(lapply(cluster.list[which(unlist(lapply(cluster.list, length)) < input$size*nrow(data) & cluster.pvalue < input$plimit & cluster.height < (input$treelimit * max(cluster.height)))], length))
                        clusters_p_heights <- cluster.height[which(unlist(lapply(cluster.list, length)) < input$size*nrow(data) & cluster.pvalue < input$plimit & cluster.height < (input$treelimit * max(cluster.height)))]
                        clusters_p_pvalues <- cluster.pvalue[which(unlist(lapply(cluster.list, length)) < input$size*nrow(data) & cluster.pvalue < input$plimit & cluster.height < (input$treelimit * max(cluster.height)))]
                        clusters_p<- data.frame(clusters_p_pvalues, clusters_p_heights)
                        
                        #scatterplot of x and y variables
                        scatter <- ggplot(clusters_p, aes(x= clusters_p_heights, y= clusters_p_pvalues))+
                                geom_point()+
                                xlab("treeheight")+
                                ylab("p-value")+
                                theme_bw()+
                                theme(legend.position = "none") 
                        
                        scatter
                })
                
                output$plot2 <- renderPlot({
                        
                        data$clusters <- assign.clusters(data = data, cluster.list = cluster.list, cluster.pvalue = cluster.pvalue, 
                                cluster.height = cluster.height, p.limit = input$plimit, max.size = input$size, max.tree.height = (input$treelimit * max(cluster.height)))
                        
                        clust_tree <- hclust(as.dist(d_combi), method = "single")
                        dend <- as.dendrogram(clust_tree)
                        
                        dendrogram.layout<- create_layout(den_to_igraph(dend, even = FALSE), layout="dendrogram")
                        dendrogram.layout<- dendrogram.layout[,5:9]
                        cluster_output<- data.frame(pvalues= cluster.pvalue, clustersize= unlist(lapply(cluster.list, length)), heights= cluster.height)
                        # keep only significant clusters
                        cluster_output<- subset(cluster_output, pvalues< input$plimit)
                        # merge with dendrogram layout
                        dendrogram.layout.2<- merge(dendrogram.layout, cluster_output, by.x= c("layout.y", "members"), by.y= c("heights", "clustersize"), all.x =T)
                        
                        # add cluster- number to layout dataframe
                        label_dend<- as.numeric(labels(dend))
                        cluster<- as.factor(data$clusters)[as.numeric(labels(dend))]
                        labeldf<- data.frame(label_dend, cluster)
                        dendrogram.layout.2<- merge(dendrogram.layout.2, labeldf, by.x= "label", by.y= "label_dend", all.x = T)
                        
                        ggraph(dend, 'dendrogram') +
                                geom_edge_elbow(edge_width = 0.4, edge_colour = "darkgrey") +
                                geom_node_point(aes(x = layout.x, y = layout.y), col="black", size = 1.5, data= subset(dendrogram.layout.2, pvalues< input$plimit))+
                                geom_node_point(aes(x = layout.x, y = layout.y, col = cluster), size = 2, shape = 17, data= subset(dendrogram.layout.2, leaf == TRUE)) +
                                scale_color_manual(values = cols)+
                                geom_hline(yintercept = (input$treelimit * max(cluster.height)), size = 0.3, linetype="dashed")+
                                ylab("Combined dissimilarity (d_combi)")+
                                #coord_cartesian(ylim = c(0, 650))+
                                theme_bw()+
                                theme(axis.line.x = element_blank(),
                                        axis.title.x = element_blank(),
                                        axis.text.x = element_blank(),
                                        axis.ticks.x = element_blank(),
                                        legend.position = "none",
                                        text = element_text(size = 10),
                                        axis.text.y = element_text(size = 8),
                                        panel.grid.major = element_blank(), 
                                        panel.grid.minor = element_blank(),
                                        panel.background = element_blank(),
                                        axis.line.y = element_line(size = 0.3))
                        
                })
                
                output$plottime<- renderPlot({
                        data$clusters <- assign.clusters(data = data, cluster.list = cluster.list, cluster.pvalue = cluster.pvalue, 
                                cluster.height = cluster.height, p.limit = input$plimit, max.size = input$size, max.tree.height = (input$treelimit * max(cluster.height)))
                        
                        p<- ggplot(data, aes(x = data[,time.col], fill = as.factor(clusters)))+
                                geom_histogram(binwidth = 30) +
                                scale_fill_manual(values = colors, name= " ")+
                                scale_y_continuous(breaks= pretty_breaks())+
                                xlab("Date onset disease")+
                                ylab("Number of cases")+
                                theme_bw()+
                                theme(axis.title = element_text(size = 10),
                                        axis.text = element_text(size = 8),
                                        legend.position = "none")
                        p
                })
                
                output$plot3 <- renderLeaflet({
                        data$clusters <- assign.clusters(data = data, cluster.list = cluster.list, cluster.pvalue = cluster.pvalue, 
                                cluster.height = cluster.height, p.limit = input$plimit, max.size = input$size, max.tree.height = (input$treelimit * max(cluster.height)))
                        
                        map_plot<- function(data) {
                                leaflet(data = data) %>% addTiles() %>%
                                        addCircleMarkers(~data[,x.col], ~data[,y.col], popup = ~as.character(id.col),
                                                color = ~cols[as.factor(clusters)],
                                                stroke = FALSE, fillOpacity = 0.5
                                        )
                        }
                        
                        map<- map_plot(data = data)
                        map
                        
                })
                
                output$plotgen<- renderPlot({
                        ddata$clusters <- assign.clusters(data = data, cluster.list = cluster.list, cluster.pvalue = cluster.pvalue, 
                                cluster.height = cluster.height, p.limit = input$plimit, max.size = input$size, max.tree.height = (input$treelimit * max(cluster.height)))
                        
                        # preprocessing of gen data -> ml tree
                        gen.dat<- as.matrix(data[,gen.col])
                        rownames(gen.dat)<- data[,1]
                        gen.dat <- phyDat(gen.dat, return.index = TRUE)
                        nj.tree <- NJ(dist.ml(gen.dat, model = "JC69"))
                        ml.tree <- pml(tree = nj.tree, data = gen.dat)
                        ml.tree <- optim.pml(ml.tree) 
                        
                        phylo_plot <- function(data, clusters, cols, tree) {
                                ggtree(tree)+ 
                                        geom_treescale(fontsize = 4, offset = 2.5)+ 
                                        geom_tippoint(color= cols[as.factor(data[, clusters])], shape= 20, size=1.5)+
                                        theme(plot.margin = unit(c(0,0,0,0), "lines"))
                        }
                        
                        phylo<- phylo_plot(data = data, clusters= "clusters", cols=cols, tree= ml.tree$tree)
                        phylo
                })
                
                output$plot4 <- renderPlot({
                        
                        data$clusters <- assign.clusters(data = data, cluster.list = cluster.list, cluster.pvalue = cluster.pvalue, 
                                cluster.height = cluster.height, p.limit = input$plimit, max.size = input$size, max.tree.height = (input$treelimit * max(cluster.height)))
                        
                        df_dist <- data.frame(d= numeric(0), dim= character(0), cluster= numeric(0) )
                        for (i in sort(unique(data[,"clusters"]))) {
                                id_temp <- which(data[,"clusters"]==i)
                                offdiag <- !diag(length(id_temp))
                                n.elem <- length(id_temp)*(length(id_temp)-1)
                                df_temp <- data.frame(d = c(d_time[id_temp, id_temp][offdiag], d_geo[id_temp, id_temp][offdiag], d_gen[id_temp, id_temp][offdiag], d_combi[id_temp, id_temp][offdiag]^(1/3)),
                                        dim = c(rep("time", n.elem), rep("geo", n.elem), rep("gen", n.elem), rep("combi", n.elem)),
                                        cluster = c(rep(i, 4*n.elem)))
                                df_dist <- rbind(df_dist, df_temp)
                        }
                        
                        medians<- tapply(df_dist$d[df_dist$dim=="combi"], df_dist$cluster[df_dist$dim=="combi"], median)
                        
                        # Plot interpatient distance per cluster and dimension
                        ggplot(df_dist, aes(x = dim, y= d))+
                                geom_boxplot(aes(x = dim, y= d, fill = as.factor(cluster)), alpha= 0.7, notch = T)+
                                scale_fill_manual(values = cols)+
                                ylab("Interpatient distance")+
                                facet_grid(.~cluster)+
                                theme_bw()+
                                theme(axis.text.x= element_text(angle = 90),
                                        legend.position= "none",
                                        text = element_text(size = 10))
                        
                })
                
                output$plot5 <- renderPlotly({
                        
                        data$clusters <- assign.clusters(data = data, cluster.list = cluster.list, cluster.pvalue = cluster.pvalue, 
                                cluster.height = cluster.height, p.limit = input$plimit, max.size = input$size, max.tree.height = (input$treelimit * max(cluster.height)))
                        
                        df_corr <- data.frame(d_time= numeric(0), d_geo= numeric(0), d_gen= numeric(0), d_combi= numeric(0), cluster= numeric(0) )
                        for (i in sort(unique(data[,"clusters"]))) {
                                id_temp <- which(data[,"clusters"]==i)
                                offdiag <- !diag(length(id_temp))
                                n.elem <- length(id_temp)*(length(id_temp)-1)
                                df_temp <- data.frame(d_time = c(d_time[id_temp, id_temp][offdiag]), d_geo = c(d_geo[id_temp, id_temp][offdiag]),
                                        d_gen = c(d_gen[id_temp, id_temp][offdiag]), d_combi = c(d_combi[id_temp, id_temp][offdiag]),
                                        cluster = c(rep(i, n.elem)))
                                df_corr <- rbind(df_corr, df_temp)
                        }
                        
                        df_corr$cluster<- as.factor(df_corr$cluster)
                        
                        #test significance of correlations
                        cor.mtest <- function(mat, ...) {
                                mat <- as.matrix(mat)
                                n <- ncol(mat)
                                p.mat<- matrix(NA, n, n)
                                diag(p.mat) <- 0
                                for (i in 1:(n - 1)) {
                                        for (j in (i + 1):n) {
                                                tmp <- cor.test(mat[, i], mat[, j], method = "spearman", exact = F, continuity = T)
                                                p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
                                        }
                                }
                                colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
                                p.mat
                        } 
                        
                        
                        t <- list(
                                family = "helvetica",
                                size = 10)
                        
                        pltList <- list()
                        par(mfrow=c(ceiling(sqrt(max(data[,"clusters"]))), ceiling(sqrt(max(data[,"clusters"])))))
                        for (i in 1: max(data[,"clusters"])){
                                pltName <- paste( 'p', i, sep = '' )
                                p.mat.cluster<- cor.mtest(subset(df_corr, cluster == i)[,1:4])
                                p.mat.cluster[lower.tri(p.mat.cluster)]<- NA 
                                c.mat.cluster<- cor(subset(df_corr, cluster == i)[,1:4], method = "spearman")
                                c.mat.cluster[lower.tri(c.mat.cluster)]<- NA
                                x<- c(rep(c("time", "place", "gen", "combi"), 4))
                                y<- c(rep(c("time", "place", "gen", "combi"), each= 4))
                                temp<- data.frame(x= x, y= y, r= c(c.mat.cluster), p=c(p.mat.cluster))
                                p<- plot_ly(temp,
                                        x = ~x, y = ~y, z = ~r, 
                                        type = "heatmap",
                                        colors = colorRamp(c(cols[i+1], "white", cols[i+1])),
                                        zmin = -1, zmax = 1,
                                        text = ~paste('r(spearman): ', round(r, 3), '\np-value: ', round(p,3)),
                                        hoverinfo= "text"
                                )
                                pltList[[ pltName ]]<- p %>% layout(xaxis = list(autorange = "reversed", title = ""), 
                                        yaxis = list(title = "") , font = t) %>% hide_colorbar()
                                
                        }
                        
                        
                        subplot(pltList, nrows=floor(sqrt(length(pltList))), shareX = TRUE, shareY = TRUE)
                        
                })
                
                
                
                
                
        }
        
        # Run the application 
        shinyApp(ui = ui, server = server)
        
}