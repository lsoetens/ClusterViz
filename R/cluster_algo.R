#' Clustering algorithm using time, geographic space and genetic sequence information (Ypma et al 2013)
#' 
#' @author J.Backer and L.Soetens
#' @param data 
#' @param time.col 
#' @param x.col 
#' @param y.col 
#' @param gen.col 
#' @param parallel 
#'
#' @import ape
#' @import parallel
#' 
#' @return clustering algortihm returns list with clusterids, clusterheights and pvalues
#' @export
cluster.algo<- function(data, time.col, x.col, y.col, gen.col, parallel= FALSE){
        
        # pairwise absolute distances
        
        # time absolute distance matrix (D)
        data$time2<- as.Date(data[, time.col], origin = "1970-01-01")
        time.col<- grep("time2", colnames(data))
        D_time<- dist(data[, time.col], method = "manhattan", diag = TRUE, upper = TRUE)
        D_time<- as.matrix(D_time)
        
        # gen distance matrix (D)
        D_gen<- dist.gene(data[, gen.col], method = "pairwise")
        D_gen<- as.matrix(D_gen)
        
        # geo euclidean distance matrix (D)
        D_geo<- dist(data[, c(x.col, y.col)], method = "euclidean", diag = TRUE, upper = TRUE)
        D_geo<- as.matrix(D_geo)
        
        
        
        # pairwise relative dissimilaritites
        
        d_time <- diss.relative(D_time)
        d_geo <- diss.relative(D_geo)
        d_gen <- diss.relative(D_gen)
        d_combi <- d_time * d_geo * d_gen
        
        
        ############# Tree from dissimilarity matrix ###########################
        
        
        make.timed.mst <- function(d_combi2) {
                n.samples <- nrow(d_combi2)
                d_combi2[d_combi2 <= 0] <- Inf
                infectors <- rep(0, n.samples)
                infected <- c()
                
                max.duplicates <- max(table(d_combi2[d_combi2 != Inf]))
                ordered.by.distance <- order(d_combi2)
                ordered.by.distance <- ordered.by.distance[which(d_combi2[ordered.by.distance]!=Inf)]
                
                while(length(ordered.by.distance) > 0) {
                        candidates <- ordered.by.distance[1:min(length(ordered.by.distance),max.duplicates)]
                        min.dist <- min(d_combi2[candidates])
                        candidates <- candidates[d_combi2[candidates]==min.dist]
                        
                        tos <- candidates %% n.samples
                        tos[tos == 0] <- n.samples
                        froms <- ceiling(candidates/n.samples)  
                        
                        min.pos <- which.max(sapply(tos, function(i) min(d_combi2[i, d_combi2[i,] != min.dist])))
                        to <- tos[min.pos]
                        from <- froms[min.pos]
                        p <- from
                        while(p!=0 & p!=to) p <- infectors[p]
                        if(p==to) {
                                d_combi2[from, to] <- Inf
                                d_combi2[to, from] <- Inf
                        } else {
                                infectors[to] <- from
                                infected <- c(infected, candidates[min.pos])
                                ordered.by.distance <- ordered.by.distance[ordered.by.distance != candidates[min.pos]]
                                d_combi2[to, (1:n.samples)[-from]] <- Inf
                                d_combi2[from, to] <- Inf
                        }
                        ordered.by.distance <- ordered.by.distance[which(d_combi2[ordered.by.distance]!=Inf)]
                }
                
                d_combi2[d_combi2 == Inf] <- 0
                # make the matrix symmetric again
                d_combi2 <- d_combi2 + t(d_combi2)
                #print(sum(d_combi2))
                d_combi2[d_combi2 == 0] <- 10*max(d_combi2)
                return(d_combi2)
        }
        
        timed.mst <- make.timed.mst(d_combi)
        
        
        # clustering without timing
        clust_tree <- hclust(as.dist(d_combi), method = "single")
        
        # determine for each node (= cluster) the height and which members
        all.tree.heights <- clust_tree$height
        all.clusters <- list()
        for(i in 1:nrow(clust_tree$merge)) {
                set <- clust_tree$merge[i,]
                if(all(set < 0)) {
                        all.clusters[[i]] <- -set
                } else if(all(set > 0)) {
                        all.clusters[[i]] <- c(all.clusters[[set[1]]], all.clusters[[set[2]]])
                } else {
                        all.clusters[[i]] <- c(-set[set < 0], all.clusters[[set[set > 0]]])
                }
        }  
        
        # determine unique clusters (their height and which members)
        clusters <- list(all.clusters[[1]])
        tree.heights <- all.tree.heights[1]
        c <- 1
        for(i in 2:length(all.clusters)) {
                if(all.tree.heights[i] != tree.heights[c] || any(!(clusters[[c]] %in% all.clusters[[i]]))) {
                        c <- c+1
                        tree.heights <- c(tree.heights, all.tree.heights[i])
                }
                clusters[[c]] <- sort(all.clusters[[i]])
        }
        
        cluster.sizes <- sapply(clusters, length)
        
        
        ############# Random trees ###########################
        
        
        # do the hustle
        
        hustle.tree <- function(p, use.timing = FALSE, dgen, dtime, dgeo, treeheights, clustersizes) {
                n.samples <- nrow(dgen)
                new.order <- sample(1:n.samples)
                new.dgen <- dgen[new.order, new.order]
                new.order <- sample(1:n.samples)
                new.dtime <- dtime[new.order, new.order]
                new.order <- sample(1:n.samples)
                new.dgeo <- dgeo[new.order, new.order]
                new.dcombi <- new.dtime * new.dgeo * new.dgen
                if(use.timing) new.dcombi <- make.timed.mst(new.dcombi)
                
                new.tree <- hclust(as.dist(new.dcombi), method = "single")
                new.cut.tree <- cutree(new.tree, h = treeheights)
                new.cluster.sizes <- sapply(1:ncol(new.cut.tree), function(i) max(table(new.cut.tree[,i])))
                
                return(new.cluster.sizes >= clustersizes)
        }
        
        
        if (parallel == TRUE){
        
        # parallel hustle
        
        n.cores <- detectCores()
        cl <- makeCluster(n.cores-1)
        
        clusterExport(cl, varlist = c("make.timed.mst", "hustle.tree", "d_gen", "d_time", "d_geo", "tree.heights", "cluster.sizes"))
        
        
        res.hustle <- parSapply(cl,
                        X = 1:10000,
                        FUN = hustle.tree,
                        use.timing = FALSE,
                        dgen = d_gen, 
                        dtime = d_time, 
                        dgeo = d_geo, 
                        treeheights = tree.heights, 
                        clustersizes = cluster.sizes)
        
        
        stopCluster(cl)
        
       
        
        } else {
                
        res.hustle <- replicate(10000, hustle.tree(1, use.timing = FALSE, dgen = d_gen, dtime = d_time, dgeo = d_geo, treeheights = tree.heights, clustersizes = cluster.sizes))
                
        }
        
        
        p_values<- apply(res.hustle, MARGIN = 1, sum)/10000    
        
        res<- list(clusters = clusters, p.values = p_values, tree.heights = tree.heights)
        
}




