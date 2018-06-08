assign.clusters <- function(data, cluster.list, cluster.pvalues, cluster.height, p.limit = 0.001, max.size = 1, max.tree.height = max(cluster.height)) {
        clustersizes <- unlist(lapply(cluster.list, length))
        clusters_p <- cluster.list[which(clustersizes < max.size*nrow(data) & cluster.pvalues < p.limit & cluster.height < max.tree.height)]
        clusters_p<- clusters_p[order(sapply(clusters_p,length),decreasing=F)]
        # remove nested clusters
        clusters_unnested <- list()
        for(i in 1:length(clusters_p)) {
                c <- unlist(clusters_p[i]) %in% unlist(clusters_p[i+1:length(clusters_p)])
                if(all(c) == FALSE) {
                        clusters_unnested[[length(clusters_unnested)+1]] <- clusters_p[i]
                }
        }
        assigned.clusters <- rep(0, nrow(data))
        for (i in 1:length(clusters_unnested)) {
                assigned.clusters[unlist(clusters_unnested[i])] <- i
        }
        return(assigned.clusters)  
}
