gvarclust <- function(
  data,
  vars,
  idvar,
  beepvar,
  dayvar,
  search = TRUE,
  scale = TRUE,
  maxit.in,
  centerWithin = TRUE,
  gamma = 0.50, # Gamma used in glasso in qgraph
  verbose = TRUE,
  redo = FALSE,
  subjectNetworks = TRUE, # Or a vector of which subjects to use!
  lambda_min_kappa_fixed = 0.001,
  lambda_min_beta_fixed = 0.001,
  lambda_min_kappa = 0.05,
  lambda_min_beta = lambda_min_kappa,
  lambda_min_glasso = 0.01, 
  dir = getwd(),
  btwnNetworks = FALSE,
  mlgvclust = TRUE,
  indplot = FALSE,
  elim.groups = FALSE,
  g.ID = NULL,
  ...){
  dir.create(file.path(paste(dir, '/graphical VAR/', sep = '')), showWarnings = F)
  dir.create(file.path(dir, paste('/graphical VAR/Individual Solutions/', sep = '')), showWarnings = F)
  dir.create(file.path(dir, paste('/graphical VAR/Group Solution/', sep = '')), showWarnings = F)
  dir.create(file.path(dir, paste('/graphical VAR/Subgroup Solutions/', sep = '')), showWarnings = F)
  # Determining sample size from elongated data
    dim <- length(unique(data[,idvar]))
    # if(!is.null(g.ID)){
    #   # grups = cbind(data[,c(ID, g.ID)])
    #   grups = data[diff(c(0,ID)) != 0,g.ID]
    #   write.csv(grups, paste0(dir, '/graphical VAR/check.csv', sep = ''), row.names = FALSE)
    # }
  message('Beginning gVAR Clustering. You may want to go grab a coffee because this will take a while.')
  tryCatch({
  group.output <- graphicalVAR::mlGraphicalVAR(data = data,
                                               vars = vars,
                                               beepvar = beepvar,
                                               dayvar = dayvar,
                                               idvar = idvar,
                                               scale = scale,
                                               centerWithin = centerWithin,
                                               gamma = gamma, 
                                               verbose = verbose,
                                               subjectNetworks = subjectNetworks, 
                                               btwnNetworks = FALSE,
                                               lambda_min_kappa_fixed = lambda_min_kappa_fixed,
                                               lambda_min_beta_fixed = lambda_min_beta_fixed,
                                               lambda_min_kappa = lambda_min_kappa,
                                               lambda_min_beta = lambda_min_kappa,
                                               lambda_min_glasso = lambda_min_glasso,
                                               mlgvclust = TRUE,
                                               ...)
  }, error = function(e){})
  # Saving out the fixed effects and the between-subjects model for the sample-level model
  write.table(group.output$fixedPCC, file = paste(dir, '/graphical VAR/Group Solution/fixedPCC.csv', sep = ''), 
              sep = ',', row.names = F, col.names = F)
  write.table(group.output$fixedPDC, file = paste(dir, '/graphical VAR/Group Solution/fixedPDC.csv', sep = ''), 
              sep = ',', row.names = F, col.names = F)
  write.table(group.output$fixedbeta, file = paste(dir, '/graphical VAR/Group Solution/fixedbeta.csv', sep = ''), 
              sep = ',', row.names = F, col.names = F)
  write.table(group.output$fixedkappa, file = paste(dir, '/graphical VAR/Group Solution/fixedkappa.csv', sep = ''), 
              sep = ',', row.names = F, col.names = F)
  if(btwnNetworks == TRUE){
    write.table(group.output$betweenNet, file = paste(dir, '/graphical VAR/Group Solution/BTW Net.csv', sep = ''), 
                  sep = ',', row.names = F, col.names = F)
  }else(message('No Between Subjects networks were generated. I hope you did this intentionally!'))
      group.plot.pdf = qgraph::qgraphMixed(group.output$fixedPCC, group.output$fixedPDC,
                                           parallelEdge = TRUE, theme = 'gimme',
                                           ltyUndirected = 1, ltyDirected = 2,
                                           labels = vars, nNodes = length(vars),
                                           DoNotPlot = TRUE, layout = 'circle')
    pdf(file.path(paste(dir, '/graphical VAR/Group Solution/Group Plot.pdf', sep = '')))
    plot(group.plot.pdf)
    dev.off()
  if(indplot == FALSE){message = "No Individual Plots Generated; Change indplot to TRUE"}  
  for(i in 1:dim){
    write.table(group.output$subjectPCC[[i]], file = paste(dir, '/graphical VAR/Individual Solutions/Subj', i, 'PCC.csv', sep = ''),
                sep = ',', row.names = F, col.names = F)
    write.table(group.output$subjectPDC[[i]], file = paste(dir, '/graphical VAR/Individual Solutions/Subj', i, 'PDC.csv', sep = ''),
                sep = ',', row.names = F, col.names = F)
    write.table(group.output$subjectbeta[[i]], file = paste(dir, '/graphical VAR/Individual Solutions/Subj', i, 'beta.csv', sep = ''),
                sep = ',', row.names = F, col.names = F)
    write.table(group.output$subjectkappa[[i]], file = paste(dir, '/graphical VAR/Individual Solutions/Subj', i, 'kappa.csv', sep = ''),
                sep = ',', row.names = F, col.names = F)
    if (indplot == TRUE) {
      individual.plot.pdf = qgraph::qgraphMixed(group.output$subjectPCC[[i]], 
                                                group.output$subjectPDC[[i]], 
                                                theme = 'gimme', fade = F, 
                                                layout = 'circle', parallel = T, 
                                                ltyUndirected = 1, ltyDirected = 2,
                                                labels = vars)
    pdf(file.path(dir, paste('/graphical VAR/Individual Solutions/Subj', i,'Plot.pdf', sep = '')))
    plot(individual.plot.pdf)
    dev.off()
    }
  }
  message('Beginning Community Detection')
  # Defining N x N matrices: Similarity matrices for lagged and contemporaneous networks
    c.sim = l.sim = matrix(0, max(dim), max(dim))
    cs.sim = ls.sim = list(matrix(NA, length(vars), length(vars)))
    for(i in 1:dim){
      cs.sim[[i]] = ls.sim[[i]] = matrix(NA, length(vars), length(vars))
    }
  # Compare to groups
  # if(!is.null(g.ID)){
  #   for(i in 1:dim){
  #     for(j in 1:dim){
  #       if(grups[i] == grups[j]){c.sim[i,j] = c.sim[i,j] + 10}else{c.sim[i,j] = 0}
  #     }
  #   }
  # }
  if(elim.groups == TRUE){
    for(i in 1:dim){
      for(j in 1:length(vars)){
        for(k in 1:length(vars)){
          if(abs(group.output$subjectPCC[[i]][j,k]) > 0 & abs(group.output$fixedPCC[j,k]) > 0){
            cs.sim[[i]][j,k] = group.output$subjectPCC[[i]][j,k]}else{
              cs.sim[[i]][j,k] = 0
              }
          if(abs(group.output$subjectPDC[[i]][j,k]) > 0 & abs(group.output$fixedPDC[j,k]) > 0){
            ls.sim[[i]][j,k] = group.output$subjectPDC[[i]][j,k]}else{
              ls.sim[[i]][j,k] = 0
              }
        }
      }
    }
    for(i in 1:dim){
      for(j in 1:dim){
        c.sim[i,j] <- sum(abs(c(cs.sim[[i]])) > 0 & abs(c(cs.sim[[j]])) > 0 &
                            sign(c(cs.sim[[i]])) == sign(c(cs.sim[[j]])))
        l.sim[i,j] <- sum(abs(c(ls.sim[[i]])) > 0.01 & abs(c(ls.sim[[j]])) > 0.01 &
                            sign(c(ls.sim[[i]])) == sign(c(ls.sim[[j]])))
      }
    }
  }else if(elim.groups == FALSE || is.null(elim.groups)){
    # Generate the N x N matrices
    for(i in 1:dim){
      for(j in 1:dim){
        c.sim[i,j] <- sum(abs(c(group.output$subjectPCC[[i]])) > 0 & abs(c(group.output$subjectPCC[[j]])) > 0 &
                            sign(c(group.output$subjectPCC[[i]])) == sign(c(group.output$subjectPCC[[j]])))
        l.sim[i,j] <- sum(abs(c(group.output$subjectPDC[[i]])) > 0.01 & abs(c(group.output$subjectPDC[[j]])) > 0.01 &
                            sign(c(group.output$subjectPDC[[i]])) == sign(c(group.output$subjectPDC[[j]])))
      }
    }
  }
  

    # for(i in 1:dim){
    #   for(j in 1:dim){
    #     c.sim[i,j] <- sum(abs(c(group.output$subjectPCC[[i]])) > 0 & abs(c(group.output$subjectPCC[[j]])) > 0)
    #     cs.sim[i,j] <- sum(sign(c(group.output$subjectPCC[[i]])) == sign(c(group.output$subjectPCC[[j]])))
    #     l.sim[i,j] <- sum(abs(c(group.output$subjectPDC[[i]])) > 0.05 & abs(c(group.output$subjectPDC[[j]])) > 0.05)
    #     ls.sim[i,j] <- sum(sign(c(group.output$subjectPDC[[i]])) == sign(c(group.output$subjectPDC[[j]])))
    #   }
    # }
  # Generate similarity matrix | The sum of lagged & contemporaneous similarity matrices
    # sim.mat <- c.sim + l.sim
    sim.mat <- c.sim + l.sim
    diag(sim.mat) = 0
    sim.mat = as.matrix(sim.mat)
    write.table(sim.mat, paste(dir, '/sim.mat.test.csv', sep = ''), sep = ',')
    jpmodmax1 <- function(x, m){
      # sparsify m based on x
      m[which(m <= quantile(m[upper.tri(m, diag = FALSE)], x))] = 0
      p = igraph::cluster_walktrap(igraph::graph_from_adjacency_matrix(m))
      modval <- igraph::modularity(p)
      sk = ask = list()
      mem = p$membership
      dim = length(mem)

      for(i in 1:dim){
        sk[i] = 0
        ask[i] = 0
        for(j in 1:dim){
          if(m[i,j] == 0){NULL}else{
            if(mem[i] == mem[j]){ask[[i]] = ask[[i]] + m[i,j]}else{
              sk[[i]] = sk[[i]] + m[i,j]
            }
          }
        }
        sk[[i]] = sk[[i]]
        ask[[i]] = sk[[i]] + ask[[i]]
      }
      sk = unlist(sk); ask = unlist(ask)

      temp.1 = cbind(mem, sk, ask)
      c_list = list()
      for(i in 1:max(temp.1[,1])){
        if(temp.1[mem == i,2] == 0 && temp.1[mem == i,3] == 0){
          # c_list[[i]] = sum(temp.1[group_membership == i,2]) / 
          #   (sum(temp.1[,3]) - sum(temp.1[group_membership == i,3]))
          c_list[[i]] = 1
        }else{
          c_list[[i]] = sum(temp.1[mem == i,2]) / 
            sum(temp.1[mem == i,3])  
        }
      }
      simplified = unlist(c_list)
      conductance = 1 - (sum(simplified)/length(simplified))
      # for(i in 1:max(temp.1[,1])){
      #   if(nrow(temp.1[mem == i,]) <= 1 || is.null(nrow(temp.1[mem == i,]))){
      #     c_list[[i]] = sum(temp.1[mem == i,2]) / (sum(temp.1[,3]) - sum(temp.1[mem == i,3]))
      #   }else{
      #     c_list[[i]] = sum(temp.1[mem == i,2]) / sum(temp.1[mem == i,3])
      #   }
      # }
      
      # conductance = 1 - sum(unlist(c_list))/max(temp.1[,1])
      return(-1*conductance)
      #
      #
      #
      
      # return(-1*modval) # solvers generally want to minimze the
      # return(-1*rat) # solvers generally want to minimze the
    }
    if(search == TRUE){
    res <- nloptr::nloptr(
      x0 = 0.5,  # starting value
      m = sim.mat,  # similarity matrix
      eval_f = jpmodmax1, # objective function
      lb = .01, # upper bound
      ub = .99, # lower bound
      opts = list(
        "algorithm"="NLOPT_GN_DIRECT_L",
        "ftol_rel"=1.0e-8,
        maxeval = 1000)
    )
    sim.mat[which(sim.mat <= quantile(sim.mat[upper.tri(sim.mat, diag = FALSE)], (res$solution)))] <- 0
    write.table(res$solution, paste(dir, '/modularity_val.csv', sep = ''), sep = ',')
    }else if(search == FALSE){
      sim.mat = as.matrix(sim.mat)
      sim.mat = sim.mat - min(sim.mat, na.rm = TRUE)
    }
    diag(sim.mat) = 0
    sim.mat = as.matrix(sim.mat)
    lay.sim = sim.mat
    write.table(lay.sim, paste(dir, '/sim.mat.csv', sep = ''), sep = ',')
  # Use igraph to prep dat for community detection
    g = igraph::graph.adjacency(sim.mat, mode = "undirected", weighted = TRUE)
    weights = igraph::E(g)$weight
  # Apply walktrap for community detection of the graphical VAR model
    res = igraph::cluster_walktrap(g, weights = weights, steps = 4)
    # res$membership
    if(search == FALSE){
      write.table(igraph::modularity(res), paste(dir, '/modularity_val.csv', sep = ''), sep = ',')
    }
    e = igraph::get.edgelist(igraph::graph.adjacency(lay.sim, mode = "undirected", weighted = TRUE))
    l = qgraph::qgraph.layout.fruchtermanreingold(e, vcount = igraph::vcount(g), weights = weights/5,
                                                  area = 8*(igraph::vcount(g)^2),
                                                  repulse.rad = (igraph::vcount(g)^3.1))
    
    par(mfrow = c(1, 1), mar = c(0, 0, 0, 0))
    pdf(file.path(paste(dir,'/graphical VAR/Subgroups Plot.pdf', sep = '')))
    plot(res, g, layout = l)
    dev.off() 
  group.membership = as.matrix(cbind(unique(data[,idvar]), 
                                  res$membership), max(dim), 2)
  write.table(group.membership, file = paste(dir, '/graphical VAR/Subgroup Solutions/Group Membership.csv', sep = ''), 
                  sep = ',', row.names = F, col.names = F)
  groups = sum(colSums(table(data.frame(group.membership))) > 1)
  if(groups == 0){
    message("No subgroups detected.")
  }else{
  for(i in 1:groups){
      dir.create(file.path(dir, paste('/graphical VAR/Subgroup Solutions/Subgroup ', i,'/', sep = '')), 
                 showWarnings = F)  
    }
  for(i in 1:groups){
    message(paste('Network Estimation for Subgroup ', i, sep = ''))
    subsetting = group.membership[group.membership[,2] == i, 1]
    # sub.dat = subset(as.data.frame(data), ID %in% subsetting)
    sub.dat = subset(as.data.frame(data), select = c(vars, idvar))
    sub.dat = sub.dat[sub.dat[,ncol(sub.dat)] %in% subsetting,]
    # subgrouping <- as.data.frame(tsData(dat = sub.dat, vars = vars, idvar = idvar,
    #                      scale = scale, centerWithin = centerWithin)$dat)
    subgroup = graphicalVAR::mlGraphicalVAR(data = sub.dat, vars = vars, idvar = idvar, #beepvar = beepvar,
                              subjectNetworks = FALSE, scale = scale, #dayvar = dayvar,
                              lambda_min_kappa = lambda_min_kappa,
                              lambda_min_beta = lambda_min_beta,
                              centerWithin = centerWithin, gamma = gamma, verbose = verbose,
                              btwnNetworks = FALSE, mlgvclust = TRUE, ...)
      write.table(subgroup$fixedPCC, file = paste(dir, '/graphical VAR/Subgroup Solutions/Subgroup ', i, '/fixedPCC.csv',sep = ''),
                  sep = ',', row.names = F, col.names = F)
      write.table(subgroup$fixedPDC, file = paste(dir, '/graphical VAR/Subgroup Solutions/Subgroup ', i, '/fixedPDC.csv',sep = ''),
                  sep = ',', row.names = F, col.names = F)
      write.table(subgroup$fixedbeta, file = paste(dir, '/graphical VAR/Subgroup Solutions/Subgroup ', i, '/fixedbeta.csv',sep = ''),
                  sep = ',', row.names = F, col.names = F)
      write.table(subgroup$fixedkappa, file = paste(dir, '/graphical VAR/Subgroup Solutions/Subgroup ', i, '/fixedkappa.csv',sep = ''),
                  sep = ',', row.names = F, col.names = F)
      # write.table(btwNetwork$optnet, file = paste(dir, '/graphical VAR/Subgroup Solutions/Subgroup ', i, '/BTW Net.csv',sep = ''),
      #             sep = ',', row.names = F, col.names = F)
      subgroup.plot.pdf = qgraph::qgraphMixed(subgroup$fixedPCC, subgroup$fixedPDC,
                                              parallelEdge = TRUE, theme = 'gimme',
                                              ltyUndirected = 1, ltyDirected = 2,
                                              labels = vars, nNodes = length(vars),
                                              DoNotPlot = TRUE, layout = 'circle')
    pdf(file.path(dir, paste('graphical VAR/Subgroup Solutions/Subgroup ', i,'/Subgroup Plot.pdf', sep = '')))
    plot(subgroup.plot.pdf)
    dev.off()
  }
  }
}