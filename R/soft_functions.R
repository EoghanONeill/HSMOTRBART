
# The get internal function reports the most important information for each internal node, input: tree
get_internals = function(tree){

  save_ancestor = NULL

  if(is.null(tree) | nrow(tree$tree_matrix) == 1) {
    save_ancestor = cbind(terminal = NULL,
                          split_value = NULL)
  } else {
    which_internal = which(tree$tree_matrix[,'terminal'] == 0)
    for (k in 1:length(which_internal)){
      get_parent = as.numeric(as.character(tree$tree_matrix[which_internal[k],'parent'])) # get the subsequent parent
      get_child_left = as.numeric(as.character(tree$tree_matrix[which_internal[k], 'child_left'])) # get the subsequent parent
      get_child_right = as.numeric(as.character(tree$tree_matrix[which_internal[k], 'child_right'])) # get the subsequent parent
      get_n_left = as.numeric(as.character(tree$tree_matrix[get_child_left, 'node_size'])) # then, get the split value from the associated parent
      get_n_right = as.numeric(as.character(tree$tree_matrix[get_child_right, 'node_size'])) # then, get the split value from the associated parent
      get_split_val = as.numeric(as.character(tree$tree_matrix[which_internal[k], 'split_value'])) # then, get the split value from the associated parent
      get_split_var = as.numeric(as.character(tree$tree_matrix[which_internal[k], 'split_variable'])) # then, get the split variable from the associated parent

      #
      #   get_left = as.integer( tree$tree_matrix[get_parent, 'child_left'] == which_terminal[k] ) # then, determine whether it was a move to the left
      #   get_split_var = as.numeric(as.character(tree$tree_matrix[get_parent, 'split_variable'])) # then, get the split variable from the associated parent
      #
      #
      save_ancestor = rbind(save_ancestor,
                            cbind(internal = which_internal[k],
                                  parent = get_parent,
                                  child_left   = get_child_left,
                                  child_right  = get_child_right,
                                  n_left       = get_n_left,
                                  n_right      = get_n_right,
                                  split_value = get_split_val,
                                  var = get_split_var))
    }
    #   while (get_parent > 1){
    #     get_parent = as.numeric(as.character(tree$tree_matrix[get_parent,'parent'])) # get the subsequent parent
    #     get_split_val = as.numeric(as.character(tree$tree_matrix[get_parent, 'split_value'])) # then, get the split value from the associated parent
    #     get_left = as.integer( tree$tree_matrix[get_parent, 'child_left'] == tree$tree_matrix[which_terminal[k], 'parent'] ) # then, determine whether it was a move to the left
    #     get_split_var = as.numeric(as.character(tree$tree_matrix[get_parent, 'split_variable'])) # then, get the split variable from the associated parent
    #
    #     save_ancestor = rbind(save_ancestor,
    #                           cbind(terminal = which_terminal[k],
    #                                 parent   = get_parent,
    #                                 split_value = get_split_val,
    #                                 left     = get_left,
    #                                 var      = get_split_var))
    #   }
    # }
    save_ancestor = unique(save_ancestor) # remove duplicates
    save_ancestor = save_ancestor[order(save_ancestor[,1], save_ancestor[,2]),] # sort by terminal and ancestor
  }

  return(save_ancestor)
}

# The parent list function lists all the parents corresponding to the internal/terminal node, input: tree, observation number
parent_list = function(mtree, index){
  terminal_index = as.numeric(mtree$node_indices[index])
  parent = as.numeric(mtree$tree_matrix[terminal_index,'parent'])
  parents = cbind(terminal_index,parent)
  while(!is.na(parent)){
    parent = as.numeric(mtree$tree_matrix[parent,'parent'])
    parents = cbind(parents, parent)
  }
  return(parents)
}

# The phi matrix function constructs based on the hierarchical shrinkage method, input: tree, data and information about the internals
phi_matrix = function(mtree,int, x){

  n = nrow(x)
  phi_matrix = NULL
  if(!is.matrix(int)){
    phi_matrix = matrix(NA,n,1)
    for(i in 1:n){
      phi_matrix[i,1] = ( int['n_right']*sum(ifelse(is.element(parent_list(mtree,i),int['child_left']),1,0)) -
                            int['n_left']*sum(ifelse(is.element(parent_list(mtree,i),int['child_right']),1,0)) )/ sqrt(int['n_left']*int['n_right'])
    }
  }
  else{
    phi_matrix = matrix(NA,n,nrow(int))
    for(j in 1:nrow(int)){
      for(i in 1:n){

        phi_matrix[i,j] = (int[j,'n_right']*sum(ifelse(is.element(parent_list(mtree,i),int[j,'child_left']),1,0)) -
                             int[j,'n_left']*sum(ifelse(is.element(parent_list(mtree,i),int[j,'child_right']),1,0)) )/ sqrt(int[j,'n_left']*int[j,'n_right'])

      }
    }

  }

  return(phi_matrix)
}




# fill_tree_details

# The parent list function lists all the parents corresponding to the internal/terminal node, input: tree, observation number
parent_list_test = function(tree_matrix_test, index, node_indices_test){
  terminal_index = as.numeric(node_indices_test[index])
  parent = as.numeric(tree_matrix_test[terminal_index,'parent'])
  parents = cbind(terminal_index,parent)
  while(!is.na(parent)){
    parent = as.numeric(tree_matrix_test[parent,'parent'])
    parents = cbind(parents, parent)
  }
  return(parents)
}

# The phi matrix function constructs based on the hierarchical shrinkage method, input: tree, data and information about the internals
phi_matrix_test = function(mtree,int, x){

  details_test = fill_tree_details(mtree, x)
  node_indices_test = details_test$node_indices
  # tree_matrix_test = details_test$tree_matrix
  tree_matrix_test = mtree$tree_matrix

  # print("tree_matrix_test = ")
  # print(tree_matrix_test)
  #
  # print("mtree = ")
  # print(mtree)
  # print("node_indices_test = ")
  # print(node_indices_test)

  n = nrow(x)
  phi_matrix = NULL
  if(!is.matrix(int)){
    phi_matrix = matrix(NA,n,1)
    for(i in 1:n){
      phi_matrix[i,1] = ( int['n_right']*sum(ifelse(is.element(parent_list_test(tree_matrix_test,i, node_indices_test),int['child_left']),1,0)) -
                            int['n_left']*sum(ifelse(is.element(parent_list_test(tree_matrix_test,i, node_indices_test),int['child_right']),1,0)) )/ sqrt(int['n_left']*int['n_right'])
    }
  }
  else{
    phi_matrix = matrix(NA,n,nrow(int))
    for(j in 1:nrow(int)){
      for(i in 1:n){

        phi_matrix[i,j] = (int[j,'n_right']*sum(ifelse(is.element(parent_list_test(tree_matrix_test,i, node_indices_test),int[j,'child_left']),1,0)) -
                             int[j,'n_left']*sum(ifelse(is.element(parent_list_test(tree_matrix_test,i, node_indices_test),int[j,'child_right']),1,0)) )/ sqrt(int[j,'n_left']*int[j,'n_right'])

      }
    }

  }

  return(phi_matrix)
}

# The variable list function lists the variables that correspond to the splitting of each internal node, input; tree, number of the internal
variable_list = function(tree, internal){
  list_var = 1
  get_parent = as.numeric(as.character(tree$tree_matrix[internal,'parent'])) # get the subsequent parent
  get_split_var = as.numeric(as.character(tree$tree_matrix[get_parent, 'split_variable']))
  list_var = cbind(list_var,get_split_var)

  while(!is.na(get_parent)){
    get_parent = as.numeric(as.character(tree$tree_matrix[get_parent,'parent'])) # get the subsequent parent
    get_split_var = as.numeric(as.character(tree$tree_matrix[get_parent, 'split_variable']))
    list_var = cbind(list_var,get_split_var)
  }

  list_var = unique(as.vector(list_var))
  list_var = list_var[!is.na(list_var)]

  return(list_var)
}

# The get_brand function outputs branching step for each terminal node in a given tree (inspired by the get ancestors function)
# This includes the relevant terminal node, the parents, move to the left, the split variable and value. inputs: tree
get_branch = function(tree){

  save_ancestor = NULL

  if(is.null(tree) | nrow(tree$tree_matrix) == 1) {
    save_ancestor = cbind(terminal = NULL,
                          split_value = NULL)
  } else {
    which_terminal = which(tree$tree_matrix[,'terminal'] == 1)
    for (k in 1:length(which_terminal)){
      get_parent = as.numeric(as.character(tree$tree_matrix[which_terminal[k], 'parent'])) # get the subsequent parent
      get_split_val = as.numeric(as.character(tree$tree_matrix[get_parent, 'split_value'])) # then, get the split value from the associated parent
      get_left = as.integer( tree$tree_matrix[get_parent, 'child_left'] == which_terminal[k] ) # then, determine whether it was a move to the left
      get_split_var = as.numeric(as.character(tree$tree_matrix[get_parent, 'split_variable'])) # then, get the split variable from the associated parent


      save_ancestor = rbind(save_ancestor,
                            cbind(terminal = which_terminal[k],
                                  parent   = get_parent,
                                  split_value = get_split_val,
                                  left     = get_left,
                                  var      = get_split_var))
      while (get_parent > 1){
        get_parent = as.numeric(as.character(tree$tree_matrix[get_parent,'parent'])) # get the subsequent parent
        get_split_val = as.numeric(as.character(tree$tree_matrix[get_parent, 'split_value'])) # then, get the split value from the associated parent
        get_left = as.integer( tree$tree_matrix[get_parent, 'child_left'] == tree$tree_matrix[which_terminal[k], 'parent'] ) # then, determine whether it was a move to the left
        get_split_var = as.numeric(as.character(tree$tree_matrix[get_parent, 'split_variable'])) # then, get the split variable from the associated parent

        save_ancestor = rbind(save_ancestor,
                              cbind(terminal = which_terminal[k],
                                    parent   = get_parent,
                                    split_value = get_split_val,
                                    left     = get_left,
                                    var      = get_split_var))
      }
    }
    save_ancestor = unique(save_ancestor) # remove duplicates
    save_ancestor = save_ancestor[order(save_ancestor[,1], save_ancestor[,2]),] # sort by terminal and ancestor
  }

  return(save_ancestor)
}

# The psi function is simply the logit transformation which will be useful for the phi function, inputs: data, bandwidth, split value
psi = function(x,c,tau){
  input = (x-c)/tau # scaling of the distance
  psi = 1/(1+exp(-input)) # the logit transformation
  return(psi)
}

# The phi function is the likelihood of following a deterministic path down the tree, inputs: data, bandwidth, split value
phi = function(x,anc,tau){
  # Firstly, the probability of each consecutive branching step is calculated using the likelihood function
  prob = psi(x[anc[,'var']], anc[,'split_value'],tau)^(anc[,'left']) * (1-psi(x[anc[,'var']], anc[,'split_value'],tau))^(1-anc[,'left'])
  new.anc = cbind(terminal = anc[,'terminal'], prob = prob) # a new object with the probabilities at each branching step is constructed
  agg = aggregate(new.anc,
                  by = list(new.anc[,'terminal']),
                  FUN = prod)
  phi = agg[,'prob'] # The probabilities for each terminal node are multiplied with each other and the result is obtained
  return(phi)
}

# The design matrix function takes the element wise products of the relevant covariates in the terminal node and the phi-matrix
# The design matrix function constructs the design matrix based on the phi matrix and the variables involved in each terminal, input: tree, data, phi matrix and information about the internals
design_matrix = function(mtree,x, phi_matrix, int){

  design = matrix(1,nrow(x),1)*phi_matrix[,1]

  if(is.matrix(int)){

    for(j in 2:nrow(int)){
      new_design = x[,variable_list(mtree,int[j,'internal'])]*phi_matrix[,j]
      design = cbind(design,new_design)
    }

  }
  return(design)
}

# The condition tilde function computes the marginalized log likelihood for all nodes for a given tree in accordance to soft MOTR
# This function is inspired by the tree full conditional function and takes in the same input plus the number of trees
conditional_tilde = function(tree, X, R, sigma2, V, inv_V, nu, lambda, tau_b, ancestors, ntrees) {
  #Note that we need to add the number of trees to our input
  #Note that not all input variables are required for the function, but to keep matters congruent we leave them

  # T = ntrees
  N = nrow(X)
  p = ncol(X)

  X_node = X
  r_node = R

  # Calculation of covariance matrix
  Sigma = sigma2*diag(N) + (1/(ntrees*tau_b)) * X_node%*%t(X_node)
  # Sigma_inv = solve(sigma2*diag(N) + (1/(T*tau_b)) * X_node%*%t(X_node))
  temp_chol <- chol(Sigma)
  Sigma_inv = chol2inv(temp_chol)

  logdettemp <- log(prod(diag(temp_chol)^2))

  # log_lik= (-N/2)*log(2*pi) + (-1/2)*log(det(Sigma)) + -(1/2)*t(R)%*%Sigma_inv%*%R

  log_lik= (-N/2)*log(2*pi) + (-1/2)*logdettemp + -(1/2)*t(R)%*%Sigma_inv%*%R

  # log_lik= (-N/2)*log(2*pi) + (-1/2)*log(det(Sigma)) + -(1/2)*t(R)%*%Sigma_inv%*%R

  if(is.infinite(log_lik)){
    log_lik = -1e301
  }
  return(log_lik)
}

#The get w function determine the number of second generation internal nodes in a given tree, input: tree
get_w = function(tree){
  if(is.null(tree) | nrow(tree$tree_matrix) == 1) { # In case the tree is empty or a stump
    B = 0
  } else {
    indeces = which(tree$tree_matrix[,'terminal'] == '1') #determine which indexes have a terminal node label
    # determine the parent for each terminal node and sum the number of duplicated parents
    B = as.numeric(sum(duplicated(tree$tree_matrix[indeces,'parent'])))}
  return(B)
}

# The get number of terminals function calculates the number of terminal nodes, input: tree
get_nterminal = function(tree){
  if(is.null(tree)) { # In case the tree is empty
    L = 0
  } else {
    indeces = which(tree$tree_matrix[,'terminal'] == '1')  # determine which indexes have a terminal node label
    L = as.numeric(length(indeces))} # take the length of these indexes, to determine the number of terminal nodes/leaves
  return(L)
}

# The ratio grow function calculates ratio of transition probabilities of the grow step
# This function is in accordance with the soft BART paper, input: current tree, proposed/new tree
ratio_grow = function(curr_tree,new_tree){
  grow_ratio = get_nterminal(curr_tree)/(get_w(new_tree)+1)
  return(grow_ratio)
}

# The ratio prune function calculates ratio of transition probabilities of the prune step
# This function is in accordance with the soft BART paper, input: current tree, proposed/new tree
ratio_prune = function(curr_tree,new_tree){
  prune_ratio = get_w(new_tree)/(get_nterminal(curr_tree)-1)
  return(prune_ratio)
}

# The alpha MH function calculates the alpha of the Metropolis-Hastings step based on the type of transformation
# input: transformation type, current tree, proposed/new tree and the respective likelihoods
alpha_mh = function(l_new,l_old, curr_trees,new_trees, type){
  if(type == 'grow'){
    a = exp(l_new - l_old)*ratio_grow(curr_trees, new_trees)

  } else if(type == 'prune'){

    a = exp(l_new - l_old)*ratio_prune(curr_trees, new_trees)
  } else{
    a = exp(l_new - l_old)
  }
  return(a)
}

# The simulate beta tilde simulates the beta tilde in accordance to soft MOTR
# This function is inspired by the simulates beta function and takes in the same inputs
simulate_beta_tilde = function(tree, X, R, sigma2, inv_V, tau_b, nu, ancestors) {
  #Note that not all input variables are required for the function, but to keep matters congruent we leave them

  p = ncol(X)
  inv_V = diag(p)*inv_V
  X_node = X
  r_node = R

  # Lambda_node = chol2inv(chol(t(X_node)%*%X_node + inv_V))
  #
  # # Generate betas
  # beta_hat = rmvnorm(1,
  #                    mean = Lambda_node%*%(t(X_node)%*%r_node),
  #                    sigma = sigma2*Lambda_node)

  U = chol ( crossprod ( X_node )+ inv_V )
  IR = backsolve (U , diag ( p ))
  # btilde = crossprod ( t ( IR ))%*%( crossprod (X_node , r_node ) )
  # btilde = crossprod ( t ( IR ))%*%( crossprod (X_node , r_node ) )

  btilde = tcrossprod (  IR , crossprod( crossprod (X_node , r_node ), IR ))

  beta_hat = btilde + sqrt ( sigma2 )* IR %*% rnorm ( p )


  # Add the beta hat results to the tree matrix
  tree$tree_matrix[,'beta_hat'] = paste(beta_hat, collapse = ',')


  return(tree)
}

# The get beta hat function extracts the estimated beta hat matrix from the tree object in vector form, input: tree
get_beta_hat = function(sim_beta){
  beta_hat = as.numeric(unlist(strsplit(sim_beta$tree_matrix[[1,'beta_hat']], split=",")))
  return(beta_hat)
}

# The tau prior function calculates the prior probability for the bandwidth
# This is relevant for the second Metropolis Hastings step, input; tau value and tau rate
tau_prior= function(tau_value, tau_rate){
  tau = dexp(tau_value,tau_rate)
  return(tau)
}

# The paste betas function concatenates the estimated beta hats as a single vector, inputs: trees, number of trees
paste_betas = function(trees,ntrees){
  betas_trees = NULL # initialize vector

  for(i in 1:ntrees){ # loop though all the trees
    tree = trees[[i]]
    beta_hat = get_beta_hat(tree) # get the corresponding estimated beta hat for a single tree
    betas_trees = rbind(betas_trees,
                        cbind(betas_trees = beta_hat)) # concatenate the betas vectors vertically
  }
  return(betas_trees)
}

# The simulate tau b function simulates the full conditional that follows from a prior in Prado et al. (2021),
# inputs: beta hat vector, sigma squared and the prior parameters
simulate_tau_b = function(betas_trees,sigma2, a,b){

  # simulate tau_b for the gamma distribution, which is equivalent to simulating sigma beta for the inverse gamma
  tau_b = rgamma(1, shape = length(betas_trees)/2 + a, rate = t(betas_trees)%*%betas_trees/(2*sigma2) + b)

  return(tau_b)
}

# The test function produces the out of sample predictions with the hierarchical shrinkage method, input: test data and the model object
#' @export
test_function = function(newdata,object){

  # Get the means and standard deviation to standardize the covariates from the test data
  center = object$center_x
  scale = object$scale_x
  newdata = as.matrix(cbind(1,scale(newdata, center=center, scale=scale))) # standardize based on the centers and scales of the training data
  ntrees = object$ntrees

  # Create holder for predicted values
  n_its = object$npost

  # Initiliaze matrices to store the predictions for each observation and iteration
  preds = matrix(NA, nrow = nrow(newdata),
                 ncol = n_its)

  # Now loop through iterations and get predictions
  for(i in 1:n_its) {
    pred = numeric(nrow(newdata))

    for(j in 1:ntrees){

      # get the tree, beta vector and bandwidth of the soft motr object
      tree = object$trees[[i]][[j]]
      beta = object$beta_trees[[i]][[j]]
      int = get_internals(tree)

      # tree$node_indices = fill_tree_details(tree,newdata)$node_indices

      # get the branching information and bandwidth of the trained trees and apply to the test data
      if(!is.null(int)){

        # phi_matrix = phi_matrix(tree,int,newdata)
        # phi_matrix = phi_matrix_test(tree,int,newdata)

        int_temp <- int
        if(!is.matrix(int_temp)){
          # print("int_temp = ")
          # print(int_temp)

          int_temp <- t(int_temp)
        }


        phi_matrix = suppressWarnings(phi_app_hs_test(matrix(as.numeric(tree$tree_matrix),
                                                        nrow = nrow(tree$tree_matrix),
                                                        ncol = ncol(tree$tree_matrix)) ,
                                                 # matrix(as.numeric(tree$node_indices),
                                                 #        nrow = length(tree$node_indices),
                                                 #        ncol = 1)  ,
                                                 matrix(as.numeric(int_temp),
                                                        nrow = nrow(int_temp),
                                                        ncol = ncol(int_temp)) ,
                                                 as.matrix(newdata))
        )


        # phi_matrix = suppressWarnings(phi_app_hs(matrix(as.numeric(tree$tree_matrix),
        #                                                      nrow = nrow(tree$tree_matrix),
        #                                                      ncol = ncol(tree$tree_matrix)) ,
        #                                          tree$node_indices,
        #                                               # matrix(as.numeric(tree$node_indices),
        #                                               #        nrow = length(tree$node_indices),
        #                                               #        ncol = 1)  ,
        #                                               matrix(as.numeric(int_temp),
        #                                                      nrow = nrow(int_temp),
        #                                                      ncol = ncol(int_temp)) ,
        #                                               as.matrix(newdata))
        # )



        design = design_matrix(tree,newdata,phi_matrix,int)

        # calculate the model fit
        pred = pred + (design %*% beta)
      }
      else{
        pred = pred + rep(beta, nrow(newdata))
      }

    }

    # re-scale the predictions
    preds[,i] = object$y_mean + object$y_sd*pred

  }

  return(list(predictions = preds))

}
