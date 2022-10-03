
# include <RcppArmadillo.h>
# include <cmath>

using namespace Rcpp;


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat phi_app_hs( arma::mat treemat,
                      arma::vec node_indices,
                      arma::mat internalmat,
                      arma::mat xmat) {


  // int treenrow = treemat.n_rows;


  // arma::field<std::string> treemat(treemat1.nrow(), treemat1.ncol());
  //
  // arma::field<std::string> internalmat();

  // Rcpp::Rcout << "treemat = " << treemat << " .\n " ;
  // Rcpp::Rcout << "internalmat = " << internalmat << " .\n " ;

  //ensure that int is a matrix (maybe use as.matrix)

  int ntemp = xmat.n_rows;

  // Rcpp::Rcout << "internalmat.n_rows = " << internalmat.n_rows << " .\n " ;

  // Rcpp::Rcout << "ntemp = " << ntemp << " .\n " ;

  arma::mat phi_matrix(ntemp, internalmat.n_rows);



  for(unsigned int j=0; j < internalmat.n_rows ; j++){

    // Rcpp::Rcout << "Line 39, j = " <<  j << " .\n " ;


    double child_left = internalmat(j,2);
    double child_right = internalmat(j,3);
    double n_left = internalmat(j, 4);
    double n_right = internalmat(j, 5);

    // Rcpp::Rcout << "Line 47, j = " <<  j << " .\n " ;

    double tempprod = n_left*n_right;
    double temp_denom = std::sqrt(tempprod );

    for(int i=0; i < ntemp ; i++){

      // Rcpp::Rcout << "Line 54, i = " <<  i << " .\n " ;


      int term_ind = node_indices(i) ;
      int parent =  treemat(term_ind-1 , 3) ;

      // Rcpp::Rcout << "Line 60, i = " <<  i << " .\n " ;

      arma::uvec parents = { term_ind, parent};
      // arma::vec parents = { term_ind, parent};

      double parent_d = parent;

      // Rcpp::Rcout << "Line 67, i = " <<  i << " .\n " ;


      while(!( std::isnan(parent))){
        int parent_min1 = parent  - 1;
        // int parentint = parent_d;

        // Rcpp::Rcout << "Line 74, i = " <<  i << " .\n " ;
        // Rcpp::Rcout << "Line 74, parent_min1 = " <<  parent_min1 << " .\n " ;
        // Rcpp::Rcout << "Line 74, parent = " <<  parent << " .\n " ;


        if(std::isnan(treemat(parent_min1 , 3))){
          arma::uvec parenttemp = {  arma::datum::nan};
          // Rcpp::Rcout << "Line 81, i = " <<  i << " .\n " ;


          parents = arma::join_cols(parents, parenttemp );

          // Rcpp::Rcout << "Line 87, i = " <<  i << " .\n " ;
          double numeqleft = arma::sum( parents == child_left) ;
          double numeqright = arma::sum( parents == child_right) ;


          // Rcpp::Rcout << "Line 92, i = " <<  i << " .\n " ;

          phi_matrix(i,j) = (n_right*numeqleft - n_left*numeqright )/temp_denom;
          break;

        }else{
          parent =  treemat(parent_min1 , 3)  ;
          arma::uvec parenttemp = {  parent};
          // Rcpp::Rcout << "Line 100, i = " <<  i << " .\n " ;


          parents = arma::join_cols(parents, parenttemp );

          // Rcpp::Rcout << "Line 105, i = " <<  i << " .\n " ;
          double numeqleft = arma::sum( parents == child_left) ;
          double numeqright = arma::sum( parents == child_right) ;


          // Rcpp::Rcout << "Line 110, i = " <<  i << " .\n " ;

          phi_matrix(i,j) = (n_right*numeqleft - n_left*numeqright )/temp_denom;

        }
        // double parent_d = parent;




        // arma::field<std::string> newparents(parents.n_elem + 1);
        // newparents.subfield(0,parents.n_elem-1, 0 , 0) =parents ;
        // newparents(parents.n_elem, 0 ) = parent ;

        // newparents = {parents, parenttemp };
        // parents = newparents;


        // arma::uvec templeftvec = (arma::find(parents == child_left));
        // arma::uvec temprightvec = (arma::find(parents == child_right));

        // double numeqleft = templeftvec.n_elem ;
        // double numeqright = temprightvec.n_elem ;

        // double numeqleft = 0;
        // double numeqright = 0 ;
        // for(unsigned int k=0; k < parents.n_elem; k++){
        //   if(parents(k) == child_left){
        //     numeqleft = numeqleft +1;
        //   }
        //   if(parents(k) == child_right){
        //     numeqright = numeqleft +1;
        //   }
        // }



      }


    }

  }


  return( phi_matrix);
}










// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat phi_app_hs_test( arma::mat treemat,// arma::vec node_indices1,
                      arma::mat internalmat,
                      arma::mat xmat) {

  // arma::vec node_indices = node_indices1;

  arma::vec node_indices(xmat.n_rows, arma::fill::value(1));

  // For all but the top row, find the number of observations falling into each one
  // Therefore loop intentionally begins at 1 instead of zero
  for( int i=1; i < treemat.n_rows ; i++){
    int curr_parent =  treemat(i , 3) ;
    int split_var = treemat(curr_parent - 1 , 4) ;
    double split_val = treemat(curr_parent - 1 , 5) ;

    arma::vec splitcol = xmat.col(split_var - 1);


    // Rcpp::Rcout << "splitcol = " << splitcol << " .\n " ;
    // Rcpp::Rcout << "split_var = " << split_var << " .\n " ;
    // Rcpp::Rcout << "split_val = " << split_val << " .\n " ;

    // subtracting 1 because indexing i from zero, whereas matrix elements defined for indexing from 1
    bool left_or_right = ( (treemat(curr_parent - 1 , 1) -1)  == i);

    double curr_parent_d = curr_parent;

    // Rcpp::Rcout << "i = " << i << " .\n " ;

    // Rcpp::Rcout << "curr_parent_d = " << curr_parent_d << " .\n " ;
    // Rcpp::Rcout << "treemat(curr_parent - 1 , 1) = " << treemat(curr_parent - 1 , 1) << " .\n " ;
    // Rcpp::Rcout << "left_or_right = " << left_or_right << " .\n " ;

    if(left_or_right == true){
      // arma::uvec found_inds = (node_indices == curr_parent) ;
      // arma::uvec found_left = (splitcol < split_val) ;
      // Rcpp::Rcout << "(node_indices == curr_parent_d) = " << (node_indices == curr_parent_d) << " .\n " ;
      // Rcpp::Rcout << " (splitcol >= split_val) = " <<  (splitcol < split_val) << " .\n " ;


      arma::uvec found_inds = arma::find( (node_indices == curr_parent_d)&& (splitcol < split_val) ) ;
      // Rcpp::Rcout << "found_inds = " << found_inds << " .\n " ;

      // arma::uvec overall_found = arma::find(found_inds*found_left );
      // not sure whether to set to i or i +1
      // double temp_d = i  + 1;
      node_indices.elem(found_inds) = (i+1)*arma::ones(found_inds.n_elem);
      // for(unsigned int k=0; k < found_inds.n_elem ; i++){
      //   node_indices(found_inds(k)) = i+1 ;
      //
      // }

    }else{
      // arma::uvec found_inds = (node_indices == curr_parent) ;
      // arma::uvec found_right = (splitcol >= split_val) ;
      arma::uvec found_inds = arma::find((node_indices == curr_parent_d)&& (splitcol >= split_val) ) ;

      // arma::uvec overall_found = arma::find(found_inds*found_right );
      // not sure whether to set to i or i +1
      // node_indices.elem(found_inds) = i +1 ;

      // node_indices.elem(arma::find((node_indices == curr_parent)&& (splitcol >= split_val))) = i +1 ;
      node_indices.elem(found_inds) = (i+1)*arma::ones(found_inds.n_elem);

      // for(unsigned int k=0; k < found_inds.n_elem ; i++){
      //   node_indices(found_inds(k)) = i+1 ;
      //
      // }
    }
    // Rcpp::Rcout << "End iteration, i = " << i << " .\n " ;
    // Rcpp::Rcout << "node_indices = " << node_indices << " .\n " ;
  }

  // arma::uvec found_inds = arma::find(node_indices == 1 ) ;
  // node_indices.elem(found_inds) = 2*arma::ones(found_inds.n_elem);

  // node_indices = node_indices +1;
  // Rcpp::Rcout << "node_indices = " << node_indices << " .\n " ;



  // int treenrow = treemat.n_rows;


  // arma::field<std::string> treemat(treemat1.nrow(), treemat1.ncol());
  //
  // arma::field<std::string> internalmat();

  // Rcpp::Rcout << "treemat = " << treemat << " .\n " ;
  // Rcpp::Rcout << "internalmat = " << internalmat << " .\n " ;

  //ensure that int is a matrix (maybe use as.matrix)

  int ntemp = xmat.n_rows;

  // Rcpp::Rcout << "internalmat.n_rows = " << internalmat.n_rows << " .\n " ;

  // Rcpp::Rcout << "ntemp = " << ntemp << " .\n " ;

  arma::mat phi_matrix(ntemp, internalmat.n_rows);



  for(unsigned int j=0; j < internalmat.n_rows ; j++){

    // Rcpp::Rcout << "Line 256, j = " <<  j << " .\n " ;


    double child_left = internalmat(j,2);
    double child_right = internalmat(j,3);
    double n_left = internalmat(j, 4);
    double n_right = internalmat(j, 5);

    // Rcpp::Rcout << "Line 264, j = " <<  j << " .\n " ;

    double tempprod = n_left*n_right;
    double temp_denom = std::sqrt(tempprod );

    for(int i=0; i < ntemp ; i++){

      // Rcpp::Rcout << "Line 271, i = " <<  i << " .\n " ;
      // Rcpp::Rcout << "Line 272, node_indices = " <<  node_indices << " .\n " ;

      int term_ind = node_indices(i) ;
      double parent =  treemat(term_ind-1 , 3) ;
      // Rcpp::Rcout << "Line 272, term_ind-1 = " <<  term_ind-1 << " .\n " ;
      // Rcpp::Rcout << "Line 272, treemat = " <<  treemat << " .\n " ;
      // Rcpp::Rcout << "Line 272, parent = " <<  parent << " .\n " ;

      // Rcpp::Rcout << "Line 277, i = " <<  i << " .\n " ;

      arma::uvec parents = { term_ind, parent};
      // arma::vec parents = { term_ind, parent};

      double parent_d = parent;

      // Rcpp::Rcout << "Line 284, i = " <<  i << " .\n " ;


      while(!( std::isnan(parent))){
        int parent_min1 = parent  - 1;
        // int parentint = parent_d;

        // Rcpp::Rcout << "Line 291, i = " <<  i << " .\n " ;
        // Rcpp::Rcout << "Line 74, parent_min1 = " <<  parent_min1 << " .\n " ;
        // Rcpp::Rcout << "Line 293, parent = " <<  parent << " .\n " ;


        if(std::isnan(treemat(parent_min1 , 3))){
          arma::uvec parenttemp = {  arma::datum::nan};

          // Rcpp::Rcout << "Line 299, i = " <<  i << " .\n " ;

          parents = arma::join_cols(parents, parenttemp );

          // Rcpp::Rcout << "Line 87, i = " <<  i << " .\n " ;
          double numeqleft = arma::sum( parents == child_left) ;
          double numeqright = arma::sum( parents == child_right) ;


          // Rcpp::Rcout << "Line 92, i = " <<  i << " .\n " ;

          phi_matrix(i,j) = (n_right*numeqleft - n_left*numeqright )/temp_denom;
          break;

        }else{
          parent =  treemat(parent_min1 , 3)  ;
          arma::uvec parenttemp = {  parent};
          // Rcpp::Rcout << "Line 316, i = " <<  i << " .\n " ;


          parents = arma::join_cols(parents, parenttemp );

          // Rcpp::Rcout << "Line 105, i = " <<  i << " .\n " ;
          double numeqleft = arma::sum( parents == child_left) ;
          double numeqright = arma::sum( parents == child_right) ;


          // Rcpp::Rcout << "Line 110, i = " <<  i << " .\n " ;

          phi_matrix(i,j) = (n_right*numeqleft - n_left*numeqright )/temp_denom;

        }
        // double parent_d = parent;




        // arma::field<std::string> newparents(parents.n_elem + 1);
        // newparents.subfield(0,parents.n_elem-1, 0 , 0) =parents ;
        // newparents(parents.n_elem, 0 ) = parent ;

        // newparents = {parents, parenttemp };
        // parents = newparents;


        // arma::uvec templeftvec = (arma::find(parents == child_left));
        // arma::uvec temprightvec = (arma::find(parents == child_right));

        // double numeqleft = templeftvec.n_elem ;
        // double numeqright = temprightvec.n_elem ;

        // double numeqleft = 0;
        // double numeqright = 0 ;
        // for(unsigned int k=0; k < parents.n_elem; k++){
        //   if(parents(k) == child_left){
        //     numeqleft = numeqleft +1;
        //   }
        //   if(parents(k) == child_right){
        //     numeqright = numeqleft +1;
        //   }
        // }



      }


    }

  }


  return( phi_matrix);
}

