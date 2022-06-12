
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

          phi_matrix(i,j) = (n_right*numeqleft + n_left*numeqright )/temp_denom;
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

          phi_matrix(i,j) = (n_right*numeqleft + n_left*numeqright )/temp_denom;

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
