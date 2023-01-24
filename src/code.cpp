#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace std;

// Creating the pos function
double pos(double x, double x_i){

     double dif = (x-x_i) ;

     // Getting the positive part only
     if( dif> 0 ){
          return dif;
     } else {
          return 0.0;
     }
}

// Creating the pos function
arma::vec pos_vec(arma::vec x, double x_i){

     arma::vec dif(x.size());

     for(int i = 0; i < x.size(); i++){

          // Getting the positive part only
          if( (x(i)-x_i)>0){
               dif(i) = x(i)-x_i;
          } else {
               dif(i) = 0.0;
          }
     }

     return dif;

}

// Function to generate B
arma::mat bspline(arma::vec x,
                  arma::vec x_obs){

     arma::mat B(x.n_rows, x_obs.n_rows+2, arma::fill::ones);

     // Storing a copy
     arma::vec x_obs_copy = x_obs;
     arma::uword max_ind = x_obs.index_max();
     x_obs_copy(max_ind) = -std::numeric_limits<double>::infinity();
     double x_n_1 = x_obs_copy.max();

     // Setting the values for all columns
     B.col(1) = x;
     double x_n = max(x_obs);

     for(int i = 0; i < x_obs.n_rows; i++) {
          if((x_n-x_obs(i))!=0){
               B.col(i+2) = (pow(pos_vec(x,x_obs(i)),3) - pow(pos_vec(x,x_n),3))/(x_n-x_obs(i)) - (pow(pos_vec(x,x_n_1),3)-pow(pos_vec(x,x_n),3))/(x_n-x_n_1);
          } else {
               B.col(i+2) = arma::zeros(x.n_rows);
          }
     }
     return B;
}

// Sample beta coefficients
arma::vec beta_sample(arma::mat bspline, arma::vec y,
                      double tau, double tau_b){

     arma::mat bspline_t = bspline.t();
     arma::mat btb = bspline_t*bspline;

     arma::mat precision_diag  = arma::eye<arma::mat>(bspline.n_cols,bspline.n_cols)*(tau_b/tau);

     // Left part of precision, it also will be used in the mean and in the variance
     arma::mat inv_b_precision = inv(btb+precision_diag);

     arma::vec mvn_mean = inv_b_precision*(bspline_t*y);

     // Getting the random sample
     arma::mat sample = arma::randn<arma::mat>(bspline.n_cols,1);
     arma::mat result = arma::chol(inv_b_precision,"lower")*sample + mvn_mean;

     return result;
}

// Updating the tau parameter
double updateTau(arma::vec &y_hat,
               arma::vec& y,
               double a_tau,
               double d_tau){

     // Getting the sum of residuals square
     double tau_res_sq_sum = dot((y_hat-y),(y_hat-y));

     double tau = R::rgamma((0.5*y.size()+a_tau),1/(0.5*tau_res_sq_sum+d_tau));

     return tau;
}

//[[Rcpp::export]]
Rcpp::List mcmc_sampler(arma::vec x,
                        arma::vec x_new,
                        arma::vec y,
                        int n_post,
                        int n_burn,
                        double tau = 1.0){


          // Creating a progress bar
          const int width = 70;
          double pb = 0;


          // Getting the spline matrix
          arma::mat B = bspline(x, x);
          arma::mat B_new = bspline(x_new,x);

          // Getting MCMC samples
          int n_mcmc = n_post + n_burn;
          arma::mat y_hat_post(n_post,y.size(),arma::fill::zeros);
          arma::mat y_hat_new_post(n_post,x_new.size(),arma::fill::zeros);
          arma::mat beta_post(n_post,B.n_cols,arma::fill::zeros);
          arma::vec tau_post(n_post,arma::fill::zeros);
          int curr = 0;

          // Initialise the sampled quantities
          arma::vec beta(B.n_cols);
          arma::vec y_hat(y.size());
          arma::vec y_hat_new(x_new.size());

          // Initialising the MCMC sampler
          for(int  i=0; i < n_mcmc; i++) {

               // Initialising PB
               std::cout << "[";
               int k = 0;
               // Evaluating progress bar
               for(;k<=pb*width/n_mcmc;k++){
                        std::cout << "=";
               }

               for(; k < width;k++){
                        std:: cout << " ";
               }

               std::cout << "] " << std::setprecision(5) << (pb/n_mcmc)*100 << "%\r";
               std::cout.flush();

               // cout << "Error on beta sample" << endl;
               beta = beta_sample(B,y,tau,100);
               // cout << "Error on y_hat" << endl;
               y_hat = B*beta;
               y_hat_new = B_new*beta;
               // cout << "Error on updateTau" << endl;
               tau = updateTau(y_hat,y,0.00001,0.00001);

               // Storing only the post
               if(i > n_burn){
                    y_hat_post.row(curr) = y_hat.t();
                    y_hat_new_post.row(curr) = y_hat_new.t();
                    beta_post.row(curr) = beta.t();
                    tau_post(curr) = tau;
                    curr ++ ;

               }

               // Adding one MCMC iteration
               pb++;
          }


          // FINISHING THE PROGRESS BAR
          std::cout << "[";
          int k = 0;
          // Evaluating progress bar
          for(;k<=pb*width/n_mcmc;k++){
                  std::cout << "=";
          }

          for(; k < width;k++){
                  std:: cout << " ";
          }

          std::cout << "] " << std::setprecision(5) << 100 << "%\r";
          std::cout.flush();

          std::cout << std::endl;

          return Rcpp::List::create(y_hat_post,
                                    y_hat_new_post,
                                    beta_post,
                                    tau_post);

}



