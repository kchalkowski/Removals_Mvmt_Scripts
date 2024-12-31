#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
//define the main function, MovementRcpp
arma::vec TestPairwise(arma::vec& x1,
					   arma::vec& y1,
                       arma::vec& x2,
                       arma::vec& y2,
                       const double& cdist) {

//inputs:
//x1-vector of x coordinates for individual 1
//y1-vector of y coordinates for individual 1
//x2-vector of x coordinates for individual 2
//y2-vector of y coordinates for individual 2
//cdist is a numeric/double specifying distance defining contact
//all x/y vectors must be same length

//outputs distances in same units as projected x/y coordinates,
    //is vector of distances between the two individuals for each timestamp

double dist;
arma::vec distances(x1.size());
arma::uvec contact(x1.size());

for(std::size_t s=0; s < x1.size(); ++s) {

distances(s)=sqrt(pow((x1(s)-x2(s)),2)+pow((y1(s)-y2(s)),2));

} //loop through xy vectors closing bracket

return distances;                        
                       } //function closing bracket
