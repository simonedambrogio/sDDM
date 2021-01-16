#include <Rcpp.h>
#include <RcppParallel.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <time.h>

using namespace Rcpp;
using namespace RcppParallel;

typedef boost::mt19937                     ENG;    // Mersenne Twister
typedef boost::normal_distribution<double> DIST;   // Normal Distribution
typedef boost::variate_generator<ENG,DIST> GEN;    // Variate generator

// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppParallel)]]
struct sddm : public Worker {
      
  // input vectors/matrices to read from
  const double dp; //delta_p (personal drift)
  const double B_intercept; //Beta (starting point)
  const double B_slope; //Beta (starting point)
  const double temptation; //Beta (starting point)
  const double theta; //Boundary separation
  const double nDT; // non decision time
  const double social_strength; //M(t) = N_a(t) - N_b(t)
  const double dt; //Delta_t (time step)
  const double qu; //shapes the power function (social drift)
  const double s; //scaling parameter that influences the strength of social drift
  GEN gen;
  
  // output vector to write to
  RVector<double> vecOut;
  
  // initialize from Rcpp input and output matrices/vectors (the RMatrix/RVector class
  // can be automatically converted to from the Rcpp matrix/vector type)
  sddm(const double dp, const double B_intercept, const double B_slope, const double theta, const double nDT, const double social_strength, 
       const double temptation, const double dt, const double qu, const double s, NumericVector vecOut , GEN gen)
    : dp(dp), B_intercept(B_intercept), B_slope(B_slope), theta(theta), nDT(nDT), social_strength(social_strength), 
      temptation(temptation), dt(dt), qu(qu), s(s), gen(gen), vecOut(vecOut) {}
  
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    
    double sign; //Sign of the strength of social information 
    double ds; //Social Drift
    double B; //Starting Point
    
    // Social Drift
    if(social_strength==0){sign=0;} else {sign=social_strength/abs(social_strength);}
    ds = s * sign * pow(abs(social_strength), qu);
    //Starting Point
    B = 1/(1+exp(B_intercept*(temptation-B_slope)));
    B = B*0.9+0.05; //Avoid start too close to boundary
      
    for (std::size_t i = begin; i < end; i++) {
      vecOut[i] = 6;
      int flag = 0;
      double noise;
      double iterations = 0; //Iterations
      double idx=0;
      
      //Starting point
      double L = theta*B - (theta - theta*B);
      
      while (flag==0 && iterations < (6/dt)) {
        
        noise=gen()*sqrt(dt);
        L = L + (dp + ds)*dt + noise;
        
        if (L > theta) { flag=1; vecOut[i] =  nDT + iterations*dt; };
        if (L < -theta) { flag=1; vecOut[i] =  -nDT -iterations*dt; };
        
        iterations++;
      }
    }
  }
};


// [[Rcpp::export]]
NumericVector sddm_parallel(double dp, double B_intercept, double B_slope, double theta, double nDT, double social_strength, double temptation, double dt, double qu, double s, unsigned int N) {
  
  const double sd_n = 1;
  struct timespec time;
  clock_gettime(CLOCK_REALTIME, &time);
  ENG  eng;
  eng.seed(time.tv_nsec);
  DIST dist(0,sd_n);
  GEN  gen(eng,dist);
  
  //output vector
  NumericVector vecOut(N);
  
  // create the worker
  sddm sddm(dp, B_intercept, B_slope, theta, nDT, social_strength, temptation, dt, qu, s, vecOut, gen);
  
  // call the worker
  parallelFor(0, N, sddm);
  
  return vecOut;
}
