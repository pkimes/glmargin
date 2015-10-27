// C++ header file for interval estimation (Kpegasos)
//
// Patrick Kimes


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>

#include <armadillo>
#include <random>
#include <cmath>
#include <stdio.h>
#include <iostream>

using namespace Rcpp;
using namespace arma;


class KPEGASOS {
public:
    //public use Rcpp classes

    KPEGASOS(mat r_x, vec r_y, 
             vec r_py, double r_lam,
             int r_max_iter, double r_min_eps,
             bool r_verbose);

    void Solve();
    double get_b();
    vec get_W();
    vec boundaries();

 private:
    //const static int MAX_ITER = 100;
    //constexpr static double MIN_EPS = 0.0000001;
    
    int i, j, k;
    
    int n, p, K;

    mat x;
    vec y;
    vec py;
    double lam;
    int max_iter;
    double min_eps;
    bool verbose;

    double b;
    vec W;

    vec h_b;
    mat h_W;

    mat pie_mat;
    mat hinge_mat;
    mat slope_mat;

    double Update();

    double eps;
    int iter = 0;
    
    double LogitAplus(double pyi);
    double LogitBplus(double pyi);
    double LogitCuts(double pyi);
    vec CalcHinges(vec py);

};

