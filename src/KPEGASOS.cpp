// C++ implementation for interval estimation (Kpegasos)
//
// Patrick Kimes


//include corresponding header
#include "KPEGASOS.h"


// constructor method
KPEGASOS::KPEGASOS(mat r_x, vec r_y, 
                   vec r_py, double r_lam,
                   int r_max_iter, double r_min_eps,
                   bool r_verbose) {
    
    x = mat(r_x);
    y = vec(r_y);
    py = vec(r_py);
    py = sort(py);
    lam = r_lam;
    max_iter = r_max_iter;
    min_eps = r_min_eps;
    verbose = r_verbose;

    //dimensions of problem
    n = x.n_rows;
    p = x.n_cols;
    K = py.n_rows;

    //temporary R scope for RNG
    RNGScope scope;

    //set default values for parameters
    b = 0;
    W = Rcpp::floor(Rcpp::runif(p)+0.5)*2.0 - 1.0;

}


//getter functions
double KPEGASOS::get_b() {
    return b;
}

vec KPEGASOS::get_W() {
    return W;
}


//solve function
void KPEGASOS::Solve() {

    //copy of labels taking value in 0, 1
    mat y_01 = vec(y);
    y_01 = (y_01 + 1) / 2;


    //table of cols, py, 1-py
    mat pie_table(K, 2);
    pie_table.col(0) = py;
    pie_table.col(1) = 1.0 - pie_table.col(0);

    //K x n matrix of weights for each sample
    pie_mat = mat(K, n);
    for (i = 0; i < n; ++i) {
        pie_mat.col(i) = pie_table.col(y_01(i));
    }

    
    //table of hinges, +class, -class
    mat hinge_table(K, 2);
    hinge_table.col(0) = CalcHinges(py);
    hinge_table.col(1) = CalcHinges(1.0 - py);

    //K x n matrix of hinges for each sample
    hinge_mat = mat(K, n);
    for (i = 0; i < n; ++i) {
        hinge_mat.col(i) = hinge_table.col(y_01(i));
    }

    slope_mat = mat(K, 2);
    for (k = 0; k < K; ++k) {
        slope_mat(k, 0) = LogitBplus(py(k));
        slope_mat(K-1-k, 1) = LogitBplus(1.0 - py(k));
    }

    //update parameters iteratively
    do {
        iter++;
        eps = Update();
        //std::cout << "eps: " << eps << std::endl;
        
    } while (eps > min_eps && iter < max_iter);
    
    if (verbose) {
        std::cout << "converged in " << iter << " steps, " 
                  << "eps=" << eps << std::endl;
    }

}



//helper function for updating W,b
double KPEGASOS::Update() {

    vec W_1 = vec(W);
    double b_1 = b;

    //compute functional margin
    vec fm = y % (x * W + b);
    mat fm_mat = ones<vec>(K) * fm.t();
    umat is_pos = (hinge_mat <= fm_mat); //use hinge_mat to define hinge locations
    fm_mat = is_pos % ones<mat>(K, n);
    vec bi = fm_mat.t() * ones<mat>(K); //get sum to find total

    //update W
    W -= (1.0 / n) * W_1;
    W += (1.0 / (n*iter*lam)) * (x.t() * (bi % y));

    //update b
    b += (1.0 / (n*iter*lam)) * accu(bi % y);

    //project to 1/lam sphere
    double c = std::pow(lam, -0.5) / (accu(square(W)) + std::pow(b, 2));
    if (c > 1.0) { c = 1.0; }
    W = c*W;
    b = c*b;

    return accu(square(W - W_1)) + pow(b - b_1, 2);
}



//internals for determining shape of hinge loss
double KPEGASOS::LogitAplus(double pyi) {
    if (pyi >= 1.0 | pyi <= 0.0) {
        return 0;
    } else {
        return -pyi * log(pyi) - (1.0 - pyi) * log(1.0 - pyi);
    }
}

double KPEGASOS::LogitBplus(double pyi) {
    return pyi - 1.0;
}

double KPEGASOS::LogitCuts(double pyi) {
    return log(pyi / (1.0 - pyi));
}


//determine hinges
vec KPEGASOS::CalcHinges(vec py) {

    vec slopes(K);
    vec interc(K);
    vec hinges(K);

    py = sort(py);

    for (k = 0; k < K; ++k) {
        interc(k) = LogitAplus(py(k));
        slopes(k) = LogitBplus(py(k));
    }
    
    hinges(span(0, K-2)) = interc(span(1, K-1)) - interc(span(0, K-2));
    hinges(span(0, K-2)) %= - (slopes(span(1, K-1)) - slopes(span(0, K-2)));
    hinges(K-1) = - (interc(K-1) / slopes(K-1));

    return hinges;
}



//compute test error for new x, y
vec KPEGASOS::boundaries() {
    vec bounds = vec(K);
    for (k = 0; k < K; ++k) {
        bounds[k] = log(py[k] / (1.0 - py[k]));
    }
    return bounds;
}



//expose to R using Rcpp-modules
RCPP_MODULE(kpegasos) {
    class_<KPEGASOS>("KPEGASOS")
        .constructor<mat, vec, vec, double, int, double, bool>("x train, y train, py, lambda, max iter, min eps, verbose")
        .property("b", &KPEGASOS::get_b)
        .property("W", &KPEGASOS::get_W)
        .method("Solve", &KPEGASOS::Solve)
        .method("boundaries", &KPEGASOS::boundaries)
        ;
}
