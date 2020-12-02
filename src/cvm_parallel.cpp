#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

struct Mycvm: public Worker {
  const RMatrix<double> mat1,mat2;

  RMatrix<double> rmat;

  Mycvm(const NumericMatrix mat1,NumericMatrix mat2, NumericMatrix rmat)
    : mat1(mat1), mat2(mat2), rmat(rmat){}

  void operator()(std::size_t begin, std::size_t end) {

    for (std::size_t i = begin; i < end; i++) {

      for (std::size_t j = 0; j < mat2.nrow(); j++) {

        double p = 1.0;
        for (int k = 0; k < mat1.ncol(); k++) {

          double x1, x2;
          x1 = 1 - mat1(i,k);
          x2 = 1 - mat2(j,k);

          if(x1 < x2) p *= x1;
          else p *= x2;
        }
        rmat(i,j) = p;
      }
    }
  }
};

NumericMatrix rcpp_cvm_term(NumericMatrix mat1,NumericMatrix mat2) {

  NumericMatrix rmat(mat1.nrow(), mat2.nrow());

  NumericVector output;

  Mycvm mycvm(mat1,mat2,rmat);

  parallelFor(0, mat1.nrow(), mycvm);

  return rmat;
}

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
using namespace RcppParallel;

struct Sum : public Worker
{
  const RMatrix<double> input_data;

  double value;

  Sum(const NumericMatrix input) : input_data(input), value(0) {}

  Sum(const Sum& sum, Split) : input_data(sum.input_data), value(0) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++){
      for(std::size_t j = 0; j < input_data.ncol(); j++){
        value += input_data(i,j);
      }
    }
  }

  void join(const Sum& rhs) {
    value += rhs.value;
  }
};

double parallelMatrixSum(NumericMatrix x) {
  Sum sum(x);

  parallelReduce(0, x.nrow(), sum);

  return sum.value;
}

// [[Rcpp::depends(RcppParallel)]]
NumericVector rcpp_pobs_vec(NumericVector x){

  NumericVector y = clone(x);
  NumericVector r;

  std::sort(y.begin(), y.end());

  r = match(x,y);
  return r;
}

// [[Rcpp::export]]
NumericMatrix rcpp_pobs(NumericMatrix x){

  NumericMatrix y(x.nrow(),x.ncol());

  for(int i = 0; i < y.ncol(); i++){
    y(_,i) = rcpp_pobs_vec(x(_,i)) / (x.nrow()+1);
  }

    return y;
}

// [[Rcpp::export]]
double rcpp_cvm(NumericMatrix mat1,NumericMatrix mat2) {

  NumericMatrix rmat1,rmat2,rmat3;
  double term01,term02,term03;
  double term1,term2,term3,term4,stat;
  double n = mat1.nrow();
  double m = mat2.nrow();

  rmat1 = rcpp_cvm_term(mat1,mat1);
  rmat2 = rcpp_cvm_term(mat2,mat2);
  rmat3 = rcpp_cvm_term(mat1,mat2);

  term01 = parallelMatrixSum(rmat1);
  term02 = parallelMatrixSum(rmat2);
  term03 = parallelMatrixSum(rmat3);

  term1 = term01 / (n * n);
  term2 = term02 / (m * m);
  term3 = 2 * term03 / (n * m);
  term4 = n * m / (n + m);
  stat = term4 * (term1 + term2 - term3);
  return stat;
}


// [[Rcpp::export]]
List rcpp_coptest(NumericMatrix mat1, NumericMatrix mat2,int nperm){

  int N = mat1.nrow() + mat2.nrow();
  int D = mat1.ncol();
  NumericMatrix U(mat1.nrow(),mat1.ncol());
  NumericMatrix V(mat2.nrow(),mat2.ncol());
  double stat;

  U = rcpp_pobs(mat1);
  V = rcpp_pobs(mat2);

  stat = rcpp_cvm(U,V);

  NumericMatrix W(N,D);
  int cnt = 0;

  for(int i = 0; i < U.nrow(); i++){
    W(cnt,_) = U(i,_);
    cnt++;
  }

  for(int i = 0; i < V.nrow(); i++){
    W(cnt,_) = V(i,_);
    cnt++;
  }

  NumericVector idx;
  NumericMatrix UU(mat1.nrow(),mat1.ncol());
  NumericMatrix VV(mat2.nrow(),mat2.ncol());
  NumericVector null_stat(nperm);

  for(int i = 0; i < nperm; i++){
    idx = sample(N, N) - 1;

    int cnt = 0;
    for(int j = 0; j < UU.nrow(); j++){
      UU(j,_) = W(idx(cnt),_);
      cnt++;
    }

    for(int j = 0; j < VV.nrow(); j++){
      VV(j,_) = W(idx(cnt),_);
      cnt++;
    }

    UU = rcpp_pobs(UU);
    VV = rcpp_pobs(VV);

    null_stat(i) = rcpp_cvm(UU,VV);
  }

  int score = 0;

  for(int i = 0; i < null_stat.length(); i++){
    if(stat >= null_stat[i]){
      score++;
    }
  }

  double pval;
  pval = 1 - (double) score / (double) null_stat.length();

  return List::create(Named("stat")=stat,Named("pval")=pval,Named("null_stat")=null_stat);
}


