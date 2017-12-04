// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::plugins(cpp11)]]
#include <RcppEigen.h>
#include <Rmath.h>

#include <string>
#include <vector>
#include <algorithm> // for std::unique
#include <iterator>

using namespace Eigen;

typedef Eigen::Map<Eigen::VectorXd> MapVxd;
typedef Eigen::Map<Eigen::MatrixXd> MapMxd;
typedef Eigen::Map<Eigen::VectorXi> MapVxi;
typedef const Eigen::Ref<const Eigen::VectorXd>& VecRef;
typedef const Eigen::Ref<const Eigen::MatrixXd>& MatRef;
typedef const Eigen::Ref<const Eigen::ArrayXd> & ArrRef;
typedef const Eigen::Ref<const Eigen::ArrayXXd> & Arr2Ref;

/*
Implementing the methodology of Diao, Hanlon, and Vidyashankar 2014:  multiple testing for high-dimensional data


*/

template <typename T>
void printVec(std::vector<T> &v){
	for(T i : v){ Rcpp::Rcout << i << " ";}
	Rcpp::Rcout << std::endl;
}

template <typename T>
int partition(T *V,int *P, const int left, const int right){ //note: it is up to this function's caller to ensure that left and right are within V's range!
	const int mid = left+(right-left)/2; // avoid overflow
	const T pivot = V[mid];
	std::swap(V[left],V[mid]);
	std::swap(P[left],P[mid]);
	int i = left + 1;
	int j = right;
	while(i <= j){ // until the indices overlap
		while( i <= j && V[i] <= pivot){
			i++;
		}
		while(i <= j && V[j] > pivot){
			j--;
		}
		if(i < j){
			std::swap(V[i],V[j]);
			std::swap(P[i],P[j]);
		}
	}
	std::swap(V[i-1],V[left]); // return the pivot value to the point where it is in order
	std::swap(P[i-1],P[left]); // keeping indices consistent with values
	return i-1; // return the location of the pivot: the recursive step splits the array here
}

template<typename T>
void qisort(T *V, int *P, const int leftInx, const int rightInx, const int SIZE){ // SIZE is the length of the entire array, and is not really necessary
	if(leftInx >= rightInx){
		return;
	}
	// do the swapping and return the index defining the cutpoint
	int part = partition(V,P,leftInx,rightInx);
	qisort(V,P,leftInx,part-1,SIZE);
	qisort(V,P,part+1,rightInx,SIZE);
}

// take a given (dynamic size) Eigen vector of any supported type, sort it, and return a VectorXi map between sorted & original positions
template<typename T>
VectorXi EigenSort(Eigen::Ref<Eigen::Matrix<T,Eigen::Dynamic,1> > data){
	const int L = data.size();
	T* vv = &data(0);
	VectorXi inxs = VectorXi::LinSpaced(L,0,L-1); // initial, zero-based indices
	int* pp = &inxs(0);
	qisort(vv,pp,0,L-1,L);
	return inxs;
}

 // generate an i.i.d. sequence of Rademacher random variables
VectorXd rademachers(int N){
	ArrayXd res(N);
	for(int i=0;i<N;i++){
		res(i) = Rf_rbinom(1.,0.5);
	}
	return 2*res - 1;
}

// return a matrix of normal random variates; if P == 1 then its just a vector.
MatrixXd rnorm(unsigned int N,unsigned int P,const double mu,const double sig){
	MatrixXd res(N,P);
	for(unsigned int i=0;i<N;i++){
		for(unsigned int j=0;j<P;j++){
			res(i,j) = Rf_rnorm(mu,sig);
		}
	}
	return res;
}

// generating multivariate normal random vectors. Uses a Cholesky decomposition. 
// An alternative to .llt() is .ldlt() which works for *semi*definite matrices (is robust in that sense)
MatrixXd mvrNorm(unsigned int N,VecRef mu,MatRef Sig){
	if(mu.size() != Sig.rows() || Sig.rows() != Sig.cols()){
		Rcpp::stop("Need to pass a square covariance matrix the same size as mean vector.");
	}
	MatrixXd Z(rnorm(N,mu.size(),0.,1.));
	const MatrixXd chol( Sig.llt().matrixL() );
	return (chol*Z.transpose()).colwise() + mu;
}

double calcNorm(const int d,VecRef V){
	switch(d){
		case -1: // codes for supremum norm
			return V.lpNorm<Eigen::Infinity>();
		break;
		case 1:
			return V.lpNorm<1>();
		break;
		case 2:
			return V.norm();
		break;
		case 3:
			return V.lpNorm<3>();
		break;
		default:
			return V.norm(); // default is Euclidean length
	}
}

// the bootstrap estimation of the distribution of the d-norm for a given data matrix U
// returns a length-B vector of norms - typically we'd like to use a sample quantile of this vector to set a rejection region.
VectorXd wNorms(const int d,unsigned int B,const MapMxd& U){
	const int N = U.rows(), P = U.cols();
	ArrayXd Utilde(P);
	ArrayXd res(B);
	VectorXd vv(P);
	const ArrayXd V = (U.cwiseProduct(U)).colwise().sum().array().sqrt(); // sample covariances (score functions should be zero-mean by definition, at least under H0)
	for(unsigned int i=0;i<B;i++){
		Utilde = U.transpose()*rademachers(N);
		vv = Utilde/V;
		res(i) = calcNorm(d,vv);
	}
	return res;
}

// this is probably a stupid copy, but it's not clear the relation between map class and matrixBase::derived classes
// possibly, the third argument can be templatized since a MatRef and a MapMxd expose the same methods and operators
ArrayXXd wNorms(const std::vector<int>& d,unsigned int B,MatRef U){
	const int N = U.rows(), P = U.cols();
	ArrayXd Utilde(P);
	ArrayXXd res(B,d.size());
	VectorXd vv(P);
	const ArrayXd V = (U.cwiseProduct(U)).colwise().sum().array().sqrt(); // sample covariances (score functions should be zero-mean by definition, at least under H0)
	for(unsigned int i=0;i<B;i++){
		Utilde = U.transpose()*rademachers(N);
		vv = Utilde/V;
		for(size_t j = 0;j < d.size(); j++){
			res(i,j) = calcNorm(d[j],vv);
		}
	}
	return res;
}

/*
 Given an input matrix U, number of bootstraps B, norm value(s) d, and significance level alpha,
 calculate the sample quantile of the d-norms of the multiplier bootstrap samples
 this function is for "real data analysis" as opposed to simulations.
 This algorithm *does* "studentize" the score vectors in uu, but doesn't center them (naturally)
*/
// [[Rcpp::export]]
Rcpp::List mbTest(Rcpp::NumericMatrix uu,unsigned int B,Rcpp::IntegerVector dd){
	const MapMxd U(Rcpp::as<MapMxd>(uu));
	std::vector<int> d;
	for(int z=0;z<dd.size();z++){
		if(dd[z] >= -1){
			d.push_back(dd[z]);
		}
	}
	std::sort(d.begin(),d.end());
	std::vector<int>::iterator last = std::unique(d.begin(),d.end()); // third argument (not passed here) is a comparison function, which defaults to '<'
	d.erase(last,d.end());
	Rcpp::Rcout << "Norms used:\n";
	for(size_t i=0;i<d.size();i++) Rcpp::Rcout << d[i] << " ";
	Rcpp::Rcout << std::endl;
	const int nNorm = d.size();
	ArrayXXd WW( wNorms(d,B,U) );
	ArrayXd W = U.colwise().sum().array() / (U.cwiseProduct(U)).colwise().sum().array().sqrt(); // scaling by standard deviation
	VectorXd wd(nNorm), pvals(nNorm);
	for(int z=0;z<nNorm;z++){
		wd(z) = calcNorm(d[z],W);
		pvals(z) = (WW.col(z) > wd(z)).count()/double(B);
	}
	// associate the norms with the entries since input order doesn't necessarily match return order
	Rcpp::CharacterVector norm_names(d.size());
	std::string buf;
	for(size_t i=0;i<d.size();i++){
		buf = "<" + std::to_string(d[i]) + ">";
		norm_names(i) = buf;
	}
	Rcpp::NumericVector pv(Rcpp::wrap(pvals));
	pv.attr("names") = norm_names;
	Rcpp::NumericVector obsnorm(Rcpp::wrap(wd));
	obsnorm.attr("names") = norm_names;
	Rcpp::NumericMatrix NN(Rcpp::wrap(WW));
	//NN.attr("colnames") = norm_names;
	colnames(NN) = norm_names;
	return Rcpp::List::create(
		Rcpp::Named("obs.norm") = obsnorm,
		Rcpp::Named("p.value") = pv,
		Rcpp::Named("norms") = NN
	);
}