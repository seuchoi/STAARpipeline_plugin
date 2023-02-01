// [[Rcpp::depends(RcppArmadillo)]]


#define ARMA_64BIT_WORD 1
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// declare K
double K(double x, arma::vec egvalues);
// declare K1 (first derivative)
double K1(double x, arma::vec egvalues, double q);
// declare K2 (second derivative)
double K2(double x, arma::vec egvalues);
// declare bisection
double Bisection(arma::vec egvalues, double q, double xmin, double xmax);
// declare saddlepoint
double Saddle(double q, arma::vec egvalues);
// declare CCT_pval
double CCT_pval(arma::vec x, arma::vec weights);

// [[Rcpp::export]]
arma::vec Burden_score_SMMAT_sparse(arma::sp_mat G, arma::sp_mat Sigma_i, arma::mat Sigma_iX, arma::mat cov, arma::vec residuals, arma::mat weights_B, arma::mat weights_S, arma::mat weights_A, arma::vec mac, int mac_thres=10)
{
	int i,k;
	int j;
	double sum0 = 0.0;
	double sumw = 0.0;
	double sumx = 0.0;

	bool lower = false;
	bool logp = false;


	// int vn = G.n_rows;
	int un = G.n_cols;
	int vvn = Sigma_iX.n_cols;
	int wn = weights_B.n_cols;

	arma::vec res;
	res.zeros(3*wn);

	arma::mat tSigma_iX_G;
	tSigma_iX_G.zeros(vvn,un);


	arma::mat Cov;
	Cov.zeros(un,un);

	arma::mat Covw;
	Covw.zeros(un,un);

	arma::rowvec x = trans(residuals)*G;

	arma::vec eigenvals;
	eigenvals.zeros(un);


	arma::mat Wleft;
	Wleft.zeros(un,un);

	arma::mat Wright;
	Wright.zeros(un,un);

	double c1 = 0.0;
	double c2 = 0.0;
	double c4 = 0.0;
	double l = 0.0;

	int ii;

	tSigma_iX_G = trans(Sigma_iX)*G;
	Cov = trans(Sigma_i*G)*G - trans(tSigma_iX_G)*cov*tSigma_iX_G;
	// Cov = trans(G)*Sigma_i*G - trans(tSigma_iX_G)*cov*tSigma_iX_G;


	//ACAT
	arma::uvec id_veryrare = arma::find(mac <= mac_thres);
	arma::uvec id_common = arma::find(mac > mac_thres);


	int n0 = id_veryrare.size();
	int n1 = id_common.size();

	arma::vec pseq;
	pseq.zeros(un);

	arma::vec wseq;
	wseq.zeros(un);

	// double SSR = 0.0;
	// double SST = 0.0;

	for(k = 0; k < n1; k++)
	{
		pseq(k) = pow(x(id_common(k)),2)/Cov(id_common(k),id_common(k));
		pseq(k) = R::pchisq(pseq(k),1,lower,logp);
	}




	int n = x.size();
	for(i = 0; i < wn; i++)
	{

		// Burden
		Wright.each_row() = trans(weights_B.col(i));
		Wleft.each_col() = weights_B.col(i);

		Covw = Wleft%Cov%Wright;

		sumw = arma::accu(Covw);

		sum0 = 0.0;
		for (k = 0; k < n; k++)
		{
			sum0 = sum0 + x(k)*weights_B(k,i);
		}

		sumx = pow(sum0, 2) / sumw;
        res(i) = sum0;
        res(wn + i) = sqrt(sumw);
				res(2*wn + i) = R::pchisq(sumx,1,lower,logp);



	}


	return res;
}
