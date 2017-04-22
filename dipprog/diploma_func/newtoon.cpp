
int NEWTON(double *beta) {
	int i, j;
	double res[N], res_w[N], beta_w[N], dbeta[N], gamma;
	bool flag;
	
	l(beta,res);
	
	for(j=0; j<15; j++) {
		//cout << "\n-------\nj = " << j << "\n-------\n";
	
		if(fabs(res[0])<ens && fabs(res[1])<ens) {
			cout << endl <<"Ended by " << j << " iteration" << endl; 
			return j;
		}
		gradf(beta,dbeta,res);
	
		gamma=1.0;
		flag=false;
		while(gamma>eps_gamma) {
			for(i=0; i<N; i++) beta_w[i]=beta[i]+gamma*dbeta[i]; //+ заменил на -
				
			l(beta_w, res_w);
			
			if(norm(res_w,N)<norm(res,N)) {flag=true; break;}
			gamma/=2.0;
		}
		if(flag==false) {cout << endl << "Broken on " << j << "'s iteration" << endl; return -1;}
		for(i=0; i<N; i++) {
			beta[i]=beta_w[i];
			res[i]=res_w[i];
		}
	}
cout << endl << "not enough iteration" << endl;
return -2;
}
