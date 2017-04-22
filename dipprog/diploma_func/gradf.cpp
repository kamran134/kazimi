void gradf(double *beta, double *dbeta, double *res) {
	int i, j;
	double beta_p[N], beta_m[N], res_p[N], res_m[N], res_minus[N], h=1.e-6;
	double **pro;
	
	
	pro=(double**)malloc(N*sizeof(double*));
	for(i=0; i<N; i++) {
		pro[i]=(double*)malloc(N*sizeof(double));
	}

	for(j=0; j<N; j++) {
		for(i=0; i<N; i++) {
			beta_p[i]=beta[i];
			beta_m[i]=beta[i];
		}
		
		beta_p[j]+=h;
		beta_m[j]-=h;
		
		l(beta_p, res_p);
		l(beta_m, res_m);
		
		for(i=0; i<N; i++) pro[i][j]=(res_p[i]-res_m[i])/(2*h);
	}
	//matrix_print_f(pro, N);
	
	for(i=0; i<N; i++) res_minus[i]=-res[i];
	linear_sys(pro, dbeta, res_minus, N);
	
	for(i=0; i<N; i++) printf("dbeta: %.16f\n",dbeta[i]);
}
