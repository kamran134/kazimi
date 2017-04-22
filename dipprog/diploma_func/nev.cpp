void l(double *beta, double *res) {
	double HP;
	double y[6];
	
	y[0] = 0; 			// x(0)
	y[1] = 0; 			// y(0)
	y[2] = 0;		 	// v(0)
	y[3] = beta[0];	 	//px(0)
	y[4] = beta[1]; 	//py(0)
	y[5] = 1;			//pv(0)

	HP = HPONTR(y,0);
	ddopri5(6,fcn,0,y,10.e0,1.e-11,0.01e0,0.01e0);
	
	res[0] = y[0]-XT;
	res[1] = y[1]-YT;
	res[2] = y[2];
}
