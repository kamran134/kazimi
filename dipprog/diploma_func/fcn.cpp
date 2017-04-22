void fcn(double x, double *y, double *f) {
	double RHO, C1, C2, sinT, cosT;
	double V, Px, Py, Pv;
	
	//------------------------------
/*
	X = y[0];
	Y = y[1];
*/
	V = y[2];
	Px = y[3];
	Py = y[4];
	Pv = y[5];
	//------------------------------
	C1 = Px*V;
	C2 = Py*V-g*Pv;
	RHO = sqrt(C1*C1+C2*C2);
	cosT = C1/RHO;
	sinT = C2/RHO;
	//------------------------------
	
	f[0] = V*cosT; //x.
	f[1] = V*sinT; //y.
	f[2] = -g*sinT; //v.
	f[3] = 0; //px.
	f[4] = 0; //py.
	f[5] = -Px*cosT-Py*sinT; //pv.
}
