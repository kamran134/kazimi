double HPONTR(double *y, double t) {
	double V,Px,Py,Pv;
	double C1, C2, RHO;
	double cosT,sinT;
	//--------------------
/*
	X = y[0];
	Y = y[1];
*/
	V = y[2];
	Px = y[3];
	Py = y[4];
	Pv = y[5];
	//--------------------
	C1 = Px*V;
	C2 = Py*V-g*Pv;
	RHO = sqrt(C1*C1+C2*C2);
	cosT = C1/RHO;
	sinT = C2/RHO;
	return (Px*V*cosT + (Py*V - g*Pv)*sinT);
}
