#include <iostream>
#include <cmath>
#include <cstring>
#include <fstream>
#include <cstdio>
#include <iomanip>
#include <cstdlib>

#define EPS_DDOPRI 1e-11
#define ens 1e-10 //epsilon newton stop
#define eps_gamma 2e-15

using namespace std;

//GLOBAL VARIABLES
double g = 9.8;
double XT, YT;
int N=2;
//----------------
//OUTPUT FILES
ofstream coor_out;
ofstream HE_out;
//----------------


#include "diploma_func/utilits.cpp"

//------------------------------------------
//Правые части
#include "diploma_func/fcn.cpp"
//Функция Понтрягина
#include "diploma_func/hpontr.cpp"
//Метод Рунге-Кутты 5-ого порядка с автоматическим выбором шага
#include "diploma_func/ddopri5_brakhistakhrona.cpp"
//Решение системы линейных уравнений методом Гаусса
#include "diploma_func/gauss_m.cpp"
//Метод стрельбы, функции невязок
#include "diploma_func/nev.cpp"
//
#include "diploma_func/gradf.cpp"
//Метод Ньютона
#include "diploma_func/newtoon.cpp"

int main() {
	double beta[2] = {1, 1};
	int check;
	//string str="diploma_output/";
	
	//OPEN FILES
	coor_out.open("diploma_output/coordinates.txt");
	HE_out.open("diploma_output/HE.txt");
	//----------------------------------------------
	cout << "input XT and YT: ";
	cin >> XT >> YT;
	
	cout << "0 for ddopri 5, 1 for NEWTON: ";
	cin >> check;

	if (check==0) {
			double res[2];
			l(beta,res);
	}
	else if (check==1) {
		NEWTON(beta);
		
	}

	cout << "\n\n--------Analitic solution---------\n\n";
	//cout << "v(" << beta[2] << ") = " << sqrt(-2*g*y[2]);
return 0;
}
