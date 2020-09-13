#include <iostream>
#include <cmath>
#include <complex>
using namespace std;

int main()
{
	for (int c = 50000; c < 50001; c++) {
		for (int d = 0; d < 250; d++) {
			




				int i, j;
				complex<long double> res[4][4];
				complex<long double> ans[4][1];
				double N = c;
				double param = 45;
				double param1 = 70;
				double param2 = 0;
				complex<long double> y (-cos(param2*3.141592653589793238/180), -sin(param2*3.141592653589793238 / 180));
				complex<long double> PI = (3.14159, 0);
				int r = d;
				complex<long double> mati[4][1] = { { complex<long double>(cos(param * 3.141592653589793238 / 180)/ sqrt(N),0) },{ complex<long double>(sin(param * 3.141592653589793238 / 180)*(sqrt(N-2)) / sqrt(N),0) },{ complex<long double>(sqrt(2)*sin(param * 3.141592653589793238 / 180)/sqrt(N), 0)},{ complex<long double>(cos(param * 3.141592653589793238 / 180)*sqrt(N-1)/sqrt(N),0) } };
				complex<long double> matii[4][1] = { { (1 / N) },{ sqrt((N - 1)*(N - 2)) / N },{ sqrt(2 * (N - 1)) / N },{ sqrt(N - 1) / N } };
				complex<long double> res2[4][4];
				complex<long double> mmat1[4][4];
				complex<long double> mmat2[4][4];
				complex<long double> mmmat1[4][4];
				complex<long double> mmmat2[4][4];
				complex<long double> mat1[4][4] = { { complex<long double>(pow(real(y),2.0)-pow(imag(y),2.0),2*real(y)*imag(y)), complex<long double>((pow(real(y),2.0)-pow(imag(y),2.0)-2*real(y)+1)/N,(2*real(y)*imag(y)-2*imag(y))/N), complex<long double>(2*(pow(real(y),2.0)-pow(imag(y),2.0))-2*real(y),2*2*real(y)*imag(y)-2*imag(y)), complex<long double>(pow(real(y),2.0)-pow(imag(y),2.0)-1,2*real(y)*imag(y)) },
				{ 0,1,0,0 },{ 0,complex<long double>((real(y)-1)/N, imag(y)/N), complex<long double>(real(y),imag(y)), 0 },{ 0,0,0,1 } };
				complex<long double> mat2[4][4] = { { 1,0,0,0 },{ 4 / N,1,4,0 },{ -2 / N,0,-1,0 },{ 0,0,0,1 } };
				complex<long double> mat5[4][4] = { {1,1 / N,2,1},{0, sqrt((N -1)*(N - 2)) / N,0,0}, {0,sqrt(2 * (N - 1)) / N,sqrt(2 * (N -1)), 0},{0,sqrt(N - 1) / N,0,sqrt(N - 1)} };
				complex<long double> mat6[4][4] = { {1,2 / sqrt((N - 1)*(N -2)),-sqrt(2) / sqrt(N - 1),-1/ sqrt(N - 1)},{0,N / sqrt((N - 1)*(N -2)),0,0},{0,-1 / sqrt((N -1)*(N - 2)),1 / sqrt(2 * (N - 1)),0},{0,-1 / sqrt((N - 1)*(N -2)),0,1 / sqrt(N-1)} };
				complex<long double> mat7[4][4] = { {complex<long double>(cos(param1 * 3.141592653589793238 / 180),sin(param1 * 3.141592653589793238 / 180)),0,0,0},{0,1,0,0,},{0,0,1,0},{0,0,0,complex<long double>(cos(param1 * 3.141592653589793238 / 180),sin(param1 * 3.141592653589793238 / 180)) } };
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						complex<long double> ini (0,0);
						mmat1[i][j] = ini;
						complex<long double> sum  (0,0);
						for (int k = 0; k < 4; k++) {

							sum += mat1[i][k] * mat6[k][j];
						}
						mmat1[i][j] = sum;
					}
				}
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						mmmat1[i][j] = (0, 0);
						complex<long double> sum  (0, 0);
						for (int k = 0; k < 4; k++) {

							sum += mat5[i][k] * mmat1[k][j];
						}
						mmmat1[i][j] = sum;
					}
				}
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						mmat2[i][j] = (0,0);
						complex<long double> sum  (0,0);
						for (int k = 0; k < 4; k++) {

							sum += mat2[i][k] * mat6[k][j];
						}
						mmat2[i][j] = sum;
					}
				}
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						mmmat2[i][j] = (0, 0);
						complex<long double> sum  (0,0);
						for (int k = 0; k < 4; k++) {

							sum += mat5[i][k] * mmat2[k][j];
						}
						mmmat2[i][j] = sum;
					}
				}
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						res[i][j] = (0, 0);
						complex<long double> sum  (0,0);
						for (int k = 0; k < 4; k++) {

							sum += mmmat2[i][k] * mmmat1[k][j];
						}
						res[i][j] = sum;
					}
				}
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 4; j++) {
						res2[i][j] = (0, 0);
						complex<long double> sum  (0, 0);
						for (int k = 0; k < 4; k++) {

							sum += mat7[i][k] * res[k][j];
						}
						res2[i][j] = sum;
					}
				}
				complex<long double> res1[4][4] = { { 1,0,0,0 },{ 0,1,0,0 },{ 0,0,1,0 },{ 0,0,0,1 } };
				for (int l = 0; l < r; l++) {
					for (int i = 0; i < 4; i++) {
						for (int j = 0; j < 4; j++) {
							complex<long double> sum  (0,0);
							for (int k = 0; k < 4; k++) {
								sum += res1[i][k] * res2[k][j];
							}
							res1[i][j] = sum;
						}
					}
				}
				for (int i = 0; i < 4; i++) {
					for (int j = 0; j < 1; j++) {
						ans[i][j] = (0,0);
						for (int k = 0; k < 4; k++) {
							ans[i][j] += res1[i][k] * mati[k][j];
						}
					}
				}
				
				cout << "(" << d <<  "," << pow(abs(ans[0][0]),2.0) + pow(abs(ans[0][2]),2.0) << ")" <<  " ";
				
				cout << "\n";
			
		}
	}
	return 0;
}


