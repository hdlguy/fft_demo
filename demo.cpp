#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "fft.hpp"

using std::cout;
using std::endl;
using std::vector;

vector<double> fftshift(vector<double> invec, int n)
{
    vector<double> tempvec(n);
}


int main(int argc, char **argv) {
    const int N  = 4096;
    const double Fs = 327.68e6;
    const double Fc = (+32.0/N)*Fs;
    const double Wc = 2*M_PI*Fc/Fs;
    const int B = 8;
    const double A = pow(2.0,(B-1))-1.0;
    cout << "A =" << A << "\n";

	vector<double> real_vec(N);
	vector<double> imag_vec(N);

    for(int i=0; i<N; i++) {
        real_vec[i] = round(A*cos(Wc*i));
        imag_vec[i] = round(A*sin(Wc*i));
    }

    //for(int i=0; i<16; i++) cout << real_vec[i] << "\n";

	Fft::transform(real_vec, imag_vec);
    vector<double> powvec(N);
    for(int i=0; i<N; i++) powvec[i] = pow(real_vec[i],2.0) + pow(imag_vec[i],2.0);


    FILE* gnuplotPipe = popen ("gnuplot -persistent", "w");
    fprintf(gnuplotPipe, "plot '-' ps 0.4 lt 6\n");
    for (int i=0; i<N; i++) fprintf(gnuplotPipe, "%lf\n", 10*log10(powvec[i]));
    fprintf(gnuplotPipe, "e");
    pclose(gnuplotPipe);

    cout << "test complete\n";
	return EXIT_SUCCESS;
}

