#include <cmath>
#include <cstdlib>
#include <iostream>
#include <vector>
#include "fft.hpp"

using std::cout;
using std::endl;
using std::vector;

void window(vector<double> *vec) 
{
    // Hamming.
    //const double a = 0.54;
    //const double b = 1.0 - a;
    //vector<double>::size_type n = vec->size();
    //for(int i=0; i<n; i++) (*vec)[i] *= a - b*cos(2.0*M_PI*i/(n-1));

    // Hanning
    vector<double>::size_type n = vec->size();
    for(int i=0; i<n; i++) (*vec)[i] *= 0.5*(1.0-cos(2.0*M_PI*i/(n-1)));
}


int main(int argc, char **argv) {
    const int N  = 4096;
    const int P = 16;
    const double Fs = 327.68e6;
    const double Fc = (+500.1/N)*Fs;
    const double Wc = 2*M_PI*Fc/Fs;
    const int B = 10;
    const double A = pow(2.0,(B-1))-1.0;
    cout << "A = " << A << "\n";

    vector<double> real_vec(N);
    vector<double> imag_vec(N);

    // synthesize quantized sinusoids
    for(int i=0; i<N; i++) {
        real_vec[i] = round(A*cos(Wc*i));
        imag_vec[i] = round(A*sin(Wc*i));
    }

    // window the data for spectral analysis
    window(&real_vec);
    window(&imag_vec);

    // zero pad the vectors.
    vector<double> real_vec_pad(N*P, 0.0);
    vector<double> imag_vec_pad(N*P, 0.0);
    for (int i=0; i<N; i++) real_vec_pad[i] = real_vec[i];
    for (int i=0; i<N; i++) imag_vec_pad[i] = imag_vec[i];

    // plot time domain data
    FILE* gnuplotPipe = popen ("gnuplot -persistent", "w");
    fprintf(gnuplotPipe, "plot '-' ps 0.1 lt 6\n");
    for (int i=0; i<N*P; i++) fprintf(gnuplotPipe, "%lf\n", real_vec_pad[i]);
    fprintf(gnuplotPipe, "e");
    pclose(gnuplotPipe);

    // commpute and convert to power
    Fft::transform(real_vec_pad, imag_vec_pad);
    vector<double> powvec(N*P);
    for(int i=0; i<N*P; i++) powvec[i] = pow(real_vec_pad[i]/A,2.0) + pow(imag_vec_pad[i]/A,2.0);

    // Plot frequency domain data in dB
    gnuplotPipe = popen ("gnuplot -persistent", "w");
    fprintf(gnuplotPipe, "plot '-' ps 0.1 lt 6\n");
    for (int i=0; i<N*P; i++) fprintf(gnuplotPipe, "%lf\n", 10*log10(powvec[i]));
    fprintf(gnuplotPipe, "e");
    pclose(gnuplotPipe);

    return EXIT_SUCCESS;
}

