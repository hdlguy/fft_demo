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
    const int Nsig  = 256;
    const int Nt    =  32;
    const int Off   =   5;
    const int Nfft  =  (int)pow(2.0, ceil(log2(2.0*Nt-1))); // 2^(ceil(log2(2*Nt-1)));


    // random signal
    vector<double> s_real(Nsig);
    vector<double> s_imag(Nsig);
    for(int i=0; i<Nsig; i++) {
        s_real[i] = ((double)rand())/((double)RAND_MAX) - 0.5;
        s_imag[i] = ((double)rand())/((double)RAND_MAX) - 0.5;
    }

    // take a portion of the signal to be the convolution template.
    vector<double> h_real(Nt);
    vector<double> h_imag(Nt);
    for(int i=0; i<Nt; i++) {
        h_real[i] = s_real[i+Off];
        h_imag[i] = s_imag[i+Off];
    }

/*
    // zero pad the vectors.
    vector<double> real_vec_pad(Nfft, 0.0);
    vector<double> imag_vec_pad(Nfft, 0.0);
    for (int i=0; i<N; i++) real_vec_pad[i] = real_vec[i];
    for (int i=0; i<N; i++) imag_vec_pad[i] = imag_vec[i];

    // plot time domain data
    FILE* gnuplotPipe = popen ("gnuplot -persistent", "w");
    fprintf(gnuplotPipe, "plot '-' ps 0.1 lt 6\n");
    for (int i=0; i<N*P; i++) fprintf(gnuplotPipe, "%lf\n", real_vec_pad[i]);
    fprintf(gnuplotPipe, "e");
    pclose(gnuplotPipe);


    // Plot frequency domain data in dB
    gnuplotPipe = popen ("gnuplot -persistent", "w");
    fprintf(gnuplotPipe, "plot '-' ps 0.1 lt 6\n");
    for (int i=0; i<N*P; i++) fprintf(gnuplotPipe, "%lf\n", 10*log10(powvec[i]));
    fprintf(gnuplotPipe, "e");
    pclose(gnuplotPipe);
*/

    return EXIT_SUCCESS;
}

