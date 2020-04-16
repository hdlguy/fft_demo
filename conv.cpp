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
    vector<double>::size_type n = vec->size();
    for(int i=0; i<n; i++) (*vec)[i] *= 0.5*(1.0-cos(2.0*M_PI*i/(n-1))); // Hanning
}


int main(int argc, char **argv) {
    const int Nsig  = 256;
    const int Nt    =  32;
    const int Off   =   5;
    const int Nfft  =  (int)pow(2.0, ceil(log2(2.0*Nt-1))); // 2^(ceil(log2(2*Nt-1)));
    cout << "Nfft = " << Nfft << "\n";


    // random signal
    vector<double> s_real(Nsig);
    vector<double> s_imag(Nsig);
    for(int i=0; i<Nsig; i++) {
        s_real[i] = ((double)rand())/((double)RAND_MAX) - 0.5;
        s_imag[i] = ((double)rand())/((double)RAND_MAX) - 0.5;
    }

    // take the reverse order complex conjugate of a portion of the signal to be the convolution template.
    vector<double> h_real(Nt);
    vector<double> h_imag(Nt);
    for(int i=0; i<Nt; i++) {
        h_imag[i] = -s_imag[Nt-1-i+Off];
        h_real[i] = +s_real[Nt-1-i+Off];
    }

    vector<double> h_pad_imag(Nfft, 0.0);
    vector<double> h_pad_real(Nfft, 0.0);
    vector<double> s_pad_imag(Nfft, 0.0);
    vector<double> s_pad_real(Nfft, 0.0);
    for(int i=0; i<Nt; i++) {
        h_pad_imag[i] = h_imag[i];
        h_pad_real[i] = h_real[i];
        s_pad_imag[i] = s_imag[i];
        s_pad_real[i] = s_real[i];
    }
    
    vector<double> H_imag(Nfft); H_imag = h_pad_imag;
    vector<double> H_real(Nfft); H_real = h_pad_real;
    Fft::transformRadix2(H_real, H_imag);

    vector<double> S_imag(Nfft); S_imag = s_pad_imag;
    vector<double> S_real(Nfft); S_real = s_pad_real;
    Fft::transformRadix2(S_real, S_imag);

    vector<double> Y_imag(Nfft);
    vector<double> Y_real(Nfft);
    for(int i=0; i<Nfft; i++) {
        Y_imag[i] = (H_real[i]*S_real[i] - H_imag[i]*S_imag[i]);
        Y_real[i] = (H_real[i]*S_imag[i] + H_imag[i]*S_real[i]);
    }

    vector<double> y_imag(Nfft); y_imag = Y_imag;
    vector<double> y_real(Nfft); y_real = Y_real;
    Fft::transformRadix2(y_imag, y_real);  // real <-> imag swapped to implement inverse fft (without scaling)

    // plot time domain data
    FILE* gnuplotPipe = popen ("gnuplot --persist", "w");
    fprintf(gnuplotPipe, "plot '-' ps 0.1 lt 6\n");
    for (int i=0; i<Nfft; i++) fprintf(gnuplotPipe, "%lf\n", h_pad_real[i]);
    fprintf(gnuplotPipe, "e");
    pclose(gnuplotPipe);


    return EXIT_SUCCESS;
}

