//#pragma once

typedef double fft_float;



	 //* Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 //* The vector can have any length. This is a wrapper function.
	void transform          (fft_float *real, fft_float *imag, int n);
    void inverseTransform   (fft_float *real, fft_float *imag, int n);
	
	 //* Computes the discrete Fourier transform (DFT) of the given complex vector, storing the result back into the vector.
	 //* The vector's length must be a power of 2. Uses the Cooley-Tukey decimation-in-time radix-2 algorithm.
	void transformRadix2(fft_float *real, fft_float *imag, int n);
	
