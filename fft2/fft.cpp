#include <cmath>
#include <cstdint>
//#include <vector>
#include "fft.hpp"

using std::size_t;
//using std::vector;


static size_t reverseBits(size_t x, unsigned int n); // Private function prototype

void Fft::transform        (fft_float *real, fft_float *imag, int n) { 
    transformRadix2 (real, imag, n); 
}

void Fft::inverseTransform (fft_float *real, fft_float *imag, int n) { 
    transformRadix2 (imag, real, n); 
    for (int i=0; i<n; i++) {
        real[i] /= n;
        imag[i] /= n;
    }
}

void Fft::transformRadix2(fft_float *real, fft_float *imag, int n) {
	// Compute levels = floor(log2(n))
	unsigned int levels;
	{
		size_t temp = n;
		levels = 0;
		while (temp > 1) {
			levels++;
			temp >>= 1;
		}
		if (1u << levels != n)
			throw "Length is not a power of 2";
	}
	
	// Trignometric tables
	fft_float cosTable[n/2];
	fft_float sinTable[n/2];
	for (size_t i = 0; i < n / 2; i++) {
		cosTable[i] = cos(2 * M_PI * i / n);
		sinTable[i] = sin(2 * M_PI * i / n);
	}
	
	// Bit-reversed addressing permutation
	for (size_t i = 0; i < n; i++) {
		size_t j = reverseBits(i, levels);
		if (j > i) {
			fft_float temp = real[i];
			real[i] = real[j];
			real[j] = temp;
			temp = imag[i];
			imag[i] = imag[j];
			imag[j] = temp;
		}
	}
	
	// Cooley-Tukey decimation-in-time radix-2 FFT
	for (size_t size = 2; size <= n; size *= 2) {
		size_t halfsize = size / 2;
		size_t tablestep = n / size;
		for (size_t i = 0; i < n; i += size) {
			for (size_t j = i, k = 0; j < i + halfsize; j++, k += tablestep) {
				fft_float tpre =  real[j+halfsize] * cosTable[k] + imag[j+halfsize] * sinTable[k];
				fft_float tpim = -real[j+halfsize] * sinTable[k] + imag[j+halfsize] * cosTable[k];
				real[j + halfsize] = real[j] - tpre;
				imag[j + halfsize] = imag[j] - tpim;
				real[j] += tpre;
				imag[j] += tpim;
			}
		}
		if (size == n)  // Prevent overflow in 'size *= 2'
			break;
	}
}


static size_t reverseBits(size_t x, unsigned int n) {
	size_t result = 0;
	unsigned int i;
	for (i = 0; i < n; i++, x >>= 1)
		result = (result << 1) | (x & 1);
	return result;
}
