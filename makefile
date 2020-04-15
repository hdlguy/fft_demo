all: demo conv

conv: conv.cpp fft.cpp fft.hpp
	g++ -std=gnu++11 -o conv fft.cpp conv.cpp

demo: demo.cpp fft.cpp fft.hpp
	g++ -std=gnu++11 -o demo fft.cpp demo.cpp

clean: 
	rm demo conv
