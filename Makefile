makeT: SO6.cpp Z2.cpp main.cpp
	g++ main.cpp SO6.cpp Z2.cpp pattern.cpp -pthread -O3 -o main.out -fopenmp
