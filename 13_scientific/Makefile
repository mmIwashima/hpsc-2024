10_cavity: 10_cavity.cpp
	g++ 10_cavity.cpp -march=native -O3
10_cavity_omp: 10_cavity_omp.cpp
	g++ 10_cavity.cpp -march=native -O3 -fopenmp
10_cavity_acc: 10_cavity_acc.cpp
	nvc++ 10_cavity_acc.cpp -O2 -acc -gpu=managed -Minfo
10_cavity_mpi: 10_cavity_mpi.cpp
	mpicxx 10_cavity_mpi.cpp -march=native -O3	
