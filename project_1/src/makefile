app.x: main_test.o 
	g++ -o app.x main_test.o vmc.o wavefunc.o gaussian_noninter.o naive_mh.o importance.o

wavefunc.o: ../include/wavefunc.h wavefunc.cpp
	clang wavefunc.cpp -c --std=c++11

gaussian_noninter.o: ../include/gaussian_noninter.h gaussian_noninter.cpp
	clang gaussian_noninter.cpp -c --std=c++11 ~/armadillo/armadillo-8.400.0/include

vmc.o: ../include/vmc.h vmc.cpp ../include/wavefunc.h wavefunc.cpp
	clang vmc.cpp -std=c++11 -c -I ~/armadillo/armadillo-8.400.0/include

importance.o: importance.cpp ../include/vmc.h vmc.cpp 
	clang importance.cpp -c -std=c++11

naive_mh.o: naive_mh.cpp ../include/vmc.h vmc.cpp 
	clang naive_mh.cpp -c -std=c++11

main_test.o: main_test.cpp vmc.o gaussian_noninter.o wavefunc.o gaussian_noninter.o  importance.o naive_mh.o
	clang main_test.cpp -std=c++11 -I ~/armadillo/armadillo-8.400.0/include -c

