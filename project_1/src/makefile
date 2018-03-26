app_GIN.x: main_GIN.o
	g++ -o app_GIN.x main_GIN.o vmc.o wavefunc.o gaussian_inter_numeric.o naive_mh.o importance.o 

app.x: main_prod.o 
	g++ -o app.x main_prod.o vmc.o wavefunc.o gaussian_noninter_numeric.o gaussian_noninter_analytic.o naive_mh.o importance.o -O3 -larmadillo

wavefunc.o: ../include/wavefunc.h wavefunc.cpp
	clang wavefunc.cpp -c -fPIC --std=c++11

gaussian_noninter_numeric.o: ../include/gaussian_noninter_numeric.h gaussian_noninter_numeric.cpp
	clang gaussian_noninter_numeric.cpp -c -O3 -fPIC --std=c++11 -I ~/armadillo/armadillo-8.400.0/include

gaussian_noninter_analytic.o: ../include/gaussian_noninter_analytic.h gaussian_noninter_analytic.cpp
	clang gaussian_noninter_analytic.cpp -c -O3 -fPIC --std=c++11 -I ~/armadillo/armadillo-8.400.0/include

vmc.o: ../include/vmc.h vmc.cpp ../include/wavefunc.h wavefunc.cpp
	clang vmc.cpp -std=c++11 -c -O3 -fPIC -I ~/armadillo/armadillo-8.400.0/include

importance.o: importance.cpp ../include/vmc.h vmc.cpp 
	clang importance.cpp -c -O3 -fPIC -std=c++11

naive_mh.o: naive_mh.cpp ../include/vmc.h vmc.cpp 
	clang naive_mh.cpp -c -fPIC -std=c++11

main_prod.o: main_prod.cpp vmc.o gaussian_noninter_numeric.o gaussian_noninter_analytic.o wavefunc.o importance.o naive_mh.o
	clang main_prod.cpp -O3 -fPIC -std=c++11 -I ~/armadillo/armadillo-8.400.0/include -c

main_GIN.o: main_GIN.cpp vmc.o gaussian_inter_numeric.o wavefunc.o importance.o naive_mh.o
	clang main_GIN.cpp -fPIC -std=c++11 -I ~/armadillo/armadillo-8.400.0/include -c

clean:
	\rm *.o
	\rm ../data/*
