FLAGS = -O3 -fPIC --std=c++11 

app_b.x: main_b.o 
	g++ -o app_b.x main_b.o vmc.o wavefunc.o gaussian_noninter_numeric.o gaussian_noninter_analytic.o naive_mh.o -O3

app_c.x: main_c.o 
	g++ -o app_c.x main_c.o vmc.o wavefunc.o gaussian_noninter_numeric.o gaussian_noninter_analytic.o importance.o -O3 -larmadillo

app_f.x: main_f.o
	g++ -o app_f.x main_f.o vmc.o wavefunc.o gaussian_noninter_analytic.o naive_mh.o

wavefunc.o: ../include/wavefunc.h wavefunc.cpp
	clang wavefunc.cpp -c $(FLAGS) 

gaussian_noninter_numeric.o: ../include/gaussian_noninter_numeric.h gaussian_noninter_numeric.cpp
	clang gaussian_noninter_numeric.cpp -c $(FLAGS) 

gaussian_noninter_analytic.o: ../include/gaussian_noninter_analytic.h gaussian_noninter_analytic.cpp
	clang gaussian_noninter_analytic.cpp -c $(FLAGS) 

vmc.o: ../include/vmc.h vmc.cpp ../include/wavefunc.h wavefunc.cpp
	clang vmc.cpp -std=c++11 -c $(FLAGS) 

importance.o: importance.cpp ../include/vmc.h vmc.cpp 
	clang importance.cpp -c $(FLAGS)

naive_mh.o: naive_mh.cpp ../include/vmc.h vmc.cpp 
	clang naive_mh.cpp -c $(FLAGS)

main_b.o: main_b.cpp vmc.o gaussian_noninter_numeric.o gaussian_noninter_analytic.o wavefunc.o naive_mh.o
	clang main_b.cpp $(FLAGS)  -c

main_c.o: main_c.cpp vmc.o gaussian_noninter_numeric.o gaussian_noninter_analytic.o wavefunc.o importance.o
	clang main_c.cpp $(FLAGS)  -c

main_f.o: main_f.cpp vmc.o gaussian_noninter_analytic wavefunc.o naive_mh.o
	clang main_f.cpp --std=c++11 -c

clean:
	rm *.o

