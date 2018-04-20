FLAGS =  -g -fPIC --std=c++11 -larmadillo
depend_all = wavefunc.o vmc.o

gin = gaussian_inter_numeric.o
gia = gaussian_inter_analytic.o
gnin = gaussian_noninter_numeric.o
gnia = gaussian_noninter_analytic.o

mh = naive_mh.o
imp = importance.o

wavefunc.o: ../include/wavefunc.h wavefunc.cpp
	clang wavefunc.cpp -c $(FLAGS) 

gaussian_inter_numeric.o: ../include/gaussian_inter_numeric.h gaussian_inter_numeric.cpp
	clang gaussian_inter_numeric.cpp -c $(FLAGS)

gaussian_noninter_numeric.o: ../include/gaussian_noninter_numeric.h gaussian_noninter_numeric.cpp
	clang gaussian_noninter_numeric.cpp -c $(FLAGS) 

gaussian_noninter_analytic.o: ../include/gaussian_noninter_analytic.h gaussian_noninter_analytic.cpp
	clang gaussian_noninter_analytic.cpp -c $(FLAGS) 

gaussian_inter_analytic.o: ../include/gaussian_inter_analytic.h gaussian_inter_analytic.cpp
	clang gaussian_inter_analytic.cpp -c $(FLAGS) 

vmc.o: ../include/vmc.h vmc.cpp ../include/wavefunc.h wavefunc.cpp
	clang vmc.cpp -std=c++11 -c $(FLAGS) 

importance.o: importance.cpp ../include/vmc.h vmc.cpp 
	clang importance.cpp -c $(FLAGS)

naive_mh.o: naive_mh.cpp ../include/vmc.h vmc.cpp 
	clang naive_mh.cpp -c $(FLAGS)

main_b.o: main_b.cpp 
	clang main_b.cpp $(FLAGS) -c

app_b.x: main_b.o $(gnin) $(gnia) $(mh) $(depend_all)
	g++ -o app_b.x main_b.o $(FLAGS) $(gnin) $(gnia) $(mh) $(depend_all)

main_f.o: main_f.cpp $(gia) $(mh) $(depend_all)
	clang main_f.cpp $(FLAGS) -c

app_f.x: main_f.o $(gnia) $(mh) $(depend_all)
	g++ -o app_f.x main_f.o $(FLAGS) $(gnia) $(mh) $(depend_all)

main_test_2.o: main_test_2.cpp 
	clang main_test_2.cpp $(FLAGS) -c

app_test2.o: main_test_2.o $(gia) $(mh) $(depend_all)
	g++ -o test2.x main_test_2.o $(gia) $(mh) $(depend_all) $(FLAGS)
	
clean:
	rm -f *.x

