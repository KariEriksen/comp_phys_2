FLAGS =  -g -fPIC --std=c++11 -larmadillo -O3 -llapack -lblas
depend_all = wavefunc.o vmc.o nqs.o

mh = naive_mh.o
imp = importance.o
gb = gibbs.o

# ==========================================================
# Class Object files
wavefunc.o: ../include/wavefunc.h wavefunc.cpp
	clang wavefunc.cpp -c $(FLAGS) 

nqs.o: ../include/nqs.h nqs.cpp 
	clang nqs.cpp -c $(FLAGS)

vmc.o: ../include/vmc.h vmc.cpp ../include/wavefunc.h wavefunc.cpp
	clang vmc.cpp -std=c++11 -c $(FLAGS) 

importance.o: importance.cpp ../include/vmc.h vmc.cpp 
	clang importance.cpp -c $(FLAGS)

naive_mh.o: naive_mh.cpp ../include/vmc.h vmc.cpp 
	clang naive_mh.cpp -c $(FLAGS)

gibbs.o: gibbs.cpp ../include/vmc.h vmc.cpp 
	clang gibbs.cpp -c $(FLAGS)

# ==========================================================
# Main Object files
main_b.o: main_b.cpp 
	clang main_b.cpp $(FLAGS) -c

main_c.o: main_c.cpp
	clang main_c.cpp $(FLAGS) -c

main_f.o: main_f.cpp
	clang main_f.cpp $(FLAGS) -c
# ==========================================================
# app.x files

app_b.x: main_b.o $(depend_all) $(mh)
	g++ -o app_b.x main_b.o $(depend_all) $(mh) $(FLAGS)  

app_c.x: main_c.o $(depend_all) $(imp)
	g++ -o app_c.x main_c.o $(depend_all) $(imp) $(FLAGS)  

app_f.x: main_f.o $(depend_all) $(gb)
	g++ -o app_f.x main_f.o $(depend_all) $(gb) $(FLAGS)  

clean:
	rm -f *.o
dataclean_b:
	rm ../data/b_data/iteration*
	rm ../data/b_data/time*

dataclean_c:
	rm ../data/c_data/iteration*
	rm ../data/c_data/time*
