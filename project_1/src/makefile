main_test.o: main_test.cpp vmc.h
	g++ main_test.cpp -std=c++11 -c

vmc.o: vmc.h vmc.cpp
	g++ vmc.cpp -std=c++11 -c

app.x: main_test.o vmc.o
	g++ -o app.x main_test.o vmc.o 
