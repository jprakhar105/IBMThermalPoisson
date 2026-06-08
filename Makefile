all:
	g++ -O3 -g -static -std=c++17 -Wall main.cpp mesh.cpp flowFunctions.cpp matrixAssembly.cpp writeOutput.cpp linSolve.cpp sparse.cpp -o myprog.exe -lpsapi

clean:
	del /Q myprog.exe
