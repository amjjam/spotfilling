CPPFLAGS=-Wall -g -I ../submodules/include/
CPP=g++

build: runDGCPM

runDGCPM: runDGCPM.o spotfilling.o
	$(CPP) -o $@ $^ -I ../submodules/include \
	-L ../submodules/lib -lDGCPM -laTime -lkp -lgfortran -lz -lrt

clean:
	- rm -f runDGCPM.o spotfilling.o

