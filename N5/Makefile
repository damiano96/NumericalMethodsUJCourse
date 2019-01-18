CC = g++
ARG1 = 1
ARG2 = 2
ARG3 = 3
ARG4 = 4

n5.x:
	$(CC) N5.cpp -o N5.x

run2: n5.x
	clear
	./n5.x $(ARG2) > wyniki/jacobi_method.dat

run3: n5.x
	clear
	./n5.x $(ARG3) > wyniki/gauss_method.dat

run4: n5.x
	clear
	./n5.x $(ARG4) > wyniki/sor_method.dat

run1: n5.x
	clear
	./n5.x $(ARG1) > wyniki/richardson.dat

.PHONY: clean
clean:
	rm -f *.o  *~ *.a *.so *.x core core* a.out *.swp
