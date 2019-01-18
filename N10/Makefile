CC = g++

N10.x:
	$(CC) N10.cpp -o N10.x

.PHONY: run
run: N10.x
	clear
	./N10.x > output.dat

.PHONY: clean
clean:
	rm -f *.o  *~ *.a *.so *.x core core* a.out *.swp
