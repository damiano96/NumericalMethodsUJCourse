CC = g++

n3.x:
	$(CC) n3.cpp -o n3.x

.PHONY: run
run: n3.x
	clear
	./n3.x > output.dat

.PHONY: clean
clean:
	rm -f *.o  *~ *.a *.so *.x core core* a.out *.swp
