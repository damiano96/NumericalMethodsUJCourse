CC = g++

n2.x:
	$(CC) n2.cpp -o n2.x

.PHONY: run
run: n2.x
	clear
	./n2.x > output.dat

.PHONY: clean
clean:
	rm -f *.o  *~ *.a *.so *.x core core* a.out *.swp
