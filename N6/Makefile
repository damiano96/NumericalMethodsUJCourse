CC = g++

n6.x:
	$(CC) n6.cpp -o n6.x

.PHONY: run
run: n6.x
	clear
	./n6.x > output.dat

.PHONY: clean
clean:
	rm -f *.o  *~ *.a *.so *.x core core* a.out *.swp
