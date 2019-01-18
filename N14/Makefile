CC = g++

n14.x:
	$(CC) n14.cpp -o n14.x

.PHONY: run
run: n14.x
	clear
	./n14.x > output.dat

.PHONY: clean
clean:
	rm -f *.o  *~ *.a *.so *.x core core* a.out *.swp
