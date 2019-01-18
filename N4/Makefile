CC = g++

n4.x:
	$(CC) n4.cpp -o n4.x

.PHONY: run
run: n4.x
	clear
	./n4.x > dane.dat

.PHONY: clean
clean:
	rm -f *.o  *~ *.a *.so *.x core core* a.out *.swp
