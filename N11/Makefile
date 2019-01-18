CC = g++

n11.x:
	$(CC) n11.cpp -o n11.x

.PHONY: run
run: n11.x
	clear
	./n11.x > dane.dat

.PHONY: clean
clean:
	rm -f *.o  *~ *.a *.so *.x core core* a.out *.swp
