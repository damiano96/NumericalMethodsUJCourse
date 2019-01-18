CC = g++

n1.x:
	$(CC) n1.cpp -o n1.x

.PHONY: run
run: n1.x
	clear
	./n1.x > output.txt

.PHONY: clean
clean:
	rm -f *.o  *~ *.a *.so *.x core core* a.out *.swp
