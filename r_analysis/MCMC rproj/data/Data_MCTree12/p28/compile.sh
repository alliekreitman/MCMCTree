/usr/bin/gcc -Wall -I/usr/include -I/usr/lib -I/usr/local/include -c $1.c -o $1.o
/usr/bin/gcc -L/usr/local/lib $1.o -lgsl -lgslcblas -o $1.e
