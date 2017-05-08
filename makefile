
all: build

build:
	mpiCC filtru.cpp -o filtru

clean:
	rm filtru
