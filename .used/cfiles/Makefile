CFLAGS = -Wall -Wextra -c

all: datagen verify transmission

datagen: datagen.o prop.o prop.h
	gcc -o ../otherexecs/datagen datagen.o prop.o

verify: verify.o prop.o prop.h
	gcc -o ../otherexecs/verify verify.o prop.o

transmission: transmission.o prop.o prop.h
	gcc -o ../otherexecs/transmission transmission.o prop.o

datagen.o: datagen.c prop.h
	gcc ${CFLAGS} datagen.c

verify.o: verify.c prop.h
	gcc ${CFLAGS} verify.c

transmission.o: transmission.c prop.h
	gcc ${CFLAGS} transmission.c

prop.o: prop.h
	gcc ${CFLAGS} prop.c
