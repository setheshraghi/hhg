all:
	+$(MAKE) -C .used/cfiles

clean:
	rm .used/cfiles/*.o
