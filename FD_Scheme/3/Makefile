all:run

exe:main.c
	gcc -o $@ $< -lm

run:exe
	./exe
	./visual_gnuplot	

.PHONY:clean
clean:
	-rm res exe
