all:run

exe:main.c
	gcc -o $@ $< -lm

run:exe
	./exe
	matlab -nodesktop -r show_res	

.PHONY:clean
clean:
	-rm *.csv exe
