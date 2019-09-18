all: wsn

ROOT_DIR=.
CC=/usr/local/bin/mpicc
SRC=./src
FLAGS=-fopenmp -I./src --std=c99
.PHONY: run
run: wsn 
	mpirun -n 7 --hostfile ./host ./wsn --mca orte_base_help_aggregate 0
wsn: wsn.c
	$(CC) $< -o $@ $(FLAGS)

test_ctr: $(SRC)/ctr.c 
	gcc $^ -o $@
clean: 
	rm -f wsn test_aes test_ctr
