
all:  run
install_dir=./WSN_INSTALL
X_SIZE?=2
Y_SIZE?=3
build: wsn.c
	mkdir -p ./build/$(install_dir)
	cd build && cmake .. -DX_SIZE=$(X_SIZE) -DY_SIZE=$(Y_SIZE) -DCMAKE_INSTALL_PREFIX=$(install_dir) && make install
run: build
	mpirun -n "$$(( $(X_SIZE) * $(Y_SIZE) + 1 ))" --oversubscribe ./build/$(install_dir)/bin/wsn --mca orte_base_help_aggregate 0
	# mpirun -n "$$(( $(X_SIZE) * $(Y_SIZE) + 1 ))"  ./build/$(install_dir)/bin/wsn --mca orte_base_help_aggregate 0

test: test.c cipher.c
	CPATH=$(CPATH):/home/luhanc/scratch/fit3143/fit3143_assignment2/code/src/tiny-AES-c LIBRARY_PATH=$(LIBRARY_PATH):/home/luhanc/scratch/fit3143/fit3143_assignment2/code/src/tiny-AES-c mpicc -g $^ -o $@ -fopenmp -std=c11 

debug: build
	mpirun -n "$$(( $(X_SIZE) * $(Y_SIZE) + 1 ))" --oversubscribe xterm -e lldb ./wsn.c  ./build/$(install_dir)/bin/wsn 
clean:
	rm -rf build *.mpiP test a.out omp_aes

