
all: test
install_dir=./WSN_INSTALL
X_SIZE?=2
Y_SIZE?=2
build: 
	mkdir -p ./build/$(install_dir)
	cd build && cmake .. -DX_SIZE=$(X_SIZE) -DY_SIZE=$(Y_SIZE) -DCMAKE_INSTALL_PREFIX=$(install_dir) && make install
run: build
	mpirun -n "$$(( $(X_SIZE) * $(Y_SIZE) + 1 ))" --oversubscribe ./build/$(install_dir)/bin/wsn --mca btl_vader_backing_directory /tmp --mca orte_base_help_aggregate 0

test: test.c
	mpicc $< -o $@
clean:
	rm -rf build *.mpiP test

