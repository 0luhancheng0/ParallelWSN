
all: run
X_SIZE?=4
Y_SIZE?=5
DEBUG?=0
wsn: wsn.c WjCryptLib_Aes.c WjCryptLib_AesCtr.c
	mpicc -DX_SIZE=$(X_SIZE) -DY_SIZE=$(Y_SIZE) -DDEBUG=$(DEBUG) -o $@ $^ -I.
run: wsn
	mpirun -n "$$(( $(X_SIZE) * $(Y_SIZE) + 1 ))" --oversubscribe ./wsn
debug: wsn
	mpirun -n "$$(( $(X_SIZE) * $(Y_SIZE) + 1 ))" --oversubscribe xterm -e lldb -o run ./wsn.c  ./wsn 
clean:
	rm -rf wsn *.mpiP  log.txt


