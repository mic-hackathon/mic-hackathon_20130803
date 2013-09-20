all : dp dest

clean :
	rm -f dp dp-mic

dest : dp-mic
	scp dp-mic /opt/intel/composerxe/lib/mic/libiomp5.so mic1:~/

dp : dp.cpp
	icpc -openmp -vec-report -O3 -std=c++11 dp.cpp -o dp

dp-mic : dp.cpp
	icpc -openmp -vec-report -O3 -std=c++11 dp.cpp -mmic -o dp-mic

dp-mic.s : dp.cpp
	icpc -openmp -std=c++11 -O3 dp.cpp -mmic -S -o dp-mic.s

.PHONY: all clean dest
