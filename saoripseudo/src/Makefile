.PHONY: all clean

CFLAGS = -O3

TARGET = phi xeon phi_single xeon_single prefetch
ICPC = icpc
OPENMP = -openmp
MIC = -mmic
IPO = -ipo
FPMODEL = -fp-model fast=1
AVX = -xAVX

all: $(TARGET)
phi: squares_mic.cpp
	$(ICPC) $(CFLAGS) -vec-report $< -lm $(MIC) $(OPENMP) -fno-alias $(FPMODEL) -o $@
xeon: squares_mic.cpp
	$(ICPC) $(CFLAGS) -vec-report $< -lm $(OPENMP)  -fno-alias $(FPMODEL) -o $@
phi_single:squares_single.cpp
	$(ICPC) $(CFLAGS) -vec-report $< -lm $(MIC) -fno-alias $(FPMODEL) -o $@
xeon_single:squares_single.cpp
	$(ICPC) $(CFLAGS) -vec-report $< -lm -fno-alias $(FPMODEL) -o $@
prefetch: squares_prefetch.cpp
	$(ICPC) $(CFLAGS) -vec-report $< -lm $(MIC) $(OPENMP) -fno-alias $(FPMODEL) -o $@

clean:
	-rm -rf $(TARGET)
