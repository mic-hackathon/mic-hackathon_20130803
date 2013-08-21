########################################################

SRC = main.c
TRG = qn24b_offload

########################################################
$(TRG): $(SRC)
	icc -offload-attribute-target=mic -openmp -O3 $(SRC) -o $(TRG)

clean :
	rm -f *.il *~ core *.o *.s $(TRG)

all:
	$(MAKE) $(TRG)
##### End of Makefile ##################################
