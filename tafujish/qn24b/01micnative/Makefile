########################################################

SRC = main.c
TRG = qn24b_micnative

########################################################
$(TRG): $(SRC)
	icc -openmp -mmic -O3 $(SRC) -o $(TRG)

nonmic: $(SRC)
	icc -openmp -O3 $(SRC) -o $(TRG)

nonomp: $(SRC)
	icc -mmic -O3 $(SRC) -o $(TRG)

clean :
	rm -f *.il *~ core *.o *.s $(TRG)

all:
	$(MAKE) $(TRG)
##### End of Makefile ##################################
