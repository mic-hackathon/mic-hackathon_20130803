########################################################

SRC = main.c
TRG = qn24b_pthread

########################################################
$(TRG): $(SRC)
	icc -mmic -lpthread -O3 $(SRC) -o $(TRG)

clean :
	rm -f *.il *~ core *.o *.s $(TRG)

all:
	$(MAKE) $(TRG)
##### End of Makefile ##################################
