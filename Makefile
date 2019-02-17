IDIR	= source/headers
CC	= gcc
CFLAGS	= -I$(IDIR) -O3 -Wall
LIBS 	= -lm -llapack

ODIR	= source/obj
SDIR	= source

_DEPS 	= Analysis.h aux_func.h ConjGrad.h Matrix.h Minimization.h ML_SHAKE.h ML_BSHAKE.h output.h ReadIn.h Reduced.h setup.h Tensor.h Vector.h
DEPS 	= $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ 	= Analysis.o aux_func.o ConjGrad.o main.o Matrix.o Minimization.o ML_SHAKE.o ML_BSHAKE.o output.o ReadIn.o Reduced.o setup.o Tensor.o Vector.o
OBJ 	= $(patsubst %,$(ODIR)/%,$(_OBJ))

EXEC 	= num-min

$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) $(CFLAGS) $(LIBS) -c -o $@ $<

$(EXEC): $(OBJ)
	gcc $(CFLAGS) $(LIBS) -o $@ $^

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~ 
