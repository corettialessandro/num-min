IDIR	= source/headers
CC	= gcc
CFLAGS	= -I$(IDIR) -O3 -Wall
LIBS 	= -lm -llapack

ODIR	= source/obj
SDIR	= source

_DEPS 	= Analysis.h aux_func.h ConjGrad.h Constrained.h Matrix.h ML_SHAKE.h ML_BSHAKE.h output.h ReadIn.h Reduced.h setup.h Tensor.h Unconstrained.h Vector.h
DEPS 	= $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ 	= Analysis.o aux_func.o ConjGrad.o Constrained.o main.o Matrix.o ML_SHAKE.o ML_BSHAKE.o output.o ReadIn.o Reduced.o setup.o Tensor.o Unconstrained.o Vector.o
OBJ 	= $(patsubst %,$(ODIR)/%,$(_OBJ))

EXEC 	= num-min

$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

$(EXEC): $(OBJ)
	gcc $(CFLAGS) $(LIBS) -o $@ $^

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~ 
