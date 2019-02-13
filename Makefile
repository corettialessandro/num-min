IDIR	= source/headers
CC	= gcc
CFLAGS	= -I$(IDIR) -O3 -Wall
LIBS 	= -lm -llapack

ODIR	= source/obj
SDIR	= source

_DEPS 	= aux_func.h ConjGrad.h Matrix.h ML_SHAKE.h ML_BSHAKE.h output.h ReadIn.h Tensor.h Vector.h
DEPS 	= $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ 	= aux_func.o ConjGrad.o main.o Matrix.o ML_SHAKE.o ML_BSHAKE.o output.o ReadIn.o Tensor.o Vector.o
OBJ 	= $(patsubst %,$(ODIR)/%,$(_OBJ))

EXEC 	= NumMin

$(ODIR)/%.o: $(SDIR)/%.c $(DEPS)
	$(CC) $(CFLAGS) $(LIBS) -c -o $@ $<

$(EXEC): $(OBJ)
	gcc $(CFLAGS) $(LIBS) -o $@ $^

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core $(IDIR)/*~ 
