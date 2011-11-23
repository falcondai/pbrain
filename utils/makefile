CFLAGS= -g -DTEST -std=gnu89
OBJS= getstringvalue.o svdcmp.o  \
gaussj.o nrutil.o getreallist.o
OBJS1= fakepoints.o fit-transforms.o

LIBS1 = -lm
LIBS= -lm -ltermcap -lcurses
PROG= procrustes
PROG1= fakepoints

$(PROG): $(PROG).o $(OBJS)
	$(CC) $(LDFLAGS) $(PROG).o $(OBJS) $(LIBS) -o $(PROG)

$(PROG1): $(OBJS1)
	$(CC) $(LDFLAGS)  $(OBJS1) $(LIBS1) -o $(PROG1)
