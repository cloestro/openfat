CC=gcc
DEBUG=-g
LDFLAGS=-lm
CFLAGS=-g -Wall -DDEBUG
AR=ar rsv

progtest: progtest.o libopenfat.a
	$(CC) -o $@ $(CFLAGS) $^ $(LDFLAGS)
libopenfat.a: matrix.o utilities.o sines.o crossland.o
	$(AR) $@ $^
progtest.o: progtest.c
	$(CC) $(CFLAGS) -c $<

matrix.o: matrix.c
	$(CC) $(CFLAGS) -c $<

sines.o: sines.c
	$(CC) $(CFLAGS) -c $<
crossland.o: crossland.c
	$(CC) $(CFLAGS) -c $<

utilities.o: utilities.c
	$(CC) $(CFLAGS) -c $<
clean:
	rm *.o progtest *.a

