TARGET = gngosa
OBJS = main.o Gauss_Newton.o gauss.o
CFLAGS = -I$(INCLUDEDIR) -g -Wall -O2
LIBS = -L$(LIBDIR) -lnd_malloc -limageio -lm 
CC = gcc
.c.o:
	$(CC) $(CFLAGS) -c $<

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $(TARGET) $(OBJS) $(LIBS)
#	install $(TARGET) sss

clean:
	rm -f *.o $(TARGET) core *~

#main.o : func.h
#main.o : lpf.h
#func.o : func.h
