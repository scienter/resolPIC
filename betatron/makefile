EXEC = show
CC = mpicc
OBJS = main.o findparam.o parameterSetting.o shareData.o solveSpectrum.o
INCL = parameter.h constants.h

CFLAGS = -I/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include
LDFLAGS = -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi
LIBS = -lm -lhdf5
$(EXEC):$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -o $(EXEC)
#	$(CC)           $(OBJS)            $(LIBS) -o $(EXEC)
$(OBJS):$(INCL)
clean:
	@rm -f *.o *~ $(EXEC)
