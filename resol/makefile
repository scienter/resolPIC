EXEC = show
CC = mpicc

OBJS = main.o boundary.o hdfUtil.o clean.o saveFile.o saveDumpHDF.o restoreDumpHDF.o reCreateParticle.o reHighCreateField.o reLowCreateField.o fieldShareX.o fieldShareY.o changePosition.o
INCL = mesh.h
#CFLAGS= -I/home/scienter/gsl/include
#LDFLAGS= -L/home/scienter/gsl/lib
#INCDIR = /home/scienter/gsl/include /home/scienter/gsl/lib
CFLAGS = -I/mnt/d/Git/gsl2.5/include -I/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include
LDFLAGS = -L/mnt/d/Git/gsl2.5/lib -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi
LIBS = -lm -lhdf5 -lgsl -lgslcblas
$(EXEC):$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -o $(EXEC)
$(OBJS):$(INCL)
clean:
	@rm -f *.o *~ $(EXEC)
