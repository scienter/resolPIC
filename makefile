EXEC = show
CC = mpicc
#CC = h5pcc
OBJS = main.o parameterSetting.o findparam.o boundary.o clean.o loadPlasma.o restoreDumpHDF.o loadMovingPlasma.o utility.o saveDumpHDF.o saveFile.o saveFieldHDF.o saveParticleHDF.o saveDensityHDF.o fieldShareX.o fieldShareY.o fieldShareZ.o fieldSolve.o pml.o loadLaser.o interpolation.o particlePush.o ionization.o updateCurrent.o movingDomain.o rearrangeParticles.o particleShareX.o particleShareY.o particleShareZ.o removeEdge.o particlePeriod.o loadBeam.o shotLaser.o plasmaLens.o

#-------- for Beam server ----------# 
CFLAGS = -I/opt/gsl/2.6/include -I/opt/fftw/3.3.8/include -I/opt/hdf5/1.10.3/include 
LDFLAGS = -L/opt/gsl/2.6/lib -L/opt/fftw/3.3.8/lib -L/opt/hdf5/1.10.3/lib

#-------- for PAL ----------#
#CFLAGS = -I/home/scienter/gsl2.6/include -I/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include -I/home/scienter/fftw/include
#LDFLAGS = -L/home/scienter/gsl2.6/lib -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -L/home/scienter/fftw/lib

#-------- for home computer ----------#
#CFLAGS= -I/opt/gsl/2.6/include -I/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include -I/home/scienter/fftw/include
#LDFLAGS= -L/opt/gsl/2.6/lib -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi -L/home/scienter/fftw/lib

#CFLAGS= -I/opt/gsl/2.6/include -I/usr/lib/x86_64-linux-gnu/hdf5/openmpi/include
#LDFLAGS= -L/opt/gsl/2.6/lib -L/usr/lib/x86_64-linux-gnu/hdf5/openmpi


INCL = constants.h mesh.h particle.h laser.h plasma.h
#LIBS = -lm -lgsl -lgslcblas -lhdf5 -lz
LIBS = -lm -lgsl -lgslcblas -lhdf5 -lz -lfftw3

$(EXEC):$(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -o $(EXEC)
#	$(CC)           $(OBJS)            $(LIBS) -o $(EXEC)
$(OBJS):$(INCL)
clean:
	@rm -f *.o *~ $(EXEC)
