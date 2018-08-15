# Makefile created by myfmkmf.pl at 09:46:00 on Jan 5, 2012
#
# Modification of fmkmf.pl by Rick D. Saylor, NOAA/ARL/ATDD
#
FC=gfortran 
CFLAGS=-O3 
LFLAGS=  

# create executable
canacc: GlobalData.o EnvironData.o PhysChemData.o Output.o Utils.o VertTransport.o Initialize.o CanopyPhysics.o Main.o 
	 $(FC) -o ../bin/canacc $(LFLAGS) GlobalData.o EnvironData.o PhysChemData.o Output.o Utils.o VertTransport.o Initialize.o CanopyPhysics.o Main.o

# compilation dependencies

GlobalData.o:./GlobalData.f90  
	 $(FC) $(CFLAGS) -c ./GlobalData.f90 

PhysChemData.o:./PhysChemData.f90 GlobalData.o
	 $(FC) $(CFLAGS) -c ./PhysChemData.f90 

Utils.o:./Utils.f90 GlobalData.o PhysChemData.o
	 $(FC) $(CFLAGS) -c ./Utils.f90 

Output.o:./Output.f90 GlobalData.o Utils.o
	 $(FC) $(CFLAGS) -c ./Output.f90 

VertTransport.o:./VertTransport.f90 GlobalData.o Utils.o
	 $(FC) $(CFLAGS) -c ./VertTransport.f90 

CanopyPhysics.o:./CanopyPhysics.f90 GlobalData.o PhysChemData.o Utils.o
	 $(FC) $(CFLAGS) -c ./CanopyPhysics.f90 

EnvironData.o:./EnvironData.f90 GlobalData.o PhysChemData.o Utils.o
	 $(FC) $(CFLAGS) -c ./EnvironData.f90 

Initialize.o:./Initialize.f90 GlobalData.o EnvironData.o PhysChemData.o Utils.o Output.o 
	 $(FC) $(CFLAGS) -c ./Initialize.f90 

Main.o:Main.f90 GlobalData.o Initialize.o EnvironData.o CanopyPhysics.o Output.o Utils.o 
	 $(FC) $(CFLAGS) -c Main.f90 

clean:
	rm -f -r *.o *.mod 

veryclean:
	rm -f -r *.o *.mod ../bin/canacc
