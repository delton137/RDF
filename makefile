

COMP = ifort 
FFLAGS =  -O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE_SOURCE  

all: RDF_master.x clean

debug: FFLAGS += -traceback -debug
debug: RDF_master.x clean

RDF_master.x: global_vars.o ACF.o m_RDF.o m_timer.o RDF_master.o
	$(COMP) $(FFLAGS) global_vars.o ACF.o RDF_master.o m_RDF.o m_timer.o libxdrf.a -o RDF_master.x

global_vars.o:
	$(COMP) $(FFLAGS) -c global_vars.f90

ACF.o:
	$(COMP) $(FFLAGS) -c ACF.f90

m_RDF.o:
	$(COMP) $(FFLAGS) -c m_RDF.f90

m_timer.o:
	$(COMP) $(FFLAGS) -c m_timer.f90

RDF_master.o:
	$(COMP) $(FFLAGS) -c RDF_master.f90

clean:
	rm -rf *.o
	rm -rf *.mod
