# Makefile for ORNL MJD scan

CC = gcc
CPP = g++
CFLAGS = -O3 -Wall -std=gnu99
RM = rm -f

# common files and headers
mk_header_files = MJDSort.h
all_programs = PZcal skim skim_opt CTcal PSAcal
#               siginspect


All: $(all_programs)

PZcal:  $(mk_header_files) PZcal.c decode_runfile_header.c ep_util.c pz_fitter.c
	$(CC) $(CFLAGS) -o $@ PZcal.c decode_runfile_header.c ep_util.c pz_fitter.c -lm

skim:   $(mk_header_files) skim.c decode_runfile_header.c ep_util.c pz_fitter.c
	$(CC) $(CFLAGS) -o $@ skim.c decode_runfile_header.c ep_util.c pz_fitter.c -lm

skim_opt:   $(mk_header_files) skim.c decode_runfile_header.c ep_util.c pz_fitter.c
	$(CC) $(CFLAGS) -DOPT -o $@ skim.c decode_runfile_header.c ep_util.c pz_fitter.c -lm

CTcal:  $(mk_header_files) CTcal.c ep_util.c
	$(CC) $(CFLAGS) -o $@ CTcal.c decode_runfile_header.c ep_util.c -lm

PSAcal: $(mk_header_files) PSAcal.c  ep_util.c
	$(CC) $(CFLAGS) -o $@ PSAcal.c decode_runfile_header.c ep_util.c -lm

siginspect: $(mk_header_files) siginspect.c decode_runfile_header.c ep_util.c
	$(CC) $(CFLAGS) -o $@ siginspect.c decode_runfile_header.c ep_util.c -lm

FORCE:

clean: 
	$(RM) *.o core* *[~%] *.trace
	$(RM) $(all_programs) siginspect
very-clean: 
	$(RM) *.o core* *[~%] *.trace *.rms run*.pdt* *.dcl* DS* run_data_for_plots.dat
	$(RM) $(all_programs) siginspect
