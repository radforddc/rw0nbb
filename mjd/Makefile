# Makefile for ORNL MJD scan

CC = gcc
CPP = g++
CFLAGS = -O3 -Wall -std=c99
RM = rm -f

# common files and headers
mk_source_files = sortrun.c decode_runfile_header.c eventbuild.c ep_util.c
mk_header_files = MJDSort.h runBits.h
all_programs = pulser_tag_init presort PZcal skim skim_opt CTcal PSAcal \
               sortrun rundiff rundiff2 siginspect threshold threshold2 \
               sumrms run_check prebuild \
               inlcheck compress decompress deadtime enr_rms

extras = noise coinc 


All: $(all_programs)
extras: $(extras)


PZcal: $(mk_header_files) PZcal.c decode_runfile_header.c ep_util.c pz_fitter.c
	$(CC) $(CFLAGS) -o $@ PZcal.c decode_runfile_header.c ep_util.c pz_fitter.c -lm

skim:  $(mk_header_files) skim.c decode_runfile_header.c ep_util.c pz_fitter.c
	$(CC) $(CFLAGS) -o $@ skim.c decode_runfile_header.c ep_util.c pz_fitter.c -lm

skim_opt:   $(mk_header_files) skim.c decode_runfile_header.c ep_util.c pz_fitter.c
	$(CC) $(CFLAGS) -DOPT -o $@ skim.c decode_runfile_header.c ep_util.c pz_fitter.c -lm

CTcal:  $(mk_header_files) CTcal.c decode_runfile_header.c ep_util.c
	$(CC) $(CFLAGS) -o $@ CTcal.c decode_runfile_header.c ep_util.c -lm

PSAcal: $(mk_header_files) PSAcal.c  decode_runfile_header.c ep_util.c
	$(CC) $(CFLAGS) -o $@ PSAcal.c decode_runfile_header.c ep_util.c -lm

sortrun: $(mk_source_files) $(mk_header_files) eventprocess.c ep_finalize.c pz_fitter.c
	$(CC) $(CFLAGS) -D DO_PSA -o $@ $(mk_source_files) eventprocess.c ep_finalize.c pz_fitter.c -lm

noise: $(mk_source_files) $(mk_header_files) eventprocess.c ep_finalize.c
	$(CC) $(CFLAGS) -D QUIET -o $@ $(mk_source_files) eventprocess.c ep_finalize.c -lm

rundiff: decode_runfile_header.c rundiff.c $(mk_header_files)
	$(CC) $(CFLAGS) -o $@ rundiff.c decode_runfile_header.c -lm

rundiff2: decode_runfile_header.c rundiff2.c $(mk_header_files)
	$(CC) $(CFLAGS) -o $@ rundiff2.c decode_runfile_header.c -lm

pulser_tag_init: $(mk_source_files) ep_pulser_tag_init.c $(mk_header_files)
	$(CC) $(CFLAGS) -o $@ $(mk_source_files) ep_pulser_tag_init.c -lm

siginspect:  $(mk_header_files) siginspect.c decode_runfile_header.c ep_util.c
	$(CC) $(CFLAGS) -o $@ siginspect.c decode_runfile_header.c ep_util.c -lm

presort: $(mk_source_files) $(mk_header_files) presort.c eventprocess.c ep_finalize.c
	$(CC) $(CFLAGS) -D PRESORT -o $@ presort.c decode_runfile_header.c eventbuild.c \
                 eventprocess.c ep_util.c ep_finalize.c -lm

prebuild: $(mk_header_files)  presort.c decode_runfile_header.c ep_util.c prebuild.c
	$(CC) $(CFLAGS) -o $@ presort.c decode_runfile_header.c ep_util.c prebuild.c -lm

sumrms: sumrms.c
	$(CC) $(CFLAGS) -o $@ sumrms.c

threshold: $(mk_header_files) sortrun.c decode_runfile_header.c ep_util.c eb_thresh.c
	$(CC) $(CFLAGS) -o $@ sortrun.c decode_runfile_header.c ep_util.c eb_thresh.c -lm

threshold2: $(mk_header_files) sortrun.c decode_runfile_header.c ep_util.c eb_thresh2.c
	$(CC) $(CFLAGS) -o $@ sortrun.c decode_runfile_header.c ep_util.c eb_thresh2.c -lm

run_check: $(mk_header_files) sortrun.c decode_runfile_header.c ep_util.c run_check.c
	$(CC) $(CFLAGS) -o $@ sortrun.c decode_runfile_header.c ep_util.c run_check.c -lm

inlcheck: $(mk_source_files) ep_inlcheck.c $(mk_header_files)
	$(CC) $(CFLAGS) -o $@ $(mk_source_files) ep_inlcheck.c -lm

compress: $(mk_header_files)  presort.c decode_runfile_header.c ep_util.c compress.c
	$(CC) $(CFLAGS) -o $@ presort.c decode_runfile_header.c ep_util.c compress.c -lm

compress_dptc: $(mk_header_files)  presort.c decode_runfile_header.c ep_util.c compress_dptc.c
	$(CC) $(CFLAGS) -o $@ presort.c decode_runfile_header.c ep_util.c compress_dptc.c dptc_pack.c -lm

decompress: $(mk_header_files)  presort.c decode_runfile_header.c ep_util.c decompress.c
	$(CC) $(CFLAGS) -o $@ presort.c decode_runfile_header.c ep_util.c decompress.c -lm

deadtime: $(mk_source_files) $(mk_header_files) eventprocess.c ep_finalize_DT_thresh.c
	$(CC) $(CFLAGS) -o $@ $(mk_source_files) eventprocess.c ep_finalize_DT_thresh.c -lm

enr_rms: enr_rms.c
	$(CC) $(CFLAGS) -o $@ enr_rms.c

coinc: $(mk_source_files) $(mk_header_files) ep_coinc.c pz_fitter.c
	$(CC) $(CFLAGS) -o $@ $(mk_source_files) ep_coinc.c pz_fitter.c -lm


FORCE:

clean: 
	$(RM) *.o core* *[~%] *.trace
	$(RM) $(all_programs) $(extras)
very-clean: 
	$(RM) *.o core* *[~%] *.trace *.rms run*.pdt* *.dcl* DS* run_data_for_plots.dat TFChannels.txt
	$(RM) $(all_programs) $(extras)
