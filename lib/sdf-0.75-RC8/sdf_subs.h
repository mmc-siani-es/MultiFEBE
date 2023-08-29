/*
 Simple Data Format
 http://solarmuri.ssl.berkeley.edu/~fisher/public/software/SDF
 Copyright (C) 2006,2007 University of California

 This is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation;
 either version 2.1 of the License, or (at your option) any later version.

 This software is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 See the GNU Lesser General Public License for more details.

 To view the GNU Lesser General Public License visit
 http://www.gnu.org/copyleft/lesser.html
 or write to the Free Software Foundation, Inc.,
 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* Next line is to keep msvc compiler from complaining about ANSI C functions
 * that they claim are deprecated */

# define _CRT_SECURE_NO_DEPRECATE 1

/* Next line is to tell IBM SP compiler to use 64-bit integers in file i/o */
# define _LARGE_FILES 1

/* next test to see if INTF_64 compiler flag was set with gcc */
#ifdef INTF_64
#    define INTF_64 1
#else
#    undef INTF_64
#endif

/* next test to see if INTF_CLEN_64 compiler flag was set with gcc */

/*
#ifdef INTF_CLEN_64
#    define INTF_CLEN_64 1
#else
#    undef INTF_CLEN_64
#endif
*/


#ifndef HINITSZ
#    define HINITSZ 2000
#endif /* HINITSZ */

/* To hard-wire setting default fortran integers to be 64-bits, 
   uncomment the next line */
/*  #define INTF_64 1 */

/* To hard-wire length of fortran string dimension variable to be 64-bits, 
   uncomment the next line */
/* #define INTF_CLEN_64 1 */
# include "./strlenfmt.h"
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <ctype.h>
# include <errno.h>

/* Note that in mingw, WIN32 is defined already as an environment variable,
 * but in msvc, it is only defined within the Makefile.  In other words, if you
 * use the cl compiler alone, WIN32 won't be set and things will go down the
 * unix path.  But if you use nmake things will work OK. */

/* following defines library for "ftruncate" or "SetEndOfFIle" (in MS win): */

#ifdef WIN32
#include <windows.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <io.h>
#else
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif

/* START of function prototypes: */

#ifdef WIN32
	typedef struct _stati64 sdf_stat;
	typedef fpos_t pos; /* in mswin (xp at least) fpos_t is 8-byte int */
	pos ftello(FILE *fd); /* fake ftello */
	int fseeko(FILE *fd, pos offset, int whence); /*fake fseeko */
	char *strtok_r(char *s1, char *sep, char **last); /* fake strtok_r */
	/* actual fake (!) definitions for ftello, fseeko, strtok_r for mswin
	 * occurs toward end of sdf_subs.c*/
         
        /* ftruncate now replaced by uniform syntax sdf_file_truncate, at end
         of sdf_subs.c file */
#else
	typedef struct stat sdf_stat;
	typedef off_t pos; /* file position type in modern *nix */
#endif

typedef int i4;
typedef unsigned int u4;
typedef float f4;
typedef double f8;
typedef long long i8;
typedef short i2; 
#ifdef INTF_64
	typedef i8 intf;
#else
	typedef int intf; 
#endif

#ifdef INTF_CLEN_64
	typedef size_t intf_clen;
#else
	typedef int intf_clen; 
#endif
typedef pos intf8; /* this should always correspond to 64-bit Fortran integer */

/* the following typdef for structure type data_id must occur *after* 
 * typedef for pos and i4 but *before* the function prototype definition in 
 * data_to_string */

struct data_header
{
	i4 order;
	char *label;
	char datatype;
	i4 nbpw;
	i4 ndim;
	pos *dims;
}; 

typedef struct data_header data_id;

int sdf_file_truncate(FILE *fd, pos fs);
int gen_stat(char *fname, sdf_stat *fileinfo);
int data_to_string(data_id * id, FILE *fp);
int data_to_buff(data_id *id, char *buff, pos *offset);
data_id * string_to_data(FILE *fp);
i4 sdf_write(char *fname, data_id *id, void * data);
void sdf_write_f77__(char *fnamef77, char *labf77, char *dtf77, intf *nbpwf77, 
   intf *ndimf77, intf8 *dimsf77, void *dataf77,  intf_clen len_fnamef77, 
   intf_clen len_labf77, intf_clen len_dtf77, intf_clen len_dataf77 );
void sdf_write_f77_(char *fnamef77, char *labf77, char *dtf77, intf *nbpwf77, 
   intf *ndimf77, intf8 *dimsf77, void *dataf77,  intf_clen len_fnamef77, 
   intf_clen len_labf77, intf_clen len_dtf77, intf_clen len_dataf77 );
void sdf_write_f77(char *fnamef77, char *labf77, char *dtf77, intf *nbpwf77, 
   intf *ndimf77, intf8 *dimsf77, void *dataf77,  intf_clen len_fnamef77, 
   intf_clen len_labf77, intf_clen len_dtf77, intf_clen len_dataf77 );
void sdf_fl2ds_f77__(char * sdf_fnamef77, char *fnamef77, char *labelf77,
      intf_clen len_sdf_fnamef77, intf_clen len_fnamef77, intf_clen
      len_labelf77);
void sdf_fl2ds_f77_(char * sdf_fnamef77, char *fnamef77, char *labelf77,
      intf_clen len_sdf_fnamef77, intf_clen len_fnamef77, intf_clen
      len_labelf77);
void sdf_fl2ds_f77(char * sdf_fnamef77, char *fnamef77, char *labelf77,
      intf_clen len_sdf_fnamef77, intf_clen len_fnamef77, intf_clen
      len_labelf77);
i4 sdf_fl2ds(char * sdf_fname, char *fname, char *label);
void sdf_ds2fl_f77__(char *sdf_fnamef77, intf *orderf77, char *fnamef77,
      intf_clen len_sdf_fnamef77, intf_clen len_fnamef77);
void sdf_ds2fl_f77_(char *sdf_fnamef77, intf *orderf77, char *fnamef77,
      intf_clen len_sdf_fnamef77, intf_clen len_fnamef77);
void sdf_ds2fl_f77(char *sdf_fnamef77, intf *orderf77, char *fnamef77,
      intf_clen len_sdf_fnamef77, intf_clen len_fnamef77);
i4 sdf_ds2fl(char *sdf_fname, i4 order, char *fname);
void sdf_read_f77__(char *fnamef77, intf *iorder_f77, char *labf77, 
     char *dtf77, intf *nbpwf77, intf *ndimf77, intf8 *dimsf77, void *dataf77, 
     intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77, 
     intf_clen len_dataf77);
void sdf_read_f77_(char *fnamef77, intf *iorder_f77, char *labf77, 
     char *dtf77, intf *nbpwf77, intf *ndimf77, intf8 *dimsf77, void *dataf77, 
     intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77, 
     intf_clen len_dataf77);
void sdf_read_f77(char *fnamef77, intf *iorder_f77, char *labf77, 
     char *dtf77, intf *nbpwf77, intf *ndimf77, intf8 *dimsf77, void *dataf77, 
     intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77, 
     intf_clen len_dataf77);
void sdf_delete_f77__(char *fnamef77, intf* idelete, intf_clen len_fnamef77);
void sdf_delete_f77_(char *fnamef77, intf* idelete, intf_clen len_fnamef77);
void sdf_delete_f77(char *fnamef77, intf* idelete, intf_clen len_fnamef77);
void sdf_insert_f77__(char *fnamef77, intf* insert, char *labf77, char *dtf77, 
    intf *nbpwf77, intf *ndimf77, intf8 *dimsf77, void *dataf77, 
    intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77, 
    intf_clen len_dataf77 );
void sdf_insert_f77_(char *fnamef77, intf* insert, char *labf77, char *dtf77, 
    intf *nbpwf77, intf *ndimf77, intf8 *dimsf77, void *dataf77, 
    intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77, 
    intf_clen len_dataf77 );
void sdf_insert_f77(char *fnamef77, intf* insert, char *labf77, char *dtf77, 
    intf *nbpwf77, intf *ndimf77, intf8 *dimsf77, void *dataf77, 
    intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77, 
    intf_clen len_dataf77 );
void sdf_replace_f77__(char *fnamef77, intf* replace, char *labf77, 
       char *dtf77, intf *nbpwf77, intf *ndimf77, intf8 *dimsf77, 
       void *dataf77,  intf_clen len_fnamef77, intf_clen len_labf77, 
       intf_clen len_dtf77, intf_clen len_dataf77 );
void sdf_replace_f77_(char *fnamef77, intf* replace, char *labf77, 
       char *dtf77, intf *nbpwf77, intf *ndimf77, intf8 *dimsf77, 
       void *dataf77,  intf_clen len_fnamef77, intf_clen len_labf77, 
       intf_clen len_dtf77, intf_clen len_dataf77 );
void sdf_replace_f77(char *fnamef77, intf* replace, char *labf77, 
       char *dtf77, intf *nbpwf77, intf *ndimf77, intf8 *dimsf77, 
       void *dataf77,  intf_clen len_fnamef77, intf_clen len_labf77, 
       intf_clen len_dtf77, intf_clen len_dataf77 );
void sdf_replace_id_f77__(char *fnamef77, intf * replace, char *labf77,
       char *dtf77, intf  *nbpwf77, intf  *ndimf77, intf8  *dimsf77,
       intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77);
void sdf_replace_id_f77_(char *fnamef77, intf * replace, char *labf77,
       char *dtf77, intf  *nbpwf77, intf  *ndimf77, intf8  *dimsf77,
       intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77);
void sdf_replace_id_f77(char *fnamef77, intf * replace, char *labf77,
       char *dtf77, intf  *nbpwf77, intf  *ndimf77, intf8  *dimsf77,
       intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77);
i4 sdf_replace_id(char *fname, i4 replace, data_id * id_new);
data_id *sdf_read(char *fname, i4 order, void **data);
void sdf_details_f77__(char *fnamef77, intf *iorder, char *labf77, char *dtf77, 
     intf *nbpwf77, intf *ndimf77, intf8 *dimsf77,
     intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77);
void sdf_details_f77_(char *fnamef77, intf *iorder, char *labf77, char *dtf77, 
     intf *nbpwf77, intf *ndimf77, intf8 *dimsf77,
     intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77);
void sdf_details_f77(char *fnamef77, intf *iorder, char *labf77, char *dtf77, 
     intf *nbpwf77, intf *ndimf77, intf8 *dimsf77,
     intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77);
data_id *sdf_details(char *fname, i4 order);
i4 sdf_query(char *fname);
i4 sdf_ndat(char *fname);
void sdf_labmatch_f77__(char *fname_f77, intf * ndat,
     char *labtest_f77, intf *match, intf_clen len_fnamef77, intf_clen 
     len_labtestf77);
void sdf_labmatch_f77_(char *fname_f77, intf * ndat,
     char *labtest_f77, intf *match, intf_clen len_fnamef77, 
     intf_clen len_labtestf77);
void sdf_labmatch_f77(char *fname_f77, intf * ndat,
     char *labtest_f77, intf *match, intf_clen len_fnamef77, 
     intf_clen len_labtestf77);
i4 * sdf_labmatch(char *fname, i4 *ndat, char * labtest);
i4 sdf_sizes(char *fname, i4 *ndat, char ** datatypes, pos ** datasizes,
	i4 **nbpw);
void sdf_sizes_f77__(char *fnamef77, intf *ndatf77, char dtypesf77[][1],
   intf8 * datasizesf77, intf *nbpwf77, intf_clen len_fnamef77, 
   intf_clen len_dtypesf77);
void sdf_sizes_f77_(char *fnamef77, intf *ndatf77, char dtypesf77[][1],
   intf8 * datasizesf77, intf *nbpwf77, intf_clen len_fnamef77, 
   intf_clen len_dtypesf77);
void sdf_sizes_f77(char *fnamef77, intf *ndatf77, char dtypesf77[][1],
   intf8 * datasizesf77, intf *nbpwf77, intf_clen len_fnamef77, 
   intf_clen len_dtypesf77);
void  sdf_query_f77__(char *fnamef77, intf * ndat, intf_clen len_fname_f77);
void  sdf_query_f77_(char *fnamef77, intf * ndat, intf_clen len_fname_f77);
void  sdf_query_f77(char *fnamef77, intf * ndat, intf_clen len_fname_f77);
void  sdf_ndat_f77__(char *fnamef77, intf * ndat, intf_clen len_fname_f77);
void  sdf_ndat_f77_(char *fnamef77, intf * ndat, intf_clen len_fname_f77);
void  sdf_ndat_f77(char *fnamef77, intf * ndat, intf_clen len_fname_f77);
pos data_size(data_id * id);
i4 cs2fs(char *cs, char *fs, intf_clen len_fs);
i4 fs2cs(char **cs, char *fs, intf_clen len_fs);
void sdf_rm_f77__(char *fnamef77, intf_clen len_fnamef77);
void sdf_rm_f77_(char *fnamef77, intf_clen len_fnamef77);
void sdf_rm_f77(char *fnamef77, intf_clen len_fnamef77);
i4 sdf_rm(char *fname);
int output_int64 (char * outs, pos l1);
pos atopos(char *s);
void rm_blanks(char *s);
i4 move_dn(FILE *fp, pos shiftdn, pos startloc, pos endloc);
i4 move_up(FILE *fp, pos shiftup, pos startloc, pos endloc);
i4 sdf_delete(char *fname, i4 idelete);
i4 sdf_insert(char *fname, i4 insert, data_id * id_new, void * data_new);
i4 sdf_replace(char *fname, i4 replace, data_id * id_new, void * data_new);
void sdf_housekeeping(char *fname, i4 *norder, pos *hdrposout,
                pos *dataposout, pos *hdrsizeout);
void test_sizes();
i4 sdf_rm2cm(data_id *id, void *data);
i4 sdf_transpose(i4 * indorder, i4 * directions, data_id *id, void *data);
void sdf_transpose_f77__(intf *indorderf77, intf *directionsf77, char * dtf77,
     intf *nbpwf77, intf *ndimf77, intf8* dimsf77, void *dataf77, 
     intf_clen len_dtf77);
void sdf_transpose_f77_(intf *indorderf77, intf *directionsf77, char * dtf77,
     intf *nbpwf77, intf *ndimf77, intf8* dimsf77, void *dataf77, 
     intf_clen len_dtf77);
void sdf_transpose_f77(intf *indorderf77, intf *directionsf77, char * dtf77,
     intf *nbpwf77, intf *ndimf77, intf8* dimsf77, void *dataf77, 
     intf_clen len_dtf77);
i4 readmask(u4 * mask, pos sizem, pos loc);
i4 writemask(u4 *mask, pos sizem, pos loc, i4 value);
i4 rindcalc(i4 ndim, pos *dims, pos *indices, pos memloc_r);
void rindcalc_f77(intf *ndimf77, intf8 *dimsf77, intf8 *indf77, 
      intf8 * memloc_rf77);
void rindcalc_f77_(intf *ndimf77, intf8 *dimsf77, intf8 *indf77, 
      intf8 * memloc_rf77);
void rindcalc_f77__(intf *ndimf77, intf8 *dimsf77, intf8 *indf77, 
      intf8 * memloc_rf77);
i4 cindcalc(i4 ndim, pos *dims, pos *indices, pos memloc_c);
void cindcalc_f77(intf *ndimf77, intf8 *dimsf77, intf8 *indf77,
      intf8 * memloc_cf77);
void cindcalc_f77_(intf *ndimf77, intf8 *dimsf77, intf8 *indf77,
      intf8 * memloc_cf77);
void cindcalc_f77__(intf *ndimf77, intf8 *dimsf77, intf8 *indf77,
      intf8 * memloc_cf77);
i4 memcalcr(pos *dims, pos *indices, i4 ndim, pos *memloc_r);
void memcalcr_f77(intf8 *dimsf77, intf8 *indf77, intf *ndimf77, 
      intf8 *memloc_rf77);
void memcalcr_f77_(intf8 *dimsf77, intf8 *indf77, intf *ndimf77, 
      intf8 *memloc_rf77);
void memcalcr_f77__(intf8 *dimsf77, intf8 *indf77, intf *ndimf77, 
      intf8 *memloc_rf77);
i4 memcalcc(pos *dims, pos *indices, i4 ndim, pos *memloc_c);
void memcalcc_f77(intf8 *dimsf77, intf8 *indf77, intf *ndimf77,
      intf8 *memloc_cf77);
void memcalcc_f77_(intf8 *dimsf77, intf8 *indf77, intf *ndimf77,
      intf8 *memloc_cf77);
void memcalcc_f77__(intf8 *dimsf77, intf8 *indf77, intf *ndimf77,
      intf8 *memloc_cf77);
i4 byteswap (unsigned char *arr, pos arrsize, i4 nbpw);
i4 is_big_endian (void);
void ***** sdf_mk_5d(data_id *id);
void ***** sdf_1d_to_5d(data_id *id, void *data);
void **** sdf_mk_4d(data_id *id);
void **** sdf_1d_to_4d(data_id *id, void *data);
void *** sdf_1d_to_3d(data_id *id, void *data);
void *** sdf_mk_3d(data_id *id);
void ** sdf_1d_to_2d(data_id *id, void *data);
void ** sdf_mk_2d(data_id *id);
void sdf_free_5d(data_id *id, void ***** data5d);
void sdf_free_4d(data_id *id, void **** data4d);
void sdf_free_3d(data_id *id, void *** data3d);
void sdf_free_2d(data_id *id, void ** data2d);
data_id *sdf_cp_id(data_id *id_0);
data_id *sdf_create_id(i4 iorder, char *label, char dtype, i4 nbpw, 
	i4 ndim, pos *dims);
void sdf_free_id(data_id *id);
void *sdf_malloc(pos nmemb);
void sdf_free(void *ptr);
int i4cmp(const void *vp, const void *vq);
i4 sdf_wb(char *fname, pos *datapos, pos nelem, i4 nbpw, FILE *fp, void *data);
i4 sdf_rb(char *fname, pos *datapos, pos nelem, i4 nbpw, FILE *fp, void *data);
void sdf_wb_f77__(char *fnamef77, intf8 *dataposf77, intf8 *nelemf77, 
	intf *nbpwf77, FILE **fp, void *dataf77,  intf_clen len_fnamef77, 
        intf_clen len_dataf77 );
void sdf_wb_f77_(char *fnamef77, intf8 *dataposf77, intf8 *nelemf77, 
	intf *nbpwf77, FILE **FP, void *dataf77,  intf_clen len_fnamef77, 
        intf_clen len_dataf77 );
void sdf_wb_f77(char *fnamef77, intf8 *dataposf77, intf8 *nelemf77, 
	intf *nbpwf77, FILE **fp, void *dataf77,  intf_clen len_fnamef77, 
        intf_clen len_dataf77 );
void sdf_rb_f77__(char *fnamef77, intf8 *dataposf77, intf8 *nelemf77, 
	intf *nbpwf77, FILE **fp, void *dataf77,  intf_clen len_fnamef77, 
        intf_clen len_dataf77 );
void sdf_rb_f77_(char *fnamef77, intf8 *dataposf77, intf8 *nelemf77, 
	intf *nbpwf77, FILE **fp, void *dataf77,  intf_clen len_fnamef77, 
        intf_clen len_dataf77 );
void sdf_rb_f77(char *fnamef77, intf8 *dataposf77, intf8 *nelemf77, 
	intf *nbpwf77, FILE **fp, void *dataf77,  intf_clen len_fnamef77, 
        intf_clen len_dataf77 );
void fread_f77__(void *data, intf * nbpwf77, intf8 *nelemf77, FILE **fp);
void fread_f77_(void *data, intf * nbpwf77, intf8 *nelemf77, FILE **fp);
void fread_f77(void *data, intf * nbpwf77, intf8 *nelemf77, FILE **fp);
void fwrite_f77__(void *data, intf * nbpwf77, intf8 *nelemf77, FILE **fp);
void fwrite_f77_(void *data, intf * nbpwf77, intf8 *nelemf77, FILE **fp);
void fwrite_f77(void *data, intf * nbpwf77, intf8 *nelemf77, FILE **fp);
void fopen_f77__(char *fnamef77, char * modef77, FILE **fp,
               intf_clen len_fnamef77, intf_clen len_modef77);
void fopen_f77_(char *fnamef77, char * modef77, FILE **fp,
               intf_clen len_fnamef77, intf_clen len_modef77);
void fopen_f77(char * fnamef77, char * modef77, FILE **fp,
               intf_clen len_fnamef77, intf_clen len_modef77);
void fclose_f77__(FILE **fp);
void fclose_f77_(FILE **fp);
void fclose_f77(FILE **fp);
void fseeko_f77__(FILE **fp, intf8 *filepos);
void fseeko_f77_(FILE **fp, intf8 *filepos);
void fseeko_f77(FILE **fp, intf8 *filepos);
void ftello_f77__(FILE **fp, intf8 *filepos);
void ftello_f77_(FILE **fp, intf8 *filepos);
void ftello_f77(FILE **fp, intf8 *filepos);
pos file_sz(char *fname);
void file_sz_f77__(char *fnamef77, intf8 *fsize, intf_clen len_fnamef77);
void file_sz_f77_(char *fnamef77, intf8 *fsize, intf_clen len_fnamef77);
void file_sz_f77(char *fnamef77, intf8 *fsize, intf_clen len_fnamef77);
void sdf_ftc_f77__(FILE **fp, intf8 *datasizef77);
void sdf_ftc_f77_(FILE **fp, intf8 *datasizef77);
void sdf_ftc_f77(FILE **fp, intf8 *datasizef77);
void ibe_f77__(intf * ibef77);
void ibe_f77_(intf * ibef77);
void ibe_f77(intf * ibef77);
void byteswap_f77__(void *data, intf8 *nelemf77, intf * nbpwf77);
void byteswap_f77_(void *data, intf8 *nelemf77, intf * nbpwf77);
void byteswap_f77(void *data, intf8 *nelemf77, intf * nbpwf77);
void fflush_f77__(FILE **fp);
void fflush_f77_(FILE **fp);
void fflush_f77(FILE **fp);
void ffstdout_f77__();
void ffstdout_f77_();
void ffstdout_f77();

/* END OF FUNCTION PROTOTYPES */

