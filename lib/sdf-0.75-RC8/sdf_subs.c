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


#include <sdf_subs.h>
void sdf_write_f77__(char *fnamef77, char *labf77, char *dtf77, intf *nbpwf77, 
       intf *ndimf77, intf8 *dimsf77, void *dataf77,  intf_clen len_fnamef77, 
       intf_clen len_labf77, intf_clen len_dtf77, intf_clen len_dataf77 )
{
sdf_write_f77( fnamef77, labf77, dtf77, nbpwf77, 
      ndimf77, dimsf77, dataf77,  len_fnamef77, 
       len_labf77, len_dtf77, len_dataf77 );
return;
}
void sdf_write_f77_(char *fnamef77, char *labf77, char *dtf77, intf *nbpwf77, 
       intf *ndimf77, intf8 *dimsf77, void *dataf77,  intf_clen len_fnamef77, 
       intf_clen len_labf77, intf_clen len_dtf77, intf_clen len_dataf77 )
{
sdf_write_f77( fnamef77, labf77, dtf77, nbpwf77, 
      ndimf77, dimsf77, dataf77,  len_fnamef77, 
       len_labf77, len_dtf77, len_dataf77 );
return;
}
void sdf_write_f77(char *fnamef77, char *labf77, char *dtf77, intf *nbpwf77, 
       intf *ndimf77, intf8 *dimsf77, void *dataf77,  intf_clen len_fnamef77, 
       intf_clen len_labf77, intf_clen len_dtf77, intf_clen len_dataf77 )
{


/*
NOTE -- Fortran character
strings are not passed as a simple pointer, as one might naively guess.  
There's a 2nd hidden integer argument related to the character string length.
Think maybe I got this problem solved now....
*/


data_id *id;
i4 nlab,i;
char *fname;
char *datatype;
char *label;


/* 1st, print out the sizes of f77 char arrays (Debug) */

/*
printf("sizeof(len_dtf77) = %d\n",(int) sizeof(len_dtf77));
printf("len_fnamef77 = %d\n",(int) len_fnamef77);
printf("len_labf77 = %d\n",(int) len_labf77);
printf("len_dtf77 = %d\n",(int) len_dtf77);
printf("len_dataf77 = %d\n",(int) len_dataf77);
*/

/* Convert fortran strings fnamef77,labf77,dtf77 to C-strings fname,label,
   datatype */

fs2cs(&fname,fnamef77,len_fnamef77);
fs2cs(&label,labf77,len_labf77);
fs2cs(&datatype,dtf77,len_dtf77);

/* Debug - print out C-strings fname,label and character variable datatype */

/*
printf("fname = %s\n",fname);
printf("label = %s\n",label);
printf("datatype = %c\n",*datatype);
*/

id=malloc(sizeof(data_id));
nlab=strlen(label);
id->label=calloc(nlab+1,1);


for (i=0;i<nlab;i++)
{
        *(id->label+i)=*(label+i);
}


id->nbpw=(i4) *nbpwf77;
id->datatype=*datatype;
id->order=0; /* value shouldn't matter */
id->ndim=(i4) *ndimf77;
id->dims=malloc(id->ndim * sizeof(pos));


/* I wonder if the f77 integer values of the dimensions shouldn't be
   assumed to be of regular integer type rather than (pos).  
   They could be recast as (pos) inside this fn.  OK, now implemented. */


for (i=0;i<id->ndim;i++)
{
        *(id->dims+i)= (pos) *(dimsf77+id->ndim-1-i); /* flip dims for f77 */
}
sdf_write(fname, id, dataf77);
free(label);
free(fname);
free(datatype);
free(id->label);
free(id->dims);
free(id);
}

void sdf_fl2ds_f77__(char * sdf_fnamef77, char *fnamef77, char *labelf77,
      intf_clen len_sdf_fnamef77, intf_clen len_fnamef77, intf_clen
      len_labelf77)
{
sdf_fl2ds_f77(sdf_fnamef77,fnamef77,labelf77,len_sdf_fnamef77,len_fnamef77,
       len_labelf77);
return;
}

void sdf_fl2ds_f77_(char * sdf_fnamef77, char *fnamef77, char *labelf77,
      intf_clen len_sdf_fnamef77, intf_clen len_fnamef77, intf_clen
      len_labelf77)
{
sdf_fl2ds_f77(sdf_fnamef77,fnamef77,labelf77,len_sdf_fnamef77,len_fnamef77,
       len_labelf77);
return;
}

void sdf_fl2ds_f77(char * sdf_fnamef77, char *fnamef77, char *labelf77,
      intf_clen len_sdf_fnamef77, intf_clen len_fnamef77, intf_clen
      len_labelf77)
{

/*  Provide a fortran-callable version of sdf_fl2ds, which writes an existing
     file into an SDF file as a dataset */

char *sdf_fname;
char *fname;
char *label;

/* convert fortran string arguments to C-string arguments: */

/* debug
printf("len_sdf_fnamef77 = %d\n",len_sdf_fnamef77);
printf("len_fnamef77 = %d\n",len_fnamef77);
printf("len_labelf77 = %d\n",len_labelf77);
*/

fs2cs(&sdf_fname,sdf_fnamef77,len_sdf_fnamef77);
fs2cs(&fname,fnamef77,len_fnamef77);
fs2cs(&label,labelf77,len_labelf77);

/* Call the C-version sdf_fl2ds: */

sdf_fl2ds(sdf_fname,fname,label);

/* free the C-string arguments */

free(sdf_fname);
free(fname);
free(label);

return;
}

i4 sdf_fl2ds(char * sdf_fname, char *fname, char *label)
{

/* Write the file fname into the sdf file sdf_fname as a dataset.  So that
the conversion can occur for files larger than 2GB on a 32-bit machine, 
will do the file transfer in pieces. */

FILE *fp,*fn;
data_id * id=NULL;
pos hdrpos, datapos, hdrtmp, datatmp, datasize, hdrsize;
pos hdrsztmp, oldhdrsize, startloc,dims0,bigbuf,nbuf,fag_end,bufsize;
char *sdfstr="SDF format";
char *filler,*labeltmp;
char testid[300];
char datatype;
void *data;
i4 i,ise,ibe, next_order, next_tmp,nbpw,ndim;
pos hdrinit=HINITSZ, safety=100;
int fctest, fotest,fntest;
/* fotest is a book-keeping variable that is 1 if file is open, 0 otherwise */
fotest=0;


/* Make certain (pos) is defined as 64bit int, otherwise, exit */
test_sizes(); /* NEW WAY */

ibe = is_big_endian ();
ise = 0;
if (ibe == 0) ise = 1;    /* set flag for byteswapping if small endian */

/* Open up file to be read in as read-only, make sure it opens OK */

fn=fopen(fname,"rb");
if(fn == NULL)
{
    printf("sdf_fl2ds: input file %s cannot be read, exiting\n",fname);
    printf("sdf_fl2ds: errno = %d\n",errno);
    fflush(stdout);
    fotest=0;
    exit(1);
}

/* Initially, open the file as read-only */


fp=fopen(sdf_fname,"rb");


if (fp == NULL) /* does file sdf_fname exist? */
{
        /* if file doesnt exist, must create it and initialize it */


/*
        printf("sdf_fl2ds:  File %s does not exist; open and initialize\n",
              sdf_fname);
*/
        fp=fopen(sdf_fname,"wb");
        fotest=1;
        /* write the identifying string 1st thing in the file */
        fwrite(sdfstr,sizeof(char),strlen(sdfstr)+1,fp);
        /* now calculate where initial header info for data will be,
         * after leaving enough space to write in the positions of
         * the header and data parts of the file */
        hdrpos = (pos) ftello(fp)+(pos) 2* sizeof(pos)+(pos) sizeof(i4)
                 +(pos) sizeof(pos);
        hdrtmp=hdrpos;
        /* byteswap the hdr pos if necessary, and write it to file */
        if(ise) byteswap((void *)&hdrtmp,(pos)1,sizeof(pos));
        fwrite(&hdrtmp, sizeof(pos), 1, fp);
        /* now write out 2KB of zeroes, after which the data part will start */
        hdrsize=hdrinit;
        filler=calloc(hdrsize,sizeof(char));
        fwrite(filler,sizeof(char),hdrsize,fp);
        free(filler);
        /* now write the position of the data part to the file */
        datapos = (pos) ftello(fp);
        datatmp=datapos;
        if(ise) byteswap((void *)&datatmp,(pos)1,sizeof(pos));
        fseeko(fp, hdrpos-(pos) sizeof(pos)-(pos)(sizeof(i4))-
                (pos)(sizeof(pos)) ,SEEK_SET); 
        fwrite(&datatmp, sizeof(pos), 1, fp);
        /* now write the order number for next data to be written (0) */
        next_order=(i4) 0;
        next_tmp=next_order;
        /* byte swapping really not necessary since it's 0, but too confusing
         * to not do it */
        if(ise) byteswap((void *)&next_tmp,(pos)1,sizeof(i4));
        fwrite(&next_tmp, sizeof(i4),1,fp);
        hdrsztmp=hdrsize;
        if(ise) byteswap((void *)&hdrsztmp,(pos)1,sizeof(pos));
        fwrite(&hdrsztmp, sizeof(pos),1,fp);
        /* file is now initialized; close it */
        fclose(fp);
        fotest=0;
}
else
{
        fotest=1;
        fgets(testid,strlen(sdfstr)+1,fp);
        /* Now, test to make sure the SDF file identifier string is there */
        if(strncmp(sdfstr,testid,strlen(sdfstr)))
        {
                fclose(fp);
                fotest=0;
                printf("sdf_fl2ds: File %s is not an SDF file, exiting\n",
                        sdf_fname);
                fflush(stdout);
                exit(1);
        }
}


if (fotest) 
{
        fclose(fp); 
        fotest=0;
}


/* re-open file with both read/write permissions */


fp=fopen(sdf_fname,"rb+");
fotest=1;


if(fp == NULL)
{
        fotest=0;
        printf("sdf_fl2ds:  Can't open SDF file %s\n",sdf_fname);
        printf("sdf_fl2ds: errno = %d\n",errno);
        fflush(stdout);
        exit(1);
}



/* Now check to see if label is a NULL string or not.  If so, id->label will
   be set equal to fname */

if(label[0] == '\0')
{
    labeltmp=fname;
}
else
{
    labeltmp=label;
}


/* Now, read the locations for next writes of header info and data */

fseeko(fp, (pos) strlen(sdfstr)+1, SEEK_SET);
fread(&hdrpos,sizeof(pos),1,fp);
if(ise) byteswap((void *)&hdrpos,(pos)1,sizeof(pos));


fread(&datapos,sizeof(pos),1,fp);
if(ise) byteswap((void *)&datapos,(pos)1,sizeof(pos));


fread(&next_order,sizeof(i4),1,fp);
if(ise) byteswap((void *)&next_order,(pos)1,sizeof(i4));


fread(&hdrsize,sizeof(pos),1,fp);
if(ise) byteswap((void *)&hdrsize,(pos)1,sizeof(pos));

/* Now set all of the members of the id structure to 
   correspond to file fname: */

datatype='b';
ndim=(i4)1;
nbpw=(i4)1;
dims0=file_sz(fname);
id=sdf_create_id(next_order,labeltmp,datatype,nbpw,ndim,&dims0);

/* Check to see if header almost full; if so, print warning and increment,
   and record new value at its proper location in the file */


if( hdrpos >= (pos) (hdrsize - safety) )
{
/* now increase size of header area, and move the rest of the file up,
   adjust the current value of datapos accordingly */


 /* For now, comment out this next print statement
 printf("sdf_fl2ds: header area almost full, incrementing by hdrinit.\n");  
 fflush(stdout);
 */

 oldhdrsize=hdrsize;
 hdrsize+=hdrinit;
 hdrsztmp=hdrsize;
 fseeko(fp,(pos)(strlen(sdfstr)+1+2*sizeof(pos)+sizeof(i4)),SEEK_SET);
 if(ise) byteswap((void *)&hdrsztmp,(pos)1,sizeof(pos));
 fwrite(&hdrsztmp, sizeof(pos),1,fp);
 startloc=(pos) strlen(sdfstr)+1+sizeof(pos)+oldhdrsize;
 move_up(fp,hdrinit,startloc,datapos);
 fflush(fp);
 datapos += hdrinit;
 sdf_file_truncate(fp, (pos) datapos);
}

/* check that id->datatype is printable, otherwise in big trouble and must
 * exit.  Actually is really only supposed to be f,i,c, or b acc. to sdf
 * manifesto, but we'll let folks slide here and use anything as long
 * as it's not blank. */


if (!isgraph(id->datatype))
{
        printf("sdf_fl2ds: error - id->datatype = 0x%02x not printable\n",
                id->datatype);
        fflush(stdout);
        fclose(fp);
        fclose(fn);
        fotest=0;
        exit(1);
}


/* Now make sure that id->label is 'blank-free': */


rm_blanks(id->label);


/* Now, output the header information for the data to the file header area */


fseeko(fp,hdrpos,SEEK_SET);
data_to_string(id,fp);
hdrpos=ftello(fp);

datasize=data_size(id);  

/* Set up buffer sizes, etc. here */

bigbuf=1048576; /* 1MB buffersize.  Could be changed */
nbuf=datasize/bigbuf;
fag_end=datasize % bigbuf;
bufsize = (nbuf > 0) ? bigbuf : fag_end;
data=(char *)malloc(bufsize);

/* position file pointers: */

fseeko(fp,datapos,SEEK_SET);
fseeko(fn,(pos)0,SEEK_SET);

/* copy the data from fname to sdf_fname, one buffer-size at a time. 
   It is done this way to avoid limitations of 2GB memory on 32-bit arch. */

for (i=0;i<=nbuf;i++)
{
   bufsize = (nbuf > i) ? bigbuf : fag_end;
   fread(data,bufsize,1,fn);
   fwrite(data,bufsize,1,fp);
}

free(data); /* free the buffer now that file has been copied */

datapos=(pos) ftello(fp);

/* Finish up: increment next_order, and update hdrpos and datapos in the file */

next_order++;
hdrtmp=hdrpos;
datatmp=datapos;
next_tmp=next_order;
hdrsztmp=hdrsize;


if(ise) byteswap((void *)&hdrtmp,(pos)1,sizeof(pos));
if(ise) byteswap((void *)&datatmp,(pos)1,sizeof(pos));
if(ise) byteswap((void *)&next_tmp,(pos)1,sizeof(i4));
if(ise) byteswap((void *)&hdrsztmp,(pos)1,sizeof(pos));


fseeko(fp, (pos) strlen(sdfstr)+1, SEEK_SET);


fwrite(&hdrtmp,sizeof(pos),1,fp);
fwrite(&datatmp,sizeof(pos),1,fp);
fwrite(&next_tmp,sizeof(i4),1,fp);
fwrite(&hdrsztmp,sizeof(pos),1,fp);


fflush(fp);
if(fotest)
{
        fctest=fclose(fp);
        fotest=0;
}
if(fctest != 0)
{
        printf("sdf_fl2ds:  failure to close file %s\n",sdf_fname);
        printf("sdf_fl2ds: errno = %d\n",errno);
        fotest=1;
}
 
fntest=fclose(fn);
if(fntest != 0)
{
        printf("sdf_fl2ds:  failure to close input file %s\n",fname);
        printf("sdf_fl2ds: errno = %d\n",errno);
}

/* We're done */


return 0;
}

void sdf_ds2fl_f77__(char *sdf_fnamef77, intf *orderf77, char *fnamef77,
      intf_clen len_sdf_fnamef77, intf_clen len_fnamef77)
{
sdf_ds2fl_f77(sdf_fnamef77,orderf77,fnamef77,len_sdf_fnamef77,len_fnamef77);
return;
}
void sdf_ds2fl_f77_(char *sdf_fnamef77, intf *orderf77, char *fnamef77,
      intf_clen len_sdf_fnamef77, intf_clen len_fnamef77)
{
sdf_ds2fl_f77(sdf_fnamef77,orderf77,fnamef77,len_sdf_fnamef77,len_fnamef77);
return;
}
void sdf_ds2fl_f77(char *sdf_fnamef77, intf *orderf77, char *fnamef77,
      intf_clen len_sdf_fnamef77, intf_clen len_fnamef77)
{
/* The fortran callable version of sdf_ds2fl, ie sdf_ds2fl_f77, which
   copies datasets in an sdf file to a new file, fnamef77.  If fnamef77
   is null or all blank, the name of the created file is id->label */

char *sdf_fname;
char *fname;

/* convert fortran string arguments to C-string arguments: */

/* debug
printf("len_sdf_fnamef77 = %d\n",len_sdf_fnamef77);
printf("len_fnamef77 = %d\n",len_fnamef77);
*/

fs2cs(&sdf_fname,sdf_fnamef77,len_sdf_fnamef77);
fs2cs(&fname,fnamef77,len_fnamef77);
/* Call the C-version sdf_ds2fl: */

sdf_ds2fl(sdf_fname,(i4) *orderf77,fname);

/* free the C-string arguments */

free(sdf_fname);
free(fname);

/* we're done */
return;
}

i4 sdf_ds2fl(char *sdf_fname, i4 order, char *fname)
{


/* read in the dataset corresponding to order, and write it as file fname. 
   If fname is null, will use id->label from the order'th dataset 
   as assumed filename. This function is dangerous, as it can easily 
   over-write existing files if you are not careful! */


FILE *fp;
FILE *fn;
data_id *id=NULL;
void *data=NULL;
pos datasize=(pos)0, curpos;
char *sdfstr="SDF format";
char temp[301], testid[300];
char *fnametmp;
i4 ise,ibe, iorder, ndatasets,i;
pos hdrsize;
pos bigbuf,fag_end,bufsize,nbuf;


/* Make certain (pos) is defined as 64bit int, otherwise, exit */


test_sizes(); /* New way */


ibe = is_big_endian ();
ise = 0;
if (ibe == 0) ise = 1;    /* set flag for byteswapping if small endian */


/* Open the sdf file as read-only */


fp=fopen(sdf_fname,"rb");

if (fp == NULL) /* does file sdf_fname exist? */
{
        printf("sdf_ds2fl: file %s does not exist, exiting\n",sdf_fname);
        fflush(stdout);
        exit(1);
}
else
{
        fgets(testid,strlen(sdfstr)+1,fp);
        /* Now, test to make sure the SDF file identifier string is there */
        if(strncmp(sdfstr,testid,strlen(sdfstr)))
        {
                fclose(fp);
                printf("sdf_ds2fl: File %s is not an SDF file, exiting\n",
                        sdf_fname);
                fflush(stdout);
                exit(1);
        }
}

/* read the number of datasets in the file: */


fseeko(fp, strlen(sdfstr)+1 + 2*sizeof(pos) , SEEK_SET);
fread(&ndatasets,sizeof(i4),1,fp);
if(ise) byteswap((void *)&ndatasets,(pos)1,sizeof(i4));
/*
 Next print statement possibly important but nuke for now, as it's not an error
*/

/*
printf("sdf_ds2fl: existing datasets = %d\n",ndatasets);
fflush(stdout);
*/

/*
printf("sdf_ds2fl: desired dataset order = %d\n",order);
*/

if(order > (ndatasets-1))
{
        printf("sdf_ds2fl:  Desired dataset order %d > ndatasets-1 = %d\n",
             order, ndatasets-1);
        fflush(stdout);
        exit(1);
}
/* read in header size */
fread(&hdrsize,sizeof(pos),1,fp);
if(ise) byteswap((void *)&hdrsize,(pos)1,sizeof(pos));
/* debug */

/*
output_int64(temp,hdrsize);
printf("sdf_ds2fl: hdrsize = %s\n",temp);
fflush(stdout);
*/


/* read the string data of all the datasets, print it out */
curpos = (pos) strlen(sdfstr)+1+sizeof(pos)+hdrsize;
/*
        printf("DEBUG - sdf_ds2fl: curpos = %qu\n",curpos);
*/
for (iorder=0; iorder <= order; iorder++)
{
/*
        fgets(temp,300,fp);
        printf("sdf_ds2fl: iorder = %d , header string = %s",iorder,temp);
*/


        id=string_to_data(fp);
        datasize=data_size(id);



        if(iorder < order)
        {
                if(id->datatype != 'c') /* test for complex var */
                {
                        curpos += (pos) datasize*(id->nbpw);
                }
                else
                {
                        curpos += (pos) datasize*2*(id->nbpw);
                        /* kludge for complex variables */
                }
                free(id->label);
                free(id->dims);
                free(id);
        }
        else
        {
        output_int64(temp,datasize);
        /* next prints possibly important, but nuke for now, not error */
        /*
        printf("sdf_ds2fl: dataset %d, size %selements, descriptor: ",
               iorder,temp);
        fflush(stdout);
        */


        /*
        data_to_string(id,stdout);
        printf("\n");
        fflush(stdout);
        */

        }
}


/* if fname is NULL, use id->label as the filename */

if(fname[0] == '\0')
{
   fnametmp=id->label;
}
else
{
   fnametmp=fname;
}

/* Open the file for which the dataset will be written */
fn=fopen(fnametmp,"wb");

if (fn == NULL) /* Problem opening file for writing? */
{
        printf("sdf_ds2fl: output file %s will not open, exiting\n",fnametmp);
        fflush(stdout);
        exit(1);
}

if (id->datatype != 'c') /* test for complex var, mpy no. bytes by 2 if so */
{
        datasize*=(pos)(id->nbpw);
}
else
{
        datasize*=(pos)2 * (pos)(id->nbpw);
}

/* Set up buffer sizes, etc. here */

bigbuf=1048576; /* 1MB buffersize.  Could be changed */
nbuf=datasize/bigbuf;
fag_end=datasize % bigbuf;
bufsize = (nbuf > 0) ? bigbuf : fag_end;
data=(char *)malloc(bufsize);

/* position file pointers: */

fseeko(fp,curpos,SEEK_SET);
fseeko(fn,(pos)0,SEEK_SET);

/* copy the data from sdf_fname to fnametmp, one buffer-size at a time.
   It is done this way to avoid limitations of 2GB memory on 32-bit arch. */

for (i=0;i<=nbuf;i++)
{
   bufsize = (nbuf > i) ? bigbuf : fag_end;
   fread(data,bufsize,1,fp);
   fwrite(data,bufsize,1,fn);
}

free(data); /* free the buffer now that file has been copied */
fclose(fn); /* close output file */

sdf_free_id(id);
fclose(fp); /* close sdf file */
return 0;
}

void sdf_read_f77__(char *fnamef77, intf *iorder_f77, char *labf77, 
     char *dtf77, intf *nbpwf77, intf *ndimf77, intf8 *dimsf77, void *dataf77, 
     intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77, 
     intf_clen len_dataf77)
{
sdf_read_f77(fnamef77, iorder_f77, labf77, dtf77, 
     nbpwf77, ndimf77, dimsf77, dataf77, 
     len_fnamef77, len_labf77, len_dtf77, len_dataf77);
return;
}


void sdf_read_f77_(char *fnamef77, intf *iorder_f77, char *labf77, char *dtf77, 
     intf *nbpwf77, intf *ndimf77, intf8 *dimsf77, void *dataf77, 
     intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77, 
     intf_clen len_dataf77)
{
sdf_read_f77(fnamef77, iorder_f77, labf77, dtf77, 
     nbpwf77, ndimf77, dimsf77, dataf77, 
     len_fnamef77, len_labf77, len_dtf77, len_dataf77);
return;
}


void sdf_read_f77(char *fnamef77, intf *iorder_f77, char *labf77, char *dtf77, 
     intf *nbpwf77, intf *ndimf77, intf8 *dimsf77, void *dataf77, 
     intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77, 
     intf_clen len_dataf77)


{
data_id *id=NULL;
char *fname;
FILE *fp;
pos datasize=(pos)0, curpos;
pos frtest;
char *sdfstr="SDF format";
char temp[301], testid[300];
i4 i, ise, ibe, iorder, order, ndatasets;
pos hdrsize;


/* Make certain (pos) is defined as 64bit int, otherwise, exit */
test_sizes(); 


/* DEBUG - print out lengths of fortran stings (hidden args) */
/*
printf("len_fnamef77 = %d\n",len_fnamef77);
printf("len_labf77 = %d\n",len_labf77);
printf("len_dtf77 = %d\n",len_dtf77);
printf("len_dataf77 = %d\n",len_dataf77);
*/


/* Convert fnamef77 to C-string fname */
fs2cs(&fname,fnamef77,len_fnamef77);
order = (i4) *iorder_f77;
/* DEBUG - print out C-string version of filename fname */
/*
printf("fname = %s\n",fname);
*/


ibe = is_big_endian ();
ise = 0;
if (ibe == 0) ise = 1;    /* set flag for byteswapping if small endian */


/* Initially, open the file as read-only */


fp=fopen(fname,"rb");


if (fp == NULL) /* does file fname exist? */
{
        printf("sdf_read_f77: file %s does not exist, exiting\n",fname);
        fflush(stdout);
        exit(1);
}
else
{
        fgets(testid,strlen(sdfstr)+1,fp);
        /* Now, test to make sure the SDF file identifier string is there */
        if(strncmp(sdfstr,testid,strlen(sdfstr)))
        {
                fclose(fp);
                printf("sdf_read_f77: File %s is not an SDF file, exiting\n",
                        fname);
                fflush(stdout);
                exit(1);
        }
}
/* read the number of datasets in the file: */


fseeko(fp, strlen(sdfstr)+1 + 2*sizeof(pos) , SEEK_SET);
fread(&ndatasets,sizeof(i4),1,fp);
if(ise) byteswap((void *)&ndatasets,(pos)1,sizeof(i4));
/*  Next print statement possibly important, but nuke for now */


/*
printf("sdf_read_f77: existing datasets = %d\n",ndatasets);
fflush(stdout);
*/


/*
printf("sdf_read_f77: desired dataset order = %d\n",order);
*/
if(order > (ndatasets-1))
{
        printf("sdf_read_f77:  Desired dataset order %d > ndatasets-1 = %d\n",
             order, ndatasets-1);
        fflush(stdout);
        exit(1);
}
if(order < 0)
{
        printf("sdf_read_f77:  Desired dataset order %d < 0 \n",
             order);
        fflush(stdout);
        exit(1);
}
/* read in header size */
fread(&hdrsize,sizeof(pos),1,fp);
if(ise) byteswap((void *)&hdrsize,(pos)1,sizeof(pos));
/* debug: */

/*
output_int64(temp,hdrsize);
printf("sdf_read_f77: hdrsize = %s\n",temp);
fflush(stdout);
*/


/* read the string data of all the datasets, print it out */
curpos = (pos) strlen(sdfstr)+1+sizeof(pos)+hdrsize;
/*
        printf("DEBUG - sdf_read_f77: curpos = %qu\n",curpos);
*/
for (iorder=0; iorder <= order; iorder++)
{
/*
        fgets(temp,300,fp);
        printf("iorder = %d , header string = %s",iorder,temp);
*/


        id=string_to_data(fp);
        datasize=data_size(id);



        if(iorder < order)
        {
                if(id->datatype != 'c') /* test for complex var */
                {
                        curpos += (pos) datasize*(id->nbpw);
                }
                else
                {
                        curpos += (pos) datasize*2*(id->nbpw);
                        /* kludge for complex variables */
                }
                free(id->label);
                free(id->dims);
                free(id);
        }
        else
        {


/*
        output_int64(temp,datasize);
*/


        /*
         Next print statements possibly important, but nuke for now
        printf("sdf_read_f77: dataset %d, size %selements, descriptor: ",
               iorder,temp);
        fflush(stdout);
        data_to_string(id,stdout);
        printf("\n");
        fflush(stdout);
        */
        }
}


fseeko(fp,curpos,SEEK_SET);                /* reposition file pointer  */


if (id->datatype != 'c') /* test for complex var */
{
        frtest=(pos) fread(dataf77, id->nbpw , datasize,fp); /* read the data */
	if(frtest < datasize)
	{
           /* test for unsuccessful read: */
           printf("sdf_read_f77: read error\n");
           printf("sdf_read_f77: errno = %d\n",errno);
           output_int64(temp,frtest);
           printf("sdf_read_f77: frtest = %s\n",temp);
           output_int64(temp,datasize);
           printf("sdf_read_f77: datasize =%s\n",temp);
           fflush(stdout);
           exit(1);
	}
        if(ise) byteswap(dataf77,datasize,id->nbpw); /* byteswap if neeeded */
}
else
{
        /* kludge for reading complex variable data */
        frtest=(pos) fread(dataf77, id->nbpw,2*datasize,fp);/* read the data */
	if(frtest < (pos)2 * datasize)
	{
           /* test for unsuccessful read: */
           printf("sdf_read_f77: read error\n");
           printf("sdf_read_f77: errno = %d\n",errno);
           output_int64(temp,frtest);
           printf("sdf_read_f77: frtest = %s\n",temp);
           output_int64(temp,(pos)2*datasize);
           printf("sdf_read_f77(complex): 2*datasize =%s\n",temp);
           fflush(stdout);
           exit(1);
	}
        if(ise) byteswap(dataf77,2*datasize,id->nbpw); /* byteswap if neeeded */
}

/* now convert id->label into fortran string labf77 */
cs2fs(id->label,labf77,len_labf77);


*dtf77 = id->datatype;
*nbpwf77 = (intf) id->nbpw;
*ndimf77 = (intf) id->ndim;
for (i=0; i < id->ndim; i++)
{
        /* flip order for f77 */
        *(dimsf77+i) = (intf8) *(id->dims+id->ndim-1 -i); 
}
datasize=data_size(id);


free(fname);
free(id->label);
free(id->dims);
free(id);
fclose(fp);


}


void sdf_delete_f77__(char *fnamef77, intf* idelete, intf_clen len_fnamef77)
{
sdf_delete_f77(fnamef77, idelete, len_fnamef77);
return;
}
void sdf_delete_f77_(char *fnamef77, intf* idelete, intf_clen len_fnamef77)
{
sdf_delete_f77(fnamef77, idelete, len_fnamef77);
return;
}
void sdf_delete_f77(char *fnamef77, intf* idelete, intf_clen len_fnamef77)
{


/*
NOTE -- Fortran character
strings are not passed as a simple pointer, as one might naively guess.  
There's a 2nd hidden integer argument related to the character string length.
Think maybe I got this problem solved now....
*/


char *fname;


/* 1st, print out the sizes of f77 char arrays (debug) */


/*
printf("len_fnamef77 = %d\n",len_fnamef77);
printf("len_labf77 = %d\n",len_labf77);
printf("len_dtf77 = %d\n",len_dtf77);
printf("len_dataf77 = %d\n",len_dataf77);
*/


/* Convert fortran strings fnamef77,labf77,dtf77 to C-strings fname,label,
   datatype */


fs2cs(&fname,fnamef77,len_fnamef77);


/* Debug - print out C-strings fname,label and character variable datatype */


/*
printf("fname = %s\n",fname);
printf("label = %s\n",label);
printf("datatype = %c\n",*datatype);
*/



/* I wonder if the f77 integer values of the dimensions shouldn't be
   assumed to be of regular integer type rather than (pos).  
   They could be recast as (pos) inside this fn.  OK, now implemented. */


sdf_delete(fname, (i4) *idelete);
free(fname);
}


void sdf_insert_f77__(char *fnamef77, intf* insert, char *labf77, 
       char *dtf77, intf *nbpwf77, intf *ndimf77, intf8 *dimsf77, 
       void *dataf77,  intf_clen len_fnamef77, intf_clen len_labf77, 
       intf_clen len_dtf77, intf_clen len_dataf77 )
{
sdf_insert_f77(fnamef77, insert, labf77, 
       dtf77, nbpwf77, ndimf77, dimsf77, 
       dataf77,  len_fnamef77, len_labf77, len_dtf77, 
       len_dataf77 );
return;
}
void sdf_insert_f77_(char *fnamef77, intf* insert, char *labf77, 
       char *dtf77, intf *nbpwf77, intf *ndimf77, intf8 *dimsf77, 
       void *dataf77,  intf_clen len_fnamef77, intf_clen len_labf77, 
       intf_clen len_dtf77, intf_clen len_dataf77 )
{
sdf_insert_f77(fnamef77, insert, labf77, 
       dtf77, nbpwf77, ndimf77, dimsf77, 
       dataf77,  len_fnamef77, len_labf77, len_dtf77, 
       len_dataf77 );
return;
}


void sdf_insert_f77(char *fnamef77, intf* insert, char *labf77, 
       char *dtf77, intf *nbpwf77, intf *ndimf77, intf8 *dimsf77, 
       void *dataf77,  intf_clen len_fnamef77, intf_clen len_labf77, 
       intf_clen len_dtf77, intf_clen len_dataf77 )
{


/*
NOTE -- Fortran character
strings are not passed as a simple pointer, as one might naively guess.  
There's a 2nd hidden integer argument related to the character string length.
Think maybe I got this problem solved now....
*/


data_id *id_new;
i4 nlab,i;
char *fname;
char *datatype;
char *label;


/* 1st, print out the sizes of f77 char arrays (debug) */


/*
printf("len_fnamef77 = %d\n",len_fnamef77);
printf("len_labf77 = %d\n",len_labf77);
printf("len_dtf77 = %d\n",len_dtf77);
printf("len_dataf77 = %d\n",len_dataf77);
*/


/* Convert fortran strings fnamef77,labf77,dtf77 to C-strings fname,label,
   datatype */


fs2cs(&fname,fnamef77,len_fnamef77);
fs2cs(&label,labf77,len_labf77);
fs2cs(&datatype,dtf77,len_dtf77);


/* Debug - print out C-strings fname,label and character variable datatype */


/*
printf("fname = %s\n",fname);
printf("label = %s\n",label);
printf("datatype = %c\n",*datatype);
*/


id_new=malloc(sizeof(data_id));
nlab=strlen(label);
id_new->label=calloc(nlab+1,1);


for (i=0;i<nlab;i++)
{
        *(id_new->label+i)=*(label+i);
}


id_new->nbpw=(i4) *nbpwf77;
id_new->datatype=*datatype;
id_new->order=0; /* value shouldn't matter */
id_new->ndim=(i4) *ndimf77;
id_new->dims=malloc(id_new->ndim * sizeof(pos));


/* I wonder if the f77 integer values of the dimensions shouldn't be
   assumed to be of regular integer type rather than (pos).  
   They could be recast as (pos) inside this fn.  OK, now implemented. */


for (i=0;i<id_new->ndim;i++)
{
        *(id_new->dims+i)= (pos) *(dimsf77+id_new->ndim-1-i); /*flip dims f77*/
}
sdf_insert(fname, (i4) *insert, id_new, dataf77);
free(label);
free(fname);
free(datatype);
free(id_new->label);
free(id_new->dims);
free(id_new);
}


void sdf_replace_f77__(char *fnamef77, intf * replace, char *labf77, 
       char *dtf77, intf  *nbpwf77, intf  *ndimf77, intf8  *dimsf77, 
       void *dataf77,  intf_clen len_fnamef77, intf_clen len_labf77, 
       intf_clen len_dtf77, intf_clen len_dataf77 )
{
sdf_replace_f77(fnamef77, replace, labf77, 
       dtf77, nbpwf77, ndimf77, dimsf77, 
       dataf77, len_fnamef77, len_labf77, len_dtf77, len_dataf77 );
return;
}
void sdf_replace_f77_(char *fnamef77, intf * replace, char *labf77, 
       char *dtf77, intf  *nbpwf77, intf  *ndimf77, intf8  *dimsf77, 
       void *dataf77,  intf_clen len_fnamef77, intf_clen len_labf77, 
       intf_clen len_dtf77, intf_clen len_dataf77 )
{
sdf_replace_f77(fnamef77, replace, labf77, 
       dtf77, nbpwf77, ndimf77, dimsf77, 
       dataf77, len_fnamef77, len_labf77, len_dtf77, len_dataf77 );
return;
}
void sdf_replace_f77(char *fnamef77, intf * replace, char *labf77, 
       char *dtf77, intf  *nbpwf77, intf  *ndimf77, intf8 *dimsf77, 
       void *dataf77,  intf_clen len_fnamef77, intf_clen len_labf77, 
       intf_clen len_dtf77, intf_clen len_dataf77 )
{


/*
NOTE -- Fortran character
strings are not passed as a simple pointer, as one might naively guess.  
There's a 2nd hidden integer argument related to the character string length.
Think maybe I got this problem solved now....
*/


data_id *id_new;
i4 nlab,i;
char *fname;
char *datatype;
char *label;


/* 1st, print out the sizes of f77 char arrays (debug) */


/*
printf("len_fnamef77 = %d\n",len_fnamef77);
printf("len_labf77 = %d\n",len_labf77);
printf("len_dtf77 = %d\n",len_dtf77);
printf("len_dataf77 = %d\n",len_dataf77);
*/


/* Convert fortran strings fnamef77,labf77,dtf77 to C-strings fname,label,
   datatype */


fs2cs(&fname,fnamef77,len_fnamef77);
fs2cs(&label,labf77,len_labf77);
fs2cs(&datatype,dtf77,len_dtf77);


/* Debug - print out C-strings fname,label and character variable datatype */


/*
printf("fname = %s\n",fname);
printf("label = %s\n",label);
printf("datatype = %c\n",*datatype);
*/


id_new=malloc(sizeof(data_id));
nlab=strlen(label);
id_new->label=calloc(nlab+1,1);


for (i=0;i<nlab;i++)
{
        *(id_new->label+i)=*(label+i);
}


id_new->nbpw=(i4) *nbpwf77;
id_new->datatype=*datatype;
id_new->order=0; /* value shouldn't matter */
id_new->ndim= (i4) *ndimf77;
id_new->dims=malloc(id_new->ndim * sizeof(pos));


/* I wonder if the f77 integer values of the dimensions shouldn't be
   assumed to be of regular integer type rather than (pos).  
   They could be recast as (pos) inside this fn.  OK, now implemented. */


for (i=0;i<id_new->ndim;i++)
{
        *(id_new->dims+i)= (pos) *(dimsf77+id_new->ndim-1-i); /*flip dims f77*/
}
sdf_replace(fname, (i4) *replace, id_new, dataf77);
free(label);
free(fname);
free(datatype);
free(id_new->label);
free(id_new->dims);
free(id_new);

return;
}

void sdf_replace_id_f77__(char *fnamef77, intf * replace, char *labf77,
       char *dtf77, intf  *nbpwf77, intf  *ndimf77, intf8  *dimsf77,
       intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77)
{
sdf_replace_id_f77(fnamef77, replace, labf77, dtf77, nbpwf77, ndimf77, dimsf77,
       len_fnamef77, len_labf77, len_dtf77);
return;
}

void sdf_replace_id_f77_(char *fnamef77, intf * replace, char *labf77,
       char *dtf77, intf  *nbpwf77, intf  *ndimf77, intf8  *dimsf77,
       intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77)
{
sdf_replace_id_f77(fnamef77, replace, labf77, dtf77, nbpwf77, ndimf77, dimsf77,
       len_fnamef77, len_labf77, len_dtf77);
return;
}

void sdf_replace_id_f77(char *fnamef77, intf * replace, char *labf77,
       char *dtf77, intf  *nbpwf77, intf  *ndimf77, intf8 *dimsf77,
       intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77)
{

data_id *id_new;
i4 nlab,i;
char *fname;
char *datatype;
char *label;


/* 1st, print out the sizes of f77 char arrays (debug) */

/*
printf("len_fnamef77 = %d\n",len_fnamef77);
printf("len_labf77 = %d\n",len_labf77);
printf("len_dtf77 = %d\n",len_dtf77);
*/


/* Convert fortran strings fnamef77,labf77,dtf77 to C-strings fname,label,
   datatype */


fs2cs(&fname,fnamef77,len_fnamef77);
fs2cs(&label,labf77,len_labf77);
fs2cs(&datatype,dtf77,len_dtf77);


/* Debug - print out C-strings fname,label and character variable datatype */


/*
printf("fname = %s\n",fname);
printf("label = %s\n",label);
printf("datatype = %c\n",*datatype);
*/


id_new=malloc(sizeof(data_id));
nlab=strlen(label);
id_new->label=calloc(nlab+1,1);


for (i=0;i<nlab;i++)
{
        *(id_new->label+i)=*(label+i);
}


id_new->nbpw=(i4) *nbpwf77;
id_new->datatype=*datatype;
id_new->order=0; /* value shouldn't matter */
id_new->ndim= (i4) *ndimf77;
id_new->dims=malloc(id_new->ndim * sizeof(pos));

/* I wonder if the f77 integer values of the dimensions shouldn't be
   assumed to be of regular integer type rather than (pos).
   They could be recast as (pos) inside this fn.  OK, now implemented. */

for (i=0;i<id_new->ndim;i++)
{
        *(id_new->dims+i)= (pos) *(dimsf77+id_new->ndim-1-i); /*flip dims f77*/
}
sdf_replace_id(fname, (i4) *replace, id_new);
free(label);
free(fname);
free(datatype);
sdf_free_id(id_new);
return;
}

i4 sdf_replace_id(char *fname, i4 replace, data_id * id_new)
{


/* Replace id descriptor string in dataset "replace" with that from
descriptor id_new in file fname.  If the no. of bytes corresponding to
the descriptor in id_new is not equal to those of the original id, no
change is made, and an error message is printed. */


FILE *fp;
data_id *id;
pos hdrpos, datapos, hdrtmp, datatmp, curpos, off, offset;
pos filepos;
char *sdfstr="SDF format";
char *buff;
char testid[300];
i4 ise,ibe, next_order, next_tmp, bufsize, iorder;
pos hdrinit=HINITSZ, safety=100;
pos hdrsize, oldhdrsize, hdrsztmp, startloc;
pos datasize, datasize_old, bytesize, bytesize_old;
int fotest;


/* Make certain (pos) is defined as 64bit int, otherwise, exit */


test_sizes();
fotest=0;


ibe = is_big_endian ();
ise = 0;
if (ibe == 0) ise = 1;    /* set flag for byteswapping if small endian */


/* Initially, open the file as read-only */


fp=fopen(fname,"rb");


if (fp == NULL) /* does file fname exist? */
{
        /* if not, print msg and exit */
        printf("sdf_replace_id: file %s does not exist; exiting\n",fname);
        fflush(stdout);
        exit(1);
}
else
{
        fotest=1;
        fgets(testid,strlen(sdfstr)+1,fp);
        /* Now, test to make sure the SDF file identifier string is there */
        /* If not, print msg and exit */
        if(strncmp(sdfstr,testid,strlen(sdfstr)))
        {
                fclose(fp);
                fotest=0;
                printf("sdf_replace_id: File %s is not an SDF file, exiting\n",
                        fname);
                fflush(stdout);
                exit(1);
        }
}


if (fotest)
{
        /*close file if it is open */
        fclose(fp);
        fotest=0;
}


/* re-open file with both read/write permissions */


fp=fopen(fname,"rb+");


/* Now, read the locations for next writes of header info and data */


fseeko(fp, (pos) strlen(sdfstr)+1, SEEK_SET);
fread(&hdrpos,sizeof(pos),1,fp);
if(ise) byteswap((void *)&hdrpos,(pos)1,sizeof(pos));


fread(&datapos,sizeof(pos),1,fp);
if(ise) byteswap((void *)&datapos,(pos)1,sizeof(pos));


fread(&next_order,sizeof(i4),1,fp);
if(ise) byteswap((void *)&next_order,(pos)1,sizeof(i4));


fread(&hdrsize,sizeof(pos),1,fp);
if(ise) byteswap((void *)&hdrsize,(pos)1,sizeof(pos));


if(( replace >= next_order) || (replace < 0))
{
        printf("sdf_replace_id: can't replace id for dataset %d in %s\n",
            replace,fname);
        fflush(stdout);
        fclose(fp);
        exit(1);
}


if( hdrpos >= (pos) (hdrsize - safety) )
{
/* now increase size of header area, and move the rest of the file up,
   adjust the current value of datapos accordingly */
 printf("sdf_replace_id: header area almost full, incrementing by hdrinit.\n");
 fflush(stdout);
 oldhdrsize=hdrsize;
 hdrsize+=hdrinit;
 hdrsztmp=hdrsize;
 fseeko(fp,(pos)(strlen(sdfstr)+1+2*sizeof(pos)+sizeof(i4)),SEEK_SET);
 if(ise) byteswap((void *)&hdrsztmp,(pos)1,sizeof(pos));
 fwrite(&hdrsztmp, sizeof(pos),1,fp);
 startloc=(pos) strlen(sdfstr)+1+sizeof(pos)+oldhdrsize;
 move_up(fp,hdrinit,startloc,datapos);
 fflush(fp);
 datapos += hdrinit;
 sdf_file_truncate(fp, (pos) datapos);
}


bufsize=hdrsize-sizeof(pos)-sizeof(i4); 
buff=(char *)calloc(bufsize*sizeof(char),1);


off=(pos)0;
curpos = (pos) strlen(sdfstr)+1+sizeof(pos)+hdrsize;
/* ensure filepos is at first string descriptor if file resized as above */ 
filepos=strlen(sdfstr)+1 + 2*sizeof(pos) + sizeof(i4) + sizeof(pos);
fseeko(fp,filepos,SEEK_SET);

for (iorder=0; iorder < next_order; iorder++) 
{
        id=(data_id *) string_to_data(fp);
        filepos=ftello(fp);
        if(iorder == replace) 
        {
                datasize_old=data_size(id);
                bytesize_old=datasize_old* (pos)(id->nbpw);
                if(id->datatype == 'c') /* test for complex variable */
                {
                   bytesize_old*=(pos)2;
                }
                id_new->order=replace;
                datasize=data_size(id_new);
                bytesize=datasize* (pos) (id_new->nbpw);
                if(id_new->datatype == 'c') /* test for complex variable */
                {
                   bytesize*=(pos)2;
                }
                /* Next, test to make sure no. of bytes for the data 
                   described by new header is the same as that of old header.  
                   If this is not true, the edited file will be corrupt!  
                   Must then exit without doing anything */
                if(bytesize != bytesize_old)
                {
                   printf("sdf_replace_id: new datasize incompatible w/data\n");
                   printf("sdf_replace_id: no change made, exiting\n");
                   fflush(stdout);
                   fclose(fp);
                   exit(1);
/*
                   Decided to delete the return statement and replace it
                   with the above exit.  If datasize is changed from original
                   value, the file would become corrupted.  This should 
                   be regarded as a fatal user error.  
                   Must force user to fix their code.
*/

                }
                /* create new string for new dataset */
                fseeko(fp,filepos,SEEK_SET); 
                data_to_buff(id_new,buff+off,&offset);
                off+=offset;

        }

        else

        {

                /* do this stuff only if dataset is not being replaced */

                data_to_buff(id,buff+off,&offset);
                off+=offset;

                /* end of stuff to be done if datset is not being replaced */

        }


        sdf_free_id(id);


}
        /* remember to free id_new in calling program */


hdrpos = strlen(sdfstr)+1 + 2*sizeof(pos) + sizeof(i4) + sizeof(pos);
fseeko(fp, hdrpos, SEEK_SET);
/*  Gonna change this to get rid of \0 at end */
/* fprintf(fp,"%s%c",buff,'\0'); */
fprintf(fp,"%s",buff);
hdrpos= (pos) ftello(fp);

hdrtmp=hdrpos;
datatmp=datapos;
next_tmp=next_order;
hdrsztmp=hdrsize;


if(ise) byteswap((void *)&hdrtmp,(pos)1,sizeof(pos));
if(ise) byteswap((void *)&datatmp,(pos)1,sizeof(pos));
if(ise) byteswap((void *)&next_tmp,(pos)1,sizeof(i4));
if(ise) byteswap((void *)&hdrsztmp,(pos)1,sizeof(pos));


fseeko(fp, (pos) strlen(sdfstr)+1, SEEK_SET);


fwrite(&hdrtmp,sizeof(pos),1,fp);
fwrite(&datatmp,sizeof(pos),1,fp);
fwrite(&next_tmp,sizeof(i4),1,fp);
fwrite(&hdrsztmp,sizeof(pos),1,fp);


fclose(fp);
free(buff);


/* We're done */


return 0;
}

i4 sdf_write(char *fname, data_id *id, void *data)
{


/* Write the data id header and the data into file fname */


/* File structure: (assuming the id-string is "SDF format"):


 * 1.  Bytes 0-10 - file identifier string, including null char terminator 
 * 2.  Bytes 11-18 - 64-bit integer, containing file position of next byte
 *     in the header.  Byte 11 is also strlen(sdfstr)+1 .
 * 3.  Bytes 19-26 - 64-bit integer, containing file position of next byte
       in the data portion of the file (which starts at byte 2019).
       Byte 19 is also strlen(sdfstr)+1+sizeof(pos) .
 * 4.  Bytes 27-30 - 32-bit integer, containing the number of datasets 
       currently written in the file.  Byte 27 is also
       strlen(sdfstr)+1 +2*sizeof(pos) .
 * 5.  Bytes 31-38 - 64-bit integer, containing the size of the file header,
       'hdrsize', counting from byte 19.  Initially this is set to HINITSZ, 
       whose default value is 2000, but can be increased in blocks of 
       HINITSZ bytes as needed.  HINITSZ can be set to something other than
       2000 with a compiler flag, e.g. -DHINITSZ=5000, or some other value.
       In the following discussion, a default value of HINITSZ=2000 has been
       assumed.
 * 6.  Bytes 39-2018 - A series of linefeed terminated strings, containing
       descriptions of the data to be written.  These strings generated by
       the data_to_string function.  Byte 39 is also
       strlen(sdfstr)+1 + 2*sizeof(pos) + sizeof(i4) *sizeof(pos) .
 * 7.  Bytes 2019-eof - the data, with each dataset written squentially,
       in order.  Byte 2019 is also equal to
       strlen(sdfstr)+1+sizeof(pos)+hdrsize .  If hdrsize is larger than
       HINITSZ, the data will start correspondingly later in the file
       
       The integers in the header and the data are all written in large-endian
       byte order, if the data isn't simply byte data. */


FILE *fp;
pos hdrpos, datapos, hdrtmp, datatmp, datasize, hdrsize;
pos hdrsztmp, oldhdrsize, startloc;
char *sdfstr="SDF format";
char *filler;
char testid[300];
i4 i,ise,ibe, next_order, next_tmp;
pos hdrinit=HINITSZ, safety=100;
int fctest, fotest;
/* fotest is a book-keeping variable that is 1 if file is open, 0 otherwise */
fotest=0;


/* Make certain (pos) is defined as 64bit int, otherwise, exit */


test_sizes(); /* NEW WAY */


/*  OLD WAY
if(sizeof(pos) != 8)
{
        printf("sdf_write: error building sdf_write: sizeof(pos) = %d not 8\n",
            sizeof(pos));
        fflush(stdout);
        exit(1);
}
*/


ibe = is_big_endian ();
ise = 0;
if (ibe == 0) ise = 1;    /* set flag for byteswapping if small endian */

/* test for 0 or negative id->npbw; exit with error message if so */

if(id->nbpw <= 0)
        {
           printf("sdf_write: id->nbpw = %d <= 0; fatal\n",(int)(id->nbpw));
           fflush(stdout);
           exit(1);
        }


/* Initially, open the file as read-only */


fp=fopen(fname,"rb");


if (fp == NULL) /* does file fname exist? */
{
        /* if file doesnt exist, must create it and initialize it */


/*
        printf("sdf_write:  File %s does not exist; open and initialize\n",
              fname);
*/
        fp=fopen(fname,"wb");
        fotest=1;
        /* write the identifying string 1st thing in the file */
        fwrite(sdfstr,sizeof(char),strlen(sdfstr)+1,fp);
        /* now calculate where initial header info for data will be,
         * after leaving enough space to write in the positions of
         * the header and data parts of the file */
        hdrpos = (pos) ftello(fp)+(pos) 2* sizeof(pos)+(pos) sizeof(i4)
                 +(pos) sizeof(pos);
        hdrtmp=hdrpos;
        /* byteswap the hdr pos if necessary, and write it to file */
        if(ise) byteswap((void *)&hdrtmp,(pos)1,sizeof(pos));
        fwrite(&hdrtmp, sizeof(pos), 1, fp);
        /* now write out 2KB of zeroes, after which the data part will start */
        hdrsize=hdrinit;
        filler=calloc(hdrsize,sizeof(char));
        fwrite(filler,sizeof(char),hdrsize,fp);
        free(filler);
        /* now write the position of the data part to the file */
        datapos = (pos) ftello(fp);
        datatmp=datapos;
        if(ise) byteswap((void *)&datatmp,(pos)1,sizeof(pos));
        fseeko(fp, hdrpos-(pos) sizeof(pos)-(pos)(sizeof(i4))-
                (pos)(sizeof(pos)) ,SEEK_SET); 
        fwrite(&datatmp, sizeof(pos), 1, fp);
        /* now write the order number for next data to be written (0) */
        next_order=(i4) 0;
        next_tmp=next_order;
        /* byte swapping really not necessary since it's 0, but too confusing
         * to not do it */
        if(ise) byteswap((void *)&next_tmp,(pos)1,sizeof(i4));
        fwrite(&next_tmp, sizeof(i4),1,fp);
        hdrsztmp=hdrsize;
        if(ise) byteswap((void *)&hdrsztmp,(pos)1,sizeof(pos));
        fwrite(&hdrsztmp, sizeof(pos),1,fp);
        /* file is now initialized; close it */
        fclose(fp);
        fotest=0;
}
else
{
        fotest=1;
        fgets(testid,strlen(sdfstr)+1,fp);
        /* Now, test to make sure the SDF file identifier string is there */
        if(strncmp(sdfstr,testid,strlen(sdfstr)))
        {
                fclose(fp);
                fotest=0;
                printf("sdf_write: File %s is not an SDF file, exiting\n",
                        fname);
                fflush(stdout);
                exit(1);
        }
}


if (fotest) 
{
        fclose(fp); 
        fotest=0;
}


/* re-open file with both read/write permissions */


fp=fopen(fname,"rb+");
fotest=1;


if(fp == NULL)
{
        fotest=0;
        printf("sdf_write:  Can't open file %s\n",fname);
        printf("sdf_write: errno = %d\n",errno);
        fflush(stdout);
        exit(1);
}


/* Now, read the locations for next writes of header info and data */


fseeko(fp, (pos) strlen(sdfstr)+1, SEEK_SET);
fread(&hdrpos,sizeof(pos),1,fp);
if(ise) byteswap((void *)&hdrpos,(pos)1,sizeof(pos));


fread(&datapos,sizeof(pos),1,fp);
if(ise) byteswap((void *)&datapos,(pos)1,sizeof(pos));


fread(&next_order,sizeof(i4),1,fp);
if(ise) byteswap((void *)&next_order,(pos)1,sizeof(i4));


fread(&hdrsize,sizeof(pos),1,fp);
if(ise) byteswap((void *)&hdrsize,(pos)1,sizeof(pos));


/* make sure order value written out is actually the order of
 * data written to file */


id->order=next_order;


/* Check to see if header almost full; if so, print warning and increment,
   and record new value at its proper location in the file */


if( hdrpos >= (pos) (hdrsize - safety) )
{
/* now increase size of header area, and move the rest of the file up,
   adjust the current value of datapos accordingly */


 /* For now, comment out this next print statement
 printf("sdf_write: header area almost full, incrementing by hdrinit.\n");  
 fflush(stdout);
 */
 oldhdrsize=hdrsize;
 hdrsize+=hdrinit;
 hdrsztmp=hdrsize;
/*
 fseeko(fp,hdrpos-(pos)(sizeof(pos)),SEEK_SET);
*/
 fseeko(fp,(pos)(strlen(sdfstr)+1+2*sizeof(pos)+sizeof(i4)),SEEK_SET);
 if(ise) byteswap((void *)&hdrsztmp,(pos)1,sizeof(pos));
 fwrite(&hdrsztmp, sizeof(pos),1,fp);
 startloc=(pos) strlen(sdfstr)+1+sizeof(pos)+oldhdrsize;
 move_up(fp,hdrinit,startloc,datapos);
 fflush(fp);
/*
 fn=fileno(fp);
*/
 datapos += hdrinit;
 sdf_file_truncate(fp, (pos) datapos);
/*
 ftruncate(fn, (pos) datapos);
*/
}


/* check that id->datatype is printable, otherwise in big trouble and must
 * exit.  Actually is really only supposed to be f,i, or b acc. to sdf
 * manifesto, but we'll let folks slide here and use anything as long
 * as it's not blank. */


if (!isgraph(id->datatype))
{
        printf("sdf_write: error - id->datatype = 0x%02x not printable\n",
                id->datatype);
        fflush(stdout);
        fclose(fp);
        fotest=0;
        exit(1);
}


/* Now make sure that id->label is 'blank-free': */


rm_blanks(id->label);


/* Now, output the header information for the data to the file header area */


fseeko(fp,hdrpos,SEEK_SET);
data_to_string(id,fp);
hdrpos=ftello(fp);


/* Now, figure out the size of the data file from the info in the id struc: 
 * (Basically, just multiply all the dimensions together) 
 *
 * Note that this was done before the data_size fn was written, so prob.
 * should replace the next 5 lines of code with call to data_size.
*/



datasize=(pos) *(id->dims+0);  /* there's always at least one value */
for (i=1; i<id->ndim; i++)     /* sure hope zero trip loops work */
{
        datasize *=  (pos)( *(id->dims+i) );
}
if(id->datatype != 'c') /* test for complex var */
{
        if(ise) byteswap((void *)data, (pos) datasize, id->nbpw);
}
else
{
        /* kludge for complex var */
        if(ise) byteswap((void *)data, (pos) 2*datasize, id->nbpw);
}


/* Finally, write out the data: */


fseeko(fp,datapos,SEEK_SET);


if (id->datatype != 'c') /* test for complex variable type */
{
        /* Is next stmt 64-bit compatible? 

        Well, kind of.  arguments to fread, fwrite are of type size_t,
        which on a 32-bit machine is a 32-bit integer, while on a 64-bit
        machine is a 64-bit integer, at least in linux.  datasize is a
        64-bit integer.  In principle, one could get in trouble on a 32-bit
        machine, but since memory is limited to 2GB, datasize can't exceed
        that anyway, so the conversion from a (pos) integer to a size_t 
        integer should be OK. */

        fwrite(data,id->nbpw,datasize,fp); 
        /* for non-destructive write, need to btye-swap back! */
        if(ise) byteswap((void *)data, (pos) datasize, id->nbpw);
}
else
{
        fwrite(data,id->nbpw,(pos)2*datasize,fp); /*kludge for complex var */
        /* for non-destructive, need to btye-swap back! */
        if(ise) byteswap((void *)data, (pos) 2*datasize, id->nbpw);
}


datapos=(pos) ftello(fp);


/* Finish up: increment next_order, and update hdrpos and datapos in the file */


next_order++;
hdrtmp=hdrpos;
datatmp=datapos;
next_tmp=next_order;
hdrsztmp=hdrsize;


if(ise) byteswap((void *)&hdrtmp,(pos)1,sizeof(pos));
if(ise) byteswap((void *)&datatmp,(pos)1,sizeof(pos));
if(ise) byteswap((void *)&next_tmp,(pos)1,sizeof(i4));
if(ise) byteswap((void *)&hdrsztmp,(pos)1,sizeof(pos));


fseeko(fp, (pos) strlen(sdfstr)+1, SEEK_SET);


fwrite(&hdrtmp,sizeof(pos),1,fp);
fwrite(&datatmp,sizeof(pos),1,fp);
fwrite(&next_tmp,sizeof(i4),1,fp);
fwrite(&hdrsztmp,sizeof(pos),1,fp);


fflush(fp);
if(fotest)
{
        fctest=fclose(fp);
        fotest=0;
}
if(fctest != 0)
{
        printf("sdf_write:  failure to close file %s\n",fname);
        printf("sdf_write: errno = %d\n",errno);
        fotest=1;
}
 
/* We're done */


return 0;
}


i4 sdf_query(char *fname)
{



FILE *fp;
pos datasize;
char *sdfstr="SDF format";
char temp[301], testid[300];
i4 ise,ibe, iorder, ndatasets;
pos hdrsize;
data_id *id;


/* Make certain (pos) is defined as 64bit int, otherwise, exit */


test_sizes(); /* New Way */


/*  OLD WAY
if(sizeof(pos) != 8)
{
        printf("sdf_query: error building sdf_query: sizeof(pos) = %d not 8\n",
            sizeof(pos));
        exit(1);
}
*/


ibe = is_big_endian ();
ise = 0;
if (ibe == 0) ise = 1;    /* set flag for byteswapping if small endian */


/* Initially, open the file as read-only */


fp=fopen(fname,"rb");


if (fp == NULL) /* does file fname exist? */
{
        printf("sdf_query: file %s does not exist, exiting\n",fname);
        fflush(stdout);
        exit(1);
}
else
{
        fgets(testid,strlen(sdfstr)+1,fp);
        /* Now, test to make sure the SDF file identifier string is there */
        if(strncmp(sdfstr,testid,strlen(sdfstr)))
        {
                fclose(fp);
                printf("sdf_query: File %s is not an SDF file, exiting\n",
                        fname);
                fflush(stdout);
                exit(1);
        }
}
 
/* read the number of datasets in the file, and the header size: */


fseeko(fp, strlen(sdfstr)+1 + 2*sizeof(pos) , SEEK_SET);
fread(&ndatasets,sizeof(i4),1,fp);
if(ise) byteswap((void *)&ndatasets,(pos)1,sizeof(i4));
fread(&hdrsize,sizeof(pos),1,fp);
if(ise) byteswap((void *)&hdrsize,(pos)1,sizeof(pos));
/*
printf("ndatasets = %d\n",ndatasets);
*/


/* read the string data of all the datasets, print it out */


for (iorder=0; iorder<ndatasets; iorder++)
{
/*
        fgets(temp,300,fp);
        printf("iorder = %d , header string = %s",iorder,temp);
*/
        id=(data_id *)string_to_data(fp);
/*
        printf("sdf_query debug: id->label =%s",id->label);
*/
        datasize=data_size(id);
        output_int64(temp,datasize);
        printf("sdf_query: dataset %d, size %selements, descriptor: ",
               iorder,temp);
        fflush(stdout);
/* OLD WAY
        #if defined (__SVR4) && defined (__sun) 
        printf("sdf_query: dataset no. %d, datasize = %llu , descriptor = ",
             iorder, datasize);
        #elif defined WIN32
        printf("sdf_query: dataset no. %d, datasize = %I64u , descriptor = ",
             iorder, datasize);
        #else
        printf("sdf_query: dataset no. %d, datasize = %llu , descriptor = ",
             iorder, datasize);
        #endif
*/
        data_to_string(id,stdout);
        fflush(stdout);
        free(id->label);
        free(id->dims);
        free(id);
}


fclose(fp);
return ndatasets;
}

i4 sdf_ndat(char *fname)
{



FILE *fp;
pos datasize;
char *sdfstr="SDF format";
char temp[301], testid[300];
i4 ise,ibe, iorder, ndatasets;
pos hdrsize;
data_id *id;


/* Make certain (pos) is defined as 64bit int, otherwise, exit */


test_sizes(); /* Test definitions, New Way */

ibe = is_big_endian ();
ise = 0;
if (ibe == 0) ise = 1;    /* set flag for byteswapping if small endian */


/* Initially, open the file as read-only */


fp=fopen(fname,"rb");


if (fp == NULL) /* does file fname exist? */
{
        printf("sdf_ndat: file %s does not exist, exiting\n",fname);
        fflush(stdout);
        exit(1);
}
else
{
        fgets(testid,strlen(sdfstr)+1,fp);
        /* Now, test to make sure the SDF file identifier string is there */
        if(strncmp(sdfstr,testid,strlen(sdfstr)))
        {
                fclose(fp);
                printf("sdf_ndat: File %s is not an SDF file, exiting\n",
                        fname);
                fflush(stdout);
                exit(1);
        }
}
 
/* read the number of datasets in the file, and the header size: */


fseeko(fp, strlen(sdfstr)+1 + 2*sizeof(pos) , SEEK_SET);
fread(&ndatasets,sizeof(i4),1,fp);
if(ise) byteswap((void *)&ndatasets,(pos)1,sizeof(i4));
/*
printf("ndatasets = %d\n",ndatasets);
*/




fclose(fp);
return ndatasets;
}

void sdf_housekeeping(char *fname, i4 *norder, pos *hdrposout, 
                pos *dataposout, pos *hdrsizeout)
{


/* return the values of norder, hdrpos, datapos, and hdrsize to the
 * the caller.  Note they are passed by reference, not value.
*/


FILE *fp;
pos hdrpos, datapos;
char *sdfstr="SDF format";
char testid[300];
/* char temp[300]; */
i4 ise,ibe, ndatasets;
pos hdrsize;


/* Make certain (pos) is defined as 64bit int, otherwise, exit */


test_sizes(); 


ibe = is_big_endian ();
ise = 0;
if (ibe == 0) ise = 1;    /* set flag for byteswapping if small endian */


/* Initially, open the file as read-only */


fp=fopen(fname,"rb");


if (fp == NULL) /* does file fname exist? */
{
        printf("sdf_housekeeping: file %s does not exist, exiting\n",fname);
        fflush(stdout);
        exit(1);
}
else
{
        fgets(testid,strlen(sdfstr)+1,fp);
        /* Now, test to make sure the SDF file identifier string is there */
        if(strncmp(sdfstr,testid,strlen(sdfstr)))
        {
                fclose(fp);
                printf("sdf_housekeeping: File %s is not an SDF file, exiting\n"
                        , fname);
                fflush(stdout);
                exit(1);
        }
}
/* read hdrpos,datapos,ndatasets,hdrsize: */
fseeko(fp, strlen(sdfstr)+1, SEEK_SET);


fread(&hdrpos,sizeof(pos),1,fp);

/* debug
output_int64(temp,hdrpos);
printf("sdf_housekeeping - before byteswap: hdrpos = %s\n",temp);
printf("sdf_housekeeping - before byteswap, ise = %d\n",ise);
fflush(stdout);
 end debug*/

if(ise) byteswap((void *)&hdrpos,(pos)1,sizeof(pos));

/* debug
output_int64(temp,hdrpos);
printf("sdf_housekeeping - after byteswap: hdrpos = %s\n",temp);
fflush(stdout);
end debug */


fread(&datapos,sizeof(pos),1,fp);
if(ise) byteswap((void *)&datapos,(pos)1,sizeof(pos));


fread(&ndatasets,sizeof(i4),1,fp);
if(ise) byteswap((void *)&ndatasets,(pos)1,sizeof(i4));


fread(&hdrsize,sizeof(pos),1,fp);
if(ise) byteswap((void *)&hdrsize,(pos)1,sizeof(pos));


/*
 printf("sdf_housekeeping: existing datasets = %d\n",ndatasets);
 fflush(stdout);
*/


/* debug */
/*
output_int64(temp,hdrsize);
printf("sdf_housekeeping: hdrsize = %s\n",temp);
fflush(stdout);
*/


/*
 * Now set output variables to the ones read in from file:
*/
*norder=ndatasets;
*hdrposout = hdrpos;
*dataposout = datapos;
*hdrsizeout = hdrsize;


fclose(fp);
return;
}


void sdf_labmatch_f77__(char *fname_f77, intf * ndat_f77,
     char *labtest_f77, intf *match, intf_clen len_fnamef77, 
     intf_clen len_labtestf77)
{
sdf_labmatch_f77(fname_f77, ndat_f77, labtest_f77, match, len_fnamef77,
    len_labtestf77);
return;
}
void sdf_labmatch_f77_(char *fname_f77, intf * ndat_f77,
     char *labtest_f77, intf *match, intf_clen len_fnamef77, 
     intf_clen len_labtestf77)
{
sdf_labmatch_f77(fname_f77, ndat_f77, labtest_f77, match, len_fnamef77,
    len_labtestf77);
return;
}


void sdf_labmatch_f77(char *fname_f77, intf * ndat_f77,
     char *labtest_f77, intf *match, intf_clen len_fnamef77, 
     intf_clen len_labtestf77)
{
i4 i;
i4 ndat_temp;
i4 *matchind;
char *fname, *labtest;


/* DEBUG - print out lengths of fortran stings (hidden args) */
/*
printf("len_fnamef77 = %d\n",len_fnamef77);
printf("len_labf77 = %d\n",len_labf77);
printf("len_dtf77 = %d\n",len_dtf77);
*/


/* Convert fnamef77 to C-string fname */
fs2cs(&fname,fname_f77,len_fnamef77);
fs2cs(&labtest,labtest_f77,len_labtestf77);
/* DEBUG - print out C-string version of filename fname */
/*
printf("fname = %s\n",fname);
*/
matchind = sdf_labmatch(fname,&ndat_temp,labtest);
*ndat_f77= (intf) ndat_temp;
for (i=0;i<ndat_temp;i++)
{
        *(match + i) = (intf) *(matchind+i);
}


free(fname);
free(labtest);
free(matchind);
return;
}


i4 * sdf_labmatch(char *fname, i4 *ndat, char * labtest)
{
/*
Returns a pointer to an array of dataset indices for which 
id->label matches labtest.  
ndat is determined independently, and its value set in *ndat.  Note that in
the caller, ndat must be allocated before call to sdf_labmatch, or
alternatively its address can be passed to sdf_labmatch.


This function allocates and returns the space for the integer array of
matched indices, and is always of size ndat.  It will be filled with 
values of -1 after the end of matched indices reached.  Be sure to free
the returned index array in caller when finished with it.


*/


FILE *fp;
char *sdfstr="SDF format";
char testid[300];
i4 ise,ibe, iorder, ndatasets;
i4 match, nclabtest, nclab, ncidlab, *matchind;
pos hdrsize;
data_id *id;


/* Make certain (pos) is defined as 64bit int, otherwise, exit */


test_sizes(); /* New way */


/*  OLD WAY
if(sizeof(pos) != 8)
{
        printf("sdf_sizes: error building sdf_sizes: sizeof(pos) = %d not 8\n",
            sizeof(pos));
        exit(1);
}
*/


ibe = is_big_endian ();
ise = 0;
if (ibe == 0) ise = 1;    /* set flag for byteswapping if small endian */


/* Initially, open the file as read-only */


fp=fopen(fname,"rb");


if (fp == NULL) /* does file fname exist? */
{
        printf("sdf_labmatch: file %s does not exist, exiting\n",fname);
        fflush(stdout);
        exit(1);
}
else
{
        fgets(testid,strlen(sdfstr)+1,fp);
        /* Now, test to make sure the SDF file identifier string is there */
        if(strncmp(sdfstr,testid,strlen(sdfstr)))
        {
                fclose(fp);
                printf("sdf_sizes: File %s is not an SDF file, exiting\n",
                        fname);
                fflush(stdout);
                exit(1);
        }
}
 
/* read the number of datasets in the file: */


fseeko(fp, strlen(sdfstr)+1 + 2*sizeof(pos) , SEEK_SET);
fread(&ndatasets,sizeof(i4),1,fp);
if(ise) byteswap((void *)&ndatasets,(pos)1,sizeof(i4));
fread(&hdrsize,sizeof(pos),1,fp);
if(ise) byteswap((void *)&hdrsize,(pos)1,sizeof(pos));


/* set *ndat to ndatasets */
*ndat = ndatasets;


/*
 DEBUG - printf("sdf_size: *ndat = %d\n",*ndat);
*/


/* Now, allocate the memory for matchind */


matchind=malloc(ndatasets*sizeof(i4));
for (iorder=0; iorder<ndatasets; iorder++)
{
        *(matchind+iorder)=(i4) (-1);
}


match=0;
nclabtest=strlen(labtest);
for (iorder=0; iorder<ndatasets; iorder++)
{
        id=(data_id *)string_to_data(fp);
/*
          printf("sdf_query debug: id->label =%s",id->label);
*/
        ncidlab=strlen(id->label);
        nclab = (nclabtest > ncidlab) ? nclabtest : ncidlab;
        if(!strncmp(id->label,labtest,nclab))
        {
                *(matchind+match)=iorder;
                match++;
        }
        
        free(id->label);
        free(id->dims);
        free(id);
}


fclose(fp);
return matchind;
}
i4 sdf_sizes(char *fname, i4 *ndat, char ** datatypes, pos ** datasizes,
        i4 **nbpw)
{
/*
Returns *ndat, the number of datasets, the char array **datatypes,
the pos array **datasizes, and the i4 array ** nbpw, reflecting the 
*ndat values of datatypes, datasizes, and nbpw values from each dataset
in the sdf file.  Note **datatypes is a true character array, not a string!
After calling sdf_sizes, calling routine responsible for freeing malloc'ed
arrays **datatypes, **datasizes, **nbpw .
Note also that ndat must be allocated BEFORE this function is called.
*/


FILE *fp;
char *sdfstr="SDF format";
char testid[300];
i4 ise,ibe, iorder, ndatasets;
pos hdrsize;
data_id *id;


/* Make certain (pos) is defined as 64bit int, otherwise, exit */


test_sizes(); /* New way */


/*  OLD WAY
if(sizeof(pos) != 8)
{
        printf("sdf_sizes: error building sdf_sizes: sizeof(pos) = %d not 8\n",
            sizeof(pos));
        exit(1);
}
*/


ibe = is_big_endian ();
ise = 0;
if (ibe == 0) ise = 1;    /* set flag for byteswapping if small endian */


/* Initially, open the file as read-only */


fp=fopen(fname,"rb");


if (fp == NULL) /* does file fname exist? */
{
        printf("sdf_sizes: file %s does not exist, exiting\n",fname);
        fflush(stdout);
        exit(1);
}
else
{
        fgets(testid,strlen(sdfstr)+1,fp);
        /* Now, test to make sure the SDF file identifier string is there */
        if(strncmp(sdfstr,testid,strlen(sdfstr)))
        {
                fclose(fp);
                printf("sdf_sizes: File %s is not an SDF file, exiting\n",
                        fname);
                fflush(stdout);
                exit(1);
        }
}
 
/* read the number of datasets in the file: */


fseeko(fp, strlen(sdfstr)+1 + 2*sizeof(pos) , SEEK_SET);
fread(&ndatasets,sizeof(i4),1,fp);
if(ise) byteswap((void *)&ndatasets,(pos)1,sizeof(i4));
fread(&hdrsize,sizeof(pos),1,fp);
if(ise) byteswap((void *)&hdrsize,(pos)1,sizeof(pos));


/* set *ndat to ndatasets */
*ndat = ndatasets;


/*
 DEBUG - printf("sdf_size: *ndat = %d\n",*ndat);
*/


/* Now, allocate the memory for *datatypes, *datasizes, *nbpw: */


*datatypes=malloc(ndatasets*sizeof(char));
*datasizes=malloc(ndatasets*sizeof(pos));
*nbpw=malloc(ndatasets*sizeof(i4));


for (iorder=0; iorder<ndatasets; iorder++)
{
/*
        fgets(temp,300,fp);
        printf("iorder = %d , header string = %s",iorder,temp);
*/
        id=(data_id *)string_to_data(fp);
/*
          printf("sdf_query debug: id->label =%s",id->label);
*/
        *(*datasizes+iorder)=data_size(id);
        *(*datatypes+iorder)=id->datatype;
        *(*nbpw+iorder)=id->nbpw;
        
        free(id->label);
        free(id->dims);
        free(id);
}


fclose(fp);
return 0;
}


void sdf_sizes_f77__(char *fnamef77, intf *ndatf77, char dtypesf77[][1],
intf8 * datasizesf77, intf *nbpwf77, intf_clen len_fnamef77, 
intf_clen len_dtypesf77)
{
sdf_sizes_f77(fnamef77, ndatf77, dtypesf77,
datasizesf77, nbpwf77, len_fnamef77, len_dtypesf77);
return;
}
void sdf_sizes_f77_(char *fnamef77, intf *ndatf77, char dtypesf77[][1],
intf8 *datasizesf77, intf *nbpwf77, intf_clen len_fnamef77, 
intf_clen len_dtypesf77)
{
sdf_sizes_f77(fnamef77, ndatf77, dtypesf77,
datasizesf77, nbpwf77, len_fnamef77, len_dtypesf77);
return;
}


void sdf_sizes_f77(char *fnamef77, intf *ndatf77, char dtypesf77[][1],
intf8 * datasizesf77, intf *nbpwf77, intf_clen len_fnamef77, 
intf_clen len_dtypesf77)
{


/*
Fortran 77 call-able version of sdf_sizes - returns size information on
each dataset in an sdf file 
*/
char *dtypes;
pos *dsizes;
i4 *nbpw;
i4 *ndat;
i4 i;
char *fname;
ndat=malloc(sizeof(i4)); /*note ndat must be allocated before sdf_sizes called*/
fs2cs(&fname,fnamef77,len_fnamef77); /*convert fname from fstring to cstring */
sdf_sizes(fname, ndat, &dtypes, &dsizes, &nbpw); /* get dataset sizes */


/*
 DEBUG - printf("sdf_sizes_f77: *ndat = %d\n",*ndat);
*/


for (i=0;i<*ndat;i++)
{
        *(nbpwf77+i)= (intf) *(nbpw+i);
        *(datasizesf77+i)=(intf8) *(dsizes+i); /* converted to intf for f77 */
        dtypesf77[i][0]=(char)*(dtypes+i);
}
*ndatf77= (intf) *ndat;


/* free all the temporary junk created in this function */


free(ndat);
free(fname);
free(dtypes);
free(nbpw);
free(dsizes);
return;
}
void sdf_rm_f77__(char *fnamef77, intf_clen len_fnamef77)
{
sdf_rm_f77(fnamef77, len_fnamef77);
return;
}
void sdf_rm_f77_(char *fnamef77, intf_clen len_fnamef77)
{
sdf_rm_f77(fnamef77, len_fnamef77);
return;
}
void sdf_rm_f77(char *fnamef77, intf_clen len_fnamef77)
{
/* Fortran callable function to remove sdf files */
char *fname;
fs2cs(&fname,fnamef77,len_fnamef77);
sdf_rm(fname);
free(fname);
return;
}


i4 sdf_rm(char *fname)
/*
Function to remove the file fname, but only if it's an sdf file.
*/
{
FILE *fp;
char *sdfstr="SDF format";
char testid[300];


fp=fopen(fname,"rb");


if (fp == NULL) /* does file fname exist? */
{
/*  This print statement is annoying, since this isn't typically
an error.  Will comment it out for now.

        printf("sdf_rm: file %s cannot be opened for deletion, returning\n",
          fname);
        fflush(stdout);
*/
        return 1;
}


fgets(testid,strlen(sdfstr)+1,fp);
/* Now, test to make sure the SDF file identifier string is there */
if(strncmp(sdfstr,testid,strlen(sdfstr)))
{
        fclose(fp);
        printf("sdf_rm: File %s is not an SDF file, cannot remove\n",
                fname);
        fflush(stdout);
        exit(1);
}
fclose(fp);
remove(fname);
return 0;
}


void sdf_query_f77__(char *fnamef77, intf *ndat, intf_clen len_fname_f77)
{
   sdf_query_f77(fnamef77, ndat, len_fname_f77);
   return;
}
void sdf_query_f77_(char *fnamef77, intf *ndat, intf_clen len_fname_f77)
{
   sdf_query_f77(fnamef77, ndat, len_fname_f77);
   return;
}
void  sdf_query_f77(char *fnamef77, intf * ndat, intf_clen len_fname_f77)
{


/* This is the fortran callable version of sdf_query. *ndat is
the number of datasets in the file. */


char *fname;

/* Make certain (pos) is defined as 64bit int, otherwise, exit */


test_sizes(); /* New way */


/* OLD WAY
if(sizeof(pos) != 8)
{
   printf("sdf_query_f77: error building, sizeof(pos) = %d not 8\n",
            sizeof(pos));
   exit(1);
}
*/


/*
OK, now convert the Fortran string containing the name of the file to the
C-string containing the name of the file, and call sdf_query:
*/


fs2cs(&fname, fnamef77, len_fname_f77);
*ndat=sdf_query(fname);
free(fname);
return;
}

void sdf_ndat_f77__(char *fnamef77, intf *ndat, intf_clen len_fname_f77)
{
   sdf_ndat_f77(fnamef77, ndat, len_fname_f77);
   return;
}
void sdf_ndat_f77_(char *fnamef77, intf *ndat, intf_clen len_fname_f77)
{
   sdf_ndat_f77(fnamef77, ndat, len_fname_f77);
   return;
}
void  sdf_ndat_f77(char *fnamef77, intf * ndat, intf_clen len_fname_f77)
{


/* This is the fortran callable version of sdf_query. *ndat is
the number of datasets in the file. */


char *fname;

/* Make certain (pos) is defined as 64bit int, otherwise, exit */


test_sizes(); /* New way */


/*
OK, now convert the Fortran string containing the name of the file to the
C-string containing the name of the file, and call sdf_query:
*/


fs2cs(&fname, fnamef77, len_fname_f77);
*ndat=sdf_ndat(fname);
free(fname);
return;
}

data_id *sdf_read(char *fname, i4 order, void **data)
{


/* read in the dataset corresponding to order, and return idnew */


FILE *fp;
data_id *id=NULL;
pos datasize=(pos)0, curpos;
char *sdfstr="SDF format";
char temp[301], testid[300];
i4 ise,ibe, iorder, ndatasets;
pos frtest;
pos hdrsize;


/* Make certain (pos) is defined as 64bit int, otherwise, exit */


test_sizes(); /* New way */


/* OLD WAY
if(sizeof(pos) != 8)
{
        printf("sdf_read: error building sdf_read: sizeof(pos) = %d not 8\n",
            sizeof(pos));
        exit(1);
}
*/


ibe = is_big_endian ();
ise = 0;
if (ibe == 0) ise = 1;    /* set flag for byteswapping if small endian */


/* Initially, open the file as read-only */


fp=fopen(fname,"rb");


if (fp == NULL) /* does file fname exist? */
{
        printf("sdf_read: file %s does not exist, exiting\n",fname);
        fflush(stdout);
        exit(1);
}
else
{
        fgets(testid,strlen(sdfstr)+1,fp);
        /* Now, test to make sure the SDF file identifier string is there */
        if(strncmp(sdfstr,testid,strlen(sdfstr)))
        {
                fclose(fp);
                printf("sdf_read: File %s is not an SDF file, exiting\n",
                        fname);
                fflush(stdout);
                exit(1);
        }
}
/* read the number of datasets in the file: */


fseeko(fp, strlen(sdfstr)+1 + 2*sizeof(pos) , SEEK_SET);
fread(&ndatasets,sizeof(i4),1,fp);
if(ise) byteswap((void *)&ndatasets,(pos)1,sizeof(i4));
/*
 Next print statement possibly important but nuke for now
*/


/*
printf("sdf_read: existing datasets = %d\n",ndatasets);
fflush(stdout);
*/


/*
printf("sdf_read: desired dataset order = %d\n",order);
*/
if(order > (ndatasets-1))
{
        printf("sdf_read:  Desired dataset order %d > ndatasets-1 = %d\n",
             order, ndatasets-1);
        fflush(stdout);
        exit(1);
}
if(order < 0)
{
        printf("sdf_read:  Desired dataset order %d < 0 \n",
             order);
        fflush(stdout);
        exit(1);
}
/* read in header size */
fread(&hdrsize,sizeof(pos),1,fp);
if(ise) byteswap((void *)&hdrsize,(pos)1,sizeof(pos));
/* debug */
/*
output_int64(temp,hdrsize);
printf("sdf_read: hdrsize = %s\n",temp);
fflush(stdout);
*/


/* read the string data of all the datasets, print it out */
curpos = (pos) strlen(sdfstr)+1+sizeof(pos)+hdrsize;
/*
        printf("DEBUG - sdf_read: curpos = %qu\n",curpos);
*/
for (iorder=0; iorder <= order; iorder++)
{
/*
        fgets(temp,300,fp);
        printf("iorder = %d , header string = %s",iorder,temp);
*/


        id=string_to_data(fp);
        datasize=data_size(id);



        if(iorder < order)
        {
                if(id->datatype != 'c') /* test for complex var */
                {
                        curpos += (pos) datasize*(id->nbpw);
                }
                else
                {
                        curpos += (pos) datasize*2*(id->nbpw);
                        /* kludge for complex variables */
                }
                free(id->label);
                free(id->dims);
                free(id);
        }
        else
        {
        /* NEW WAY */
        output_int64(temp,datasize);
        /* next print statements possibly important, but nuke for now */
        /*
        printf("sdf_read: dataset %d, size %selements, descriptor: ",
               iorder,temp);
        fflush(stdout);
        */


        /* OLD WAY 
        #if defined (__SVR4) && defined (__sun) 
        printf("sdf_read: dataset no. %d, datasize = %llu , descriptor = ",
             iorder, datasize);
        #elif defined WIN32
        printf("sdf_read: dataset no. %d, datasize = %I64u , descriptor = ",
             iorder, datasize);
        #else
        printf("sdf_read: dataset no. %d, datasize = %llu , descriptor = ",
             iorder, datasize);
        #endif
        */
        /*
        data_to_string(id,stdout);
        printf("\n");
        fflush(stdout);
        */
        }
}


fseeko(fp,curpos,SEEK_SET);                /* reposition file pointer  */


if (id->datatype != 'c') /* test for complex var */
{
        *data=malloc(datasize* (pos) (id->nbpw));         /* allocate data */
        frtest=(pos) fread(*data, id->nbpw , datasize,fp); /* read the data */
	if(frtest < datasize)
	{
           /* test for unsuccessful read: */
           printf("sdf_read: read error\n");
           printf("sdf_read: errno = %d\n",errno);
           output_int64(temp,frtest);
           printf("sdf_read: frtest = %s\n",temp);
           output_int64(temp,datasize);
           printf("sdf_read: datasize =%s\n",temp);
           fflush(stdout);
           exit(1);
	}
        if(ise) byteswap(*data,datasize,id->nbpw); /* byteswap if neeeded */
}
else
{
        /* kludge for reading complex variable data */
        *data=malloc(datasize*(pos)2*(id->nbpw)); 
        frtest=(pos) fread(*data, id->nbpw , 2*datasize,fp); /* read the data */
	if(frtest < 2*datasize)
	{
           /* test for unsuccessful read: */
           printf("sdf_read: read error\n");
           printf("sdf_read: errno = %d\n",errno);
           output_int64(temp,frtest);
           printf("sdf_read: frtest = %s\n",temp);
           output_int64(temp,(pos)2*datasize);
           printf("sdf_read(complex): 2*datasize =%s\n",temp);
           fflush(stdout);
           exit(1);
	}
        if(ise) byteswap(*data,2*datasize,id->nbpw); /* byteswap if neeeded */
}


fclose(fp);
return id;


}


void sdf_details_f77__(char *fnamef77, intf *iorder, char *labf77, char *dtf77, 
     intf *nbpwf77, intf *ndimf77, intf8 *dimsf77,
     intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77)
{
sdf_details_f77(fnamef77, iorder, labf77, dtf77, 
     nbpwf77, ndimf77, dimsf77,
     len_fnamef77, len_labf77, len_dtf77);
return;
}
void sdf_details_f77_(char *fnamef77, intf *iorder, char *labf77, char *dtf77, 
     intf *nbpwf77, intf *ndimf77, intf8 *dimsf77,
     intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77)
{
sdf_details_f77(fnamef77, iorder, labf77, dtf77, 
     nbpwf77, ndimf77, dimsf77,
     len_fnamef77, len_labf77, len_dtf77);
return;
}


void sdf_details_f77(char *fnamef77, intf *iorder, char *labf77, char *dtf77, 
     intf *nbpwf77, intf *ndimf77, intf8 *dimsf77,
     intf_clen len_fnamef77, intf_clen len_labf77, intf_clen len_dtf77)


{
data_id *id;
i4 i;
char *fname;


/* DEBUG - print out lengths of fortran stings (hidden args) */

/*
printf("len_fnamef77 = %d\n",len_fnamef77);
printf("len_labf77 = %d\n",len_labf77);
printf("len_dtf77 = %d\n",len_dtf77);
*/

/* Convert fnamef77 to C-string fname */
fs2cs(&fname,fnamef77,len_fnamef77);
/* DEBUG - print out C-string version of filename fname */
/*
printf("fname = %s\n",fname);
*/


id=(data_id *) sdf_details(fname, (i4) *iorder);
/*
 id=(data_id *) sdf_read(fname, *iorder, (void **)&dataf77);
*/

/* now convert id->label into fortran string labf77 */
cs2fs(id->label,labf77,len_labf77);


*dtf77 = id->datatype;
*nbpwf77 = (intf) id->nbpw;
*ndimf77 = (intf) id->ndim;
for (i=0; i < id->ndim; i++)
{
        /* flip order for f77 */
        *(dimsf77+i) = (intf8) *(id->dims+id->ndim-1 -i); 
}


/* rather than go through the next intermediate stage, could call sdf_read with
 * dataf77 directly, maybe?  ...tried that, didn't work, but prob. needs
 * further investigation ... */


free(fname);
free(id->label);
free(id->dims);
free(id);
return;
}
data_id *sdf_details(char *fname, i4 order)
{


/* return the id structure corresponding to dataset iorder */


FILE *fp;
data_id *id=NULL;
pos datasize;
char *sdfstr="SDF format";
char temp[301], testid[300];
i4 ise,ibe, iorder, ndatasets;
pos hdrsize;


/* Make certain (pos) is defined as 64bit int, otherwise, exit */


test_sizes(); 


ibe = is_big_endian ();
ise = 0;
if (ibe == 0) ise = 1;    /* set flag for byteswapping if small endian */


/* Initially, open the file as read-only */


fp=fopen(fname,"rb");


if (fp == NULL) /* does file fname exist? */
{
        printf("sdf_details: file %s does not exist, exiting\n",fname);
        fflush(stdout);
        exit(1);
}
else
{
        fgets(testid,strlen(sdfstr)+1,fp);
        /* Now, test to make sure the SDF file identifier string is there */
        if(strncmp(sdfstr,testid,strlen(sdfstr)))
        {
                fclose(fp);
                printf("sdf_details: File %s is not an SDF file, exiting\n",
                        fname);
                fflush(stdout);
                exit(1);
        }
}
/* read the number of datasets in the file: */


fseeko(fp, strlen(sdfstr)+1 + 2*sizeof(pos) , SEEK_SET);
fread(&ndatasets,sizeof(i4),1,fp);
if(ise) byteswap((void *)&ndatasets,(pos)1,sizeof(i4));
/* Nuke the copious output
printf("sdf_details: existing datasets = %d\n",ndatasets);
fflush(stdout);
*/


/*
printf("sdf_details: desired dataset order = %d\n",order);
*/
if(order > (ndatasets-1))
{
        printf("sdf_details:  Desired dataset order %d > ndatasets-1 = %d\n",
             order, ndatasets-1);
        fflush(stdout);
        exit(1);
}
if(order < 0)
{
        printf("sdf_details:  Desired dataset order %d < 0 \n",
             order);
        fflush(stdout);
        exit(1);
}
/* read in header size */
fread(&hdrsize,sizeof(pos),1,fp);
if(ise) byteswap((void *)&hdrsize,(pos)1,sizeof(pos));
/* debug */
/*
output_int64(temp,hdrsize);
printf("sdf_details: hdrsize = %s\n",temp);
fflush(stdout);
*/


/* read the string data of all the datasets, print it out */
/* Now nuking all the printout starting in version 0.75 */


for (iorder=0; iorder <= order; iorder++)
{
/*
        fgets(temp,300,fp);
        printf("iorder = %d , header string = %s",iorder,temp);
*/


        id=string_to_data(fp);
        datasize=data_size(id);



        if(iorder < order)


        {
                free(id->label);
                free(id->dims);
                free(id);
        }


        else


        {
        output_int64(temp,datasize);
/*  Nuke the copious output for now:
        printf("sdf_details: dataset %d, size %selements, descriptor: ",
               iorder,temp);
        fflush(stdout);
        data_to_string(id,stdout);
        printf("\n");
        fflush(stdout);
*/
        }
}


fclose(fp);
return id;
}


data_id * string_to_data( FILE *fp )
{


/* Read in header string, extract information, return it as a data_id struc */


        char temp[301];
        int ipos,i,ii;
        char *p, *last;
        char * delim = " ";
        data_id *idnew;


        idnew=malloc(sizeof(data_id)); /* create space for structure idnew */
        fgets(temp, 300, fp);          /* read in header string */


/*
        printf("DEBUG - string_to_data: temp = %s\n",temp);
*/


        p=(char *)strtok_r(temp,delim,&last);  /* get 1st token from string */


        /* loop over remaining tokens in the string and put them into idnew */


        i=0;
        ipos=0;
        idnew->order=atol(p);
        while ((p = (char *)strtok_r(NULL,delim,&last)) != NULL)
        {
                i++;
/*
                printf(" string_to_data DEBUG - i = %u, token = %s\n",i,p);
                printf("strlen(p) = %d\n",strlen(p));
*/
                if(i == 1) 
                {
                        idnew->label=calloc(strlen(p)+1,1);
/*
                        strncpy(idnew->label,p,strlen(p));
*/
                        for (ii=0;ii<=strlen(p);ii++)
                        {
                                *(idnew->label+ii)=*(p+ii);
                        }
                }
                if(i == 2) idnew->datatype=*p;
                if(i == 3) idnew->nbpw=atol(p);
                if(i == 4) idnew->ndim=atol(p);


                /* now create space for the individual dimension values */


                if (i == 5) idnew->dims=malloc(idnew->ndim*sizeof(pos));


                if (i >= 5)
                {
                        if(ipos < idnew->ndim)
                        {
                                /* read in the array dimensions */
                                *(idnew->dims+ipos)= (pos) atopos(p);   
                                ipos++;
                        }
                }
        }
        return idnew;
}
int data_to_buff(data_id *id, char *buff, pos *offset)
{


/* take information from the id structure, and write it as string
to the buffer buff; return the new offset.  This function very similar
to data_to_string function, which writes the string directly to file */


        char *na = "NA";
        /* char *p, *label; */   /* not used anymore */
        pos ipos;
        i4 i,blank_label,ic,lablen;
        ipos=0;
        sprintf((buff+ipos),"%u ",(u4) id->order);
        ipos=(pos) strlen(buff);


        /* Here, need to check that label is not a blank or null string; 
        if so, set to "NA". */


        blank_label=1;
        lablen=strlen(id->label);
        for (ic=0; ic < lablen; ic++)
        {
                if( isprint(*(id->label+ic)) && ! isspace(*(id->label+ic)) )
                     blank_label=0;
        }
/*
        if(blank_label) id->label="NA";
*/
        if(blank_label) strncpy(id->label,na,2);


        sprintf((buff+ipos),"%s ",id->label);


/*
        printf("inside data_to_buff, id->label = %s\n",id->label);
*/


        ipos=(pos) strlen(buff);
        sprintf((buff+ipos),"%c ",id->datatype);
        ipos=(pos) strlen(buff);
        sprintf((buff+ipos),"%u ",(u4) id->nbpw);
        ipos=(pos) strlen(buff);
        sprintf((buff+ipos),"%u ",(u4) id->ndim);
        for (i=0;i<id->ndim;i++)
        {
                ipos=(pos) strlen(buff);
                output_int64((buff+ipos),*(id->dims+i)); /* new way */
        }


        ipos=(pos) strlen(buff);
        *(buff+ipos)='\n'; /* write linefeed to buffer */
        ipos++;
/*
        printf("data_to_buff: ipos = %d\n",ipos);
*/
        *offset=ipos;
        return 0;
}
int data_to_string(data_id *id, FILE *fp)
{


/* take information from the id structure, and write it to header as string */


        char temp[300];
        char *na = "NA";
        /* char *p, *label; */  /* not used anymore */
        i4 ipos,i,blank_label,ic,lablen;
        ipos=0;
        sprintf((temp+ipos),"%u ",(u4) id->order);
        ipos=strlen(temp);


        /* Here, need to check that label is not a blank or null string; 
        if so, set to "NA". */


        blank_label=1;
        lablen=strlen(id->label);
        for (ic=0; ic < lablen; ic++)
        {
                if( isprint(*(id->label+ic)) && ! isspace(*(id->label+ic)) )
                     blank_label=0;
        }
/*
        if(blank_label) id->label="NA";
*/
        if(blank_label) strncpy(id->label,na,2);


        sprintf((temp+ipos),"%s ",id->label);


/*
        printf("inside data_to_string, id->label = %s\n",id->label);
*/


        ipos=strlen(temp);
        sprintf((temp+ipos),"%c ",id->datatype);
        ipos=strlen(temp);
        sprintf((temp+ipos),"%u ",(u4) id->nbpw);
        ipos=strlen(temp);
        sprintf((temp+ipos),"%u ",(u4) id->ndim);
        for (i=0;i<id->ndim;i++)
        {
                ipos=strlen(temp);
                output_int64((temp+ipos),*(id->dims+i)); /* new way */


                /*  OLD WAY
                #if defined (__SVR4) && defined (__sun) 
                   sprintf((temp+ipos),"%llu ",*(id->dims + i));
                #elif defined WIN32
                   sprintf((temp+ipos),"%I64u ",*(id->dims + i));
                #else
                   sprintf((temp+ipos),"%llu ",*(id->dims + i));
                #endif
                */
        }
/*
        fprintf(fp, "%s", temp); 
*/
        /* write out the string to the file */
        /* add a line feed */
/*
        fputc((unsigned char)'\n',fp); 
*/


        fprintf(fp, "%s\n", temp); /* write string and linefeed to file */
        return 0;
}


int output_int64(char *buff, pos ill)


/* Function to print a 64-bit integer into the string buff,
   since there seems to be no platform independent format stmt for 64 bit 
   integers.  Adds a trailing blank to the string. User must make sure
   buff is dimensioned or allocated big enough in calling program.  */


{
char work[100];
int neg=0;
int rem,i;
int counter=0;


if(ill < 0)
{
        ill=-ill;
        neg=1;
}


while(ill > 0)
{
        /* digit by digit, get the chars for each decimal place into work */
        rem= (int) (ill % (pos)10);
        ill /= (pos)10;
        sprintf(&work[counter],"%d",rem);
        counter++;
}


if(neg)
{
        work[counter]='-'; /* add - sign for neg. integers */
        counter++;
}


if(counter == 0)
{
        /* make sure that if integer = 0, that at least a 0 gets printed */
        work[counter]='0';
        counter++;
}


counter++; /* add one extra character for the trailing blank */


for (i=0;i<counter-1;i++)
{
        *(buff+counter-i-2)=work[i]; /* digits and minus sign if it exists */
}
        *(buff+counter-1)=' '; /* trailing blank */
        *(buff+counter)='\0'; /* null terminator */


return 0;
}


pos atopos(char *s)
/* Try to make a platform-independent version of atoll, which doesn't work
on Loraine's MAC laptop running OS 10.3.9 .  Returns the value of the
64-bit integer corresponding to the ascii string s. */
{
i4 nchar;
i4 i;
pos value=0;
i4 neg=0;
nchar=strlen(s);
/* should put in an upper limit on nchar */
if(nchar > 19)
{
        printf("atopos: string too long to convert to pos, > 19 chars\n");
        exit(1);
}
for(i=0;i<nchar;i++)
{
        if((s[i] >= '0') && (s[i] <= '9')) 
        {
                value=10*value+(s[i]-'0');
                /* build up the integer value digit by digit */
        }
        if(s[i] == '-') neg=1;
}
if(neg) value=-value;
return value;
}


pos data_size(data_id * id)


/* compute the size of the dataset (in elements) from the id structure */


{
   pos datasize;
   i4 i;
   datasize=(pos) *(id->dims+0);  /* there's always at least one value */
   for (i=1; i<id->ndim; i++)     /* sure hope zero trip loops work */
   {
        datasize *=  (pos)( *(id->dims+i) );
   }
   return datasize;
}


i4 cs2fs(char *cs, char *fs, intf_clen len_fs)
{
/*
  Convert the C string *cs to the Fortran string *fs.  Here it is assumed
  the C string *cs already exists, and the fortran string *fs has already
  been allocated by the fortran calling program, and has length len_fs.
*/
i4 i,csl;
if (len_fs <= 0)
{
        printf("cs2fs:  Bad value of len_fs = %d\n",(int)len_fs);
        fflush(stdout);
        exit(1);
}
csl=strlen(cs);
if ( (csl < 0) || (csl > len_fs))
{
        printf("cs2fs:  Bad value of csl = %d\n",csl);
        fflush(stdout);
        exit(1);
}
for (i=0;i<csl;i++)
{
        *(fs+i)=*(cs+i);
}
for (i=csl;i<len_fs;i++)
{
        *(fs+i)=' ';
}
return 0;
}


i4 fs2cs(char **cs, char *fs, intf_clen len_fs)
{
/*
   Convert the fortran string *fs to the c-style string *cs.  Fortran strings
   are right-filled with blanks; C strings are null terminated.
   Start from end of fortran string, find first non-blank character, then
   copy from the start to that character into the c-style string, add null
   character.  The parameter len_fs has been acquired from the fortran
   call to a C function already.  This malloc's *cs as a new string, so
   must free cs in the calling function when done.

   Note that if len_fs == 0 then this will just create a null C-string.

*/
i4 i,csl;
if (len_fs < 0)
{
        printf("fs2cs:  Bad value of len_fs = %d\n",(int)len_fs);
        fflush(stdout);
        exit(1);
}
/* csl=len_fs;  -- results in bug when converting all blank fortran strings */
csl=0; /* Hopefully this fixes the problem */
if (len_fs > 0)
{
   for (i=len_fs-1; i >= 0; i--) /* this should be a while loop - cleaner */
   {
           if(*(fs+i) != ' ') 
           {
                   csl=i+1;
                   i=0;
           }
   }
}

/*
 printf("string fs has nontrivial length of %d\n",csl);
*/
if(csl <=0)
{
    printf("fs2cs in sdf: fortran input string blank, exiting\n");
    exit(1);
    
}
*cs=malloc(csl+1);
for (i=0;i<csl;i++)
{
        *(*cs+i)=*(fs+i);
}
*(*cs+csl)='\0';
return 0;
}


void rm_blanks(char *s)
{
/* remove any blanks from a character string.  This is needed for a blank-free
 * string in id->label in sdf_write */
char *work;
i4 i,ls,inb;
ls=strlen(s);
work=malloc(ls+1);
inb=0;
for (i=0;i<ls;i++)
{
        if(*(s+i) != ' ')
        {
                *(work+inb)=*(s+i);
                inb++;
        }
}
*(work+inb)='\0';
for (i=0; i <= inb; i++)
{
        *(s+i)=*(work+i);
}
free(work);
return;
}


i4 move_up(FILE *fp, pos shiftup, pos startloc, pos endloc)
/*


This function shifts the contents of a file upward in the file
by an amount shiftup bytes, starting at position startloc, and ending
at endloc (counted before the shift is done).  endloc is assumed to be the
position at the end of the file before this process starts.
The file (FILE *fp) is assumed to have already been opened.


*/
{
pos bigbuf = 1048576; /* a 1 MB max buffer size, could be adjusted */
pos bufsize;
pos sz_shift;
pos nbuf, lower_loc;
pos fag_end;
char *buf;
pos i;
sz_shift=(endloc-startloc);
nbuf=sz_shift/bigbuf;
fag_end=sz_shift % bigbuf;
bufsize = (nbuf > 0) ? bigbuf : fag_end;
buf = malloc(bufsize);
lower_loc=endloc;
for (i=0;i <= nbuf;i++)
{
        bufsize = (nbuf > i) ? bigbuf : fag_end;
        lower_loc -= bufsize;
        fseeko(fp,lower_loc,SEEK_SET);
/* debug
        printf("lower_loc = %llu\n",lower_loc);
*/
        fread(buf,bufsize,1,fp);
        fseeko(fp,lower_loc+shiftup,SEEK_SET);
        fwrite(buf,bufsize,1,fp);
/* debug 
        location=ftello(fp);
        printf("location = %llu\n",location);
        fflush(stdout);
*/
}


free(buf);
return 0;
}


i4 move_dn(FILE *fp, pos shiftdn, pos startloc, pos endloc)
/*


This function shifts the contents of a file downward in the file
by an amount shiftdn bytes, starting at position startloc, and ending
at endloc (counted before the shift is done).  endloc is assumed to be the
position at the end of the file before this process starts.
The file (FILE *fp) is assumed to have already been opened, and remains open.


*/
{
pos bigbuf = 1048576; /* a 1 MB max buffer size, could be adjusted */


pos bufsize;
pos sz_shift;
pos nbuf, lower_loc;
pos fag_end;
char *buf;
pos i;


/* debug 
printf("shiftdn = %llu, startloc = %llu, endloc = %llu\n",shiftdn,startloc,
          endloc);
*/


sz_shift=(endloc-startloc-shiftdn); /* note this is diff. than in move_up! */
nbuf=sz_shift/bigbuf;
fag_end=sz_shift % bigbuf;
bufsize = (nbuf > 0) ? bigbuf : fag_end;


if(bufsize > 0)
{


   /* debug
   printf("move_dn: bufsize = %llu\n",bufsize);
   */


   buf = malloc(bufsize);
   lower_loc=startloc;
   for (i=0;i <= nbuf;i++)
   {
        bufsize = (nbuf > i) ? bigbuf : fag_end;


   /* debug
        printf("i = %llu\n",i);
        printf("lower_loc = %llu\n",lower_loc);
        printf("bufsize = %llu\n", bufsize);
        printf("shiftdn = %llu\n", shiftdn);
        fflush(stdout);
  */


        fseeko(fp,lower_loc+shiftdn,SEEK_SET);
        fread(buf,bufsize,1,fp);
        fseeko(fp,lower_loc,SEEK_SET);
        fwrite(buf,bufsize,1,fp);
        lower_loc+=bufsize;
   }
   free(buf);


}


/*
If after shifting down, you want to close the file and have it be
shorter, do the following after the call to shift_dn:


fflush(fp);
fn=fileno(fp);
ftruncate(fn, (pos)endloc - (pos)shiftdn )


*/
return 0;
}


i4 sdf_delete(char *fname, i4 idelete)
{


/* Delete dataset no. idelete from the file fname */


FILE *fp;
data_id *id;
pos hdrpos, datapos, hdrtmp, datatmp, datasize, curpos, off, offset;
pos size_del, byte_del=(pos)0, filepos;
char *sdfstr="SDF format";
char *buff;
char testid[300];
i4 ise,ibe, next_order, next_tmp, bufsize, iorder;
pos hdrsize, hdrsztmp;
int fotest;


/* Make certain (pos) is defined as 64bit int, otherwise, exit */


test_sizes();
fotest=0;



ibe = is_big_endian ();
ise = 0;
if (ibe == 0) ise = 1;    /* set flag for byteswapping if small endian */


/* Initially, open the file as read-only */


fp=fopen(fname,"rb");


if (fp == NULL) /* does file fname exist? */
{
        /* if not, print msg and exit */
        printf("sdf_delete: file %s does not exist; exiting\n",fname);
        fflush(stdout);
        exit(1);
}
else
{
        fotest=1;
        fgets(testid,strlen(sdfstr)+1,fp);
        /* Now, test to make sure the SDF file identifier string is there */
        /* If not, print msg and exit */
        if(strncmp(sdfstr,testid,strlen(sdfstr)))
        {
                fclose(fp);
                printf("sdf_delete: File %s is not an SDF file, exiting\n",
                        fname);
                fflush(stdout);
                fotest=0;
                exit(1);
        }
}


if (fotest)
{
        /*close file if it is open */
        fclose(fp);
        fotest=0;
}


/* re-open file with both read/write permissions */


fp=fopen(fname,"rb+");


/* Now, read the locations for next writes of header info and data */


fseeko(fp, (pos) strlen(sdfstr)+1, SEEK_SET);
fread(&hdrpos,sizeof(pos),1,fp);
if(ise) byteswap((void *)&hdrpos,(pos)1,sizeof(pos));


fread(&datapos,sizeof(pos),1,fp);
if(ise) byteswap((void *)&datapos,(pos)1,sizeof(pos));


fread(&next_order,sizeof(i4),1,fp);
if(ise) byteswap((void *)&next_order,(pos)1,sizeof(i4));


fread(&hdrsize,sizeof(pos),1,fp);
if(ise) byteswap((void *)&hdrsize,(pos)1,sizeof(pos));


if((idelete >= next_order) || (idelete < 0))
{
        printf("sdf_delete: can't delete non-existent dataset no. %d in %s\n",
            idelete,fname);
        fflush(stdout);
        fclose(fp);
        exit(1);
}
if(next_order == 1)
{
        fclose(fp);
        fotest=0;
        sdf_rm(fname);
        return 0;
}


bufsize=hdrsize-sizeof(pos)-sizeof(i4); 
buff=(char *)calloc(bufsize*sizeof(char),1);


off=(pos)0;
curpos = (pos) strlen(sdfstr)+1+sizeof(pos)+hdrsize;
for (iorder=0; iorder < next_order; iorder++)
{
        id=(data_id *) string_to_data(fp);
        datasize=data_size(id);
        filepos=ftello(fp);
        if(iorder == idelete)
        {
                size_del = data_size(id);
                byte_del = size_del * id->nbpw;
                if(id->datatype == 'c')
                {
                        byte_del *= (pos)2;
                }
                /* Remove the offending dataset, shorten file, reposition */


                /* debug
                printf("byte_del = %d, curpos = %d, datapos = %d\n",
                    (i4)byte_del,(i4)curpos,(i4)datapos);
                */


                move_dn(fp, byte_del, curpos, datapos);
                fflush(fp);
/*
                fn=fileno(fp);
                ftruncate(fn, (pos) datapos - (pos) byte_del );
*/
                sdf_file_truncate(fp, (pos) datapos - (pos) byte_del );
                fseeko(fp,filepos,SEEK_SET);


        }
        if(iorder >= idelete) /* right limit for this? */
        {
                id->order--;
        }
        if(iorder != idelete)
        {
                data_to_buff(id,buff+off,&offset);
                off+=offset;
        }
        if(id->datatype != 'c') /* test for complex var */
        {
                curpos += (pos) datasize*(id->nbpw);
        }
        else
        {
                curpos += (pos) datasize*2*(id->nbpw);
                /* kludge for complex variables */
        }
        free(id->label);
        free(id->dims);
        free(id);
}
hdrpos = strlen(sdfstr)+1 + 2*sizeof(pos) + sizeof(i4) + sizeof(pos);
fseeko(fp, hdrpos, SEEK_SET);
/* Gonna delete the \0 here, caused a bug */
/* fprintf(fp,"%s%c",buff,'\0'); */
fprintf(fp,"%s",buff);
hdrpos= (pos) ftello(fp);
datapos -= byte_del;
next_order--;


hdrtmp=hdrpos;
datatmp=datapos;
next_tmp=next_order;
hdrsztmp=hdrsize;


if(ise) byteswap((void *)&hdrtmp,(pos)1,sizeof(pos));
if(ise) byteswap((void *)&datatmp,(pos)1,sizeof(pos));
if(ise) byteswap((void *)&next_tmp,(pos)1,sizeof(i4));
if(ise) byteswap((void *)&hdrsztmp,(pos)1,sizeof(pos));


fseeko(fp, (pos) strlen(sdfstr)+1, SEEK_SET);


fwrite(&hdrtmp,sizeof(pos),1,fp);
fwrite(&datatmp,sizeof(pos),1,fp);
fwrite(&next_tmp,sizeof(i4),1,fp);
fwrite(&hdrsztmp,sizeof(pos),1,fp);


fclose(fp);
free(buff);


/* We're done */


return 0;
}


i4 sdf_insert(char *fname, i4 insert, data_id * id_new, void * data_new)
{


/* Insert new dataset after no. "insert" from the file fname */


FILE *fp;
data_id *id;
pos hdrpos, datapos, hdrtmp, datatmp, datasize, curpos, off, offset;
pos size_ins, byte_ins, filepos;
char *sdfstr="SDF format";
char *buff;
char testid[300];
i4 ise,ibe, next_order, next_tmp, bufsize, iorder;
pos hdrinit=HINITSZ, safety=100;
pos hdrsize, oldhdrsize, hdrsztmp, startloc;
int fotest;
/* Make certain (pos) is defined as 64bit int, otherwise, exit */


test_sizes();
fotest=0;



ibe = is_big_endian ();
ise = 0;
if (ibe == 0) ise = 1;    /* set flag for byteswapping if small endian */

/* test for 0 or negative id_new->npbw; exit with error message if so */

if(id_new->nbpw <= 0)
        {
           printf("sdf_insert: id_new->nbpw = %d <= 0; fatal\n",
                (int)(id_new->nbpw));
           fflush(stdout);
           exit(1);
        }


/* Initially, open the file as read-only */


fp=fopen(fname,"rb");


if (fp == NULL) /* does file fname exist? */
{
        /* if not, print msg and exit */
        printf("sdf_insert: file %s does not exist; exiting\n",fname);
        fflush(stdout);
        exit(1);
}
else
{
        fotest=1;
        fgets(testid,strlen(sdfstr)+1,fp);
        /* Now, test to make sure the SDF file identifier string is there */
        /* If not, print msg and exit */
        if(strncmp(sdfstr,testid,strlen(sdfstr)))
        {
                fclose(fp);
                fotest=0;
                printf("sdf_insert: File %s is not an SDF file, exiting\n",
                        fname);
                fflush(stdout);
                exit(1);
        }
}


if (fotest)
{
        /*close file if it is open */
        fclose(fp);
        fotest=0;
}


/* re-open file with both read/write permissions */


fp=fopen(fname,"rb+");


/* Now, read the locations for next writes of header info and data */


fseeko(fp, (pos) strlen(sdfstr)+1, SEEK_SET);
fread(&hdrpos,sizeof(pos),1,fp);
if(ise) byteswap((void *)&hdrpos,(pos)1,sizeof(pos));


fread(&datapos,sizeof(pos),1,fp);
if(ise) byteswap((void *)&datapos,(pos)1,sizeof(pos));


fread(&next_order,sizeof(i4),1,fp);
if(ise) byteswap((void *)&next_order,(pos)1,sizeof(i4));


fread(&hdrsize,sizeof(pos),1,fp);
if(ise) byteswap((void *)&hdrsize,(pos)1,sizeof(pos));


if((insert > next_order) || (insert < 0))
{
        printf("sdf_insert: can't insert before dataset %d in %s\n",
            insert,fname);
        fflush(stdout);
        fclose(fp);
        exit(1);
}


if( hdrpos >= (pos) (hdrsize - safety) )
{
/* now increase size of header area, and move the rest of the file up,
   adjust the current value of datapos accordingly */
 printf("sdf_insert: header area almost full, incrementing by hdrinit.\n");
 fflush(stdout);
 oldhdrsize=hdrsize;
 hdrsize+=hdrinit;
 hdrsztmp=hdrsize;
/*
 fseeko(fp,hdrpos-(pos)(sizeof(pos)),SEEK_SET);
*/
 fseeko(fp,(pos)(strlen(sdfstr)+1+2*sizeof(pos)+sizeof(i4)),SEEK_SET);
 if(ise) byteswap((void *)&hdrsztmp,(pos)1,sizeof(pos));
 fwrite(&hdrsztmp, sizeof(pos),1,fp);
 startloc=(pos) strlen(sdfstr)+1+sizeof(pos)+oldhdrsize;
 move_up(fp,hdrinit,startloc,datapos);
 fflush(fp);
/*
 fn=fileno(fp);
*/
 datapos += hdrinit;
/*
  ftruncate(fn, (pos) datapos);
*/
 sdf_file_truncate(fp, (pos) datapos);
}


bufsize=hdrsize-sizeof(pos)-sizeof(i4); 
buff=(char *)calloc(bufsize*sizeof(char),1);


off=(pos)0;
curpos = (pos) strlen(sdfstr)+1+sizeof(pos)+hdrsize;

/* ensure filepos is at first string descriptor if file resized as above */ 
filepos=strlen(sdfstr)+1 + 2*sizeof(pos) + sizeof(i4) + sizeof(pos);
fseeko(fp,filepos,SEEK_SET);

for (iorder=0; iorder <= next_order; iorder++) 
/* loop must include next_order for case of insertion at eof */
{
        filepos=ftello(fp);
        if(iorder == insert) 
        {
                id_new->order=insert;
                size_ins = data_size(id_new);
                byte_ins = size_ins * id_new->nbpw;
                if(id_new->datatype == 'c')
                {
                        byte_ins *= (pos)2;
                }
                /* Insert the new dataset, lengthen file, reposition */


                /* DEBUG
                printf("byte_ins = %d, curpos = %d, datapos = %d\n",
                    (i4)byte_ins,(i4)curpos,(i4)datapos);
                */
                move_up(fp, byte_ins, curpos, datapos);
                fflush(fp);
/*
                fn=fileno(fp);
                ftruncate(fn, (pos) datapos + (pos) byte_ins );
*/
                sdf_file_truncate(fp, (pos) datapos + (pos) byte_ins );
                fseeko(fp,curpos,SEEK_SET);
                /* test for complex variable type */
                if (id_new->datatype != 'c') 
                {
                        if(ise) byteswap((void *)data_new, 
                              (pos) size_ins, id_new->nbpw);
                        fwrite(data_new,id_new->nbpw,size_ins,fp); 
                        /* is this 64-bit compatible? */
                        /* for non-destructive, should btye-swap back! */
                        if(ise) byteswap((void *)data_new, 
                              (pos) size_ins, id_new->nbpw);
                }
                else
                {
                        if(ise) byteswap((void *)data_new, 
                               (pos) 2*size_ins, id_new->nbpw);
                        fwrite(data_new,id_new->nbpw,(pos)2*size_ins,fp); 
                        /*kludge for complex var */
                        /* for non-destructive, should btye-swap back! */
                        if(ise) byteswap((void *)data_new, 
                               (pos) 2*size_ins, id_new->nbpw);
                }
                curpos=ftello(fp);  /* SHOULD THIS BE HERE? */
                /* insert new string for new dataset */
                fseeko(fp,filepos,SEEK_SET);
                data_to_buff(id_new,buff+off,&offset);
                off+=offset;


        }
        if(iorder < next_order)
        {
                id=(data_id *) string_to_data(fp);
                datasize=data_size(id);
                if(iorder >= insert) /* right limit for this? */
                {
                        id->order++;
                }


                data_to_buff(id,buff+off,&offset);
                off+=offset;
                if(id->datatype != 'c') /* test for complex var */
                {
                        curpos += (pos) datasize*(id->nbpw);
                }
                else
                {
                        curpos += (pos) datasize*2*(id->nbpw);
                        /* kludge for complex variables */
                }


                free(id->label);
                free(id->dims);
                free(id);
        }


        /* make sure you free id_new in calling program */
}


hdrpos = strlen(sdfstr)+1 + 2*sizeof(pos) + sizeof(i4) + sizeof(pos);
fseeko(fp, hdrpos, SEEK_SET);
/* Gonna delete the \0 here */
/* fprintf(fp,"%s%c",buff,'\0'); */
fprintf(fp,"%s",buff);
hdrpos= (pos) ftello(fp);
datapos = curpos;
next_order++;


hdrtmp=hdrpos;
datatmp=datapos;
next_tmp=next_order;
hdrsztmp=hdrsize;


if(ise) byteswap((void *)&hdrtmp,(pos)1,sizeof(pos));
if(ise) byteswap((void *)&datatmp,(pos)1,sizeof(pos));
if(ise) byteswap((void *)&next_tmp,(pos)1,sizeof(i4));
if(ise) byteswap((void *)&hdrsztmp,(pos)1,sizeof(pos));


fseeko(fp, (pos) strlen(sdfstr)+1, SEEK_SET);


fwrite(&hdrtmp,sizeof(pos),1,fp);
fwrite(&datatmp,sizeof(pos),1,fp);
fwrite(&next_tmp,sizeof(i4),1,fp);
fwrite(&hdrsztmp,sizeof(pos),1,fp);


fclose(fp);
free(buff);


/* We're done */


return 0;
}


i4 sdf_replace(char *fname, i4 replace, data_id * id_new, void * data_new)
{


/* Replace dataset "replace" with descriptor id_new and data_new in file fname*/


FILE *fp;
data_id *id;
pos hdrpos, datapos, hdrtmp, datatmp, datasize, curpos, off, offset;
pos size_new, byte_new, byte_old, shiftup, shiftdn, filepos;
char *sdfstr="SDF format";
char *buff;
char testid[300];
i4 ise,ibe, next_order, next_tmp, bufsize, iorder;
pos hdrinit=HINITSZ, safety=100;
pos hdrsize, oldhdrsize, hdrsztmp, startloc;
int fotest;


/* Make certain (pos) is defined as 64bit int, otherwise, exit */


test_sizes();
fotest=0;


ibe = is_big_endian ();
ise = 0;
if (ibe == 0) ise = 1;    /* set flag for byteswapping if small endian */

/* test for 0 or negative id_new->npbw; exit with error message if so */

if(id_new->nbpw <= 0)
        {
           printf("sdf_replace: id_new->nbpw = %d <= 0; fatal\n",
                (int)(id_new->nbpw));
           fflush(stdout);
           exit(1);
        }


/* Initially, open the file as read-only */


fp=fopen(fname,"rb");


if (fp == NULL) /* does file fname exist? */
{
        /* if not, print msg and exit */
        printf("sdf_replace: file %s does not exist; exiting\n",fname);
        fflush(stdout);
        exit(1);
}
else
{
        fotest=1;
        fgets(testid,strlen(sdfstr)+1,fp);
        /* Now, test to make sure the SDF file identifier string is there */
        /* If not, print msg and exit */
        if(strncmp(sdfstr,testid,strlen(sdfstr)))
        {
                fclose(fp);
                fotest=0;
                printf("sdf_replace: File %s is not an SDF file, exiting\n",
                        fname);
                fflush(stdout);
                exit(1);
        }
}


if (fotest)
{
        /*close file if it is open */
        fclose(fp);
        fotest=0;
}


/* re-open file with both read/write permissions */


fp=fopen(fname,"rb+");


/* Now, read the locations for next writes of header info and data */


fseeko(fp, (pos) strlen(sdfstr)+1, SEEK_SET);
fread(&hdrpos,sizeof(pos),1,fp);
if(ise) byteswap((void *)&hdrpos,(pos)1,sizeof(pos));


fread(&datapos,sizeof(pos),1,fp);
if(ise) byteswap((void *)&datapos,(pos)1,sizeof(pos));


fread(&next_order,sizeof(i4),1,fp);
if(ise) byteswap((void *)&next_order,(pos)1,sizeof(i4));


fread(&hdrsize,sizeof(pos),1,fp);
if(ise) byteswap((void *)&hdrsize,(pos)1,sizeof(pos));


if(( replace >= next_order) || (replace < 0))
{
        printf("sdf_replace: can't replace dataset %d in %s\n",
            replace,fname);
        fflush(stdout);
        fclose(fp);
        exit(1);
}


if( hdrpos >= (pos) (hdrsize - safety) )
{
/* now increase size of header area, and move the rest of the file up,
   adjust the current value of datapos accordingly */
 printf("sdf_replace: header area almost full, incrementing by hdrinit.\n");
 fflush(stdout);
 oldhdrsize=hdrsize;
 hdrsize+=hdrinit;
 hdrsztmp=hdrsize;
/*
 fseeko(fp,hdrpos-(pos)(sizeof(pos)),SEEK_SET);
*/
 fseeko(fp,(pos)(strlen(sdfstr)+1+2*sizeof(pos)+sizeof(i4)),SEEK_SET);
 if(ise) byteswap((void *)&hdrsztmp,(pos)1,sizeof(pos));
 fwrite(&hdrsztmp, sizeof(pos),1,fp);
 startloc=(pos) strlen(sdfstr)+1+sizeof(pos)+oldhdrsize;
 move_up(fp,hdrinit,startloc,datapos);
 fflush(fp);
/*
 fn=fileno(fp);
*/
 datapos += hdrinit;
/*
  ftruncate(fn, (pos) datapos);
*/
 sdf_file_truncate(fp, (pos) datapos);
}


bufsize=hdrsize-sizeof(pos)-sizeof(i4); 
buff=(char *)calloc(bufsize*sizeof(char),1);


off=(pos)0;
curpos = (pos) strlen(sdfstr)+1+sizeof(pos)+hdrsize;
/* ensure filepos is at first string descriptor if file resized as above */ 
filepos=strlen(sdfstr)+1 + 2*sizeof(pos) + sizeof(i4) + sizeof(pos);
fseeko(fp,filepos,SEEK_SET);

for (iorder=0; iorder < next_order; iorder++) 
{
        id=(data_id *) string_to_data(fp);
        filepos=ftello(fp);
        datasize=data_size(id); /* size of old dataset */
        if(iorder == replace) 
        {
                id_new->order=replace;
                size_new = data_size(id_new); /* size of new dataset */
                byte_new = size_new * id_new->nbpw;
                if(id_new->datatype == 'c')
                {
                        byte_new *= (pos)2;
                }
                byte_old = datasize*id->nbpw;
                if(id->datatype == 'c')
                {
                        byte_old*= (pos)2;
                }
                shiftup=byte_new - byte_old;
                if(shiftup < 0) 
                {
                        /* if true, shrink file size */
                        shiftdn = -shiftup;
                        move_dn(fp,shiftdn,curpos,datapos);
                        fflush(fp);
/*
                        fn=fileno(fp);
                        ftruncate(fn,(pos)datapos - (pos)shiftdn );
*/
                        sdf_file_truncate(fp, (pos) datapos - (pos) shiftdn );
                }
                else if (shiftup > 0)
                {
                        /* if true, expand file size */
                        move_up(fp,shiftup,curpos,datapos);
                        fflush(fp);
/*
                        fn=fileno(fp);
                        ftruncate(fn,(pos)datapos + (pos)shiftup );
*/
                        sdf_file_truncate(fp, (pos) datapos + (pos) shiftup );
                }


                /* now must actually write the new data to the file */


                /* DEBUG
                printf("byte_ins = %d, curpos = %d, datapos = %d\n",
                    (i4)byte_ins,(i4)curpos,(i4)datapos);
                */


                fseeko(fp,curpos,SEEK_SET); /* position file for write */


                /* test for complex variable type */
                if (id_new->datatype != 'c') 
                {
                        if(ise) byteswap((void *)data_new, 
                              (pos) size_new, id_new->nbpw);
                        fwrite(data_new,id_new->nbpw,size_new,fp); 
                        /* is this 64-bit compatible? */
                        /* for non-destructive, should btye-swap back! */
                        if(ise) byteswap((void *)data_new, 
                              (pos) size_new, id_new->nbpw);
                }
                else
                {
                        if(ise) byteswap((void *)data_new, 
                               (pos) 2*size_new, id_new->nbpw);
                        fwrite(data_new,id_new->nbpw,(pos)2*size_new,fp); 
                        /*kludge for complex var */
                        /* for non-destructive, should btye-swap back! */
                        if(ise) byteswap((void *)data_new, 
                               (pos) 2*size_new, id_new->nbpw);
                }
                curpos=ftello(fp);  /* SHOULD THIS BE HERE? */
                /* create new string for new dataset */
                fseeko(fp,filepos,SEEK_SET); 
                data_to_buff(id_new,buff+off,&offset);
                off+=offset;


        }


        else


        {


                /* do this stuff only if dataset is not being replaced */
                data_to_buff(id,buff+off,&offset);
                off+=offset;
                if(id->datatype != 'c') /* test for complex var */
                {
                        curpos += (pos) datasize*(id->nbpw);
                }
                else
                {
                        curpos += (pos) datasize*2*(id->nbpw);
                        /* kludge for complex variables */
                }
                /* end of stuff to be done if datset is not being replaced */


        }


        free(id->label);
        free(id->dims);
        free(id);


}
        /* remember to free id_new in calling program */


hdrpos = strlen(sdfstr)+1 + 2*sizeof(pos) + sizeof(i4) + sizeof(pos);
fseeko(fp, hdrpos, SEEK_SET);
/*  Gonna change this to get rid of \0 at end */
/* fprintf(fp,"%s%c",buff,'\0'); */
fprintf(fp,"%s",buff);
hdrpos= (pos) ftello(fp);
datapos = curpos;


hdrtmp=hdrpos;
datatmp=datapos;
next_tmp=next_order;
hdrsztmp=hdrsize;


if(ise) byteswap((void *)&hdrtmp,(pos)1,sizeof(pos));
if(ise) byteswap((void *)&datatmp,(pos)1,sizeof(pos));
if(ise) byteswap((void *)&next_tmp,(pos)1,sizeof(i4));
if(ise) byteswap((void *)&hdrsztmp,(pos)1,sizeof(pos));


fseeko(fp, (pos) strlen(sdfstr)+1, SEEK_SET);


fwrite(&hdrtmp,sizeof(pos),1,fp);
fwrite(&datatmp,sizeof(pos),1,fp);
fwrite(&next_tmp,sizeof(i4),1,fp);
fwrite(&hdrsztmp,sizeof(pos),1,fp);


fclose(fp);
free(buff);


/* We're done */


return 0;
}


void sdf_transpose_f77__(intf *indorderf77, intf *directionsf77, char * dtf77,
     intf *nbpwf77, intf *ndimf77, intf8* dimsf77, void *dataf77, 
     intf_clen len_dtf77)
{
sdf_transpose_f77(indorderf77, directionsf77, dtf77, nbpwf77, ndimf77, dimsf77,
        dataf77, len_dtf77);
return;
}


void sdf_transpose_f77_(intf *indorderf77, intf *directionsf77, char * dtf77,
     intf *nbpwf77, intf *ndimf77, intf8* dimsf77, void *dataf77, 
     intf_clen len_dtf77)
{
sdf_transpose_f77(indorderf77, directionsf77, dtf77, nbpwf77, ndimf77, dimsf77,
        dataf77, len_dtf77);
return;
}
void sdf_transpose_f77(intf *indorderf77, intf *directionsf77, char * dtf77,
     intf *nbpwf77, intf *ndimf77, intf8* dimsf77, void *dataf77, 
     intf_clen len_dtf77)
{
char na[3]="NA";
data_id *id;
i4 *indorder;
i4 *directions;
i4 i,nlab;
id=malloc(sizeof(data_id));
id->ndim=(i4) *ndimf77;
id->nbpw=(i4) *nbpwf77;
id->datatype= (char) *dtf77;
nlab=strlen(na);
id->label=calloc(nlab+1,1);
strncpy(id->label,na,2);
id->dims=malloc(id->ndim*sizeof(pos));
indorder=malloc(id->ndim*sizeof(i4));
directions=malloc(id->ndim*sizeof(i4));


for (i=0;i<id->ndim;i++)
{
        *(id->dims+i)= (pos) *(dimsf77+id->ndim-1-i); /* flip dims for f77 */


        /*indorder starts at 1 in f77  - must flip order for f77. 
         * Note also that definitions of indorder for C routine
         * are ndim-1- defs in fortran caller, +1 (because indorder starts at
         * 1 in fortran version).
         * I know this all seems baffling, but it gives correct results*/


        *(indorder+i)= (i4) id->ndim-1 - *(indorderf77+id->ndim-1-i)+1; 
        *(directions+i)= (i4) *(directionsf77+id->ndim-1-i); /*flip order f77 */
}


sdf_transpose(indorder, directions, id, dataf77);


/* must also change the dimsf77 array on output, since id->dims gets changed */
for (i=0; i < id->ndim; i++)
{
        *(dimsf77+i) = (intf8) *(id->dims+id->ndim-1 -i); /* flip dims for f77 */
}


free(id->label);
free(id->dims);
free(id);
free(directions);
free(indorder);
return;
}


i4 sdf_transpose(i4 * indorder, i4 *directions, data_id *id, void *data)
/*


Use vacancy cycle tracking technique to convert multi-d array from original
index and dimension order to a permuted order, determined by the values 
pointed to by indorder.  The corresponding values of the directions
array determines whether or not the order of the new indices occurs
in their original order, or are reversed, depending on the sign of
the corresponding directions element.  After the function
is called, data will have been shuffled in place, and the order of id->dims
will be changed, to be consistent with indorder and the
sdf notation for array dimensions.
The vacancy cycle tracking technique of Ding (2001) 
[Chris H. Q. Ding, "An Optimal Index Reshuffle Algorithm for Multidimensional
Arrays and its Applications for Parallel Architectures", IEEE Transactions 
on Parallel and Distributed Systems, vol. 12, No. 3, pp 306-315, 2001]
allows the transformation
to be done in place, rather than copying from one array to a new array.


The speed of the reshuffle part of this procedure could be greatly speeded up
if we knew that the last index, or last set of indices, were the same as the
original.  In that case, the overall loop indexed by ir would be much smaller,
and the size of the memory transfers (the size of tmp) would be correspondingly
bigger.  The factor in each case would be the value of dim[last] (or for
several (q) dimensions together), the factor would be dim[last-q-1] x ... x
dim[last-2]*dim[last-1] - not sure that last bit is correct, but you get
the idea.  Roughing this out, which would involve testing indorder against
[0,...n-1] (ie newindorder, assuming it has no errors)
and then seeing how many contiguous last elements there were.
The dimension of the shuffle would also be reduced, so the values of ndim
passed into rindcalc and memcalcr would have to be changed as well.


Index reversal is done as a second step after
the re-shuffling has been accomplished.
It is possible the two steps could be combined into one, but that is too
complicated for me to contemplate.


*/
{
pos datasize,ir;
pos *indices_r, *indices;
/* pos memloc_c; */   /* not used here anymore */ 
/* pos memloc_r; */   /* not used here anymore */
pos istart,ilast,inext;
pos nbytes, nelem_cycle;
pos *newdims;
i4 *newindorder;
/* i4 *newdirections; */
u4 *mask;
i4 i,reverse,badindex,trivial;
char *tmp;
test_sizes();
/* define the original order of dimensions and indices: */
datasize=data_size(id);
nbytes=(pos)id->nbpw;
if(id->datatype == 'c') nbytes*= (pos)2;
tmp=malloc(nbytes);
indices_r=malloc(id->ndim*sizeof(pos));
indices=malloc(id->ndim*sizeof(pos));
newindorder=malloc(id->ndim*sizeof(i4));
/* newdirections=malloc(id->ndim*sizeof(i4)); */ /*shuffled directions array */
/* make sure the shuffled dimensions are correctly specified: */
for (i=0;i<id->ndim;i++)
{
        *(newindorder+i)=*(indorder+i);

                        
}
/* sort copy of indorder array */
qsort((void *)newindorder,id->ndim,sizeof(i4),i4cmp); 
badindex=0;
trivial=1;
for (i=0;i<id->ndim;i++)
{
        if(*(newindorder+i) != (i4)i) 
        {
                badindex=1;
                printf("sdf_transpose: newindorder[%d] = %d\n",
                        i,*(newindorder+i));
                fflush(stdout);
        }
        if(*(indorder+i) != (i4)i) trivial=0;
}


/* if the shuffled order of indorder and the original dim orders don't match,
   we're screwed, so bail and return without doing anything */


if(badindex)
{
        printf("sdf_transpose:  error in indorder, no transpose done\n");
        fflush(stdout);
        free(newindorder);
        /* free(newdirections); */
        free(tmp);
        free(indices_r);
        free(indices);
        return 1;
}


if(trivial == 0)


/* START OF RESHUFFLE: */


{
        newdims=malloc(id->ndim*sizeof(pos));
        mask=(u4 *)calloc(datasize/(sizeof(u4)*(pos)8)+(pos)1,sizeof(u4));
/* First, get array of new dimensions in the desired order: */
        for (i=0;i<id->ndim;i++)
        {
                newdims[i]=id->dims[indorder[i]];
        /* newdirections array is not currently used for anything 
           But I am leaving this in here as a reminder for future work */
        /* newdirections[i]=directions[indorder[i]]; */
        }
/* Loop over each gridpoint, doing the vacancy cycle if not done already */
        for (ir=0;ir<datasize;ir++)
        {
                while(readmask(mask,datasize,ir) == 0)
                {
                        istart=ir;
                        ilast=istart;
                        rindcalc(id->ndim,newdims,indices_r,ir); 
                        for (i=0;i<id->ndim;i++)
                        {
                                indices[indorder[i]]=indices_r[i];
                        }
                        memcalcr(id->dims,indices,id->ndim,&inext);
                        memcpy((void *)tmp,(void *)((char *)data+inext*nbytes),
                               nbytes);
                        nelem_cycle=(pos)0;
                        writemask(mask,datasize,istart,1);
                        while(inext != istart)
                        {
                                ilast=inext;
                                rindcalc(id->ndim,newdims,indices_r,ilast); 
                                for (i=0; i< id->ndim; i++)
                                {
                                        indices[indorder[i]]=indices_r[i];
                                }
                                memcalcr(id->dims,indices,id->ndim,&inext);
                                memcpy((void *)((char *)data+ilast*nbytes),
                                         (void *)((char *)data+inext*nbytes), 
                                           nbytes);
                                nelem_cycle++;
                                writemask(mask,datasize,ilast,1);
                        }
                        /*
                         printf("nelem_cycle = %d\n",(i4)nelem_cycle);
                        */
                        memcpy((void *)((char *)data+istart*nbytes),
                                (void *)tmp,nbytes);
                }
        }
        
/* flip order of id->dims */
        for (i=0; i < id->ndim; i++)
        {
                *(id->dims+i)=(pos) *(newdims+i);
        }
        free(newdims);
        free(mask);
/*
END OF RESHUFFLE
*/
}


reverse=0;
for (i=0;i<id->ndim;i++)
{
        if(directions[i] < 0) reverse=1; 
        /* if(newdirections[i] < 0) reverse=1; */
}
/*
Next block of code does the index reversal, if any of the elements of
the directions array are negative:


BEGIN INDEX REVERSAL
*/


if(reverse)
{
        /* re-allocate and reset bit mask: */
        mask=(u4 *)calloc(datasize/(sizeof(u4)*(pos)8)+(pos)1,sizeof(u4));
        for (ir=0;ir<datasize;ir++)
        {
                while(readmask(mask,datasize,ir) == 0)
                {
                        rindcalc(id->ndim,id->dims,indices_r,ir);
                        for (i=0;i < id->ndim; i++)
                        {
                                indices[i] = (directions[i] < 0) 
                                ? id->dims[i]-(pos)1-(pos)indices_r[i] 
                                : indices_r[i]; 
                        }
                        memcpy((void *)tmp,(void *)((char *)data+ir*nbytes),
                                nbytes);
                        memcalcr(id->dims,indices,id->ndim,&inext);
                        memcpy((void *)((char *)data+ir*nbytes),
                                (void *)((char *)data+inext*nbytes), nbytes);
                        memcpy((void *)((char *)data+inext*nbytes),
                                (void *)tmp,nbytes);
                        writemask(mask,datasize,ir,1);
                        writemask(mask,datasize,inext,1);
                }
        }
        free(mask);
}
/* END INDEX REVERSAL */


/* clean up and return */
free(newindorder);
/* free(newdirections); */
free(tmp);
free(indices);
free(indices_r);
return 0;
}

void rindcalc_f77__(intf *ndimf77, intf8 *dimsf77, intf8 *indf77, 
      intf8 * memloc_rf77)
{
rindcalc_f77(ndimf77, dimsf77, indf77, memloc_rf77);
return;
}

void rindcalc_f77_(intf *ndimf77, intf8 *dimsf77, intf8 *indf77, 
      intf8 * memloc_rf77)
{
rindcalc_f77(ndimf77, dimsf77, indf77, memloc_rf77);
return;
}

void rindcalc_f77(intf *ndimf77, intf8 *dimsf77, intf8 *indf77, 
      intf8 * memloc_rf77)
{
/* Fortran callable version of rindcalc.  Returns an array of indices into
indf77, given memory location *memloc_rf77 and dimensions dims. 
Note that rindcalc is not really applicable to fortran arrays, should use
cindcalc_f77 (column major indexing) for use with fortran arrays.  */

intf8 bigq;
intf8 temp;
intf i;
bigq = *memloc_rf77 - (intf8)1; /* subtract 1 from fortran 1-d index */
for (i=*ndimf77-1;i>=0;i--)
{
        temp = bigq / *(dimsf77+i);
        /* add one to C-style index: */
        *(indf77+i) = (bigq - temp* (*(dimsf77+i))) + (intf8)1;
        bigq = temp;
}
return;
}

i4 rindcalc(i4 ndim, pos *dims, pos *indices, pos memloc_r)
{
/*  This function returns an array of indices corresponding to the memory
    location memloc_r (assuming row major) and dimensions dims.
    Assumes memory for pointers dims and indices are allocated in the
    caller.
*/
pos bigq;
pos temp;
i4 i;
bigq = memloc_r;
for (i=ndim-1;i>=0;i--)
{
/* new code (somewhat faster, no modulus operation): */
        temp = bigq / *(dims+i);
        *(indices+i) = bigq - temp* (*(dims+i));
        bigq = temp;

/* old code:
        *(indices+i) = bigq % *(dims+i);
        bigq /= *(dims+i);
*/

}
return 0;
}

void cindcalc_f77__(intf *ndimf77, intf8 *dimsf77, intf8 *indf77,
      intf8 * memloc_cf77)
{
cindcalc_f77(ndimf77, dimsf77, indf77, memloc_cf77);
return;
}

void cindcalc_f77_(intf *ndimf77, intf8 *dimsf77, intf8 *indf77,
      intf8 * memloc_cf77)
{
cindcalc_f77(ndimf77, dimsf77, indf77, memloc_cf77);
return;
}

void cindcalc_f77(intf *ndimf77, intf8 *dimsf77, intf8 *indf77,
      intf8 * memloc_cf77)
{
/* Code to convert from a 1-d index to a multi-dimensional index, using
   column major indexing.  This is the version to use from fortran.  Assumes
   array indices go from 1 to dims(i), where i is the ith index. */

intf8 bigm;
intf8 temp;
intf i;
bigm = *memloc_cf77-(intf8)1;
for (i=0;i<*ndimf77;i++)
{
        temp = bigm / *(dimsf77+i);
        *(indf77+i) = (bigm - temp* (*(dimsf77+i))) + (intf8)1;
        bigm = temp;
}
return;
}

i4 cindcalc(i4 ndim, pos *dims, pos *indices, pos memloc_c)
{
/*  This function returns an array of indices corresponding to the memory
    location memloc_c (assuming column major) and dimensions dims.
    Assumes memory for pointers dims and indices are allocated in the
    caller.
*/
pos bigm;
pos temp;
i4 i;
bigm = memloc_c;
for (i=0;i<ndim;i++)
{
/* new code: */
        temp = bigm / *(dims+i);
        *(indices+i) = bigm - temp* (*(dims+i));
        bigm = temp;

/* old code:
        *(indices+i) = bigm % *(dims+i);
        bigm /= *(dims+i);
*/

}
return 0;
}

void memcalcr_f77__(intf8 *dimsf77, intf8 *indf77, intf *ndimf77, 
      intf8 *memloc_rf77)
{
memcalcr_f77(dimsf77,indf77,ndimf77,memloc_rf77);
return;
}

void memcalcr_f77_(intf8 *dimsf77, intf8 *indf77, intf *ndimf77, 
      intf8 *memloc_rf77)
{
memcalcr_f77(dimsf77,indf77,ndimf77,memloc_rf77);
return;
}

void memcalcr_f77(intf8 *dimsf77, intf8 *indf77, intf *ndimf77, 
      intf8 *memloc_rf77)
{
/* Fortran callable version of memcalcr.  Generally, this is not applicable
   to fortran arrays, will probably want to use memcalcc_f77. */

/* For sake of efficiency, make the (dangerous) assumption that intf8 and
   pos integer arrays are equivalent */
pos bigq; 
/*
pos bigm;
*/
intf i;
bigq =(intf8)0;
for (i=0;i<*ndimf77;i++)
{
        bigq=*(dimsf77+i)*bigq+(*(indf77+i) - (intf8)1 );
}
*memloc_rf77=bigq+(intf8)1;
return;
}

i4 memcalcr(pos *dims, pos *indices, i4 ndim, pos *memloc_r)
{
/*
Given dimensions dims and indices indices, allocated in caller, this function
calculates the row major equivalent 1-d location in memory
(memloc_r) for a multidimensional array element.  Code to compute the
column major location memloc_c as well is here, but commented out to speed
things up.
*/
pos bigq; 
/*
pos bigm;
*/
i4 i;
bigq =(pos)0;
/*
 bigm =(pos)0;
*/
for (i=0;i<ndim;i++)
{
        bigq=*(dims+i)*bigq+*(indices+i);
        /*
         bigm=*(dims+ndim-1-i)*bigm+*(indices+ndim-1-i);
        */
}
/*
 *memloc_c=bigm;
*/
*memloc_r=bigq;
return 0;
}

void memcalcc_f77__(intf8 *dimsf77, intf8 *indf77, intf *ndimf77,
      intf8 *memloc_cf77)
{
memcalcc_f77(dimsf77, indf77, ndimf77, memloc_cf77);
return;
}

void memcalcc_f77_(intf8 *dimsf77, intf8 *indf77, intf *ndimf77,
      intf8 *memloc_cf77)
{
memcalcc_f77(dimsf77, indf77, ndimf77, memloc_cf77);
return;
}

void memcalcc_f77(intf8 *dimsf77, intf8 *indf77, intf *ndimf77,
      intf8 *memloc_cf77)
{
/* Fortran callable version of memcalcc.  Converts multi-dimensional index
   into an equivalent 1-d array index.  This is the function relevant for
   fortran arrays.  Must adjust indices, 1-d array locations for fortran
   index notation, c.f. C notation in memcalcc. */

intf8 bigm;
intf i;
bigm =(intf8)0;
for (i=0;i<*ndimf77;i++)
{
        bigm=*(dimsf77+*ndimf77-1-i)*bigm+(*(indf77+*ndimf77-1-i) - (intf8)1);
}
*memloc_cf77= (intf8) bigm + (intf8)1;
return;

}

i4 memcalcc(pos *dims, pos *indices, i4 ndim, pos *memloc_c)
{
/*
Given dimensions dims and indices indices, allocated in caller, this function
calculates the column major equivalent 1-d location in memory
(memloc_c) for a multidimensional array element.  Row-major analogues are
included, but commented out.
*/

/*
pos bigq; 
*/
pos bigm;
i4 i;
/*
bigq =(pos)0;
*/
bigm =(pos)0;
for (i=0;i<ndim;i++)
{
        /*
        bigq=*(dims+i)*bigq+*(indices+i);
        */
        bigm=*(dims+ndim-1-i)*bigm+*(indices+ndim-1-i);
}
*memloc_c=bigm;
/*
*memloc_r=bigq;
*/
return 0;
}


i4 readmask(u4 * mask, pos sizem, pos loc)


/* reads the single bit corresponding to location loc for a mask of
   size sizem.  Assumes user has allocated enough memory for mask in caller,
   and which is assumed to be divided into an array of 32 bit unsigned integers.
   Here is an example:


   mask=calloc(sizem/(sizeof(u4)*(pos)8)+(pos)1,sizeof(u4));


   */
{
/* i4 result; */
pos intindex;
u4 bitindex, baseint;
if((loc >= sizem) || (loc < (pos)0))
{
   printf("readmask: loc is out of range\n");
   exit(1);
}
/*
intindex = (pos) (loc >> 5);
bitindex = (u4) (loc - (intindex << 5));
*/
intindex = (pos)(loc / 32);
bitindex = (u4)(loc - intindex * 32);
/* intindex= (pos)(loc / (sizeof(u4)*(pos)8)); */
/* bitindex= (u4)(loc % (sizeof(u4)*(pos)8));  */
/* bitindex= (u4)(loc - intindex*(sizeof(u4)*(pos)8)); */
baseint=*(mask+intindex);
return ((baseint & ((u4)1 << bitindex)) ? 1 : 0);
/* return result; */
}

i4 writemask(u4 *mask, pos sizem, pos loc, i4 value)
{
/* Write out a single bit according to value; if value=0, a 0 will be written,
   for any other value, a 1 will be written.  This is done for a bit at
   location loc, corresponding to a total bit mask length of sizem.  
   It is assumed
   the user has allocated enough space for mask in caller,
   which is assumed to be an array of unsigned 32 bit integers 
   Here is an example:


   mask=calloc(sizem/(sizeof(u4)*(pos)8)+(pos)1,sizeof(u4));


   */



pos intindex;
/* i4 result; */
u4 shiftone;
u4 baseint,bitindex;
if((loc >= sizem) || (loc < (pos)0))
{
   printf("writemask: loc is out of range\n");
   exit(1);
}
intindex = (pos)(loc / 32);
bitindex = (u4)(loc - intindex * 32);

/* intindex = (pos) (loc >> 5); */
/* bitindex = (u4) (loc - (intindex << 5)); */

/* intindex= (pos)(loc / (sizeof(u4)*(pos)8)); */
/* bitindex= (u4)(loc % (sizeof(u4)*(pos)8)); */

baseint=*(mask+intindex);
shiftone = ((u4)1 << bitindex);
*(mask+intindex) = (value ? (baseint | (shiftone)) : 
        (baseint & ~shiftone)); /* is there a simpler way? */
/* *(mask+intindex)=result; */
return 0;

}

i4 byteswap (unsigned char *arr, pos arrsize, i4 nbpw)
/* Pretty simple:  arr is input array, which is byte-swapped in place,
   nbpw is the number of bytes per word, and arrsize is the size of the array
   (in units of nbpw bytes).  It is assumed that arr has
   already have been correctly defined and allocated in the calling program. 

   This version of byteswap was written by Jack Vernetti at SSL, UC Berkeley 
   It is much faster than the original version I wrote. */
{
    unsigned char tmpb ;
    unsigned char *next_bytes, *bptr, *optr, *lastp ;

    if (nbpw < 2)
    {
       return (0) ; /* no point in byteswapping an array of bytes */
    }

/* Initialize pointers */
    lastp = arr + (nbpw * arrsize) ;
    next_bytes = arr ;
    optr = next_bytes ;
    bptr = optr + nbpw - 1 ;
    next_bytes += nbpw ;
/* Main loop: */
    while (optr < lastp)
    {
/* swap byte values with those of the conjugate location */
	tmpb = *optr ;
	*optr = *bptr ;
	*bptr = tmpb ;
/* increment and decrement pointers to next bytes for swapping */
	optr++ ;
	bptr-- ;
	if (bptr <= optr)
/* If all the bytes have been swapped for this element, move to next element */
        {
	    optr = next_bytes ;
	    bptr = optr + nbpw - 1 ;
	    next_bytes += nbpw;
        }
    }

/* we're done */
    return (0) ;

}

i4 is_big_endian (void)
/* This function returns 1 if it is a big endian machine, 0 if small endian */

/* NOTE:  to make the test below bulletproof, should probably put fakeword
   and realword into a union.  But so far, has worked on all platforms as is. 
   Therefore, will not fix it till I find out it's broken somewhere */
{
  const unsigned char fakeword[4] = { 0xFF, 0x00, 0xFF, 0x00 };
  i4 *realword;
  realword=(i4 *)fakeword;
  if (*realword == 0xFF00FF00)
    {
      return (1);
    }
  else if (*realword == 0x00FF00FF)
    {
      return (0);
    }
  else
    {
      printf("is_big_endian: fatal - value of *realword = 0x%08x\n",*realword);
      fflush(stdout);
      exit(1);
    }
  return (-1); /* will never get here.  Just makes SGI origin compiler happy */
}

void ***** sdf_1d_to_5d(data_id *id, void *data)
{


/* This function returns a 5d dynamically allocated array, based on the
   values of nbpw, ndim, and dims in the id structure.  The values of the
   array are set by the 1-d array data read in from sdf_read */


        char *****data5d;
        i4 sizeof_elem;
        i4 sizeof_ptr;
        pos i,j,k,l;
        sizeof_elem=id->nbpw; /* compute size of data element */
        if(id->datatype == 'c')
        {
                sizeof_elem*=2;
        }
        sizeof_ptr=sizeof(char *); /* compute size of pointer  */
        test_sizes();
        if(id->ndim != 5)
        {
                printf("sdf_mk_5d: id->ndim != 5; returning NULL\n");
                fflush(stdout);
                return NULL;
        }
        data5d=sdf_malloc(id->dims[0]*sizeof_ptr);
        for (i=0;i<id->dims[0];i++)
        {
                data5d[i]=malloc(id->dims[1]*sizeof_ptr);


                for (j=0;j<id->dims[1];j++)
                {
                        data5d[i][j]=malloc(id->dims[2]*sizeof_ptr);
                        for (k=0;k<id->dims[2];k++)
                        {
                                data5d[i][j][k]=malloc(id->dims[3]*sizeof_ptr);
                        }
                }
        }


        data5d[0][0][0][0]=data;


        /* compute addresses for each start of each column */
        for (i=1;i<id->dims[0];i++)
        {
                data5d[i][0][0][0]=data5d[i-1][0][0][0]+id->dims[4]*id->dims[3]*
                        id->dims[2]*id->dims[1]*sizeof_elem;
        }


        for (i=0;i<id->dims[0];i++)
        {
                for (j=1;j<id->dims[1];j++)
                {
                        data5d[i][j][0][0]=data5d[i][j-1][0][0]+id->dims[4]*
                                id->dims[3]*id->dims[2]*sizeof_elem;
                }
        }
        
        for (i=0;i<id->dims[0];i++)
        {
                for (j=0;j<id->dims[1];j++)
                {
                        for (k=1;k<id->dims[2];k++)
                        {
                                data5d[i][j][k][0]=data5d[i][j][k-1][0] +
                                        id->dims[4]*id->dims[3]*sizeof_elem;
                        }
                }
        }
        
        for (i=0;i<id->dims[0];i++)
        {
                for (j=0;j<id->dims[1];j++)
                {
                        for (k=0;k<id->dims[3];k++)
                        {
                                for (l=1;l<id->dims[4];l++)
                                {
                                        data5d[i][j][k][l]=data5d[i][j][k][l-1]
                                        + id->dims[4]*sizeof_elem;
                                }
                        }
                }
        }


        return (void *****)data5d;
}
void **** sdf_1d_to_4d(data_id *id, void *data)
{


/* This function returns a 4d dynamically allocated array, using the
   array values from the 1d array data returned e.g. from sdf_read */


        char ****data4d;
        i4 sizeof_elem;
        i4 sizeof_ptr;
        pos i,j,k;
        sizeof_elem=id->nbpw; /* compute size of data element */
        if(id->datatype == 'c')
        {
                sizeof_elem*=2;
        }
        sizeof_ptr=sizeof(char *); /* compute size of pointer  */
        test_sizes();
        if(id->ndim != 4)
        {
                printf("sdf_1d_to_4d: id->ndim != 4; returning NULL\n");
                fflush(stdout);
                return NULL;
        }
        data4d=sdf_malloc(id->dims[0]*sizeof_ptr);
        for (i=0;i<id->dims[0];i++)
        {
                data4d[i]=malloc(id->dims[1]*sizeof_ptr);


                for (j=0;j<id->dims[1];j++)
                {
                        data4d[i][j]=malloc(id->dims[2]*sizeof_ptr);
                }
        }


        data4d[0][0][0]=data;


        /* compute addresses for each start of each column */
        for (i=1;i<id->dims[0];i++)
        {
                data4d[i][0][0]=data4d[i-1][0][0]+id->dims[3]*id->dims[2]*
                        id->dims[1]*sizeof_elem;
        }
        for (i=0;i<id->dims[0];i++)
        {
                for (j=1;j<id->dims[1];j++)
                {
                        data4d[i][j][0]=data4d[i][j-1][0]+id->dims[3]*
                                id->dims[2]*sizeof_elem;
                }
        }
        
        for (i=0;i<id->dims[0];i++)
        {
                for (j=0;j<id->dims[1];j++)
                {
                        for (k=1;k<id->dims[2];k++)
                        {
                                data4d[i][j][k]=data4d[i][j][k-1] +
                                        id->dims[3]*sizeof_elem;
                        }
                }
        }


        return (void ****)data4d;
}


void *** sdf_1d_to_3d(data_id *id, void *data)
{


/* This function returns a 3-d dynamically allocated array, whose
   elements also coincide in memory with the elements of the 1-d array 
   returned by a call to sdf_read */


        char ***data3d;
        i4 sizeof_elem;
        i4 sizeof_ptr;
        pos i,j;
        sizeof_elem=id->nbpw; /* compute size of data element */
        if(id->datatype == 'c')
        {
                sizeof_elem*=2;
        }
        sizeof_ptr=sizeof(char *); /* compute size of pointer  */
        test_sizes();
        if(id->ndim != 3)
        {
                printf("sdf_1d_to_3d: id->ndim != 3; returning NULL\n");
                fflush(stdout);
                return NULL;
        }
        data3d=sdf_malloc(id->dims[0]*sizeof_ptr);
        for (i=0;i<id->dims[0];i++)
        {
                data3d[i]=malloc(id->dims[1]*sizeof_ptr);
        }


        data3d[0][0]=data;  /* here is where the data from sdf_read is */


        /* compute addresses for each start of each column */
        for (i=1;i<id->dims[0];i++)
        {
                data3d[i][0]=data3d[i-1][0]+id->dims[2]*id->dims[1]*sizeof_elem;
        }
        for (i=0;i<id->dims[0];i++)
        {
                for (j=1;j<id->dims[1];j++)
                {
                        data3d[i][j]=data3d[i][j-1]+id->dims[2]*sizeof_elem;
                }
        }
        return (void ***) data3d;
}


void ** sdf_1d_to_2d(data_id *id, void *data)
{


/* This function returns a 2-d dynamically allocated array, whose
   elements also coincide in memory with the elements of the 1-d array 
   returned by a call to sdf_read */


        char **data2d;
        i4 sizeof_elem;
        i4 sizeof_ptr;
        pos i;
        sizeof_elem=id->nbpw; /* compute size of data element */
        if(id->datatype == 'c')
        {
                sizeof_elem*=2;
        }
        sizeof_ptr=sizeof(char *); /* compute size of pointer  */
        test_sizes();
        if(id->ndim != 2)
        {
                printf("sdf_1d_to_2d: id->ndim != 2; returning NULL\n");
                fflush(stdout);
                return NULL;
        }
        data2d=sdf_malloc(id->dims[0]*sizeof_ptr);


        data2d[0]=data;  /* here is where the data from sdf_read is */


        /* compute addresses for each start of each column */
        for (i=1;i<id->dims[0];i++)
        {
                data2d[i]=data2d[i-1]+id->dims[1]*sizeof_elem;
        }
        return (void **) data2d;
}


void ** sdf_mk_2d(data_id *id)
{


/* This function returns a 2-d dynamically allocated array, whose
   elements have npbw, ndim, and dims from the id structure */


        char **data2d;
        i4 sizeof_elem;
        i4 sizeof_ptr;
        pos i;
        sizeof_elem=id->nbpw; /* compute size of data element */
        if(id->datatype == 'c')
        {
                sizeof_elem*=2;
        }
        sizeof_ptr=sizeof(char *); /* compute size of pointer  */
        test_sizes();
        if(id->ndim != 2)
        {
                printf("sdf_mk_2d: id->ndim != 2; returning NULL\n");
                fflush(stdout);
                return NULL;
        }
        data2d=sdf_malloc(id->dims[0]*sizeof_ptr);


        data2d[0]=sdf_malloc(id->dims[0]*id->dims[1]*sizeof_elem);


        /* compute addresses for each start of each column */
        for (i=1;i<id->dims[0];i++)
        {
                data2d[i]=data2d[i-1]+id->dims[1]*sizeof_elem;
        }
        return (void **) data2d;
}
void *** sdf_mk_3d(data_id *id)
{


/* This function returns a 3d dynamically allocated array, based on the
   values of nbpw, ndim, and dims in the id structure. */


        char ***data3d;
        i4 sizeof_elem;
        i4 sizeof_ptr;
        pos i,j;
        sizeof_elem=id->nbpw; /* compute size of data element */
        if(id->datatype == 'c')
        {
                sizeof_elem*=2;
        }
        sizeof_ptr=sizeof(char *); /* compute size of pointer  */
        test_sizes();
        if(id->ndim != 3)
        {
                printf("sdf_mk_3d: id->ndim != 3; returning NULL\n");
                fflush(stdout);
                return NULL;
        }
        data3d=sdf_malloc(id->dims[0]*sizeof_ptr);
        for (i=0;i<id->dims[0];i++)
        {
                data3d[i]=malloc(id->dims[1]*sizeof_ptr);
        }


        data3d[0][0]=sdf_malloc(id->dims[0]*id->dims[1]*id->dims[2]*
                sizeof_elem);


        /* compute addresses for each start of each column */
        for (i=1;i<id->dims[0];i++)
        {
                data3d[i][0]=data3d[i-1][0]+id->dims[2]*id->dims[1]*sizeof_elem;
        }
        for (i=0;i<id->dims[0];i++)
        {
                for (j=1;j<id->dims[1];j++)
                {
                        data3d[i][j]=data3d[i][j-1]+id->dims[2]*sizeof_elem;
                }
        }
        return (void ***)data3d;
}


void **** sdf_mk_4d(data_id *id)
{


/* This function returns a 4d dynamically allocated array, based on the
   values of nbpw, ndim, and dims in the id structure. */


        char ****data4d;
        i4 sizeof_elem;
        i4 sizeof_ptr;
        pos i,j,k;
        sizeof_elem=id->nbpw; /* compute size of data element */
        if(id->datatype == 'c')
        {
                sizeof_elem*=2;
        }
        sizeof_ptr=sizeof(char *); /* compute size of pointer  */
        test_sizes();
        if(id->ndim != 4)
        {
                printf("sdf_mk_4d: id->ndim != 4; returning NULL\n");
                fflush(stdout);
                return NULL;
        }
        data4d=sdf_malloc(id->dims[0]*sizeof_ptr);
        for (i=0;i<id->dims[0];i++)
        {
                data4d[i]=malloc(id->dims[1]*sizeof_ptr);


                for (j=0;j<id->dims[1];j++)
                {
                        data4d[i][j]=malloc(id->dims[2]*sizeof_ptr);
                }
        }


        data4d[0][0][0]=sdf_malloc(id->dims[0]*id->dims[1]*id->dims[2]*
                id->dims[3]*sizeof_elem);


        /* compute addresses for each start of each column */
        for (i=1;i<id->dims[0];i++)
        {
                data4d[i][0][0]=data4d[i-1][0][0]+id->dims[3]*id->dims[2]*
                        id->dims[1]*sizeof_elem;
        }
        for (i=0;i<id->dims[0];i++)
        {
                for (j=1;j<id->dims[1];j++)
                {
                        data4d[i][j][0]=data4d[i][j-1][0]+id->dims[3]*
                                id->dims[2]*sizeof_elem;
                }
        }
        
        for (i=0;i<id->dims[0];i++)
        {
                for (j=0;j<id->dims[1];j++)
                {
                        for (k=1;k<id->dims[2];k++)
                        {
                                data4d[i][j][k]=data4d[i][j][k-1] +
                                        id->dims[3]*sizeof_elem;
                        }
                }
        }


        return (void ****)data4d;
}


void ***** sdf_mk_5d(data_id *id)
{


/* This function returns a 5d dynamically allocated array, based on the
   values of nbpw, ndim, and dims in the id structure. */


        char *****data5d;
        i4 sizeof_elem;
        i4 sizeof_ptr;
        pos i,j,k,l;
        sizeof_elem=id->nbpw; /* compute size of data element */
        if(id->datatype == 'c')
        {
                sizeof_elem*=2;
        }
        sizeof_ptr=sizeof(char *); /* compute size of pointer  */
        test_sizes();
        if(id->ndim != 5)
        {
                printf("sdf_mk_5d: id->ndim != 5; returning NULL\n");
                fflush(stdout);
                return NULL;
        }
        data5d=sdf_malloc(id->dims[0]*sizeof_ptr);
        for (i=0;i<id->dims[0];i++)
        {
                data5d[i]=malloc(id->dims[1]*sizeof_ptr);


                for (j=0;j<id->dims[1];j++)
                {
                        data5d[i][j]=malloc(id->dims[2]*sizeof_ptr);
                        for (k=0;k<id->dims[2];k++)
                        {
                                data5d[i][j][k]=malloc(id->dims[3]*sizeof_ptr);
                        }
                }
        }


        data5d[0][0][0][0]=sdf_malloc(id->dims[0]*id->dims[1]*id->dims[2]*
                id->dims[3]*id->dims[4]*sizeof_elem);


        /* compute addresses for each start of each column */
        for (i=1;i<id->dims[0];i++)
        {
                data5d[i][0][0][0]=data5d[i-1][0][0][0]+id->dims[4]*id->dims[3]*
                        id->dims[2]*id->dims[1]*sizeof_elem;
        }


        for (i=0;i<id->dims[0];i++)
        {
                for (j=1;j<id->dims[1];j++)
                {
                        data5d[i][j][0][0]=data5d[i][j-1][0][0]+id->dims[4]*
                                id->dims[3]*id->dims[2]*sizeof_elem;
                }
        }
        
        for (i=0;i<id->dims[0];i++)
        {
                for (j=0;j<id->dims[1];j++)
                {
                        for (k=1;k<id->dims[2];k++)
                        {
                                data5d[i][j][k][0]=data5d[i][j][k-1][0] +
                                        id->dims[4]*id->dims[3]*sizeof_elem;
                        }
                }
        }
        
        for (i=0;i<id->dims[0];i++)
        {
                for (j=0;j<id->dims[1];j++)
                {
                        for (k=0;k<id->dims[2];k++)
                        {
                                for (l=1;l<id->dims[3];l++)
                                {
                                        data5d[i][j][k][l]=data5d[i][j][k][l-1]
                                        + id->dims[4]*sizeof_elem;
                                }
                        }
                }
        }


        return (void *****)data5d;
}


void sdf_free_5d(data_id *id, void ***** data5d)
{
        pos i,j,k;
        sdf_free(data5d[0][0][0][0]);
        for (i=0;i<id->dims[0];i++)
        {
                for (j=0;j<id->dims[1];j++)
                {
                        for(k=0;k<id->dims[2];k++)
                        {
                                sdf_free(data5d[i][j][k]);
                        }


                        sdf_free(data5d[i][j]);
                }


                sdf_free(data5d[i]);


        }
        sdf_free(data5d);
}


void sdf_free_4d(data_id *id, void **** data4d)
{
        pos i,j;
        sdf_free(data4d[0][0][0]);
        for (i=0;i<id->dims[0];i++)
        {
                for (j=0;j<id->dims[1];j++)
                {
                        sdf_free(data4d[i][j]);
                }


                sdf_free(data4d[i]);


        }
        sdf_free(data4d);
}


void    sdf_free_3d(data_id *id, void *** data3d)
{


/* This function frees a 3d dynamically allocated array.  It arguments are 
   the id structure used to create the array and the array itself. */


        pos i;
        sdf_free(data3d[0][0]);
        for (i=0;i<id->dims[0];i++) 
        {
                sdf_free(data3d[i]);
        }
        sdf_free(data3d);
        return;
}


void    sdf_free_2d(data_id *id, void ** data2d)
{


/* This function frees a 2d dynamically allocated array.  It arguments are 
   the id structure used to create the array and the array itself. */


        sdf_free(data2d[0]);
        sdf_free(data2d);
        return;
}


data_id * sdf_cp_id(data_id *id_0)
{


/* function that returns a copy of the id structure that is its argument */


        data_id *id;
        i4 i,nlab;
        id=sdf_malloc(sizeof(data_id));
        id->order=id_0->order;
        id->nbpw=id_0->nbpw;
        id->ndim=id_0->ndim;
        id->datatype=id_0->datatype;
        nlab=strlen(id_0->label)+1;
        id->label=malloc(nlab*sizeof(char));
        id->dims=malloc(id->ndim*sizeof(pos));
        for (i=0;i<nlab;i++)
        {
                *(id->label+i)=*(id_0->label+i);        
        }
        for (i=0;i<id->ndim;i++) 
        {
                id->dims[i]=id_0->dims[i];
        }
        return id;
}
void sdf_free_id(data_id *id)
{
/* free the members of id, then free id itself.  It uses sdf_free so it
 * can be called from outside of sdf. */


        sdf_free(id->label);
        sdf_free(id->dims);
        sdf_free(id);


        return;
}
data_id *sdf_create_id(i4 order, char *label, char dtype, i4 nbpw, 
                i4 ndim, pos *dims)


/* Returns pointer to a data_id type structure, whose members are assigned
 * to values from the calling argument list.  The caller program needs to
 * free the structure with sdf_free_id once they are done using it, to avoid
 * a memory leak. */


{
        data_id * idnew;
        i4 i, nchar;


        idnew=sdf_malloc(sizeof(data_id));
        idnew->order=order;
        nchar=strlen(label)+1;
        idnew->label=sdf_malloc((pos)nchar);
        for (i=0;i<nchar-1;i++)
        {
                *(idnew->label+i) = *(label+i);
        }
        *(idnew->label+nchar-1) = '\0';
        idnew->datatype=dtype;
        idnew->nbpw=nbpw;
        idnew->ndim=ndim;
        idnew->dims=sdf_malloc((pos)ndim*sizeof(pos));
        for (i=0;i<ndim;i++)
        {
                *(idnew->dims+i)=*(dims+i);
        }
        return idnew;
}


void *sdf_malloc(pos nmemb)
{
        /* purpose of sdf_malloc and sdf_free are to provide self-consistent
         * versions of malloc and free which can be called from outside the
         * sdf library (like in a calling program) to avoid inconsistencies
         * between versions of malloc/free from within sdf and versions that
         * might use different versions of these functions in libraries
         * outside of sdf.
         */


        void *ptr;
        ptr=malloc(nmemb);
        if(ptr == NULL)
        {
                printf("sdf_malloc: memory allocation failed\n");
                fflush(stdout);
                exit(1);
        }
        return ptr;
}
void sdf_free(void *ptr)
{
        free(ptr);
        return;
}



i4 sdf_wb(char *fname, pos *datapos, pos nelem, i4 nbpw, FILE *fp, void *data)


{


/* Just do a simple, large-endian binary write, starting at byte *datapos, doing
   nelem words of length nbpw.  If data is complex, user must mpy nelem by
   2 (not adjusted in this function).  If fp != NULL, assume file already
   open */


FILE *fpc;
char temp[300];
int fctest, fsr;
int ise,ibe;
fctest=0;


/* Make certain (pos) is defined as 64bit int, otherwise, exit */


test_sizes(); /* NEW WAY */
fpc=NULL;


ibe = is_big_endian ();
ise = 0;
if (ibe == 0) ise = 1;    /* set flag for byteswapping if small endian */


/* if file fp not already open, then open file with both read/write 
   permissions */


if(fp == NULL) 
{
        fpc=fopen(fname,"rb+");
}
else
{
        fpc=fp;
}


if(fpc == NULL)
{
        fpc=fopen(fname,"wb");
        if(fpc == NULL)
        {
                printf("sdf_wb:  Can't open file %s\n",fname);
                printf("sdf_wb: errno = %d\n",errno);
                fflush(stdout);
                exit(1);
        }
}



/* Write out the data: */


fsr=fseeko(fpc,*datapos,SEEK_SET);
if(fsr != 0)
{
        output_int64(temp,*datapos);
        printf("sdf_wb: unable to seek to position %s\n",temp);
        printf("sdf_wb: errno = %d\n",errno);
        fflush(stdout);
        exit(1);
}


if(ise) byteswap((void *)data, (pos) nelem, (i4)nbpw);
fwrite(data,nbpw,nelem,fpc);
if(ise) byteswap((void *)data, (pos) nelem, (i4)nbpw);


*datapos=(pos) ftello(fpc);


/* Finish up */


fflush(fpc);


/* Rick Devore said he thought he wanted to truncate the file */


sdf_file_truncate(fpc, *datapos);


if(fp == NULL)
{
        fctest=fclose(fpc);
}
if(fp == NULL && fctest != 0)
{
        printf("sdf_wb:  failure to close file %s\n",fname);
        printf("sdf_wb: errno = %d\n",errno);
}
 
/* We're done */


return 0;
}


i4 sdf_rb(char *fname, pos *datapos, pos nelem, i4 nbpw, FILE *fp, void *data)


{


/* Just do a simple, binary read, starting at byte *datapos, doing
   nelem words of length nbpw.  If data is complex, user must mpy nelem by
   2 (not adjusted in this function) */


FILE *fpc;
char temp[300];
int fctest, fsr;
int ise,ibe;
pos frtest;
fctest=0;

/* Make certain (pos) is defined as 64bit int, otherwise, exit */


test_sizes(); /* NEW WAY */
fpc=NULL;


ibe = is_big_endian ();
ise = 0;
if (ibe == 0) ise = 1;    /* set flag for byteswapping if small endian */


/* open file with read permission */


if(fp == NULL)
{
        fpc=fopen(fname,"rb");
}
else
{
        fpc=fp;
}


if(fpc == NULL)
{
        printf("sdf_rb:  Can't open file %s\n",fname);
        printf("sdf_rb: errno = %d\n",errno);
        fflush(stdout);
        exit(1);
}



/* Read the data: */


fsr=fseeko(fpc,*datapos,SEEK_SET);
        /* test for unsuccessful seek */
if(fsr != 0)
{
        output_int64(temp,*datapos);
        printf("sdf_rb: unable to seek to position %s\n",temp);
        printf("sdf_rb: errno = %d\n",errno);
        fflush(stdout);
        exit(1);
}


frtest=(pos) fread(data,nbpw,nelem,fpc);
if(frtest < nelem)
{
        /* test for unsuccessful read: */
        printf("sdf_rb: read error: attempt to read past EOF?\n");
        printf("sdf_rb: errno = %d\n",errno); 
        output_int64(temp,frtest);
        printf("sdf_rb: frtest = %s\n",temp);
        output_int64(temp,nelem);
        printf("sdf_rb: nelem = %s\n",temp);
        fflush(stdout);
        exit(1);
}
if(ise) byteswap((void *)data, (pos) nelem, (i4)nbpw);


*datapos=(pos) ftello(fpc);


/* Finish up */


fflush(fpc);
if(fp == NULL)
{
        fctest=fclose(fpc);
}
if(fp == NULL && fctest != 0)
{
        printf("sdf_rb:  failure to close file %s\n",fname);
        printf("sdf_rb: errno = %d\n",errno);
}
 
/* We're done */


return 0;
}


void sdf_wb_f77__(char *fnamef77, intf8 *dataposf77, intf8 *nelemf77, 
        intf *nbpwf77, FILE **fp, void *dataf77,  intf_clen len_fnamef77, 
        intf_clen len_dataf77 )
{
sdf_wb_f77(fnamef77, dataposf77, nelemf77, nbpwf77, fp, dataf77,
    len_fnamef77, len_dataf77);
return;
}


void sdf_wb_f77_(char *fnamef77, intf8 *dataposf77, intf8 *nelemf77, 
        intf *nbpwf77, FILE **fp, void *dataf77,  intf_clen len_fnamef77, 
        intf_clen len_dataf77 )
{
sdf_wb_f77(fnamef77, dataposf77, nelemf77, nbpwf77, fp, dataf77,
    len_fnamef77, len_dataf77);
return;
}


void sdf_wb_f77(char *fnamef77, intf8 *dataposf77, intf8 *nelemf77, 
        intf *nbpwf77, FILE **fp, void *dataf77,  intf_clen len_fnamef77, 
        intf_clen len_dataf77 )
{


char *fname;
i4 nbpw;
pos datapos,nelem;


/* 1st, print out the sizes of f77 char arrays (debug) */


/*
printf("len_fnamef77 = %d\n",len_fnamef77);
printf("len_dataf77 = %d\n",len_dataf77);
*/


/* Convert fortran string fnamef77 to C-strings fname */


fs2cs(&fname,fnamef77,len_fnamef77);


/* Debug - print out C-string fname */


/*
printf("fname = %s\n",fname);
*/


nbpw=(i4) *nbpwf77;
nelem = (pos) *nelemf77;
datapos= (pos) *dataposf77;


sdf_wb(fname, &datapos, nelem, nbpw, *fp, (void *)dataf77);
*dataposf77 = (intf8) datapos;
free(fname);
}


void sdf_rb_f77__(char *fnamef77, intf8 *dataposf77, intf8 *nelemf77, 
        intf *nbpwf77, FILE **fp, void *dataf77,  intf_clen len_fnamef77, 
        intf_clen len_dataf77 )
{
sdf_rb_f77(fnamef77, dataposf77, nelemf77, nbpwf77, fp, dataf77,
    len_fnamef77, len_dataf77);
return;
}


void sdf_rb_f77_(char *fnamef77, intf8 *dataposf77, intf8 *nelemf77, 
        intf *nbpwf77, FILE **fp, void *dataf77,  intf_clen len_fnamef77, 
        intf_clen len_dataf77 )
{
sdf_rb_f77(fnamef77, dataposf77, nelemf77, nbpwf77, fp, dataf77,
    len_fnamef77, len_dataf77);
return;
}


void sdf_rb_f77(char *fnamef77, intf8 *dataposf77, intf8 *nelemf77, 
        intf *nbpwf77, FILE **fp, void *dataf77,  intf_clen len_fnamef77, 
        intf_clen len_dataf77 )
{


char *fname;
i4 nbpw;
pos datapos,nelem;


/* 1st, print out the sizes of f77 char arrays (debug) */


/*
printf("len_fnamef77 = %d\n",len_fnamef77);
printf("len_dataf77 = %d\n",len_dataf77);
*/


/* Convert fortran string fnamef77 to C-string fname */


fs2cs(&fname,fnamef77,len_fnamef77);


/* Debug - print out C-string fname */


/*
printf("fname = %s\n",fname);
*/


nbpw=(i4) *nbpwf77;
nelem = (pos) *nelemf77;
datapos= (pos) *dataposf77;


sdf_rb(fname, &datapos, nelem, nbpw, *fp, (void *)dataf77);
*dataposf77 = (intf8) datapos;
free(fname);
}


void fopen_f77__(char *fnamef77, char * modef77, FILE **fp,
               intf_clen len_fnamef77, intf_clen len_modef77)
{
fopen_f77(fnamef77, modef77, fp, len_fnamef77, len_modef77);
return ;
}
void fopen_f77_(char *fnamef77, char * modef77, FILE **fp,
               intf_clen len_fnamef77, intf_clen len_modef77)
{
fopen_f77(fnamef77, modef77, fp, len_fnamef77, len_modef77);
return ;
}
void fopen_f77(char * fnamef77, char * modef77, FILE **fp,
               intf_clen len_fnamef77, intf_clen len_modef77)
{


/* Provide a fortran callable version of the C fn fopen.  In the fortran
   calling program, fp should be declared as an integer*8 or (kind=8). */


char *fname;
char *mode;


fs2cs(&fname,fnamef77,len_fnamef77);
fs2cs(&mode,modef77,len_modef77);


*fp=fopen(fname,mode);


free(fname);
free(mode);


}


void fclose_f77__(FILE **fp)
{
fclose_f77(fp);
return ;
}
void fclose_f77_(FILE **fp)
{
fclose_f77(fp);
return ;
}
void fclose_f77(FILE **fp)
{


/* provide a fortran callable version of the C fn fclose */


fclose(*fp);


}


void fflush_f77__(FILE **fp)
{
fflush_f77(fp);
return ;
}
void fflush_f77_(FILE **fp)
{
fflush_f77(fp);
return ;
}
void fflush_f77(FILE **fp)
{


/* provide a fortran callable version of the C fn fflush */


fflush(*fp);


}


void fwrite_f77__(void *data, intf * nbpwf77, intf8 *nelemf77, FILE **fp)
{
fwrite_f77(data,nbpwf77,nelemf77,fp);
return ;
}
void fwrite_f77_(void *data, intf * nbpwf77, intf8 *nelemf77, FILE **fp)
{
fwrite_f77(data,nbpwf77,nelemf77,fp);
return ;
}
void fwrite_f77(void *data, intf * nbpwf77, intf8 * nelemf77, FILE **fp)
{


/* Fortran callable version of the fwrite C function.  *fp assumed open */


pos nelem, nbpw;
nelem=(pos) *nelemf77;
nbpw=(pos) *nbpwf77;


fwrite(data, nbpw, nelem, *fp);


}


void fread_f77__(void *data, intf * nbpwf77, intf8 *nelemf77, FILE **fp)
{
fread_f77(data,nbpwf77,nelemf77,fp);
return ;
}
void fread_f77_(void *data, intf * nbpwf77, intf8 *nelemf77, FILE **fp)
{
fread_f77(data,nbpwf77,nelemf77,fp);
return ;
}
void fread_f77(void *data, intf * nbpwf77, intf8 * nelemf77, FILE **fp)
{


/* Fortran callable version of the fread C function.  *fp assumed open */


pos nelem, nbpw;
nelem=(pos) *nelemf77;
nbpw=(pos) *nbpwf77;


fread(data, nbpw, nelem, *fp);


}


void fseeko_f77__(FILE **fp, intf8 *filepos)
{
fseeko_f77(fp, filepos);
return ;
}
void fseeko_f77_(FILE **fp, intf8 *filepos)
{
fseeko_f77(fp, filepos);
return ;
}
void fseeko_f77(FILE **fp, intf8 *filepos)
{


/* provide a fortran callable version of the C fn fseeko */


int result;
char temp[300];
result = fseeko(*fp, (pos) *filepos, SEEK_SET);
if (result != 0)
{
        output_int64(temp,(pos)*filepos);
        printf("fseeko_f77: seek failed for filepos = %s \n", temp);
        fflush(stdout);
        exit(1);
        /* return; */
}
return;

}



void ftello_f77__(FILE **fp, intf8 *filepos)
{
ftello_f77(fp, filepos);
return ;
}
void ftello_f77_(FILE **fp, intf8 *filepos)
{
ftello_f77(fp, filepos);
return ;
}
void ftello_f77(FILE **fp, intf8 *filepos)
{


/* provide a fortran callable version of the C fn ftello */


*filepos = ftello(*fp);


}


pos file_sz(char *fname)

/* Function to return the size of the file.  Uses generic stat function,
   gen_stat. */

{
sdf_stat fileinfo;
int result;
pos filesz;
result=gen_stat(fname,&fileinfo);
if(result == 0)
{
        filesz=(pos) fileinfo.st_size;
        return filesz;
}
else
{
        printf("file_sz: bad return from generic stat = %d\n",result);
        printf("file_sz: file %s non-existent or permission problem?\n",fname);
        fflush(stdout);
        exit(1);
        /* return (pos)0; */
}

/* Code cannot get here, but this makes SGI Origin compiler happy */
return (pos)0;

}


void file_sz_f77__(char *fnamef77, intf8 *fsize, intf_clen len_fnamef77)
{
file_sz_f77(fnamef77,fsize,len_fnamef77);
return;
}


void file_sz_f77_(char *fnamef77, intf8 *fsize, intf_clen len_fnamef77)
{
file_sz_f77(fnamef77,fsize,len_fnamef77);
return;
}


void file_sz_f77(char *fnamef77, intf8 *fsize, intf_clen len_fnamef77)

/* Fortran callable function to return size of file. */

{
char *fname;
fs2cs(&fname,fnamef77,len_fnamef77);
*fsize = (intf8) file_sz(fname);
free(fname);
return;
}


void byteswap_f77__(void *data, intf8 *nelemf77, intf * nbpwf77)
{
byteswap_f77(data,nelemf77,nbpwf77);
return ;
}
void byteswap_f77_(void *data, intf8 *nelemf77, intf * nbpwf77)
{
byteswap_f77(data,nelemf77,nbpwf77);
return ;
}
void byteswap_f77(void *data, intf8 *nelemf77, intf * nbpwf77)
{


/* Fortran callable version of the C byteswap function from this library. */


pos nelem;
i4 nbpw;
nelem=(pos) *nelemf77;
nbpw=(i4) *nbpwf77;


byteswap((void *)data, nelem, nbpw);


}


void ibe_f77__(intf * ibef77)
{
ibe_f77(ibef77);
return ;
}
void ibe_f77_(intf * ibef77)
{
ibe_f77(ibef77);
return ;
}
void ibe_f77(intf * ibef77)
{


/* Fortran callable version of the C is_big_endian()  function from 
   this library. */


*ibef77=(intf) is_big_endian();


}


void sdf_ftc_f77__(FILE **fp, intf8 *datasizef77)
{
sdf_ftc_f77(fp,datasizef77);
return ;
}
void sdf_ftc_f77_(FILE **fp, intf8 *datasizef77)
{
sdf_ftc_f77(fp,datasizef77);
return ;
}
void sdf_ftc_f77(FILE **fp, intf8 *datasizef77)
{


/* provides a fortran callable version of the C fn sdf_file_truncate */


pos datasize;
datasize= (pos) *datasizef77;


sdf_file_truncate(*fp, datasize);


}


void ffstdout_f77__()
{
ffstdout_f77();
return;
}


void ffstdout_f77_()
{
ffstdout_f77();
return;
}


void ffstdout_f77()
{
        /* This function just provides a fortran callable version of
         * fflush(stdout) */
        fflush(stdout);
        return;
}


int i4cmp(const void *vp, const void *vq)
{
/* comparison fn for 2 i4 integers in qsort. */ 


        const i4 *p = vp;
        const i4 *q = vq;
        int result=0;
        /*
         i4 diff = *p - *q;
        */
        if(*p < *q) result=(i4)(-1);
        if(*p > *q) result=(i4)1;
        if(*p == *q) result=(i4)0;
        return result;
}


void test_sizes()
{


/*
test to make sure all defined variables types have correct size 
*/


int i2size;
int i4size;
int u4size;
int i8size;
int f4size;
int f8size;
int possize;
int intfsize;
int intf8size;
# ifdef INTF_CLEN_64
  size_t intf_clensize;
# else
  int intf_clensize;
# endif
i2size=sizeof(i2);
i4size=sizeof(i4);
u4size=sizeof(u4);
i8size=sizeof(i8);
f4size=sizeof(f4);
f8size=sizeof(f8);
possize=sizeof(pos);
intfsize=sizeof(intf);
intf8size=sizeof(intf8);
intf_clensize=sizeof(intf_clensize);
if (i2size != 2)
{
        printf("test_sizes: sizeof(i2) = %d != 2\n",i2size);
        fflush(stdout);
        exit(1);
}
if (i4size != 4)
{
        printf("test_sizes: sizeof(i4) = %d != 4\n",i4size);
        fflush(stdout);
        exit(1);
}
if (u4size != 4)
{
        printf("test_sizes: sizeof(u4) = %d != 4\n",u4size);
        fflush(stdout);
        exit(1);
}
if (f4size != 4)
{
        printf("test_sizes: sizeof(f4) = %d != 4\n",f4size);
        fflush(stdout);
        exit(1);
}
if (i8size != 8)
{
        printf("test_sizes: sizeof(i8) = %d != 8\n",i8size);
        fflush(stdout);
        exit(1);
}
if (f8size != 8)
{
        printf("test_sizes: sizeof(f8) = %d != 8\n",f8size);
        fflush(stdout);
        exit(1);
}
if (possize != 8)
{
        printf("test_sizes: sizeof(pos) = %d != 8\n",possize);
        fflush(stdout);
        exit(1);
}
if (intf8size != 8)
{
        printf("test_sizes: sizeof(intf8size) = %d != 8\n",intf8size);
        fflush(stdout);
        exit(1);
}

#ifdef INTF_CLEN_64
if (intf_clensize != 8)
{
        printf("test_sizes: sizeof(intf_clensize) = %d != 8\n",(int) 
        intf_clensize);
        fflush(stdout);
        exit(1);
}
#else
if (intf_clensize != 4)
{
        printf("test_sizes: sizeof(intf_clensize) = %d != 4\n",intf_clensize);
        fflush(stdout);
        exit(1);
}
#endif

#ifdef INTF_64
if (intfsize != 8)
{
        printf("test_sizes: sizeof(intf) = %d != 8\n",intfsize);
        printf("NOTE: above test is for default 64-bit fortran integer\n");
        printf("If you want a default fortran integer size, you must\n");
        printf("undefine INTF_64 in sdf_subs.h\n");

        fflush(stdout);
        exit(1);
}
#else
if (intfsize != sizeof(double)/2)
{
        printf("test_sizes: sizeof(intf) = %d != sizeof(double)/2\n",intfsize);
        printf("NOTE: above test is the fortran standard for the size of an\n");
        printf("integer.  If you want a 64-bit default integer, you must \n");
        printf("define INTF_64 1 in sdf_subs.h\n");


        /* NOTE: above test is the fortran standard for the size of an
        integer.  If you use something else, you'll need to change this test
        accordingly! */

        fflush(stdout);
        exit(1);
}
#endif


}

/*Original line 6952: #ifdef WIN32*/
#if defined(WIN32) && !defined(__MINGW32__)


/* define strtok_r in terms of strtok; 
   fseeko, ftello for MS windoze, using 
   fsetpos, fgetpos, which are
   supposedly large-file compatible, unlike fseek and ftell. */
char *strtok_r(char* s1, char* sep, char** last)
{
        /* supposedly strtok itself is threadsafe in Windoze, so just use it */
        return strtok(s1, sep);
}


int fseeko(FILE *fd, pos offset, int whence)
{
        /* note -- whence implicitly assumed = SEEK_SET
         * therefore, this isn't going to work to e.g. get filesize
         * where you want to set whence to SEEK_END 
         * But this could be fixed if you invoke the _stati64 function
         * to find filesize, then compute offset relative to SEEK_END
         */


        return fsetpos(fd, &offset);
}
pos ftello(FILE *fd)
{
        pos offset;
        fgetpos(fd, &offset);
        return offset;
}
int gen_stat(char *fname, sdf_stat *fileinfo)
{
        return _stati64(fname, fileinfo);
}
int sdf_file_truncate(FILE *fd, pos size)
{
/* Windoze version of 64-bit compatible file_truncate */
        int fn;
        pos sizetest;
        HANDLE fh; /* a windoze file handle definition */
        fseeko(fd,size,SEEK_SET);
        sizetest=ftello(fd);
        if(sizetest != size)
        {
                char temp[300];
                printf("sdf_file_truncate win: ftello value != fseeko value\n");
                output_int64(temp,sizetest);
                printf("sdf_file_truncate win: sizetest = %s\n",temp);
                fflush(stdout);
                fclose(fd);
                exit(1);
        }
        fn=fileno(fd);
        fflush(fd);
        fh= (HANDLE) _get_osfhandle(fn); /* get windoze file handle */
        if(SetEndOfFile(fh) ==0) /* this sets 64-bit compatible file size */
        {
                char buf[100];
                sprintf(buf, "SetEndOfFile failed in sdf_file_truncate: %d\n",
                                (int)GetLastError() );
                printf("%s\n",buf);
                fflush(stdout);
                return 1;
        }
        fflush(fd);
        return 0;
}
#else
int gen_stat(char *fname, sdf_stat *fileinfo)
{
        return stat(fname, fileinfo);
}
int sdf_file_truncate(FILE *fd, pos size)
/* Unix-like version of sdf_file_truncate using ftruncate */
{
        int fn;
        fn=fileno(fd);
        fflush(fd);
        if(ftruncate(fn,size))
        {
                printf("problem ftruncating the file in sdf_file_truncate\n");
                fflush(stdout);
                fclose(fd);
                exit(1);
        }
        fflush(fd);
        return 0;
}
#endif
