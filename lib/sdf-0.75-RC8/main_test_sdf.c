#include <sdf_subs.h>
int main(int argc, char ** argv);

int main(int argc, char ** argv)
{
/*
This test C program writes and reads some test datasets into an sdf file.
It illustrates how one can use sdf i/o functions from C.

This program takes one argument from the command line, which is 
expected to be "ro" (read only); "wo" (write only); "qo"
(query only), or "rw" (read & write).  Default with no argument is rw.

*/

/* junk.sdf will be the name of the sdf file created */

char * fname = "junk.sdf";
char * label ;
/*
	char * label = "";
*/
char rwopts[100]="";
char temp[300];
f4 *data;
f8 *ddata;
char *cdata;
float deltat;
char dtype;
/* some test string data */
char * stringdata="Flatulence simulation with ridget hurling";
data_id *id, *idnew;
pos datasize;
pos dims[20];
i4 i, j, ndata, nchar, order, nbpw, ndim;
i4 ro,wo,rw;
#if defined (__SVR4) && defined (__sun) 
  printf("Solaris OS detected\n");
#elif  defined __WIN32
  printf("Windows OS detected\n");
#else
  printf("Other Unix-like OS detected\n");
#endif
fflush(stdout);
test_sizes(); /*test to make sure variable types correctly defined */
if(sizeof(pos) !=8)
{
	printf("main_test_sdf.c was not compiled correctly.\n");
	fflush(stdout);
	exit(1);
}
ro=0;wo=0;rw=1;
if(argc >= 2) strncpy(rwopts,argv[1],99);

if(!strncmp("ro",rwopts,2))
{
	rw=0;
	ro=1;
	wo=0;
}
if(!strncmp("wo",rwopts,2))
{
	rw=0;
	wo=1;
	ro=0;
}
if(!strncmp("qo",rwopts,2))
{
	rw=0;
	ro=0;
	wo=0;
}
	
if( rw == 1 || wo == 1)

{
/*  CASE 0: output scalar float deltat to the sdf file: */
	deltat=5e-4;
	data=sdf_malloc(sizeof(f4));
	*data=deltat;
	order=0;
	label="deltat(s)";
	dtype='f';
	ndim=1;
	nbpw=4;
	dims[0]=(pos)1;
	id=sdf_create_id(order,label,dtype,nbpw,ndim,dims);
	sdf_write(fname, id,(void *)data);
	sdf_free_id(id);
	sdf_free(data);

/* CASE 1: output a string to the sdf file: */

	nchar=strlen(stringdata);
	cdata=sdf_malloc(nchar);
	strncpy((char *)cdata,stringdata,nchar);
	order=0;
	label="stringdata";
	dtype='b';
	nbpw=1;
	ndim=1;
	dims[0]=(pos)nchar;
	id=sdf_create_id(order,label,dtype,nbpw,ndim,dims);
	sdf_write(fname, id, (void *)cdata);
	sdf_free_id(id);
	sdf_free(cdata);

/* CASE 2:  output a singly dimensioned floating array */

	order=0;
	label="lambda";
	dtype='f';
	nbpw=4;
	ndim=1;
	dims[0]=(pos)101;
	id=sdf_create_id(order,label,dtype,nbpw,ndim,dims);
	data=sdf_malloc(*(id->dims+0)*id->nbpw);
	for (i=0;i<*(id->dims+0);i++)
	{
		*(data+i)=(float)i;
	}

	sdf_write(fname,id,(void *)data);
	sdf_free_id(id);
	sdf_free(data);

/* CASE 3:  output a doubly dimensioned double precision array */

	order=0;
	label="vdotn";
	dtype='f';
	nbpw=8;
	ndim=2;
	dims[0]=12;
	dims[1]=13;
	id=sdf_create_id(order,label,dtype,nbpw,ndim,dims);

/* Note that in the following, the 2d array ddata is statically allocated.  
   In this case, the 2d array ddata will be output fine with sdf if 
   ddata is referenced as a simple pointer as in this sdf_write call. 

   However, for 2d arrays which are dynamically allocated using
   nested malloc calls, this will not work correctly, as the memory locations
   for ddata will not necessarily be contiguous, which sdf_write requires.  
   Instead, the best
   way to allocate the 2d array is as a 1-d array, but then calculate the
   standard row-major memory locations as if it were a 2d array.  This is
   the best way to define allocatable multi-dimensional arrays for use with
   sdf software in C.  The statements below that are commented out show
   how this would be done with ddata treated as a dynamically allocated array.

   Note in this example that ddata gets defined within its own code block,
   so it doesn't conflict with the ddata that gets used later on when reading
   the data back in. */

	{
	/* uncomment next line for dyn. alloc. */
	/*
	ddata=malloc(nx*ny*id->nbpw); 
	*/
	double ddata[12][13]; /* comment this stmt out for dyn. mem. alloc. */
		for (i=0;i<*(id->dims+0);i++)
		{
			for(j=0;j<*(id->dims+1);j++)
			{
				/* uncomment next stmt for dyn. alloc case */
      			/* *(ddata+i*(*(id->dims+1))+j)=(f8)100*i+(f8)j; */

				/* comment out next stmt for dyn. alloc case */
				ddata[i][j]=(f8)100*i+(f8)j; 
			}
		}			
		sdf_write(fname,id,(void *)ddata);
	}
	sdf_free_id(id);
	/* uncomment next stmt for dyn. alloc. case */
	/*
	sdf_free(ddata); 
	*/
}

/* CASE 4:  Query the output file */

        ndata=sdf_ndat(fname);
        printf("ndata (from sdf_ndat) = %d\n",ndata);
        ndata=sdf_query (fname);
        printf("ndata (from sdf_query) = %d\n",ndata);

if((rw == 1) || (ro == 1))
{

/* CASE 5: read the 2nd dataset (dataset 1): */

	idnew = sdf_read(fname,1,(void **)&cdata);
        printf("idnew->order = %u\n",idnew->order);
        printf("idnew->label = %s\n",idnew->label);
        printf("idnew->datatype = %c\n",idnew->datatype);
        printf("idnew->nbpw = %u\n",idnew->nbpw);
        printf("idnew->ndim = %u\n",idnew->ndim);
        for (i=0;i<idnew->ndim;i++)
        {
		/* output_int64 is a function that converts a 64 bit
                integer into a string.  Necessary because there's no
		platform indep. format stmt for 64-bit integers */

		output_int64(temp,idnew->dims[i]);
		printf("idnew->dims[%u] = %s\n",i,temp);
        }
        datasize=data_size(idnew);
	printf("dataset 2 = cdata = ");
	for (i=0;i<datasize;i++)
	{
		if(isprint(*cdata+i)) printf("%c",*(cdata+i));
	}
	printf("\n");
	sdf_free_id(idnew);
	sdf_free(cdata);

/* CASE 6:  Read the 0th dataset: (a single float) */

	idnew = sdf_read(fname,0,(void **)&data);
        printf("idnew->order = %u\n",idnew->order);
        printf("idnew->label = %s\n",idnew->label);
        printf("idnew->datatype = %c\n",idnew->datatype);
        printf("idnew->nbpw = %u\n",idnew->nbpw);
        printf("idnew->ndim = %u\n",idnew->ndim);
        for (i=0;i<idnew->ndim;i++)
        {
		output_int64(temp,idnew->dims[i]);
		printf("idnew->dims[%u] = %s\n",i,temp);
        }
        datasize=data_size(idnew);
	for (i=0;i<datasize;i++)
	{
		printf(" 0th dataset: data = %g",*(data+i));
	}
	printf("\n");
	sdf_free_id(idnew);
	sdf_free(data);

/* CASE 7:  Read the 2nd dataset: (a 1-d array of floats) */

	idnew = sdf_read(fname,2,(void **)&data);
        printf("idnew->order = %u\n",idnew->order);
        printf("idnew->label = %s\n",idnew->label);
        printf("idnew->datatype = %c\n",idnew->datatype);
        printf("idnew->nbpw = %u\n",idnew->nbpw);
        printf("idnew->ndim = %u\n",idnew->ndim);
        for (i=0;i<idnew->ndim;i++)
        {
		output_int64(temp,idnew->dims[i]);
		printf("idnew->dims[%u] = %s\n",i,temp);
        }
        datasize=data_size(idnew);
	printf("dataset 2: data = ");
	for (i=0;i<datasize;i++)
	{
		printf("%g ",*(data+i));
	}
	printf("\n");
	sdf_free_id(idnew);
	sdf_free(data);

/* CASE 8:  Read the 3rd dataset: (a 2d array of doubles) */

	idnew = sdf_read(fname,3,(void **)&ddata);
        printf("idnew->order = %u\n",idnew->order);
        printf("idnew->label = %s\n",idnew->label);
        printf("idnew->datatype = %c\n",idnew->datatype);
        printf("idnew->nbpw = %u\n",idnew->nbpw);
        printf("idnew->ndim = %u\n",idnew->ndim);
        for (i=0;i<idnew->ndim;i++)
        {
		output_int64(temp,idnew->dims[i]);
		printf("idnew->dims[%u] = %s\n",i,temp);
        }
        datasize=data_size(idnew);
	printf("ddata = \n");
	for (i=0;i<datasize;i++)
	{
		printf("%g ",*(ddata+i));
	}
	printf("\n");
	sdf_free_id(idnew);
	sdf_free(ddata);
}
return 0;
}
