#include <sdf_subs.h>
/* function prototypes not in sdf_subs.h */
int print_data(data_id *id, char *data, pos nelem);

int main(int argc, char **argv)
{
data_id * id;
char temp[301],ans[100];
char stuff[100];
char fname[100];
char labmatch[100];
i4 i,ii,iorder,ise,ibe,ndat;
i4 *matchind;
pos datasize, hdrpos, datapos, hdrsize, nelem;
void *data;
test_sizes();
printf("\nsdf_browse: version 0.75\n");
ise=1;
ibe=is_big_endian();
if(ibe) ise=0;
if(ise)
{
	printf("\nsmall endian platform, byteswapping implemented\n\n");
        fflush(stdout);
}
else
{
	printf("\nlarge endian platform, no byteswapping needed\n\n");
        fflush(stdout);
}

if(argc != 2)
{
	printf("SDF browser syntax: %s <filename>\n",argv[0]);
	fflush(stdout);
	return 0;
}
strncpy (fname, argv[1], 99);

/* debugging the housekeeping data */

printf("Now checking SDF housekeeping data in file %s\n\n",fname);
sdf_housekeeping(fname,&ndat,&hdrpos,&datapos,&hdrsize);
printf("ndat = %d ",ndat);
output_int64(temp,hdrpos);
printf("hdrpos = %s\n",temp);
output_int64(temp,datapos);
printf("datapos = %s",temp);
output_int64(temp,hdrsize);
printf("hdrsize = %s\n\n",temp);
fflush(stdout);

top:  printf("Querying the SDF file: %s\n",fname);
fflush(stdout);

ndat=sdf_query(fname);

printf("Browse datasets by (o)rder or by (l)abel [ or (q)uit ]?\n");
fflush(stdout);
fscanf(stdin,"%2s",ans);

if(!strncmp(ans,"q",2))
{
	printf("quit entered, exiting SDF file browser\n");
	fflush(stdout);
	return 0;
}

if(!strncmp(ans,"l",2))
{
	printf("Enter string for matching with SDF labels:\n");
	fflush(stdout);
	fscanf(stdin,"%99s",labmatch); 
	
	matchind=sdf_labmatch(fname,&ndat,labmatch);
	
	if(matchind[0] == -1)
	{
		printf("no match to label %s found\n",labmatch);
		fflush(stdout);
		goto top;
	}
	else
	{
		for (i=0;i<1;i++) /* just first match done */
		{
			if(matchind[i] != -1)
			{
				printf("matchind = %d\n",matchind[i]);
		
				id=sdf_read(fname,matchind[i],(void **)&data);
	                        datasize=data_size(id);
	                        output_int64(temp,datasize);
	                        printf("datasize = %s\n",temp);
				fflush(stdout);
				printf("id->order = %d\n",id->order);
				printf("id->label = %s\n",id->label);
				printf("id->datatype = %c\n",id->datatype);
				printf("id->nbpw = %d\n",id->nbpw);
				printf("id->ndim = %d\n",id->ndim);
				for (ii=0;ii< id->ndim ;ii++)
				{
					output_int64(temp,*(id->dims+ii));	
					printf("id->dims[%d] = %s\n",ii,temp);
				}
				fflush(stdout);
				printf("Enter no. elements to print:\n");
				fflush(stdout);
				fscanf(stdin,"%99s",stuff);
				nelem=atopos(stuff);
                                if(nelem < (pos) 0)
                                {
                                  nelem = (pos) 0;
                                }
	/*
				output_int64(temp,nelem);
				printf("nelem = %s\n",temp);
				fflush(stdout);
	*/
	
				print_data(id,(void *)data,nelem); 
				/* 
				printf("got back from print_data OK\n");
				fflush(stdout);
                                */
				sdf_free(data);
				sdf_free_id(id);
				sdf_free(matchind);
				/*
				printf("finished freeing stuff OK\n");
				fflush(stdout);
				*/
				goto top;
			}
		}
	}
}

if(!strncmp(ans,"o",2))
{
	printf("Enter desired dataset order:\n");
	fflush(stdout);
	fscanf(stdin,"%d",&iorder);
	if((iorder < 0) || (iorder >= ndat))
	{
		printf("dataset order = %d does not exist in file %s\n",
			iorder,fname);
		fflush(stdout);
		goto top;
	}
	printf("selected dataset order = %d\n",iorder);
	fflush(stdout);
	id=sdf_read(fname,iorder,(void **)&data);
	datasize=data_size(id);
	output_int64(temp,datasize);
	printf("datasize = %s\n",temp);
	fflush(stdout);
	printf("id->order = %d\n",id->order);
	printf("id->label = %s\n",id->label);
	printf("id->datatype = %c\n",id->datatype);
	printf("id->nbpw = %d\n",id->nbpw);
	printf("id->ndim = %d\n",id->ndim);
	for (ii=0;ii< id->ndim ;ii++)
	{
		output_int64(temp,*(id->dims+ii));	
		printf("id->dims[%d] = %s\n",ii,temp);
	}
	fflush(stdout);
	printf("Enter no. elements to print:\n");
	fflush(stdout);
	fscanf(stdin,"%99s",stuff);
	nelem=atopos(stuff);
/*
			output_int64(temp,nelem);
			printf("nelem = %s\n",temp);
			fflush(stdout);
*/
	
	print_data(id,(void *)data,nelem); 
	sdf_free(data);
	sdf_free_id(id);
	goto top;

}

printf("Unrecognized response, exiting SDF file browser\n");
fflush(stdout);
return 1;

}

int print_data(data_id *id, char *data, pos nelem)
{
/*
    Print nelem elements of a chosen dataset.

    Create a new variable for use in a complicated case statement for
    printing out the data.  Call it idat, an i4 integer, whose values
    will be:
    idat = -1:  undefined/invalid
    idat = 0: byte data (print as hex)
    idat = 1: byte data (print as characters)
    idat = 2: short int data
    idat = 3: int data
    idat = 4: long long data
    idat = 5: float data
    idat = 6: complex data
    idat = 7: double data
    idat = 8: double complex data
*/
	pos datasize,nelem_max,i;
	char ans[20], temp[300];
	i4 idat, printcyes;
	f4 *f4arr;
	i8 *i8arr;
	i4 *i4arr;
	i2 *i2arr;
	f8 *f8arr;
	f4 *c4arr;
	f8 *c8arr;
	datasize=data_size(id);
	nelem_max = (nelem < datasize) ? nelem : datasize;
        printcyes=0;
        idat = -1;
        if((id->datatype == 'b') && (id->nbpw == 1))
	{
		printf("Print byte data as chars? (y or n) (if n then hex)\n");
		fflush(stdout);
		fscanf(stdin,"%2s",ans);
		if(!strncmp(ans,"y",2)) printcyes=1;
	}
        if((id->datatype == 'b') && (id->nbpw == 1) && (!printcyes)) idat = 0;
        if((id->datatype == 'b') && (id->nbpw == 1) && (printcyes)) idat=1;
        if((id->datatype == 'i') && (id->nbpw == 2)) idat = 2;
        if((id->datatype == 'i') && (id->nbpw == 4)) idat = 3;
        if((id->datatype == 'i') && (id->nbpw == 8)) idat = 4;
        if((id->datatype == 'f') && (id->nbpw == 4)) idat = 5;
        if((id->datatype == 'c') && (id->nbpw == 4)) idat = 6;
        if((id->datatype == 'f') && (id->nbpw == 8)) idat = 7;
        if((id->datatype == 'c') && (id->nbpw == 8)) idat = 8;
        
	switch (idat)
	{
		case -1:
			printf("undefined or invalid datatype, exiting\n");
			fflush(stdout);
			exit(1);
		case 0:
			for (i=0;i<nelem_max;i++)
			{
				printf("%02x:",*(data+i));
			}
			printf("\n");
			fflush(stdout);
			break;
		case 1:
			for (i=0;i<nelem_max;i++)
			{
				printf("%c",*(data+i));
			}
			printf("\n");
			fflush(stdout);
			break;

		case 2:
			
			i2arr = sdf_malloc(nelem_max* (pos) id->nbpw);
			memcpy((void *)i2arr, (void *)data, nelem_max*
				(pos) id->nbpw );
			for (i=0;i<nelem_max;i++)
			{
				printf("%d ",*(i2arr+i));
			}
			printf("\n");
			fflush(stdout);
			sdf_free(i2arr);
			break;

		case 3:
			
			i4arr = sdf_malloc(nelem_max* (pos) id->nbpw);
			memcpy((void *)i4arr, (void *)data, nelem_max*
				(pos) id->nbpw );
			for (i=0;i<nelem_max;i++)
			{
				printf("%d ",*(i4arr+i));
			}
			printf("\n");
			fflush(stdout);
			sdf_free(i4arr);
			break;

		case 4:
			
			i8arr = sdf_malloc(nelem_max* (pos) id->nbpw);
			memcpy((void *)i8arr, (void *)data, nelem_max*
				(pos) id->nbpw );
			for (i=0;i<nelem_max;i++)
			{
				output_int64(temp,*(i8arr+i));
				printf("%s",temp);
			}
			printf("\n");
			fflush(stdout);
			sdf_free(i8arr);
			break;

		case 5:
			
			f4arr = sdf_malloc(nelem_max* (pos) id->nbpw);
			memcpy((void *)f4arr, (void *)data, nelem_max*
				(pos) id->nbpw );
			for (i=0;i<nelem_max;i++)
			{
				printf("%g ",*(f4arr+i));
			}
			printf("\n");
			fflush(stdout);
			sdf_free(f4arr);
			break;

		case 6:
			
			c4arr = sdf_malloc(nelem_max* (pos) id->nbpw *(pos)2);
			memcpy((void *)c4arr, (void *)data, nelem_max*
				(pos) id->nbpw *(pos) 2);
			for (i=0;i<2*nelem_max;i+=2)
			{
				printf("(%g, %g)",*(c4arr+i),*(c4arr+i+1));
			}
			printf("\n");
			fflush(stdout);
			sdf_free(c4arr);
			break;

		case 7:
			
			f8arr = sdf_malloc(nelem_max* (pos) id->nbpw);
			memcpy((void *)f8arr, (void *)data, nelem_max*
				(pos) id->nbpw );
			for (i=0;i<nelem_max;i++)
			{
				printf("%g ",*(f8arr+i));
			}
			printf("\n");
			fflush(stdout);
			sdf_free(f8arr);
			break;

		case 8:
			
			c8arr = sdf_malloc(nelem_max* (pos) id->nbpw *(pos)2);
			memcpy((void *)c8arr, (void *)data, nelem_max*
				(pos) id->nbpw *(pos) 2);
			for (i=0;i<2*nelem_max;i+=2)
			{
				printf("(%g, %g)",*(c8arr+i),*(c8arr+i+1));
			}
			printf("\n");
			fflush(stdout);
			sdf_free(c8arr);
			break;

	}
	printf("(c)ontinue or (q)uit?\n");
	fflush(stdout);
	fscanf(stdin,"%2s",ans);
	if(!strncmp(ans,"c",2))
	{
		return 0;
	}
	else
	{
		printf("exiting SDF file browser\n");
		exit(0);
	}
return 0;
}
