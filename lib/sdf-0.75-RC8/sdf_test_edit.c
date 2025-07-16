#include <sdf_subs.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/* function prototypes not in sdf_subs.h */
int main(int argc, char **argv);

int main(int argc, char **argv)
{

/* This simple test program will delete, insert, or replace datasets in an
sdf file.  It is intended only to provide simple tests of sdf_delete,
sdf_insert, and sdf_replace done from a C calling program.  Here, for 
insertion/replacement, the new dataset is a 101
length floating point array equal to sin (pi x), with x ranging from 0 to 1. */

data_id * id;
char temp[301],ans[100];
char fname[100];
char labmatch[100];
char *labelnew = "sinx";
char dtype='f';
i4 i,iorder,ise,ibe,ndat;
i4 nbpw, ndim;
i4 *matchind;
pos hdrpos, datapos, hdrsize;
pos dims[20];
float pi=3.141592653;
float datanew[101];
for (i=0;i<101;i++)
{
	datanew[i]=sin(0.01*pi*(float)i);
}
test_sizes();
nbpw=4;
ndim=1;
dims[0]=(pos)101;
ise=1;
ibe=is_big_endian();
if(ibe) ise=0;
if(ise)
{
	printf("\nsmall endian platform detected, byteswapping enabled\n\n");
        fflush(stdout);
}
else
{
	printf("\nlarge endian platform detected, no byteswapping needed\n\n");
        fflush(stdout);
}

if(argc != 2)
{
	printf("SDF test editor syntax: %s <filename>\n",argv[0]);
	fflush(stdout);
	return 0;
}
strncpy (fname, argv[1], 99);

/* checking the housekeeping data */

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

printf("Edit datasets by (o)rder or by (l)abel [ or (q)uit ]?\n");
fflush(stdout);
fscanf(stdin,"%2s",ans);

if(!strncmp(ans,"q",2))
{
	printf("quit entered, exiting SDF test file editor\n");
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
	/* add logic to ask about insert, delete, or replace */
	printf("(i)nsert, (d)elete, (r)eplace, or (q)uit?\n");
	fflush(stdout);
	fscanf(stdin,"%2s",ans);
	if(!strncmp(ans,"d",2))
	{
		sdf_delete(fname,matchind[0]);
	}
	if(!strncmp(ans,"i",2))
	{
		id=sdf_create_id(matchind[0],labelnew,dtype,nbpw,ndim,dims);
		sdf_insert(fname,matchind[0],id,datanew);
		sdf_free_id(id);
	}
	if(!strncmp(ans,"r",2))
	{
		id=sdf_create_id(matchind[0],labelnew,dtype,nbpw,ndim,dims);
		sdf_replace(fname,matchind[0],id,datanew);
		sdf_free_id(id);
	}
	if(!strncmp(ans,"q",2))
	{
		printf("quit entered, exiting sdf test editor\n");
		fflush(stdout);
		return 0;
	}
	{
		for (i=0;i<ndat;i++)
		{
			if(matchind[i] != -1)
			{
	/*
				output_int64(temp,nelem);
				printf("nelem = %s\n",temp);
				fflush(stdout);
	*/
	
				sdf_free(matchind);
				/* printf("finished freeing stuff OK\n"); */
				fflush(stdout);
			        /* bug - can only be in this loop for i=0 */
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
	if((iorder < 0) || (iorder >= ndat+1))
	{
		printf("dataset order = %d does not exist in file %s\n",
			iorder,fname);
		fflush(stdout);
		goto top;
	}
	printf("selected dataset order = %d\n",iorder);
	fflush(stdout);
	/* logic for edit */
	printf("(i)nsert, (d)elete, (r)eplace, or (q)uit?\n");
	fflush(stdout);
	fscanf(stdin,"%2s",ans);
	if(!strncmp(ans,"d",2))
	{
		sdf_delete(fname,iorder);
	}
	if(!strncmp(ans,"i",2))
	{
		id=sdf_create_id(iorder,labelnew,dtype,nbpw,ndim,dims);
		sdf_insert(fname,iorder,id,datanew);
		sdf_free_id(id);
	}
	if(!strncmp(ans,"r",2))
	{
		id=sdf_create_id(iorder,labelnew,dtype,nbpw,ndim,dims);
		sdf_replace(fname,iorder,id,datanew);
		sdf_free_id(id);
	}
	if(!strncmp(ans,"q",2))
	{
		printf("quit entered, exiting sdf test editor\n");
		fflush(stdout);
		return 0;
	}
/*
			output_int64(temp,nelem);
			printf("nelem = %s\n",temp);
			fflush(stdout);
*/
	
	goto top;

}

printf("Unrecognized response, exiting test SDF file editor\n");
fflush(stdout);
return 1;

}

