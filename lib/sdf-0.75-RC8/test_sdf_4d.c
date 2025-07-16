#include <sdf_subs.h>
#include <math.h>

int main(int argc, char **argv);
int main(int argc, char **argv)
{

/* This is a test program illustrating how one can make an allocatable 4D
   floating point array in C and use sdf_write and sdf_read to write from,
   and read to, such arrays.  Also illustrates the use of sdf_transpose
   to convert from row-major to column major indexing */

	pos size,i,j,k,l;
        char temp[300];
        float pi=3.141592653;
	i4 szflt;
	data_id *id, *idcp, *id_in;
        float ****outdata;
        float ****indata;
	float *data;
        pos n1,n2,n3,n4;
	pos dims[20];
	i4 indorder[20];
	i4 reverse[20];
	char *outfile = "testout4d.sdf";
	szflt=sizeof(float);
        
	if(sizeof(pos) != 8)
	{
		printf("test_sdf_4d: sizeof(pos) = %d != 8, returning\n"
			,(int)sizeof(pos));
		fflush(stdout);
		return 1;
	}
	sdf_rm(outfile);
        n1=(pos)63;
        n2=(pos)64;
        n3=(pos)65;
	n4=(pos)66;
        dims[0]=n1;
	dims[1]=n2;
	dims[2]=n3;
	dims[3]=n4;
	id=sdf_create_id(0,"sinxyzt",'f',szflt,4,dims);

	/* make copy of id structure; will be necessary for freeing 4d array */
	idcp=sdf_cp_id(id);

        size=data_size(id);
	output_int64(temp,size);
        printf("size = %s\n",temp);
        fflush(stdout);
        outdata = (float ****)sdf_mk_4d(id);

	for (i=0;i<n1;i++)
	{
		for (j=0;j<n2;j++)
		{
			for (k=0;k<n3;k++)
			{
				for (l=0;l<n4;l++)
				{
					outdata[i][j][k][l]=
					sin(pi*(float)i/((float)(n1-1))) *
					sin(pi*(float)j/((float)(n2-1))) *
					sin(pi*(float)k/((float)(n3-1))) *
					sin(pi*(float)l/((float)(n4-1)));
				}
			}
		}
	}


/* If you want the array to look the same way in IDL as it does in C, 
   transpose it.  Arrays are stored in row-major order in C, 
   and column major in IDL. */

	indorder[0]=3;
	indorder[1]=2;
	indorder[2]=1;
	indorder[3]=0;
	reverse[0]=1;
	reverse[1]=1;
	reverse[2]=1;
	reverse[3]=1;
	sdf_transpose(indorder,reverse,id,outdata[0][0][0]);

/* transpose completed.  Note sdf_transpose does the transpose "in place". */

/* Loop to write out three of copies of this array as sdf datasets */

        for (i=0;i<3;i++)
	{
                sdf_write(outfile,id,(void *) outdata[0][0][0]);
		output_int64(temp,i);
                printf("test_sdf_4d: i = %s\n",temp);
		fflush(stdout);
	}

/* must use original structure, not transpose structure, to free 4d array,
   since the indexable array outdata didn't have its dimensions changed after 
   transpose */

	sdf_free_4d(idcp, (void ****) outdata);
	sdf_free_id(idcp);
        sdf_free_id(id);

/* OK, now let's read the dataset 1 in, and see if we can manipulate it like
   a 4d array.  First, read in the 1th dataset as though it were a 1d array: */

	id_in=sdf_read(outfile,1,(void **)&data);
/* Now convert it to a 4-d, indexable array: */
	indata=(float ****)sdf_1d_to_4d(id_in,(void *)data);

/* define n1, n2, n3, n4  in terms of id_in structure just for readability: */

	n1=id_in->dims[0];
	n2=id_in->dims[1];
	n3=id_in->dims[2];
	n4=id_in->dims[3];

	printf("test_sdf_4d: gonna print out one column of 4d array:\n\n");
	for (j=0;j<n2;j++)
	{
		printf("%g ",indata[n1/2][j][n3/2][n4/2]);
	}
	printf("\n");
	fflush(stdout);

/* clean up all the allocated stuff associated with the 4d array: */

	sdf_free_4d(id_in, (void ****)indata);
        sdf_free_id(id_in);
	return 0;
}
