#include <sdf_subs.h>
#include <math.h>

int main(int argc, char **argv);
int main(int argc, char **argv)
{

/* This is a test program illustrating how one can make an allocatable 5D
   floating point array in C and use sdf_write and sdf_read to write from,
   and read to, such arrays.  Also illustrates the use of sdf_transpose
   to convert from row-major to column major indexing */

	pos size,i,j,k,l,m;
        char temp[300];
        float pi=3.141592653;
	i4 szflt;
	data_id *id, *idcp, *id_in;
        float *****outdata;
        float *****indata;
	float *data;
        pos n1,n2,n3,n4,n5;
	pos dims[20];
	i4 indorder[20];
	i4 reverse[20];
	char *outfile = "testout5d.sdf";
	szflt=sizeof(float);
        
	if(sizeof(pos) != 8)
	{
		printf("test_sdf_5d: sizeof(pos) = %d != 8, returning\n"
			,(int)sizeof(pos));
		fflush(stdout);
		return 1;
	}
	sdf_rm(outfile);
        n1=(pos)10;
        n2=(pos)11;
        n3=(pos)12;
	n4=(pos)13;
	n5=(pos)14;
        dims[0]=n1;
	dims[1]=n2;
	dims[2]=n3;
	dims[3]=n4;
	dims[4]=n5;
	id=sdf_create_id(0,"sinxyzt",'f',szflt,5,dims);

	/* make copy of id structure; will be necessary for freeing 5d array */
	idcp=sdf_cp_id(id);

        size=data_size(id);
	output_int64(temp,size);
        printf("size = %s\n",temp);
        fflush(stdout);
        outdata = (float *****)sdf_mk_5d(id);

	for (i=0;i<n1;i++)
	{
		for (j=0;j<n2;j++)
		{
			for (k=0;k<n3;k++)
			{
				for (l=0;l<n4;l++)
				{
					for (m=0;m<n5;m++)
					{
					outdata[i][j][k][l][m]=
					sin(pi*(float)i/((float)(n1-1))) *
					sin(pi*(float)j/((float)(n2-1))) *
					sin(pi*(float)k/((float)(n3-1))) *
					sin(pi*(float)l/((float)(n4-1))) *
					sin(pi*(float)m/((float)(n5-1)));
					}
				}
			}
		}
	}


/* If you want the array to look the same way in IDL as it does in C, 
   transpose it.  Arrays are stored in row-major order in C, 
   and column major in IDL. */

	indorder[0]=4;
	indorder[1]=3;
	indorder[2]=2;
	indorder[3]=1;
	indorder[4]=0;
	reverse[0]=1;
	reverse[1]=1;
	reverse[2]=1;
	reverse[3]=1;
	reverse[4]=1;
	sdf_transpose(indorder,reverse,id,outdata[0][0][0][0]);

/* transpose completed.  Note sdf_transpose does the transpose "in place". */

/* Loop to write out three of copies of this array as sdf datasets */

        for (i=0;i<3;i++)
	{
                sdf_write(outfile,id,(void *) outdata[0][0][0][0]);
		output_int64(temp,i);
                printf("test_sdf_5d: i = %s\n",temp);
		fflush(stdout);
	}

/* must use original structure, not transpose structure, to free 4d array,
   since the indexable array outdata didn't have its dimensions changed after 
   transpose */

	sdf_free_5d(idcp, (void *****) outdata);
	sdf_free_id(idcp);
        sdf_free_id(id);

/* OK, now let's read the dataset 1 in, and see if we can manipulate it like
   a 4d array.  First, read in the 1th dataset as though it were a 1d array: */

	id_in=sdf_read(outfile,1,(void **)&data);
/* Now convert it to a 5-d, indexable array: */
	indata=(f4 *****)sdf_1d_to_5d(id_in,(void *)data);

/* define n1, n2, n3, n4, n5  in terms of id_in structure just 
   for readability: */

	n1=id_in->dims[0];
	n2=id_in->dims[1];
	n3=id_in->dims[2];
	n4=id_in->dims[3];
	n5=id_in->dims[4];

	printf("test_sdf_5d: gonna print out one column of 5d array:\n\n");
	for (j=0;j<n2;j++)
	{
		printf("%g ",indata[n1/2][j][n3/2][n4/2][n5/2]);
	}
	printf("\n");
	fflush(stdout);

/* clean up all the allocated stuff associated with the 5d array: */

	sdf_free_5d(id_in, (void *****)indata);
        sdf_free_id(id_in);
	return 0;
}
