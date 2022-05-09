#include <sdf_subs.h>
#include <math.h>

int main(int argc, char **argv);
int main(int argc, char **argv)
{

/* This is a test program illustrating how one can make an allocatable 3D
   floating point array in C and use sdf_write and sdf_read to write from,
   and read to, such arrays.  Also illustrates the use of sdf_transpose
   to convert from row-major to column major indexing */

	pos size,i,j,k;
        char temp[300];
        float pi=3.141592653;
	i4 szflt;
	data_id *id, *idcp, *id_in;
        float ***outdata;
        float ***indata;
	float *data;
        pos nx,ny,nz;
	pos dims[20];
	i4 indorder[20];
	i4 reverse[20];
	char *outfile = "testout3d.sdf";
	szflt=sizeof(float);
        
	if(sizeof(pos) != 8)
	{
		printf("test_sdf_large: sizeof(pos) = %d != 8, returning\n"
			,(int)sizeof(pos));
		fflush(stdout);
		return 1;
	}
	sdf_rm(outfile);
        nx=(pos)64;
        ny=(pos)128;
        nz=(pos)256;
        dims[0]=nx;
	dims[1]=ny;
	dims[2]=nz;
	id=sdf_create_id(0,"sinxyz",'f',szflt,3,dims);

	/* make copy of id structure; will be necessary for freeing 3d array */
	idcp=sdf_cp_id(id);

        size=data_size(id);
	output_int64(temp,size);
        printf("size = %s\n",temp);
        fflush(stdout);
        outdata = (float ***)sdf_mk_3d(id);

	for (i=0;i<nx;i++)
	{
		for (j=0;j<ny;j++)
		{
			for (k=0;k<nz;k++)
			{
				outdata[i][j][k]=
					sin(pi*(float)i/((float)(nx-1))) *
					sin(pi*(float)j/((float)(ny-1))) *
					sin(pi*(float)k/((float)(nz-1)));
			}
		}
	}


/* If you want the array to look the same way in IDL as it does in C, 
   transpose it.  Arrays are stored in row-major order in C, 
   and column major in IDL. */

	indorder[0]=2;
	indorder[1]=1;
	indorder[2]=0;
	reverse[0]=1;
	reverse[1]=1;
	reverse[2]=1;
	sdf_transpose(indorder,reverse,id,outdata[0][0]); 

/* transpose completed.  Note sdf_transpose does the transpose "in place". */

/* Loop to write out three of copies of this array as sdf datasets */

        for (i=0;i<3;i++)
	{
                sdf_write(outfile,id,(void *) outdata[0][0]);
		output_int64(temp,i);
                printf("test_sdf_3d: i = %s\n",temp);
		fflush(stdout);
	}

/* must use original structure, not transpose structure, to free 3d array,
   since the indexable array outdata didn't have its dimensions changed after 
   transpose */

	sdf_free_3d(idcp, (void ***) outdata);
	sdf_free_id(idcp);
        sdf_free_id(id);

/* OK, now let's read the dataset 0 in, and see if we can manipulate it like
   a 3d array.  First, read in the 0th dataset as though it were a 1d array: */

	id_in=sdf_read(outfile,0,(void **)&data);
/* Now convert it to a 3-d, indexable array: */
	indata=(float ***)sdf_1d_to_3d(id_in,(void *)data);

/* define nx, ny, nz  in terms of id_in structure just for code readability: */

	nx=id_in->dims[0];
	ny=id_in->dims[1];
	nz=id_in->dims[2];

	printf("test_sdf_3d: gonna print out one column of 3d array:\n\n");
	for (j=0;j<ny;j++)
	{
		printf("%g ",indata[nx/2][j][nz/2]);
	}
	printf("\n");
	fflush(stdout);

/* clean up all the allocated stuff associated with the 3d array: */

	sdf_free_3d(id_in, (void ***)indata);
        sdf_free_id(id_in);
	return 0;
}
