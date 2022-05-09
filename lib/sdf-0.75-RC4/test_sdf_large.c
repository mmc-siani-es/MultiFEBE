#include <sdf_subs.h>
#include <math.h>

int main(int argc, char **argv);
int main(int argc, char **argv)
{

/* This is a test program, illustrating the use of sdf with an allocatable,
   2D floating point array in C.  It also illustrates the use of 
   sdf_transpose to switch from row major to column major indexing before 
   the array is written out to disk.

   Most importantly, this test program creates a 5GB file, to test whether
   large file support is working correctly.  Expect it to take about 10-20
   minutes to run to completion, with occasional long delays in execution as
   the OS tries to keep up with the disk I/O. Make sure you have enough disk
   space to run this before starting it, and make sure you are not on a FAT32
   partition! */

	pos size,i,j;
        char temp[300];
	char ans[100];
        float pi=3.141592653;
	i4 szflt;
	data_id *id, *idcp;
        float **indata;
        pos nx,ny;
	pos dims[20];
	i4 indorder[20];
	i4 reverse[20];
	char *outfile = "testoutlarge.sdf";

	szflt=sizeof(float); /* compute sizeof float */
	if(sizeof(pos) != 8)
	{
		printf("test_sdf_large: sizeof(pos) = %d != 8, returning\n"
			,(int)sizeof(pos));
		fflush(stdout);
		return 1;
	}
        nx=(pos)8192;
        ny=(pos)4096;
	dims[0]=nx;
	dims[1]=ny;
        /* Print out warning */
        printf("WARNING: this program will try to write a 5GB file!!\n");
        printf("(q) to quit, (c) to continue, then hit Enter:\n");
        fflush(stdout);
        fscanf(stdin,"%2s",ans);

        if(!strncmp(ans,"q",2))
        {
               printf("quit entered, exiting xlarge\n");
               fflush(stdout);
               return 0;
        }

	sdf_rm(outfile);

	id=sdf_create_id(0,"sinxy",'f',szflt,2,dims);
	idcp=sdf_cp_id(id);
/*      Next statement uses sdf_mk_2d to create a 2D indexable float array */
	indata = (float **)sdf_mk_2d(id);
	for (i=0;i<nx;i++)
	{
		for (j=0;j<ny;j++)
		{
			indata[i][j]=sin(pi*(float)i/((float)(nx-1))) *
				sin(pi*(float)j/((float)(ny-1)));
		}
	}
        size=data_size(id);
	output_int64(temp,size);
        printf("size = %s\n",temp);
        fflush(stdout);

/* If you want the array to look the same way in IDL as it does in C, 
   transpose it.  Arrays are stored in row-major order in C, 
   and column major in IDL. */

	indorder[0]=1;
	indorder[1]=0;
	reverse[0]=1;
	reverse[1]=1;
	sdf_transpose(indorder,reverse,id,indata[0]);

/* transpose completed.  Note sdf_transpose does the transpose "in place". */

/* Loop to write out a bunch of copies of this array as sdf datasets */

        for (i=0;i<38;i++)
	{
/*              Next stmt writes out the dynamically allocated 2d array */
                sdf_write(outfile,id,(void *) indata[0]);
		output_int64(temp,i);
                printf("test_sdf_large: i = %s\n",temp);
		fflush(stdout);
	}
/*      next statement frees the 2d array */
	sdf_free_2d(idcp, (void **)indata);

/*      Now free the id structures */
        sdf_free_id(idcp);
        sdf_free_id(id);
	return 0;
}
