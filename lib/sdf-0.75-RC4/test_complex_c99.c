#include <sdf_subs.h>
#include <complex.h>

int main (int argc, void **argv);
int main (int argc, void **argv)
{

/* Note that in C, unlike Fortran, the default complex variable is 
   double precision, not single precision.  Thus to read a single precision
   complex number, must declare it as complex float instead of just complex.

   This test program reads a complex array from sdf_example_file.sdf */

/* Note this will only work with a C99 compliant compiler.  See 
   test_complex_c89.c for an example using C89. */

pos nx,ny,i,j;
complex float* data;
complex float **cdata;
complex float **cdataj;
complex float **cdataa;
data_id *id, *ida;
/* read in the complex array as dataset number 6 from test_example_file.sdf: */
id=sdf_read("sdf_example_file.sdf",6,(void **)&data);
ida=sdf_create_id(0,"cdataa",id->datatype,id->nbpw,id->ndim,id->dims);

/* convert the 1-d complex array to a 2-d indexable array */
cdata=(complex float **)sdf_1d_to_2d(id,(void *)data);

/* create a new 2-d complex array that will be complex conjugate */
cdataj=(complex float **)sdf_mk_2d(id);

/* cdataa will be the array times its complex conj., which in this case
   should be just 1.0 for all array elements */
cdataa=(complex float **)sdf_mk_2d(id);
nx=id->dims[0];
ny=id->dims[1];
for (i=0;i<nx;i++)
{
	for(j=0;j<ny;j++)
	{
		/* calc. compl. conj. and array times compl. conj: */

                /* if using std complex variable in C instead of single 
                   precision, use conj instead of conjf */
		cdataj[i][j]=conjf(cdata[i][j]); 
		cdataa[i][j]=cdata[i][j]*cdataj[i][j];
	}
}
sdf_rm("complex99.sdf");
/* output the array times its complex conjugate */
sdf_write("complex99.sdf",ida,cdataa[0]);

sdf_free_2d(id,(void **)cdata);
sdf_free_2d(id,(void **)cdataj);
sdf_free_2d(ida,(void **)cdataa);
sdf_free_id(ida);
sdf_free_id(id);

return 0;
}
