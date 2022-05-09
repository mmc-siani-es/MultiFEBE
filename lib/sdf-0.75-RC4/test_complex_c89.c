#include <sdf_subs.h>

int main (int argc, char **argv);
int main (int argc, char **argv)
{

/* In this program, we assume one does not have a C99 compliant compiler.
   Instead of reading in a 2-d complex float array, we will read in the array as
   a 3d float array, in which the last dimension is 2, and the value of the
   3rd index is 0 or 1 depending on whether one accesses the real or imaginary
   part.  This is kind of consistent with the convention used by fftw.

   This test program illustrates how one can read in a complex array, perform
   complex arithmetic, and write out a complex array even without having
   complex data types defined.  This is important for portability, as some C
   compilers (e.g. MSVC) are not C99 compliant and don't have complex variables
   defined.  This test program also illustrates the
   use of sdf_1d_to_3d and sdf_mk_3d functions to create indexable dynamically
   allocated arrays.  Note that because we can create dynamically allocated
   arrays with up to 5 dimensions, this allows us to have up to 4
   dimensions for complex arrays.
*/

pos nx,ny,i,j;
pos dims[20];
f4* data;
f4 ***cdata;
f4 ***cdataj;
f4 ***cdataa;
data_id *id=NULL, *idc=NULL, *ida=NULL;
if(sizeof(pos) !=8)
{
	printf("test_complex_c89: sizeof(pos) = %d != 8\n",(int)sizeof(pos));
	fflush(stdout);
	return 1;
}
/* read 2d complex array as dataset number 6 from sdf_example_file.sdf: */
id=sdf_read("sdf_example_file.sdf",6,(void **)&data);
for (i=0;i<id->ndim;i++)
{
	dims[i]=id->dims[i];
}
dims[id->ndim]=(pos)2;
if(id->datatype == 'c')
{
	/* idc is created to treat cdata as a 3d float array */
	idc=sdf_create_id(0,"cdata",'f',id->nbpw,(id->ndim)+1,dims);
}

/* for ida, must pretend like it's really complex, since dataa will be written
   out as though it were really a complex variable.  
   Use values from id not idc. */
ida=sdf_create_id(0,"cdataa",id->datatype,id->nbpw,id->ndim,id->dims);

/* convert 1d array read in with sdf_read to a 3d indexable float array, 
   last dim is real/imag part of each complex number */
cdata=(f4 ***)sdf_1d_to_3d(idc,(void *)data);

/* create a new 3d array that will be complex conjugate */
cdataj=(f4 ***)sdf_mk_3d(idc);

/* cdataa will be the orig. array times its complex conj., which in this case
   should be just 1.0 for all array elements */
cdataa=(f4 ***)sdf_mk_3d(idc);

nx=id->dims[0];
ny=id->dims[1];
for (i=0;i<nx;i++)
{
	for(j=0;j<ny;j++)
	{
		/* calc. compl. conj: */
		cdataj[i][j][0]=cdata[i][j][0]; /* real part */
		cdataj[i][j][1]=-cdata[i][j][1]; /* imag part */

		/* calc. cdata times its complex conjugate */
		cdataa[i][j][0]=cdata[i][j][0]*cdataj[i][j][0]
				-cdata[i][j][1]*cdataj[i][j][1]; /* real part */
		cdataa[i][j][1]=cdata[i][j][0]*cdataj[i][j][1]
				+cdata[i][j][1]*cdataj[i][j][0]; /* imag part */
	}
}
/* get rid of pre-existing copies of complex.sdf */
sdf_rm("complex89.sdf");

/* output the array times its complex conjugate.  Note that even though
  complex variables are not defined in C89, the array will be written
  out in the correct form because ida is set up to treat cdataa as complex */
sdf_write("complex89.sdf",ida,cdataa[0][0]);

/* free up the dynamically allocated arrays that were created by the program */
sdf_free_3d(id,(void ***)cdata);
sdf_free_3d(id,(void ***)cdataj);
sdf_free_3d(ida,(void ***)cdataa);

/* free up the id structures too */
sdf_free_id(ida);
sdf_free_id(idc);
sdf_free_id(id);

return 0;
}
