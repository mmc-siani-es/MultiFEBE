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
