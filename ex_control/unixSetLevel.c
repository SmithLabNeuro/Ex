/*unixSetLevel.c****************************************************************

mex unixSetLevel.c -lcomedi
(the -lrt isn't needed, although I don't know why not)

Set the level high or low and leave it that way

18 is juice, 19-22 are the digital outputs
use a 1 to set the value high, and 0 to set it low

e.g.,
unixSetLevel(19,1)
unixSetLevel(19,0)

***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <comedilib.h>
#include "mex.h"

#define SUBDEV 2 /* digital device (ports A/B/C) on the PCI-DAS1602 */	

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   comedi_t *it;
   int flag;
   int chan = 19; /* default to 3rd bit on port C (16-23) - the juice line */
   int digval = 0;
   int stype;

   if (nrhs < 2) {
       mexPrintf("*** unixSetLevel.c: Error, must pass a channel and a value\n");
       return;
   }

   chan = (int)*mxGetPr(prhs[0]);
   digval = (int)*mxGetPr(prhs[1]);
   
   if (digval != 0 && digval != 1){
     mexPrintf("*** unixSetLevel.c: Error, must pass a value of 0 or 1\n");
     return;
   }
   
   it = comedi_open("/dev/comedi0");
   if (it == NULL)
     {
       mexPrintf("*** unixSetLevel.c: Error opening /dev/comedi0\n");
       return;
     }
   
   stype = comedi_get_subdevice_type(it,SUBDEV);
   if(stype!=COMEDI_SUBD_DIO){
     mexPrintf("*** unixSetLevel.c: Error, subdevice %d is not configurable for digital I/O\n",SUBDEV);
     return;
   }
   
   flag=comedi_dio_config(it,SUBDEV,chan,COMEDI_OUTPUT);
   if (flag==-1) {
     mexPrintf("*** unixSetLevel.c: Error configuring channel %d for output\n",chan);
     return;
   }

   flag = comedi_dio_write(it,SUBDEV,chan,digval);
   if (flag==-1) {
     mexPrintf("*** unixSetLevel.c: Error setting channel %d to %d\n",chan,digval);
     return;
   }
   
   flag = comedi_close(it);
   if (flag==-1) {
     mexPrintf("*** unixSetLevel.c: Error closing /dev/comedi0\n");
     return;
   }
   
   return;
}
