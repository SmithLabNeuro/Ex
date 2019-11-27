/*unixSendPulse.C****************************************************************

mex unixSendPulse.c -lcomedi
(the -lrt isn't needed, although I don't know why not)

if no input is passed, sends a '1' for 20 ms on channel 18 (port C2)

16 is the strobe!!!!!

Can accept as input the channel and duration (in milliseconds)

e.g., 
unixSendPulse(18)
unixSendPulse(18,100)

***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <comedilib.h>
#include "mex.h"

#define SUBDEV 2 /* digital device (ports A/B/C) on the PCI-DAS1602 */	

void waitMS(float ms)
{
  struct timespec start,current;
  double elapsed=0.0;
  double waitSec=ms/1000;

  clock_gettime(CLOCK_REALTIME, &start);

  while (elapsed < waitSec) {
    clock_gettime(CLOCK_REALTIME, &current);
    elapsed = current.tv_sec + current.tv_nsec/1E9;
    elapsed -= start.tv_sec + start.tv_nsec/1E9;
  }
  /*mexPrintf("CLOCK_REALTIME: %f\n", elapsed);*/

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   float duration = 20.0;
   comedi_t *it;
   int flag;
   int chan = 18; /* default to 3rd bit on port C (16-23) - the juice line */
   int stype;

   if (nrhs > 0) {
     chan = (int)*mxGetPr(prhs[0]);
     if (nrhs > 1) {
       duration = *mxGetPr(prhs[1]);
     }
   }

   it = comedi_open("/dev/comedi0");
   if (it == NULL)
     {
       mexPrintf("*** unixSendPulse.c: Error opening /dev/comedi0\n");
       return;
     }
   
   stype = comedi_get_subdevice_type(it,SUBDEV);
   if(stype!=COMEDI_SUBD_DIO){
     mexPrintf("*** unixSendPulse.c: Error, subdevice %d is not configurable for digital I/O\n",SUBDEV);
     return;
   }
   
   flag=comedi_dio_config(it,SUBDEV,chan,COMEDI_OUTPUT);
   if (flag==-1) {
     mexPrintf("*** unixSendPulse.c: Error configuring channel %d for output\n",chan);
     return;
   }

   flag = comedi_dio_write(it,SUBDEV,chan,1);
   if (flag==-1) {
     mexPrintf("*** unixSendPulse.c: Error sending setting channel %d to high\n",chan);
     return;
   }
   
   /* wait duration ms before turning the pulse off */  
   waitMS(duration);
   
   flag = comedi_dio_write(it,SUBDEV,chan,0);
   if (flag==-1) {
     mexPrintf("*** unixSendPulse.c: Error sending setting channel %d to low\n",chan);
     return;
   }

   flag = comedi_close(it);
   if (flag==-1) {
     mexPrintf("*** unixSendPulse.c: Error closing /dev/comedi0\n");
     return;
   }
   
   return;
}
