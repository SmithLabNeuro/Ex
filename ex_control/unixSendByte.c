/*unixSendByte.C****************************************************************
 *
 * mex unixSendByte.c -lcomedi
 * (the -lrt isn't needed, although I don't know why not)
 *
 * Sends a 16-bit digital code to a measurement computing board,
 * specifically the PCIM-DAS1602.  Uses ports A & B for the
 * 16 bits, and sends a strobe on the first bit of port C.
 *
 *
 * 2022.01.26 - MAS - Set port bits to zero after write
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include <comedilib.h>
#include "mex.h"

#define SUBDEV 2 /* digital device (ports A/B/C) on the PCI-DAS1602 */

#define LOWPORT 0 /* 0-7 is port A */
#define HIGHPORT 8 /* 8-15 is port B */
#define STROBECHAN 16 /* first line of port C */

#define WRITEMASK 0xFFFF /* 16 1's - whole port open for writing */
#define LOWPORTMASK 0xFF00 /* 8 1's and 8 0's - low port open */
#define HIGHPORTMASK 0xFF /* 8 0's and 8 1's - high port open */

/* These variables set the lag between bit flips (ms) */
/* 20-50 uS is the range needed for the Ripple system (30 kHz) */
/* faster than that will likely cause dropped codes */
#define WAIT1_MICROSEC 0.050 /* between setting bits and the strobe */
#define WAIT2_MICROSEC 0.050 /* between setting the strobe to 1 and back to 0 */

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
    int i;
    double* doubleData;
    mxChar * charData;
    comedi_t *it;
    int flag;
    int stype;
    unsigned int byteToSend;
    
    it = comedi_open("/dev/comedi0");
    if (it == NULL)
    {
        mexPrintf("*** unixSendByte.c: Error opening /dev/comedi0\n");
        return;
    }
    
    stype = comedi_get_subdevice_type(it,SUBDEV);
    if(stype!=COMEDI_SUBD_DIO){
        mexPrintf("*** unixSendByte.c: Error, subdevice %d is not configurable for digital I/O\n",SUBDEV);
        return;
    }
    
    /* Configure the 3 ports for output */
    flag=comedi_dio_config(it,SUBDEV,LOWPORT,COMEDI_OUTPUT);
    if (flag==-1) {
        mexPrintf("*** unixSendByte.c: Error configuring low port for output\n");
        return;
    }
    flag=comedi_dio_config(it,SUBDEV,HIGHPORT,COMEDI_OUTPUT);
    if (flag==-1) {
        mexPrintf("*** unixSendByte.c: Error configuring high port for output\n");
        return;
    }
    flag=comedi_dio_config(it,SUBDEV,STROBECHAN,COMEDI_OUTPUT);
    if (flag==-1) {
        mexPrintf("*** unixSendByte.c: Error configuring strobe port for output\n");
        return;
    }
    
    if (nrhs > 0)
    {
        switch (mxGetClassID(prhs[0]))
        {
            case mxDOUBLE_CLASS:
                doubleData = mxGetPr(prhs[0]);
                
                for (i = 0; i < mxGetNumberOfElements(prhs[0]); i++)
                {
                    /* OLD CODE FROM WINDOWS USING CBW32.LIB */
                    /*ULStat = cbDOut(boardNum, lowPort, ((int)doubleData[i]) % 256);
                     * ULStat = cbDOut(boardNum, highPort, ((int)doubleData[i]) / 256);*/
                    
                    byteToSend = (int)doubleData[i];
                    flag=comedi_dio_bitfield(it,SUBDEV,WRITEMASK,&byteToSend);
                    
                    if (flag==-1) {
                        mexPrintf("*** unixSendByte.c: Error sending byte\n");
                        return;
                    }
                    
                    waitMS(WAIT1_MICROSEC);
                    flag = comedi_dio_write(it,SUBDEV,STROBECHAN,1);
                    if (flag==-1) {
                        mexPrintf("*** unixSendByte.c: Error turning strobe on\n");
                        return;
                    }
                    waitMS(WAIT2_MICROSEC);
                    flag = comedi_dio_write(it,SUBDEV,STROBECHAN,0);
                    if (flag==-1) {
                        mexPrintf("*** unixSendByte.c: Error turning strobe off\n");
                        return;
                    }
                }
                /* After we've written all the data, set all the port bits to zero */
                byteToSend = 0;
                flag=comedi_dio_bitfield(it,SUBDEV,WRITEMASK,&byteToSend);
                    
                break;
            case mxCHAR_CLASS:
                charData = (mxChar*)mxGetData(prhs[0]);
                
                for (i = 0; i < mxGetNumberOfElements(prhs[0]); i++)
                {
                    /* OLD CODE FROM WINDOWS USING CBW32.LIB */
                    /*ULStat = cbDOut(boardNum, lowPort, charData[i]); */
                    byteToSend = (int)charData[i];
                    flag=comedi_dio_bitfield(it,SUBDEV,WRITEMASK,&byteToSend);
                    if (flag==-1) {
                        mexPrintf("*** unixSendByte.c: Error sending byte\n");
                        return;
                    }
                    
                    waitMS(WAIT1_MICROSEC);
                    flag = comedi_dio_write(it,SUBDEV,STROBECHAN,1);
                    if (flag==-1) {
                        mexPrintf("*** unixSendByte.c: Error turning strobe on\n");
                        return;
                    }
                    waitMS(WAIT2_MICROSEC);
                    flag = comedi_dio_write(it,SUBDEV,STROBECHAN,0);
                    if (flag==-1) {
                        mexPrintf("*** unixSendByte.c: Error turning strobe off\n");
                        return;
                    }
                }
                /* After we've written all the data, set all the port bits to zero */
                byteToSend = 0;
                flag=comedi_dio_bitfield(it,SUBDEV,WRITEMASK,&byteToSend);                
                
                break;
            default:
                mexPrintf("Sorry, this data type cannot be transmitted.\n");
                mexEvalString("drawnow;");
        }
    }
    
    flag = comedi_close(it);
    if (flag==-1) {
        mexPrintf("*** unixSendByte.c: Error closing /dev/comedi0\n");
        return;
    }
    
    return;
    
}

