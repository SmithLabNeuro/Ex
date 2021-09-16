/*unixGetAnalogInput.C****************************************************************

mex unixGetAnalogInput.c -lcomedi

e.g., 
chInX = 0; % often the x-channel for the eye tracker
chInY = 1; % often the y-channel for the eye tracker
posX=unixGetAnalogInput(chInX)
posY=unixGetAnalogInput(chInY)

Returns the analog voltage on the input channel, scaled to +/- 5V

2021.03.19 - updated version to handle changes in the I/O card driver

***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <comedilib.h>
#include "mex.h"

#define SUBDEV 0 /* analog input (16 chan) on the PCIe-DAS1602/16 */	
#define VOLTOFFSET 10
#define MVPERVOLT 1000
#define RANGE 0
#define AREF AREF_GROUND
#define SETTLETIME 1000 /* in nanoseconds, 1000 means 1 microsecond */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    comedi_t *it;
    int flag=0;
    lsampl_t analogval;
    int maxdata;
    double voltsp;
    int stype;
    int analogCh;
    double * out;
    static comedi_range fixed_range =
    {
	.min = -1.0 * VOLTOFFSET,
	.max = 1.0 * VOLTOFFSET,
	.unit = UNIT_volt
    };
    comedi_range *rng;
    
    if (nrhs < 1) {
       mexPrintf("*** unixGetAnalogInput.c: Error, must pass the analog channel on the PCIe-DAS1602-16 DAQ being grabbed\n");
       return;
    }
    else
    {
        analogCh = (int)*mxGetPr(prhs[0]);
    }
    
    it=comedi_open("/dev/comedi0");
    if (it == NULL)
    {
        mexPrintf("*** unixGetAnalogInput.c: Error opening /dev/comedi0\n");
        return;
    }
    
    stype = comedi_get_subdevice_type(it,SUBDEV);
    if(stype!=COMEDI_SUBD_AI){
        mexPrintf("*** unixGetAnalogInput.c: Error, subdevice %d is not configurable for analog input\n",SUBDEV);
        return;
    }
    
    maxdata=comedi_get_maxdata(it,SUBDEV,analogCh);
    
    /* currently any number less than 1 microsecond in SETTLETIME gets rounded up to a microsecond */
    comedi_data_read_delayed(it,SUBDEV,analogCh,RANGE,AREF,&analogval,SETTLETIME);
    
    rng = comedi_get_range(it,SUBDEV,analogCh,RANGE);
    if (rng->unit == UNIT_none)
    {
	/* Replace incorrect range reported by old cb_pcimdas driver. */
	rng = &fixed_range;
    }
    voltsp=comedi_to_phys(analogval,rng,maxdata);

    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    out = mxGetPr(plhs[0]);
    
    out[0] = voltsp*MVPERVOLT; /* returns signal in mV */
    
    flag = comedi_close(it);
    if (flag==-1) {
        mexPrintf("*** unixGetPupil.c: Error closing /dev/comedi0\n");
        return;
    }
    
}
