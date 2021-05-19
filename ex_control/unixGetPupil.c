/*unixGetPupil.C****************************************************************

mex unixGetPupil.c -lcomedi

e.g., 
pos=unixGetPupil

Returns the current pupil diameter in terms of voltage, scaled to +/- 5V

2021.03.19 - updated version to handle changes in the I/O card driver

***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <comedilib.h>
#include "mex.h"

#define SUBDEV 0 /* analog input (16 chan) on the PCI-DAS1602 */	
#define VOLTOFFSET 10
#define MVPERVOLT 1000
#define PUPILCHAN 2
#define RANGE 0
#define AREF AREF_GROUND
#define SETTLETIME 1000 /* in nanoseconds, 1000 means 1 microsecond */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    comedi_t *it;
    int flag=0;
    lsampl_t pupilval;
    int maxdata;
    double voltsp;
    int stype;
    double * out;
    static comedi_range fixed_range =
    {
	.min = -1.0 * VOLTOFFSET,
	.max = 1.0 * VOLTOFFSET,
	.unit = UNIT_volt
    };
    comedi_range *rng;
    
    it=comedi_open("/dev/comedi0");
    if (it == NULL)
    {
        mexPrintf("*** unixGetPupil.c: Error opening /dev/comedi0\n");
        return;
    }
    
    stype = comedi_get_subdevice_type(it,SUBDEV);
    if(stype!=COMEDI_SUBD_AI){
        mexPrintf("*** unixGetPupil.c: Error, subdevice %d is not configurable for analog input\n",SUBDEV);
        return;
    }
    
    maxdata=comedi_get_maxdata(it,SUBDEV,PUPILCHAN);
    
    /* currently any number less than 1 microsecond in SETTLETIME gets rounded up to a microsecond */
    comedi_data_read_delayed(it,SUBDEV,PUPILCHAN,RANGE,AREF,&pupilval,SETTLETIME);
    
    rng = comedi_get_range(it,SUBDEV,PUPILCHAN,RANGE);
    if (rng->unit == UNIT_none)
    {
	/* Replace incorrect range reported by old cb_pcimdas driver. */
	rng = &fixed_range;
    }
    voltsp=comedi_to_phys(pupilval,rng,maxdata);

    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    out = mxGetPr(plhs[0]);
    
    out[0] = voltsp*MVPERVOLT; /* returns signal in mV */
    
    flag = comedi_close(it);
    if (flag==-1) {
        mexPrintf("*** unixGetPupil.c: Error closing /dev/comedi0\n");
        return;
    }
    
}
