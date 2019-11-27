/*unixGetPupil.C****************************************************************

mex unixGetPupil.c -lcomedi

e.g., 
pos=unixGetPupil

Returns the current pupil diameter in terms of voltage, scaled to +/- 5V
 
***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <comedilib.h>
#include "mex.h"

#define SUBDEV 0 /* analog input (16 chan) on the PCI-DAS1602 */	
#define VOLTSCALE 20
#define VOLTOFFSET 10
#define MVPERVOLT 1000
#define EYEPUPILCHAN 2
#define RANGE 0
#define AREF AREF_GROUND
#define SETTLETIME 1000 /* in nanoseconds, 1000 means 1 microsecond */

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    comedi_t *it;
    int flag=0;
    lsampl_t eyex;
    int maxdata;
    double voltsx;
    int stype;
    double * out;
    
    
    it=comedi_open("/dev/comedi0");
    if (it == NULL)
    {
        mexPrintf("*** unixGetEyes.c: Error opening /dev/comedi0\n");
        return;
    }
    
    stype = comedi_get_subdevice_type(it,SUBDEV);
    if(stype!=COMEDI_SUBD_AI){
        mexPrintf("*** unixGetEyes.c: Error, subdevice %d is not configurable for analog input\n",SUBDEV);
        return;
    }
    
    maxdata=comedi_get_maxdata(it,SUBDEV,EYEPUPILCHAN);
    
    /*comedi_data_read(it,SUBDEV,EYEXCHAN,RANGE,AREF,&eyex);
    comedi_data_read(it,SUBDEV,EYEYCHAN,RANGE,AREF,&eyey);*/
    
    /* currently any number less than 1 microsecond in SETTLETIME gets rounded up to a microsecond */
    comedi_data_read_delayed(it,SUBDEV,EYEPUPILCHAN,RANGE,AREF,&eyex,SETTLETIME);

    /*int comedi_data_read_delayed(	comedi_t * device,
 	unsigned int subdevice,
 	unsigned int channel,
 	unsigned int range,
 	unsigned int aref,
 	lsampl_t * data,
 	unsigned int nanosec);
 
 int  comedi_data_read  (comedi_t  *  device,  unsigned  int  subdevice,
       unsigned int channel, unsigned int range, unsigned int aref, lsampl_t *
       data);*/
    
    /*mexPrintf("Eyes: %d %d\n",eyex,eyey);*/
    
    voltsx=comedi_to_phys(eyex,comedi_get_range(it,SUBDEV,EYEPUPILCHAN,RANGE),maxdata);
    
    /*mexPrintf("Volts: %g %g\n",voltsx,voltsy);*/

    plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    out = mxGetPr(plhs[0]);
    
    out[0] = (voltsx*VOLTSCALE-VOLTOFFSET)*MVPERVOLT; /* returns signal in mV */
    /*out[0] = voltsx;
    out[1] = voltsy;*/
    
    flag = comedi_close(it);
    if (flag==-1) {
        mexPrintf("*** unixGetEyes.c: Error closing /dev/comedi0\n");
        return;
    }
    
}
