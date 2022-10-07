/* MATLABUDP2.C
 *
 *	MATLABUDP2.c contains a few c-routines
 *	to be called from MATLAB so that machines can chat via
 *	ethernet and the UDP/IP protocols.
 *
 *  This code was originally written by Ben Heasly in Josh Gold's lab at U Penn (circa 2006)
 *
 *  Slight changes were made over the years in the Smith Lab, and the name
 *  was changed to MATLABUDP2 by Ryan Williamson when he modified it to allow
 *  communication over more than one UDP port at a time. 
 *
 */

#include "matlabUDP2.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    char *command=NULL;
    int buf_len;
    int socketindex;
    /* If no arguments given, print usage string */
    if(nrhs < 1) {
        mexPrintf("matlabUDP usage:\n socketindex = matlabUDP('open', (string)localIP, (string)remoteIP, (int)port);%% should return small int\n matlabUDP('send', socketindex, (string)message);\n messageIsAvailable = matlabUDP('check', socketindex);\n message = matlabUDP('receive', socketindex);\n socketIsOpen = matlabUDP('close', socketindex);%% should return -1\n");
        return;
    }
    
    /* First argument is command string... get and convert to char * */
    if(mxGetM(prhs[0]) == 1 && mxGetN(prhs[0]) >= 1 && mxIsChar(prhs[0])) {
        buf_len =  mxGetN(prhs[0]) + 1;
        command = mxCalloc(buf_len, sizeof(char));
        if(mxGetString(prhs[0], command, buf_len))
            mexWarnMsgTxt("matlabUDP: Not enough heap space. String (command) is truncated.");
    } else {
        mexErrMsgTxt("matlabUDP: First argument should be a string (command).");
    }
    
    /* case on command string... */
    if(!strncmp(command, "open", 4)) {
        
        /* done with command */
        mxFree(command);
        
        /* Modify this code to eliminate socket in a sensible way*/
        /* register exit routine to free socket */
        /*if(mexAtExit(mat_UDP_close) != 0 ) {
            mat_UDP_close();
            mexErrMsgTxt("matlabUDP: failed to register exit routine, mat_UDP_close.");
        }*/
        
        /* only open a fresh socket if */
        /*  PORT arg is a number, and */
        /*  IP addr args are short strings e.g. "111.222.333.444" */
        if(nrhs==4 && mxIsNumeric(prhs[3])
        && mxIsChar(prhs[2]) && mxGetN(prhs[2])<=15
        && mxIsChar(prhs[1]) && mxGetN(prhs[1])<=15){

            /* close old socket? */
            /*if(mat_UDP_sockfd>=0)
               mat_UDP_close();*/
            
            /* format args for socket opener function */
            char local[16],remote[16];
            mxGetString(prhs[1],local,16);
            mxGetString(prhs[2],remote,16);

            /* openerup */
            mexPrintf("matlabUDP opening socket\n");
            socketindex = mat_UDP_open(local, remote, (int)mxGetScalar(prhs[3]));
            
        }
        
        /* always return socket index */

        /* What is msCreateNumericArray???*/
        /* build me a return value worthy of MATLAB */
        if(!(plhs[0] = mxCreateDoubleScalar(socketindex)))
            mexErrMsgTxt("matlabUDP: mxCreateNumericArray failed.");
        
        
    } else if(!strncmp(command, "receive", 4)) {
        
        /* done with command */
        mxFree(command);
        
        int i;
        mwSize dims[2];
        unsigned short *outPtr;
        socketindex = (int)mxGetScalar(prhs[1]);
        dims[0] = 1;
        
        if(nlhs<=1){
            
            if(mat_UDP_sockfd[socketindex]<0){

                /* socket closed so zero bytes are read */
                i = 0;

            } else {
                
                /* read new bytes from socket */
                mat_UDP_read(socketindex, &mat_UDP_messBuff[socketindex][0], MAX_NUM_BYTES); /* sets mat_UDP_numBytes */
                i = mat_UDP_numBytes[socketindex];

            }

            /* always provide at least an empty return value */
//             mexPrintf("num Bytes: %d\n", i);
            dims[1] = i;
            const mwSize finalDims[2] = { dims[0],  dims[1]};
            if(!(plhs[0] = mxCreateCharArray(2, finalDims)))
                mexErrMsgTxt("matlabUDP: mxCreateCharArray failed.");
            
            /* fill in report with any new bytes */
            outPtr = (unsigned short *) mxGetData(plhs[0]);
            for(i--; i>=0; i--){
                *(outPtr + i) = mat_UDP_messBuff[socketindex][i];
            }

        }

        
    } else if(!strncmp(command, "send", 4)) {
        
        /* done with command */
        mxFree(command);
        
        socketindex = (int)mxGetScalar(prhs[1]);
        
        if(mat_UDP_sockfd[socketindex]<0){
            /* warn that no message was not sent */
            mexWarnMsgTxt("matlabUDP: Message not sent.  No socket is open.");
            

            
        } else {
            
            /* only send message if message arg is a 1-by-N char array */
            if(nrhs==3 && (mxIsChar(prhs[2]) || mxIsUint8(prhs[2])) && mxGetM(prhs[2])==1 && mxGetN(prhs[2])>0){
                
                /* format ye string and send forth */
                mxClassID classID = mxGetClassID(prhs[2]);
                switch (classID) {
                    case mxUINT8_CLASS:
                    {
                        buf_len = mxGetN(prhs[2]);
                        unsigned char *uintArrayMsg;
                        uintArrayMsg = (unsigned char*)mxGetUint8s(prhs[2]);
                        memcpy(&mat_UDP_messBuff[socketindex][0], uintArrayMsg, mxGetN(prhs[2])*sizeof(unsigned char));
                        break;
                    }
                    case mxCHAR_CLASS:
                    {
                        mxGetString(prhs[2],&mat_UDP_messBuff[socketindex][0],mxGetN(prhs[2])+1);
                        break;
                    }
                }
                /*mexPrintf("sending...\n");*/
                mat_UDP_send(socketindex, &mat_UDP_messBuff[socketindex][0], mxGetN(prhs[2]));
                
            }else{
                
                /* warn that no message was not sent */
                mexWarnMsgTxt("matlabUDP: Message not sent.  Must be 1-by-N char array.");
                
            }
        }

        
    } else if(!strncmp(command, "check", 4)) {
        
        /* done with command */
        mxFree(command);
        
        socketindex = (int)mxGetScalar(prhs[1]);
        
        /* always provide a return value */
        /* if socket is closed, && will short-circuit and skip the actual socket check */
        if(!(plhs[0] = mxCreateDoubleScalar( (double) (mat_UDP_sockfd[socketindex]>=0) && mat_UDP_check(socketindex) )))
            mexErrMsgTxt("matlabUDP: mxCreateNumericArray failed.");
        
        
    } else if(!strncmp(command, "close", 4)) {
        
        
        /* done with command */
        mxFree(command);

        socketindex = (int)mxGetScalar(prhs[1]);
        
        /* only try to close if socket is open */
        if(mat_UDP_sockfd[socketindex] >= 0)
            mat_UDP_close(socketindex);
        
        /* always return socket index */
        if(nlhs==1){
            if(!(plhs[0] = mxCreateDoubleScalar((double)mat_UDP_sockfd[socketindex])))
                mexErrMsgTxt("matlabUDP: mxCreateNumericArray failed.");
        }
    } else if(!strncmp(command, "all_close", 4)) {


        /* done with command */
        mxFree(command);

        /* only try to close if socket is open */
        for(socketindex = 0; socketindex<MAX_SOCKETS;socketindex++){
            if(mat_UDP_sockfd[socketindex] >= 0)
                mat_UDP_close(socketindex);
        }  

    } else {
        
        /* done with command */
        mxFree(command);
        
        mexWarnMsgTxt("matlabUDP: Unknown command option");
    }
}

/* initialize UDP socket */
int mat_UDP_open (char localIP[], char remoteIP[], int port){
    
    int sockind = nextSocket();
    
    
    mat_UDP_REMOTE_addr[sockind].sin_family = AF_INET;	/* host byte order */
    mat_UDP_REMOTE_addr[sockind].sin_port = htons(port);	/* short, network byte order */
    mat_UDP_REMOTE_addr[sockind].sin_addr.s_addr = inet_addr(remoteIP);
    memset(&(mat_UDP_REMOTE_addr[sockind].sin_zero), '\0', 8);/* zero the rest of the struct */
    
    mat_UDP_LOCAL_addr[sockind].sin_family = AF_INET;         /* host byte order */
    mat_UDP_LOCAL_addr[sockind].sin_port = htons(port);     /* short, network byte order */
    mat_UDP_LOCAL_addr[sockind].sin_addr.s_addr = inet_addr(localIP);
    memset(&(mat_UDP_LOCAL_addr[sockind].sin_zero), '\0', 8); /* zero the rest of the struct */
    
    /* mexPrintf("localIP = <%s>\n",inet_ntoa(mat_UDP_LOCAL_addr.sin_addr)); */
    /* mexPrintf("remoteIP = <%s>\n",inet_ntoa(mat_UDP_REMOTE_addr.sin_addr)); */
    /* mexPrintf("ports = <%i>,<%i>\n",mat_UDP_LOCAL_addr.sin_port,mat_UDP_REMOTE_addr.sin_port  ); */
    
    if ((mat_UDP_sockfd[sockind]=socket(AF_INET, SOCK_DGRAM, 0)) == -1) {
        mexErrMsgTxt("Couldn't create UDP socket.");
    }
    
    /* mexPrintf("sockFD = %i\n",mat_UDP_sockfd); */
    
    if (bind(mat_UDP_sockfd[sockind], (struct sockaddr *)&(mat_UDP_LOCAL_addr[sockind]), mat_UDP_addr_len) == -1){
        mexErrMsgTxt("Couldn't bind socket.  Maybe invalid local address.");
    }
    
    return sockind;
}


/* send a string to MATLAB */
void mat_UDP_send (int sockind, char mBuff[], int mLen){
    
    /*     const char drPhil[] = {"you're a loser"}; */
    /*     int callyourwifeabitch = 11; */
    /*     int youreabigfatgooneybird = 0; */
    /* */
    /*     youreabigfatgooneybird=sendto(mat_UDP_sockfd, drPhil, callyourwifeabitch, MSG_DONTWAIT,(struct sockaddr *)&mat_UDP_REMOTE_addr, mat_UDP_addr_len); */
    /*     mexPrintf("loser=<%s>, loserLen=%i, retVal=%i\n",drPhil,callyourwifeabitch,youreabigfatgooneybird); */
    
    if ((mLen=sendto(mat_UDP_sockfd[sockind], mBuff, mLen, MSG_DONTWAIT,
    (struct sockaddr *)&(mat_UDP_REMOTE_addr[sockind]), mat_UDP_addr_len)) == -1)
        mexWarnMsgTxt("Couldn't send string.  Are computers connected??");
}


/* is a return message available? */
int mat_UDP_check (int sockind){
    static struct timeval timout[MAX_SOCKETS];
    static fd_set readfds[MAX_SOCKETS];
    FD_ZERO(&(readfds[sockind]));
    FD_SET(mat_UDP_sockfd[sockind],&(readfds[sockind]));
    select(mat_UDP_sockfd[sockind]+1,&(readfds[sockind]),NULL,NULL,&(timout[sockind]));
    return(FD_ISSET(mat_UDP_sockfd[sockind],&(readfds[sockind])));
}


/* read any available message */
void mat_UDP_read (int sockind, char mBuff[], int messUpToLen){
    
    if ((mat_UDP_numBytes[sockind]=recvfrom(mat_UDP_sockfd[sockind],mBuff, messUpToLen, MSG_DONTWAIT,
    (struct sockaddr *)&(mat_UDP_REMOTE_addr[sockind]), &mat_UDP_addr_len)) <0 )
        mat_UDP_numBytes[sockind]=0;
    
}

/*overload this function so we can close all sockets in one command*/
/* cleanup UDP socket */
void mat_UDP_close (int sockind){
    if(mat_UDP_sockfd[sockind]>=0){
        mexPrintf("matlabUDP closing socket\n");
        close(mat_UDP_sockfd[sockind]);
        mat_UDP_sockfd[sockind]=-1;
    }
}

int nextSocket (void){
    int n;
    for(n = 0; n<MAX_SOCKETS; n++){
        if(mat_UDP_sockfd[n]==-1)
            return n;    
    }
}
