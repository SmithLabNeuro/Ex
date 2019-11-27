/* MATLABUDP.H
 *
 *	Header file for MATLABUDP.c, which contains a few c-routines
 *	to be called from MATLAB so that dotsX machines can chat via
 *	ethernet and the UDP/IP protocols.
 *
 *  This is as close as possible to the code in matlabUDP.h on the
 *  REX machine.  The only real difference is that we have to implement
 *  the mexFunction interface for MATLAB.
 *
 *
 *	BSH 20 Jan 2006
*/

#ifndef MATLABUDP_H_
#define MATLABUDP_H_

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <sys/time.h>

#include "mex.h"

#define MAX_NUM_BYTES 2000
#define MAX_SOCKETS   4 /*Modify initialization of mat_UDP_sockfd when changing MAX_SOCKETS*/

/*Globals for UDP socket*/
static int                      mat_UDP_sockfd[MAX_SOCKETS]={-1,-1,-1,-1},      /*descriptor of UDP socket*/
                				mat_UDP_addr_len        =sizeof(struct sockaddr),
                                mat_UDP_numBytes[MAX_SOCKETS];       /*length of return message*/
                                
static char mat_UDP_messBuff[MAX_SOCKETS][MAX_NUM_BYTES];            /*used by send and receive*/

static struct sockaddr_in       mat_UDP_LOCAL_addr[MAX_SOCKETS],     /*holds LOCAL IP address */
                                mat_UDP_REMOTE_addr[MAX_SOCKETS];	/*holds REMOTE IP address*/


/*functions for exchanging strings with remote machines */
int	mat_UDP_open	(char*, char*, int);            /*initialize UDP socket*/
void	mat_UDP_send	(int,char*, int);                   /*send a string to MATLAB*/
int     mat_UDP_check	(int);                         /*is a return message available?*/
void	mat_UDP_read	(int,char*, int);                   /*read any available message*/
void	mat_UDP_close	(int);                         /*cleanup UDP socket*/
/*void	mat_UDP_send	(char*, int);                   /*send a string to MATLAB*/
/*int     mat_UDP_check	(void);                         /*is a return message available?*/
/*void	mat_UDP_read	(char*, int);                   /*read any available message*/
/*void	mat_UDP_close	(void);
int     nextSocket      (void);                     /* Checks for next empty socket*/

void mexFunction(
    int           nlhs,           /* number of expected outputs */
    mxArray       *plhs[],        /* array of pointers to output arguments */
    int           nrhs,           /* number of inputs */
    const mxArray *prhs[]         /* array of pointers to input arguments */
    );


#endif /* MATLABUDP_H_ */
