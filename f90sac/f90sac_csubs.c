/*==============================================================================
 -------------------------------------------------------------------------------
 
   C Source Code File
 
 -------------------------------------------------------------------------------
 ===============================================================================
 
   PROGRAM : f90sac
   VERSION : 4.0
   CVS: $Revision: 1.5 $ $Date: 2008/10/17 10:48:03 $
 
   (C) James Wookey
   Department of Earth Sciences, University of Bristol
   Wills Memorial Building, Queen's Road, Bristol, BR8 1RJ, UK
   j.wookey@bristol.ac.uk
 
 -------------------------------------------------------------------------------
 
    The module provides data structures and functions for reading,
    writing and handling SAC files in Fortran 90/95.
 
    Please report bugs/problems to email address above
  
    NOTE: This version of the code assumes IO filestream 99 is available 
          for reading and writing. 
 
 -------------------------------------------------------------------------------
 
   This software is distributed under the term of the BSD free software license.
 
   Copyright:
      (c) 2003-2008, James Wookey
 
   All rights reserved.
 
    * Redistribution and use in source and binary forms, with or without
      modification, are permitted provided that the following conditions are
      met:
         
    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
         
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
      
    * Neither the name of the copyright holder nor the names of its
      contributors may be used to endorse or promote products derived from
      this software without specific prior written permission.
 
 
    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
    "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
    LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
    A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT
    OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 
 -----------------------------------------------------------------------------*/

//
// These routines are the low level routines used to read and write SAC files.
// They are not designed to be used directly, but by the routines in f90sac.90 
// Use them in other ways at your peril!
//

#include <stdio.h>
#include <stdlib.h>

// Global file descriptor, used for reading and writing traces
   FILE *pfd ;

//=============================================================================
void f90sac_c_openr_(char pfname[], int pfnlen) 
//-----------------------------------------------------------------------------
//
//    Open a file for reading in C   
//
{  
	char fname[256] ;
	int i;

	for(i=0;i<pfnlen;++i) fname[i] = pfname[i] ;
	for(--i;fname[i]==' ' && i>0;--i) ;  /* get rid of trailing blanks */
	fname[i+1]='\0' ; 
	
   pfd = fopen(fname,"r") ;
	if (pfd==NULL) {
		printf("F90SAC: Error: Specified file ( %s ) could not be opened for reading\n",fname) ;
		exit(1);
	}

	return ;
}
//=============================================================================

//=============================================================================
void f90sac_c_openw_(char pfname[], int pfnlen) 
//-----------------------------------------------------------------------------
//
//    Open a file for writing in C   
//
{  
	char fname[256] ;
	int i;

	for(i=0;i<pfnlen;++i) fname[i] = pfname[i] ;
	for(--i;fname[i]==' ' && i>0;--i) ;  /* get rid of trailing blanks */
	fname[i+1]='\0' ; 
	
   pfd = fopen(fname,"w") ;
	if (pfd==NULL) {
		printf("F90SAC: Error: Specified file ( %s ) could not be opened for writing\n",fname) ;
		exit(1);
	}

	return ;
}
//=============================================================================

//=============================================================================
void f90sac_c_openrw_(char pfname[], int pfnlen) 
//-----------------------------------------------------------------------------
//
//    Open a file for reading/overwriting in C, used for writing header   
//
{  
	char fname[256] ;
	int i;

	for(i=0;i<pfnlen;++i) fname[i] = pfname[i] ;
	for(--i;fname[i]==' ' && i>0;--i) ;  /* get rid of trailing blanks */
	fname[i+1]='\0' ; 
	
   pfd = fopen(fname,"r+") ;
	if (pfd==NULL) {
		printf("F90SAC: Error: Specified file ( %s ) could not be opened for writing\n",fname) ;
		exit(1);
	}

	return ;
}
//=============================================================================

//=============================================================================
void f90sac_c_whd_(float rhd[], int ihd[], char chd[], int chdlen) 
//-----------------------------------------------------------------------------
//
// Read the header.    
//
{  
	fwrite(rhd,4,70,pfd) ;
	fwrite(ihd,4,40,pfd) ;
	fwrite(chd,1,192,pfd) ;

	return ;
}
//=============================================================================

//=============================================================================
void f90sac_c_rhd_(float rhd[], int ihd[], char chd[], int chdlen) 
//-----------------------------------------------------------------------------
//
// Write the header.    
//
{  
	fread(rhd,4,70,pfd) ;
	fread(ihd,4,40,pfd) ;
	fread(chd,1,192,pfd) ;

	return ;
}
//=============================================================================

//=============================================================================
void f90sac_c_wtr_(int *npts, float tr[]) 
//-----------------------------------------------------------------------------
//
// Read the trace.    
//
{  
	fwrite(tr,4,*npts,pfd) ;
	
	return ;
}
//=============================================================================

//=============================================================================
void f90sac_c_rtr_(int *npts, float tr[]) 
//-----------------------------------------------------------------------------
//
// Read the trace.    
//
{  
	fread(tr,4,*npts,pfd) ;
	
	return ;
}
//=============================================================================

//=============================================================================
void f90sac_c_close_() 
//-----------------------------------------------------------------------------
//
// Close the file.    
//
{  
   int status ;
	status = fclose(pfd) ;
	
	return ;
}
//=============================================================================
