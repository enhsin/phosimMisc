/* -------------------------------------------- */
/*                  ImSim                       */
/*            Electron to ADC code              */
/* -------------------------------------------- */
/* Written by John R. Peterson (Purdue)         */
/* for the ImSim group                          */
/* as part of the LSST Project                  */
/* -------------------------------------------- */ 
/* WARNING:   This code is not fully validated  */
/* and not ready for full release.  Please      */
/* treat results with caution.                  */
/* -------------------------------------------- */ 

#include <cstdio>
#include <string>
#include <fitsio.h>
#include <fitsio2.h>



int main(int argc,char **argv) {

 char inputfilename[100], comment[10], val[100], command[100];
 std::string keyname, value;
 fitsfile *foptr;
 int status=0;

 sprintf(inputfilename,argv[1]);
 comment[0]='\0';
 sprintf(command,"gunzip %s.gz",argv[1]);
 system(command);
 if (fits_open_file(&foptr,inputfilename,READWRITE,&status)) {  //only work with uncompressed file in readwrite mode
   printf("Error opening %s\n",inputfilename);
   return 1;
 }
 keyname="CTYPE1";
 value="RA--TAN";
 sprintf(val,value.c_str());
 fits_update_key(foptr,TSTRING,keyname.c_str(),val,comment,&status);
 if (status != 0) {
     printf("FITS Header Key Error: %d\n", status);
     return 1;
 }
 keyname="RADESYS";
 value="ICRS";
 sprintf(val,value.c_str());
 fits_update_key(foptr,TSTRING,keyname.c_str(),val,comment,&status);
 if (status != 0) {
     printf("FITS Header Key Error: %d\n", status);
     return 1;
 }
 fits_close_file(foptr,&status);
 sprintf(command,"gzip %s",argv[1]);
 system(command);
 return 0;

}
