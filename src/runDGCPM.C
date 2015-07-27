/******************************************************************************
 * This program is used to run DGCPM with a input Kp file and other           *
 * parameters.                                                                *
 ******************************************************************************/

/*=============================================================================
  runDGCPM [-s yr mo dy hr ] [-e yr mo dy hr] [-so yr mo dy hr] 
  [-dt float] [-T float] [-o <file> ] [-samples <file> ] 
  [-filling|-f <fMax> <tauClosed> <tauOpen>] [-saturation <A> <B>]
  <ifile1> [<ifile2>.. ]

  Runs the DGCPM model and writes the output to a file.

  -s yr mo dy hr - the start time of the run in year month day hour
     UT. If not specified the default start time is the first time in
     the Kp input files
  -e yr mo dy hr - the end time of the run in year month day hour UT. If 
     not specified then use -T for runtime
  -so yr mo dy hr - the start time for writing output. This is useful for 
     pre-conditioning the run for a few days before generating output.
  -dt <float> - time, in seconds, between writing images to file. If
     not specified the default is 15 minutes (900 s).
  -T Duration of the run in seconds. Ignored if -e is specified. If
     neither -e or -T are specified then run goes to the last time 
     specified in the Kp input files.
  -o <ofile> - the file to write the output to. If not specified the
     default is output.dat
  -filling|-f <fMax> <tauClosed> <tauOpen> Set the parameters fMax, tauClosed, 
     and tauOpen in the filling function. Units are particles/m^2/s for fMax,
     and days for tauClosed and tauOpen. Default values if -f is not
     used are 2e12, 10, and 1. 
  -saturation <A> <B> Use saturation function neq=10^(A+B*l) with those values
     for parametes A and B. Default values are 3.9043 and -0.3145
  -samples <file> If specified then don't write out images but instead write
     out samples at the locations specified in this file. The file contains 
     pairs of L-shell and magnetic longitude (in degrees). The output will
     be for the L,MLT locations corresponding to those locations.  
  The following are inputs for the spot
  -sStart int - The time to turn on the spot, in seconds after start time.
  -sStop int - The time to turn off the spot, in seconds after stop time.
  -sT float - the co-latitude of the center of the spot, in degrees
  -sP float - the local time of the center of the spot, in degrees 
     east from midnight. 
  -sR float - the radius of the spot in kilometers at the surface of the Earth.
  -sF float - the amplification factor of the spot. fMax and dSat in filling
     formula are increased by this factor in the spot.
  <ifiles> - input Kp files in WDC format. Can be specified multiple
     times and the files are added in the order they appear on the
     command line. Make sure they are specified in increasing time
     order.
  ============================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <zlib.h>
#include <string.h>
#include <libgen.h>
#include <time.h>

#include "../submodules/include/dgcpm.H"
#include "../submodules/include/sample.H"
#include "../submodules/include/kp.H"

#include "../include/spotfilling.H"

void parseArgs(int argc, char *argv[]);
void printTime(aTime &t);
void writeState(aTime &t, gzFile fp, DGCPM &m);
aTime &writeSamples(aTime &t, DGCPM &m);

std::vector<std::string> iFiles;
std::string oFile;
std::string samplesIFile;
double dt=900;
double T=-1;
aTime tStart,tStop,tOut;

int filling=0;
float fMax,tauClosed,tauOpen;

int saturation=0;
float saturationA,saturationB;

int ePotModel=EPOT_SOJKA;

SAMPLE *samples=NULL;

// Parameters related to the spot
double sStartDt=1e31,sStopDt=-1e31;
double sT=30,sP=315,sR=1000,sF=10;

int main(int argc, char *argv[]){
  tStart.set(0);
  tStop.set(0);
  tOut.set(0);

  parseArgs(argc,argv);
  
  // Load the Kp data
  KPS kp(iFiles);
  
  // Determine start and end times
  if(tStart.get()<1)
    tStart=kp[0].getTime();
  if(tStop<tStart){
    if(T<0)
      tStop=kp[kp.size()-1].getTime();
    else{
      tStop=tStart;
      tStop+=T;
    }
  }
  if(tOut.get()<1)
    tOut=tStart;
  
  // Set initial pointer in kp
  int iKp=kp.find(tStart);
  aTime tKp=kp[iKp].getTime();

   // Create the parameters array
  float par[1]={kp[iKp].getKp()};
  
  // Create the DGCPM model
  DGCPM m;
  m.setEPot(ePotModel,par);

  // If a different filling function was specified then create it here
  // and attach it.
  SPOTFILLING *f;
  if(filling==1){
    f=new SPOTFILLING(fMax,tauClosed,tauOpen);
    m.setFilling(f);
    aTime sStart,sStop;
    sStart=tStart;
    sStart+=sStartDt;
    sStop=tStart;
    sStop+=sStopDt;
    f->setSpot(sStart,sStop,sT,sP,sR,sF);
  }

  // If a different saturation function was specified then create it
  // here and attach it to the filling function
  SATURATION *s;
  if(saturation==1){
    s=new SATURATION(saturationA,saturationB);
    f->setSaturation(s);
  }

  // If doing samples create the samples object
  aTime tWriteSample=tStop;
  tWriteSample+=1;
  if(samplesIFile.size()>0){
    samples=new SAMPLE(samplesIFile,tOut,dt,oFile);
    tWriteSample=samples->getTime();
  }

  // If not doing samples then do density images
  aTime tWriteState=tStop;
  tWriteState+=1;
  gzFile oFp;
  if(samples==NULL){
    if(oFile.size()==0)
      oFile="output.dat";
    
    oFp=gzopen(oFile.c_str(),"w9");;
    m.writeHeader(oFp);
    tWriteState=tOut;
  }


  // Loop over time
  aTime t=tStart;
  aTime tNext=tStart;
  aTime tFilling=tStart;
  for(;tNext<=tStop;){
    // Set the time for the filling function
    f->setTime(t);
    tFilling+=300;

    if(tNext-t>0){
      std::cout << tNext-t << std::endl;
      m.advance(tNext-t);
      t=tNext;
    }
    
    printTime(t);

    if(t>=tKp){
      std::cout << "Kp " << kp[iKp].getKp() << std::endl;
      par[0]=kp[iKp].getKp();
      m.setEPot(ePotModel,par);
      iKp++;
      if(iKp>=kp.size()){
	tKp=tStop;
	tKp+=1;
      }
      else
	tKp=kp[iKp].getTime();
    }
    
    if(t>=tWriteState){
      std::cout << "Writing state" << std::endl;
      writeState(t,oFp,m);
      tWriteState+=dt;
    }
    
    if(t>=tWriteSample){
      std::cout << "Writing sample" << std::endl;
      tWriteSample=writeSamples(t,m);
    }
    
    tNext=tWriteSample;
    if(tWriteState<tNext)
      tNext=tWriteState;
    if(tKp<tNext)
      tNext=tKp;
    if(tFilling<tNext)
      tNext=tFilling;
  }
  
  if(filling==1)
    delete f;

  if(samples!=NULL)
    delete samples;

  if(samples==NULL)
    gzclose(oFp);

  return 0;
}


/*============================================================================
  parseArgs - parse command line arguments.
  ============================================================================*/
void parseArgs(int argc, char *argv[]){
  int i;

  for(i=1;i<argc;i++)
    if(strcmp(argv[i],"-h")==0||strcmp(argv[i],"-help")==0||
       strcmp(argv[i],"--help")==0){
      std::cout << "runDGCPM [-s yr mo dy hr ] [-e yr mo dy hr] "
		<< "[-so yr mo dy hr]" << std::endl; 
      std::cout << "[-dt float] [-T float] [-o <file> ] [-samples <file> ]" 
		<< std::endl;
      std::cout << "[-filling|-f <fMax> <tauClosed> <tauOpen>] "
		<< "[-saturation <A> <B>]" << std::endl;
      std::cout << "<ifile1> [<ifile2>.. ]" << std::endl;
      std::cout << "" << std::endl;
      std::cout << "Runs the DGCPM model and writes the output to a file." 
		<< std::endl;
      std::cout << "" << std::endl;
      std::cout << "-s yr mo dy hr - the start time of the run in year " 
		<< "month day hour" << std::endl;
      std::cout << "   UT. If not specified the default start time is the "
		<< "first time in" << std::endl;
      std::cout << "   the Kp input files" << std::endl;
      std::cout << "-e yr mo dy hr - the end time of the run in year "
		<<"month day hour UT. If" << std::endl;
      std::cout << "   not specified then use -T for runtime" << std::endl;
      std::cout << "-so yr mo dy hr - the start time for writing output. "
		<< "This is useful for " << std::endl;
      std::cout << "   pre-conditioning the run for a few days before "
		<< "generating output." << std::endl;
      std::cout << "-dt <float> - time, in seconds, between writing images " 
		<< "to file. If" << std::endl;
      std::cout << "   not specified the default is 15 minutes (900 s)." 
		<< std::endl;
      std::cout << "-T Duration of the run in seconds. Ignored if -e "
		<<"is specified. If" << std::endl;
      std::cout << "   neither -e or -T are specified then run goes to "
		<<"the last time" << std::endl;
      std::cout << "   specified in the Kp input files." << std::endl;
      std::cout << "-o <ofile> - the file to write the output to. If not "
		<< "specified the" << std::endl;
      std::cout << "   default is output.dat" << std::endl;
      std::cout << "  -filling|-f <fMax> <tauClosed> <tauOpen> Set the "
		<< "parameters fMax, tauClosed," << std::endl;
      std::cout << "     and tauOpen in the filling function. Units are "
		<< "particles/m^2/s for fMax," << std::endl;
      std::cout << "     and days for tauClosed and tauOpen. Default "
		<< "values if -f is not" << std::endl;
      std::cout << "     used are 2e12, 10, and 1." << std::endl;
      std::cout << "  -saturation <A> <B> Use saturation function "
		<< "neq=10^(A+B*l) with those values" << std::endl;
      std::cout << "     for parametes A and B. Default values are "
		<< "3.9043 and -0.3145" << std::endl;
      std::cout << "-samples <file> If specified then don't write out "
		<< "images but instead write" << std::endl;
      std::cout << "   out samples at the locations specified in this file. "
		<< "The file contains" << std::endl; 
      std::cout << "   pairs of L-shell and magnetic longitude (in degrees). "
		<< "The output will" << std::endl;
      std::cout << "   be for the L,MLT locations corresponding to those "
		<< "locations." << std::endl;
      std::cout << "<ifiles> - input Kp files in WDC format. Can be "
		<< "specified multiple" << std::endl;
      std::cout << "   times and the files are added in the order they "
		<< "appear on the" << std::endl;
      std::cout << "   command line. Make sure they are specified in "
		<< "increasing time" << std::endl;
      std::cout << "   order." << std::endl;
      std::cout << "The following are inputs for the spot" << std::endl;
      std::cout << "-sStart int - The time to turn on the spot, in seconds "
		<< "after start time." << std::endl;
      std::cout << "-sStop int - The time to turn off the spot, in seconds "
		<< "after stop time." << std::endl;
      std::cout << "-sT float - the co-latitude of the center of the spot, "
		<< "in degrees" << std::endl;
      std::cout << "-sP float - the local time of the center of the spot, "
		<< "in degrees" << std::endl; 
      std::cout << "   east from midnight. " << std::endl;
      std::cout << "-sR float - the radius of the spot in kilometers at the "
		<< "surface of the Earth." << std::endl;
      std::cout << "-sF float - the amplification factor of the spot. fMax "
		<< "and dSat in filling" << std::endl;
      std::cout << "   formula are increased by this factor in the spot." 
		<< std::endl;
      exit(0);
    }
  
  for(i=1;i<argc;i++){
    if(strcmp(argv[i],"-s")==0){
      int yr,mo,dy,hr;
      i++;
      yr=atoi(argv[i]);
      i++;
      mo=atoi(argv[i]);
      i++;
      dy=atoi(argv[i]);
      i++;
      hr=atoi(argv[i]);
      tStart.set(yr,mo,dy,hr);
    }
    else if(strcmp(argv[i],"-e")==0){
      int yr,mo,dy,hr;
      i++;
      yr=atoi(argv[i]);
      i++;
      mo=atoi(argv[i]);
      i++;
      dy=atoi(argv[i]);
      i++;
      hr=atoi(argv[i]);
      tStop.set(yr,mo,dy,hr);
    }
    else if(strcmp(argv[i],"-so")==0){
      int yr,mo,dy,hr;
      i++;
      yr=atoi(argv[i]);
      i++;
      mo=atoi(argv[i]);
      i++;
      dy=atoi(argv[i]);
      i++;
      hr=atoi(argv[i]);
      tOut.set(yr,mo,dy,hr);
    }
    else if(strcmp(argv[i],"-dt")==0){
      i++;
      dt=atof(argv[i]);
    }
    else if(strcmp(argv[i],"-T")==0){
      i++;
      T=atof(argv[i]);
    }
    else if(strcmp(argv[i],"-o")==0){
      i++;
      oFile=argv[i];
    }
    else if(strcmp(argv[i],"-filling")==0||strcmp(argv[i],"-f")==0){
      filling=1;
      i++;
      fMax=atof(argv[i]);
      i++;
      tauClosed=atof(argv[i])*86400;
      i++;
      tauOpen=atof(argv[i])*86400;
    }
    else if(strcmp(argv[i],"-saturation")==0){
      saturation=1;
      i++;
      saturationA=atof(argv[i]);
      i++;
      saturationB=atof(argv[i]);
    }
    else if(strcmp(argv[i],"-samples")==0){
      i++;
      samplesIFile=std::string(argv[i]);
    }
    else if(strcmp(argv[i],"-sStart")==0){
      i++;
      sStartDt=atof(argv[i]);
    }
    else if(strcmp(argv[i],"-sStop")==0){
      i++;
      sStopDt=atof(argv[i]);
    }
    else if(strcmp(argv[i],"-sT")==0){
      i++;
      sT=atof(argv[i]);
    }
    else if(strcmp(argv[i],"-sP")==0){
      i++;
      sP=atof(argv[i]);
    }
    else if(strcmp(argv[i],"-sR")==0){
      i++;
      sR=atof(argv[i]);
    }
    else if(strcmp(argv[i],"-sF")==0){
      i++;
      sF=atof(argv[i]);
    }
    else if(argv[i][0]=='-'){
      std::cout << "Error: unknown option: " << argv[i] << std::endl;
      exit(1);
    }
    else
      iFiles.push_back(argv[i]);
  }

  if(iFiles.size()==0){
    std::cout << "No input Kp files specified." << std::endl;
    exit(1);
  }

  if(saturation==1&&filling==0){
    std::cout << "Must use custom filling model in order to use custom "
	      << "saturation model." << std::endl;
    exit(1);
  }
}


/*=============================================================================
  void printTime(aTime &t) - print the time
  ============================================================================*/
void printTime(aTime &t){
  int yr,mo,dy,hr,mn,se;
  t.get(yr,mo,dy,hr,mn,se);
  std::cout << yr << "/" << mo << "/" << dy << " " << hr << ":" << mn << ":"
	    << se << std::endl;
}


/*=============================================================================
  void writeState(aTime &t, DGCPM &m) - write the state of the model to file
  
  aTime &t - the current time to associate with the sate written
  DGCPM &m - the model whose state is to be written
  ============================================================================*/
void writeState(aTime &t, gzFile fp, DGCPM &m){
  int yr,mo,dy,hr,mn,se;
  t.get(yr,mo,dy,hr,mn,se);

  gzwrite(fp,&yr,sizeof(int));
  gzwrite(fp,&mo,sizeof(int));
  gzwrite(fp,&dy,sizeof(int));
  gzwrite(fp,&hr,sizeof(int));
  gzwrite(fp,&mn,sizeof(int));
  gzwrite(fp,&se,sizeof(int));
  
  m.writeState(fp);
}


/*=============================================================================
  aTime &writeSamples(aTime &t, DGCPM &m) - 
  ============================================================================*/
aTime &writeSamples(aTime &t, DGCPM &m){
  (*samples)(m);
  ++(*samples);
  return samples->getTime();
}
