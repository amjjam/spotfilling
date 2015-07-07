#include "../include/spotfilling.H"

/*=============================================================================
  SPOTFILLING(float fmax=2e12, float tauclosed=86400, float
  tauopen=86400) - constructor
  ============================================================================*/
SPOTFILLING::SPOTFILLING(float fmax, float tauclosed, float tauopen):
  FILLING(fmax,tauclosed,tauopen){  
}


/*=============================================================================
  ~SPOTFILLING() - destructor
  ============================================================================*/
SPOTFILLING::~SPOTFILLING(){

}


/*=============================================================================
  void setSpot(double tStart, double tEnd, float t, float p, float r,
  float f) - sets the spot information

  double tStart - the time at which the spot is activated
  double tEnd - the time at which the spot is deactivated
  float t - theta of the center of the spot
  float p - phi of the center of the spot (local time in degrees)
  float r - radius of spot in km at the surface of the Earth
  float f - factor. saturation and fmax are multiplied by this in the spot. 
  ============================================================================*/
void SPOTFILLING::setSpot(double tStart, double tEnd, float t, float p, 
			  float r, float f){
  SPOTFILLING::tStart=tStart;
  SPOTFILLING::tEdn=tEnd;
  SPOTFILLING::tCenter=tCenter;
  SPOTFILLING::pCenter=pCenter;
  SPOTFILLING::R=r;
  SPOTFILLING::f=f;
}


/*=============================================================================
  void setTime(double time) - set the time (used to determine whether
  the spot is on or off).
  ============================================================================*/
void SPOTFILLING::setTime(double time){
  t=time;
}


/*=============================================================================
  void filling(std::vector<float> &vR, std::vector<float> &vT,
  std::vector<float> &vP, GRID &mGridN, GRID &mGridDen, GRID
  &mGridVol, GRID &mGridOc, GRID &mGridBi, float dt) - alternate
  filling function. Calls the default filling function and then does
  special filling in the spot.
  ============================================================================*/
void SPOTFILLING::filling(std::vector<float> &vR, std::vector<float> &vT, 
			  std::vector<float> &vP, GRID &mGridN, GRID &mGridDen,
			  GRID &mGridVol, GRID &mGridOc, GRID &mGridBi, 
			  float dt){
  FILLING::filling(vR,vT,vP,mGridN,mGridDen,mGridVol,mGridOc,mGridBi,dt);
  
  if(tStart<=t&&t<=tEnd){
    // Convert latitude into radius
    float dSat=(*saturation)(1/(sin(tCenter/180*M_PI)^2));
    float sSat=f*dSat;
    float sFMax=f*fMax;
    int iT,nT=vT.size();
    int iP,nP=vP.size();
    float r;
    float dT,dP;
    for(iT=0;iT<nT;iT++)
      for(iP=0;iP<nP;iP++){
	// Compute radial distance from center
	dT=(vT[iT]-tCenter)/180*M_PI*RE;
	dP=vP[iP]-pCenter;
	if(dP>180)
	  dP-=360;
	if(dp<-180)
	  dp+=180;
	dP=dP/180*M_PI*RE*sin(vT[iT]/180*M_PI);
	r=sqrt(dT*dT+dP*dP);
	// Compute a Gaussian based on the radial distance. 
	if(r<R){
	  flux=(dSat-mGridDen[j][i])/dSat*ff;
	  mGridN[j][i]+=flux*dt/mGridBi[j][i];
	  mGridDen[j][i]=mGridN[j][i]/mGridVol[j][i];
	}
      }
  }
}
