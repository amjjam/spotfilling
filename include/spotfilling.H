/******************************************************************************
 * This is class SPOTFILLING. It is an alternative to the default filling     *
 * function in DGCPM. It calls the default filling function and then, during  *
 * a time-period, it fills a spot much more rapidly. The rapid filling is     *
 * supposed to simulate the increased ionization from a substorm.             *
 ******************************************************************************/

#include "../submodules/include/dgcpm.H"
#include "../submodules/include/aTime.H"

class SPOTFILLING: public FILLING{
public:
  SPOTFILLING(float fmax=2e12, float tauclosed=86400, float tauopen=86400);
  ~SPOTFILLING();
  void setSpot(aTime tstart, aTime tend, float t, float p, float r, float f);
  void setTime(aTime time);
  void filling(std::vector<float> &vR, std::vector<float> &vT, 
	       std::vector<float> &vP, GRID &mGridN, GRID &mGridDen,
	       GRID &mGridVol, GRID &mGridOc, GRID &mGridBi, float dt);
private:
  aTime t,tStart,tEnd;
  float tCenter,pCenter,R,f;
};
