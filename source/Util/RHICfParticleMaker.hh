#ifndef RHICfParticleMaker_hh
#define RHICfParticleMaker_hh

#include <vector>

#include "RHICfOptContainer.hh"
#include "StRHICfEventDst.h"
#include "StRHICfEvent.h"
#include "StRHICfDetHit.h"
#include "StRHICfDetPoint.h"

#include "StRHICfSimDst.h"
#include "StRHICfSimEvent.h"
#include "StRHICfSimRHICfHit.h"
#include "StRHICfSimRHICfPoint.h"

struct MiniParticle
{    
    int particleNum;
    std::vector<int> label;
    std::vector<int> towerIdx;
    std::vector<int> type;
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
    std::vector<double> pt;
    std::vector<double> xf;
    std::vector<double> px;
    std::vector<double> py;
    std::vector<double> pz;
    std::vector<double> e;
    std::vector<double> m;
};

class RHICfParticleMaker
{
    public:
        RHICfParticleMaker();
        ~RHICfParticleMaker();

        void Init();

        void SetEventDst(StRHICfEventDst* eventDst);
        void SetSimDst(StRHICfSimDst* simDst);

        MiniParticle GetMiniParticle();

        int GetGlobalPosition(double& localX, double& localY, int towerIdx, int fillNum, int method, int ref);
        double BeamCenterPosition(int fillNum, int xy, int method, int ref);
        double BeamCenterPositionErr(int fillNum, int xy, int method, int ref);

    private:
        MiniParticle GetGamma();
        MiniParticle GetPi0();
        MiniParticle GetNeutron();
        MiniParticle GetLambda0();


        RHICfOptContainer* mOptContainer;
        bool mIsSimData;

        StRHICfEventDst* mEventDst;
        StRHICfEvent* mEvent;
        StRHICfDetHit* mDetHit;
        StRHICfDetPoint* mDetPoint;

        StRHICfSimDst* mSimDst;
        StRHICfSimEvent* mSimEvent;
        StRHICfSimRHICfHit* mSimRHICfHit;
        StRHICfSimRHICfPoint* mSimRHICfPoint;

        int mBeamCenterMethod;
        int mBeamCenterRef;

        // const double distZatIP = 18000.0; // [mm] Distance of the beam from interaction point to RHICf detector.
        const double distZatIP = 17812.; // [mm] Distance of the beam from interaction point to RHICf detector.
		const double geoCenterTL = 20.0; // [mm] geometrical large tower center
		const double geoCenterTS = 10.0; // [mm] geometrical small tower center
		const double distTStoTL = 47.4; // [mm] Distance between center of TL and cernter of TS


};

#endif
