#ifndef RHICfParticleMaker_hh
#define RHICfParticleMaker_hh

#include <vector>

#include "RHICfOptContainer.hh"
#include "StRHICfEventDst.h"
#include "StRHICfEvent.h"
#include "StRHICfDetPoint.h"
#include "StRHICfParticle.h"

struct MiniParticle
{    
    int particleNum;
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
        MiniParticle GetMiniParticle();

    private:
        MiniParticle GetGamma();
        MiniParticle GetPi0();
        MiniParticle GetNeutron();
        MiniParticle GetLambda0();

        RHICfOptContainer* mOptContainer;

        StRHICfEventDst* mEventDst;
        StRHICfEvent* mEvent;
        StRHICfDetPoint* mDetPoint;
        StRHICfParticle* mParticle;

};

#endif
