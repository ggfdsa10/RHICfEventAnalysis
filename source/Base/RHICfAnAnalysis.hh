#ifndef RHICfAnAnalysis_hh
#define RHICfAnAnalysis_hh

#include "RHICfOptContainer.hh"
#include "RHICfEventDstReader.hh"
#include "RHICfTableMaker.hh"
#include "RHICfParticleMaker.hh"

#include "RHICfFigureDrawing.hh"
#include "RHICfDLECondition.hh"
#include "RHICfMassFitting.hh"
#include "RHICfBinning.hh"
#include "RHICfDilutionFactor.hh"
#include "RHICfPolarization.hh"

#include "StRHICfEventDst.h"
#include "StRHICfEvent.h"
#include "StRHICfTPCTrack.h"
#include "StRHICfBTof.h"
#include "StRHICfBBC.h"
#include "StRHICfRPS.h"

class RHICfAnAnalysis : public RHICfOptContainer, RHICfTableMaker
{
    public:
        RHICfAnAnalysis();
        ~RHICfAnAnalysis();

        int Init();
        int Calculate();
        int Finish();

    private:
        int InitData();
        int MassCalculation();
        int BinningCalculation();
        int DilutionAndPolCalculation();

        // Utilized classes
        RHICfEventDstReader* mEventReader;
        RHICfParticleMaker* mParticleMaker;

        // Calculation classes
        RHICfFigureDrawing* mFigureDrawing;
        RHICfDLECondition* mDLECondition;
        RHICfMassFitting* mMassFitting;
        RHICfBinning* mBinning;
        RHICfDilutionFactor* mDilution;
        RHICfPolarization* mPolarization;

        // StRHICfEventDst 
        StRHICfEventDst* mEventDst;
        StRHICfEvent* mEvent;
        StRHICfBTof* mBTof;
        StRHICfBBC* mBBC;
        StRHICfRPS* mRPS;


};

#endif
