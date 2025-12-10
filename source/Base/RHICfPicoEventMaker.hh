#ifndef RHICfPicoEventMaker_hh
#define RHICfPicoEventMaker_hh

#include <iostream>
#include <fstream>
#include <map>

#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TStyle.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TObjArray.h"

#include "RHICfOptContainer.hh"
#include "RHICfEventDstReader.hh"
#include "RHICfTableMaker.hh"
#include "RHICfParticleMaker.hh"
#include "RHICfDLECondition.hh"

#include "StRHICfEventDst.h"
#include "StRHICfEvent.h"
#include "StRHICfDetHit.h"
#include "StRHICfDetPoint.h"
#include "StRHICfBBC.h"
#include "StRHICfRPS.h"
#include "StRHICfZDC.h"

class RHICfPicoEventMaker : public RHICfOptContainer, RHICfTableMaker
{
    public:
        RHICfPicoEventMaker();
        ~RHICfPicoEventMaker();

        int Init();
        int Calculate();
        int Finish();

        void MakeRHICfStream(){mIsRHICfStream = true;}
        void MakePhysicsStream(){mIsPhysicsStream = true;}
        void MakeRPSEvent(){mIsRPSEvent = true;}
        void MakeRHICfParticleEvent(){mIsRHICfParticleEvent = true;}

    private:
        void RHICfStream();
        void PhysicsStream();
        void RPSEvent();
        void RHICfParticleEvent();

        // Utilized classes
        RHICfEventDstReader* mEventReader;
        RHICfParticleMaker* mParticleMaker;
        RHICfDLECondition* mDLECondition;

        // Options
        TString mOutputPath;
        bool mIsRHICfStream;
        bool mIsPhysicsStream;
        bool mIsRPSEvent;
        bool mIsRHICfParticleEvent;

        // StRHICfMiniSimDst 
        StRHICfEventDst* mEventDst;
        StRHICfEvent* mEvent;
        StRHICfBBC* mBBC;
        StRHICfDetHit* mDetHit;
        StRHICfDetPoint* mDetPoint;
        StRHICfRPS* mRPS;
        StRHICfZDC* mZDC;

};

#endif
