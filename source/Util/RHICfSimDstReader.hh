#ifndef RHICfSimDstReader_hh
#define RHICfSimDstReader_hh

#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TSystem.h"
#include "TChain.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TString.h"
#include "TVector3.h"

#include "RHICfOptContainer.hh"
#include "RHICfParticleMaker.hh"
#include "RHICfDLECondition.hh"

#include "StRHICfSimDst.h"
#include "StRHICfSimEvent.h"
#include "StRHICfSimTrack.h"
#include "StRHICfSimBTof.h"
#include "StRHICfSimBBC.h"
#include "StRHICfSimZDC.h"
#include "StRHICfSimRHICfHit.h"
#include "StRHICfSimRHICfPoint.h"

#include "StRHICfMiniSimDst.h"

class RHICfSimDstReader
{
    public:
        RHICfSimDstReader();
        ~RHICfSimDstReader();

        void Init();
        void Make();

        Int_t GetEventNum();
        StRHICfMiniSimDst* GetMiniSimDst(int idx);

    private:
        bool FindMiniSimDst();
        void MakeMiniSimDst();

        void SaveTruthSimTracks();

        RHICfOptContainer* mOptContainer;
        RHICfParticleMaker* mParticleMaker;
        RHICfDLECondition* mDLECondition;

        Int_t mEventNum;

        TChain* mChain;
        StRHICfSimDst* mSimDst;
        StRHICfSimEvent* mSimEvent;
        StRHICfSimTrack* mSimTrack;
        StRHICfSimBTof* mSimBTof;
        StRHICfSimBBC* mSimBBC;
        StRHICfSimZDC* mSimZDC;
        StRHICfSimRHICfHit* mSimRHICfHit;
        StRHICfSimRHICfPoint* mSimRHICfPoint;

        TFile* mMiniDstFile;
        TTree* mMiniDstTree;
        StRHICfMiniSimDst* mMiniSimDst;

        TString mModelName;
        TString mInputDataPath;
        TString mDataPath;
        TString mDataList;
        TString mMiniDstName;
};

#endif
