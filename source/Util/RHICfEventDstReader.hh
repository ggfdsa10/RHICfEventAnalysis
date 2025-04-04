#ifndef RHICfEventDstReader_hh
#define RHICfEventDstReader_hh

#include <iostream>
#include <fstream>

#include "TSystem.h"
#include "TString.h"
#include "TChain.h"
#include "TObjArray.h"

#include "RHICfOptContainer.hh"
#include "StRHICfEventDst.h"

class RHICfEventDstReader
{
    public:
        RHICfEventDstReader();
        ~RHICfEventDstReader();

        void Init();

        Int_t GetEventNum();
        StRHICfEventDst* GetEventDst(int idx);

    private:
        RHICfOptContainer* mOptContainer;

        Int_t mEventNum;

        TChain* mChain;
		StRHICfEventDst* mRHICfEventDst;


};

#endif
