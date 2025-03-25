#ifndef StRHICfEventDst_HH
#define StRHICfEventDst_HH

#include "TChain.h"
#include "TTree.h"
#include "TObject.h"
#include "TClonesArray.h"

#include "StRHICfParameters.h"
#include "StRHICfEvent.h"
#include "StRHICfTPCTrack.h"
#include "StRHICfBTof.h"
#include "StRHICfBBC.h"
#include "StRHICfVPD.h"
#include "StRHICfZDC.h"
#include "StRHICfParticle.h"
#include "StRHICfDetPoint.h"
#include "StRHICfFMS.h"
#include "StRHICfRPS.h"

using namespace std;

class StRHICfEventDst : public TObject
{
    public:
        StRHICfEventDst();
        virtual ~StRHICfEventDst();

        void Clear(Option_t* option = "");

        Int_t Init();
        Int_t CreateDstArray(TTree* tree);
        Int_t ReadDstArray(TTree* tree);
        Int_t ReadDstArray(TChain* chain);

        void OffTPCTrack();
        void OffFMS();
        void OffRPS();

        StRHICfEvent* GetEvent();
        StRHICfTPCTrack* GetTPCTrack(int idx);
        StRHICfBTof* GetBTof(int idx);
        StRHICfBBC* GetBBC();
        StRHICfVPD* GetVPD();
        StRHICfZDC* GetZDC();
        StRHICfParticle* GetRHICfParticle(int idx);
        StRHICfDetPoint* GetRHICfDetPoint(int idx);
        StRHICfFMS* GetFMS(int idx);
        StRHICfRPS* GetRPS(int idx);

        Int_t GetTPCTrackNum();
        Int_t GetBTofNum();
        Int_t GetRHICfParticleNum();
        Int_t GetRHICfDetPointNum();
        Int_t GetFMSNum();
        Int_t GetRPSNum();
    
    private:
        StRHICfEvent* mEvent = nullptr;
        TClonesArray* mTPCTrack = nullptr;
        TClonesArray* mBTof = nullptr;
        StRHICfBBC* mBBC = nullptr;
        StRHICfVPD* mVPD = nullptr;
        StRHICfZDC* mZDC = nullptr;
        TClonesArray* mRHICfParticle = nullptr;
        TClonesArray* mRHICfDetPoint = nullptr;
        TClonesArray* mFMS = nullptr;
        TClonesArray* mRPS = nullptr;

        bool mOffTPCTrack;
        bool mOffFMS;
        bool mOffRPS;

    ClassDef(StRHICfEventDst, 1)
};

#endif