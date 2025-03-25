#include "StRHICfEventDst.h"

ClassImp(StRHICfEventDst)

StRHICfEventDst::StRHICfEventDst()
: mOffTPCTrack(false), mOffFMS(false), mOffRPS(false)
{
}

StRHICfEventDst::~StRHICfEventDst()
{
}

void StRHICfEventDst::Clear(Option_t* option)
{
    if(mEvent){mEvent -> Clear();}
    if(mTPCTrack && !mOffTPCTrack){mTPCTrack -> Clear("C");}
    if(mBTof){mBTof -> Clear("C");}
    if(mBBC){mBBC -> Clear();}
    if(mVPD){mVPD -> Clear();}
    if(mZDC){mZDC -> Clear();}
    if(mRHICfParticle){mRHICfParticle -> Clear("C");}
    if(mRHICfDetPoint){mRHICfDetPoint -> Clear("C");}
    if(mFMS && !mOffFMS){mFMS -> Clear("C");}
    if(mRPS && !mOffRPS){mRPS -> Clear("C");}
}

Int_t StRHICfEventDst::Init()
{
    if(!mEvent){mEvent = new StRHICfEvent();}
    if(!mTPCTrack && !mOffTPCTrack){mTPCTrack = new TClonesArray("StRHICfTPCTrack");}
    if(!mBTof){mBTof = new TClonesArray("StRHICfBTof");}
    if(!mBBC){mBBC = new StRHICfBBC();}
    if(!mVPD){mVPD = new StRHICfVPD();}
    if(!mZDC){mZDC = new StRHICfZDC();}
    if(!mRHICfParticle){mRHICfParticle = new TClonesArray("StRHICfParticle");}
    if(!mRHICfDetPoint){mRHICfDetPoint = new TClonesArray("StRHICfDetPoint");}
    if(!mFMS && !mOffFMS){mFMS = new TClonesArray("StRHICfFMS");}
    if(!mRPS && !mOffRPS){mRPS = new TClonesArray("StRHICfRPS");}

    return 1;
}

Int_t StRHICfEventDst::CreateDstArray(TTree* tree)
{
    if(mEvent){tree -> Branch("StRHICfEvent", &mEvent);}
    if(mTPCTrack){tree -> Branch("StRHICfTPCTrack", &mTPCTrack);}
    if(mBTof){tree -> Branch("StRHICfBTof", &mBTof);}
    if(mBBC){tree -> Branch("StRHICfBBC", &mBBC);}
    if(mVPD){tree -> Branch("StRHICfVPD", &mVPD);}
    if(mZDC){tree -> Branch("StRHICfZDC", &mZDC);}
    if(mRHICfParticle){tree -> Branch("StRHICfParticle", &mRHICfParticle);}
    if(mRHICfDetPoint){tree -> Branch("StRHICfDetPoint", &mRHICfDetPoint);}
    if(mFMS){tree -> Branch("StRHICfFMS", &mFMS);}
    if(mRPS){tree -> Branch("StRHICfRPS", &mRPS);}

    return 1;
}

Int_t StRHICfEventDst::ReadDstArray(TTree* tree)
{
    if(mEvent){tree -> SetBranchAddress("StRHICfEvent", &mEvent);}
    if(mTPCTrack){tree -> SetBranchAddress("StRHICfTPCTrack", &mTPCTrack);}
    if(mBTof){tree -> SetBranchAddress("StRHICfBTof", &mBTof);}
    if(mBBC){tree -> SetBranchAddress("StRHICfBBC", &mBBC);}
    if(mVPD){tree -> SetBranchAddress("StRHICfVPD", &mVPD);}
    if(mZDC){tree -> SetBranchAddress("StRHICfZDC", &mZDC);}
    if(mRHICfParticle){tree -> SetBranchAddress("StRHICfParticle", &mRHICfParticle);}
    if(mRHICfDetPoint){tree -> SetBranchAddress("StRHICfDetPoint", &mRHICfDetPoint);}
    if(mFMS){tree -> SetBranchAddress("StRHICfFMS", &mFMS);}
    if(mRPS){tree -> SetBranchAddress("StRHICfRPS", &mRPS);}

    return 1;
}

Int_t StRHICfEventDst::ReadDstArray(TChain* chain)
{
    if(mEvent){chain -> SetBranchAddress("StRHICfEvent", &mEvent);}
    if(mTPCTrack){chain -> SetBranchAddress("StRHICfTPCTrack", &mTPCTrack);}
    if(mBTof){chain -> SetBranchAddress("StRHICfBTof", &mBTof);}
    if(mBBC){chain -> SetBranchAddress("StRHICfBBC", &mBBC);}
    if(mVPD){chain -> SetBranchAddress("StRHICfVPD", &mVPD);}
    if(mZDC){chain -> SetBranchAddress("StRHICfZDC", &mZDC);}
    if(mRHICfParticle){chain -> SetBranchAddress("StRHICfParticle", &mRHICfParticle);}
    if(mRHICfDetPoint){chain -> SetBranchAddress("StRHICfDetPoint", &mRHICfDetPoint);}
    if(mFMS){chain -> SetBranchAddress("StRHICfFMS", &mFMS);}
    if(mRPS){chain -> SetBranchAddress("StRHICfRPS", &mRPS);}

    return 1;
}

void StRHICfEventDst::OffTPCTrack(){mOffTPCTrack =true;}
void StRHICfEventDst::OffFMS(){mOffFMS =true;}
void StRHICfEventDst::OffRPS(){mOffRPS =true;}

StRHICfEvent* StRHICfEventDst::GetEvent(){return mEvent;}
StRHICfTPCTrack* StRHICfEventDst::GetTPCTrack(int idx){return (StRHICfTPCTrack*)mTPCTrack->ConstructedAt(idx);}
StRHICfBTof* StRHICfEventDst::GetBTof(int idx){return (StRHICfBTof*)mBTof->ConstructedAt(idx);}
StRHICfBBC* StRHICfEventDst::GetBBC(){return mBBC;}
StRHICfVPD* StRHICfEventDst::GetVPD(){return mVPD;}
StRHICfZDC* StRHICfEventDst::GetZDC(){return mZDC;}
StRHICfParticle* StRHICfEventDst::GetRHICfParticle(int idx){return (StRHICfParticle*)mRHICfParticle->ConstructedAt(idx);}
StRHICfDetPoint* StRHICfEventDst::GetRHICfDetPoint(int idx){return (StRHICfDetPoint*)mRHICfDetPoint->ConstructedAt(idx);}
StRHICfFMS* StRHICfEventDst::GetFMS(int idx){return (StRHICfFMS*)mFMS->ConstructedAt(idx);}
StRHICfRPS* StRHICfEventDst::GetRPS(int idx){return (StRHICfRPS*)mRPS->ConstructedAt(idx);}

Int_t StRHICfEventDst::GetTPCTrackNum(){return mTPCTrack->GetEntries();}
Int_t StRHICfEventDst::GetBTofNum(){return mBTof->GetEntries();}
Int_t StRHICfEventDst::GetRHICfParticleNum(){return mRHICfParticle->GetEntries();}
Int_t StRHICfEventDst::GetRHICfDetPointNum(){return mRHICfDetPoint->GetEntries();}
Int_t StRHICfEventDst::GetFMSNum(){return mFMS->GetEntries();}
Int_t StRHICfEventDst::GetRPSNum(){return mRPS->GetEntries();}
