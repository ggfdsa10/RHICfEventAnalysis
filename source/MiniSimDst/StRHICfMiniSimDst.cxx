#include "StRHICfMiniSimDst.h"

ClassImp(StRHICfMiniSimDst)

StRHICfMiniSimDst::StRHICfMiniSimDst()
{
    Clear();
}

StRHICfMiniSimDst::~StRHICfMiniSimDst()
{
}

void StRHICfMiniSimDst::Clear(Option_t *option)
{
    mProcessID = -1;
    mEventType = -1;
    mRHICfRunType = -1;
    mGeneratorIdx = -1;
    mDLEIdx = -1;

    mBTofMult = -1;
    memset(mBBCSumADC, 0., sizeof(mBBCSumADC));
    
    mTrkIdx.clear();
    mTrkparentId.clear();
    mTrkPid.clear();
    mTrkEnergy.clear();
    mTrkP[0].clear();
    mTrkP[1].clear();
    mTrkP[2].clear();
    mTrkVtxStart[0].clear();
    mTrkVtxStart[1].clear();
    mTrkVtxStart[2].clear();
    mTrkVtxEnd[0].clear();
    mTrkVtxEnd[1].clear();
    mTrkVtxEnd[2].clear();
    mTrkIsPrimary.clear();
    mTrkIsPropagate.clear();
    mTrkIsRHICfHit.clear();
    mTrkIsFinal.clear();

    mTowerIdx.clear();
    mPID.clear();
    mPos[0].clear();
    mPos[1].clear();
    mMomentum[0].clear();
    mMomentum[1].clear();
    mMomentum[2].clear();
    mMomentum[3].clear();
    mPi0Type.clear();
    mMass.clear();
    mRHICfTruthId[0].clear();
    mRHICfTruthId[1].clear();
    mRHICfTruthIncidentE[0].clear();
    mRHICfTruthIncidentE[1].clear();
    mRHICfTruthIncidentPosX[0].clear();
    mRHICfTruthIncidentPosX[1].clear();
    mRHICfTruthIncidentPosY[0].clear();
    mRHICfTruthIncidentPosY[1].clear();

    memset(mZDCPmtPhotonNum, 0., sizeof(mZDCPmtPhotonNum));
    memset(mZDCPmtdE, 0., sizeof(mZDCPmtdE));
    memset(mZDCSMDXdE, 0., sizeof(mZDCSMDXdE));
    memset(mZDCSMDYdE, 0., sizeof(mZDCSMDYdE));
    mZDCTruthId.clear();
    mZDCTruthIncidentE.clear();
    mZDCTruthIncidentPosX.clear();
    mZDCTruthIncidentPosY.clear();

    memset(mRHICfPlateE, 0., sizeof(mRHICfPlateE));
}

void StRHICfMiniSimDst::SetProcessID(int id){mProcessID = id;}
void StRHICfMiniSimDst::SetEventType(int type){mEventType = type;}
void StRHICfMiniSimDst::SetRHICfRunType(int type){mRHICfRunType = type;}
void StRHICfMiniSimDst::SetGeneratorIdx(int idx){mGeneratorIdx = idx;}
void StRHICfMiniSimDst::SetDLEIdx(int idx){mDLEIdx = idx;}

void StRHICfMiniSimDst::SetBTofMult(int mult){mBTofMult = mult;}
void StRHICfMiniSimDst::SetBBCSumADC(int ew, int sl, int adc){mBBCSumADC[ew][sl] = adc;}

void StRHICfMiniSimDst::SetSimTrack(int id, int parentId, int pid, double e, double px, double py, double pz, double vx, double vy, double vz, double ex, double ey, double ez, bool isPrimary, bool isPropagate, bool isRHICfHit, bool isFinal)
{
    mTrkIdx.push_back(id);
    mTrkparentId.push_back(parentId);
    mTrkPid.push_back(pid);
    mTrkEnergy.push_back(e);
    mTrkP[0].push_back(px);
    mTrkP[1].push_back(py);
    mTrkP[2].push_back(pz);
    mTrkVtxStart[0].push_back(vx);
    mTrkVtxStart[1].push_back(vy);
    mTrkVtxStart[2].push_back(vz);
    mTrkVtxEnd[0].push_back(ex);
    mTrkVtxEnd[1].push_back(ey);
    mTrkVtxEnd[2].push_back(ez);
    mTrkIsPrimary.push_back(isPrimary);
    mTrkIsPropagate.push_back(isPropagate);
    mTrkIsRHICfHit.push_back(isRHICfHit);
    mTrkIsFinal.push_back(isFinal);
}

void StRHICfMiniSimDst::SetRHICfParticle(int towerIdx, int pid, double x, double y, double px, double py, double pz, double e, int pi0Type, double mass)
{
    mTowerIdx.push_back(towerIdx);
    mPID.push_back(pid);
    mPos[0].push_back(x);
    mPos[1].push_back(y);
    mMomentum[0].push_back(px);
    mMomentum[1].push_back(py);
    mMomentum[2].push_back(pz);
    mMomentum[3].push_back(e);
    mPi0Type.push_back(pi0Type);
    mMass.push_back(mass);
}

void StRHICfMiniSimDst::SetRHICfTruthId(int tower, int idx, double posX, double posY, double e)
{
    mRHICfTruthId[tower].push_back(idx);
    mRHICfTruthIncidentE[tower].push_back(e);
    mRHICfTruthIncidentPosX[tower].push_back(posX);
    mRHICfTruthIncidentPosY[tower].push_back(posY);
}

void StRHICfMiniSimDst::SetZDCPmtPhotonNum(int pmt, int num){mZDCPmtPhotonNum[pmt] = num;}
void StRHICfMiniSimDst::SetZDCPmtdE(int idx, double e){mZDCPmtdE[idx] = e;}
void StRHICfMiniSimDst::SetZDCSMDdE(int xy, int smd, double e)
{
    if(xy == 0){mZDCSMDXdE[smd] = e;}
    if(xy == 1){mZDCSMDYdE[smd] = e;}
}
void StRHICfMiniSimDst::SetZDCTruthId(int idx, double posX, double posY, double e)
{
    mZDCTruthId.push_back(idx);
    mZDCTruthIncidentE.push_back(e);
    mZDCTruthIncidentPosX.push_back(posX);
    mZDCTruthIncidentPosY.push_back(posY);
}

void StRHICfMiniSimDst::SetRHICfPlatedE(int tower, int plate, double e){mRHICfPlateE[tower][plate] = e;}

Int_t StRHICfMiniSimDst::GetProcessID(){return mProcessID;}
Int_t StRHICfMiniSimDst::GetEventType(){return mEventType;}
Int_t StRHICfMiniSimDst::GetRHICfRunType(){return mRHICfRunType;}
Int_t StRHICfMiniSimDst::GetGeneratorIdx(){return mGeneratorIdx;}
Int_t StRHICfMiniSimDst::GetDLEIdx(){return mDLEIdx;}

Int_t StRHICfMiniSimDst::GetBTofMult(){return mBTofMult;}
Int_t StRHICfMiniSimDst::GetBBCSumADC(int ew, int sl){return mBBCSumADC[ew][sl];}

Int_t StRHICfMiniSimDst::GetSimTrackNum(){return mTrkIdx.size();}
Int_t StRHICfMiniSimDst::GetSimTrackId(int idx){return mTrkIdx[idx];}
Int_t StRHICfMiniSimDst::GetSimTrackParentId(int idx){return mTrkparentId[idx];}
Int_t StRHICfMiniSimDst::GetSimTrackPid(int idx){return mTrkPid[idx];}
Double_t StRHICfMiniSimDst::GetSimTrackEnergy(int idx){return mTrkEnergy[idx];}
Double_t StRHICfMiniSimDst::GetSimTrackPx(int idx){return mTrkP[0][idx];}
Double_t StRHICfMiniSimDst::GetSimTrackPy(int idx){return mTrkP[1][idx];}
Double_t StRHICfMiniSimDst::GetSimTrackPz(int idx){return mTrkP[2][idx];}
Double_t StRHICfMiniSimDst::GetSimTrackVxStart(int idx){return mTrkVtxStart[0][idx];}
Double_t StRHICfMiniSimDst::GetSimTrackVyStart(int idx){return mTrkVtxStart[1][idx];}
Double_t StRHICfMiniSimDst::GetSimTrackVzStart(int idx){return mTrkVtxStart[2][idx];}
Double_t StRHICfMiniSimDst::GetSimTrackVxEnd(int idx){return mTrkVtxEnd[0][idx];}
Double_t StRHICfMiniSimDst::GetSimTrackVyEnd(int idx){return mTrkVtxEnd[1][idx];}
Double_t StRHICfMiniSimDst::GetSimTrackVzEnd(int idx){return mTrkVtxEnd[2][idx];}
Bool_t StRHICfMiniSimDst::IsSimTrackPrimary(int idx){return mTrkIsPrimary[idx];}
Bool_t StRHICfMiniSimDst::IsSimTrackPropagate(int idx){return mTrkIsPropagate[idx];}
Bool_t StRHICfMiniSimDst::IsSimTrackRHICfHit(int idx){return mTrkIsRHICfHit[idx];}
Bool_t StRHICfMiniSimDst::IsSimTrackFinal(int idx){return mTrkIsFinal[idx];}

Int_t StRHICfMiniSimDst::GetParticleNum(){return mPID.size();}
Int_t StRHICfMiniSimDst::GetTowerIdx(int idx){return mTowerIdx[idx];}
Int_t StRHICfMiniSimDst::GetPID(int idx){return mPID[idx];}
Double_t StRHICfMiniSimDst::GetPosX(int idx){return mPos[0][idx];}
Double_t StRHICfMiniSimDst::GetPosY(int idx){return mPos[1][idx];}
Double_t StRHICfMiniSimDst::GetPx(int idx){return mMomentum[0][idx];}
Double_t StRHICfMiniSimDst::GetPy(int idx){return mMomentum[1][idx];}
Double_t StRHICfMiniSimDst::GetPz(int idx){return mMomentum[2][idx];}
Double_t StRHICfMiniSimDst::GetEnergy(int idx){return mMomentum[3][idx];}
Int_t StRHICfMiniSimDst::GetPi0Type(int idx){return mPi0Type[idx];}
Double_t StRHICfMiniSimDst::GetMass(int idx){return mMass[idx];}

Int_t StRHICfMiniSimDst::GetRHICfTruthNum(int tower){return mRHICfTruthId[tower].size();}
Int_t StRHICfMiniSimDst::GetRHICfTruthId(int tower, int idx){return mRHICfTruthId[tower][idx];}
Double_t StRHICfMiniSimDst::GetRHICfTruthIncidentEnergy(int tower, int idx){return mRHICfTruthIncidentE[tower][idx];}
Double_t StRHICfMiniSimDst::GetRHICfTruthIncidentPos(int tower, int idx, int xy)
{
    if(xy == 0){return  mRHICfTruthIncidentPosX[tower][idx];}
    if(xy == 1){return  mRHICfTruthIncidentPosY[tower][idx];}
    return -999.;
}

Int_t StRHICfMiniSimDst::GetZDCPmtPhotonNum(int idx){return mZDCPmtPhotonNum[idx];}
Double_t StRHICfMiniSimDst::GetZDCPmtdE(int idx){return mZDCPmtdE[idx];}
Double_t StRHICfMiniSimDst::GetZDCSMDdE(int xy, int idx)
{
    if(xy == 0){return mZDCSMDXdE[idx];}
    if(xy == 1){return mZDCSMDYdE[idx];}
    return -999.;
}
Int_t StRHICfMiniSimDst::GetZDCTruthNum(){return mZDCTruthId.size();}
Int_t StRHICfMiniSimDst::GetZDCTruthId(int idx){return mZDCTruthId[idx];}
Double_t StRHICfMiniSimDst::GetZDCTruthIncidentEnergy(int idx){return mZDCTruthIncidentE[idx];}
Double_t StRHICfMiniSimDst::GetZDCTruthIncidentPos(int idx, int xy)
{
    if(xy == 0){return  mZDCTruthIncidentPosX[idx];}
    if(xy == 1){return  mZDCTruthIncidentPosY[idx];}
    return -999.;
}

Double_t StRHICfMiniSimDst::GetRHICfPlatedE(int tower, int plate){return mRHICfPlateE[tower][plate];}