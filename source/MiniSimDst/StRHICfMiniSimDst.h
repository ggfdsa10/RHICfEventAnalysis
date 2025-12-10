#ifndef StRHICfMiniSimDst_HH
#define StRHICfMiniSimDst_HH

#include <vector>
#include "TObject.h"

using namespace std;

class StRHICfMiniSimDst : public TObject
{
    public: 
        StRHICfMiniSimDst();
        ~StRHICfMiniSimDst();

        void Clear(Option_t *option = "");

        // Event header
        void SetProcessID(int id);
        void SetEventType(int type);
        void SetRHICfRunType(int type);
        void SetGeneratorIdx(int idx);
        void SetDLEIdx(int idx);

        void SetBTofMult(int mult);
        void SetBBCSumADC(int ew, int sl, int adc);

        void SetSimTrack(int id, int parentId, int pid, double e, double px, double py, double pz, double vx, double vy, double vz, double ex, double ey, double ez, bool isPrimary, bool isPropagate, bool isRHICfHit, bool isFinal);

        void SetRHICfParticle(int towerIdx, int pid, double x, double y, double px, double py, double pz, double e, int pi0Type=-1, double mass=-1.);
        void SetRHICfTruthId(int tower, int idx, double posX, double posY, double e);

        void SetZDCPmtPhotonNum(int pmt, int num);
        void SetZDCPmtdE(int idx, double e);
        void SetZDCSMDdE(int xy, int smd, double e);
        void SetZDCTruthId(int idx, double posX, double posY, double e);

        void SetRHICfPlatedE(int tower, int plate, double e);

        Int_t GetProcessID();
        Int_t GetEventType();
        Int_t GetRHICfRunType();
        Int_t GetGeneratorIdx();
        Int_t GetDLEIdx();

        Int_t GetBTofMult();
        Int_t GetBBCSumADC(int ew, int sl);

        Int_t GetSimTrackNum();
        Int_t GetSimTrackId(int idx);
        Int_t GetSimTrackParentId(int idx);
        Int_t GetSimTrackPid(int idx);
        Double_t GetSimTrackEnergy(int idx);
        Double_t GetSimTrackPx(int idx);
        Double_t GetSimTrackPy(int idx);
        Double_t GetSimTrackPz(int idx);
        Double_t GetSimTrackVxStart(int idx);
        Double_t GetSimTrackVyStart(int idx);
        Double_t GetSimTrackVzStart(int idx);
        Double_t GetSimTrackVxEnd(int idx);
        Double_t GetSimTrackVyEnd(int idx);
        Double_t GetSimTrackVzEnd(int idx);
        Bool_t IsSimTrackPrimary(int idx);
        Bool_t IsSimTrackPropagate(int idx);
        Bool_t IsSimTrackRHICfHit(int idx);
        Bool_t IsSimTrackFinal(int idx);
        
        Int_t GetParticleNum();
        Int_t GetTowerIdx(int idx=0);
        Int_t GetPID(int idx=0);
        Double_t GetPosX(int idx=0);
        Double_t GetPosY(int idx=0);
        Double_t GetPx(int idx=0);
        Double_t GetPy(int idx=0);
        Double_t GetPz(int idx=0);
        Double_t GetEnergy(int idx=0);
        Int_t GetPi0Type(int idx=0);
        Double_t GetMass(int idx=0);

        Int_t GetRHICfTruthNum(int tower);
        Int_t GetRHICfTruthId(int tower, int idx);
        Double_t GetRHICfTruthIncidentEnergy(int tower, int idx);
        Double_t GetRHICfTruthIncidentPos(int tower, int idx, int xy);

        Int_t GetZDCPmtPhotonNum(int idx);
        Double_t GetZDCPmtdE(int idx);
        Double_t GetZDCSMDdE(int xy, int idx);

        Int_t GetZDCTruthNum();
        Int_t GetZDCTruthId(int idx);
        Double_t GetZDCTruthIncidentEnergy(int idx);
        Double_t GetZDCTruthIncidentPos(int idx, int xy);

        Double_t GetRHICfPlatedE(int tower, int plate);

    private:
        Int_t mProcessID;
        Char_t mEventType;
        Char_t mRHICfRunType;
        Char_t mGeneratorIdx;
        Char_t mDLEIdx;

        // STAR detector data
        Int_t mBBCSumADC[2][2]; // [east, west][small, large];
        Int_t mBTofMult;

        // Truth track data
        vector<int> mTrkIdx;
        vector<int> mTrkparentId;
        vector<int> mTrkPid;
        vector<double> mTrkEnergy;
        vector<double> mTrkP[3];
        vector<double> mTrkVtxStart[3];
        vector<double> mTrkVtxEnd[3];
        vector<bool> mTrkIsPrimary;
        vector<bool> mTrkIsPropagate;
        vector<bool> mTrkIsRHICfHit;
        vector<bool> mTrkIsFinal;

        // RHICf data
        vector<char> mTowerIdx;
        vector<char> mPID; // [0 == pi0, 1 == neutron]
        vector<double> mPos[2];
        vector<double> mMomentum[4];
        vector<char> mPi0Type;
        vector<double> mMass;
        vector<int> mRHICfTruthId[2];
        vector<double> mRHICfTruthIncidentE[2];
        vector<double> mRHICfTruthIncidentPosX[2];
        vector<double> mRHICfTruthIncidentPosY[2];
        
        Int_t mZDCPmtPhotonNum[3];
        Double_t mZDCPmtdE[3];
        Double_t mZDCSMDXdE[8];
        Double_t mZDCSMDYdE[7];
        vector<int> mZDCTruthId;
        vector<double> mZDCTruthIncidentE;
        vector<double> mZDCTruthIncidentPosX;
        vector<double> mZDCTruthIncidentPosY;
        
        Double_t mRHICfPlateE[2][16];

    ClassDef(StRHICfMiniSimDst,1)
};

#endif
