#ifndef RHICfOptContainer_hh
#define RHICfOptContainer_hh

#include <iostream>
#include "TString.h"

enum parameters{
    kTLRun = 0,
    kTSRun = 1,
    kTOPRun = 2,
    kALLRun = -1,

    kPhoton = 0,
    kHadron = 1,

    kSDLE = 0,
    kDDLE = 1,
    kNDLE = 2,
    kALLDLE = 3,

    kGammaRun = 0,
    kPi0Run = 1,
    kNeutronRun = 2,
    kLambda0Run = 3,

    kType1 = 0,
    kType2 = 1,

    kNULL = 0
};

static const double kRHICfPosZ = 17800.; // [mm]
static const int kFillNum = 5;
static const int kRunNum = 3; 
static const int kTypeNum = 2;
static const int kDLENum = 4;

using namespace std;

class RHICfOptContainer
{
    public:
        static RHICfOptContainer* GetOptContainer();

        RHICfOptContainer();
        ~RHICfOptContainer();

        virtual int Init();

        // ========== Set function for options ==========
        void SetRunType(TString type);
        void SetInputDataPath(TString path);
        void SetInputDataList(TString listFile);
        void SetExecuteEventNum(int event);
        void SetExecuteFileNum(int num);

        void ForceCalculateMass();
        void ForceCalculateBinning();
        void ForceCalculateDilution();
        void ForceCalculateAsymmetry();
        void ForceCalculateSystematicError();

        void ForceDefaultBinning();

        // specific particle calculation option
        void CalculateGamma();
        void CalculatePi0();
        void CalculateNeutron();
        void CalculateLambda0();

        // ========== On-Off data options for RHICfEventDst
        void SetOffDetPoint();
        void SetOffParticle();
        void SetOffTPCTrack();
        void SetOffBTof();
        void SetOffBBC();
        void SetOffVPD();
        void SetOffZDC();
        void SetOffFMS();
        void SetOffRPS();

        // ========== Get function for options ==========
        int GetRunType();
        int GetFillNumIdx(int fillNum);
        int GetFillNumToRunIdx(int fillNum);
        int GetFillToRunIdx(int fillIdx);

        TString GetRunTypeName(int runIdx=-1);
        TString GetDLEName(int dleIdx);

        TString GetInputDataPath();
        TString GetInputDataList();
        TString GetTablePath();
        TString GetDataPath();
        TString GetFigurePath();
        int GetExecuteEventNum();
        int GetExecuteFileNum();

        bool IsForceCalculateMass();
        bool IsForceCalculateBinning();
        bool IsForceCalculateDilution();
        bool IsForceCalculateAsymmetry();
        bool IsForceCalculateSystematicError();

        bool IsForceDefaultBinning();

        bool GetOffDetPoint();
        bool GetOffParticle();
        bool GetOffTPCTrack();
        bool GetOffBTof();
        bool GetOffBBC();
        bool GetOffVPD();
        bool GetOffZDC();
        bool GetOffFMS();
        bool GetOffRPS();

        TString GetParticleRunName();
        Int_t GetParticleRunIdx();

    private:
        void InitOptions();
        void FindDirPath();
        void PrintOptions();

        static RHICfOptContainer* mOptContainer;
        
        int mRunType;
        TString mInputDataPath;
        TString mInputDataList;
        TString mTablePath;
        TString mDataPath;
        TString mFigurePath;
        int mExecuteEventNum;
        int mExecuteFileNum;

        bool mGammaCalc;
        bool mPi0Calc;
        bool mNeutronCalc;
        bool mLambda0Calc;

        bool mForceCalculateMass;
        bool mForceCalculateBinning;
        bool mForceCalculateDilution;
        bool mForceCalculateAsymmetry;
        bool mForceCalculateSystematicError;
        
        bool mForceDefaultBinning;

        bool mOffDetPoint;
        bool mOffParticle;
        bool mOffTPCTrack;
        bool mOffBTof;
        bool mOffBBC;
        bool mOffVPD;
        bool mOffZDC;
        bool mOffFMS;
        bool mOffRPS;

};

#endif
