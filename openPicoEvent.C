const int kRunNum = 3;
const TString kRunName[kRunNum] = {"TL", "TS", "TOP"};
const TString kPicoEventPath = "/gpfs01/star/pwg/slee5/RHICfPicoEvent";

TFile* mFiles[kRunNum];
TTree* mTrees[kRunNum];

// Run info
Int_t mRunNum;
Int_t mFillNum;
Int_t mRHICfRunType;

// RHICf trigger 
Bool_t mRHICfShowerTrig;
Bool_t mRHICfPi0Trig;
Bool_t mRHICfHighEMTrig;

// STAR trigger 
Bool_t mRHICfTrig;
Bool_t mRHICfDiffractiveTrig;
Bool_t mRHICfVPDMC30Trig;
Bool_t mRHICfTPCTrig;
Bool_t mRPSDTTrig;
Bool_t mRPETTrig;
Bool_t mRPCPTnoBBCLTrig;
Bool_t mRPCPT2Trig;
Bool_t mRPCPT2noBBCLTrig;
Bool_t mVPDMB30Trig;
Bool_t mVPDMB100Trig;
Bool_t mBHT3Trig;
Bool_t mTofHighMultTrig;
Bool_t mEPDTrig;
Bool_t mBBCTrig;
Bool_t mBBCTacTrig;
Bool_t mZDCTrig;
Bool_t mZDCTacTrig;
Bool_t mVPDTrig;
Bool_t mZeroBiasTrig;

// star detectors
Int_t mTOFMult;
Int_t mBBCSmallADC[2][16]; // [east, west][BBC Small tower pmt index]
Int_t mBBCSmallTDC[2][16]; // [east, west][BBC Small tower pmt index]
Int_t mBBCSmallTAC[2][16]; // [east, west][BBC Small tower pmt index]
Int_t mBBCLargeADC[2][8];  // [east, west][BBC Large tower pmt index]
Int_t mBBCLargeTDC[2][8];  // [east, west][BBC Large tower pmt index]
Int_t mBBCLargeTAC[2][8];  // [east, west][BBC Large tower pmt index]

Int_t mZDCADC[2][3]; // [east, west][ZDC pmt index]
Int_t mZDCTDC[2][3]; // [east, west][ZDC pmt index]
Int_t mZDCSumADC[2][2]; // [east, west][Attenuated, Unattenuated]
Int_t mZDCSMD[2][2][7];  //[east, west][vertical, horizontal][smd pmt index]

// rhicf detector 
Double_t mL20[2]; // [TS, TL]
Double_t mL90[2]; // [TS, TL]
Double_t mPlateEnergy[2][16]; // [TS, TL][plate layers]
Int_t mGSObarMaxLayer[2][2]; // [TS, TL][1st Max layer, 2nd Max layer]
Int_t mResultHitNum[2]; // [TS, TL]
Int_t mEvalHitNum[2][4][2]; //[TS, TL][GSOBar layers][x, y]
Double_t mSingleHit[2][4][2][2]; // [TS, TL][GSOBar layers][x, y][pos, height]
Double_t mMultiHit[2][4][2][2][2]; // [TS, TL][GSOBar layers][x, y][multi order][pos, height]

// rhicf point
Int_t mPointNum;
Int_t mPointTowerIdx[4]; // [num] (0 = TS, 1 = TL)
Int_t mPointPID[4]; // [num] (0 = photon, 1 = hadron)
Double_t mPointPos[4][2]; // [num][x, y]
Double_t mPointEnergy[4][2]; // [num][photon, hadron]

void openPicoEvent()
{
    for(int run=0; run<kRunNum; run++){
        mFiles[run] = new TFile(Form("%s/RHICfPicoEvent_%s.root", kPicoEventPath.Data(), kRunName[run].Data()), "read");
        mTrees[run] = (TTree*)mFiles[run] -> Get("event");
        mTrees[run] -> SetBranchAddress("RunNum", &mRunNum);
        mTrees[run] -> SetBranchAddress("FillNum", &mFillNum);
        mTrees[run] -> SetBranchAddress("RHICfShowerTrig", &mRHICfShowerTrig);
        mTrees[run] -> SetBranchAddress("RHICfPi0Trig", &mRHICfPi0Trig);
        mTrees[run] -> SetBranchAddress("RHICfHighEMTrig", &mRHICfHighEMTrig);
        mTrees[run] -> SetBranchAddress("RHICfRunType", &mRHICfRunType);
        mTrees[run] -> SetBranchAddress("RHICfTrig", &mRHICfTrig);
        mTrees[run] -> SetBranchAddress("RHICfDiffractiveTrig", &mRHICfDiffractiveTrig);
        mTrees[run] -> SetBranchAddress("RHICfVPDMC30Trig", &mRHICfVPDMC30Trig);
        mTrees[run] -> SetBranchAddress("RHICfTPCTrig", &mRHICfTPCTrig);
        mTrees[run] -> SetBranchAddress("RPSDTTrig", &mRPSDTTrig);
        mTrees[run] -> SetBranchAddress("RPETTTrig", &mRPETTrig);
        mTrees[run] -> SetBranchAddress("RPPTnoBBCLTTrig", &mRPCPTnoBBCLTrig);
        mTrees[run] -> SetBranchAddress("RPPT2TTrig", &mRPCPT2Trig);
        mTrees[run] -> SetBranchAddress("RPSPT2noBBCLTrig", &mRPCPT2noBBCLTrig);
        mTrees[run] -> SetBranchAddress("VPDMB30Trig", &mVPDMB30Trig);
        mTrees[run] -> SetBranchAddress("VPDMB100Trig", &mVPDMB100Trig);
        mTrees[run] -> SetBranchAddress("BHT3Trig", &mBHT3Trig);
        mTrees[run] -> SetBranchAddress("TofHighMultTrig", &mTofHighMultTrig);
        mTrees[run] -> SetBranchAddress("EPDTrig", &mEPDTrig);
        mTrees[run] -> SetBranchAddress("BBCTrig", &mBBCTrig);
        mTrees[run] -> SetBranchAddress("BBCTacTrig", &mBBCTacTrig);
        mTrees[run] -> SetBranchAddress("ZDCTrig", &mZDCTrig);
        mTrees[run] -> SetBranchAddress("ZDCTacTrig", &mZDCTacTrig);
        mTrees[run] -> SetBranchAddress("VPDTrig", &mVPDTrig);
        mTrees[run] -> SetBranchAddress("ZeroBiasTrig", &mZeroBiasTrig);
        mTrees[run] -> SetBranchAddress("TOFMult", &mTOFMult);
        mTrees[run] -> SetBranchAddress("BBCSmallADC", &mBBCSmallADC);
        mTrees[run] -> SetBranchAddress("BBCSmallTDC", &mBBCSmallTDC);
        mTrees[run] -> SetBranchAddress("BBCSmallTAC", &mBBCSmallTAC);
        mTrees[run] -> SetBranchAddress("BBCLargeADC", &mBBCLargeADC);
        mTrees[run] -> SetBranchAddress("BBCLargeTDC", &mBBCLargeTDC);
        mTrees[run] -> SetBranchAddress("BBCLargeTAC", &mBBCLargeTAC);
        mTrees[run] -> SetBranchAddress("ZDCADC", &mZDCADC);
        mTrees[run] -> SetBranchAddress("ZDCTDC", &mZDCTDC);
        mTrees[run] -> SetBranchAddress("ZDCSumADC", &mZDCSumADC);
        mTrees[run] -> SetBranchAddress("ZDCSMD", &mZDCSMD);
        mTrees[run] -> SetBranchAddress("L20", &mL20);
        mTrees[run] -> SetBranchAddress("L90", &mL90);
        mTrees[run] -> SetBranchAddress("PlateEnergy", &mPlateEnergy);
        mTrees[run] -> SetBranchAddress("GSObarMaxLayer", &mGSObarMaxLayer);
        mTrees[run] -> SetBranchAddress("ResultHitNum", &mResultHitNum);
        mTrees[run] -> SetBranchAddress("EvalHitNum", &mEvalHitNum);
        mTrees[run] -> SetBranchAddress("SingleHit", &mSingleHit);
        mTrees[run] -> SetBranchAddress("MultiHit", &mMultiHit);
        mTrees[run] -> SetBranchAddress("PointNum", &mPointNum);
        mTrees[run] -> SetBranchAddress("PointTowerIdx", &mPointTowerIdx);
        mTrees[run] -> SetBranchAddress("PointPID", &mPointPID);
        mTrees[run] -> SetBranchAddress("PointPos", &mPointPos);
        mTrees[run] -> SetBranchAddress("PointEnergy", &mPointEnergy);

        int eventNum = mTrees[run] -> GetEntries();

        cout << "========================== RHICf " << kRunName[run] << " run : " << eventNum << " events " << endl;

        for(int event=0; event<eventNum; event++){
            if(event%1000 == 0){cout << " event: " << event << " / " << eventNum << endl;}
            mTrees[run] -> GetEntry(event);

            if(!mRHICfTrig){continue;}

            if(!mRHICfHighEMTrig){continue;}

            // ====================== events =========================
            bool isNeutronLikeEvent = false;
            for(int p=0; p<mPointNum; p++){
                int pid = mPointPID[p]; // 0 = photon, 1 = hadron
                int towerIdx = mPointTowerIdx[p]; // 0 = TS, 1 = TL
                double l90 = mL90[towerIdx];
                double energy = mPointEnergy[p][0]; // [GeV]

                if(l90 > 20 && pid == 1 && energy > 20){
                    isNeutronLikeEvent = true;
                }
            }
            if(!isNeutronLikeEvent){continue;}

            int zdcEastSum = 0;
            int zdcWestSum = 0;
            for(int i=0; i<3; i++){
                zdcEastSum += mZDCADC[0][i];
                zdcWestSum += mZDCADC[1][i];
            }

        }

    }


}
