#include "RHICfPicoEventMaker.hh"

RHICfPicoEventMaker::RHICfPicoEventMaker() 
: mOutputPath("/gpfs01/star/pwg/slee5/RHICfPicoEvent")
{
    mIsRHICfStream = false;
    mIsPhysicsStream = false;
    mIsRPSEvent = false;
    mIsRHICfParticleEvent = false;
}

RHICfPicoEventMaker::~RHICfPicoEventMaker()
{
}

int RHICfPicoEventMaker::Init()
{
    if(mIsPhysicsStream){GetOptContainer()->SetTestMode();} // Not used particle maker
    if(mIsRHICfStream){GetOptContainer()->SetTestMode();} // Not used particle maker
    if(!RHICfOptContainer::Init()){return 0;}
    if(!RHICfTableMaker::Init()){return 0;}

    // Utilized classes
    mEventReader = new RHICfEventDstReader();
    mParticleMaker = new RHICfParticleMaker();
    mDLECondition = new RHICfDLECondition();

    mEventReader -> Init();
    mParticleMaker -> Init(); 
    mDLECondition -> Init();  

    return 1;
}

int RHICfPicoEventMaker::Calculate()
{
    cout << "RHICfPicoEventMaker::Calculate() -- Main calculation has started. " << endl;

    if(mIsRHICfStream){RHICfStream();}
    if(mIsPhysicsStream){PhysicsStream();}
    if(mIsRPSEvent){RPSEvent();}
    if(mIsRHICfParticleEvent){RHICfParticleEvent();}

    cout << "RHICfPicoEventMaker::Calculate() -- Done. " << endl;
    return 1;
}

int RHICfPicoEventMaker::Finish()
{
    return 1;
}

void RHICfPicoEventMaker::RHICfStream()
{
    int runType = GetOptContainer()->GetRunType();
    TString runName = GetOptContainer()->GetRunTypeName(runType);
    TFile* file = new TFile(Form("%s/RHICfPicoEvent_%s.root", mOutputPath.Data(), runName.Data()), "recreate");
    TTree* tree = new TTree("event", "event");

    Int_t outRunNum;
    Int_t outFillNum;
    Int_t outRHICfRunType;

    // RHICf standalone trigger
    Bool_t outRHICfShowerTrig;
    Bool_t outRHICfPi0Trig;
    Bool_t outRHICfHighEMTrig;

    // STAR trigger
    Bool_t outRHICfTrig;
    Bool_t outRHICfDiffractiveTrig;
    Bool_t outRHICfVPDMC30Trig;
    Bool_t outRHICfTPCTrig;
    Bool_t outRPSDTTrig;
    Bool_t outRPETTrig;
    Bool_t outRPCPTnoBBCLTrig;
    Bool_t outRPCPT2Trig;
    Bool_t outRPCPT2noBBCLTrig;
    Bool_t outVPDMB30Trig;
    Bool_t outVPDMB100Trig;
    Bool_t outBHT3Trig;
    Bool_t outTofHighMultTrig;
    Bool_t outEPDTrig;
    Bool_t outBBCTrig;
    Bool_t outBBCTacTrig;
    Bool_t outZDCTrig;
    Bool_t outZDCTacTrig;
    Bool_t outVPDTrig;
    Bool_t outZeroBiasTrig;

    Int_t outTOFMult;
    Int_t outBBCSmallADC[kBeamSideNum][kBBCSmallPmtNum];
    Int_t outBBCSmallTDC[kBeamSideNum][kBBCSmallPmtNum];
    Int_t outBBCSmallTAC[kBeamSideNum][kBBCSmallPmtNum];
    Int_t outBBCLargeADC[kBeamSideNum][kBBCLargePmtNum];
    Int_t outBBCLargeTDC[kBeamSideNum][kBBCLargePmtNum];
    Int_t outBBCLargeTAC[kBeamSideNum][kBBCLargePmtNum];

    Int_t outZDCADC[kBeamSideNum][kZDCPmtNum];
    Int_t outZDCTDC[kBeamSideNum][kZDCPmtNum];
    Int_t outZDCSumADC[kBeamSideNum][2]; // [Attenuated, Unattenuated]
    Int_t outZDCSMD[kBeamSideNum][2][kSMDNum];  //[east, west][vertical, horizontal][smd]

    Double_t outL20[kTowerNum];
    Double_t outL90[kTowerNum];
    Double_t outPlateEnergy[kTowerNum][kPlateNum];
    Int_t outGSObarMaxLayer[kTowerNum][2]; // [1st Max layer, 2nd Max layer]
    Int_t outResultHitNum[kTowerNum];
    Int_t outEvalHitNum[kTowerNum][kLayerNum][kXYNum];
    Double_t outSingleHit[kTowerNum][kLayerNum][kXYNum][2]; // [pos, height]
    Double_t outMultiHit[kTowerNum][kLayerNum][kXYNum][2][2]; // [multi order][pos, height]

    Int_t outPointNum;
    Int_t outPointTowerIdx[4]; // [num]
    Int_t outPointPID[4]; // [num]
    Double_t outPointPos[4][2]; // [num][x, y]
    Double_t outPointEnergy[4][2]; // [num][photon, hadron]

    tree -> Branch("RunNum", &outRunNum, "RunNum/I");
    tree -> Branch("FillNum", &outFillNum, "FillNum/I");
    tree -> Branch("RHICfRunType", &outRHICfRunType, "RHICfRunType/I");

    tree -> Branch("RHICfShowerTrig", &outRHICfShowerTrig, "RHICfShowerTrig/O");
    tree -> Branch("RHICfPi0Trig", &outRHICfPi0Trig, "RHICfPi0Trig/O");
    tree -> Branch("RHICfHighEMTrig", &outRHICfHighEMTrig, "RHICfHighEMTrig/O");

    tree -> Branch("RHICfTrig", &outRHICfTrig, "RHICfTrig/O");
    tree -> Branch("RHICfDiffractiveTrig", &outRHICfDiffractiveTrig, "RHICfDiffractiveTrig/O");
    tree -> Branch("RHICfVPDMC30Trig", &outRHICfVPDMC30Trig, "RHICfVPDMC30Trig/O");
    tree -> Branch("RHICfTPCTrig", &outRHICfTPCTrig, "RHICfTPCTrig/O");
    tree -> Branch("RPSDTTrig", &outRPSDTTrig, "RPSDTTrig/O");
    tree -> Branch("RPETTTrig", &outRPETTrig, "RPETTTrig/O");
    tree -> Branch("RPPTnoBBCLTTrig", &outRPCPTnoBBCLTrig, "RPPTnoBBCLTTrig/O");
    tree -> Branch("RPPT2TTrig", &outRPCPT2Trig, "RPPT2TTrig/O");
    tree -> Branch("RPSPT2noBBCLTrig", &outRPCPT2noBBCLTrig, "RPSPT2noBBCLTrig/O");
    tree -> Branch("VPDMB30Trig", &outVPDMB30Trig, "VPDMB30Trig/O");
    tree -> Branch("VPDMB100Trig", &outVPDMB100Trig, "VPDMB100Trig/O");
    tree -> Branch("BHT3Trig", &outBHT3Trig, "BHT3Trig/O");
    tree -> Branch("TofHighMultTrig", &outTofHighMultTrig, "TofHighMultTrig/O");
    tree -> Branch("EPDTrig", &outEPDTrig, "EPDTrig/O");
    tree -> Branch("BBCTrig", &outBBCTrig, "BBCTrig/O");
    tree -> Branch("BBCTacTrig", &outBBCTacTrig, "BBCTacTrig/O");
    tree -> Branch("ZDCTrig", &outZDCTrig, "ZDCTrig/O");
    tree -> Branch("ZDCTacTrig", &outZDCTacTrig, "ZDCTacTrig/O");
    tree -> Branch("VPDTrig", &outVPDTrig, "VPDTrig/O");
    tree -> Branch("ZeroBiasTrig", &outZeroBiasTrig, "ZeroBiasTrig/O");

    tree -> Branch("TOFMult", &outTOFMult, "TOFMult/I");
    tree -> Branch("BBCSmallADC", outBBCSmallADC, "BBCSmallADC[2][16]/I");
    tree -> Branch("BBCSmallTDC", outBBCSmallTDC, "BBCSmallTDC[2][16]/I");
    tree -> Branch("BBCSmallTAC", outBBCSmallTAC, "BBCSmallTAC[2][16]/I");
    tree -> Branch("BBCLargeADC", outBBCLargeADC, "BBCLargeADC[2][8]/I");
    tree -> Branch("BBCLargeTDC", outBBCLargeTDC, "BBCLargeTDC[2][8]/I");
    tree -> Branch("BBCLargeTAC", outBBCLargeTAC, "BBCLargeTAC[2][8]/I");

    tree -> Branch("ZDCADC", outZDCADC, "ZDCADC[2][3]/I");
    tree -> Branch("ZDCTDC", outZDCTDC, "ZDCTDC[2][3]/I");
    tree -> Branch("ZDCSumADC", outZDCSumADC, "ZDCSumADC[2][2]/I");
    tree -> Branch("ZDCSMD", outZDCSMD, "ZDCSMD[2][2][7]/I");

    tree -> Branch("L20", outL20, "L20[2]/D");
    tree -> Branch("L90", outL90, "L90[2]/D");
    tree -> Branch("PlateEnergy", outPlateEnergy, "PlateEnergy[2][16]/D");
    tree -> Branch("GSObarMaxLayer", outGSObarMaxLayer, "GSObarMaxLayer[2][2]/I");
    tree -> Branch("ResultHitNum", outResultHitNum, "ResultHitNum[2]/I");
    tree -> Branch("EvalHitNum", outEvalHitNum, "EvalHitNum[2][4][2]/I");
    tree -> Branch("SingleHit", outSingleHit, "SingleHit[2][4][2][2]/D");
    tree -> Branch("MultiHit", outMultiHit, "MultiHit[2][4][2][2][2]/D");

    tree -> Branch("PointNum", &outPointNum, "PointNum/I");
    tree -> Branch("PointTowerIdx", outPointTowerIdx, "PointTowerIdx[PointNum]/I");
    tree -> Branch("PointPID", outPointPID, "PointPID[PointNum]/I");
    tree -> Branch("PointPos", outPointPos, "outPointPos[PointNum][2]/D");
    tree -> Branch("PointEnergy", outPointEnergy, "PointEnergy[PointNum][2]/D");


    int eventNum = mEventReader->GetEventNum();
    for(int event=0; event<eventNum; event++){
        if(event%10000 == 0){cout << "RHICfPicoEventMaker::ParticleTest() -- Event: " << event << " / " << eventNum << endl;}
        mEventDst = mEventReader -> GetEventDst(event);

        fill_n(&outBBCSmallADC[0][0], kBeamSideNum*kBBCSmallPmtNum, -1);
        fill_n(&outBBCSmallTDC[0][0], kBeamSideNum*kBBCSmallPmtNum, -1);
        fill_n(&outBBCSmallTAC[0][0], kBeamSideNum*kBBCSmallPmtNum, -1);
        fill_n(&outBBCLargeADC[0][0], kBeamSideNum*kBBCLargePmtNum, -1);
        fill_n(&outBBCLargeTDC[0][0], kBeamSideNum*kBBCLargePmtNum, -1);
        fill_n(&outBBCLargeTAC[0][0], kBeamSideNum*kBBCLargePmtNum, -1);
        fill_n(&outZDCADC[0][0], kBeamSideNum*kZDCPmtNum, -1);
        fill_n(&outZDCTDC[0][0], kBeamSideNum*kZDCPmtNum, -1);
        fill_n(&outZDCSumADC[0][0], kBeamSideNum*2, -1);
        fill_n(&outZDCSMD[0][0][0], kBeamSideNum*2*kSMDNum, -1);

        fill_n(&outL20[0], kTowerNum, -1);
        fill_n(&outL90[0], kTowerNum, -1);
        fill_n(&outPlateEnergy[0][0], kTowerNum*kPlateNum, -1);
        fill_n(&outGSObarMaxLayer[0][0], kTowerNum*2, -1);
        fill_n(&outResultHitNum[0], kTowerNum, -1);
        fill_n(&outEvalHitNum[0][0][0], kTowerNum*kLayerNum*kXYNum, -1);
        fill_n(&outSingleHit[0][0][0][0], kTowerNum*kLayerNum*kXYNum*2, -1);
        fill_n(&outMultiHit[0][0][0][0][0], kTowerNum*kLayerNum*kXYNum*2*2, -1);

        outPointNum = 0;

        fill_n(&outPointTowerIdx[0], 4, -1);
        fill_n(&outPointPID[0], 4, -1);
        fill_n(&outPointPos[0][0], 4*2, -1);
        fill_n(&outPointEnergy[0][0], 4*2, -1);

        mEvent = mEventDst -> GetEvent();
        outFillNum = mEvent -> GetFillNumber();
        outRunNum = mEvent -> GetRunNumber();
        outRHICfRunType = mEvent -> GetRHICfRunType();

        outRHICfShowerTrig = mEvent -> GetRHICfShowerTrig();
        outRHICfPi0Trig = mEvent -> GetRHICfPi0Trig();
        outRHICfHighEMTrig = mEvent -> GetRHICfHighEMTrig();

        outRHICfTrig = mEvent -> GetRHICfTrig();
        outRHICfDiffractiveTrig = mEvent -> GetRHICfDiffractiveTrig();
        outRHICfVPDMC30Trig = mEvent -> GetRHICfVPDMB30Trig();
        outRHICfTPCTrig = mEvent -> GetRHICfTPCTrig();
        outRPSDTTrig = mEvent -> GetRPSDTTrig();
        outRPETTrig = mEvent -> GetRPETTrig();
        outRPCPTnoBBCLTrig = mEvent -> GetRPCPTnoBBCLTrig();
        outRPCPT2Trig = mEvent -> GetRPCPT2Trig();
        outRPCPT2noBBCLTrig = mEvent -> GetRPCPT2noBBCLTrig();
        outVPDMB30Trig = mEvent -> GetVPDMB30Trig();
        outVPDMB100Trig = mEvent -> GetVPDMB100Trig();
        outBHT3Trig = mEvent -> GetBHT3Trig();
        outTofHighMultTrig = mEvent -> GetTofHighMultTrig();
        outEPDTrig = mEvent -> GetEPDTrig();
        outBBCTrig = mEvent -> GetBBCTrig();
        outBBCTacTrig = mEvent -> GetBBCTacTrig();
        outZDCTrig = mEvent -> GetZDCTrig();
        outZDCTacTrig = mEvent -> GetZDCTacTrig();
        outVPDTrig = mEvent -> GetVPDTrig();
        outZeroBiasTrig = mEvent -> GetZeroBiasTrig();

        outTOFMult = mEvent -> GetBTofMult();
        mBBC = mEventDst -> GetBBC();
        mZDC = mEventDst -> GetZDC();

        for(int ew=0; ew<kBeamSideNum; ew++){
            int beamSide = (ew==0)? kBeamEast : kBeamWest;

            for(int i=0; i<kBBCSmallPmtNum; i++){
                outBBCSmallADC[beamSide][i] = mBBC -> GetBBCSmallADC(beamSide, i);
                outBBCSmallTDC[beamSide][i] = mBBC -> GetBBCSmallTDC(beamSide, i);
                outBBCSmallTAC[beamSide][i] = mBBC -> GetBBCSmallTAC(beamSide, i);
            }
            for(int i=0; i<kBBCLargePmtNum; i++){
                outBBCLargeADC[beamSide][i] = mBBC -> GetBBCLargeADC(beamSide, i);
                outBBCLargeTDC[beamSide][i] = mBBC -> GetBBCLargeTDC(beamSide, i);
                outBBCLargeTAC[beamSide][i] = mBBC -> GetBBCLargeTAC(beamSide, i);
            }

            for(int i=0; i<kZDCPmtNum; i++){
                outZDCADC[beamSide][i] = mZDC -> GetZDCPmtADC(beamSide, i);
                outZDCTDC[beamSide][i] = mZDC -> GetZDCPmtTDC(beamSide, i);
            }
            outZDCSumADC[beamSide][0] = mZDC -> GetZDCPmtAttenuatedSumADC(beamSide);
            outZDCSumADC[beamSide][1] = mZDC -> GetZDCPmtUnAttenuatedSumADC(beamSide);
        }

        for(int ew=0; ew<kBeamSideNum; ew++){
            for(int ixy=0; ixy<2; ixy++){
                for(int smd=0; smd<7; smd++){
                    outZDCSMD[ew][ixy][smd] = mZDC -> GetZDCSmdADC(ew, ixy, smd);
                }
            }
        }

        // RHICf detector hit information 
        mDetHit = mEventDst -> GetRHICfDetHit();

        for(int it=0; it<kTowerNum; it++){
            outL20[it] = mDetHit -> GetL20(it);
            outL90[it] = mDetHit -> GetL90(it);
            outResultHitNum[it] = mDetHit -> GetResultHitNum(it);

            for(int ip=0; ip<kPlateNum; ip++){
                outPlateEnergy[it][ip] = mDetHit -> GetPlateEnergy(it, ip);
            }

            for(int io=0; io<2; io++){
                outGSObarMaxLayer[it][io] = mDetHit -> GetGSOBarMaxLayer(it, io);
            }

            for(int il=0; il<kLayerNum; il++){
                for(int ixy=0; ixy<kXYNum; ixy++){
                    outEvalHitNum[it][il][ixy] = mDetHit -> GetEvalHitNum(it, il, ixy);
                    outSingleHit[it][il][ixy][0] = mDetHit -> GetSingleHitPos(it, il, ixy);
                    outSingleHit[it][il][ixy][1] = mDetHit -> GetSingleHitHeight(it, il, ixy);

                    for(int io=0; io<2; io++){
                        outMultiHit[it][il][ixy][io][0] = mDetHit -> GetMultiHitPos(it, il, ixy, io);
                        outMultiHit[it][il][ixy][io][1] = mDetHit -> GetMultiHitHeight(it, il, ixy, io);
                    }
                }
            }
        }

        // RHICf Point information
        outPointNum = mEventDst -> GetRHICfDetPointNum();
        if(outPointNum > 3){
            outPointNum = 0;
            continue;
        }
        for(int i=0; i<outPointNum; i++){
            mDetPoint = mEventDst -> GetRHICfDetPoint(i);

            outPointTowerIdx[i] = mDetPoint -> GetTowerIdx();
            outPointPID[i] = mDetPoint -> GetPID();
            outPointPos[i][0] = mDetPoint -> GetPointPos(0);
            outPointPos[i][1] = mDetPoint -> GetPointPos(1);
            outPointEnergy[i][0] = mDetPoint -> GetPointEnergy(0);
            outPointEnergy[i][1] = mDetPoint -> GetPointEnergy(1);
        }

        tree -> Fill();
    }

    file -> cd();
    tree -> Write();
    file -> Close();
}

void RHICfPicoEventMaker::PhysicsStream()
{
    int runType = GetOptContainer()->GetRunType();
    TString runName = GetOptContainer()->GetRunTypeName(runType);
    TFile* file = new TFile(Form("%s/RHICfPicoPhysics_%s.root", mOutputPath.Data(), runName.Data()), "recreate");
    TTree* tree = new TTree("event", "event");

    Int_t outRunNum;
    Int_t outFillNum;
    Int_t outRHICfRunType;

    Bool_t outRPSDTTrig;
    Bool_t outRPETTrig;
    Bool_t outRPCPTnoBBCLTrig;
    Bool_t outRPCPT2Trig;
    Bool_t outRPCPT2noBBCLTrig;
    Bool_t outVPDMB30Trig;
    Bool_t outVPDMB100Trig;
    Bool_t outBHT3Trig;
    Bool_t outTofHighMultTrig;
    Bool_t outEPDTrig;
    Bool_t outBBCTrig;
    Bool_t outBBCTacTrig;
    Bool_t outZDCTrig;
    Bool_t outZDCTacTrig;
    Bool_t outVPDTrig;
    Bool_t outZeroBiasTrig;

    Int_t outTOFMult;

    Int_t outBBCSmallADC[kBeamSideNum][kBBCSmallPmtNum];
    Int_t outBBCSmallTDC[kBeamSideNum][kBBCSmallPmtNum];
    Int_t outBBCSmallTAC[kBeamSideNum][kBBCSmallPmtNum];
    Int_t outBBCLargeADC[kBeamSideNum][kBBCLargePmtNum];
    Int_t outBBCLargeTDC[kBeamSideNum][kBBCLargePmtNum];
    Int_t outBBCLargeTAC[kBeamSideNum][kBBCLargePmtNum];

    Int_t outZDCADC[kBeamSideNum][kZDCPmtNum];
    Int_t outZDCTDC[kBeamSideNum][kZDCPmtNum];
    Int_t outZDCSumADC[kBeamSideNum][2]; // [Attenuated, Unattenuated]
    Int_t outZDCSMD[kBeamSideNum][2][kSMDNum];  //[east, west][vertical, horizontal][smd]

    tree -> Branch("RunNum", &outRunNum, "RunNum/I");
    tree -> Branch("FillNum", &outFillNum, "FillNum/I");
    tree -> Branch("RHICfRunType", &outRHICfRunType, "RHICfRunType/I");

    tree -> Branch("RPSDTTrig", &outRPSDTTrig, "RPSDTTrig/O");
    tree -> Branch("RPETTTrig", &outRPETTrig, "RPETTTrig/O");
    tree -> Branch("RPPTnoBBCLTTrig", &outRPCPTnoBBCLTrig, "RPPTnoBBCLTTrig/O");
    tree -> Branch("RPPT2TTrig", &outRPCPT2Trig, "RPPT2TTrig/O");
    tree -> Branch("RPSPT2noBBCLTrig", &outRPCPT2noBBCLTrig, "RPSPT2noBBCLTrig/O");
    tree -> Branch("VPDMB30Trig", &outVPDMB30Trig, "VPDMB30Trig/O");
    tree -> Branch("VPDMB100Trig", &outVPDMB100Trig, "VPDMB100Trig/O");
    tree -> Branch("BHT3Trig", &outBHT3Trig, "BHT3Trig/O");
    tree -> Branch("TofHighMultTrig", &outTofHighMultTrig, "TofHighMultTrig/O");
    tree -> Branch("EPDTrig", &outEPDTrig, "EPDTrig/O");
    tree -> Branch("BBCTrig", &outBBCTrig, "BBCTrig/O");
    tree -> Branch("BBCTacTrig", &outBBCTacTrig, "BBCTacTrig/O");
    tree -> Branch("ZDCTrig", &outZDCTrig, "ZDCTrig/O");
    tree -> Branch("ZDCTacTrig", &outZDCTacTrig, "ZDCTacTrig/O");
    tree -> Branch("VPDTrig", &outVPDTrig, "VPDTrig/O");
    tree -> Branch("ZeroBiasTrig", &outZeroBiasTrig, "ZeroBiasTrig/O");

    tree -> Branch("TOFMult", &outTOFMult, "TOFMult/I");
    tree -> Branch("BBCSmallADC", outBBCSmallADC, "BBCSmallADC[2][16]/I");
    tree -> Branch("BBCSmallTDC", outBBCSmallTDC, "BBCSmallTDC[2][16]/I");
    tree -> Branch("BBCSmallTAC", outBBCSmallTAC, "BBCSmallTAC[2][16]/I");
    tree -> Branch("BBCLargeADC", outBBCLargeADC, "BBCLargeADC[2][8]/I");
    tree -> Branch("BBCLargeTDC", outBBCLargeTDC, "BBCLargeTDC[2][8]/I");
    tree -> Branch("BBCLargeTAC", outBBCLargeTAC, "BBCLargeTAC[2][8]/I");

    tree -> Branch("ZDCADC", outZDCADC, "ZDCADC[2][3]/I");
    tree -> Branch("ZDCTDC", outZDCTDC, "ZDCTDC[2][3]/I");
    tree -> Branch("ZDCSumADC", outZDCSumADC, "ZDCSumADC[2][2]/I");
    tree -> Branch("ZDCSMD", outZDCSMD, "ZDCSMD[2][2][7]/I");

    int eventNum = mEventReader->GetEventNum();
    for(int event=0; event<eventNum; event++){
        if(event%10000 == 0){cout << "RHICfPicoEventMaker::ParticleTest() -- Event: " << event << " / " << eventNum << endl;}
        mEventDst = mEventReader -> GetEventDst(event);

        fill_n(&outBBCSmallADC[0][0], kBeamSideNum*kBBCSmallPmtNum, -1);
        fill_n(&outBBCSmallTDC[0][0], kBeamSideNum*kBBCSmallPmtNum, -1);
        fill_n(&outBBCSmallTAC[0][0], kBeamSideNum*kBBCSmallPmtNum, -1);
        fill_n(&outBBCLargeADC[0][0], kBeamSideNum*kBBCLargePmtNum, -1);
        fill_n(&outBBCLargeTDC[0][0], kBeamSideNum*kBBCLargePmtNum, -1);
        fill_n(&outBBCLargeTAC[0][0], kBeamSideNum*kBBCLargePmtNum, -1);
        fill_n(&outZDCADC[0][0], kBeamSideNum*kZDCPmtNum, -1);
        fill_n(&outZDCTDC[0][0], kBeamSideNum*kZDCPmtNum, -1);
        fill_n(&outZDCSumADC[0][0], kBeamSideNum*2, -1);
        fill_n(&outZDCSMD[0][0][0], kBeamSideNum*2*kSMDNum, -1);

        mEvent = mEventDst -> GetEvent();
        outFillNum = mEvent -> GetFillNumber();
        outRunNum = mEvent -> GetRunNumber();
        outRHICfRunType = mEvent -> GetRHICfRunType();

        outRPSDTTrig = mEvent -> GetRPSDTTrig();
        outRPETTrig = mEvent -> GetRPETTrig();
        outRPCPTnoBBCLTrig = mEvent -> GetRPCPTnoBBCLTrig();
        outRPCPT2Trig = mEvent -> GetRPCPT2Trig();
        outRPCPT2noBBCLTrig = mEvent -> GetRPCPT2noBBCLTrig();
        outVPDMB30Trig = mEvent -> GetVPDMB30Trig();
        outVPDMB100Trig = mEvent -> GetVPDMB100Trig();
        outBHT3Trig = mEvent -> GetBHT3Trig();
        outTofHighMultTrig = mEvent -> GetTofHighMultTrig();
        outEPDTrig = mEvent -> GetEPDTrig();
        outBBCTrig = mEvent -> GetBBCTrig();
        outBBCTacTrig = mEvent -> GetBBCTacTrig();
        outZDCTrig = mEvent -> GetZDCTrig();
        outZDCTacTrig = mEvent -> GetZDCTacTrig();
        outVPDTrig = mEvent -> GetVPDTrig();
        outZeroBiasTrig = mEvent -> GetZeroBiasTrig();

        outTOFMult = mEvent -> GetBTofMult();
        mBBC = mEventDst -> GetBBC();
        mZDC = mEventDst -> GetZDC();

        for(int ew=0; ew<kBeamSideNum; ew++){
            int beamSide = (ew==0)? kBeamEast : kBeamWest;

            for(int i=0; i<kBBCSmallPmtNum; i++){
                outBBCSmallADC[beamSide][i] = mBBC -> GetBBCSmallADC(beamSide, i);
                outBBCSmallTDC[beamSide][i] = mBBC -> GetBBCSmallTDC(beamSide, i);
                outBBCSmallTAC[beamSide][i] = mBBC -> GetBBCSmallTAC(beamSide, i);
            }
            for(int i=0; i<kBBCLargePmtNum; i++){
                outBBCLargeADC[beamSide][i] = mBBC -> GetBBCLargeADC(beamSide, i);
                outBBCLargeTDC[beamSide][i] = mBBC -> GetBBCLargeTDC(beamSide, i);
                outBBCLargeTAC[beamSide][i] = mBBC -> GetBBCLargeTAC(beamSide, i);
            }

            for(int i=0; i<kZDCPmtNum; i++){
                outZDCADC[beamSide][i] = mZDC -> GetZDCPmtADC(beamSide, i);
                outZDCTDC[beamSide][i] = mZDC -> GetZDCPmtTDC(beamSide, i);
            }
            outZDCSumADC[beamSide][0] = mZDC -> GetZDCPmtAttenuatedSumADC(beamSide);
            outZDCSumADC[beamSide][1] = mZDC -> GetZDCPmtUnAttenuatedSumADC(beamSide);
        }

        for(int ew=0; ew<kBeamSideNum; ew++){
            for(int ixy=0; ixy<2; ixy++){
                for(int smd=0; smd<7; smd++){
                    outZDCSMD[ew][ixy][smd] = mZDC -> GetZDCSmdADC(ew, ixy, smd);
                }
            }
        }

        tree -> Fill();
    }

    file -> cd();
    tree -> Write();
    file -> Close();
}

void RHICfPicoEventMaker::RPSEvent()
{
    const int rpBranchNum = 4;
    TString rpBranchName[rpBranchNum] = {"EU", "ED" , "WU", "WD"};
    TString ewName[2] = {"East", "West"};
    const double protonMass = 0.938; // [GeV/c^2]

    int colorIdx[4] = {1, 2, 4, 8};
    const int kFillNum = 5;
    const int kRunTypeNum = 3;
    const int kRunNum = 65;

    const TString fillLabel[kFillNum] = {"21142", "21145", "21148", "21149", "21150"};

    static const int runIdxArr[kRunNum] = {18175022, 18175023, 18175024, 18175025, 18175026, 
                                          18175027, 18175029, 18175030, 18176011, 18176012, 
                                          18176014, 18176016, 18176017, 18176018, 18176019, 
                                          18176020, 18176021, 18176033, 18176034, 18176035, 
                                          18176040, 18176042, 18176043, 18177001, 18177002, 
                                          18177003, 18177005, 18177011, 18177012, 18177014, 
                                          18177015, 18177016, 18177017, 18177018, 18177019, 
                                          18177020, 18177024, 18177025, 18177026, 18177027, 
                                          18177028, 18177029, 18177031, 18177032, 18177034, 
                                          18177036, 18177043, 18177046, 18177047, 18177049, 
                                          18177050, 18177052, 18178002, 18178003, 18178004, 
                                          18178005, 18178006, 18178007, 18178008, 18178009, 
                                          18178011, 18178012, 18178015, 18178016, 18178017
                                         };
                                      
    static const int NoRPSRunIdx[1] = {17};

    TLatex* globalLatex = new TLatex();

    // =========== histogram ==============
    TH1D* hRPSEvent[3][3]; // [runType, fill, run][rhicfEvent, rps point > 0, good trk event]
    for(int i=0; i<3; i++){
        hRPSEvent[0][i] = new TH1D(Form("hRPSEvent_runType_%i", i), "", kRunTypeNum, 0, kRunTypeNum);
        hRPSEvent[1][i] = new TH1D(Form("hRPSEvent_Fill_%i", i), "", kFillNum, 0, kFillNum);
        hRPSEvent[2][i] = new TH1D(Form("hRPSEvent_Run_%i", i), "", kRunNum, 0, kRunNum);
        hRPSEvent[0][i] -> SetLineColor(colorIdx[i]);
        hRPSEvent[1][i] -> SetLineColor(colorIdx[i]);
        hRPSEvent[2][i] -> SetLineColor(colorIdx[i]);
        hRPSEvent[0][i] -> SetStats(0);
        hRPSEvent[1][i] -> SetStats(0);
        hRPSEvent[2][i] -> SetStats(0);
    }

    TH1D* hPlane_type[3]; // [all, local, global]
    for(int i=0; i<3; i++){
        hPlane_type[i] = new TH1D(Form("plane_type%i", i), "", 10, 0, 10);
        hPlane_type[i] -> SetLineColor(colorIdx[i]);
        hPlane_type[i] -> SetStats(0);
    }

    TH1D* hEnergyQA[rpBranchNum][3]; // [all, usedplane cut, type cut]
    for(int i=0; i<rpBranchNum; i++){
        hEnergyQA[i][0] = new TH1D(Form("energy_%s_%i", rpBranchName[i].Data(), 0), "", 60., 0., 320.);
        hEnergyQA[i][1] = new TH1D(Form("energy_%s_%i", rpBranchName[i].Data(), 1), "", 60., 0., 320.);
        hEnergyQA[i][2] = new TH1D(Form("energy_%s_%i", rpBranchName[i].Data(), 2), "", 60., 0., 320.);
        hEnergyQA[i][0] -> SetLineColor(colorIdx[0]);
        hEnergyQA[i][1] -> SetLineColor(colorIdx[1]);
        hEnergyQA[i][2] -> SetLineColor(colorIdx[2]);
        hEnergyQA[i][0] -> SetStats(0);
        hEnergyQA[i][1] -> SetStats(0);
        hEnergyQA[i][2] -> SetStats(0);
    }

    const int ewSide = 2;
    TH1D* hThetaX[ewSide][2];  // [ew][all , pi0]
    TH1D* hThetaY[ewSide][2];  // [ew][all , pi0]
    TH2D* hThetaXY[ewSide][2]; // [ew][all , pi0]
    TH2D* hPxy[ewSide][2];  // [ew][all , pi0]
    TH1D* hXi[ewSide][2];  // [ew][all , pi0] 
    TH1D* hEnergy[ewSide][2];  // [ew][all , pi0]
    TH1D* hSumE[ewSide][2];  // [ew][all, Fiducial cut]

    TH2D* hRPS_Pi0[ewSide][2];  // [ew][all, Fiducial cut]

    for(int i=0; i<ewSide; i++){
        for(int j=0; j<2; j++){
            hThetaX[i][j] = new TH1D(Form("hThetaX%s_%i", ewName[i].Data(), j), "", 80, -6.5, 6.5);
            hThetaY[i][j] = new TH1D(Form("hThetaY%s_%i", ewName[i].Data(), j), "", 80, -6.5, 6.5);
            hThetaX[i][j] -> SetLineColor(colorIdx[j]);
            hThetaY[i][j] -> SetLineColor(colorIdx[j]);
            hThetaX[i][j] -> SetStats(0);
            hThetaY[i][j] -> SetStats(0);

            hThetaXY[i][j] = new TH2D(Form("hThetaXY%s_%i", ewName[i].Data(), j), "", 80, -6.5, 6.5, 80, -6.5, 6.5);
            hThetaXY[i][j] -> SetStats(0);

            hPxy[i][j] = new TH2D(Form("hPxy%s", ewName[i].Data()), "", 80, -1., 1., 80, -1., 1.);
            hPxy[i][j] -> SetStats(0);

            hXi[i][j] = new TH1D(Form("hXi%s_%i", ewName[i].Data(), j), "", 80, 0., 1.);
            hXi[i][j] -> SetLineColor(colorIdx[j]);
            hXi[i][j] -> SetStats(0);

            hEnergy[i][j] = new TH1D(Form("hEnergy%s_%i", ewName[i].Data(), j), "", 80., 50., 450.);
            hEnergy[i][j] -> SetLineColor(colorIdx[j]);
            hEnergy[i][j] -> SetStats(0);

            hSumE[i][j] = new TH1D(Form("hEnergy%s_%i", ewName[i].Data(), j), "", 80., 50., 450.);
            hSumE[i][j] -> SetLineColor(colorIdx[j]);
            hSumE[i][j] -> SetStats(0);
        }
    }


    TFile* file = new TFile("../RHICfPi0RPS.root","recreate");
    TTree* tree = new TTree("event", "event");
    
    int outRunIdx;
    int outFillIdx;
    int outRHICfRunIdx;
    bool outIsRpsSDTrig;
    bool outIsRpsETrig;
    bool outIsRHICfDiffTrig;
    int outTOFMult;
    int outBBCSumADC[2][2];
    int outZDCADC[2][3]; 

    int outRPSNum;
    int outRPSBranch[5];
    double outRPSXi[5];
    double outRPSTheta[5][2];
    double outRPSTrkPos[5][2][3];
    double outRPSP[5][4];

    bool outIsPi0;
    double outPi0E;
    double outPi0Mass;

    tree -> Branch("runIdx", &outRunIdx, "runIdx/I");
    tree -> Branch("fillIdx", &outFillIdx, "fillIdx/I");
    tree -> Branch("RHICfRunIdx", &outRHICfRunIdx, "RHICfRunIdx/I");
    tree -> Branch("isRpsSDTrig", &outIsRpsSDTrig, "isRpsSDTrig/O");
    tree -> Branch("isRpsETrig", &outIsRpsETrig, "isRpsETrig/O");
    tree -> Branch("isRHICfDiffTrig", &outIsRHICfDiffTrig, "isRHICfDiffTrig/O");

    tree -> Branch("tofMult", &outTOFMult, "tofMult/I");
    tree -> Branch("bbcSumADC", outBBCSumADC, "bbcSumADC[2][2]/I");
    tree -> Branch("zdcADC", outZDCADC, "zdcADC[2][3]/I");

    tree -> Branch("rpsNum", &outRPSNum, "rpsNum/I");
    tree -> Branch("rpsBranch", outRPSBranch, "rpsBranch[5]/I");
    tree -> Branch("rpsXi", outRPSXi, "rpsXi[5]/D");
    tree -> Branch("rpsTheta", outRPSTheta, "rpsTheta[5][2]/D");
    tree -> Branch("rpsTrkPos", outRPSTrkPos, "rpsTrkPos[5][2][3]/D");
    tree -> Branch("rpsP", outRPSP, "rpsP[5][4]/D");
    tree -> Branch("isPi0", &outIsPi0, "isPi0/O");
    tree -> Branch("pi0E", &outPi0E, "pi0E/D");
    tree -> Branch("pi0Mass", &outPi0Mass, "pi0Mass/D");

    int eventNum = mEventReader->GetEventNum();
    for(int event=0; event<eventNum; event++){
        if(event%10000 == 0){cout << "RHICfPicoEventMaker::ParticleTest() -- Event: " << event << " / " << eventNum << endl;}
        mEventDst = mEventReader -> GetEventDst(event);

        outRPSNum = 0;
        memset(outRPSBranch, 0., sizeof(outRPSBranch));
        memset(outRPSXi, 0., sizeof(outRPSXi));
        memset(outRPSTheta, 0., sizeof(outRPSTheta));
        memset(outRPSTrkPos, 0., sizeof(outRPSTrkPos));
        memset(outRPSP, 0., sizeof(outRPSP));
        outIsPi0 = false;
        outPi0E = 0.;
        outPi0Mass = 0.;

        mEvent = mEventDst -> GetEvent();
        int fillNum = mEvent -> GetFillNumber();
        int runNum = mEvent -> GetRunNumber();
        int runType = mEvent -> GetRHICfRunType();

        outRunIdx = -1;
        for(int i=0; i<runNum; i++){
            if(runIdxArr[i] == runNum){
                outRunIdx = i;
                break;
            }
        }

        outTOFMult = mEvent -> GetBTofMult();
        mBBC = mEventDst -> GetBBC();

        outBBCSumADC[0][0] = mBBC -> GetBBCSmallSumADC(0);
        outBBCSumADC[0][1] = mBBC -> GetBBCLargeSumADC(0);
        outBBCSumADC[1][0] = mBBC -> GetBBCSmallSumADC(1);
        outBBCSumADC[1][1] = mBBC -> GetBBCLargeSumADC(1);

        mZDC = mEventDst -> GetZDC();
        outZDCADC[0][0] = mZDC -> GetZDCPmtSumADC(0);
        outZDCADC[1][0] = mZDC -> GetZDCPmtSumADC(1);
        outZDCADC[0][1] = mZDC -> GetZDCPmtAttenuatedSumADC(0);
        outZDCADC[1][1] = mZDC -> GetZDCPmtAttenuatedSumADC(1);
        outZDCADC[0][2] = mZDC -> GetZDCPmtUnAttenuatedSumADC(0);
        outZDCADC[1][2] = mZDC -> GetZDCPmtUnAttenuatedSumADC(1);

        outFillIdx = RHICfOptContainer::GetFillNumIdx(fillNum);
        outRHICfRunIdx = RHICfOptContainer::GetFillNumToRunIdx(fillNum);

        outIsRHICfDiffTrig = mEvent -> GetRHICfDiffractiveTrig();
        outIsRpsSDTrig = mEvent -> GetRPSDTTrig();
        outIsRpsETrig = mEvent -> GetRPETTrig();

        double beamEnergy = mEvent -> GetBeamEnergy();
        double beamP = sqrt(beamEnergy*beamEnergy - protonMass*protonMass);

        // ========================== RHICf Pi0 ============================
        bool isPi0Event = false;
        double pi0Energy = 0.;
        double pi0Mass = 0.;
        mParticleMaker -> SetEventDst(mEventDst);
        MiniParticle particles = mParticleMaker -> GetMiniParticle();
        int particleNum = particles.particleNum;
        for(int i=0; i<particleNum; i++){
            double m = -999.;
            if(particles.m.size() != 0){m = particles.m[i];}
            if(135.-3.*9. <= m && m <= 135.+3.*9.){
                isPi0Event = true;
                pi0Energy = particles.e[i];
                pi0Mass = m;
            }
        }

        if(isPi0Event){
            outIsPi0 = true;
            outPi0E = pi0Energy;
            outPi0Mass = pi0Mass;
        }

        // ========================= Roman pots ============================
        int rpsNum = mEventDst -> GetRPSNum();
        bool isGoodRPSEvent = false;
        for(int i=0; i<rpsNum; i++){
            mRPS = mEventDst -> GetRPS(i);

            int trkBranch = mRPS -> GetTrkBranch(); // EU=0, ED=1, WU=2, WD=3 
            int plane = mRPS -> GetUsedPlane();     // used plane
            int trkType = mRPS -> GetTrkType();     // 0 == reco. only one point (local) , 1 == reco. using two point (global) , 2 == undefined

            // position
            double hitPosX1  = mRPS -> GetRPSHitPosX(0);
            double hitPosY1  = mRPS -> GetRPSHitPosY(0);
            double hitPosZ1  = mRPS -> GetRPSHitPosZ(0);
            double hitPosX2  = mRPS -> GetRPSHitPosX(1);
            double hitPosY2  = mRPS -> GetRPSHitPosY(1);
            double hitPosZ2  = mRPS -> GetRPSHitPosZ(1);

            // momentum
            double px = mRPS -> GetPx(); // [GeV/c]
            double py = mRPS -> GetPy(); // [GeV/c]
            double pz = mRPS -> GetPz(); // [GeV/c]
            double pt = sqrt(px*px + py*py); // [GeV/c]
            double p = sqrt(px*px + py*py + pz*pz); // [GeV/c]
            double e = sqrt(p*p + protonMass*protonMass) - protonMass; // kinetic energy [GeV]

            double xi = (beamP - p)/beamP;

            // angle
            double thetaX = atan(px/fabs(pz)) * 1000.; // [mrad]
            double thetaY = atan(py/fabs(pz)) * 1000.; // [mrad]
            double thetaT = atan(pt/fabs(pz)) * 1000.; // [mrad] 

            // ===================== RPS Kinematics Fill ========================
            if(plane >= 5 && trkType == 1){
                isGoodRPSEvent = true;
                outRPSNum = rpsNum;
                outRPSBranch[i] = trkBranch;
                outRPSXi[i] = xi;
                outRPSTheta[i][0] = thetaX;
                outRPSTheta[i][1] = thetaY;
                outRPSP[i][0] = px;
                outRPSP[i][1] = py;
                outRPSP[i][2] = pz;
                outRPSP[i][3] = e;

                outRPSTrkPos[i][0][0] = hitPosX1;
                outRPSTrkPos[i][0][1] = hitPosY1;
                outRPSTrkPos[i][0][2] = hitPosZ1;
                outRPSTrkPos[i][1][0] = hitPosX2;
                outRPSTrkPos[i][1][1] = hitPosY2;
                outRPSTrkPos[i][1][2] = hitPosZ2;
            }
        }

        if(outRPSNum != 0 || isPi0Event){
            tree -> Fill();
        }
    }

    file -> cd();
    tree -> Write();
    file -> Close();
}

void RHICfPicoEventMaker::RHICfParticleEvent()
{
    int runType = GetOptContainer() -> GetRunType();
    TString runName = GetOptContainer() -> GetRunTypeName(runType);
    TString particleName = GetOptContainer() -> GetParticleRunName();
    TFile* file = new TFile(Form("%s/RHICfPicoEvent_%s_%s.root", mOutputPath.Data(), particleName.Data(), runName.Data()), "recreate");
    TTree* tree = new TTree("event", "event");

    Int_t outRunNum;
    Int_t outFillNum;
    Int_t outRHICfRunType;

    // STAR trigger
    Bool_t outRHICfTrig;
    Bool_t outRHICfDiffractiveTrig;
    Bool_t outRHICfVPDMC30Trig;
    Bool_t outRHICfTPCTrig;
    Bool_t outRPSDTTrig;
    Bool_t outRPETTrig;
    Bool_t outRPCPTnoBBCLTrig;
    Bool_t outRPCPT2Trig;
    Bool_t outRPCPT2noBBCLTrig;
    Bool_t outVPDMB30Trig;
    Bool_t outVPDMB100Trig;
    Bool_t outBHT3Trig;
    Bool_t outTofHighMultTrig;
    Bool_t outEPDTrig;
    Bool_t outBBCTrig;
    Bool_t outBBCTacTrig;
    Bool_t outZDCTrig;
    Bool_t outZDCTacTrig;
    Bool_t outVPDTrig;
    Bool_t outZeroBiasTrig;

    Int_t outDLEIdx;
    Int_t outBTofMult;
    Int_t outBBCSumADC[kBeamSideNum][2];
    
    Int_t RHICfParticleNum;
    Int_t outPi0Type[4];
    Int_t outTowerIdx[4];
    Double_t outPos[4][2];
    Double_t outMomentum[4][4];
    Double_t outMass[4];

    Double_t outRHICfPlatedE[2][16];

    Int_t outZDCADC[kBeamSideNum][kZDCPmtNum];
    Int_t outZDCTDC[kBeamSideNum][kZDCPmtNum];
    Int_t outZDCSumADC[kBeamSideNum][2]; // [Attenuated, Unattenuated]
    Int_t outZDCSMD[kBeamSideNum][2][kSMDNum];  //[east, west][vertical, horizontal][smd]

    tree -> Branch("RunNum", &outRunNum, "RunNum/I");
    tree -> Branch("FillNum", &outFillNum, "FillNum/I");
    tree -> Branch("RHICfRunType", &outRHICfRunType, "RHICfRunType/I");

    tree -> Branch("RHICfTrig", &outRHICfTrig, "RHICfTrig/O");
    tree -> Branch("RHICfDiffractiveTrig", &outRHICfDiffractiveTrig, "RHICfDiffractiveTrig/O");
    tree -> Branch("RHICfVPDMC30Trig", &outRHICfVPDMC30Trig, "RHICfVPDMC30Trig/O");
    tree -> Branch("RHICfTPCTrig", &outRHICfTPCTrig, "RHICfTPCTrig/O");
    tree -> Branch("RPSDTTrig", &outRPSDTTrig, "RPSDTTrig/O");
    tree -> Branch("RPETTTrig", &outRPETTrig, "RPETTTrig/O");
    tree -> Branch("RPPTnoBBCLTTrig", &outRPCPTnoBBCLTrig, "RPPTnoBBCLTTrig/O");
    tree -> Branch("RPPT2TTrig", &outRPCPT2Trig, "RPPT2TTrig/O");
    tree -> Branch("RPSPT2noBBCLTrig", &outRPCPT2noBBCLTrig, "RPSPT2noBBCLTrig/O");
    tree -> Branch("VPDMB30Trig", &outVPDMB30Trig, "VPDMB30Trig/O");
    tree -> Branch("VPDMB100Trig", &outVPDMB100Trig, "VPDMB100Trig/O");
    tree -> Branch("BHT3Trig", &outBHT3Trig, "BHT3Trig/O");
    tree -> Branch("TofHighMultTrig", &outTofHighMultTrig, "TofHighMultTrig/O");
    tree -> Branch("EPDTrig", &outEPDTrig, "EPDTrig/O");
    tree -> Branch("BBCTrig", &outBBCTrig, "BBCTrig/O");
    tree -> Branch("BBCTacTrig", &outBBCTacTrig, "BBCTacTrig/O");
    tree -> Branch("ZDCTrig", &outZDCTrig, "ZDCTrig/O");
    tree -> Branch("ZDCTacTrig", &outZDCTacTrig, "ZDCTacTrig/O");
    tree -> Branch("VPDTrig", &outVPDTrig, "VPDTrig/O");
    tree -> Branch("ZeroBiasTrig", &outZeroBiasTrig, "ZeroBiasTrig/O");

    tree -> Branch("DLEIdx", &outDLEIdx, "DLEIdx/I");
    tree -> Branch("BTofMult", &outBTofMult, "BTofMult/I");
    tree -> Branch("BBCSumADC", outBBCSumADC, "BBCSumADC[2][2]/I");

    tree -> Branch("particleNum", &RHICfParticleNum, "RHICfParticleNum/I");
    tree -> Branch("Pi0Type", outPi0Type, "Pi0Type[4]/I");
    tree -> Branch("TowerIdx", outTowerIdx, "TowerIdx[4]/I");
    tree -> Branch("Position", outPos, "Position[4][2]/D");
    tree -> Branch("Momentum", outMomentum, "Momentum[4][4]/D");
    tree -> Branch("Mass", outMass, "Mass[4]/D");

    tree -> Branch("RHICfPlatedE", outRHICfPlatedE, "outRHICfPlatedE[2][16]/D");
    
    tree -> Branch("ZDCADC", outZDCADC, "ZDCADC[2][3]/I");
    tree -> Branch("ZDCTDC", outZDCTDC, "ZDCTDC[2][3]/I");
    tree -> Branch("ZDCSumADC", outZDCSumADC, "ZDCSumADC[2][2]/I");
    tree -> Branch("ZDCSMD", outZDCSMD, "ZDCSMD[2][2][8]/I");

    int eventNum = mEventReader->GetEventNum();
    for(int event=0; event<eventNum; event++){
        if(event%10000 == 0){cout << "RHICfPicoEventMaker::RHICfParticleEvent() -- Event: " << event << " / " << eventNum << endl;}
        mEventDst = mEventReader -> GetEventDst(event);
        mParticleMaker -> SetEventDst(mEventDst);

        // RHICf particles
        MiniParticle particles = mParticleMaker -> GetMiniParticle();
        RHICfParticleNum = particles.particleNum;
        if(RHICfParticleNum == 0 || RHICfParticleNum > 3){continue;}
        if(GetOptContainer() -> GetParticleRunIdx() == kNeutron && RHICfParticleNum != 1){continue;}
        
        mEvent = mEventDst -> GetEvent();

        outFillNum = mEvent -> GetFillNumber();
        outRunNum = mEvent -> GetRunNumber();
        outRHICfRunType = mEvent -> GetRHICfRunType();

        outRHICfTrig = mEvent -> GetRHICfTrig();
        outRHICfDiffractiveTrig = mEvent -> GetRHICfDiffractiveTrig();
        outRHICfVPDMC30Trig = mEvent -> GetRHICfVPDMB30Trig();
        outRHICfTPCTrig = mEvent -> GetRHICfTPCTrig();
        outRPSDTTrig = mEvent -> GetRPSDTTrig();
        outRPETTrig = mEvent -> GetRPETTrig();
        outRPCPTnoBBCLTrig = mEvent -> GetRPCPTnoBBCLTrig();
        outRPCPT2Trig = mEvent -> GetRPCPT2Trig();
        outRPCPT2noBBCLTrig = mEvent -> GetRPCPT2noBBCLTrig();
        outVPDMB30Trig = mEvent -> GetVPDMB30Trig();
        outVPDMB100Trig = mEvent -> GetVPDMB100Trig();
        outBHT3Trig = mEvent -> GetBHT3Trig();
        outTofHighMultTrig = mEvent -> GetTofHighMultTrig();
        outEPDTrig = mEvent -> GetEPDTrig();
        outBBCTrig = mEvent -> GetBBCTrig();
        outBBCTacTrig = mEvent -> GetBBCTacTrig();
        outZDCTrig = mEvent -> GetZDCTrig();
        outZDCTacTrig = mEvent -> GetZDCTacTrig();
        outVPDTrig = mEvent -> GetVPDTrig();
        outZeroBiasTrig = mEvent -> GetZeroBiasTrig();

        outBTofMult = mEvent -> GetBTofMult();
        mBBC = mEventDst -> GetBBC();
        mZDC = mEventDst -> GetZDC();

        outBBCSumADC[0][0] = mBBC -> GetBBCSmallSumADC(0);
        outBBCSumADC[0][1] = mBBC -> GetBBCLargeSumADC(0);
        outBBCSumADC[1][0] = mBBC -> GetBBCSmallSumADC(1);
        outBBCSumADC[1][1] = mBBC -> GetBBCLargeSumADC(1);

        mDLECondition -> SetBTofMult(outBTofMult);
        mDLECondition -> SetBBCADC(outBBCSumADC[0][0], outBBCSumADC[1][0], outBBCSumADC[0][1], outBBCSumADC[1][1]);
        outDLEIdx = mDLECondition -> GetDLEIdx();

        for(int ew=0; ew<kBeamSideNum; ew++){
            int beamSide = (ew==0)? kBeamEast : kBeamWest;
            for(int i=0; i<kZDCPmtNum; i++){
                outZDCADC[beamSide][i] = mZDC -> GetZDCPmtADC(beamSide, i);
                outZDCTDC[beamSide][i] = mZDC -> GetZDCPmtTDC(beamSide, i);
            }
            outZDCSumADC[beamSide][0] = mZDC -> GetZDCPmtAttenuatedSumADC(beamSide);
            outZDCSumADC[beamSide][1] = mZDC -> GetZDCPmtUnAttenuatedSumADC(beamSide);
        }

        for(int ew=0; ew<kBeamSideNum; ew++){
            for(int ixy=0; ixy<2; ixy++){
                for(int smd=0; smd<kSMDNum; smd++){
                    outZDCSMD[ew][ixy][smd] = mZDC -> GetZDCSmdADC(ew, ixy, smd);
                }
            }
        }

        for(int i=0; i<RHICfParticleNum; i++){
            outPi0Type[i] = particles.type[i]-1;
            outTowerIdx[i] = particles.towerIdx[i];
            outPos[i][0] = particles.x[i];
            outPos[i][1] = particles.y[i];
            outMomentum[i][0] = particles.px[i];
            outMomentum[i][1] = particles.py[i];
            outMomentum[i][2] = particles.pz[i];
            outMomentum[i][3] = particles.e[i];
            outMass[i] = particles.m[i];
        }

        mDetHit = mEventDst -> GetRHICfDetHit();
        for(int it=0; it<kTowerNum; it++){
            for(int ip=0; ip<kPlateNum; ip++){
                outRHICfPlatedE[it][ip] = mDetHit -> GetPlateEnergy(it, ip);
            }
        }


        tree -> Fill();
    }

    file -> cd();
    tree -> Write();
    file -> Close();
}