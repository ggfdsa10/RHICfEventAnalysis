#include "RHICfAnAnalysis.hh"

RHICfAnAnalysis::RHICfAnAnalysis() 
{
}

RHICfAnAnalysis::~RHICfAnAnalysis()
{
}

int RHICfAnAnalysis::Init()
{
    RHICfOptContainer::SetOffTPCTrack();
    RHICfOptContainer::SetOffBTof();
    RHICfOptContainer::SetOffVPD();
    RHICfOptContainer::SetOffZDC();
    RHICfOptContainer::SetOffFMS();
    if(RHICfOptContainer::GetParticleRunIdx()==kGammaRun){RHICfOptContainer::SetOffParticle();}
    if(RHICfOptContainer::GetParticleRunIdx()==kPi0Run){RHICfOptContainer::SetOffDetPoint();}
    if(RHICfOptContainer::GetParticleRunIdx()==kNeutronRun){RHICfOptContainer::SetOffDetPoint();}
    if(RHICfOptContainer::GetParticleRunIdx()==kLambda0Run){RHICfOptContainer::SetOffParticle();}

    if(!RHICfOptContainer::Init()){return 0;}
    if(!RHICfTableMaker::Init()){return 0;}

    // Utilized classes
    mEventReader = new RHICfEventDstReader();
    mParticleMaker = new RHICfParticleMaker();

    // Calculation classes
    mFigureDrawing = new RHICfFigureDrawing();
    mDLECondition = new RHICfDLECondition();
    mMassFitting = new RHICfMassFitting();
    mBinning = new RHICfBinning();
    mDilution = new RHICfDilutionFactor();
    mAsymmetry = new RHICfAsymmetry();

    mFigureDrawing -> SetMassFitting(mMassFitting);
    mFigureDrawing -> SetBinning(mBinning);
    mFigureDrawing -> SetDilution(mDilution);
    mMassFitting -> SetBinning(mBinning);
    mDilution -> SetBinning(mBinning);
    mAsymmetry -> SetBinning(mBinning);
    mAsymmetry -> SetDilution(mDilution);

    mEventReader -> Init();
    mParticleMaker -> Init();   
    mFigureDrawing -> Init();
    mDLECondition -> Init();
    mMassFitting -> Init();
    mBinning -> Init();
    mDilution -> Init();
    mAsymmetry -> Init();

    return 1;
}

int RHICfAnAnalysis::Calculate()
{
    MassCalculation();
    BinningCalculation();
    if(mMassFitting->GetMassTableFlag() != kExistTable){Calculate();}

    AnCalculation();

    return 1;
}

int RHICfAnAnalysis::Finish()
{
    return 1;
}

int RHICfAnAnalysis::MassCalculation()
{
    int particleRunIdx = GetOptContainer() -> GetParticleRunIdx();
    if(particleRunIdx == kNeutronRun || particleRunIdx == kGammaRun){return 1;}
    if(mMassFitting->GetMassTableFlag() == kExistTable){return 1;}

    mMassFitting -> InitHist();

    int eventNum = mEventReader->GetEventNum();
    for(int event=0; event<eventNum; event++){
        if(event%10000 == 0){cout << "RHICfAnAnalysis::MassCalculation() -- Event: " << event << " / " << eventNum << endl;}
        mEventDst = mEventReader -> GetEventDst(event);
        mEvent = mEventDst -> GetEvent();
        mBBC = mEventDst -> GetBBC();

        mParticleMaker -> SetEventDst(mEventDst);
        int runType = mEvent -> GetRHICfRunType();

        // STAR BTOF and BBC detectors 
        int btofMult = mEvent -> GetBTofMult();
        int bbcSmallEast = mBBC -> GetBBCSmallSumADC(kBeamEast);
        int bbcSmallWest = mBBC -> GetBBCSmallSumADC(kBeamWest);
        int bbcLargeEast = mBBC -> GetBBCLargeSumADC(kBeamEast);
        int bbcLargeWest = mBBC -> GetBBCLargeSumADC(kBeamWest);

        // DLE condition
        mDLECondition -> SetBTofMult(btofMult);
        mDLECondition -> SetBBCADC(bbcSmallEast, bbcSmallWest, bbcLargeEast, bbcLargeWest);
        int DLEIdx = mDLECondition -> GetDLEIdx();

        // RHICf particles
        MiniParticle particles = mParticleMaker -> GetMiniParticle();
        int particleNum = particles.particleNum;
        for(int i=0; i<particleNum; i++){
            int typeIdx = particles.type[i]-1;
            double m = particles.m[i];
            double particlePt = particles.pt[i];
            double particleXf = particles.xf[i];

            if(DLEIdx < kALLDLE){
                mMassFitting -> FillMass(runType, typeIdx, DLEIdx, m);
            }
            mMassFitting -> FillMass(runType, typeIdx, kALLDLE, m);

            int ptNum = mBinning -> GetPtBinNum(runType, typeIdx, DLEIdx);
            int xfNum = mBinning -> GetXfBinNum(runType, typeIdx, DLEIdx);
            for(int pt=0; pt<ptNum; pt++){
                double ptLowerBoundary = mBinning -> GetPtBinBoundary(runType, typeIdx, DLEIdx, pt);
                double ptUpperBoundary = mBinning -> GetPtBinBoundary(runType, typeIdx, DLEIdx, pt+1);

                for(int xf=0; xf<xfNum; xf++){
                    double xfLowerBoundary = mBinning -> GetXfBinBoundary(runType, typeIdx, DLEIdx, xf);
                    double xfUpperBoundary = mBinning -> GetXfBinBoundary(runType, typeIdx, DLEIdx, xf+1);
                    
                    if(ptLowerBoundary <= particlePt && particlePt < ptUpperBoundary){
                        if(xfLowerBoundary <= particleXf && particleXf < xfUpperBoundary){
                            if(DLEIdx < kALLDLE){
                                mMassFitting -> FillMass(runType, typeIdx, DLEIdx, pt, xf, m);
                            }
                        }
                    }
                }
            }

            // All DLE conditions
            int ptNumAllDLE = mBinning -> GetPtBinNum(runType, typeIdx, kALLDLE);
            int xfNumAllDLE = mBinning -> GetXfBinNum(runType, typeIdx, kALLDLE);
            for(int pt=0; pt<ptNumAllDLE; pt++){
                double ptLowerBoundary = mBinning -> GetPtBinBoundary(runType, typeIdx, kALLDLE, pt);
                double ptUpperBoundary = mBinning -> GetPtBinBoundary(runType, typeIdx, kALLDLE, pt+1);

                for(int xf=0; xf<xfNumAllDLE; xf++){
                    double xfLowerBoundary = mBinning -> GetXfBinBoundary(runType, typeIdx, kALLDLE, xf);
                    double xfUpperBoundary = mBinning -> GetXfBinBoundary(runType, typeIdx, kALLDLE, xf+1);
                    
                    if(ptLowerBoundary <= particlePt && particlePt < ptUpperBoundary){
                        if(xfLowerBoundary <= particleXf && particleXf < xfUpperBoundary){
                            mMassFitting -> FillMass(runType, typeIdx, kALLDLE, pt, xf, m);
                        }
                    }
                }
            }
        }
    }

    mMassFitting -> Fitting();
    mMassFitting -> SaveMassData();
    mFigureDrawing -> DrawMassHist();

    return 1;
}

int RHICfAnAnalysis::BinningCalculation()
{
    if(mBinning->GetBinningTableFlag() == kExistTable){return 1;}

    int particleRunIdx = GetOptContainer() -> GetParticleRunIdx();
    int eventNum = mEventReader->GetEventNum();
    for(int event=0; event<eventNum; event++){
        if(event%10000 == 0){cout << "RHICfAnAnalysis::BinningCalculation() -- Event: " << event << " / " << eventNum << endl;}
        mEventDst = mEventReader -> GetEventDst(event);
        mEvent = mEventDst -> GetEvent();
        mBBC = mEventDst -> GetBBC();

        mParticleMaker -> SetEventDst(mEventDst);
        int runType = mEvent -> GetRHICfRunType();

        // STAR BTOF and BBC detectors 
        int btofMult = mEvent -> GetBTofMult();
        int bbcSmallEast = mBBC -> GetBBCSmallSumADC(kBeamEast);
        int bbcSmallWest = mBBC -> GetBBCSmallSumADC(kBeamWest);
        int bbcLargeEast = mBBC -> GetBBCLargeSumADC(kBeamEast);
        int bbcLargeWest = mBBC -> GetBBCLargeSumADC(kBeamWest);

        // DLE condition
        mDLECondition -> SetBTofMult(btofMult);
        mDLECondition -> SetBBCADC(bbcSmallEast, bbcSmallWest, bbcLargeEast, bbcLargeWest);
        int dleIdx = mDLECondition -> GetDLEIdx();

        // RHICf particles
        MiniParticle particles = mParticleMaker -> GetMiniParticle();
        int particleNum = particles.particleNum;
        for(int i=0; i<particleNum; i++){
            int typeIdx = particles.type[i]-1;
            double pt = particles.pt[i];
            double xf = particles.xf[i];
            double m = -999.;
            if(particles.m.size() != 0){m = particles.m[i];}

            double massLowerBoundary = mMassFitting -> GetMassLowerBoundary(runType, typeIdx, dleIdx);
            double massUpperBoundary = mMassFitting -> GetMassUpperBoundary(runType, typeIdx, dleIdx);
            if(massLowerBoundary <= m && m < massUpperBoundary || particleRunIdx == kNeutron || particleRunIdx == kGammaRun){
                if(dleIdx < kALLDLE){
                    mBinning -> FillKinematics(runType, typeIdx, dleIdx, pt, xf);
                }
                mBinning -> FillKinematics(runType, typeIdx, kALLDLE,  pt, xf);
            }
        }
    }

    mBinning -> Binning();
    mBinning -> SaveBinningData();
    mFigureDrawing -> DrawBinningHist();

    return 1;
}

int RHICfAnAnalysis::AnCalculation()
{
    mDilution -> InitHist();
    mAsymmetry -> InitData();

    int particleRunIdx = GetOptContainer() -> GetParticleRunIdx();
    int eventNum = mEventReader->GetEventNum();
    for(int event=0; event<eventNum; event++){
        if(event%10000 == 0){cout << "RHICfAnAnalysis::AnCalculation() -- Event: " << event << " / " << eventNum << endl;}
        mEventDst = mEventReader -> GetEventDst(event);
        mEvent = mEventDst -> GetEvent();
        mBBC = mEventDst -> GetBBC();

        mParticleMaker -> SetEventDst(mEventDst);
        int runType = mEvent -> GetRHICfRunType();
        int fillNum = mEvent -> GetFillNumber();
        int fillIdx = GetOptContainer() -> GetFillNumIdx(fillNum);
        int spinBit = mEvent -> GetSpinBit();
        if(spinBit > 10 || spinBit <= 0){continue;}

        bool blueSpinUp = false;
        bool yellowSpinUp = false;
        if(spinBit == 5 || spinBit == 6){blueSpinUp = true;}
        if(spinBit == 5 || spinBit == 9){yellowSpinUp = true;}

        // STAR BTOF and BBC detectors 
        int btofMult = mEvent -> GetBTofMult();
        int bbcSmallEast = mBBC -> GetBBCSmallSumADC(kBeamEast);
        int bbcSmallWest = mBBC -> GetBBCSmallSumADC(kBeamWest);
        int bbcLargeEast = mBBC -> GetBBCLargeSumADC(kBeamEast);
        int bbcLargeWest = mBBC -> GetBBCLargeSumADC(kBeamWest);

        // DLE condition
        mDLECondition -> SetBTofMult(btofMult);
        mDLECondition -> SetBBCADC(bbcSmallEast, bbcSmallWest, bbcLargeEast, bbcLargeWest);
        int dleIdx = mDLECondition -> GetDLEIdx();

        // RHICf particles
        MiniParticle particles = mParticleMaker -> GetMiniParticle();
        int particleNum = particles.particleNum;
        for(int i=0; i<particleNum; i++){
            int typeIdx = particles.type[i]-1;
            double x = particles.x[i];
            double y = particles.y[i];
            double particlePt = particles.pt[i];
            double particleXf = particles.xf[i];
            double phiAngle = (180./TMath::Pi())*TMath::ATan2(y, x);
            double m = -999.;
            if(particles.m.size() != 0){m = particles.m[i];}

            // pT and xF selection
            int ptNum = mBinning -> GetPtBinNum(runType, typeIdx, dleIdx);
            int xfNum = mBinning -> GetXfBinNum(runType, typeIdx, dleIdx);
            for(int pt=0; pt<ptNum; pt++){
                double ptLowerBoundary = mBinning -> GetPtBinBoundary(runType, typeIdx, dleIdx, pt);
                double ptUpperBoundary = mBinning -> GetPtBinBoundary(runType, typeIdx, dleIdx, pt+1);

                for(int xf=0; xf<xfNum; xf++){
                    double xfLowerBoundary = mBinning -> GetXfBinBoundary(runType, typeIdx, dleIdx, xf);
                    double xfUpperBoundary = mBinning -> GetXfBinBoundary(runType, typeIdx, dleIdx, xf+1);
                    
                    if(ptLowerBoundary <= particlePt && particlePt < ptUpperBoundary){
                        if(xfLowerBoundary <= particleXf && particleXf < xfUpperBoundary){
                            // signal particle mass selection
                            double massLowerBoundary = mMassFitting -> GetMassLowerBoundary(runType, typeIdx, dleIdx, pt, xf);
                            double massUpperBoundary = mMassFitting -> GetMassUpperBoundary(runType, typeIdx, dleIdx, pt, xf);
                            if(massLowerBoundary <= m && m < massUpperBoundary || particleRunIdx == kNeutron || particleRunIdx == kGammaRun){
                                mDilution -> FillAngle(runType, typeIdx, dleIdx, pt, xf, phiAngle);

                                cout << phiAngle << endl;
                                if(0. <= phiAngle && phiAngle < 180.){
                                    if(blueSpinUp == true){mAsymmetry -> FillPolarization(fillIdx, typeIdx, dleIdx, pt, xf, false);} // left 
                                    else{mAsymmetry -> FillPolarization(fillIdx, typeIdx, dleIdx, pt, xf, true);} // right 
                                }
                                else{
                                    if(blueSpinUp == true){mAsymmetry -> FillPolarization(fillIdx, typeIdx, dleIdx, pt, xf, true);} // right
                                    else{mAsymmetry -> FillPolarization(fillIdx, typeIdx, dleIdx, pt, xf, false);} // left
                                }
                            }

                            // BackGround particle mass selection
                            double massBkgLowerBoundary = mMassFitting -> GetMassBkgLowerBoundary(runType, typeIdx, dleIdx, pt, xf);
                            double massBkgUpperBoundary = mMassFitting -> GetMassBkgUpperBoundary(runType, typeIdx, dleIdx, pt, xf);
                            if((0. < m && m < massBkgLowerBoundary) || m > massBkgUpperBoundary){
                                if(0. <= phiAngle && phiAngle < 180.){
                                    if(blueSpinUp == true){mAsymmetry -> FillBkgPolarization(fillIdx, typeIdx, dleIdx, pt, xf, false);} // left 
                                    else{mAsymmetry -> FillBkgPolarization(fillIdx, typeIdx, dleIdx, pt, xf, true);} // right 
                                }
                                else{
                                    if(blueSpinUp == true){mAsymmetry -> FillBkgPolarization(fillIdx, typeIdx, dleIdx, pt, xf, true);} // right
                                    else{mAsymmetry -> FillBkgPolarization(fillIdx, typeIdx, dleIdx, pt, xf, false);} // left
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    mDilution -> Calculate();
    mDilution -> SaveDilutionData();
    mFigureDrawing -> DrawDilutionHist();    

    mAsymmetry -> Calculate();
    mAsymmetry -> SaveAsymmetryData();

}