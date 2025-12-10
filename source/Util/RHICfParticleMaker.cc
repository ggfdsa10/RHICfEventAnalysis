#include "RHICfParticleMaker.hh"

RHICfParticleMaker::RHICfParticleMaker() 
{
    mOptContainer = RHICfOptContainer::GetOptContainer();
    mIsSimData = false;
}

RHICfParticleMaker::~RHICfParticleMaker()
{
}

void RHICfParticleMaker::Init()
{
    mBeamCenterMethod = mOptContainer -> GetBeamCenterMethod();
    mBeamCenterRef = mOptContainer -> GetBeamCenterRef();
}

void RHICfParticleMaker::SetEventDst(StRHICfEventDst* eventDst)
{
    mEventDst = eventDst;
    mIsSimData = false;
}

void RHICfParticleMaker::SetSimDst(StRHICfSimDst* simDst)
{
    mSimDst = simDst;
    mIsSimData = true;
}

MiniParticle RHICfParticleMaker::GetMiniParticle()
{
    if(mOptContainer->GetParticleRunIdx() == kGammaRun){return GetGamma();}
    if(mOptContainer->GetParticleRunIdx() == kPi0Run){return GetPi0();}
    if(mOptContainer->GetParticleRunIdx() == kNeutronRun){return GetNeutron();}
    if(mOptContainer->GetParticleRunIdx() == kLambda0Run){return GetLambda0();}
    
    MiniParticle null;
    return null;
}

MiniParticle RHICfParticleMaker::GetGamma()
{
    MiniParticle particle;
    mEvent = mEventDst -> GetEvent();
    mDetHit = mEventDst -> GetRHICfDetHit();
    int pointNum = mEventDst -> GetRHICfDetPointNum();
    int gammaNum = 0;
    for(int i=0; i<pointNum; i++){
        mDetPoint = mEventDst -> GetRHICfDetPoint(i);

        int pid = mDetPoint -> GetPID();
        if(pid != kPhoton){continue;}

        double x = mDetPoint->GetPointPos(0);
        double y = mDetPoint->GetPointPos(1);
        double e = mDetPoint->GetPointEnergy(0);
        double r = std::sqrt(x*x + y*y + kRHICfPosZ*kRHICfPosZ);
        double px = x/r*e;
        double py = y/r*e;
        double pz = kRHICfPosZ/r*e;
        double pt = std::sqrt(px*px + py*py);

        double beamEnergy = mEvent->GetBeamEnergy();
        double xf = pz/beamEnergy;

        particle.label.push_back(0);
        particle.towerIdx.push_back(mDetPoint->GetTowerIdx());
        particle.x.push_back(x);
        particle.y.push_back(y);
        particle.e.push_back(e);
        particle.pt.push_back(pt);
        particle.xf.push_back(xf);
        particle.px.push_back(px);
        particle.py.push_back(py);
        particle.pz.push_back(pz);
        

        gammaNum++;
    }
    particle.particleNum = gammaNum;

    return particle;
}

MiniParticle RHICfParticleMaker::GetPi0()
{
    MiniParticle particle;
    particle.particleNum = 0;

    int pointNum = 0;
    int fillNum = 0;
    int runType = kALLRun;
    
    if(!mIsSimData){
        mEvent = mEventDst -> GetEvent();
        mDetHit = mEventDst -> GetRHICfDetHit();
        pointNum = mEventDst -> GetRHICfDetPointNum();
        fillNum = mEvent -> GetFillNumber();
        runType = mOptContainer -> GetFillNumToRunIdx(fillNum);
    }
    if(mIsSimData){
        mSimEvent = mSimDst -> GetSimEvent();
        mSimRHICfHit = mSimDst -> GetSimRHICfHit();
        pointNum = mSimDst -> GetSimRHICfPointNum();
        runType = mSimEvent -> GetRHICfRunType();
        
        if(runType == kTLRun){fillNum = 21142;}
        if(runType == kTSRun){fillNum = 21148;}
        if(runType == kTOPRun){fillNum = 21149;}
    }
    if(pointNum != 2){return particle;}

    int geoLabel = -1;
    int pi0Type = -1;
    double gammaX[2];
    double gammaY[2];
    double gammaE[2];
    double gammaTowerIdx[2];

    memset(gammaX, 0., sizeof(gammaX));
    memset(gammaY, 0., sizeof(gammaY));
    memset(gammaE, 0., sizeof(gammaE));
    memset(gammaTowerIdx, 0., sizeof(gammaTowerIdx));

    int gammaNum = 0;
    for(int i=0; i<pointNum; i++){
        int pid = -1.;
        int towerIdx = -1.;
        double x = -1.;
        double y = -1.;
        double energy = -1.;
        double l90 = -1.;

        if(!mIsSimData){
            mDetPoint = mEventDst -> GetRHICfDetPoint(i);
            pid = mDetPoint -> GetPID();
            towerIdx = mDetPoint->GetTowerIdx();
            x = mDetPoint->GetPointPos(0);
            y = mDetPoint->GetPointPos(1);
            energy = mDetPoint->GetPointEnergy(0);
            l90 = mDetHit -> GetL90(towerIdx);
        }
        if(mIsSimData){
            mSimRHICfPoint = mSimDst -> GetSimRHICfPoint(i);
            pid = mSimRHICfPoint -> GetPID();
            towerIdx = mSimRHICfPoint -> GetTowerIdx();
            x = mSimRHICfPoint -> GetPointPos(0);
            y = mSimRHICfPoint -> GetPointPos(1);
            energy = mSimRHICfPoint -> GetPointEnergy(0);
            l90 = mSimRHICfHit -> GetL90(towerIdx);
        }

        if(pid != kPhoton){continue;}
        if(energy < 2.){continue;}
        if(l90 < 8. || l90 > 18.){continue;}
        if(fabs(x - y) < 0.00001){continue;}
        double towerSize = (towerIdx == kTS)? 20. : 40.;

        if(2. > x || x > towerSize-2.){continue;}
        if(2. > y || y > towerSize-2.){continue;}

        if(4. < x && x < towerSize-4.){
            if(4. < y && y < towerSize-4.){
                geoLabel = 2;
            }
        }

        int isCut = GetGlobalPosition(x, y, towerIdx, fillNum, mBeamCenterMethod, mBeamCenterRef);
        if(isCut == 0){continue;}

        gammaX[gammaNum] = x;
        gammaY[gammaNum] = y;
        gammaE[gammaNum] = energy;
        gammaTowerIdx[gammaNum] = towerIdx;
        gammaNum++;
    }

    if(gammaNum != 2){return particle;}

    if(gammaTowerIdx[0] == gammaTowerIdx[1]){
        if(!mIsSimData){
            if(!mEvent->GetRHICfHighEMTrig()){return particle;}
        }
        if(mIsSimData){
            // if(!mSimEvent->IsHighEMTrigger()){return particle;}
        }
        if(gammaTowerIdx[0] == kTS){return particle;}
        if(fabs(gammaX[0]-gammaX[1]) < 0.000001){return particle;}
        if(fabs(gammaY[0]-gammaY[1]) < 0.000001){return particle;}
        pi0Type = 2;
        if(mIsSimData){
            // New simulated energy correction 
            // for(int i=0; i<2; i++){
            //     if(runType == kTLRun){
            //         if(gammaTowerIdx[i] == kTL){gammaE[i] = gammaE[i]*0.97 + 1.5*0.975;}
            //     }
            //     if(runType == kTSRun){
            //         if(gammaTowerIdx[i] == kTL){gammaE[i] = gammaE[i]*0.964 + 2.072*0.969;}
            //     }
            //     if(runType == kTOPRun){
            //         if(gammaTowerIdx[i] == kTL){gammaE[i] = gammaE[i]*0.965 + 1.623*0.970;}
            //     }
            // }
        }
        // New energy correction for Mass peak 
        if(!mIsSimData){
            for(int i=0; i<2; i++){
                if(gammaE[i] < 140. || gammaE[i] > 20.){
                    gammaE[i] = gammaE[i] + 2.*exp(-0.05*(gammaE[i])-20) - gammaE[i]*0.001 + 0.5;
                }
            }
        }
    }
    if(gammaTowerIdx[0] != gammaTowerIdx[1]){
        if(!mIsSimData){
            if(!mEvent->GetRHICfPi0Trig()){return particle;}
        }
        if(mIsSimData){
            // if(!mSimEvent->IsType1Pi0Trigger()){return particle;}
        }

        pi0Type = 1;
        if(mIsSimData){
            // New simulated energy correction 
            // for(int i=0; i<2; i++){
            //     if(runType == kTLRun){
            //         if(gammaTowerIdx[i] == kTS){gammaE[i] = gammaE[i]*0.964 + 0.9*0.964;}
            //         if(gammaTowerIdx[i] == kTL){gammaE[i] = gammaE[i]*0.979 + 0.0408*0.979;}
            //     }
            //     if(runType == kTSRun){
            //         if(gammaTowerIdx[i] == kTS){gammaE[i] = gammaE[i]*0.975 + 0.7232*0.975;}
            //         if(gammaTowerIdx[i] == kTL){gammaE[i] = gammaE[i]*0.967 + 0.519*0.967;}
            //     }
            //     if(runType == kTOPRun){
            //         if(gammaTowerIdx[i] == kTS){gammaE[i] = gammaE[i]*0.965 + 1.619*0.965;}
            //         if(gammaTowerIdx[i] == kTL){gammaE[i] = gammaE[i]*0.967 + 0.5*0.967;}
            //     }
            // }
        }
        // New energy correction for Mass peak
        if(!mIsSimData){
            for(int i=0; i<2; i++){
                if(gammaE[i] < 140. || gammaE[i] > 20.){
                    if(gammaTowerIdx[i] == kTS){
                        gammaE[i] = gammaE[i] + 2.*exp(-0.05*(gammaE[i] - 13.)) - gammaE[i]*0.001 + 0.15;
                    }
                    if(gammaTowerIdx[i] == kTL){
                        gammaE[i] = gammaE[i] + 1.8*exp(-0.05*(gammaE[i])) - gammaE[i]*0.001 + 0.15;
                    }
                }
            }
        }
    }

    if(gammaE[0] < 20 || gammaE[1] < 20){return particle;}

    double x = gammaX[0] + gammaX[1];
    double y = gammaY[0] + gammaY[1];

    double z = 17.812; // [m]
    double r1 = sqrt(gammaX[0]*gammaX[0] + gammaY[0]*gammaY[0] + z*z);
    double r2 = sqrt(gammaX[1]*gammaX[1] + gammaY[1]*gammaY[1] + z*z);

    double px = (gammaX[0]/r1)*gammaE[0] + (gammaX[1]/r2)*gammaE[1]; 
    double py =(gammaY[0]/r1)*gammaE[0] + (gammaY[1]/r2)*gammaE[1]; 
    double pz = (z/r1)*gammaE[0] + (z/r2)*gammaE[1]; 
    double e = gammaE[0]+gammaE[1];
    double pt = sqrt(px*px + py*py);
    double xf = pz/255.;
    double m = sqrt(e*e - px*px - py*py - pz*pz)*1000.;

    particle.label.push_back(geoLabel);
    particle.towerIdx.push_back(-1);
    particle.type.push_back(pi0Type);
    particle.x.push_back(x);
    particle.y.push_back(y);
    particle.e.push_back(e);
    particle.pt.push_back(pt);
    particle.xf.push_back(xf);
    particle.px.push_back(px);
    particle.py.push_back(py);
    particle.pz.push_back(pz);
    particle.m.push_back(m);
    particle.particleNum = 1;

    return particle;
}

MiniParticle RHICfParticleMaker::GetNeutron()
{

    MiniParticle particle;
    particle.particleNum = 0;
    int pointNum = 0;
    int fillNum = 0;
    int runType = kALLRun;
    
    if(!mIsSimData){
        mEvent = mEventDst -> GetEvent();
        mDetHit = mEventDst -> GetRHICfDetHit();
        pointNum = mEventDst -> GetRHICfDetPointNum();
        fillNum = mEvent -> GetFillNumber();
        runType = mOptContainer -> GetFillNumToRunIdx(fillNum);
    }
    if(mIsSimData){
        mSimEvent = mSimDst -> GetSimEvent();
        mSimRHICfHit = mSimDst -> GetSimRHICfHit();
        pointNum = mSimDst -> GetSimRHICfPointNum();
        runType = mSimEvent -> GetRHICfRunType();
        
        if(runType == kTLRun){fillNum = 21142;}
        if(runType == kTSRun){fillNum = 21148;}
        if(runType == kTOPRun){fillNum = 21149;}
    }
    if(pointNum != 1){return particle;}

    if(!mIsSimData){
        if(!mEvent->GetRHICfShowerTrig()){return particle;}
    }
    if(mIsSimData){
        if(!mSimEvent->IsShowerTrigger()){return particle;}
    }

    int geoLabel = -1;
    int pid = -1.;
    int towerIdx = -1.;
    double x = -1.;
    double y = -1.;
    double energy = -1.;
    double l90 = -1.;

    if(!mIsSimData){
        mDetPoint = mEventDst -> GetRHICfDetPoint(0);
        pid = mDetPoint -> GetPID();
        towerIdx = mDetPoint->GetTowerIdx();
        x = mDetPoint->GetPointPos(0);
        y = mDetPoint->GetPointPos(1);
        energy = mDetPoint->GetPointEnergy(0);
        l90 = mDetHit -> GetL90(towerIdx);
    }
    if(mIsSimData){
        mSimRHICfPoint = mSimDst -> GetSimRHICfPoint(0);
        pid = mSimRHICfPoint -> GetPID();
        towerIdx = mSimRHICfPoint -> GetTowerIdx();
        x = mSimRHICfPoint -> GetPointPos(0);
        y = mSimRHICfPoint -> GetPointPos(1);
        energy = mSimRHICfPoint -> GetPointEnergy(0);
        l90 = mSimRHICfHit -> GetL90(towerIdx);
    }

    if(pid != kHadron){return particle;}
    if(energy < 2.){return particle;}
    if(l90 < 20.){return particle;}
    if(fabs(x - y) < 0.00001){return particle;}
    double towerSize = (towerIdx == kTS)? 20. : 40.;

    if(2. > x || x > towerSize-2.){return particle;}
    if(2. > y || y > towerSize-2.){return particle;}

    if(4. < x && x < towerSize-4.){
        if(4. < y && y < towerSize-4.){
            geoLabel = 2;
        }
    }

    int isCut = GetGlobalPosition(x, y, towerIdx, fillNum, mBeamCenterMethod, mBeamCenterRef);
    if(isCut == 0){return particle;}

    // double z = 17.8; // RHICf detector surface distance
    // double r = sqrt(x*x* + y*y + z*z);

    particle.label.push_back(geoLabel);
    particle.towerIdx.push_back(towerIdx);
    particle.type.push_back(-1);
    particle.x.push_back(x);
    particle.y.push_back(y);
    particle.e.push_back(energy);
    particle.pt.push_back(-1);
    particle.xf.push_back(-1);
    particle.px.push_back(-1);
    particle.py.push_back(-1);
    particle.pz.push_back(-1);
    particle.m.push_back(-1.);
    particle.particleNum = 1;

    return particle;  
}

MiniParticle RHICfParticleMaker::GetLambda0()
{
    MiniParticle particle;
    return particle;  
}

int RHICfParticleMaker::GetGlobalPosition(double& localX, double& localY, int towerIdx, int fillNum, int method, int ref)
{
    // for convert Collision Point system, define the system in TS origin
    Int_t TowerIdxByRunType = -999;
    Double_t distOriginToOrigin = distTStoTL - sqrt(2.)*geoCenterTL + sqrt(2.)*geoCenterTS;

    double beamCenterX = 0.;
    double beamCenterY = 0.;
    // if(!mIsSimData){ test
        beamCenterX = BeamCenterPosition(fillNum, 0, method, ref);
        beamCenterY = BeamCenterPosition(fillNum, 1, method, ref);
    // }

    int runType = mOptContainer -> GetFillNumToRunIdx(fillNum);
    if(runType == kTLRun){
        TowerIdxByRunType = 1;
        beamCenterY = beamCenterY + sqrt(2.)*geoCenterTL;
    }
    if(runType == kTSRun || runType == kTOPRun){
        TowerIdxByRunType = 0;
        beamCenterY = beamCenterY + sqrt(2.)*geoCenterTS;
        if(runType == kTOPRun && mIsSimData){beamCenterY -= 21.6;} // [mm]
    }

    double posX = localX;
    double posY = localY;

    // the Hit Position convert to one tower coordinate
    if(TowerIdxByRunType == 0 && towerIdx == 1){
        posX += distOriginToOrigin/sqrt(2.);
        posY += distOriginToOrigin/sqrt(2.);
    }
    if(TowerIdxByRunType == 1 && towerIdx == 0){
        posX -= distOriginToOrigin/sqrt(2.);
        posY -= distOriginToOrigin/sqrt(2.);
    }

    // Hit position rotate to -45 degree.
    double tmpPosX = (posX - posY)/sqrt(2.);
    double tmpPosY = (posX + posY)/sqrt(2.);

    double type2TOPGeoCut = 58. + sqrt(2)*geoCenterTS;
    if(runType == kTOPRun && tmpPosY > type2TOPGeoCut){return 0;} // cut for TOP run 

    // convert to beam center and IP system [m]
    localX = (tmpPosX - beamCenterX)*0.001;
    localY = (tmpPosY - beamCenterY)*0.001;

    return 1;
}

double RHICfParticleMaker::BeamCenterPosition(int fillNum, int xy, int method, int ref)
{
    if(fillNum == 21142){
        if(method==kBeamCenterScan){
            if(xy==0){return 0.;}
            else if(xy==1){return 2.21;}
        }
        else if(method==kBeamCenterHit){
            if(xy==0){return 0.66;}
            else if(xy==1){return 1.37;}
        }

    }
    if(fillNum == 21145){
        if(method==kBeamCenterScan){
            if(xy==0){return 0.;}
            else if(xy==1){return 2.31;}
        }
        else if(method==kBeamCenterHit){
            if(xy==0){return 0.22;}
            else if(xy==1){return 1.55;}
        }
    }

    if(fillNum == 21148){
        if(method==kBeamCenterScan){
            if(xy==0){return 0.;}
            else if(xy==1){return 2.36;}
        }
        else if(method==kBeamCenterHit){
            if(xy==0){return 0.22;}
            else if(xy==1){return 1.37;}
        }
    }
    if(fillNum == 21150){
        if(method==kBeamCenterScan){
            if(xy==0){return 0.;}
            else if(xy==1){return 2.33;}
        }
        else if(method==kBeamCenterHit){
            if(xy==0){return 0.22;}
            else if(xy==1){return 1.55;}
        }
    }
    if(fillNum == 21149){
        if(ref==kBeamCenterRef21148){ // refer to fill 21148
            if(method==kBeamCenterScan){
                if(xy==0){return 0.;}
                else if(xy==1){return -21.67;}
            }
            else if(method==kBeamCenterHit){
                if(xy==0){return 0.22;}
                else if(xy==1){return -22.66;}
            }
        }
        if(ref==kBeamCenterRef21150){ // refer to fill 21150
            if(method==kBeamCenterScan){
                if(xy==0){return 0.;}
                else if(xy==1){return -21.71;}
            }
            else if(method==kBeamCenterHit){
                if(xy==0){return 0.22;}
                else if(xy==1){return -22.48;}
            }
        }
    }
    return -999.;
}
