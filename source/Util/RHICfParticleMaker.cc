#include "RHICfParticleMaker.hh"

RHICfParticleMaker::RHICfParticleMaker() 
{
    mOptContainer = RHICfOptContainer::GetOptContainer();
}

RHICfParticleMaker::~RHICfParticleMaker()
{
}

void RHICfParticleMaker::Init()
{
}

void RHICfParticleMaker::SetEventDst(StRHICfEventDst* eventDst)
{
    mEventDst = eventDst;
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
    int particleNum = mEventDst -> GetRHICfParticleNum();
    int pi0Num = 0;
    for(int i=0; i<particleNum; i++){
        mParticle = mEventDst -> GetRHICfParticle(i);
        int pid = mParticle -> GetPID();
        if(pid != kPi0){continue;}

        int pi0Type = mParticle -> GetPi0Type();
        int towerIdx = mParticle -> GetTowerIdx();
        double x = mParticle -> GetRHICfHitPosX();
        double y = mParticle -> GetRHICfHitPosY();
        double e = mParticle -> GetEnergy();
        double pt = mParticle -> GetPt();
        double xf = mParticle -> GetXf();
        double px = mParticle -> GetPx();
        double py = mParticle -> GetPy();
        double pz = mParticle -> GetPz();
        double m = mParticle -> GetMass();

        if(pi0Type == 2 && towerIdx == kTS){continue;}

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

        pi0Num++;
    }
    particle.particleNum = pi0Num;

    return particle;
}

MiniParticle RHICfParticleMaker::GetNeutron()
{
    MiniParticle particle;
    int particleNum = mEventDst -> GetRHICfParticleNum();
    int neutronNum = 0;
    for(int i=0; i<particleNum; i++){
        mParticle = mEventDst -> GetRHICfParticle(i);
        int pid = mParticle -> GetPID();
        if(pid != kNeutron){continue;}

        int towerIdx = mParticle -> GetTowerIdx();
        int typeIdx = mParticle -> GetTowerIdx();
        double x = mParticle -> GetRHICfHitPosX();
        double y = mParticle -> GetRHICfHitPosY();
        double e = mParticle -> GetEnergy();
        double pt = mParticle -> GetPt();
        double xf = mParticle -> GetXf();
        double px = mParticle -> GetPx();
        double py = mParticle -> GetPy();
        double pz = mParticle -> GetPz();

        particle.towerIdx.push_back(towerIdx);
        particle.type.push_back(typeIdx+1);
        particle.x.push_back(x);
        particle.y.push_back(y);
        particle.e.push_back(e);
        particle.pt.push_back(pt);
        particle.xf.push_back(xf);
        particle.px.push_back(px);
        particle.py.push_back(py);
        particle.pz.push_back(pz);

        neutronNum++;
    }
    particle.particleNum = neutronNum;

    return particle;  
}

MiniParticle RHICfParticleMaker::GetLambda0()
{
    MiniParticle particle;
    return particle;  
}
