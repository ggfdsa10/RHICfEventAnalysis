#include "RHICfSimDstReader.hh"

#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TObjArray.h"
#include "TObjString.h"

RHICfSimDstReader::RHICfSimDstReader()
{
    mOptContainer = RHICfOptContainer::GetOptContainer();
}

RHICfSimDstReader::~RHICfSimDstReader()
{
}

void RHICfSimDstReader::Init()
{
	int runType = mOptContainer->GetRunType();
    int modelIdx = mOptContainer->GetSimModelIdx();
    mDataPath = mOptContainer->GetDataPath();
	mInputDataPath = mOptContainer->GetSimInputDataPath();
	mModelName = mOptContainer->GetSimModelName(modelIdx);
    mDataList = mDataPath+"/RHICfSimDst_"+mModelName+"_"+mOptContainer->GetRunTypeName(runType)+".list";

    TString tag = mOptContainer->GetParticleRunName();
    TString tag2 = mOptContainer->GetTag();
    if(tag != ""){tag = "_" + tag;}
    if(tag2 != ""){tag = tag + "_" + tag2;}

	if(mOptContainer->GetInputDataList() != ""){
        mDataList = mOptContainer->GetInputDataList();
        mMiniDstName = "MiniSimDst_"+mModelName+tag;
    }
	else{
		gSystem->Exec(Form(" > %s", mDataList.Data()));
		if(runType == kTSRun){
            gSystem->Exec(Form("ls %s/TS/StarSim_%s_TS* > %s", mInputDataPath.Data(), mModelName.Data(), mDataList.Data()));
            mMiniDstName = "MiniSimDst_"+mModelName+"_TS"+tag;
		}
		else if(runType == kTLRun){
			gSystem->Exec(Form("ls %s/TL/StarSim_%s_TL* > %s", mInputDataPath.Data(), mModelName.Data(), mDataList.Data()));
            mMiniDstName = "MiniSimDst_"+mModelName+"_TL"+tag;
		}
		else if(runType == kTOPRun){
			gSystem->Exec(Form("ls %s/TOP/StarSim_%s_TOP* > %s", mInputDataPath.Data(), mModelName.Data(), mDataList.Data()));
            mMiniDstName = "MiniSimDst_"+mModelName+"_TOP"+tag;
		}
        else if(runType == kALLRun){
            cout << "RHICfSimDstReader::Init() -- RHICfSimDstReader should be seperated by runs.."<< endl;
        }
	}

    bool findFlag = FindMiniSimDst();
    if(findFlag){
        cout << "RHICfSimDstReader::Init() -- Find a " << mMiniDstName << ".root" << endl;

        mChain = new TChain("MiniSimDst");
        mChain -> Add(Form("%s/%s.root", mDataPath.Data(), mMiniDstName.Data()));
        mMiniSimDst = 0;

        mChain -> SetBranchAddress("MiniSimDst", &mMiniSimDst);
        mEventNum = mChain -> GetEntries();
    }
    else{
        cout << "RHICfSimDstReader::Init() -- Cannot find a " << mMiniDstName << ".root" << endl;
        mMiniSimDst = 0;
        mEventNum = 0;
    }
}

void RHICfSimDstReader::Make()
{
    if(FindMiniSimDst()){cout << "RHICfSimDstReader::Make() -- Duplicated file -> " << mMiniDstName << ".root... recreate file." << endl;}
    else{cout << "RHICfSimDstReader::Make() -- Create file: " << mMiniDstName << ".root" << endl;}

    mMiniDstFile = new TFile(Form("%s/%s.root", mDataPath.Data(), mMiniDstName.Data()), "recreate");
    mMiniDstTree = new TTree("MiniSimDst", "MiniSimDst");
    mMiniSimDst = new StRHICfMiniSimDst();
    mMiniDstTree -> Branch("MiniSimDst", &mMiniSimDst);

    std::ifstream inputStream(mDataList.Data());
    std::string file;
    int fileNum = 0;

    mChain = new TChain("StRHICfSimDst");
    cout << "RHICfSimDstReader::Make() -- " << mModelName << " Data list reading...." << endl;

    while(getline(inputStream, file)){
        if(fileNum == mOptContainer->GetExecuteFileNum() && mOptContainer->GetExecuteFileNum() > 0){break;}
        if(file.find("reco.rhicfsim.RHICfSimDst.root") != std::string::npos){
            TString fileName = file;
            mChain -> Add(fileName);
            if(fileNum%10 == 0){cout << mModelName << " processing " << fileNum << " th file... Entries --> "  << mChain -> GetEntries() << endl;}
            fileNum++;
        }
    }

    mSimDst = new StRHICfSimDst();
    mSimDst -> ReadDstArray(mChain);

    mEventNum = mChain -> GetEntries();
    cout << "RHICfSimDstReader::Make() -- Number of RHICfEventDst file: " << fileNum << " has initialized.. with eventNum: " << mEventNum << endl;
    if(mOptContainer->GetExecuteEventNum() < mEventNum && mOptContainer->GetExecuteEventNum() > 0){
        mEventNum = mOptContainer->GetExecuteEventNum();
        cout << "RHICfSimDstReader::Make() -- Force to total event number : " << mEventNum << endl;
    }

    MakeMiniSimDst();

    mMiniDstFile -> cd();
    mMiniDstTree -> Write();
    mMiniDstFile -> Close();

    cout << "RHICfSimDstReader::Make() -- " << mModelName << " MiniSimDst has been created." << endl;
    Init();
}

Int_t RHICfSimDstReader::GetEventNum()
{
	return mEventNum;
}

StRHICfMiniSimDst* RHICfSimDstReader::GetMiniSimDst(int idx)
{
	if(idx < mEventNum){
		mChain -> GetEntry(idx);
		return mMiniSimDst;
	}
	return 0;
}

bool RHICfSimDstReader::FindMiniSimDst()
{
    // Find a directory
    TList *listOfDirs;
    TObject *objDir;

    TSystemDirectory dir("dir", mDataPath);
    listOfDirs = dir.GetListOfFiles();
    TIter next(listOfDirs);

    while((objDir = next())){
        TSystemFile* dirPtr = dynamic_cast<TSystemFile*>(objDir);
        if(dirPtr && !dirPtr->IsDirectory()){
            TString dirName = dirPtr->GetName();
            if(dirName.Index(mMiniDstName) != -1){
                return true;
            }
        }
    }
    return false;
}

void RHICfSimDstReader::MakeMiniSimDst()
{
	cout << "RHICfSimDstReader::MakeMiniSimDst() -- making of MiniSimDst has started." << endl;

    mParticleMaker = new RHICfParticleMaker();
    mDLECondition = new RHICfDLECondition();

    mParticleMaker -> Init(); 
    mDLECondition -> Init();

	int totalEventNum = GetEventNum();
    for(int event=0; event<totalEventNum; event++){
        if(event%1000 == 0){cout << "RHICfSimDstReader::MakeMiniSimDst() -- Event: " << event << " / " << totalEventNum << endl;}
        mChain -> GetEntry(event);

        mMiniSimDst -> Clear();

        mParticleMaker -> SetSimDst(mSimDst);
        MiniParticle particles = mParticleMaker -> GetMiniParticle();
        int particleNum = particles.particleNum;
        if(particleNum == 0){continue;}
        mSimEvent = mSimDst -> GetSimEvent();
        int modelIdx = mSimEvent -> GetGeneratorIdx();
        int processID = mSimEvent -> GetProcessId();
        int eventType = mDLECondition -> GetSimEventType(modelIdx, processID);
        mMiniSimDst -> SetProcessID(processID);
        mMiniSimDst -> SetEventType(eventType);
        mMiniSimDst -> SetRHICfRunType(mSimEvent->GetRHICfRunType());
        mMiniSimDst -> SetGeneratorIdx(modelIdx);

        bool isShower = mSimEvent -> IsShowerTrigger();
        bool isPi0Trig = mSimEvent -> IsType1Pi0Trigger();
        bool isEMTrig = mSimEvent -> IsHighEMTrigger();

        if(!isShower && !isPi0Trig && !isEMTrig){
            cout << " particleNum " << particleNum << endl;
        }
        
        // STAR Det
        mSimBBC = mSimDst -> GetSimBBC();
        int btofMult = mSimDst -> GetSimBTofNum();
        int bbcSmallEast = mSimBBC -> GetSmallSum(kBeamEast);
        int bbcSmallWest = mSimBBC -> GetSmallSum(kBeamWest);
        int bbcLargeEast = mSimBBC -> GetLargeSum(kBeamEast);
        int bbcLargeWest = mSimBBC -> GetLargeSum(kBeamWest);

        mMiniSimDst -> SetBTofMult(btofMult);
        mMiniSimDst -> SetBBCSumADC(rEast, 0, bbcSmallEast);
        mMiniSimDst -> SetBBCSumADC(rEast, 1, bbcLargeEast);
        mMiniSimDst -> SetBBCSumADC(rWest, 0, bbcSmallWest);
        mMiniSimDst -> SetBBCSumADC(rWest, 1, bbcLargeWest);

        // DLE condition
        mDLECondition -> SetBTofMult(btofMult);
        mDLECondition -> SetBBCADC(bbcSmallEast, bbcSmallWest, bbcLargeEast, bbcLargeWest);
        mMiniSimDst -> SetDLEIdx(mDLECondition -> GetDLEIdx());

        // RHICf Particles
        for(int i=0; i<particleNum; i++){
            int towerIdx = particles.towerIdx[i];
            int pi0Type = particles.type[i]-1;
            double x = particles.x[i];
            double y = particles.y[i];
            double px = particles.px[i];
            double py = particles.py[i];
            double pz = particles.pz[i];
            double energy = particles.e[i];
            double mass = particles.m[i];

            int pid = -1;
            if(towerIdx == -1){pid = 0;}
            else{pid = 1;}

            mMiniSimDst -> SetRHICfParticle(towerIdx, pid, x, y, px, py, pz, energy, pi0Type, mass);
        }

        // ZDC
        mSimZDC = mSimDst -> GetSimZDC();
        for(int pmt=0; pmt<rZDCPMTNum; pmt++){
            mMiniSimDst -> SetZDCPmtPhotonNum(pmt, mSimZDC->GetPmtPhotonNum(pmt));
            mMiniSimDst -> SetZDCPmtdE(pmt, mSimZDC->GetPmtdE(pmt));
        }
        for(int smdX=0; smdX<rSMDXNum; smdX++){
            mMiniSimDst -> SetZDCSMDdE(0, smdX, mSimZDC->GetSMDdE(0, smdX));
        }
        for(int smdY=0; smdY<rSMDYNum; smdY++){
            mMiniSimDst -> SetZDCSMDdE(1, smdY, mSimZDC->GetSMDdE(1, smdY));
        }

        // RHICf Plate Energy
        mSimRHICfHit = mSimDst -> GetSimRHICfHit();
        for(int it=0; it<2; it++){
            for(int ip=0; ip<16; ip++){
                double e = mSimRHICfHit->GetPlatedE(it, ip);
                mMiniSimDst -> SetRHICfPlatedE(it, ip, e);
            }
        }

        SaveTruthSimTracks();

        mMiniDstTree -> Fill();
    }
}

void RHICfSimDstReader::SaveTruthSimTracks()
{
    for(int it=0; it<2; it++){
        for(int trk=0; trk<mSimRHICfHit->GetSimTrkNum(it); trk++){
            int truthIdx = mSimRHICfHit -> GetSimTrkId(it, trk);
            double posX = mSimRHICfHit -> GetSimTrkIncidentPos(it, trk, 0);
            double posY = mSimRHICfHit -> GetSimTrkIncidentPos(it, trk, 1);
            double energy = mSimRHICfHit -> GetSimTrkIncidentEnergy(it, trk);
            mMiniSimDst -> SetRHICfTruthId(it, truthIdx, posX, posY, energy);
        }
    }
    for(int trk=0; trk<mSimZDC->GetSimTrkNum(); trk++){
        int truthIdx = mSimZDC -> GetSimTrkId(trk);
        double posX = mSimZDC -> GetSimTrkIncidentPos(trk, 0);
        double posY = mSimZDC -> GetSimTrkIncidentPos(trk, 1);
        double energy = mSimZDC -> GetSimTrkIncidentEnergy(trk);
        mMiniSimDst -> SetZDCTruthId(truthIdx, posX, posY, energy);
    }

    int simTrackNum = mSimDst->GetSimTrackNum(); 
    for(int i=0; i<simTrackNum; i++){
        mSimTrack = mSimDst -> GetSimTrack(i);

        bool isPrimary = mSimTrack -> IsPrimary();
        bool isPropagate = mSimTrack -> IsSimPropagate();
        bool isRHICfHit = mSimTrack -> IsRHICfHit();
        bool isFinal = mSimTrack -> IsFinal();

        int id = mSimTrack -> GetId();
        int parentId = mSimTrack -> GetParentId();
        int pid = mSimTrack -> GetPid();
        double energy = mSimTrack -> GetE();
        double px = mSimTrack -> GetPx();
        double py = mSimTrack -> GetPy();
        double pz = mSimTrack -> GetPz();
        double vx = mSimTrack -> GetVxStart();
        double vy = mSimTrack -> GetVyStart();
        double vz = mSimTrack -> GetVzStart();
        double ex = mSimTrack -> GetVxEnd();
        double ey = mSimTrack -> GetVyEnd();
        double ez = mSimTrack -> GetVzEnd();

        mMiniSimDst -> SetSimTrack(id, parentId, pid, energy, px, py, pz, vx, vy, vz, ex, ey, ez, isPrimary, isPropagate, isRHICfHit, isFinal);
    }
}