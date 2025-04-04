#include "RHICfEventDstReader.hh"

RHICfEventDstReader::RHICfEventDstReader() 
{
	mOptContainer = RHICfOptContainer::GetOptContainer();
}

RHICfEventDstReader::~RHICfEventDstReader()
{
}

void RHICfEventDstReader::Init()
{
	int runType = mOptContainer->GetRunType();
	TString inputDataPath = mOptContainer->GetInputDataPath();
    TString inputList = mOptContainer->GetDataPath()+"/RHICfEventDstFile.list";
	if(mOptContainer->GetInputDataList() != ""){inputList = mOptContainer->GetInputDataList();}
	else{
		gSystem->Exec(Form(" > %s", inputList.Data()));
		if(runType == kTSRun || runType == kALLRun){
			gSystem->Exec(Form("ls %s/TS/*RHICfEventDst.root >> %s", inputDataPath.Data(), inputList.Data()));
		}
		if(runType == kTLRun || runType == kALLRun){
			gSystem->Exec(Form("ls %s/TL/*RHICfEventDst.root >> %s", inputDataPath.Data(), inputList.Data()));
		}
		if(runType == kTOPRun || runType == kALLRun){
			gSystem->Exec(Form("ls %s/TOP/*RHICfEventDst.root >> %s", inputDataPath.Data(), inputList.Data()));
		}
	}

    std::ifstream inputStream(inputList.Data());
	std::string file;
	int fileNum = 0;

	mChain = new TChain("RHICfEventDst");

	cout << "RHICfEventDstReader::Init() -- Data list reading...." << endl;

	while(getline(inputStream, file)){
        if(file.find("RHICfEventDst.root") != std::string::npos){
            TString fileName = file;
            mChain -> Add(fileName);
            fileNum++;
        }
		if(fileNum > mOptContainer->GetExecuteFileNum()){break;}
    }

	mRHICfEventDst = new StRHICfEventDst();
	if(mOptContainer->GetOffDetPoint()){mChain -> SetBranchStatus("StRHICfDetPoint*", 0);}
	if(mOptContainer->GetOffParticle()){mChain -> SetBranchStatus("StRHICfParticle*", 0);}
	if(mOptContainer->GetOffTPCTrack() && mChain->GetBranchStatus("StRHICfTPCTrack")){mChain -> SetBranchStatus("StRHICfTPCTrack*", 0);}
    if(mOptContainer->GetOffBTof() && mChain->GetBranchStatus("StRHICfBTof")){mChain -> SetBranchStatus("StRHICfBTof*", 0);}
    if(mOptContainer->GetOffBBC() && mChain->GetBranchStatus("StRHICfBBC")){mChain -> SetBranchStatus("StRHICfBBC*", 0);}
    if(mOptContainer->GetOffVPD() && mChain->GetBranchStatus("StRHICfVPD")){mChain -> SetBranchStatus("StRHICfVPD*", 0);}
    if(mOptContainer->GetOffZDC() && mChain->GetBranchStatus("StRHICfZDC")){mChain -> SetBranchStatus("StRHICfZDC*", 0);}
    if(mOptContainer->GetOffFMS() && mChain->GetBranchStatus("StRHICfFMS")){mChain -> SetBranchStatus("StRHICfFMS*", 0);}
    if(mOptContainer->GetOffRPS() && mChain->GetBranchStatus("StRHICfRPS")){mChain -> SetBranchStatus("StRHICfRPS*", 0);}
	mRHICfEventDst -> ReadDstArray(mChain);

	mEventNum = mChain -> GetEntries();
	cout << "RHICfEventDstReader::Init() -- Number of RHICfEventDst file: " << fileNum << " has initialized.. with eventNum: " << mEventNum << endl;
	if(mOptContainer->GetExecuteEventNum() < mEventNum && mOptContainer->GetExecuteEventNum() > 0){
		mEventNum = mOptContainer->GetExecuteEventNum();
		cout << "RHICfEventDstReader::Init() -- Force to total event number : " << mEventNum << endl;
 	}
}

Int_t RHICfEventDstReader::GetEventNum()
{
	return mEventNum;
}

StRHICfEventDst* RHICfEventDstReader::GetEventDst(int idx)
{
	if(idx < mEventNum){
		mChain -> GetEntry(idx);
		return mRHICfEventDst;
	}
	return 0;
}