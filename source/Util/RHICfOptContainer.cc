#include "RHICfOptContainer.hh"

#include "TROOT.h"
#include "TSystem.h"
#include "TSystemFile.h"
#include "TSystemDirectory.h"
#include "TObjArray.h"
#include "TObjString.h"

RHICfOptContainer* RHICfOptContainer::mOptContainer = nullptr;
RHICfOptContainer* RHICfOptContainer::GetOptContainer()
{
    if (mOptContainer != nullptr){return mOptContainer;}
    return new RHICfOptContainer();
}

RHICfOptContainer::RHICfOptContainer() 
{
    mOptContainer = this;
    InitOptions();
    FindDirPath();
}

RHICfOptContainer::~RHICfOptContainer()
{
}

int RHICfOptContainer::Init()
{
    int sumParticleOpt = int(mGammaCalc)+int(mPi0Calc)+int(mNeutronCalc)+int(mLambda0Calc);
    if(sumParticleOpt != 1){
        cout << "RHICfOptContainer::Init() -- Error: set to only one particle calculation option!" << endl;
        return 0;
    }
    PrintOptions();

    return 1;
}

void RHICfOptContainer::SetRunType(TString type)
{
    type.ToUpper();
    if(type.Index("TS") != -1){mRunType = kTSRun;}
    else if(type.Index("TL") != -1){mRunType = kTLRun;}
    else if(type.Index("TOP") != -1){mRunType = kTOPRun;}
    else{mRunType = kALLRun;}
}

void RHICfOptContainer::SetInputDataPath(TString path){mInputDataPath = path;}
void RHICfOptContainer::SetInputDataList(TString listFile){mInputDataList = listFile;}
void RHICfOptContainer::SetExecuteEventNum(int event){mExecuteEventNum = event;}
void RHICfOptContainer::SetExecuteFileNum(int num){mExecuteFileNum = num;}

void RHICfOptContainer::ForceCalculateMass(){mForceCalculateMass = true;}
void RHICfOptContainer::ForceCalculateBinning(){mForceCalculateBinning = true;}
void RHICfOptContainer::ForceCalculateDilution(){mForceCalculateDilution = true;}
void RHICfOptContainer::ForceCalculateAsymmetry(){mForceCalculateAsymmetry = true;}
void RHICfOptContainer::ForceCalculateSystematicError(){mForceCalculateSystematicError = true;}

void RHICfOptContainer::ForceDefaultBinning(){mForceDefaultBinning = true;}

void RHICfOptContainer::CalculateGamma(){mGammaCalc = true;}
void RHICfOptContainer::CalculatePi0(){mPi0Calc = true;}
void RHICfOptContainer::CalculateNeutron(){mNeutronCalc = true;}
void RHICfOptContainer::CalculateLambda0(){mLambda0Calc = true;}

void RHICfOptContainer::SetOffDetPoint(){mOffDetPoint = true;}
void RHICfOptContainer::SetOffParticle(){mOffParticle = true;}
void RHICfOptContainer::SetOffTPCTrack(){mOffTPCTrack = true;}
void RHICfOptContainer::SetOffBTof(){mOffBTof = true;}
void RHICfOptContainer::SetOffBBC(){mOffBBC = true;}
void RHICfOptContainer::SetOffVPD(){mOffVPD = true;}
void RHICfOptContainer::SetOffZDC(){mOffZDC = true;}
void RHICfOptContainer::SetOffFMS(){mOffFMS = true;}
void RHICfOptContainer::SetOffRPS(){mOffRPS = true;}

int RHICfOptContainer::GetRunType(){return mRunType;}

int RHICfOptContainer::GetFillNumIdx(int fillNum)
{
    if(fillNum == 21142){return 0;}
    if(fillNum == 21145){return 1;}
    if(fillNum == 21148){return 2;}
    if(fillNum == 21149){return 3;}
    if(fillNum == 21150){return 4;}
    return -1;
}

int RHICfOptContainer::GetFillNumToRunIdx(int fillNum)
{
    if(fillNum == 21142 || fillNum == 21145){return kTLRun;}
    if(fillNum == 21148 || fillNum == 21150){return kTSRun;}
    if(fillNum == 21149){return kTOPRun;}
    return kALLRun;
}

int RHICfOptContainer::GetFillToRunIdx(int fillIdx)
{
    if(fillIdx == 0 || fillIdx == 1){return kTLRun;}
    if(fillIdx == 2 || fillIdx == 4){return kTSRun;}
    if(fillIdx == 3){return kTOPRun;}
    return kALLRun;
}

TString RHICfOptContainer::GetRunTypeName(int runIdx)
{
    if(runIdx == kTLRun){return "TL";}
    if(runIdx == kTSRun){return "TS";}
    if(runIdx == kTOPRun){return "TOP";}
    return "ALL";
}

TString RHICfOptContainer::GetDLEName(int dleIdx)
{
    if(dleIdx == kSDLE){return "SDLE";}
    if(dleIdx == kDDLE){return "DDLE";}
    if(dleIdx == kNDLE){return "NDLE";}
    return "ALL";
}

TString RHICfOptContainer::GetInputDataPath(){return mInputDataPath;}
TString RHICfOptContainer::GetInputDataList(){return mInputDataList;}
TString RHICfOptContainer::GetTablePath(){return mTablePath;}
TString RHICfOptContainer::GetDataPath(){return mDataPath;}
TString RHICfOptContainer::GetFigurePath(){return mFigurePath;}
int RHICfOptContainer::GetExecuteEventNum(){return mExecuteEventNum;}
int RHICfOptContainer::GetExecuteFileNum(){return mExecuteFileNum;}

bool RHICfOptContainer::IsForceCalculateMass(){return mForceCalculateMass;}
bool RHICfOptContainer::IsForceCalculateBinning(){return mForceCalculateBinning;}
bool RHICfOptContainer::IsForceCalculateDilution(){return mForceCalculateDilution;}
bool RHICfOptContainer::IsForceCalculateAsymmetry(){return mForceCalculateAsymmetry;}
bool RHICfOptContainer::IsForceCalculateSystematicError(){return mForceCalculateSystematicError;}

bool RHICfOptContainer::IsForceDefaultBinning(){return mForceDefaultBinning;}

bool RHICfOptContainer::GetOffDetPoint(){return mOffDetPoint;}
bool RHICfOptContainer::GetOffParticle(){return mOffParticle;}
bool RHICfOptContainer::GetOffTPCTrack(){return mOffTPCTrack;}
bool RHICfOptContainer::GetOffBTof(){return mOffBTof;}
bool RHICfOptContainer::GetOffBBC(){return mOffBBC;}
bool RHICfOptContainer::GetOffVPD(){return mOffVPD;}
bool RHICfOptContainer::GetOffZDC(){return mOffZDC;}
bool RHICfOptContainer::GetOffFMS(){return mOffFMS;}
bool RHICfOptContainer::GetOffRPS(){return mOffRPS;}

TString RHICfOptContainer::GetParticleRunName()
{
    if(mGammaCalc){return "Gamma";}
    if(mPi0Calc){return "Pi0";}
    if(mNeutronCalc){return "Neutron";}
    if(mLambda0Calc){return "Lambda0";}
    return "Non";
}

Int_t RHICfOptContainer::GetParticleRunIdx()
{
    if(mGammaCalc){return kGammaRun;}
    if(mPi0Calc){return kPi0Run;}
    if(mNeutronCalc){return kNeutronRun;}
    if(mLambda0Calc){return kLambda0Run;}
    return -1;
}

void RHICfOptContainer::InitOptions()
{
    mRunType = kALLRun;
    mInputDataPath = "/gpfs01/star/pwg/slee5/RHICfEvent"; // Default input path
    mInputDataList = "";
    mTablePath = "";
    mDataPath = "";
    mFigurePath = "";
    mExecuteEventNum = -1;
    mExecuteFileNum = -1;

    mGammaCalc = false;
    mPi0Calc = false;
    mNeutronCalc = false;
    mLambda0Calc = false;

    mForceCalculateMass = false;
    mForceCalculateBinning = false;
    mForceCalculateDilution = false;
    mForceCalculateAsymmetry = false;
    mForceCalculateSystematicError = false;
    
    mForceDefaultBinning = false;

    mOffDetPoint = false;
    mOffParticle = false;
    mOffTPCTrack = false;
    mOffBTof = false;
    mOffBBC = false;
    mOffVPD = false;
    mOffZDC = false;
    mOffFMS = false;
    mOffRPS = false;
}

void RHICfOptContainer::FindDirPath()
{
    // Find a directory
    TString currentPath = gSystem -> pwd();
    TObjArray *tokens = currentPath.Tokenize("/");
    TString directory = "";

    TList *listOfDirs;
    TObject *objDir;
    for(int i=0; i<tokens->GetEntries()-3; i++){
        TString dirPath = "";
        for(int j=0; j<tokens->GetEntries()-i; j++){
            dirPath = dirPath+"/"+ ((TObjString *) tokens -> At(j)) -> GetString();
        }
        TSystemDirectory dir("dir", dirPath);
        listOfDirs = dir.GetListOfFiles();
        TIter next(listOfDirs);

        while((objDir = next())){
            TSystemFile* dirPtr = dynamic_cast<TSystemFile*>(objDir);
            if(dirPtr && dirPtr->IsDirectory()){
                TString dirName = dirPtr->GetName();
                if(dirName.Index("tables") != -1){
                    mTablePath = dirPath+"/tables";
                    cout << "RHICfOptContainer::FindDirPath() -- Tables dir was found. " << mTablePath << endl;
                }
                if(dirName.Index("data") != -1){
                    mDataPath = dirPath+"/data";
                    cout << "RHICfOptContainer::FindDirPath() -- Data dir was found. " << mDataPath << endl;
                }
                if(dirName.Index("figure") != -1){
                    mFigurePath = dirPath+"/figure";
                    cout << "RHICfOptContainer::FindDirPath() -- Figure dir was found. " << mFigurePath << endl;
                }
            }
        }
    }
}

void RHICfOptContainer::PrintOptions()
{
    cout << "== RHICfOptContainer::PrintOptions()" << endl;
    cout << "== Calculate Run Type: " << GetRunTypeName(GetRunType()) << endl;
    cout << "== Calculate Particle Type: " << GetParticleRunName() << endl;
    cout << "== Input Data Path: " << GetInputDataPath() << endl;
    cout << "== Input Data List: " << ((GetInputDataList() == "")? "Autometic : " : GetInputDataList()) << endl; 
    cout << "== Table Data Path: " << GetTablePath() << endl;
    cout << "== Output Data Path: " << GetDataPath() << endl;
    cout << "== Figure Data Path: " << GetFigurePath() << endl;
    cout << "== Execute Event Number: " << GetExecuteEventNum() << endl;
    cout << "== Execute Input File Number: " << GetExecuteFileNum() << endl;
    cout << "== Force to Calculate Mass: " << ((IsForceCalculateMass() == true)? "TRUE" : "FALSE") << endl;
    cout << "== Force to Calculate Binning: " << ((IsForceCalculateBinning() == true)? "TRUE" : "FALSE") << endl;
    cout << "== Detector data ON-OFF:" << endl;
    cout << "     RHICfDetPoint : " << ((GetOffDetPoint() == false)? "ON" : "OFF") << endl;
    cout << "     RHICfParticle : " << ((GetOffParticle() == false)? "ON" : "OFF") << endl;
    cout << "     TPCTrack      : " << ((GetOffTPCTrack() == false)? "ON" : "OFF") << endl;
    cout << "     B-TOF         : " << ((GetOffBTof() == false)? "ON" : "OFF") << endl;
    cout << "     BBC           : " << ((GetOffBBC() == false)? "ON" : "OFF") << endl;
    cout << "     VPD           : " << ((GetOffVPD() == false)? "ON" : "OFF") << endl;
    cout << "     ZDC           : " << ((GetOffZDC() == false)? "ON" : "OFF") << endl;
    cout << "     FMS           : " << ((GetOffFMS() == false)? "ON" : "OFF") << endl;
    cout << "     RPS           : " << ((GetOffRPS() == false)? "ON" : "OFF") << endl;
}