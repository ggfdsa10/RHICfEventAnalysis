// AN for RHICf particles calculation and analysis
#include "RHICfAnAnalysis.hh"
#include "RHICfAsymmetry.hh"

using namespace std;

int main(int argc, char** argv)
{
    // RHICfAnAnalysis* anAnalysis = new RHICfAnAnalysis();
    // anAnalysis -> SetRunType("TOP");

    // // anAnalysis -> ForceCalculateMass();
    // // anAnalysis -> ForceCalculateBinning(); 
    // // anAnalysis -> ForceCalculatePolarization(); 
    // // anAnalysis -> ForceCalculateSystematicError(); 
    // // anAnalysis -> ForceDefaultBinning(); 

    // // ====== beam center method option ======
    // // SetBeamCenterMethod( method , refFill for TOPrun)
    // // method : kBeamCenterHit, kBeamCenterScan
    // // ref : kBeamCenterRef21148, kBeamCenterRef21150

    // // anAnalysis -> SetConditionName(""); 
    // anAnalysis -> SetBeamCenterMethod(kBeamCenterHit, kBeamCenterRef21150);
    // anAnalysis -> CalculatePi0();  
    // // anAnalysis -> SetExecuteFileNum();    

    // anAnalysis -> Init();
    // anAnalysis -> Calculate();
    // anAnalysis -> Finish();

    // =========================================== 
    RHICfAsymmetry* asym = new RHICfAsymmetry();
    asym -> Init();
    asym -> Calculate();
    

    cout << "Terminated." << endl;
    
    return 0;
}

