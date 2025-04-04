// AN for RHICf particles calculation and analysis

#include "RHICfAnAnalysis.hh"

using namespace std;

int main(int argc, char** argv)
{
    RHICfAnAnalysis* anAnalysis = new RHICfAnAnalysis();
    anAnalysis -> SetRunType("TS");

    anAnalysis -> ForceCalculateMass();
    // anAnalysis -> ForceCalculateBinning(); 

    anAnalysis -> CalculatePi0();  
    anAnalysis -> SetExecuteEventNum(100000);    
    anAnalysis -> SetExecuteFileNum(70);    

    anAnalysis -> Init();
    anAnalysis -> Calculate();
    anAnalysis -> Finish();


    cout << "Terminated." << endl;
    return 0;
}


