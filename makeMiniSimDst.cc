// RHICf simulation calculation and analysis
#include "RHICfSimAnalysis.hh"
#include "RHICfOptContainer.hh"

using namespace std;

int main(int argc, char** argv)
{
    RHICfOptContainer* optContainer = RHICfOptContainer::GetOptContainer();
    optContainer -> SetSimInputDataPath("/gpfs01/star/pwg/slee5/simulation");
    optContainer -> SetRunType("TS");
    // optContainer -> SetSimModel(rQGSJETIII);
    // optContainer -> SetSimModel(rSIBYLL);
    optContainer -> SetSimModel(rPythia8);
    // optContainer -> SetSimModel(rEPOSLHCR_S);
    optContainer -> CalculatePi0();  
    optContainer->SetTag("2");

    RHICfSimDstReader* reader = new RHICfSimDstReader();
    reader -> Init();
    reader -> Make();

    cout << "Terminated." << endl;
    
    return 0;
}
