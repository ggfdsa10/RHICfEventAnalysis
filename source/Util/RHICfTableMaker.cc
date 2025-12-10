#include "RHICfTableMaker.hh"

#include "TObjArray.h"
#include "TObjString.h"

RHICfTableMaker* RHICfTableMaker::mTableMaker = nullptr;
RHICfTableMaker* RHICfTableMaker::GetTableMaker()
{
    if (mTableMaker != nullptr){return mTableMaker;}
    return new RHICfTableMaker();
}

RHICfTableMaker::RHICfTableMaker() 
{
    mTableMaker = this;
}

RHICfTableMaker::~RHICfTableMaker()
{
}

int RHICfTableMaker::Init()
{
    mOptContainer = RHICfOptContainer::GetOptContainer();

    mTableName.clear();
    mTable.clear();

    return 1;
}

int RHICfTableMaker::InitTable(TString tableName)
{
    TString particleRunName = mOptContainer -> GetParticleRunName();
    TString tmpFileName = Form("%s_Results_%s", tableName.Data(), particleRunName.Data());
    TString runName = mOptContainer->GetRunTypeName(mOptContainer->GetRunType());
    TString FileName = tmpFileName + "_" + runName;
    
    int tableIdx = FindTableIdx(tableName);
    if(tableIdx == -1){
        vector<TableData> table;
        if(!ReadTable(tableName, table)){
            cout << "RHICfTableMaker::InitTable(" << tableName << ") -- Table is not found, new table will be assigned" << endl;
        }

        mTableName.push_back(FileName);
        mTable.push_back(table);

        if(FindTableIdx(tableName) != -1){return 1;}
    }
    
    return 0;
}

int RHICfTableMaker::SaveTable(TString tableName, vector<TableData> table)
{
    int tableIdx = FindTableIdx(tableName);
    if(tableIdx != -1){
        mTable[tableIdx] = table;
        MakeTable(mTableName[tableIdx], mTable[tableIdx]);
    }
    return 0;
}

double RHICfTableMaker::GetTableData(TString tableName, int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, int valueIdx)
{
    int tableIdx = FindTableIdx(tableName);
    if(tableIdx != -1){
        return FindTable(mTable[tableIdx], runIdx, typeIdx, dleIdx, ptIdx, xfIdx, valueIdx);
    }
    return -999.;
}

double RHICfTableMaker::GetTableData(TString tableName, int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, int valueIdx, int valueIdx2)
{
    int tableIdx = FindTableIdx(tableName);
    if(tableIdx != -1){
        return FindTable(mTable[tableIdx], runIdx, typeIdx, dleIdx, ptIdx, xfIdx, valueIdx, valueIdx2);
    }
    return -999.;
}

int RHICfTableMaker::FindTableIdx(TString tableName)
{
    TString particleRunName = mOptContainer -> GetParticleRunName();
    TString tmpFileName = Form("%s_Results_%s", tableName.Data(), particleRunName.Data());
    TString runName = mOptContainer->GetRunTypeName(mOptContainer->GetRunType());
    TString FileName = tmpFileName + "_" + runName;
    
    FileName.ToUpper();
    int tableListSize = mTableName.size();
    for(int i=0; i<tableListSize; i++){
        TString tmpName = mTableName[i];
        tmpName.ToUpper();
        if(tmpName.Index(FileName) != -1){return i;}
    }
    return -1;
}

double RHICfTableMaker::FindTable(vector<TableData> table, int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, int valueIdx)
{
    int tableNum = table.size();
    for(int i=0; i<tableNum; i++){
        if(table[i].runIdx != runIdx){continue;}
        if(table[i].typeIdx != typeIdx){continue;}
        if(table[i].dleIdx != dleIdx){continue;}
        if(table[i].ptIdx != ptIdx){continue;}
        if(table[i].xfIdx != xfIdx){continue;}

        int valueNum = table[i].values.size();
        if(valueIdx == -1){return double(valueNum);}
        for(int j=0; j<valueNum; j++){
            if(j != valueIdx){continue;}
            return table[i].values[j];
        }
    }
    return -999.;
}

double RHICfTableMaker::FindTable(vector<TableData> table, int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, int valueIdx, int valueIdx2)
{
    int tableNum = table.size();
    for(int i=0; i<tableNum; i++){
        if(table[i].runIdx != runIdx){continue;}
        if(table[i].typeIdx != typeIdx){continue;}
        if(table[i].dleIdx != dleIdx){continue;}
        if(table[i].ptIdx != ptIdx){continue;}
        if(table[i].xfIdx != xfIdx){continue;}

        if(table[i].values.size() == 1){break;}
        if(int(table[i].values[0]) != valueIdx){continue;}
        int valueNum = table[i].values.size();
        for(int j=0; j<valueNum; j++){
            if(j != valueIdx2+1){continue;}
            return table[i].values[j];
        }
    }
    return -999.;
}

int RHICfTableMaker::ReadTable(TString tableName, vector<TableData> &tableData)
{
    TString tablePath = mOptContainer -> GetTablePath();
    TString particleRunName = mOptContainer -> GetParticleRunName();
    TString fileName = Form("%s/%s_Results_%s", tablePath.Data(), tableName.Data(), particleRunName.Data());

    for(int run=0; run<kRunNum; run++){
        if(mOptContainer->GetRunType() != kALLRun && mOptContainer->GetRunType() != run){continue;}
        if(mTableStream.is_open()){mTableStream.close();}

        TString runName = mOptContainer->GetRunTypeName(run);
        TString runFileName = fileName + Form("_%s.dat", runName.Data());
        mTableStream.open(runFileName.Data(), ios::in);

        TObjArray *tokens;
        if(mTableStream.is_open()){
            std::string tmpLine;
            while(getline(mTableStream, tmpLine)){
                tokens = TString(tmpLine).Tokenize(" ");
                if(tokens->GetEntries() == 0){continue;}

                TableData data;
                data.runIdx = (((TObjString *) tokens -> At(0)) -> GetString()).Atoi();
                if(run != data.runIdx){continue;}

                data.typeIdx = (((TObjString *) tokens -> At(1)) -> GetString()).Atoi();
                data.dleIdx = (((TObjString *) tokens -> At(2)) -> GetString()).Atoi();
                data.ptIdx = (((TObjString *) tokens -> At(3)) -> GetString()).Atoi();
                data.xfIdx = (((TObjString *) tokens -> At(4)) -> GetString()).Atoi();

                for(int i=5; i<tokens->GetEntries(); i++){
                    double value = (((TObjString *) tokens -> At(i)) -> GetString()).Atof();
                    data.values.push_back(value);   
                }
                tableData.push_back(data);
            }
            mTableStream.close();
        }
        else{return 0;}

        cout << "RHICfTableMaker::ReadTable(" << tableName << ") -- Data table has been read, " << runFileName << endl;
    }
    return 1; 
}

int RHICfTableMaker::MakeTable(TString tableName, vector<TableData> &tableData)
{
    TString tablePath = mOptContainer -> GetTablePath();
    TString fileName = Form("%s/%s.dat", tablePath.Data(), tableName.Data());

    for(int run=0; run<kRunNum; run++){
        if(mOptContainer->GetRunType() != kALLRun && mOptContainer->GetRunType() != run){continue;}

        if(mTableStream.is_open()){mTableStream.close();}
        mTableStream.open(fileName.Data(), ios::out);

        if(mTableStream.is_open()){
            int tableNum = tableData.size();
            for(int i=0; i<tableNum; i++){
                if(tableData[i].runIdx != run){continue;}

                mTableStream << tableData[i].runIdx << " " << tableData[i].typeIdx << " " << tableData[i].dleIdx << " " << tableData[i].ptIdx << " " << tableData[i].xfIdx << " ";
                int valueNum = tableData[i].values.size();
                for(int j=0; j<valueNum; j++){
                    mTableStream << tableData[i].values[j] << " ";
                }
                mTableStream << endl;
            }
            mTableStream.close();
            cout << "RHICfTableMaker::MakeTable(" << tableName << ") -- Data table has been saved: " << endl;
        }
        else{return 0;}
    }
    return 1;
}
