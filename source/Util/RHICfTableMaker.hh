#ifndef RHICfTableMaker_hh
#define RHICfTableMaker_hh

#include <vector>
#include <iostream>
#include <fstream>

#include "RHICfOptContainer.hh"

using namespace std;

enum TableFlag
{
    kNotExist = 0,
    kExistTable = 1,
    kExistPartOfTable = 2
};

class RHICfTableMaker
{
    public:
        struct TableData
        {
            int runIdx;
            int typeIdx;
            int dleIdx;
            int ptIdx;
            int xfIdx;
            vector<double> values;
        };

        static RHICfTableMaker* GetTableMaker();

        RHICfTableMaker();
        ~RHICfTableMaker();

        virtual int Init();
        int InitTable(TString tableName);
        int SaveTable(TString tableName, vector<TableData> table);
        double GetTableData(TString tableName, int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, int valueIdx);

    private:
        int FindTableIdx(TString tableName);
        double FindTable(vector<TableData> table, int runIdx, int typeIdx, int dleIdx, int ptIdx, int xfIdx, int valueIdx);
        int ReadTable(TString tableName, vector<TableData> &tableData);
        int MakeTable(TString tableName, vector<TableData> &tableData);

        static RHICfTableMaker* mTableMaker;
        RHICfOptContainer* mOptContainer;

        std::fstream mTableStream;
        vector<TString> mTableName;
        vector<vector<TableData> > mTable;

};

#endif
