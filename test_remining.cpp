#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <time.h>
#include <cmath>
#include <cstdlib>
#include <functional> //vector<reference_wrapper>
#include <queue>
#include <list>
#include <windows.h>
#include <psapi.h>
#include <iostream>
using namespace std;

class SeqData
{
public:
    int sid;
    vector<int> IndexArray;
    vector<int> TidArray;
    vector<int> ItemArray;
    vector<int> IuArray; // new
    vector<double> UtilityArray;
    vector<double> RuArray;
    map<int, vector<int>> ItemIdxTable; //<Item,[Idx set]>
    vector<int> IdxTid;                 // size n代表Tid n、<Tid起始Idx>
};
vector<SeqData> VecDataBase;

class L1_UtilityInfo
{
public:
    vector<int> VecIndex;
    vector<int> VecTid;
    vector<int> VecIu; // 其實可查表?
    vector<double> VecUtility;
    double CaseUtility = 0.0;
};

class L2_SequenceInfo
{
public:
    int sid;
    double SeqPEU = 0;
    double SeqUt = 0;
    int SeqUtCase = 0;
    int SeqPEUCase = 0;
    int OverMinUtilCase = 0;
    int IdxOfLastLevelSeq = 0;
    vector<L1_UtilityInfo> L1_UtInfo;
    map<int, vector<pair<int, int>>> IdxPos;
};

class L3_NodeInfo
{
public:
    string pattern;
    vector<int> VecPattern;
    double SumPEU = 0;
    double SumUt = 0;
    vector<L2_SequenceInfo> L2_SeqInfo;
};

vector<L3_NodeInfo> Node_SingleItem;
map<int, int> IdxNodeSingleItem;

vector<int> ItemQueue; // 紀錄Externalfile中讀取Item的先後順序，來放入VecDataBase中
map<int, double> ExternalUt;
string FileName_Eu = "";
string FileName_DB = "";

vector<reference_wrapper<L3_NodeInfo>> PatternPath;
// vector<L3_NodeInfo *> pp;
vector<L3_NodeInfo> HUSPattern; // For remining check

int Single_ItemCounter = 0;
int I_ExtensionCounter = 0;
int S_ExtensionCounter = 0;
int TotalPatternCounter = 0;
double SumSWU = 0;
double MinUtil = 0;
double memoryMB = 0;
double maxut = 0;

void Cout_2DVec(vector<vector<int>> vector2D)
{
    for (int i = 0; i < vector2D.size(); i++)
    {
        for (int j = 0; j < vector2D[i].size(); j++)
        {
            cout << vector2D[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}
void Cout_ExternalUt(map<int, double> ExternalUt)
{
    for (auto a : ExternalUt)
    {
        cout << a.first << "," << a.second << endl;
    }
}
void Cout_VecDB(vector<SeqData> VecDataBase)
{
    for (int i = 1; i < VecDataBase.size(); i++) // 0預設填空
    {
        cout << "Sid : " << VecDataBase[i].sid << endl;
        cout << "  IndexArray : ";
        for (int j = 0; j < VecDataBase[i].IndexArray.size(); j++)
        {
            cout << VecDataBase[i].IndexArray[j] << " ";
        }
        cout << endl;
        cout << " **ItemArray : ";
        for (int j = 0; j < VecDataBase[i].ItemArray.size(); j++)
        {
            cout << VecDataBase[i].ItemArray[j] << " ";
        }
        cout << endl;
        cout << "    TidArray : ";
        for (int j = 0; j < VecDataBase[i].TidArray.size(); j++)
        {
            cout << VecDataBase[i].TidArray[j] << " ";
        }
        cout << endl;
        cout << "     IuArray : ";
        for (int j = 0; j < VecDataBase[i].IuArray.size(); j++)
        {
            cout << VecDataBase[i].IuArray[j] << " ";
        }
        cout << endl;
        cout << "UtilityArray : ";
        for (int j = 0; j < VecDataBase[i].UtilityArray.size(); j++)
        {
            cout << VecDataBase[i].UtilityArray[j] << " ";
        }
        cout << endl;
        cout << "     RuArray : ";
        for (int j = 0; j < VecDataBase[i].RuArray.size(); j++)
        {
            cout << VecDataBase[i].RuArray[j] << " ";
        }
        cout << endl;
        cout << "ItemIdxTable : ";
        for (auto a : VecDataBase[i].ItemIdxTable)
        {
            cout << a.first << "[";
            for (int j = 0; j < a.second.size(); j++)
            {
                cout << a.second[j] << " ";
            }
            cout << "]  ";
        }
        cout << endl;
        cout << "      IdxTid : ";
        for (int j = 0; j < VecDataBase[i].IdxTid.size(); j++)
        {
            cout << VecDataBase[i].IdxTid[j] << " ";
        }
        cout << endl;
    }
    cout << endl;
}

void Read_ExternalUt(string ExternalFileName)
{
    string StrLine;
    ifstream read_file(ExternalFileName);
    if (!read_file.is_open())
        cout << "*** Read External Utility File FAIL. ***" << endl;
    else
    {
        while (getline(read_file, StrLine))
        {
            for (int i = 0; i < StrLine.length(); i++)
            {
                if (StrLine[i] == ':')
                {
                    string Item = StrLine.substr(0, i);
                    double Eu = stod(StrLine.substr(i + 1, StrLine.length()));

                    ExternalUt.insert(make_pair(stoi(Item), Eu));
                }
            }
        }
    }
}
void Read_Database(string DatabaseFileName)
{
    string strline;
    ifstream ifs_readDB(DatabaseFileName);
    if (!ifs_readDB.is_open())
        cout << " Read Database FAIL. " << endl;
    else
    {
        int seq = 0;
        while (getline(ifs_readDB, strline))
        {
            seq++;
            SeqData SD;
            SD.sid = seq;
            SD.IdxTid.push_back(0); // 沒有Tid 0，先填0
            SD.IdxTid.push_back(1); // Tid 1必從Index 1開始，先填1

            int IndexOfSpace = 0;
            int TidCounter = 1;
            int Index = 1;
            int ItemForUtArray = 0;
            double SU = 0;
            for (int i = strline.size() - 1; i >= 0; i--)
            {
                if (strline[i] == ':')
                {
                    SU = stod((strline.substr(i + 1, strline.size() - i + 1)));
                    SumSWU += SU;
                    break;
                }
            }
            for (int i = 0; i < strline.length(); i++)
            {
                int IndexOfLeftBrackets;
                int IndexOfRightBrackets;
                if (strline[i] == '[')
                {
                    IndexOfLeftBrackets = i;
                    if (IndexOfSpace == 0) // 開頭第一個item
                    {
                        SD.ItemArray.push_back(stoi(strline.substr(0, IndexOfLeftBrackets)));
                        SD.ItemIdxTable[stoi(strline.substr(0, IndexOfLeftBrackets))].push_back(Index);

                        SD.IndexArray.push_back(Index);
                        Index++;

                        ItemForUtArray = stoi(strline.substr(0, IndexOfLeftBrackets));
                    }
                    else
                    {
                        SD.ItemArray.push_back(stoi(strline.substr(IndexOfSpace + 1, IndexOfLeftBrackets)));
                        SD.ItemIdxTable[stoi(strline.substr(IndexOfSpace + 1, IndexOfLeftBrackets))].push_back(Index);

                        SD.IndexArray.push_back(Index);
                        Index++;

                        ItemForUtArray = stoi(strline.substr(IndexOfSpace + 1, IndexOfLeftBrackets));
                    }
                }
                if (strline[i] == ']')
                {
                    IndexOfRightBrackets = i;
                    SD.IuArray.push_back(stod(strline.substr(IndexOfLeftBrackets + 1, IndexOfRightBrackets - IndexOfLeftBrackets - 1)));

                    SD.TidArray.push_back(TidCounter);

                    for (auto a : ExternalUt)
                    {
                        if (a.first == ItemForUtArray)
                        {
                            SD.UtilityArray.push_back(a.second * stod(strline.substr(IndexOfLeftBrackets + 1, IndexOfRightBrackets - IndexOfLeftBrackets - 1)));
                            // 避免浮點數表示成3.10862e-015等數值(實則3.10862-0.15)
                            if (SU - (a.second * stod(strline.substr(IndexOfLeftBrackets + 1, IndexOfRightBrackets - IndexOfLeftBrackets - 1))) < 0.001)
                            {
                                SD.RuArray.push_back(0);
                            }
                            else
                            {
                                SD.RuArray.push_back(SU - (a.second * stod(strline.substr(IndexOfLeftBrackets + 1, IndexOfRightBrackets - IndexOfLeftBrackets - 1))));
                            }
                            SU -= a.second * stod(strline.substr(IndexOfLeftBrackets + 1, IndexOfRightBrackets - IndexOfLeftBrackets - 1));
                        }
                    }
                }
                if (strline[i] == ' ')
                {
                    IndexOfSpace = i;
                }
                if (strline[i] == '-')
                {
                    if (strline[i + 1] == '1' && strline[i + 3] != '-')
                    {
                        SD.IdxTid.push_back(Index);
                    }
                    TidCounter++;
                }
            }
            VecDataBase.push_back(SD);
        }
    }
}

void BulidSingleItems(vector<SeqData> VecDataBase)
{
    // Node_SingleItem的順序 = ExternalUt的順序 = IdxNodeSingleItem的順序
    // 藉由IdxNodeSingleItem來更新SingleItem的資訊
    int n = 1;
    for (auto Item : ExternalUt) // 開好所有SingleItem的L3 class
    {
        if (Item.first == 0)
        {
            continue;
        }
        L3_NodeInfo L3_node;
        L3_node.pattern = to_string(Item.first);
        L3_node.VecPattern.push_back(Item.first);
        Node_SingleItem.push_back(L3_node); // Node_SingleItem經map排過序

        IdxNodeSingleItem.insert(make_pair(Item.first, n));
        n++;
    }

    for (int i = 1; i < VecDataBase.size(); i++)
    {
        for (int j = 0; j < VecDataBase[i].IndexArray.size(); j++)
        {
            // 第一個S2 or 不同sid的新S2
            if (!Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() ||
                Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo[Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() - 1].sid != VecDataBase[i].sid)
            {
                L2_SequenceInfo L2_seq;
                L2_seq.sid = VecDataBase[i].sid;
                L2_seq.SeqUtCase = VecDataBase[i].IndexArray[j];
                Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.push_back(L2_seq);

                L1_UtilityInfo L1_ut;
                L1_ut.VecIndex.push_back(VecDataBase[i].IndexArray[j]);
                L1_ut.VecTid.push_back(VecDataBase[i].TidArray[j]);
                L1_ut.VecIu.push_back(VecDataBase[i].IuArray[j]);
                L1_ut.VecUtility.push_back(VecDataBase[i].UtilityArray[j]);
                L1_ut.CaseUtility = VecDataBase[i].UtilityArray[j];

                Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo[Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() - 1].L1_UtInfo.push_back(L1_ut);
                Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo[Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() - 1].SeqUt = VecDataBase[i].UtilityArray[j];

                Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].SumUt += VecDataBase[i].UtilityArray[j];

                if (VecDataBase[i].RuArray[j] > 0)
                {
                    Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo[Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() - 1].SeqPEU = VecDataBase[i].UtilityArray[j] + VecDataBase[i].RuArray[j];
                    Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].SumPEU += VecDataBase[i].UtilityArray[j] + VecDataBase[i].RuArray[j];
                    Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo[Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() - 1].SeqPEUCase = VecDataBase[i].IndexArray[j];
                }

                if (L1_ut.CaseUtility > MinUtil)
                {
                    Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo[Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() - 1].OverMinUtilCase++;
                }
                Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo[Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() - 1].IdxPos[VecDataBase[i].IndexArray[j]].push_back({Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo[Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() - 1].L1_UtInfo.size() - 1, 0});
            }
            else // 已經存在任意L2，且sid已出現過
            {
                L1_UtilityInfo L1_ut;
                L1_ut.VecIndex.push_back(VecDataBase[i].IndexArray[j]);
                L1_ut.VecTid.push_back(VecDataBase[i].TidArray[j]);
                L1_ut.VecIu.push_back(VecDataBase[i].IuArray[j]);
                L1_ut.VecUtility.push_back(VecDataBase[i].UtilityArray[j]);
                L1_ut.CaseUtility = VecDataBase[i].UtilityArray[j];

                Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo[Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() - 1].L1_UtInfo.push_back(L1_ut); // SeqPEU、SumPEU

                // SeqPEU、SumPEU
                if (Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo[Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() - 1].SeqPEU <
                        (VecDataBase[i].UtilityArray[j] + VecDataBase[i].RuArray[j]) &&
                    VecDataBase[i].RuArray[j] > 0)
                {
                    Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].SumPEU -= Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo[Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() - 1].SeqPEU; //(扣除)取出之前的PEU
                    Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].SumPEU += VecDataBase[i].UtilityArray[j] + VecDataBase[i].RuArray[j];                                                                                                               // 改成加現在這個更新後的PEU
                    Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo[Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() - 1].SeqPEU = VecDataBase[i].UtilityArray[j] + VecDataBase[i].RuArray[j];

                    Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo[Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() - 1].SeqPEUCase = VecDataBase[i].IndexArray[j];
                }

                // SeqUt、SumUt
                if (Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo[Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() - 1].SeqUt <
                    VecDataBase[i].UtilityArray[j])
                {
                    Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].SumUt -= Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo[Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() - 1].SeqUt; //(扣除)取出之前的Ut
                    Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].SumUt += VecDataBase[i].UtilityArray[j];                                                                                                                                          // 改成加現在這個更新後的Ut
                    Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo[Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() - 1].SeqUt = VecDataBase[i].UtilityArray[j];

                    Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo[Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() - 1].SeqUtCase = VecDataBase[i].IndexArray[j];
                }

                if (L1_ut.CaseUtility > MinUtil)
                {
                    Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo[Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() - 1].OverMinUtilCase++;
                }
                Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo[Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() - 1].IdxPos[VecDataBase[i].IndexArray[j]].push_back({Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo[Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() - 1].L1_UtInfo.size() - 1, 0});
            }
        }
    }
}

void Cout_HUSPL3(L3_NodeInfo L3_node)
{
    cout << "---------- Pattern:" << L3_node.pattern << " ----------" << endl;
    /*
    for (int i = 0; i < L3_node.VecPattern.size(); i++)
    {
        cout << L3_node.VecPattern[i] << " ";
    }
    cout << endl;
    */
    cout << "=== MinUtil:" << MinUtil << " ===" << endl;
    cout << "** SumUt:" << L3_node.SumUt << " **" << endl;
    cout << "** SumPEU:" << L3_node.SumPEU << " **" << endl;

    for (int j = 0; j < L3_node.L2_SeqInfo.size(); j++)
    {
        cout << "Seq " << L3_node.L2_SeqInfo[j].sid << endl;
        cout << "SeqUt:" << L3_node.L2_SeqInfo[j].SeqUt << endl;
        cout << "SeqPEU:" << L3_node.L2_SeqInfo[j].SeqPEU << endl;
        cout << "SeqUtCase:" << L3_node.L2_SeqInfo[j].SeqUtCase << endl;
        cout << "SeqPEUCase:" << L3_node.L2_SeqInfo[j].SeqPEUCase << endl;
        cout << "IdxOfLastLevelSeq:" << L3_node.L2_SeqInfo[j].IdxOfLastLevelSeq << endl;
        cout << "|Idx Iu Ut|InstanceUt" << endl;
        // cout << "Over Minutil:" << L3_node.L2_SeqInfo[j].OverMinUtilCase << endl;
        for (int k = 0; k < L3_node.L2_SeqInfo[j].L1_UtInfo.size(); k++)
        {
            cout << "[";
            for (int h = 0; h < L3_node.L2_SeqInfo[j].L1_UtInfo[k].VecIndex.size(); h++)
            {
                cout << L3_node.L2_SeqInfo[j].L1_UtInfo[k].VecIndex[h] << " ";
            }
            cout << "]";
            cout << "[";
            for (int h = 0; h < L3_node.L2_SeqInfo[j].L1_UtInfo[k].VecIndex.size(); h++)
            {
                cout << L3_node.L2_SeqInfo[j].L1_UtInfo[k].VecIu[h] << " ";
            }
            cout << "]";
            cout << "[";
            for (int h = 0; h < L3_node.L2_SeqInfo[j].L1_UtInfo[k].VecIndex.size(); h++)
            {
                cout << L3_node.L2_SeqInfo[j].L1_UtInfo[k].VecUtility[h] << " ";
            }
            cout << "]";
            cout << L3_node.L2_SeqInfo[j].L1_UtInfo[k].CaseUtility << endl;
        }
        for (auto a : L3_node.L2_SeqInfo[j].IdxPos)
        {
            cout << "Idx:" << a.first << "->";
            for (auto pair : a.second)
            {
                cout << "(" << pair.first << "," << pair.second << ")";
            }
            cout << endl;
        }
        cout << endl;
    }
    cout << endl;
}

void DeleteLowRSUFromlist(L3_NodeInfo NodeUC, set<int> &ilist, set<int> &slist)
{
    map<int, double> S_RSUmap; //<Item,RSU>
    map<int, double> I_RSUmap; //<Item,RSU>
    for (int i = 0; i < NodeUC.L2_SeqInfo.size(); i++)
    {
        // slist
        // cout << NodeUC.L2_SeqInfo[i].L1_UtInfo[0].VecIndex.back() << endl;
        // 最後一個Tid的Item不往後掃
        if (VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[NodeUC.L2_SeqInfo[i].L1_UtInfo[0].VecIndex.back() - 1] < VecDataBase[NodeUC.L2_SeqInfo[i].sid].IdxTid.back())
        {
            // cout << VecDataBase[NodeUC.L2_SeqInfo[i].sid].IdxTid[VecDataBase[NodeUC.L2_SeqInfo[i].sid].TidArray[NodeUC.L2_SeqInfo[i].L1_UtInfo[0].VecIndex.back() - 1] + 1] << " ";
            for (int j = VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[VecDataBase[NodeUC.L2_SeqInfo[i].sid].IdxTid[VecDataBase[NodeUC.L2_SeqInfo[i].sid].TidArray[NodeUC.L2_SeqInfo[i].L1_UtInfo[0].VecIndex.back() - 1] + 1] - 1] - 1; j < VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray.size(); j++)
            {
                // cout << VecDataBase[NodeUC.L2_SeqInfo[i].sid].ItemArray[j] << " ";
                S_RSUmap[VecDataBase[NodeUC.L2_SeqInfo[i].sid].ItemArray[j]] += NodeUC.L2_SeqInfo[i].SeqPEU;
                if (S_RSUmap[VecDataBase[NodeUC.L2_SeqInfo[i].sid].ItemArray[j]] >= MinUtil)
                {
                    slist.insert(VecDataBase[NodeUC.L2_SeqInfo[i].sid].ItemArray[j]);
                }
            }
        }

        // ilist
        // cout << "seq " << NodeUC.L2_SeqInfo[i].sid << endl;
        for (int j = 0; j < NodeUC.L2_SeqInfo[i].L1_UtInfo.size(); j++)
        {
            // cout << NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex.back() << endl;
            //  cout << VecDataBase[NodeUC.L2_SeqInfo[i].sid].TidArray[NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex.back() - 1] << endl;
            for (int k = VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex.back()] - 1; k < VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray.size(); k++)
            {
                if (VecDataBase[NodeUC.L2_SeqInfo[i].sid].TidArray[NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex.back() - 1] != VecDataBase[NodeUC.L2_SeqInfo[i].sid].TidArray[k])
                {
                    break;
                }
                // cout << VecDataBase[NodeUC.L2_SeqInfo[i].sid].ItemArray[k] << " ";
                I_RSUmap[VecDataBase[NodeUC.L2_SeqInfo[i].sid].ItemArray[k]] += NodeUC.L2_SeqInfo[i].SeqPEU;
                if (I_RSUmap[VecDataBase[NodeUC.L2_SeqInfo[i].sid].ItemArray[k]] >= MinUtil)
                {
                    ilist.insert(VecDataBase[NodeUC.L2_SeqInfo[i].sid].ItemArray[k]);
                }
            }
        }
    }
}

void HUSP(L3_NodeInfo &NodeUC)
{
    if (NodeUC.SumPEU < MinUtil)
    {
        return;
    }
    set<int> ilist;
    set<int> slist;
    DeleteLowRSUFromlist(NodeUC, ilist, slist);

    for (int Iitem : ilist)
    {
        L3_NodeInfo NIF;
        NIF.SumPEU = 0;
        NIF.SumUt = 0;
        NIF.pattern = NodeUC.pattern + "-" + to_string(Iitem);
        NIF.VecPattern = NodeUC.VecPattern;
        NIF.VecPattern.push_back(Iitem);

        for (int i = 0; i < NodeUC.L2_SeqInfo.size(); i++)
        {
            for (int ItemIdx : VecDataBase[NodeUC.L2_SeqInfo[i].sid].ItemIdxTable[Iitem])
            {
                for (int j = 0; j < NodeUC.L2_SeqInfo[i].L1_UtInfo.size(); j++)
                {
                    if (ItemIdx > NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex.back() && VecDataBase[NodeUC.L2_SeqInfo[i].sid].TidArray[ItemIdx - 1] == NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecTid.back())
                    {
                        if (NIF.L2_SeqInfo.size() == 0 || NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].sid != VecDataBase[NodeUC.L2_SeqInfo[i].sid].sid)
                        {
                            L2_SequenceInfo SIF;
                            SIF.sid = VecDataBase[NodeUC.L2_SeqInfo[i].sid].sid;
                            SIF.IdxOfLastLevelSeq = i;
                            NIF.L2_SeqInfo.push_back(SIF);

                            L1_UtilityInfo UIF;
                            // Index
                            UIF.VecIndex = NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex;
                            UIF.VecIndex.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[ItemIdx - 1]);
                            // Tid
                            UIF.VecTid = NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecTid;
                            UIF.VecTid.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].TidArray[ItemIdx - 1]);

                            // Iu、Ru、Utility(讀VecDB更新資訊)
                            for (int k = 0; k < NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex.size(); k++)
                            {
                                UIF.VecIu.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].IuArray[NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex[k] - 1]);
                                UIF.VecUtility.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].UtilityArray[NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex[k] - 1]);
                            }
                            UIF.VecIu.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].IuArray[ItemIdx - 1]);
                            UIF.VecUtility.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].UtilityArray[ItemIdx - 1]);

                            // CaseUtility(拿更新後的VecUtility來加總)
                            for (int k = 0; k < UIF.VecUtility.size(); k++)
                            {
                                UIF.CaseUtility += UIF.VecUtility[k];                        
                            }

                            // SeqUt、SumUt
                            NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt = UIF.CaseUtility;
                            NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUtCase = UIF.VecIndex.back();
                            NIF.SumUt += UIF.CaseUtility;

                            // SeqPEU、SumPEU
                            if (VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1] != 0)
                            {
                                NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU = UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1];
                                NIF.SumPEU += NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU;
                                NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEUCase = VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[ItemIdx - 1];
                            }
                            else
                            {
                                NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU = 0;
                            }

                            // Record Cases of Over MinUtil.
                            if (UIF.CaseUtility > MinUtil)
                            {
                                NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].OverMinUtilCase++;
                            }

                            NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.push_back(UIF);
                            UIF.CaseUtility = 0;
                        }
                        else // 已存在任何L2
                        {
                            L1_UtilityInfo UIF;
                            // Index
                            UIF.VecIndex = NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex;
                            UIF.VecIndex.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[ItemIdx - 1]);
                            // Tid
                            UIF.VecTid = NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecTid;
                            UIF.VecTid.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].TidArray[ItemIdx - 1]);

                            // Iu、Ru、Utility(讀VecDB更新資訊)
                            for (int k = 0; k < NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex.size(); k++)
                            {
                                UIF.VecIu.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].IuArray[NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex[k] - 1]);
                                UIF.VecUtility.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].UtilityArray[NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex[k] - 1]);
                            }
                            UIF.VecIu.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].IuArray[ItemIdx - 1]);
                            UIF.VecUtility.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].UtilityArray[ItemIdx - 1]);

                            // CaseUtility(拿更新後的VecUtility來加總)
                            for (int k = 0; k < UIF.VecUtility.size(); k++)
                            {
                                UIF.CaseUtility += UIF.VecUtility[k];
                            }

                            // 投影點的Index不同
                            if (NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo[NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.size() - 1].VecIndex.back() != VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[ItemIdx - 1])
                            {
                                // acu與SeqUt比較出SeqUt
                                if (UIF.CaseUtility > NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt)
                                {
                                    NIF.SumUt += (UIF.CaseUtility - NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt);                                
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt = UIF.CaseUtility;
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUtCase = UIF.VecIndex.back();
                                }

                                // ru!=0 且 acu+ru比現存的SeqPEU大，比較出SeqPEU
                                if (VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1] != 0 &&
                                    UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1] > NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU)
                                {
                                    NIF.SumPEU += ((UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1]) - NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU);
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU = UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1];
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEUCase = VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[ItemIdx - 1];
                                }

                                // Record Cases of Over MinUtil.
                                if (UIF.CaseUtility > MinUtil)
                                {
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].OverMinUtilCase++;
                                }
                                NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.push_back(UIF);
                                UIF.CaseUtility = 0;
                            }
                            else // 投影點Index相同
                            {
                                if (UIF.CaseUtility > NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo[NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.size() - 1].CaseUtility)
                                {
                                    // Record Cases of Over MinUtil.
                                    if (UIF.CaseUtility > MinUtil && NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo[NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.size() - 1].CaseUtility <= MinUtil)
                                    {
                                        NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].OverMinUtilCase++;
                                    }
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.pop_back();
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.push_back(UIF);
                                }

                                // SeqUt、SumUt
                                if (UIF.CaseUtility > NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt)
                                {
                                    NIF.SumUt += (UIF.CaseUtility - NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt);            
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt = UIF.CaseUtility;
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUtCase = UIF.VecIndex.back();
                                }
                                // SeqPEU、SumPEU
                                if (VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1] != 0 &&
                                    UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1] > NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU)
                                {
                                    NIF.SumPEU += ((UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1]) - NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU);
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU = UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1];
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEUCase = VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[ItemIdx - 1];
                                }
                                UIF.CaseUtility = 0;
                            }
                        }
                    }
                }
            }
        }

        // IdxPos
        for (int a = 0; a < NIF.L2_SeqInfo.size(); a++)
        {
            for (int b = 0; b < NIF.L2_SeqInfo[a].L1_UtInfo.size(); b++)
            {
                for (int c = 0; c < NIF.L2_SeqInfo[a].L1_UtInfo[b].VecIndex.size(); c++)
                {
                    NIF.L2_SeqInfo[a].IdxPos[NIF.L2_SeqInfo[a].L1_UtInfo[b].VecIndex[c]].push_back({b, c});
                }
            }
        }

        if (NIF.SumUt >= MinUtil)
        {
            I_ExtensionCounter++;
            //cout << NIF.pattern << " | " << NIF.SumUt << endl;
            HUSPattern.push_back(NIF);
            //Cout_HUSPL3(NIF);
        }
        if (NIF.SumUt > maxut)
        {            
            maxut = NIF.SumUt;
        }

        //Cout_HUSPL3(NIF);
        HUSP(NIF);
    }

    // S-Extension

    for (int Sitem : slist)
    {
        L3_NodeInfo NIF;
        NIF.SumPEU = 0;
        NIF.SumUt = 0;
        NIF.pattern = NodeUC.pattern + "," + to_string(Sitem);
        NIF.VecPattern = NodeUC.VecPattern;
        NIF.VecPattern.push_back(Sitem);

        for (int i = 0; i < NodeUC.L2_SeqInfo.size(); i++)
        {
            for (int ItemIdx : VecDataBase[NodeUC.L2_SeqInfo[i].sid].ItemIdxTable[Sitem])
            {
                for (int j = 0; j < NodeUC.L2_SeqInfo[i].L1_UtInfo.size(); j++)
                {
                    if (ItemIdx > NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex.back() &&
                        VecDataBase[NodeUC.L2_SeqInfo[i].sid].TidArray[ItemIdx - 1] > VecDataBase[NodeUC.L2_SeqInfo[i].sid].TidArray[NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex.back() - 1])
                    {
                        // 首個L2
                        if (NIF.L2_SeqInfo.size() == 0 || NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].sid != VecDataBase[NodeUC.L2_SeqInfo[i].sid].sid)
                        {
                            L2_SequenceInfo SIF;
                            SIF.sid = VecDataBase[NodeUC.L2_SeqInfo[i].sid].sid;
                            SIF.IdxOfLastLevelSeq = i;
                            NIF.L2_SeqInfo.push_back(SIF);

                            // Index、Tid
                            L1_UtilityInfo UIF;
                            // Index
                            UIF.VecIndex = NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex;
                            UIF.VecIndex.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[ItemIdx - 1]);
                            // IdxPos
                            for (int k = 0; k < UIF.VecIndex.size(); k++)
                            {
                                // NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].IdxPos[UIF.VecIndex[k]].push_back({NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.size(), k});
                            }
                            // Tid
                            UIF.VecTid = NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecTid;
                            UIF.VecTid.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].TidArray[ItemIdx - 1]);

                            // Iu、Utility(讀VecDB更新資訊)
                            for (int k = 0; k < NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex.size(); k++)
                            {
                                UIF.VecIu.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].IuArray[NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex[k] - 1]);
                                UIF.VecUtility.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].UtilityArray[NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex[k] - 1]);
                            }
                            UIF.VecIu.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].IuArray[ItemIdx - 1]);
                            UIF.VecUtility.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].UtilityArray[ItemIdx - 1]);

                            // CaseUtility(拿更新後的VecUtility來加總)
                            for (int k = 0; k < UIF.VecUtility.size(); k++)
                            {
                                UIF.CaseUtility += UIF.VecUtility[k];
                            }

                            // SeqUt、SumUt
                            NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt = UIF.CaseUtility;
                            NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUtCase = UIF.VecIndex.back();
                            NIF.SumUt += UIF.CaseUtility;

                            // SeqPEU、SumPEU
                            if (VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1] != 0)
                            {
                                NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU = UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1];
                                NIF.SumPEU += NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU;
                                NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEUCase = VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[ItemIdx - 1];
                            }
                            else
                            {
                                NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU = 0;
                            }

                            // Record Cases of Over MinUtil.
                            if (UIF.CaseUtility > MinUtil)
                            {
                                NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].OverMinUtilCase++;
                            }

                            NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.push_back(UIF);
                            UIF.CaseUtility = 0;
                        }
                        else // 已存在任何L2
                        {
                            L1_UtilityInfo UIF;
                            // Index
                            UIF.VecIndex = NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex;
                            UIF.VecIndex.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[ItemIdx - 1]);
                            // Tid
                            UIF.VecTid = NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecTid;
                            UIF.VecTid.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].TidArray[ItemIdx - 1]);

                            // Iu、Ru、Utility(讀VecDB更新資訊)
                            for (int k = 0; k < NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex.size(); k++)
                            {
                                UIF.VecIu.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].IuArray[NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex[k] - 1]);
                                UIF.VecUtility.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].UtilityArray[NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex[k] - 1]);
                            }
                            UIF.VecIu.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].IuArray[ItemIdx - 1]);
                            UIF.VecUtility.push_back(VecDataBase[NodeUC.L2_SeqInfo[i].sid].UtilityArray[ItemIdx - 1]);

                            // CaseUtility(拿更新後的VecUtility來加總)
                            for (int k = 0; k < UIF.VecUtility.size(); k++)
                            {
                                UIF.CaseUtility += UIF.VecUtility[k];
                            }

                            // 投影點的Index不同
                            if (NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo[NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.size() - 1].VecIndex.back() != VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[ItemIdx - 1])
                            {
                                // acu與SeqUt比較出SeqUt
                                if (UIF.CaseUtility > NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt)
                                {
                                    NIF.SumUt += (UIF.CaseUtility - NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt);
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt = UIF.CaseUtility;
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUtCase = UIF.VecIndex.back();
                                }

                                // ru!=0 且 acu+ru比現存的SeqPEU大，比較出SeqPEU
                                if (VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1] != 0 &&
                                    UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1] > NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU)
                                {
                                    NIF.SumPEU += ((UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1]) - NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU);
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU = UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1];
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEUCase = VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[ItemIdx - 1];
                                }

                                // Record Cases of Over MinUtil.
                                if (UIF.CaseUtility > MinUtil)
                                {
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].OverMinUtilCase++;
                                }

                                NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.push_back(UIF);
                            }
                            else // 投影點Index相同
                            {
                                if (UIF.CaseUtility > NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo[NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.size() - 1].CaseUtility)
                                {
                                    // Record Cases of Over MinUtil.
                                    if (UIF.CaseUtility > MinUtil && NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo[NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.size() - 1].CaseUtility <= MinUtil)
                                    {
                                        NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].OverMinUtilCase++;
                                    }
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.pop_back();
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.push_back(UIF);
                                }
                                // SeqUt、SumUt
                                if (UIF.CaseUtility > NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt)
                                {
                                    NIF.SumUt += (UIF.CaseUtility - NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt);
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt = UIF.CaseUtility;
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUtCase = UIF.VecIndex.back(); // Record Cases of Over MinUtil.
                                }
                                // SeqPEU、SumPEU
                                if (VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1] != 0 &&
                                    UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1] > NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU)
                                {
                                    NIF.SumPEU += ((UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1]) - NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU);
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU = UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1];
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEUCase = VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[ItemIdx - 1];
                                }
                            }
                            UIF.CaseUtility = 0;
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }

        // IdxPos
        for (int a = 0; a < NIF.L2_SeqInfo.size(); a++)
        {
            for (int b = 0; b < NIF.L2_SeqInfo[a].L1_UtInfo.size(); b++)
            {
                for (int c = 0; c < NIF.L2_SeqInfo[a].L1_UtInfo[b].VecIndex.size(); c++)
                {
                    NIF.L2_SeqInfo[a].IdxPos[NIF.L2_SeqInfo[a].L1_UtInfo[b].VecIndex[c]].push_back({b, c});
                }
            }
        }

        if (NIF.SumUt >= MinUtil)
        {
            S_ExtensionCounter++;
            // cout << NIF.pattern << " | " << NIF.SumUt << endl;
            HUSPattern.push_back(NIF);
            //Cout_HUSPL3(NIF);
        }
        if (NIF.SumUt > maxut)
        {             
            maxut = NIF.SumUt;
        }
        HUSP(NIF);
    }
}

void WriteOutHUSPattern(vector<L3_NodeInfo> HUSPattern, double time)
{
    string Filename = "Remining_" + FileName_DB + "_Minutil_" + to_string((int)MinUtil) + ".txt";
    ofstream writeFile(Filename);

    writeFile << "Database:" << FileName_DB << "\n";
    writeFile << "=====================" << "\n";
    writeFile << "SWU : " << SumSWU << endl;
    writeFile << "*** MinUtil = " << MinUtil << " ***" << "\n";
    writeFile << "Single item Counter : " << Single_ItemCounter << "\n";
    writeFile << "I_Extension Counter : " << I_ExtensionCounter << "\n";
    writeFile << "S_Extension Counter : " << S_ExtensionCounter << "\n";
    writeFile << "Total Non hiding pattern number : " << Single_ItemCounter + I_ExtensionCounter + S_ExtensionCounter << "\n";
    writeFile << "Memory Usage (Working Set Size): " << memoryMB << " MB\n";
    writeFile << "Time consumption:" << time << "s" << "\n";
    writeFile << "=====================" << "\n";
    writeFile << "\n";

    if (HUSPattern.empty())
    {
        writeFile << "No HUSP." << "\n";
    }
    else
    {
        for (int i = 0; i < HUSPattern.size(); i++)
        {
            writeFile << "*** Pattern:" << HUSPattern[i].pattern << " ***\n";
            writeFile << "* SumUt:" << HUSPattern[i].SumUt << " *\n";
            for (int j = 0; j < HUSPattern[i].L2_SeqInfo.size(); j++)
            {
                writeFile << "Seq: " << HUSPattern[i].L2_SeqInfo[j].sid << "\n";
                writeFile << "SeqUt:" << HUSPattern[i].L2_SeqInfo[j].SeqUt << "\n";
                for (int k = 0; k < HUSPattern[i].L2_SeqInfo[j].L1_UtInfo.size(); k++)
                {
                    writeFile << "Index:\n";
                    writeFile << "[";
                    for (int h = 0; h < HUSPattern[i].L2_SeqInfo[j].L1_UtInfo[k].VecIndex.size(); h++)
                    {
                        writeFile << HUSPattern[i].L2_SeqInfo[j].L1_UtInfo[k].VecIndex[h] << " ";
                    }
                    writeFile << "]\n";
                    writeFile << "Iu:\n";
                    writeFile << "[";
                    for (int h = 0; h < HUSPattern[i].L2_SeqInfo[j].L1_UtInfo[k].VecIndex.size(); h++)
                    {
                        writeFile << HUSPattern[i].L2_SeqInfo[j].L1_UtInfo[k].VecIu[h] << " ";
                    }
                    writeFile << "]\n";
                    writeFile << "Utility:\n";
                    writeFile << "[";
                    for (int h = 0; h < HUSPattern[i].L2_SeqInfo[j].L1_UtInfo[k].VecIndex.size(); h++)
                    {
                        writeFile << HUSPattern[i].L2_SeqInfo[j].L1_UtInfo[k].VecUtility[h] << " ";
                    }
                    writeFile << "]";
                    writeFile << HUSPattern[i].L2_SeqInfo[j].L1_UtInfo[k].CaseUtility << "\n";
                }
            }
            writeFile << "\n";
        }
    }
}

void WriteOutOnlyHUSPattern(vector<L3_NodeInfo> HUSPattern)
{
    string Filename = "Remining_" + FileName_DB + "_Minutil_" + to_string((int)MinUtil) + ".txt";
    ofstream writeFile(Filename);

    if (HUSPattern.empty())
    {
        writeFile << "No HUSP." << "\n";
    }
    else
    {
        for (int i = 0; i < HUSPattern.size(); i++)
        {
            writeFile << HUSPattern[i].pattern << "\n";
        }
    }
}

int main()
{
    SeqData SD;
    SD.sid = 0;
    VecDataBase.push_back(SD); // VecDataBase[0]留空，item從1開始記錄

    L3_NodeInfo EmptyL3;
    Node_SingleItem.push_back(EmptyL3);

    vector<int> v(1, 0);
    ExternalUt.insert(make_pair(0, 0));

    //FileName_Eu = "jzwpaper_utb.txt";
    //FileName_DB = "jzwpaper_db.txt";
    //FileName_DB = "Adjust_jzwpaper_db.txt_Minutil_50.txt";

    FileName_Eu = "simple_utb.txt";
    //FileName_DB = "simple_db.txt";
    FileName_DB = "Adjust_simple_db.txt_Minutil_584.txt";

    //FileName_Eu = "simple2_utb.txt";
    //FileName_DB = "simple2_db.txt";
    //FileName_DB = "Adjust_simple2_db.txt_Minutil_154.txt";
    
    //FileName_Eu = "4_sign_ExternalUtility.txt";
    //FileName_DB = "Adjust_4_sign.txt_Minutil_7200.txt";
    
    //FileName_Eu = "01.bible_ExternalUtility.txt";
    //FileName_DB = "Adjust_01.bible.txt_Minutil_90000.txt";
    
    //FileName_Eu = "03.bms2_ExternalUtility.txt";
    //FileName_DB = "Adjust_03.bms2.txt_Minutil_26000.txt";

    //FileName_Eu = "05.foodmart_ExternalUtility.txt";
    //FileName_DB = "Adjust_05.foodmart.txt_Minutil_1000.txt";
    //FileName_DB = "05.foodmart.txt";
    
    //FileName_Eu = "7.RSCS_ExternalUtility.txt";
    //FileName_DB = "Adjust_7.RSCS.txt_Minutil_17500.txt";
    //FileName_DB = "7.RSCS.txt";
    
    MinUtil = 584;
    cout << endl;
    cout << "--- (Remining)Min utility = " << MinUtil << " ---" << endl;
    cout << "--- (Remining)Database : " << FileName_DB << " ---" << endl;
    cout << "--- (Remining)Eu : " << FileName_Eu << " ---" << endl;
    // cout << "========= Remining HUSP process ========" << endl;
    // cout << "--- MinUtil:" << MinUtil << " ---" << endl;
    // cout << "Start Reading External Utility..." << endl;
    Read_ExternalUt(FileName_Eu);
    // cout << "Start Reading Database..." << endl;
    Read_Database(FileName_DB);

    // cout << "Start Building Single Items..." << endl;
    BulidSingleItems(VecDataBase);

    //Cout_VecDB(VecDataBase);
    //Cout_ExternalUt(ExternalUt);

    double Threshold = 1;
    //  double MinUtil = SumSWU * Threshold;
    //  double MinUtil = 0;
    //  cout << MinUtil << "=" << SumSWU << "*" << Threshold << endl;

    double time = 0;
    clock_t start, end;
    start = clock();

    for (int i = 1; i < Node_SingleItem.size(); i++)
    {
        if (Node_SingleItem[i].SumUt >= MinUtil)
        {
            Single_ItemCounter++;
            HUSPattern.push_back(Node_SingleItem[i]);
            //cout << "Pattern: " << Node_SingleItem[i].pattern << "-->SumUt = " << Node_SingleItem[i].SumUt << endl;
        }
        if (Node_SingleItem[i].SumUt > maxut)
        {
            maxut = Node_SingleItem[i].SumUt;
        }
        // cout << "Start HUSP of " << Node_SingleItem[i].pattern << " ..." << endl;
        HUSP(Node_SingleItem[i]);
    }
    end = clock();
    cout << endl;
    cout << "===== Remining Result =====" << endl;
    cout << "SWU : " << SumSWU << endl;
    // cout << "Threshold : " << Threshold << endl;
    cout << "--- MinUtil = " << MinUtil << " ***" << endl;
    cout << "Single item Counter : " << Single_ItemCounter << endl;
    cout << "I_Extension Counter : " << I_ExtensionCounter << endl;
    cout << "S_Extension Counter : " << S_ExtensionCounter << endl;
    cout << "Total Non hiding pattern Counter : " << Single_ItemCounter + I_ExtensionCounter + S_ExtensionCounter << endl;
    cout << endl;
    cout << "MaxUt = " << maxut << endl;

    // Memory Usage
    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
    {
        memoryMB = pmc.WorkingSetSize / (1024.0 * 1024.0); // 轉成 MB
        cout << "Memory Usage (Working Set Size): " << memoryMB << " MB" << endl;
    }
    /////////////////////
    time = (double)(end - start) / CLOCKS_PER_SEC;
    cout << "Time consumption:" << time << "s" << endl;
    cout << endl;

    cout << "Start Write out HUS pattern ..." << endl;
    //WriteOutHUSPattern(HUSPattern, time);
    WriteOutOnlyHUSPattern(HUSPattern);

    cout << "========= Remining HUSP process END ========" << endl;
    cout << endl;

    return 0;
}
