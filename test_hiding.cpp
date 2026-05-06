#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <map>
#include <time.h>
#include <cmath>
#include <cstdlib>
#include <functional> // vector<reference_wrapper>
#include <queue>
#include <list>
#include <windows.h>
#include <psapi.h>
#include <memory>
#include <algorithm>

using namespace std;
/*
- 投影點夠扣、不夠扣設為1、Node和DB資料一致完成
- 當一次做a、b時，因為a的Node可能被b修改，因此Node a的資訊未必是最新的(但能保證不超過門檻)
    : 如1,2,4,4在Idx 1、4、6、9，如果只做a這顆子樹，則Idx 9的資訊與最後DB是一致的；
    但做a、b時，1,2,4,4的資訊會因為2,4的修改，而與DB不一致。
    但不影響輸出結果，頂多因為被其他Node修改(如b)，使Node a紀錄的utility比實際上多但不超過門檻
    (*尚未往上層找投影點扣，因此此版Node a的樹有超過門檻的情況發生)
*/

/* 
(2025.11.14目前要解)，2025.12.20已解
若Pattern有2個Seq:
一個Seq無法再被減(所有Iu皆為1)，仍有應減去而未減去的Utility
一個雖仍有Iu可減，但已減掉該Seq應減去的Utility
此時會導致Hiding failure
(如Pattern 1-2,4,4)

解決:若仍有Utility應減去而未減，把未減的值傳遞到下一筆序列，若整輪跑完還有未減完的，最多跑到第二輪。
*/
/* 
(2026.01.13) MDU
為了減少序列隱藏失敗的次數，新增一個最大可扣(MDU):Pattern的所有iu扣到剩1後，還有多少Utility可以扣。
利用這個最大可扣去算rut，可以使有Utility，可扣值已經極低甚至是0的序列降低rut的值。
也可以讓Utility更大的去承擔更多的rut去扣。
 */
class SeqData
{
public:
    int sid;
    vector<int> IndexArray;
    vector<int> TidArray;
    vector<int> ItemArray;
    vector<int> IuArray;
    vector<double> UtilityArray;
    vector<double> RuArray;

    map<int, vector<int>> ItemIdxTable; //<Item,[Idx set]>
    vector<int> IdxTid;                 // size n代表Tid n、<Tid起始Idx>
};
vector<SeqData> VecDataBase;

class L1_UtilityInfo
{
public:
    int VecIndex;     
    int VecIu;         
    double VecUtility; 
    double CaseUtility = 0;
    int IdxOfLastLevelIns = 0;
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

    vector<pair<int, pair<int, int>>> IdxPos; // map<index of item,position of index in L2>
};

class L3_NodeInfo
{
public:
    string pattern;
    vector<int> VecPattern;
    double SumPEU = 0;
    double SumUt = 0;
    vector<L2_SequenceInfo> L2_SeqInfo;
    int ExtensionType; // 0:I，1:S
};
vector<L3_NodeInfo> Node_SingleItem;
map<int, int> IdxNodeSingleItem; // SingleItem, Node_SingItem<>的位置

vector<int> ItemQueue; // 紀錄Externalfile中讀取Item的先後順序，來放入VecDataBase中
map<int, double> ExternalUt;
string str_EuFile = "";
string str_DBFile = "";

vector<reference_wrapper<L3_NodeInfo>> PatternPath;

int Single_ItemCounter = 0;
int I_ExtensionCounter = 0;
int S_ExtensionCounter = 0;
int TotalPatternCounter = 0;
double SumSWU = 0;
double MinUtil;
double memoryMB = 0;
double SumDiff = 0;

inline double cleanUtil(double val) {
    if (abs(val) < 1e-6) return 0.0;
    return std::round(val * 10000.0) / 10000.0;
}

void Cout_ExternalUt(map<int, double> ExternalUt)
{
    for (auto a : ExternalUt)
    {
        cout << a.first << "," << a.second << endl;
    }
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
            SD.IdxTid.push_back(0);
            SD.IdxTid.push_back(1);

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
                    if (IndexOfSpace == 0)
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
                            if (SU - (a.second * stod(strline.substr(IndexOfLeftBrackets + 1, IndexOfRightBrackets - IndexOfLeftBrackets - 1))) < 001)
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

void BulidSingleItems(vector<SeqData> &VecDataBase)
{
    int n = 1;
    for (auto Item : ExternalUt)
    {
        if (Item.first == 0)
            continue;
        L3_NodeInfo L3_node;
        L3_node.pattern = to_string(Item.first);
        L3_node.VecPattern.push_back(Item.first);
        L3_node.ExtensionType = 2;
        Node_SingleItem.push_back(L3_node);

        IdxNodeSingleItem.insert(make_pair(Item.first, n));
        n++;
    }

    for (int i = 1; i < VecDataBase.size(); i++)
    {
        for (int j = 0; j < VecDataBase[i].IndexArray.size(); j++)
        {
            if (!Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.size() ||
                Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.back().sid != VecDataBase[i].sid)
            {
                L2_SequenceInfo L2_seq;
                L2_seq.sid = VecDataBase[i].sid;
                L2_seq.SeqUtCase = VecDataBase[i].IndexArray[j];
                Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]].L2_SeqInfo.push_back(L2_seq);

                L1_UtilityInfo L1_ut;
                L1_ut.VecIndex = VecDataBase[i].IndexArray[j];
                L1_ut.VecIu = VecDataBase[i].IuArray[j];
                L1_ut.VecUtility = VecDataBase[i].UtilityArray[j];
                L1_ut.CaseUtility = VecDataBase[i].UtilityArray[j];

                auto &targetNode = Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]];
                targetNode.L2_SeqInfo.back().L1_UtInfo.push_back(L1_ut);
                targetNode.L2_SeqInfo.back().SeqUt = VecDataBase[i].UtilityArray[j];
                targetNode.SumUt += VecDataBase[i].UtilityArray[j];

                if (VecDataBase[i].RuArray[j] > 0)
                {
                    targetNode.L2_SeqInfo.back().SeqPEU = VecDataBase[i].UtilityArray[j] + VecDataBase[i].RuArray[j];
                    targetNode.SumPEU += VecDataBase[i].UtilityArray[j] + VecDataBase[i].RuArray[j];
                }

                // 建立 IdxPos
                targetNode.L2_SeqInfo.back().IdxPos.push_back({VecDataBase[i].IndexArray[j],
                                                               {targetNode.L2_SeqInfo.back().L1_UtInfo.size() - 1, 0}});
            }
            else
            {
                L1_UtilityInfo L1_ut;
                L1_ut.VecIndex = VecDataBase[i].IndexArray[j];
                L1_ut.VecIu = VecDataBase[i].IuArray[j];
                L1_ut.VecUtility = VecDataBase[i].UtilityArray[j];
                L1_ut.CaseUtility = VecDataBase[i].UtilityArray[j];

                auto &targetNode = Node_SingleItem[IdxNodeSingleItem[VecDataBase[i].ItemArray[j]]];
                targetNode.L2_SeqInfo.back().L1_UtInfo.push_back(L1_ut);

                if (targetNode.L2_SeqInfo.back().SeqPEU < (VecDataBase[i].UtilityArray[j] + VecDataBase[i].RuArray[j]) && VecDataBase[i].RuArray[j] > 0)
                {
                    targetNode.SumPEU -= targetNode.L2_SeqInfo.back().SeqPEU;
                    targetNode.SumPEU += VecDataBase[i].UtilityArray[j] + VecDataBase[i].RuArray[j];
                    targetNode.L2_SeqInfo.back().SeqPEU = VecDataBase[i].UtilityArray[j] + VecDataBase[i].RuArray[j];
                }

                if (targetNode.L2_SeqInfo.back().SeqUt < VecDataBase[i].UtilityArray[j])
                {
                    targetNode.SumUt -= targetNode.L2_SeqInfo.back().SeqUt;
                    targetNode.SumUt += VecDataBase[i].UtilityArray[j];
                    targetNode.L2_SeqInfo.back().SeqUt = VecDataBase[i].UtilityArray[j];
                    targetNode.L2_SeqInfo.back().SeqUtCase = VecDataBase[i].IndexArray[j];
                }

                // 建立 IdxPos
                targetNode.L2_SeqInfo.back().IdxPos.push_back({VecDataBase[i].IndexArray[j],
                                                               {targetNode.L2_SeqInfo.back().L1_UtInfo.size() - 1, 0}});
            }
        }
    }
}

void applyDeltaOnDB(int sid, int idx1b, int diffIu)
{
    if (diffIu <= 0)
        return;

    int idx0 = idx1b - 1;
    auto &db = VecDataBase[sid];

    int item = db.ItemArray[idx0];
    auto it_eu = ExternalUt.find(item);
    if (it_eu == ExternalUt.end())
        return;

    double eu = it_eu->second;
    double deltaU = diffIu * eu;

    db.UtilityArray[idx0] -= deltaU;
    db.UtilityArray[idx0] = cleanUtil(db.UtilityArray[idx0]);

    for (int k = 0; k < idx0; ++k)
    {
        db.RuArray[k] -= deltaU;
        db.RuArray[k] = cleanUtil(db.RuArray[k]);
    }
}

void applyDeltaAlongPath(
    vector<reference_wrapper<L3_NodeInfo>> &PatternPath,
    int topSeqIdx,
    int idx1b,
    int diffIu,
    double usedDelta,
    int sid)
{
    if (diffIu <= 0 || usedDelta <= 0)
        return;

    int curSeqIdx = topSeqIdx;

    for (int pi = (int)PatternPath.size() - 1; pi >= 0; --pi)
    {
        auto &node = PatternPath[pi].get();
        auto &seq = node.L2_SeqInfo[curSeqIdx];

        for (auto &pos_record : seq.IdxPos)
        {
            if (pos_record.first == idx1b)
            {
                int p = pos_record.second.first;
                auto &inst = seq.L1_UtInfo[p];

                inst.CaseUtility -= usedDelta;
                inst.CaseUtility = cleanUtil(inst.CaseUtility);

                if (inst.VecIndex == idx1b)
                {
                    inst.VecUtility -= usedDelta;
                    inst.VecUtility = cleanUtil(inst.VecUtility);
                    inst.VecIu -= diffIu;
                    if (inst.VecIu < 1)
                        inst.VecIu = 1;
                }
            }
        }

        node.SumUt -= seq.SeqUt;
        node.SumPEU -= seq.SeqPEU;

        seq.SeqUt = 0;
        seq.SeqPEU = 0;
        seq.SeqUtCase = -1;

        for (int x = 0; x < (int)seq.L1_UtInfo.size(); ++x)
        {
            auto &inst2 = seq.L1_UtInfo[x];

            if (inst2.CaseUtility > seq.SeqUt)
            {
                node.SumUt -= seq.SeqUt;
                node.SumUt += inst2.CaseUtility;
                seq.SeqUt = inst2.CaseUtility;
                seq.SeqUtCase = inst2.VecIndex;
            }

            int idx1b_x = inst2.VecIndex;
            int idx0_x = idx1b_x - 1;
            double ru_x = VecDataBase[seq.sid].RuArray[idx0_x];

            if (ru_x > 0)
            {
                double candPEU = inst2.CaseUtility + ru_x;
                if (candPEU > seq.SeqPEU)
                {
                    node.SumPEU -= seq.SeqPEU;
                    node.SumPEU += candPEU;
                    seq.SeqPEU = candPEU;
                }
            }
        }
        curSeqIdx = seq.IdxOfLastLevelSeq;
        if (curSeqIdx < 0)
            break;
    }
}

void UpdateSingleItem(vector<L3_NodeInfo> &Node_SingleItem, int Item)
{
    int IdxItem = IdxNodeSingleItem[Item];
    Node_SingleItem[IdxItem].SumUt = 0;
    Node_SingleItem[IdxItem].SumPEU = 0;

    for (int i = 0; i < Node_SingleItem[IdxItem].L2_SeqInfo.size(); i++)
    {
        Node_SingleItem[IdxItem].L2_SeqInfo[i].SeqUt = 0;
        Node_SingleItem[IdxItem].L2_SeqInfo[i].SeqPEU = 0;

        for (int j = 0; j < Node_SingleItem[IdxItem].L2_SeqInfo[i].L1_UtInfo.size(); j++)
        {
            auto &inst = Node_SingleItem[IdxItem].L2_SeqInfo[i].L1_UtInfo[j];
            int sid = Node_SingleItem[IdxItem].L2_SeqInfo[i].sid;
            int idx0 = inst.VecIndex - 1;

            inst.VecIu = VecDataBase[sid].IuArray[idx0];
            inst.VecUtility = VecDataBase[sid].UtilityArray[idx0];
            inst.CaseUtility = VecDataBase[sid].UtilityArray[idx0];

            if (inst.CaseUtility > Node_SingleItem[IdxItem].L2_SeqInfo[i].SeqUt)
            {
                Node_SingleItem[IdxItem].SumUt -= Node_SingleItem[IdxItem].L2_SeqInfo[i].SeqUt;
                Node_SingleItem[IdxItem].L2_SeqInfo[i].SeqUt = inst.CaseUtility;
                Node_SingleItem[IdxItem].SumUt += Node_SingleItem[IdxItem].L2_SeqInfo[i].SeqUt;
                Node_SingleItem[IdxItem].L2_SeqInfo[i].SeqUtCase = VecDataBase[sid].IndexArray[idx0];
            }

            if (VecDataBase[sid].RuArray[idx0] != 0)
            {
                if (inst.CaseUtility + VecDataBase[sid].RuArray[idx0] > Node_SingleItem[IdxItem].L2_SeqInfo[i].SeqPEU)
                {
                    Node_SingleItem[IdxItem].SumPEU -= Node_SingleItem[IdxItem].L2_SeqInfo[i].SeqPEU;
                    Node_SingleItem[IdxItem].L2_SeqInfo[i].SeqPEU = inst.CaseUtility + VecDataBase[sid].RuArray[idx0];
                    Node_SingleItem[IdxItem].SumPEU += Node_SingleItem[IdxItem].L2_SeqInfo[i].SeqPEU;
                }
            }
            else
            {
                Node_SingleItem[IdxItem].L2_SeqInfo[i].SeqPEU = 0;
            }
        }
    }
}

void Cout_HUSPL3(L3_NodeInfo L3_node)
{
    cout << "---------- Pattern:" << L3_node.pattern << " ----------" << endl;

    cout << "=== MinUtil:" << MinUtil << " ===" << endl;
    cout << "** SumUt:" << L3_node.SumUt << " **" << endl;
    cout << "** SumPEU:" << L3_node.SumPEU << " **" << endl;
    cout << "Extension Type:" << L3_node.ExtensionType << endl;

    for (int j = 0; j < L3_node.L2_SeqInfo.size(); j++)
    {
        cout << "Seq " << L3_node.L2_SeqInfo[j].sid << endl;
        cout << "SeqUt:" << L3_node.L2_SeqInfo[j].SeqUt << endl;
        cout << "SeqPEU:" << L3_node.L2_SeqInfo[j].SeqPEU << endl;
        cout << "IdxOfLastLevelSeq:" << L3_node.L2_SeqInfo[j].IdxOfLastLevelSeq << endl;
        cout << "|Idx Iu Ut|InstanceUt" << endl;
        
        for (int k = 0; k < L3_node.L2_SeqInfo[j].L1_UtInfo.size(); k++)
        {
            auto &inst = L3_node.L2_SeqInfo[j].L1_UtInfo[k];
            
            vector<int> path_idx;
            vector<int> path_iu;
            vector<double> path_ut;

            int tempIns = k;
            int tempSeq = j;

            // 從當前這層，沿著 PatternPath 一路往上爬到 Root
            for (int level = (int)PatternPath.size() - 1; level >= 0; --level) 
            {
                auto &trace_node = PatternPath[level].get();
                auto &trace_inst = trace_node.L2_SeqInfo[tempSeq].L1_UtInfo[tempIns];

                // 收集資料
                path_idx.push_back(trace_inst.VecIndex);
                path_iu.push_back(trace_inst.VecIu);
                path_ut.push_back(trace_inst.VecUtility);

                // 往上跳一層
                tempIns = trace_inst.IdxOfLastLevelIns;
                tempSeq = trace_node.L2_SeqInfo[tempSeq].IdxOfLastLevelSeq;
            }

            // 因為是從葉子往上爬，順序是反的，所以我們把它反轉回來 (變成從 Root 到 Leaf)
            reverse(path_idx.begin(), path_idx.end());
            reverse(path_iu.begin(), path_iu.end());
            reverse(path_ut.begin(), path_ut.end());

            cout << inst.IdxOfLastLevelIns << " ";
            
            cout << "[";
            for (int val : path_idx) cout << val << " ";
            cout << "] [";
            
            for (int val : path_iu) cout << val << " ";
            cout << "] [";
            
            for (double val : path_ut) cout << val << " ";
            cout << "] ";

            cout << inst.CaseUtility << endl;
        }

        // 印出群組化後的定位表 (IdxPos)
        map<int, vector<pair<int, int>>> display_map;
        for (auto &record : L3_node.L2_SeqInfo[j].IdxPos) {
            display_map[record.first].push_back(record.second);
        }

        for (auto a : display_map)
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

void DeleteLowRSUFromlist(L3_NodeInfo &NodeUC, set<int> &ilist, set<int> &slist)
{
    map<int, double> S_RSUmap;
    map<int, double> I_RSUmap;
    for (int i = 0; i < NodeUC.L2_SeqInfo.size(); i++)
    {
        if (NodeUC.L2_SeqInfo[i].L1_UtInfo.empty())
            continue;
        int sid = NodeUC.L2_SeqInfo[i].sid;

        // ---- S-Extension ----
        int first_inst_idx = NodeUC.L2_SeqInfo[i].L1_UtInfo[0].VecIndex;
        if (first_inst_idx < VecDataBase[sid].IdxTid.back())
        {
            int current_tid = VecDataBase[sid].TidArray[first_inst_idx - 1];
            int next_tid_start_idx = VecDataBase[sid].IdxTid[current_tid + 1];

            for (int j = next_tid_start_idx - 1; j < VecDataBase[sid].IndexArray.size(); j++)
            {
                S_RSUmap[VecDataBase[sid].ItemArray[j]] += NodeUC.L2_SeqInfo[i].SeqPEU;
                if (S_RSUmap[VecDataBase[sid].ItemArray[j]] >= MinUtil) 
                {
                    slist.insert(VecDataBase[sid].ItemArray[j]);
                }
            }
        }

        // ---- I-Extension ----
        for (int j = 0; j < NodeUC.L2_SeqInfo[i].L1_UtInfo.size(); j++)
        {
            int cur_idx = NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex;

            for (int k = cur_idx; k < VecDataBase[sid].IndexArray.size(); k++)
            {
                if (VecDataBase[sid].TidArray[cur_idx - 1] != VecDataBase[sid].TidArray[k])
                    break;

                I_RSUmap[VecDataBase[sid].ItemArray[k]] += NodeUC.L2_SeqInfo[i].SeqPEU;
                if (I_RSUmap[VecDataBase[sid].ItemArray[k]] >= MinUtil) 
                {
                    ilist.insert(VecDataBase[sid].ItemArray[k]);
                }
            }
        }
    }
}

double TraceBack(vector<reference_wrapper<L3_NodeInfo>> &PatternPath,
                 int TopNodeSeqIdx,
                 int LastSeqIdx,
                 int LastInsIdx,
                 double KeepReduceUt,
                 int PathIdx,
                 int ExtendItemIdx)
{
    if (KeepReduceUt <= 0)
        return 0;
    if (PathIdx > (int)PatternPath.size())
        return KeepReduceUt;

    int curIdx = (int)PatternPath.size() - PathIdx;
    if (curIdx < 0 || curIdx >= (int)PatternPath.size())
        return KeepReduceUt;

    bool fromSExt = false;
    if (curIdx + 1 < (int)PatternPath.size())
    {
        if (PatternPath[curIdx + 1].get().ExtensionType == 1)
            fromSExt = true;
    }

    auto &node = PatternPath[curIdx].get();
    auto &seq = node.L2_SeqInfo[LastSeqIdx];
    auto &inst = seq.L1_UtInfo[LastInsIdx];

    double KeepSeqUt = seq.SeqUt;
    double KeepInsUt = inst.CaseUtility;

    int idx1b = inst.VecIndex;
    int idx0 = idx1b - 1;
    int sid = seq.sid;

    if (fromSExt && idx1b > ExtendItemIdx)
    {
        if (PathIdx + 1 <= (int)PatternPath.size())
        {
            int nextSeqIdx = seq.IdxOfLastLevelSeq;
            int nextInsIdx = inst.IdxOfLastLevelIns;
            return TraceBack(PatternPath, TopNodeSeqIdx,
                             nextSeqIdx, nextInsIdx,
                             KeepReduceUt, PathIdx + 1, ExtendItemIdx);
        }
        return KeepReduceUt;
    }
    //cout << seq.sid << ": " << inst.VecIndex <<"InsUt: " << KeepInsUt << "-->SeqUt: " << KeepSeqUt << "-->CurRut: " << KeepReduceUt << endl;
    if (KeepInsUt <= (KeepSeqUt - KeepReduceUt))
    {
        //cout << seq.sid << ": " << inst.VecIndex << "--> Ut <= SeqUt - rut" << endl;
        if (PathIdx + 1 <= (int)PatternPath.size())
        {
            int nextSeqIdx = seq.IdxOfLastLevelSeq;
            int nextInsIdx = inst.IdxOfLastLevelIns;
            return TraceBack(PatternPath, TopNodeSeqIdx,
                             nextSeqIdx, nextInsIdx,
                             KeepReduceUt, PathIdx + 1, ExtendItemIdx);
        }
        return KeepReduceUt;
    }

    int item = node.VecPattern.back();
    auto it_eu = ExternalUt.find(item);
    if (it_eu == ExternalUt.end())
    {
        if (PathIdx + 1 <= (int)PatternPath.size())
        {
            int nextSeqIdx = seq.IdxOfLastLevelSeq;
            int nextInsIdx = inst.IdxOfLastLevelIns;
            return TraceBack(PatternPath, TopNodeSeqIdx,
                             nextSeqIdx, nextInsIdx,
                             KeepReduceUt, PathIdx + 1, ExtendItemIdx);
        }
        return KeepReduceUt;
    }
    double eu = it_eu->second;

    int oldIu = VecDataBase[sid].IuArray[idx0];
    if (oldIu <= 1)
    {
        if (PathIdx + 1 <= (int)PatternPath.size())
        {
            int nextSeqIdx = seq.IdxOfLastLevelSeq;
            int nextInsIdx = inst.IdxOfLastLevelIns;
            return TraceBack(PatternPath, TopNodeSeqIdx,
                             nextSeqIdx, nextInsIdx,
                             KeepReduceUt, PathIdx + 1, ExtendItemIdx);
        }
        return KeepReduceUt;
    }
    

    int maxDiffIu = oldIu - 1;
    double maxDeltaU = maxDiffIu * eu;

    double needU = KeepReduceUt;
    double usedDelta = min(needU, maxDeltaU);

    int diffIu = (int)ceil(cleanUtil(usedDelta / eu));
    if (diffIu > maxDiffIu)
        diffIu = maxDiffIu;
    usedDelta = diffIu * eu;

    //cout << seq.sid << ": " << inst.VecIndex << "--> currut: "<< usedDelta <<endl;

    int newIu = oldIu - diffIu;
    VecDataBase[sid].IuArray[idx0] = newIu;

    applyDeltaOnDB(sid, idx1b, diffIu);
    applyDeltaAlongPath(PatternPath, TopNodeSeqIdx, idx1b, diffIu, usedDelta, sid);

    KeepReduceUt -= usedDelta;
    if (KeepReduceUt <= 0)
        return 0;
    if (PathIdx + 1 <= (int)PatternPath.size())
    {
        int nextSeqIdx = seq.IdxOfLastLevelSeq;
        int nextInsIdx = inst.IdxOfLastLevelIns;
        return TraceBack(PatternPath, TopNodeSeqIdx,
                         nextSeqIdx, nextInsIdx,
                         KeepReduceUt, PathIdx + 1, ExtendItemIdx);
    }
    return KeepReduceUt;
}

void REIHUSP_hiding(vector<reference_wrapper<L3_NodeInfo>> &PatternPath)
{
    if (PatternPath.empty())
        return;

    L3_NodeInfo &leafNode = PatternPath.back().get();

    double diff = leafNode.SumUt - MinUtil + 1;
    if (diff <= 0)
        return;

    int item = leafNode.VecPattern.back();
    double eu = 0;
    auto it_eu = ExternalUt.find(item);

    if (it_eu != ExternalUt.end())
    {
        eu = it_eu->second;
    }

    vector<double> VecSeqMDU(leafNode.L2_SeqInfo.size(), 0.0);
    double TotalMDU = 0;

    if (eu > 0)
    {
        for (int i = 0; i < (int)leafNode.L2_SeqInfo.size(); ++i)
        {
            auto &seq = leafNode.L2_SeqInfo[i];
            double seqMDU = 0;
            double MaxseqMDU = 0;
            for (int j = 0; j < (int)seq.L1_UtInfo.size(); ++j)
            {
                double seqMDU = 0;
                int curSeqIdx = i;
                int curInsIdx = j;

                for (int k = (int)PatternPath.size() - 1; k >= 0; --k)
                {
                    auto &trace_node = PatternPath[k].get();
                    int curIu = trace_node.L2_SeqInfo[curSeqIdx].L1_UtInfo[curInsIdx].VecIu;
                    int currentItem = trace_node.VecPattern.back();

                    if (curIu > 1)
                    {
                        seqMDU += (curIu - 1) * ExternalUt[currentItem];
                    }

                    int nextSeqIdx = trace_node.L2_SeqInfo[curSeqIdx].IdxOfLastLevelSeq;
                    int nextInsIdx = trace_node.L2_SeqInfo[curSeqIdx].L1_UtInfo[curInsIdx].IdxOfLastLevelIns;
                    curSeqIdx = nextSeqIdx;
                    curInsIdx = nextInsIdx;
                }
                MaxseqMDU = max(MaxseqMDU, seqMDU);
            }
            VecSeqMDU[i] = MaxseqMDU;
            TotalMDU += MaxseqMDU;
        }
    }
    /*    if (diff > TotalMDU) {
        cout << "[FindPattern] " << leafNode.pattern << " Out Of Math LIMIT" 
             << " (TotalDiff: " << diff << ", TotalMDU:  " << TotalMDU << ")" << endl;
    }*/

    double curdiff = diff;
    double Unpaidrut = 0;
    int RoundCounter = 0;
    int MaxRound = 2;
    while (TotalMDU > diff && curdiff > 0 && RoundCounter < MaxRound)
    {
        RoundCounter++;
        double RutInThisRound = 0;

        for (int i = 0; i < (int)leafNode.L2_SeqInfo.size(); ++i)
        {
            if (curdiff <= 0 && Unpaidrut <= 0)
                break;

            auto &seq = leafNode.L2_SeqInfo[i];
            double KeepSeqUt = seq.SeqUt;
            double KeepSumUt = leafNode.SumUt;

            if (KeepSeqUt <= 0)
                continue;

            double AllocUt = 0;

            if (TotalMDU > 0)
            {
                AllocUt = ceil(diff * (VecSeqMDU[i] / TotalMDU));
            }

            double TargetReduce = AllocUt + Unpaidrut;
            if (TargetReduce <= 0)
                continue;

            double ActualReduced = 0;
            double KeepReduceUt = TargetReduce;

            for (int j = 0; j < (int)seq.L1_UtInfo.size(); ++j)
            {
                if (KeepReduceUt <= 0)
                    break;

                auto &inst = seq.L1_UtInfo[j];
                double curReduceUt = KeepReduceUt;

                if (inst.CaseUtility <= (KeepSeqUt - KeepReduceUt))
                    continue;

                int sid = seq.sid;
                int idx1b = inst.VecIndex;
                int idx0 = idx1b - 1;

                double localDropped = 0; //本輪實際扣的Utility
                int curIu = inst.VecIu;

                if (curIu > 1 && eu > 0)
                {
                    int maxDiffIu = curIu - 1;
                    double maxDeltaU = maxDiffIu * eu;

                    double needU = KeepReduceUt;
                    double usedDelta = min(needU, maxDeltaU);

                    int diffIu = (int)ceil(cleanUtil(usedDelta / eu));
                    if (diffIu > maxDiffIu)
                        diffIu = maxDiffIu;
                    usedDelta = diffIu * eu;
                    int newIu = curIu - diffIu;

                    /*inst.VecIu = newIu;
                    inst.VecUtility -= usedDelta;
                    inst.VecUtility = cleanUtil(inst.VecUtility);
                    inst.CaseUtility -= usedDelta;
                    inst.CaseUtility = cleanUtil(inst.CaseUtility);

                    if (idx1b == seq.SeqUtCase) {
                        seq.SeqUt -= usedDelta;
                    }*/

                    VecDataBase[sid].IuArray[idx0] = newIu;
                    applyDeltaOnDB(sid, idx1b, diffIu);
                    applyDeltaAlongPath(PatternPath, i, idx1b, diffIu, usedDelta, sid);

                    localDropped += usedDelta;
                    curReduceUt -= usedDelta;
                }

                if (curReduceUt > 0)
                {
                    double remain = TraceBack(
                        PatternPath,
                        i,
                        seq.IdxOfLastLevelSeq,
                        inst.IdxOfLastLevelIns,
                        curReduceUt,
                        2,
                        idx1b);
                    double traceAmount = curReduceUt - remain;
                    localDropped += traceAmount;
                    curReduceUt = remain;
                }

                ActualReduced += localDropped;
            }

            RutInThisRound += ActualReduced;

            if (ActualReduced < TargetReduce)
            {
                Unpaidrut = TargetReduce - ActualReduced;
            }
            else
            {
                Unpaidrut = 0;
            }
        }

        curdiff -= RutInThisRound;
        if (curdiff < 0)
            curdiff = 0;
        if (RutInThisRound == 0 && curdiff > 0)
            break;
    }
}

void SingleItem_Hiding(vector<L3_NodeInfo> &Node_SingleItem, int Item)
{
    if (IdxNodeSingleItem.find(Item) == IdxNodeSingleItem.end())
        return;
    int IdxItem = IdxNodeSingleItem[Item];

    double diff = Node_SingleItem[IdxItem].SumUt - MinUtil + 1;
    if (diff <= 0)
        return;

    double eu = 0;
    auto it_eu = ExternalUt.find(Item);
    if (it_eu == ExternalUt.end())
        return;
    eu = it_eu->second;

    vector<double> VecSeqMDU(Node_SingleItem[IdxItem].L2_SeqInfo.size(), 0.0);
    double TotalMDU = 0;

    for (int i = 0; i < (int)Node_SingleItem[IdxItem].L2_SeqInfo.size(); i++)
    {
        double MaxseqMDU = 0;
        auto &seq = Node_SingleItem[IdxItem].L2_SeqInfo[i];
        
        for (auto &inst : seq.L1_UtInfo)
        {
            double seqMDU = 0; 
            if (inst.VecIu > 1)
            {
                seqMDU = (inst.VecIu - 1) * eu; 
            }
            MaxseqMDU = max(seqMDU, MaxseqMDU);
        }
        VecSeqMDU[i] = MaxseqMDU;
        TotalMDU += MaxseqMDU;
    }
    /*if (diff > TotalMDU) {
        cout << "[FindItem] " << Item << " Out Of Math LIMIT" 
             << " (TotalDiff: " << diff << ", TotalMDU:  " << TotalMDU << ")" << endl;
    }*/

    double curdiff = diff;
    double Unpaidrut = 0;
    int RoundCounter = 0;
    int MaxRound = 2;

    while (TotalMDU > diff && curdiff > 0 && RoundCounter < MaxRound) 
    {
        RoundCounter++;
        double RutInThisRound = 0;

        for (int i = 0; i < (int)Node_SingleItem[IdxItem].L2_SeqInfo.size(); i++)
        {
            if (curdiff <= 0 && Unpaidrut <= 0)
                break;

            auto &seq = Node_SingleItem[IdxItem].L2_SeqInfo[i];
            double KeepSeqUt = seq.SeqUt;

            if (KeepSeqUt <= 0)
                continue;

            double AllocUt = 0;
            if (TotalMDU > 0)
            {
                AllocUt = ceil(cleanUtil(diff * (VecSeqMDU[i] / TotalMDU)));
            }

            double TargetReduce = AllocUt + Unpaidrut;
            if (TargetReduce <= 0)
                continue;

            double KeepReduceUt = TargetReduce; 

            for (int j = 0; j < (int)seq.L1_UtInfo.size(); j++)
            {
                auto &inst = seq.L1_UtInfo[j];

                double curReduceUt = KeepReduceUt;

                if (inst.CaseUtility <= (KeepSeqUt - curReduceUt) + 1e-6)
                    continue;

                int sid = seq.sid;
                int idx1b = inst.VecIndex;
                int idx0 = idx1b - 1;
                int curIu = inst.VecIu;

                if (curIu > 1)
                {
                    int maxDiffIu = curIu - 1;
                    double maxDeltaU = maxDiffIu * eu;

                    double needU = curReduceUt;
                    double usedDelta = min(needU, maxDeltaU);

                    int diffIu = (int)ceil(cleanUtil(usedDelta / eu));
                    if (diffIu > maxDiffIu)
                        diffIu = maxDiffIu;
                    usedDelta = diffIu * eu;
                    int newIu = curIu - diffIu;

                    inst.VecUtility -= usedDelta;
                    inst.VecUtility = cleanUtil(inst.VecUtility);
                    inst.VecIu = newIu;
                    inst.CaseUtility -= usedDelta;
                    inst.CaseUtility = cleanUtil(inst.CaseUtility);

                    VecDataBase[sid].IuArray[idx0] = newIu;
                    applyDeltaOnDB(sid, idx1b, diffIu);
                    
                }
            }

            double NewSeqUt = 0;
            for (auto &inst : seq.L1_UtInfo) {
                if (inst.CaseUtility > NewSeqUt) {
                    NewSeqUt = inst.CaseUtility;
                }
            }

            double ActualReduced = KeepSeqUt - NewSeqUt;
            seq.SeqUt = NewSeqUt; 
            
            RutInThisRound += ActualReduced;

            if (ActualReduced < TargetReduce - 1e-6)
            {
                Unpaidrut = TargetReduce - ActualReduced;
            }
            else
            {
                Unpaidrut = 0;
            }
        }

        curdiff -= RutInThisRound;
        if (curdiff < 0)
            curdiff = 0;
        if (RutInThisRound == 0 && curdiff > 0)
            break;
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

    // I-Extenesion
    for (int Iitem : ilist)
    {
        L3_NodeInfo NIF;
        NIF.SumPEU = 0;
        NIF.SumUt = 0;
        NIF.pattern = NodeUC.pattern + "-" + to_string(Iitem);
        NIF.ExtensionType = 0;
        NIF.VecPattern = NodeUC.VecPattern;
        NIF.VecPattern.push_back(Iitem);

        for (int i = 0; i < NodeUC.L2_SeqInfo.size(); i++)
        {
            for (int ItemIdx : VecDataBase[NodeUC.L2_SeqInfo[i].sid].ItemIdxTable[Iitem])
            {
                for (int j = 0; j < NodeUC.L2_SeqInfo[i].L1_UtInfo.size(); j++)
                {
                    if (ItemIdx > NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex &&
                        VecDataBase[NodeUC.L2_SeqInfo[i].sid].TidArray[ItemIdx - 1] == VecDataBase[NodeUC.L2_SeqInfo[i].sid].TidArray[NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex - 1])
                    {
                        if (NIF.L2_SeqInfo.size() == 0 || NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].sid != VecDataBase[NodeUC.L2_SeqInfo[i].sid].sid)
                        {
                            L2_SequenceInfo SIF;
                            SIF.sid = VecDataBase[NodeUC.L2_SeqInfo[i].sid].sid;
                            SIF.IdxOfLastLevelSeq = i;
                            NIF.L2_SeqInfo.push_back(SIF);

                            L1_UtilityInfo UIF;
                            UIF.IdxOfLastLevelIns = j;
                            UIF.VecIndex = VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[ItemIdx - 1];
                            UIF.VecIu = VecDataBase[NodeUC.L2_SeqInfo[i].sid].IuArray[ItemIdx - 1];
                            UIF.VecUtility = VecDataBase[NodeUC.L2_SeqInfo[i].sid].UtilityArray[ItemIdx - 1];
                            UIF.CaseUtility = NodeUC.L2_SeqInfo[i].L1_UtInfo[j].CaseUtility + UIF.VecUtility;

                            NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt = UIF.CaseUtility;
                            NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUtCase = UIF.VecIndex;
                            NIF.SumUt += UIF.CaseUtility;

                            if (VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1] != 0)
                            {
                                NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU = UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1];
                                NIF.SumPEU += NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU;
                            }
                            else
                            {
                                NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU = 0;
                            }

                            if (UIF.CaseUtility > MinUtil)
                            {
                                NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].OverMinUtilCase++;
                            }

                            NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.push_back(UIF);
                            UIF.CaseUtility = 0;
                        }
                        else
                        {
                            L1_UtilityInfo UIF;
                            UIF.IdxOfLastLevelIns = j;
                            UIF.VecIndex = VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[ItemIdx - 1];
                            UIF.VecIu = VecDataBase[NodeUC.L2_SeqInfo[i].sid].IuArray[ItemIdx - 1];
                            UIF.VecUtility = VecDataBase[NodeUC.L2_SeqInfo[i].sid].UtilityArray[ItemIdx - 1];
                            UIF.CaseUtility = NodeUC.L2_SeqInfo[i].L1_UtInfo[j].CaseUtility + UIF.VecUtility;

                            if (NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo[NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.size() - 1].VecIndex != VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[ItemIdx - 1])
                            {
                                if (UIF.CaseUtility > NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt)
                                {
                                    NIF.SumUt += (UIF.CaseUtility - NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt);
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt = UIF.CaseUtility;
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUtCase = UIF.VecIndex;
                                }

                                if (VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1] != 0 &&
                                    UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1] > NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU)
                                {
                                    NIF.SumPEU += ((UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1]) - NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU);
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU = UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1];
                                }

                                if (UIF.CaseUtility > MinUtil)
                                {
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].OverMinUtilCase++;
                                }

                                NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.push_back(UIF);
                                UIF.CaseUtility = 0;
                            }
                            else
                            {
                                if (UIF.CaseUtility > NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo[NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.size() - 1].CaseUtility)
                                {
                                    if (UIF.CaseUtility > MinUtil && NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo[NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.size() - 1].CaseUtility <= MinUtil)
                                    {
                                        NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].OverMinUtilCase++;
                                    }
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.pop_back();
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.push_back(UIF);
                                }
                                if (UIF.CaseUtility > NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt)
                                {
                                    NIF.SumUt += (UIF.CaseUtility - NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt);
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt = UIF.CaseUtility;
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUtCase = UIF.VecIndex;
                                }
                                if (VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1] != 0 &&
                                    UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1] > NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU)
                                {
                                    NIF.SumPEU += ((UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1]) - NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU);
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU = UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1];
                                }
                                UIF.CaseUtility = 0;
                            }
                        }
                    }
                }
            }
        }

        for (int a = 0; a < NIF.L2_SeqInfo.size(); a++)
        {
            for (int b = 0; b < NIF.L2_SeqInfo[a].L1_UtInfo.size(); b++)
            {
                int leaf_idx = NIF.L2_SeqInfo[a].L1_UtInfo[b].VecIndex;
                int c_pos = NIF.VecPattern.size() - 1;
                NIF.L2_SeqInfo[a].IdxPos.push_back({leaf_idx, {b, c_pos}});

                int tempIns = NIF.L2_SeqInfo[a].L1_UtInfo[b].IdxOfLastLevelIns;
                int tempSeq = NIF.L2_SeqInfo[a].IdxOfLastLevelSeq;
                int level = PatternPath.size() - 1;

                while (level >= 0)
                {
                    auto &trace_node = PatternPath[level].get();
                    int idx = trace_node.L2_SeqInfo[tempSeq].L1_UtInfo[tempIns].VecIndex;
                    int ancestor_c_pos = trace_node.VecPattern.size() - 1;
                    NIF.L2_SeqInfo[a].IdxPos.push_back({idx, {b, ancestor_c_pos}});

                    tempIns = trace_node.L2_SeqInfo[tempSeq].L1_UtInfo[tempIns].IdxOfLastLevelIns;
                    tempSeq = trace_node.L2_SeqInfo[tempSeq].IdxOfLastLevelSeq;
                    level--;
                }
            }
        }

        PatternPath.push_back(ref(NIF));

        if (NIF.SumUt >= MinUtil)
        {
            Cout_HUSPL3(NIF);
            I_ExtensionCounter++;
            REIHUSP_hiding(PatternPath);
        }
        HUSP(NIF);

        PatternPath.pop_back();
    }

    // S-Extension
    for (int Sitem : slist)
    {
        L3_NodeInfo NIF;
        NIF.SumPEU = 0;
        NIF.SumUt = 0;
        NIF.pattern = NodeUC.pattern + "," + to_string(Sitem);
        NIF.ExtensionType = 1;

        NIF.VecPattern = NodeUC.VecPattern;
        NIF.VecPattern.push_back(Sitem);

        for (int i = 0; i < NodeUC.L2_SeqInfo.size(); i++)
        {
            for (int ItemIdx : VecDataBase[NodeUC.L2_SeqInfo[i].sid].ItemIdxTable[Sitem])
            {
                for (int j = 0; j < NodeUC.L2_SeqInfo[i].L1_UtInfo.size(); j++)
                {
                    if (ItemIdx > NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex &&
                        VecDataBase[NodeUC.L2_SeqInfo[i].sid].TidArray[ItemIdx - 1] > VecDataBase[NodeUC.L2_SeqInfo[i].sid].TidArray[NodeUC.L2_SeqInfo[i].L1_UtInfo[j].VecIndex - 1])
                    {
                        if (NIF.L2_SeqInfo.size() == 0 || NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].sid != VecDataBase[NodeUC.L2_SeqInfo[i].sid].sid)
                        {
                            L2_SequenceInfo SIF;
                            SIF.sid = VecDataBase[NodeUC.L2_SeqInfo[i].sid].sid;
                            SIF.IdxOfLastLevelSeq = i;
                            NIF.L2_SeqInfo.push_back(SIF);

                            L1_UtilityInfo UIF;
                            UIF.IdxOfLastLevelIns = j;
                            UIF.VecIndex = VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[ItemIdx - 1];
                            UIF.VecIu = VecDataBase[NodeUC.L2_SeqInfo[i].sid].IuArray[ItemIdx - 1];
                            UIF.VecUtility = VecDataBase[NodeUC.L2_SeqInfo[i].sid].UtilityArray[ItemIdx - 1];
                            UIF.CaseUtility = NodeUC.L2_SeqInfo[i].L1_UtInfo[j].CaseUtility + UIF.VecUtility;

                            NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt = UIF.CaseUtility;
                            NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUtCase = UIF.VecIndex;
                            NIF.SumUt += UIF.CaseUtility;

                            if (VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1] != 0)
                            {
                                NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU = UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1];
                                NIF.SumPEU += NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU;
                            }
                            else
                            {
                                NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU = 0;
                            }

                            if (UIF.CaseUtility > MinUtil)
                            {
                                NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].OverMinUtilCase++;
                            }

                            NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.push_back(UIF);
                            UIF.CaseUtility = 0;
                        }
                        else
                        {
                            L1_UtilityInfo UIF;
                            UIF.IdxOfLastLevelIns = j;
                            UIF.VecIndex = VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[ItemIdx - 1];
                            UIF.VecIu = VecDataBase[NodeUC.L2_SeqInfo[i].sid].IuArray[ItemIdx - 1];
                            UIF.VecUtility = VecDataBase[NodeUC.L2_SeqInfo[i].sid].UtilityArray[ItemIdx - 1];
                            UIF.CaseUtility = NodeUC.L2_SeqInfo[i].L1_UtInfo[j].CaseUtility + UIF.VecUtility;

                            if (NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo[NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.size() - 1].VecIndex != VecDataBase[NodeUC.L2_SeqInfo[i].sid].IndexArray[ItemIdx - 1])
                            {
                                if (UIF.CaseUtility > NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt)
                                {
                                    NIF.SumUt += (UIF.CaseUtility - NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt);
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt = UIF.CaseUtility;
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUtCase = UIF.VecIndex;
                                }

                                if (VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1] != 0 &&
                                    UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1] > NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU)
                                {
                                    NIF.SumPEU += ((UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1]) - NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU);
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU = UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1];
                                }

                                if (UIF.CaseUtility > MinUtil)
                                {
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].OverMinUtilCase++;
                                }

                                NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.push_back(UIF);
                                UIF.CaseUtility = 0;
                            }
                            else
                            {
                                if (UIF.CaseUtility > NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo[NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.size() - 1].CaseUtility)
                                {
                                    if (UIF.CaseUtility > MinUtil && NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo[NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.size() - 1].CaseUtility <= MinUtil)
                                    {
                                        NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].OverMinUtilCase++;
                                    }
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.pop_back();
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].L1_UtInfo.push_back(UIF);
                                }
                                if (UIF.CaseUtility > NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt)
                                {
                                    NIF.SumUt += (UIF.CaseUtility - NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt);
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUt = UIF.CaseUtility;
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqUtCase = UIF.VecIndex;
                                }
                                if (VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1] != 0 &&
                                    UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1] > NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU)
                                {
                                    NIF.SumPEU += ((UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1]) - NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU);
                                    NIF.L2_SeqInfo[NIF.L2_SeqInfo.size() - 1].SeqPEU = UIF.CaseUtility + VecDataBase[NodeUC.L2_SeqInfo[i].sid].RuArray[ItemIdx - 1];
                                }
                                UIF.CaseUtility = 0;
                            }
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }

        for (int a = 0; a < NIF.L2_SeqInfo.size(); a++)
        {
            for (int b = 0; b < NIF.L2_SeqInfo[a].L1_UtInfo.size(); b++)
            {
                int leaf_idx = NIF.L2_SeqInfo[a].L1_UtInfo[b].VecIndex;
                int c_pos = NIF.VecPattern.size() - 1;
                NIF.L2_SeqInfo[a].IdxPos.push_back({leaf_idx, {b, c_pos}});

                int tempIns = NIF.L2_SeqInfo[a].L1_UtInfo[b].IdxOfLastLevelIns;
                int tempSeq = NIF.L2_SeqInfo[a].IdxOfLastLevelSeq;
                int level = PatternPath.size() - 1;

                while (level >= 0)
                {
                    auto &trace_node = PatternPath[level].get();
                    int idx = trace_node.L2_SeqInfo[tempSeq].L1_UtInfo[tempIns].VecIndex;
                    int ancestor_c_pos = trace_node.VecPattern.size() - 1;
                    NIF.L2_SeqInfo[a].IdxPos.push_back({idx, {b, ancestor_c_pos}});

                    tempIns = trace_node.L2_SeqInfo[tempSeq].L1_UtInfo[tempIns].IdxOfLastLevelIns;
                    tempSeq = trace_node.L2_SeqInfo[tempSeq].IdxOfLastLevelSeq;
                    level--;
                }
            }
        }

        PatternPath.push_back(ref(NIF));


        if (NIF.SumUt >= MinUtil)
        {
            Cout_HUSPL3(NIF);
            S_ExtensionCounter++;
            REIHUSP_hiding(PatternPath);
        }
        HUSP(NIF);

        PatternPath.pop_back();
    }
}

void DBWriteBack(vector<SeqData> VecDataBase)
{
    string Filename = "Adjust_" + str_DBFile + "_Minutil_";
    Filename.append(to_string((int)MinUtil).append(".txt"));
    ofstream writeFile(Filename);
    if (writeFile.is_open())
    {
        for (int i = 1; i < VecDataBase.size(); i++)
        {
            int Tid = 1;
            for (int j = 0; j < VecDataBase[i].IndexArray.size(); j++)
            {
                if (VecDataBase[i].TidArray[j] != Tid)
                {
                    writeFile << "-1 ";
                    Tid++;
                }
                writeFile << VecDataBase[i].ItemArray[j];
                writeFile << "[";
                writeFile << VecDataBase[i].IuArray[j];
                writeFile << "] ";
            }
            writeFile << "-1 -2  SUtility:";
            writeFile << VecDataBase[i].RuArray[0] + VecDataBase[i].UtilityArray[0] << "\n";
        }
    }
}

void OutputResult(double time)
{
    string Filename = "Output_" + str_DBFile + "_Experiment_Minutil_";
    Filename.append(to_string((int)MinUtil).append(".txt"));
    ofstream writeFile(Filename);
    if (writeFile.is_open())
    {
        writeFile << "Database:" << str_DBFile << "\n";
        writeFile << "=====================" << "\n";
        writeFile << "Single item Counter : " << Single_ItemCounter << "\n";
        writeFile << "I_Extension Counter : " << I_ExtensionCounter << "\n";
        writeFile << "S_Extension Counter : " << S_ExtensionCounter << "\n";
        writeFile << "Total pattern Counter : " << Single_ItemCounter + I_ExtensionCounter + S_ExtensionCounter << "\n";
        writeFile << "Time consumption:" << time << "s" << "\n";
        writeFile << "Memory Usage:" << memoryMB << "(MB)" << "\n";
        writeFile << "=====================" << "\n";
    }
}

int main()
{
    SeqData SD;
    SD.sid = 0;
    VecDataBase.push_back(SD);

    L3_NodeInfo EmptyL3;
    Node_SingleItem.push_back(EmptyL3);

    ExternalUt.insert(make_pair(0, 0));
    cout << endl;

    //str_EuFile = "simple_utb.txt";
    //str_DBFile = "simple_db.txt";

    str_EuFile = "jzwpaper_utb.txt";
    str_DBFile = "jzwpaper_db.txt";

    //str_EuFile = "05.foodmart_ExternalUtility.txt";
    //str_DBFile = "05.foodmart.txt";

    //str_EuFile = "01.bible_ExternalUtility.txt";
    //str_DBFile = "01.bible.txt";
    
    //str_EuFile = "03.bms2_ExternalUtility.txt";
    //str_DBFile = "03.bms2.txt";

    //str_EuFile = "4_sign_ExternalUtility.txt";
    //str_DBFile = "4_sign.txt";

    //str_EuFile = "7.RSCS_ExternalUtility.txt";
    //str_DBFile = "7.RSCS.txt";   
     
    MinUtil = 0;

    cout << "*** (Hiding)Min utility = " << MinUtil << " ***" << endl;
    cout << "*** (Hiding)Database : " << str_DBFile << " ***" << endl;
    cout << "*** (Hiding)Eu : " << str_EuFile << " ***" << endl;

    Read_ExternalUt(str_EuFile);
    Read_Database(str_DBFile);
    BulidSingleItems(VecDataBase);

    /*double Percentage = 0.12;
    double Threshold = Percentage * 0.01;
    //MinUtil = 0;
    MinUtil = SumSWU * Threshold;

    //cout << MinUtil << "=" << SumSWU << "*" << Threshold << endl;
    cout << "*** (Hiding)Min utility = " << MinUtil << " ***" << endl;
    cout << "*** (Hiding)Database : " << str_DBFile << " ***" << endl;
    cout << "*** (Hiding)Eu : " << str_EuFile << " ***" << endl;*/


    double time = 0;
    clock_t start, end;
    start = clock();

    for (int i = 1; i < Node_SingleItem.size(); i++)
    {
        UpdateSingleItem(Node_SingleItem, stoi(Node_SingleItem[i].pattern));

        if (Node_SingleItem[i].SumUt >= MinUtil)
        {
            Single_ItemCounter++;
            SingleItem_Hiding(Node_SingleItem, stoi(Node_SingleItem[i].pattern));
        }

        PatternPath.push_back(ref(Node_SingleItem[i]));
        HUSP(Node_SingleItem[i]);
        PatternPath.clear();
    }

    end = clock();

    cout << endl;
    cout << "SWU : " << SumSWU << endl;
    cout << "Database:" << str_DBFile << endl;
    cout << "*** MinUtil = " << MinUtil << " ***" << endl;
    cout << "Single item Counter : " << Single_ItemCounter << endl;
    cout << "I_Extension Counter : " << I_ExtensionCounter << endl;
    cout << "S_Extension Counter : " << S_ExtensionCounter << endl;
    cout << "Total pattern Counter : " << Single_ItemCounter + I_ExtensionCounter + S_ExtensionCounter << endl;
    cout << endl;

    time = (double)(end - start) / CLOCKS_PER_SEC;

    PROCESS_MEMORY_COUNTERS pmc;
    if (GetProcessMemoryInfo(GetCurrentProcess(), &pmc, sizeof(pmc)))
    {
        memoryMB = pmc.WorkingSetSize / (1024.0 * 1024.0);
        cout << "Memory Usage (Working Set Size): " << memoryMB << " MB" << endl;
    }

    cout << "Time consumption:" << time << "s" << endl;
    cout << endl;

    cout << "Start Write out sanitized dataset ..." << endl;
    DBWriteBack(VecDataBase);

    cout << "Start Write out hiding experiment ..." << endl;
    OutputResult(time);
    cout << endl;

    cout << "******** Hiding HUSP process END ********" << endl;
    cout << endl;
    return 0;
}