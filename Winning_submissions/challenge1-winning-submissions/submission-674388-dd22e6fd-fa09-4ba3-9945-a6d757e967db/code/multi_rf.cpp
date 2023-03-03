#include <bits/stdc++.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <dirent.h>
#define SSTR( x ) dynamic_cast< std::ostringstream & >( ( std::ostringstream() << std::dec << x ) ).str()
#define SSTRF( x ) dynamic_cast< std::ostringstream & >( ( std::ostringstream() << fixed << setprecision(2) << x ) ).str()
#define GCOST(x) ( sqrt((x)*(1 - (x)) ) )

using namespace std;

const int RESULTDIM = 5;


vector <int> FEAT_LIST = {0,1,2,3,4,5,6,7,8,9,10,11,12,13};

const int FEATURES = 14;
const int ALLFEATURES = 14;

int THREADS = 4;

int TREES = 5000;
int FEATURETRY = 3;
const int MAXLEVEL = 50;
bool USEMSE = true;
const int CLUSTERS = 2;
int MINNODEMSE = 10;
int MINNODE = 5;
const int MAXSAMPLESIZE = 3000000; //10000000;
int BORDER = 1000000000;
const double EPS = 1e-10;

int NTRAIN = 64;
int NTEST = 7;
long long SAMPLES = NTRAIN;
vector <double> FEAT(SAMPLES * FEATURES);
vector <vector<double> > RESULT(NTRAIN, vector<double>(RESULTDIM, 0));
bool training;
mutex mtxG;

vector <int> featureScore;
vector <int> featureScoreC;

double getTime() {
    timeval t;
    gettimeofday(&t,NULL);
    return 1e-6*t.tv_usec + t.tv_sec;
}
// random generator
//unsigned long long nowRand = 1;
vector <unsigned long long> nowRand;
void seedBig(unsigned long long seed, int th){
	nowRand[th] = seed;
}
unsigned long long randBig(int th){
	nowRand[th] = ((nowRand[th] * 6364136223846793005ULL + 1442695040888963407ULL) >> 1);
	return nowRand[th];
}
string int2len(int v, int l){
	string ret = SSTR(v);
	int dig = floor(log10(v + 0.1)) + 1;
	while(ret.length() < l) ret = " " + ret;
	return ret;
}
string string2len(string ret, int l){
	while(ret.length() < l) ret += " ";
	return ret;
}
string trim(const string& str,const string& whitespace=" \t\r\n"){
    size_t strBegin = str.find_first_not_of(whitespace);
    if (strBegin == string::npos) return "";
    size_t strEnd = str.find_last_not_of(whitespace);
    size_t strRange = strEnd - strBegin + 1;
    return str.substr(strBegin, strRange);
}
vector <string> splitBy(const string& text, char by){   // split string by by
	vector <string> vys;
	stringstream ss(text);
    string word;
    while(getline(ss,word,by)){
        vys.push_back(word);
    }
    return vys;
}
vector <int> selectFeatures(int th){
	set <int> temp;
	while(temp.size() != FEATURETRY){
		temp.insert(randBig(th) % FEATURES);
	}
	vector<int> result(temp.begin(), temp.end());
	//random_shuffle(result.begin(),result.end());
	return result;
}
class Node{
  public:
	int left;
	int right;
	int feature;
	double value;
	int level;
	int counts[CLUSTERS];
	int total;
	double sumX[RESULTDIM];
	Node(){
		left = -1;
		right = -1;
		feature = -1;
		value = 0;
		level = -1;
		for(int i = 0; i < CLUSTERS; i++) counts[i] = 0;
		for(int i = 0; i < RESULTDIM; i++) sumX[i] = 0;
		total = 0;
	}
	Node(int lev, const vector <int>& cou, int tot){
		left = -1;
		right = -1;
		feature = -1;
		value = 0;
		level = lev;
		copy(cou.begin(), cou.end(), counts);
		for(int i = 0; i < RESULTDIM; i++) sumX[i] = 0;
		total = tot;
	}
	Node(int lev, int tot, const vector <double>& sumx){
		left = -1;
		right = -1;
		feature = -1;
		value = 0;
		level = lev;
		total = tot;
		copy(sumx.begin(), sumx.end(), sumX);
	}
};
class Tree{
  public:
	vector <Node> node;
	vector <double> resultAtNode(int nodeId, double* source) const{
		if(node[nodeId].left == -1){
            vector <double> ret(begin(node[nodeId].sumX), end(node[nodeId].sumX));
            for (int i = 0; i < ret.size(); i++) ret[i] /= node[nodeId].total;
			return ret;
		}
		if(*(source + node[nodeId].feature) <= node[nodeId].value) return resultAtNode(node[nodeId].left, source);
		return resultAtNode(node[nodeId].right, source);
	}
	vector<double> assignResult(double* source) const{
		return resultAtNode(0, source);
	}
	vector<double> OOBresultAtNode(int nodeId, int i) const{
		if(node[nodeId].left == -1){
            vector <double> ret(begin(node[nodeId].sumX), end(node[nodeId].sumX));
            for (int j = 0; j < ret.size(); j++) ret[j] /= node[nodeId].total;
			return ret;
		}
		if(FEAT[node[nodeId].feature * SAMPLES + i] <= node[nodeId].value) return OOBresultAtNode(node[nodeId].left, i);
		return OOBresultAtNode(node[nodeId].right, i);
	}
	vector<double> assignOOBResult(int i) const{
		return OOBresultAtNode(0, i);
	}
	double ratioAtNode(int nodeId, double* source) const{
		if(node[nodeId].left == -1){
			double sum = 0;
			for(int cluster = 0; cluster < CLUSTERS; cluster++){
				sum += node[nodeId].counts[cluster] * cluster;
			}
			return sum / node[nodeId].total;
		}
		if(*(source + node[nodeId].feature) <= node[nodeId].value) return ratioAtNode(node[nodeId].left, source);
		return ratioAtNode(node[nodeId].right, source);
	}
	double assignRatio(double* source) const{
		return ratioAtNode(0, source);
	}
	double OOBratioAtNode(int nodeId, int id) const{
		if(node[nodeId].left == -1){
			double sum = 0;
			for(int cluster = 0; cluster < CLUSTERS; cluster++){
				sum += node[nodeId].counts[cluster] * cluster;
			}
			return sum / node[nodeId].total;
		}
		if(FEAT[node[nodeId].feature * SAMPLES + id] <= node[nodeId].value) return OOBratioAtNode(node[nodeId].left, id);
		return OOBratioAtNode(node[nodeId].right, id);
	}
	double assignOOBRatio(int id) const{
		return OOBratioAtNode(0, id);
	}
	/*void divideNode(int nodeIndex, vector <int>& sample, int th){
		int n = sample.size();
		int nonzero = 0;
		for(int i = 0; i < CLUSTERS; i++){
			if(node[nodeIndex].counts[i] > 0) nonzero++;
			if(nonzero > 1) break;
		}
		if(node[nodeIndex].level < MAXLEVEL - 1 && nonzero > 1 && node[nodeIndex].total > MINNODE){
	    	vector <int> feaID = selectFeatures(th);
	    	double minCost = 1e30;
	    	int bestF = -1;
	    	double bestValue = 0;
	    	int bestI = 0;
	    	vector <int> bestC1(CLUSTERS, 0);
	    	int bestTotalL = 0;
			int SORTBY = 0;
			long long SORTBYSAMPLES = 0;
			for(int f = 0; f < FEATURETRY; f++){
				SORTBY = feaID[f];
				SORTBYSAMPLES = SORTBY * SAMPLES;
				sort(sample.begin(), sample.end(), [&](int aa, int bb){return FEAT[SORTBYSAMPLES + aa] < FEAT[SORTBYSAMPLES + bb];});

				//int FSIZ=FSIZE[fi*FEATURES+SORTBY];
				//int bucket[FSIZ]={};
				//for(int j=0;j<n;j++){
				//	bucket[FEAT[SORTBYSAMPLES+sample[j]]]++;
				//}
				//int cumSum[FSIZ+1];
				//cumSum[0]=0;
				//for(int i=0;i<FSIZ;i++) cumSum[i+1]=cumSum[i]+(bucket[i]--);
				//vector <int> temp=sample;
				//for(int i=0;i<n;i++){
				//	int v=FEAT[SORTBYSAMPLES+temp[i]];
				//	sample[cumSum[v]+(bucket[v]--)]=temp[i];
				//}

				vector <int> c1(CLUSTERS, 0);
	    		int totalL = 0;
	    		for(int i = 0; i < n-1; i++){
	    			c1[int(RESULT[sample[i]])]++;
	    			totalL++;
	    			if(FEAT[SORTBYSAMPLES + sample[i+1]] > FEAT[SORTBYSAMPLES + sample[i]]){
	    			    double costL = 0.0;
						double costR = 0.0;
						for(int cl = 0; cl < CLUSTERS; cl++){
							costL += GCOST(c1[cl] / static_cast<double>(totalL));
							costR += GCOST((node[nodeIndex].counts[cl] - c1[cl]) / static_cast<double>(n - totalL));
						}
						double cost = (totalL * costL + (n - totalL) * costR) / n;
						if(cost < minCost && i >= n/BORDER && i < n - n/BORDER){
	    			    	minCost = cost;
	    			    	bestF = feaID[f];
	    			    	bestValue = FEAT[SORTBYSAMPLES + sample[i]];
							bestI = i;
	    			    	bestC1 = c1;
	    			    	bestTotalL = totalL;
	    			    }
					}
	    		}
	    	}
	    	if(bestF >= 0){
				//if(bestTotalL == node[nodeIndex].total) cerr << node[nodeIndex].total << " " << n << endl;
				mtxG.lock();
				featureScore[bestF] += n;
				featureScoreC[bestF]++;
				mtxG.unlock();
		    	vector <int> sampleLeft; sampleLeft.reserve(bestI + 1);
		    	vector <int> sampleRight; sampleRight.reserve(n - bestI - 1);
		    	SORTBYSAMPLES = bestF * SAMPLES;
				for(int i = 0; i < n; i++){
		    		if(FEAT[SORTBYSAMPLES + sample[i]] <= bestValue){
						sampleLeft.push_back(sample[i]);
					}
		    		else sampleRight.push_back(sample[i]);
		    	}
		    	//if(sampleLeft.size() != bestTotalL){
				//	cout << "!" << sampleLeft.size() << " " << bestTotalL << endl;
				//	cout << bestValue << " " << bestI << endl;
				//	for(int i = 0; i < n; i++){
				//		cout << FEAT[SORTBYSAMPLES + sample[i]] << endl;
			    //	}
			    //	cout << " ------------------------------------- " << endl;
				//}
		        node[nodeIndex].feature = bestF;
		    	node[nodeIndex].value = bestValue;
		    	node.push_back(Node(node[nodeIndex].level + 1, bestC1, bestTotalL));
		    	//if(bestTotalL <= 0) cerr << "L" << bestTotalL << endl;
		    	node[nodeIndex].left = node.size() - 1;
		    	vector <int> c2(CLUSTERS, 0);
		    	for(int i = 0; i < CLUSTERS; i++){
		    		c2[i] = node[nodeIndex].counts[i] - bestC1[i];
		    	}
		    	node.push_back(Node(node[nodeIndex].level + 1, c2, node[nodeIndex].total - bestTotalL));
		    	//if(node[nodeIndex].total - bestTotalL <= 0) cerr << "R" << node[nodeIndex].total - bestTotalL << endl;

		    	node[nodeIndex].right = node.size() - 1;
			    divideNode(node[nodeIndex].left, sampleLeft, th);
				divideNode(node[nodeIndex].right, sampleRight, th);
			}
		}
	}*/
	void divideNodeMSE(int nodeIndex, vector <int>& sample, int th){
		int n = sample.size();
		if(node[nodeIndex].level < MAXLEVEL-1 && node[nodeIndex].total > MINNODEMSE){
			vector <int> feaID = selectFeatures(th);
			double minCost = 1e30;
			int bestF = -1;
			double bestValue = 0;
			int bestI = 0;
			vector<double> bestSumXL(RESULTDIM, 0);
			int bestTotalL = 0;
			int SORTBY = 0;
			long long SORTBYSAMPLES = 0;
			for(int f = 0; f < FEATURETRY; f++){
				SORTBY = feaID[f];
				SORTBYSAMPLES = SORTBY*SAMPLES;
				sort(sample.begin(), sample.end(), [&](int aa, int bb){return FEAT[SORTBYSAMPLES+aa] < FEAT[SORTBYSAMPLES+bb];});

				//sortFast(sample);

				/*int bucket[FSIZE[SORTBY]]={};
				for(int j=0;j<n;j++){
					bucket[FEAT[SORTBYSAMPLES+sample[j]]]++;
				}
				int cumSum[FSIZE[SORTBY]+1];
				cumSum[0]=0;
				for(int i=0;i<FSIZE[SORTBY];i++) cumSum[i+1]=cumSum[i]+(bucket[i]--);
				vector <int> temp=sample;
				for(int i=0;i<n;i++){
					int v=FEAT[SORTBYSAMPLES+temp[i]];
					sample[cumSum[v]+(bucket[v]--)]=temp[i];
				}*/


				vector<double> sumXL(RESULTDIM, 0);
				int totalL = 0;
	    		for(int i = 0; i < n - 1; i++){
					for (int j = 0; j < RESULTDIM; j++) sumXL[j] += RESULT[sample[i]][j];
					totalL++;
					if(FEAT[SORTBYSAMPLES + sample[i+1]] > FEAT[SORTBYSAMPLES + sample[i]]){
						double cost = 0;
						for (int j = 0; j < RESULTDIM; j++) {
                            cost += -sumXL[j]*sumXL[j]/totalL - (node[nodeIndex].sumX[j]-sumXL[j]) * (node[nodeIndex].sumX[j]-sumXL[j]) / (n-totalL);
						}
						if(cost < minCost && i >= n/BORDER && i < n - n/BORDER){
	    			    	minCost = cost;
	    			    	bestF = SORTBY;
							bestValue = FEAT[SORTBYSAMPLES+sample[i]];
							bestI = i;
	    			    	bestSumXL = sumXL;
	    			    	bestTotalL = totalL;
	    			    }
					}
	    		}
	    	}
	    	if(bestF >= 0){
	    		mtxG.lock();
				featureScore[bestF] += n;
	    		featureScoreC[bestF]++;
	    		mtxG.unlock();
				vector <int> sampleLeft; sampleLeft.reserve(bestI + 1);
		    	vector <int> sampleRight; sampleRight.reserve(n - bestI - 1);
				SORTBYSAMPLES = bestF*SAMPLES;
				for(int i = 0; i < n; i++){
					if(FEAT[SORTBYSAMPLES + sample[i]] <= bestValue){
						sampleLeft.push_back(sample[i]);
					}
		    		else sampleRight.push_back(sample[i]);
		    	}
		        node[nodeIndex].feature = bestF;
		    	node[nodeIndex].value = bestValue;
		    	node.push_back(Node(node[nodeIndex].level+1, bestTotalL, bestSumXL));
		    	vector <double> bestSumXR(begin(node[nodeIndex].sumX), end(node[nodeIndex].sumX));
		    	for (int j = 0; j < RESULTDIM; j++) bestSumXR[j] -= bestSumXL[j];
		    	node[nodeIndex].left = node.size() - 1;
		    	node.push_back(Node(node[nodeIndex].level+1, node[nodeIndex].total-bestTotalL, bestSumXR));
		    	node[nodeIndex].right = node.size() - 1;

				divideNodeMSE(node[nodeIndex].left,sampleLeft, th);
				divideNodeMSE(node[nodeIndex].right,sampleRight, th);
			}
		}
	}
	Tree(){
	}
	int toStream(ofstream &out){
		int size = node.size();
		out.write((const char *)&size, sizeof(int));
		for(int i = 0; i < size; i++){
			out.write((const char *)&node[i], sizeof(Node));
		}
		return 0;
	}
	Tree(ifstream& in){
		int size;
		in.read((char *)&size, sizeof(int));
		node.resize(size);
		for(int i = 0; i < size; i++){
			in.read((char *)&node[i], sizeof(Node));
		}
	}
};
int RFtoFile(vector <Tree>& rf, string fileName){
	ofstream out(fileName.c_str(), ios::binary);
	int trees = rf.size();
	out.write((const char *)&trees, sizeof(int));
	for(int i = 0; i < trees; i++){
		rf[i].toStream(out);
	}
	out.close();
	return 0;
}
vector <Tree> RFfromFile(string fileName){
	vector <Tree> rf;
	ifstream in(fileName.c_str(), ios::binary);
	int trees;
	in.read((char *)&trees, sizeof(int));
	for(int i = 0; i < trees; i++){
		rf.push_back(Tree(in));
	}
	in.close();
	return rf;
}
double forestAssignResultClassification(const vector <Tree>& tree, double* source){
	double result = 0;
	for(int t = 0; t < tree.size(); t++){
		result += tree[t].assignRatio(source);
	}
	return result / tree.size();
}
vector<double> forestAssignResultMSE(const vector <Tree>& tree, double* source){
	vector<double> result(RESULTDIM, 0);
	for(int t = 0; t < tree.size(); t++){
        vector<double> temp = tree[t].assignResult(source);
		for (int j = 0; j < RESULTDIM; j++) result[j] += temp[j];
	}
	for (int j = 0; j < RESULTDIM; j++) result[j] /= tree.size();
	return result;
}
vector<double> forestAssignOOBResultMSE(const vector <Tree>& tree, const vector <int>& indices, int i){
	vector<double> result(RESULTDIM, 0);
	for(int t = 0; t < indices.size(); t++){
        vector<double> temp = tree[indices[t]].assignOOBResult(i);
		for (int j = 0; j < RESULTDIM; j++) result[j] += temp[j];
	}
	for (int j = 0; j < RESULTDIM; j++) result[j] /= indices.size();
	return result;
}
/*Tree buildTree(int n, int th){
	//int n = SAMPLES;
	int ns = min(n, MAXSAMPLESIZE);
	vector <int> sample;
	sample.resize(ns);
	Tree tree;
	tree.node.resize(1, Node());
	tree.node[0].level = 0;
	for(int i = 0; i < ns; i++){
		sample[i] = randBig(th) % n;
		tree.node[0].counts[int(RESULT[sample[i]])]++;
	}
	tree.node[0].total = ns;
	tree.divideNode(0, sample, th);
	return tree;
}*/
unordered_set<int> treeSamples;
Tree buildTreeMSE(int n, int th){
	//int n = SAMPLES;
	int ns = min(n, MAXSAMPLESIZE);
	vector <int> sample;
	sample.resize(ns);
	Tree tree;
	tree.node.resize(1, Node());
	tree.node[0].level = 0;
	//tree.node[0].sumX = 0;
	treeSamples.clear();
	for(int i = 0; i < ns; i++){
		sample[i] = randBig(th) % n;
		treeSamples.insert(sample[i]);
		for (int j = 0; j < RESULTDIM; j++)	{
            tree.node[0].sumX[j] += RESULT[sample[i]][j];
		}
	}
	tree.node[0].total = ns;
	tree.divideNodeMSE(0, sample, th);
	return tree;
}

void writeVector(ostream& out, vector<double>& data, string ending = "\n") {
    for (int i = 0; i < data.size(); i++) {
        if (i != 0) out << ",";
        out << data[i];
    }
    out << ending;
}

int main(int argc, char* argv[]){
	THREADS = 1;

	nowRand.resize(THREADS);
	for(int i = 0; i < THREADS; i++) nowRand[i] = i + 1;

	vector <double> trainMean(RESULTDIM, 0);


	ifstream in("training.csv");
	int ln = 0;
	string line;
	while (!in.eof()) {
		getline(in, line);
		if (!in.eof()) {
			vector <string> row = splitBy(line, ',');
			if (row.size() == ALLFEATURES + RESULTDIM) {
                for (int i = 0; i < FEATURES; i++) {
                    FEAT[SAMPLES * i + ln] = atof(row[FEAT_LIST[i]].c_str());
                    //if (i >= 7) FEAT[SAMPLES * i + ln] -= FEAT[SAMPLES * (i-7) + ln];
                }
                for (int i = 0; i < RESULTDIM; i++) {
                    RESULT[ln][i] = atof(row[ALLFEATURES + i].c_str());
                    trainMean[i] += RESULT[ln][i];
                }
                ln++;
			}
		}
	}
	in.close();
    for (int i = 0; i < RESULTDIM; i++) {
        trainMean[i] /= ln;
    }

    double start = getTime();


	featureScore.clear();
	featureScoreC.clear();
	featureScore.resize(FEATURES, 0);
	featureScoreC.resize(FEATURES, 0);
	ofstream pred("prediktors.txt");
    vector <Tree> randomForest(TREES);
	vector <vector <int> > sampleOOBTrees(SAMPLES);
	vector <vector <int> > sampleIOBTrees(SAMPLES);
	for(int j = 0; j < TREES; j++){
        randomForest[j] = buildTreeMSE(SAMPLES, 0);
		if(j % 1 == 0) clog << j + 1 << " trees done...\r";
		for (int i = 0; i < SAMPLES; i++) {
            if(treeSamples.count(i) == 0) sampleOOBTrees[i].push_back(j);
            else sampleIOBTrees[i].push_back(j);
		}
    }
	clog << endl;
	clog << "  Average per tree: " << (getTime() - start) / TREES << " sec." << endl;

	vector <pair <int,int> > stat;
	for(int i = 0; i < FEATURES; i++) stat.push_back(make_pair(featureScore[i], i));
	sort(stat.begin(), stat.end());
	int len = log10(stat.back().first + 0.1) + 2;
	vector <pair <int,int> > statC;
	for(int i = 0; i < FEATURES; i++) statC.push_back(make_pair(featureScoreC[i], i));
	sort(statC.begin(),statC.end());
	int lenC = log10(statC.back().first + 0.1) + 2;
	for(int i = FEATURES - 1; i >= 0; i--){
		pred << int2len(stat[i].first, len) << " " << string2len(SSTR(stat[i].second), 3) << "   |   " << int2len(statC[i].first, lenC) << " " << statC[i].second << endl;
	}
	pred.close();

	//RFtoFile(randomForest, "bindata/rf.dat");
	//vector <Tree> randomForest = RFfromFile("bindata/rf" + part + ".dat");

	vector <vector<double> > resultOOB(NTRAIN, vector<double>(RESULTDIM));
	vector <vector<double> > resultIOB(NTRAIN, vector<double>(RESULTDIM));
	double errOOB = 0;
	double errIOB = 0;
	vector<double> errOOBs(NTRAIN, 0);
	vector<double> errMeans(NTRAIN, 0);
    ofstream outoob("training_oob_result.csv");
    for (int i = 0; i < NTRAIN; i++) {
        resultOOB[i] = forestAssignOOBResultMSE(randomForest, sampleOOBTrees[i], i);
        resultIOB[i] = forestAssignOOBResultMSE(randomForest, sampleIOBTrees[i], i);
        for (int j = 0; j < RESULTDIM; j++) {
            errOOBs[i] += fabs(resultOOB[i][j] - RESULT[i][j]);
            errMeans[i] += fabs(trainMean[j] - RESULT[i][j]);
            errIOB += fabs(resultIOB[i][j] - RESULT[i][j]);
        }
        errOOB += errOOBs[i];
        writeVector(outoob, RESULT[i]);
        writeVector(outoob, resultOOB[i]);
        outoob << "\n";
    }
    errOOB /= NTRAIN;
    errIOB /= NTRAIN;
	cerr << "OOB loss: " << errOOB << endl;
	cerr << "IOB loss: " << errIOB << endl;
	outoob.close();
	outoob.open("oob.csv");
	writeVector(outoob, errMeans);
	writeVector(outoob, errOOBs);
	outoob.close();

    in.open("testing.csv");
	ln = 0;
    vector <vector<double> > feat(NTEST, vector<double>(FEATURES));
	while (!in.eof()) {
		getline(in, line);
		if (!in.eof()) {
			vector <string> row = splitBy(line, ',');
			if (row.size() == ALLFEATURES) {
                for (int i = 0; i < FEATURES; i++) {
                    feat[ln][i] = atof(row[FEAT_LIST[i]].c_str());
                    //if (i >= 7) feat[ln][i] -= feat[ln][i-7];
                }
                ln++;
			}
		}
	}
	in.close();

    ofstream out("validation_output.csv");
    out << "gene,a_i,b_i,c_i,d_i,e_i" << endl;
    vector<string> geneNames = {"Aqr", "Bach2", "Bhlhe40", "Ets1", "Fosb", "Mafk", "Stat3"};
    vector <vector<double> > result(NTEST, vector<double>(RESULTDIM));
    for (int i = 0; i < NTEST; i++) {
        result[i] = forestAssignResultMSE(randomForest, &feat[i][0]);
        out << geneNames[i] << ",";
        writeVector(out, result[i]);
        if (i == 2) {
            out.close();
            out.open("test_output.csv");
            out << "gene,a_i,b_i,c_i,d_i,e_i" << endl;
        }
    }
    out.close();



	return 0;
}


