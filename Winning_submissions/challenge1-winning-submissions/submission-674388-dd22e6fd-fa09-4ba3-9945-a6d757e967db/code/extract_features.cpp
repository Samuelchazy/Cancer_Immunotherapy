#include <bits/stdc++.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <dirent.h>
#define SSTR( x ) dynamic_cast< std::ostringstream & >( ( std::ostringstream() << std::dec << x ) ).str()


using namespace std;

double getTime() {
    timeval t;
    gettimeofday(&t,NULL);
    return 1e-6*t.tv_usec + t.tv_sec;
}
string trim(const string& str, const string& whitespace=" \t\r\n") {
    size_t strBegin = str.find_first_not_of(whitespace);
    if (strBegin == string::npos) return "";
    size_t strEnd = str.find_last_not_of(whitespace);
    size_t strRange = strEnd - strBegin + 1;
    return str.substr(strBegin, strRange);
}
vector <string> splitBy(const string& text, char by) {   // split string by by
	vector <string> vys;
	stringstream ss(text);
    string word;
    while(getline(ss,word,by)){
        vys.push_back(word);
    }
    return vys;
}

unordered_map <string, int> stateId = {
    {"progenitor", 0},
    {"effector", 1},
    {"terminal exhausted", 2},
    {"cycling", 3},
    {"other", 4}
};

struct Cell {
  public:
    int state;
    int perturbedGene;
    Cell(int st, int gn) : state(st), perturbedGene(gn) {}
};

int MINPERTURBEDCOUNT = -1;

const int CG = 80497219;
const int C = 28697;
const int G = 15077;
const int F = 14;
const int S = 5;
vector<short> rows(CG);
vector<short> cols(CG);
vector<float> vals(CG);
float matrix[G][C] = {0};
unordered_map <string, int> geneId;
vector<int> unperturbedIds;
vector<int> perturbedIds;
vector<int> perturbedIdsSp;

vector<Cell> cells;
vector<int> trainGenes;
vector<int> testGenes;

double mean(const vector<int>& rows, int gene) {
    double sum = 0;
    for (int row: rows) {
        sum += matrix[gene][row];
    }
    return rows.size() > 0 ? sum / rows.size() : 0;
}
double dev(const vector<int>& rows, const double& m, int gene) {
    double sum = 0;
    for (int row: rows){
        sum += (matrix[gene][row] - m) * (matrix[gene][row] - m);
    }
	return rows.size() > 0 ? sqrt(sum / rows.size()) : 0;
}
double skew(const vector<int>& rows, const double& m, const double& d, int gene) {
    double sum = 0;
    for (int row: rows){
        sum += (matrix[gene][row] - m) * (matrix[gene][row] - m) * (matrix[gene][row] - m);
    }
	return rows.size() > 0 && d != 0 ? (sum / rows.size()) / (d * d * d) : 0;
}
vector<int> nonzeros(vector<int>& rows, int gene) {
    vector<int> ans;
    for (int row: rows) {
        if (matrix[gene][row] != 0) ans.push_back(row);
    }
    return ans;
}

vector<double> generateFeatures(int gene) {
    vector<double> features(F, 0);
    int f = 0;
    vector<int> nonzeroRows = nonzeros(unperturbedIds, gene);
    double m = mean(unperturbedIds, gene);
    double mNnz = mean(nonzeroRows, gene);
    double d = dev(unperturbedIds, m, gene);
    double dNnz = dev(nonzeroRows, mNnz, gene);
    features[f++] = nonzeroRows.size() / (double) unperturbedIds.size();
    features[f++] = m;
    features[f++] = mNnz;
    features[f++] = d;
    features[f++] = dNnz;
    features[f++] = skew(unperturbedIds, m, d, gene);
    features[f++] = skew(nonzeroRows, mNnz, dNnz, gene);

    //vector<int> temp = nonzeros(perturbedIdsSp, gene);
    //features[14] = features[0] - temp.size() / (double) perturbedIdsSp.size();
    //features[15] = m - mean()

    nonzeroRows = nonzeros(perturbedIds, gene);
    m = mean(perturbedIds, gene);
    mNnz = mean(nonzeroRows, gene);
    d = dev(perturbedIds, m, gene);
    dNnz = dev(nonzeroRows, mNnz, gene);
    features[f++] = nonzeroRows.size() / (double) perturbedIds.size();
    features[f++] = m;
    features[f++] = mNnz;
    features[f++] = d;
    features[f++] = dNnz;
    features[f++] = skew(perturbedIds, m, d, gene);
    features[f++] = skew(nonzeroRows, mNnz, dNnz, gene);
    return features;
}

int perturbedCount;
vector<double> generateGroundTruth(int gene) {
    vector<double> gt(S, 0);
    perturbedCount = 0;
    for (Cell& cell: cells) {
        if (cell.perturbedGene == gene) {
            perturbedCount++;
            gt[cell.state] += 1;
        }
    }
    if (perturbedCount > 0) {
        for (int s = 0; s < S; s++) gt[s] /= perturbedCount;
    }
    return gt;
}


void writeVector(ofstream& out, vector<double>& data, string ending = "\n") {
    for (int i = 0; i < data.size(); i++) {
        if (i != 0) out << ",";
        out << data[i];
    }
    out << ending;
}


int main(int argc, char* argv[]) {
	double start = getTime();
	string line;

	ifstream in("var.csv");
	getline(in, line);
	int ln = 0;
	while (!in.eof()) {
		getline(in, line);
		if (!in.eof()) {
			string geneName = trim(line);
			if (geneName != "") {
                geneId[geneName] = ln++;
			}
		}
	}
	in.close();
	clog << ln << " gene names loaded" << endl;
	geneId["Unperturbed"] = -1;
    geneId["Fzd1"] = -2;
    geneId["P2rx7"] = -3;

    set<int> trainGenesSet;
	in.open("obs.csv");
    getline(in, line);
	ln = 0;
	while (!in.eof()) {
		getline(in, line);
		if (!in.eof()) {
			vector <string> row = splitBy(line, ',');
			if (row.size() == 5) {
                string state = row[2];
                string perturbedGene = row[3];
                if (stateId.count(state) == 0 || geneId.count(perturbedGene) == 0) {
                    cerr << ln << " incorrect cell info" << endl;
                    cerr << state << " " << perturbedGene << endl;
                    cerr << stateId.count(state) << " " << geneId.count(perturbedGene) << endl;
                }
                cells.push_back(Cell(stateId[state], geneId[perturbedGene]));
                if (geneId[perturbedGene] >= 0) trainGenesSet.insert(geneId[perturbedGene]);
                if (geneId[perturbedGene] == -1) unperturbedIds.push_back(ln);
                else perturbedIds.push_back(ln);
                if (perturbedGene == "Dvl1") perturbedIdsSp.push_back(ln);
                ln++;
			}
		}
	}
	in.close();
	//clog << perturbedIds.size() << " " << perturbedIdsSp.size() << endl;
	clog << ln << " cells loaded" << endl;
	trainGenes.assign(trainGenesSet.begin(), trainGenesSet.end());

	start = getTime();

	in.open("bindata/matrix.bin", ios::binary);
	in.read((char *)&rows[0], CG * sizeof(short));
	in.read((char *)&cols[0], CG * sizeof(short));
	in.read((char *)&vals[0], CG * sizeof(float));
	in.close();
    clog << "sparse matrix loaded in " << getTime() - start << " sec." << endl;

	start = getTime();
	for (int i = 0; i < CG; i++) matrix[cols[i]][rows[i]] = vals[i];
	clog << "full matrix created in " << getTime() - start << " sec." << endl;

	ofstream out("training.csv");
	out.precision(16);

    int trainSize = 0;
    for (int g: trainGenes) {
        vector <double> features = generateFeatures(g);
        vector <double> groundTruth = generateGroundTruth(g);
        if (perturbedCount > MINPERTURBEDCOUNT) {
            writeVector(out, features, ",");
            writeVector(out, groundTruth, "\n");
            trainSize++;
        }
    }
    clog << trainSize << " train samples";
    out.close();
	out.open("testing.csv");
	testGenes = {geneId["Aqr"], geneId["Bach2"], geneId["Bhlhe40"], geneId["Ets1"], geneId["Fosb"], geneId["Mafk"], geneId["Stat3"]};
    for (int g: testGenes) {
        vector <double> features = generateFeatures(g);
        writeVector(out, features, "\n");
    }


	return 0;
}

