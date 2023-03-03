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

int main(int argc, char* argv[]) {
	double start = getTime();
	ifstream in("xf.csv");
	string line;
	int ln = 0;
	vector <short> rowIndices;
	vector <short> columnIndices;
	vector <float> values;
	while (!in.eof()) {
		getline(in, line);
		if (!in.eof()) {
			vector <string> row = splitBy(line, ',');
			if (row.size() == 3) {
                rowIndices.push_back(atoi(row[0].c_str()));
                columnIndices.push_back(atoi(row[1].c_str()));
                values.push_back(atof(row[2].c_str()));
			}
			ln++;
		}
		if(ln % 100000 == 0) clog << ln << "\r";
	}
	in.close();
	clog << ln << endl;
	clog << getTime() - start << endl;

	ofstream out("bindata/matrix.bin", ios::binary);
	out.write((const char *)&rowIndices[0], rowIndices.size() * sizeof(short));
	out.write((const char *)&columnIndices[0], columnIndices.size() * sizeof(short));
	out.write((const char *)&values[0], values.size() * sizeof(float));
	out.close();

	return 0;
}

