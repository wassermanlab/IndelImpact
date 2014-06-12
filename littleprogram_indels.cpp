#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include "raven_scoring.h"
//#include <EXTERN.h>
//#include <perl.h>

using namespace std;
/*
string exec(char* cmd) {
	FILE* pipe = popen(cmd, "r");
	if (!pipe) return "ERROR";
	char buffer[128];
	string result = "";
	while (!feof(pipe)) {
		if (fgets(buffer, 128, pipe) != NULL) {
			result += buffer;
		}
	}
	pclose(pipe);
	return result;
}
*/
vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while(getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


vector<string> split(const string &s, char delim) {
    vector<string> elems;
    return split(s, delim, elems);
}

void all_scores (double* matrix, char* seq, int start, int end, double* allscarray, int length) {
	double forward_score;
	double back_score;
	for (int i = start; i <= end - length + 1; ++i) {
		forward_score = 0.0;
		back_score = 0.0;
		for (int j = 0; j < length; ++j) {
			if (TRANS[seq[i+j-1]] == 4) {
				double forward_total_score = 0;
				double backward_total_score = 0;
				for (int k = 0; k < 4; ++k) {
					forward_total_score += matrix[4*j+k];
					backward_total_score += matrix[4*(length-1-j) + k];
				}
				forward_score += forward_total_score/4;
				back_score += backward_total_score/4;
			} else {
				forward_score += matrix[4*j+TRANS[seq[i+j-1]]];
				back_score += matrix[4*(length - 1 -j) + 3-TRANS[seq[i+j-1]]];
			}
		}
		allscarray[2*i] = forward_score;
		allscarray[2*i+1] = back_score;
	}
}

int main(int argc, char** argv) {
	//my_perl = perl_alloc();
	//perl_construct(my_perl);
	//perl_parse(my_perl, NULL, argc, argv, env);
	double scratch[1024];
	char garbage[128];
	double *allscarr;
	int seq_len = strlen(argv[2]);
	if (seq_len < 4096) {
		seq_len = 4096;
	}
	allscarr = new double[2*seq_len+100];
	double pwm[1024];
	int length;
	int done;
	int num_counts;
	FILE* fp;
	vector<string> filenames = split(string(argv[1]), ',');
	vector<string> starts = split(string(argv[3]), ',');
	vector<string> ends = split(string(argv[4]), ',');
	vector<string> names = split(string(argv[5]), ',');
	vector<string>::iterator itstart = starts.begin();
	vector<string>::iterator itend = ends.begin();
	vector<string>::iterator itname = names.begin();
	string outpath = "";
	outpath += (string)argv[8];
	outpath += "/";
	outpath += (string)argv[6];
	outpath += (string)argv[7];
	outpath += (string)argv[9];
	outpath += ".txt";
	ofstream outfile(outpath.c_str());
	for (vector<string>::iterator it = filenames.begin(); it != filenames.end(); ++it) {
	done = 0;
	if (fp=fopen((char*) ((string)*it).c_str(),"r")) {
		fscanf(fp, "%s", garbage); 
      for (num_counts=0; !done && num_counts<MAXCOUNTS; ++num_counts )
      {
         if ( fscanf(fp,"%lf,%*c",scratch+num_counts) == EOF )
            done = 1;
      }
	fclose(fp);
   
      length = num_counts/4;
      for (int pos=0; pos<length; ++pos )
      {
         for (int nt=0; nt<4; ++nt )
         {
            pwm[4*pos + nt] = scratch[length*nt + pos];
         }
      }
	} else {cout << "Failed to open file" << endl;}
        
		all_scores(pwm, argv[2], atoi(((string)*itstart).c_str()), atoi(((string)*itend).c_str()), allscarr, length);
		string outstr = (string)*itname;
		for (int i = 0; i <= 2*(atoi(((string)*itend).c_str()) - length)+3; ++i) {
			outstr += ",";
			ostringstream strs;
			if (i < 2*atoi(((string)*itstart).c_str())) {
				strs << 0;
			} else strs << allscarr[i];
			outstr += strs.str();
		}
//		char namebuf[64];
//		outstr += '\t';
//		snprintf(namebuf, sizeof(namebuf),"%d",names.size());
//                outstr += namebuf;
		outstr += '\n';
		outfile << outstr;
		++itstart;
		++itend;
		++itname;
	}
	delete [] allscarr;
	outfile.close();
	
	return 0;
}
