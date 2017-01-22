#include "gene.h"
#include "editdistance.h"
#include <iostream>
#include <fstream>
#include "util.h"
#include <string.h>
#include <sstream>
#include "globalsettings.h"

Gene::Gene(string name, string chr, int start, int end){
	mName = name;
    mChr = chr;
    mStart = start;
    mEnd = end;
}

Gene::Gene(const Gene& other) {
    mName = other.mName;
    mChr = other.mChr;
    mStart = other.mStart;
    mEnd = other.mEnd;
    mExons = other.mExons;
}

Gene::Gene() {
    mName = "invalid";
    mChr = "invalid";
    mStart = 0;
    mEnd = 0;
}

bool Gene::valid() {
    return mName != "invalid" && mStart != 0 && mEnd != 0;
}

void Gene::addExon(Exon exon) {
    mExons.push_back(exon);
}

void Gene::addExon(int id, int start, int end) {
    Exon exon;
    exon.id=id;
    exon.start=start;
    exon.end=end;
    mExons.push_back(exon);
}

void Gene::print() {
    cout<<mName<<","<<mChr<<":"<<mStart<<"-"<<mEnd<<endl;
    for(int i=0;i<mExons.size();++i){
        cout<<mExons[i].id<<","<<mExons[i].start<<","<<mExons[i].end<<endl;
    }
}

Gene Gene::parse(string linestr) {
    vector<string> splitted;
    split(linestr, splitted, ",");
    if(splitted.size()<2)
        return Gene();
    string name = trim(splitted[0].substr(1, splitted[0].length()-1));

    vector<string> chrPos;
    split(splitted[1], chrPos, ":");
    if(chrPos.size()<2)
        return Gene();
    string chr = trim(chrPos[0]);

    vector<string> range;
    split(chrPos[1], range, "-");
    if(range.size()<2)
        return Gene();

    int start = atoi(trim(range[0]).c_str());
    int end = atoi(trim(range[1]).c_str());

    return Gene(name, chr, start, end);

}