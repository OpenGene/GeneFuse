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
    mReversed = false;
}

Gene::Gene(const Gene& other) {
    mName = other.mName;
    mChr = other.mChr;
    mStart = other.mStart;
    mEnd = other.mEnd;
    mExons = other.mExons;
    mReversed = other.mReversed;
}

Gene::Gene() {
    mName = "invalid";
    mChr = "invalid";
    mStart = 0;
    mEnd = 0;
    mReversed = false;
}

bool Gene::valid() {
    return mName != "invalid" && mStart != 0 && mEnd != 0;
}

void Gene::addExon(Exon exon) {
    mExons.push_back(exon);
    if(mExons.size()>1) {
        if(mExons[0].start > mExons[1].start) {
            mReversed = true;
        }
    }
}

void Gene::addExon(int id, int start, int end) {
    Exon exon;
    exon.id=id;
    exon.start=start;
    exon.end=end;
    addExon(exon);
}

void Gene::print() {
    cout<<mName<<","<<mChr<<":"<<mStart<<"-"<<mEnd;
    cout<<(mReversed?" reversed":" forward")<<endl;
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

int Gene::genePos2ChrPos(int genepos) {
    int chrpos = abs(genepos) + mStart;
    if(genepos <0)
        chrpos *= -1;
    return chrpos;
}

string Gene::pos2str(int pos) {
    int pp = abs(pos) + mStart;
    stringstream ss;
    ss<<mName<<":";
    for(int i=0; i<mExons.size(); i++) {
        if(pp >= mExons[i].start && pp <= mExons[i].end) {
            ss<<"exon:"<<mExons[i].id<<"|";
            break;
        }
        if(i>0) {
            if(mReversed) {
                if(mExons[i].end < pp && pp < mExons[i-1].start){
                    ss<<"intron:"<<(mExons[i].id-1)<<"|";
                    break;
                }
            } else {
                if(mExons[i-1].end < pp && pp < mExons[i].start){
                    ss<<"intron:"<<(mExons[i].id-1)<<"|";
                    break;
                }
            }
        }
    }
    if(pos>=0)
        ss<<"+";
    else
        ss<<"-";
    ss<<mChr<<":";
    ss<<pp;

    return ss.str();
}


void Gene::getExonIntron(int pos, bool& isExon, int& number) {
    int pp = abs(pos) + mStart;
    for(int i=0; i<mExons.size(); i++) {
        if(pp >= mExons[i].start && pp <= mExons[i].end) {
            isExon = true;
            number = mExons[i].id;
            break;
        }
        if(i>0) {
            if(mReversed) {
                if(mExons[i].end < pp && pp < mExons[i-1].start){
                    isExon = false;
                    number = (mExons[i].id-1);
                    break;
                }
            } else {
                if(mExons[i-1].end < pp && pp < mExons[i].start){
                    isExon = false;
                    number = (mExons[i].id-1);
                    break;
                }
            }
        }
    }
}