#include "fusionresult.h"
#include <sstream>
#include "editdistance.h"
#include "common.h"
#include <stdlib.h>
#include "util.h"
#include <math.h>
#include "globalsettings.h"

using namespace std;

FusionResult::FusionResult() {
    mLeftIsExon = false;
    mRightIsExon = false;
    mLeftExonOrIntronID = -1;
    mRightExonOrIntronID = -1;
}

FusionResult::~FusionResult() {
    
}

void FusionResult::addMatch(Match* m) {
    mMatches.push_back(m);
}

bool FusionResult::support(Match* m) {
    for(int i=0; i<mMatches.size(); i++) {
        if(supportSame(m, mMatches[i]))
            return true;
    }
    return false;
}

bool FusionResult::supportSame(Match* m1, Match* m2) {
    const int T=3;

    if( abs(m1->mLeftGP.position - m2->mLeftGP.position) > T)
        return false;
    if( abs(m1->mRightGP.position - m2->mRightGP.position) > T)
        return false;
    if( m1->mLeftGP.contig != m2->mLeftGP.contig)
        return false;
    if( m1->mRightGP.contig != m2->mRightGP.contig)
        return false;

    return true;
}

void FusionResult::calcFusionPoint() {
    if(mMatches.size() == 0)
        return ;

    // if we can find an exact match with 0 gap, then use it
    // else we use the mean one
    long leftTotal = 0;
    long rightTotal = 0;
    for(int i=0; i<mMatches.size(); i++) {
        Match* match = mMatches[i];
        if(match->mGap == 0){
            mLeftGP = match->mLeftGP;
            mRightGP = match->mRightGP;

            return ;
        }
        leftTotal += match->mLeftGP.position;
        rightTotal += match->mRightGP.position;
    }

    mLeftGP.contig = mMatches[0]->mLeftGP.contig;
    mLeftGP.position = leftTotal/(long)mMatches.size();
    mRightGP.contig = mMatches[0]->mRightGP.contig;
    mRightGP.position = rightTotal/(long)mMatches.size();

}
    
void FusionResult::calcUnique() {
    mUnique = 1;
    // since it is sorted, so just check every match with previous one
    for(int i=1; i<mMatches.size(); i++) {
        if(mMatches[i]->mReadBreak != mMatches[i-1]->mReadBreak ||
            mMatches[i]->mRead->length() != mMatches[i-1]->mRead->length())
            mUnique ++;
    }
}

bool FusionResult::isDeletion() {
    if(mLeftGP.contig == mRightGP.contig ) {
        if(mLeftGP.position>0 && mRightGP.position>0)
            return true;
        if(mLeftGP.position<0 && mRightGP.position<0)
            return true;
    }
    return false;
}

bool FusionResult::canBeMapped() {
    if(canBeMatched(mLeftRefExt, mRightRef))
        return true;
    if(canBeMatched(mLeftRef, mRightRefExt))
        return true;
    return false;
}

bool FusionResult::canBeMatched(string& s1, string& s2) {
    int len = s1.length();
    for(int offset = -6; offset<=6; offset++) {
        int start1 = max(0, offset);
        int start2 = max(0, -offset);
        int cmplen = len - abs(offset);
        if(start1>=s1.length() || start2>= s2.length())
            return true;
        int ed = edit_distance(s1.substr(start1, cmplen), s2.substr(start2, cmplen));

        int threshold = threshold = cmplen / 10;
        if(ed <= threshold)
            return true;
    }
    return false;
}

bool FusionResult::isQualified() {
    if(mUnique<GlobalSettings::uniqueRequirement){
        return false;
    }
    if(canBeMapped())
        return false;
    if(mLeftRef.length() <= 30 || mRightRef.length()<= 30)
        return false;
    if(dis_connected_count(mLeftRef.substr(mLeftRef.length()-10, 10)) <=2)
        return false;
    if(dis_connected_count(mRightRef.substr(0, 10)) <=2)
        return false;
    return true;
}

void FusionResult::updateInfo(vector<Fusion>& fusions) {
    mLeftGene = fusions[mLeftGP.contig].mGene;
    mRightGene = fusions[mRightGP.contig].mGene;

    stringstream ss;
    if(isDeletion())
        ss  << "Deletion: ";
    else
        ss  << "Fusion: ";
    ss << mLeftGene.pos2str(mLeftGP.position) << "___";
    ss << mRightGene.pos2str(mRightGP.position) ;
    ss << "  (total: " << mMatches.size() << ", unique:" << mUnique <<")";
    mTitle = ss.str();

    mLeftPos = mLeftGene.pos2str(mLeftGP.position);
    mRightPos = mRightGene.pos2str(mRightGP.position);

    mLeftGene.getExonIntron(mLeftGP.position, mLeftIsExon, mLeftExonOrIntronID);
    mRightGene.getExonIntron(mRightGP.position, mRightIsExon, mRightExonOrIntronID);

}

void FusionResult::makeReference(string& refL, string& refR) {
    int longestLeft = 0;
    int longestRight = 0;

    for(int i=0; i<mMatches.size(); i++) {
        if(mMatches[i]->mReadBreak + 1 > longestLeft)
            longestLeft = mMatches[i]->mReadBreak + 1;
        if(mMatches[i]->mRead->length() - (mMatches[i]->mReadBreak + 1) > longestRight)
            longestRight = mMatches[i]->mRead->length() - (mMatches[i]->mReadBreak + 1);
    }

    mLeftRef = getRefSeq(refL, mLeftGP.position - longestLeft + 1, mLeftGP.position);
    mRightRef = getRefSeq(refR, mRightGP.position, mRightGP.position + longestRight - 1);

    mLeftRefExt = getRefSeq(refL, mLeftGP.position, mLeftGP.position + longestRight - 1);
    mRightRefExt = getRefSeq(refR, mRightGP.position - longestLeft + 1, mRightGP.position);
}

string FusionResult::getRefSeq(string& ref, int start, int end) {
    // check start and end are in same strand
    if( (start>=0 && end<=0) || (start<=0 && end>=0) ) {
        return "";
    }

    // check the overflow
    if(abs(start)>=ref.length() || abs(end)>=ref.length())
        return "";

    int len = abs(end - start) + 1;

    if(start <0) {
        Sequence seq(ref.substr(-end, len));
        Sequence rcseq = ~seq;
        return rcseq.mStr;
    } else {
        return ref.substr(start, len);
    }
}

void FusionResult::adjustFusionBreak() {
    for(int i=0; i<mMatches.size(); i++) {
        int smallestED = 0xFFFF;
        int shift = 0;
        for(int s=-3; s<=3; s++) {
            int leftED = 0;
            int rightED = 0;
            int ed = calcED(mMatches[i], s, leftED, rightED);
            if(ed < smallestED) {
                smallestED = ed;
                shift = s;
                mMatches[i]->mLeftDistance = leftED;
                mMatches[i]->mRightDistance = rightED;
            }
        }
        mMatches[i]->mReadBreak += shift;
        mMatches[i]->mLeftGP.position += shift;
        mMatches[i]->mRightGP.position += shift;
    }
}

int FusionResult::calcED(Match* m, int shift, int& leftED, int& rightED) {
    int readBreak = m->mReadBreak + shift;
    string seq = m->mRead->mSeq.mStr;
    int leftLen = readBreak + 1;
    int rightLen = seq.length() - leftLen;
    string leftSeq = seq.substr(0, leftLen);
    string rightSeq = seq.substr(leftLen, rightLen);

    // use the sequence near the break point to adjust
    int leftComp = min(20, min((int)leftSeq.length(), (int)mLeftRef.length()));
    int rightComp = min(20, min((int)rightSeq.length(), (int)mRightRef.length()));
    int leftPartED = edit_distance(leftSeq.substr(leftSeq.length() - leftComp, leftComp), mLeftRef.substr(mLeftRef.length() - leftComp, leftComp));
    int rightPartED = edit_distance(rightSeq.substr(0, rightComp), mRightRef.substr(0, rightComp));
    int totalED = leftPartED + rightPartED;

    // recalculate the left and right edit distance
    leftComp = min(leftLen, (int)mLeftRef.length());
    rightComp = min(rightLen, (int)mRightRef.length());
    leftED = edit_distance(leftSeq.substr(leftSeq.length() - leftComp, leftComp), mLeftRef.substr(mLeftRef.length() - leftComp, leftComp));
    rightED = edit_distance(rightSeq.substr(0, rightComp), mRightRef.substr(0, rightComp));

    return totalED;
}

void FusionResult::print(vector<Fusion>& fusions) {
    cout << endl << "#" << mTitle << endl;
    for(int i=0; i<mMatches.size(); i++) {
        cout << ">" << i+1 << ", ";
        mMatches[i]->print();
    }
}


void FusionResult::printFusionProteinHTML(ofstream& file) {
    calcLeftExonIntronNumber();
    calcRightExonIntronNumber();
    float leftSize = mLeftExonNum + mLeftIntronNum;
    float rightSize = mRightExonNum + mRightIntronNum;
    int leftPercent = round(leftSize * 100 / (leftSize + rightSize));
    int rightPercent = 100 - leftPercent;
    if(leftPercent==0)
        leftPercent = 1;
    if(rightPercent==0)
        rightPercent=1;
    file << "<table width='100%' class='protein_table'>\n";
    file << "<tr>";
    file << "<td width='" << int2str(leftPercent) << "%'>";
    file << mLeftGene.mName;
    file << "</td>";
    file << "<td width='" << int2str(rightPercent) << "%'>";
    file << mRightGene.mName;
    file << "</td>";
    file << "</tr>";
    file << "<tr>";
    file << "<td class='protein_left' width='" << int2str(leftPercent) << "%'>";
    printLeftProteinHTML(file);
    file << "</td>";
    file << "<td class='protein_right' width='" << int2str(leftPercent) << "%'>";
    printRightProteinHTML(file);
    file << "</td>";
    file << "</tr>";
    file << "</table>";
}

void FusionResult::printLeftProteinHTML(ofstream& file) {
    float totalStep = mLeftExonNum + mLeftIntronNum;
    int exon = 1;
    int intron = 1;
    int step = 1;
    float stepPercent = 100.0/totalStep;
    float halfStepPercent = stepPercent * 0.5;
    bool forward = isLeftProteinForward();
    if(!forward){
        exon = mLeftGene.mExons.size();
        intron = exon - 1;
        step = -1;
    }
    file << "<table width='100%' class='protein_table'>\n";
    file << "<tr>";

    float printExon = 0;
    float printIntron = 0;

    while(printExon < mLeftExonNum || printIntron < mLeftIntronNum) {
        if(printExon < mLeftExonNum) {
            float percent = stepPercent;
            // last one is a half exon
            if(printExon+1.0 > mLeftExonNum)
                percent = halfStepPercent;
            printExonIntronTD(file, true, forward, exon, percent, "exon_left");
            printExon += 1.0;
            exon += step;
        }
        if(printIntron < mLeftIntronNum) {
            float percent = stepPercent;
            // last one is a half intron
            if(printIntron+1.0 > mLeftIntronNum)
                percent = halfStepPercent;
            printExonIntronTD(file, false, forward, intron, percent, "intron_left");
            printIntron += 1.0;
            intron += step;
        }
    }

    file << "</tr>";
    file << "</table>";
}

void FusionResult::printRightProteinHTML(ofstream& file) {
    float totalStep = mRightExonNum + mRightIntronNum;
    int exon = mRightExonOrIntronID;
    int intron = mRightExonOrIntronID;
    int step = 1;
    float stepPercent = 100.0/totalStep;
    float halfStepPercent = stepPercent * 0.5;
    bool forward = isRightProteinForward();
    if(!forward){
        step = -1;
    }
    file << "<table width='100%' class='protein_table'>\n";
    file << "<tr>";

    float printExon = 0;
    float printIntron = 0;

    // print the first half intron
    if(!mRightIsExon) {
        printExonIntronTD(file, false, forward, intron, halfStepPercent, "intron_right");
        printIntron += 0.5;
        intron +=  step;
        if(forward)
            exon += step;
    }

    while(printExon < mRightExonNum || printIntron < mRightIntronNum) {
        if(printExon < mRightExonNum) {
            float percent = stepPercent;
            if(mRightIsExon && printExon == 0.0)
                percent = halfStepPercent;
            printExonIntronTD(file, true, forward, exon, percent, "exon_right");
            if(mRightIsExon && printExon == 0.0)
                printExon += 0.5;
            else
                printExon += 1.0;
            exon += step;
        }
        if(printIntron < mRightIntronNum) {
            float percent = stepPercent;
            printExonIntronTD(file, false, forward, intron, percent, "intron_right");
            printIntron += 1.0;
            intron += step;
        }
    }

    file << "</tr>";
    file << "</table>";
}

void FusionResult::printExonIntronTD(ofstream& file, bool isExon, bool forward, int number, float percent, string style) {
    int intPercent = (int)percent;
    if(intPercent<=0)
        intPercent = 1;
    file << "<td class='"<<style<<"' width='" << int2str(intPercent) << "%'>";
    if(isExon)
        file << "E" << int2str(number);
    else {
        if(forward) 
            file << "→";
        else
            file << "←";
    }
    file << "</td>";
}

void FusionResult::calcLeftExonIntronNumber() {
    int totalExon = mLeftGene.mExons.size();
    int totalIntron = totalExon - 1;
    if(isLeftProteinForward()) {
        if(mLeftIsExon) {
            mLeftExonNum = mLeftExonOrIntronID - 0.5;
            mLeftIntronNum = mLeftExonOrIntronID - 1;
        } else {
            mLeftExonNum = mLeftExonOrIntronID;
            mLeftIntronNum = mLeftExonOrIntronID - 0.5;
        }
    } else {
        if(mLeftIsExon) {
            mLeftExonNum = totalExon - mLeftExonOrIntronID + 0.5;
            mLeftIntronNum = totalIntron - mLeftExonOrIntronID + 1;
        } else {
            mLeftExonNum = totalExon - mLeftExonOrIntronID;
            mLeftIntronNum = totalIntron - mLeftExonOrIntronID + 0.5;
        }
    }
}

void FusionResult::calcRightExonIntronNumber() {
    int totalExon = mRightGene.mExons.size();
    int totalIntron = totalExon - 1;
    if(isRightProteinForward()) {
        if(mRightIsExon) {
            mRightExonNum = totalExon - mRightExonOrIntronID + 0.5;
            mRightIntronNum = totalIntron - mRightExonOrIntronID + 1;
        } else {
            mRightExonNum = totalExon - mRightExonOrIntronID;
            mRightIntronNum = totalIntron - mRightExonOrIntronID + 0.5;
        }
    } else {
        if(mRightIsExon) {
            mRightExonNum = mRightExonOrIntronID - 0.5;
            mRightIntronNum = mRightExonOrIntronID - 1;
        } else {
            mRightExonNum = mRightExonOrIntronID;
            mRightIntronNum = mRightExonOrIntronID - 0.5;
        }
    }
}

bool FusionResult::isLeftProteinForward() {
    if(mLeftGene.isReversed()) {
        return mLeftGP.position < 0;
    } else {
        return mLeftGP.position > 0;
    }
}

bool FusionResult::isRightProteinForward() {
    if(mRightGene.isReversed()) {
        return mRightGP.position < 0;
    } else {
        return mRightGP.position > 0;
    }
}