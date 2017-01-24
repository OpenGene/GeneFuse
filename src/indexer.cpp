#include "indexer.h"
#include "util.h"

Indexer::Indexer(string refFile, vector<Fusion>& fusions) {
    mRefFile = refFile;
    mFusions = fusions;
    mReference = new FastaReader(refFile);
    mReference->readAll();
    mUniquePos = 0;
    mDupePos = 0;
}

Indexer::~Indexer() {
    delete mReference;
    mReference = NULL;
}

void Indexer::makeIndex() {
    if(mReference == NULL)
        return ;

    map<string, string> ref = mReference->mAllContigs;
    for(int ctg=0; ctg<mFusions.size(); ctg++){
        Gene gene = mFusions[ctg].mGene;
        string chr = gene.mChr;
        if(ref.count(chr) == 0) {
            if(ref.count("chr" + chr) >0)
                chr = "chr" + chr;
            else if(ref.count(replace(chr, "chr", "")) >0)
                chr = replace(chr, "chr", "");
            else{
                mFusionSeq.push_back("");
                continue;
            }
        }
        string s = ref[chr].substr(gene.mStart, gene.mEnd - gene.mStart);
        str2upper(s);
        mFusionSeq.push_back(s);
        indexContig(ctg, s);
    }
}

void Indexer::indexContig(int ctg, string seq) {
    for(int i=0; i<seq.length() - KMER; ++i) {
        long kmer = makeKmer(seq, i);
        //cout << kmer << "\t" << seq.substr(i, KMER) << endl;
        if(kmer < 0)
            continue;
        GenePos site;
        site.contig = ctg;
        site.position = i;
        // this is a dupe
        if(mKmerPos.count(kmer) >0 ){
            GenePos gp = mKmerPos[kmer];
            // already marked as a dupe
            if(gp.contig < 0) {
                mDupeList[gp.position].push_back(site);
            } else {
                // else make a new dupe entry
                vector<GenePos> gps;
                gps.push_back(gp);
                gps.push_back(site);
                mDupeList.push_back(gps);
                // and mark it as a dupe
                mKmerPos[kmer].contig = -1;
                mKmerPos[kmer].position = mDupeList.size() -1;
                mUniquePos--;
                mDupePos++;
            }
        } else {
            mKmerPos[kmer]=site;
            mUniquePos++;
        }
    }
}

vector<GenePos> Indexer::mapRead(Read* r) {
    string seq = r->mSeq.mStr;
    const int step = 1;
    for(int i=0; i<seq.length() - KMER; i += step) {
        long kmer = makeKmer(seq, i);
        if(kmer < 0)
            continue;
        if(mKmerPos.count(kmer) <0 )
            continue;
    }
}

long Indexer::makeKmer(string & seq, int pos) {
    long kmer = 0;
    for(int i=0;i<KMER;i++){
        switch(seq[pos+i]){
            case 'A':
                kmer += 0;
                break;
            case 'T':
                kmer += 1;
                break;
            case 'C':
                kmer += 2;
                break;
            case 'G':
                kmer += 3;
                break;
            default:
                return -1;
        }
        // not the tail
        if(i<KMER-1)
            kmer = kmer << 2;
    }
    return kmer;
}

long Indexer::gp2long(const GenePos& gp){
    long ret = gp.position;
    return (ret<<16) + gp.contig;
}

GenePos Indexer::long2gp(const long val){
    GenePos gp;
    gp.contig = (val & 0xFFFF);
    gp.position = val >> 16;
    return gp;
}

void Indexer::printStat() {
    cout<<"mUniquePos:"<<mUniquePos<<endl;
    cout<<"mDupePos:"<<mDupePos<<endl;
}

bool Indexer::test() {
    GenePos gp;
    gp.contig = 5;
    gp.position = 35792;
    long val = gp2long(gp);
    GenePos gp2 = long2gp(val);
    long val2 = gp2long(gp2);
    cout << val << "," << val2 <<endl;
    cout << gp2.contig << ", " << gp2.position << endl;
    return val == val2;
}