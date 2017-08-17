#include "matcher.h"
#include "util.h"

// we use 512M memory
const int BLOOM_FILTER_LENGTH = (1<<29);
const int KMER = 16;

Matcher::Matcher(FastaReader* ref, vector<Sequence>& seqs) {
    mReference = ref;
    mUniquePos = 0;
    mDupePos = 0;
    initBloomFilter(seqs);
    makeIndex();
}

Matcher::~Matcher() {
    delete mBloomFilterArray;
    mBloomFilterArray=NULL;
}

void Matcher::initBloomFilter(vector<Sequence>& seqs) {
    mBloomFilterArray = new unsigned char[BLOOM_FILTER_LENGTH];
    memset(mBloomFilterArray, 0, BLOOM_FILTER_LENGTH);
    for(int s=0;s<seqs.size();s++) {
        Sequence seq = seqs[s];
        Sequence rSeq = ~seq;
        initBloomFilterWithSeq(seq);
        initBloomFilterWithSeq(rSeq);
    }
}

void Matcher::initBloomFilterWithSeq(Sequence& seq) {
    string str = seq.mStr;
    bool valid = false;
    for(int i=0;i<str.length()-KMER+1;++i){
        unsigned int kmer = makeKmer(str, i, valid);
        //cout << kmer << "\t" << seq.substr(i, KMER) << endl;
        if(!valid)
            continue;
        // set bloom filter
        mBloomFilterArray[kmer>>3] |= (1<<(kmer & 0x07));
    }
}

void Matcher::makeIndex() {
    if(mReference == NULL)
        return ;

    map<string, string> ref = mReference->mAllContigs;

    map<string, string>::iterator iter;
    int ctg = 0;
    for(iter = ref.begin(); iter!=ref.end(); iter++){
        string ctgName = iter->first;
        string s = iter->second;
        mContigNames.push_back(ctgName);
        str2upper(s);
        //index forward
        indexContig(ctg, s, 0);
        //index reverse complement
        //Sequence rseq = ~(Sequence(s));
        //indexContig(ctg, rseq.mStr, -s.length()+1);
        ctg++;
    }
}

void Matcher::indexContig(int ctg, string seq, int start) {
    unsigned int kmer = 0;
    bool valid = false;
    for(int i=0; i<seq.length() - KMER; ++i) {
        if(valid){
            char base = seq[i+KMER-1];
            int num = base2num(base);
            if(num<0){
                valid = false;
                continue;
            } else {
                kmer = (kmer<<2) | num;
            }
        } else {
            kmer = makeKmer(seq, i, valid);
            if(!valid)
                continue;
        }

        // check bloom filter
        if( (mBloomFilterArray[kmer>>3] & (1<<(kmer & 0x07))) == 0)
            continue;

        GenePos site;
        site.contig = ctg;
        site.position = i+start;

        if(mKmerPositions.count(kmer) == 0) {
            mKmerPositions[kmer] = vector<GenePos>();
        } 

        mKmerPositions[kmer].push_back(site);

    }
}

MatchResult* Matcher::match(Sequence& sequence) {
    Sequence rcseq = ~sequence;

    MatchResult* mc = mapToIndex(sequence);
    if(mc!=NULL)
        mc->reversed=false;
    MatchResult* rcmc = mapToIndex(rcseq);
    if(rcmc!=NULL)
        rcmc->reversed=true;

    if(mc==NULL)
        return rcmc;
    else if(rcmc==NULL)
        return mc;
    else {
        if(mc->mismatches.size() <= rcmc->mismatches.size()){
            delete rcmc;
            return mc;
        } else {
            delete mc;
            return rcmc;
        }
    }
}

MatchResult* Matcher::mapToIndex(Sequence& sequence) {
    map<long, int> kmerStat;
    kmerStat[0]=0;
    string seq = sequence.mStr;
    const int step = 1;
    int seqlen = seq.length();
    // first pass, we only want to find if this seq can be partially aligned to the target
    bool valid = false;
    for(int i=0; i< seqlen - KMER + 1; i += step) {
        unsigned int kmer = makeKmer(seq, i, valid);
        if(!valid)
            continue;
        // no match
        if(mKmerPositions.count(kmer) <=0 ){
            kmerStat[0]++;
            continue;
        }

        for(int g=0; g<mKmerPositions[kmer].size();g++) {
            long gplong = gp2long(shift(mKmerPositions[kmer][g], i));
            if(kmerStat.count(gplong)==0)
                kmerStat[gplong] = 1;
            else
                kmerStat[gplong] += 1;
        }
    }
    // get top N
    const int TOP = 5;
    long topgp[TOP] ={0};
    int topcount[TOP] = {0};
    map<long, int>::iterator iter;
    for(iter = kmerStat.begin(); iter!=kmerStat.end(); iter++){
        long gp = iter->first;
        int count = iter->second;
        // no need to update the top N
        if(gp == 0 || count <= topcount[TOP-1])
            continue;
        // update the last one first
        topgp[TOP-1]=gp;
        topcount[TOP-1]=count;
        // compare with the rest ones
        for(int t=TOP-2;t>=0;t--){
            if(count > topcount[t]) {
                topcount[t+1] = topcount[t];
                topgp[t+1] = topgp[t];
                topcount[t] = count;
                topgp[t] = gp;
            }
        }
    }

    for(int t=0;t<TOP;t++){
        unsigned char* mask = new unsigned char[seqlen];
        memset(mask, 0, sizeof(unsigned char)*seqlen);

        // make the mask
        bool valid = false;
        for(int i=0; i< seqlen - KMER + 1; i += step) {
            long kmer = makeKmer(seq, i, valid);
            if(!valid || mKmerPositions.count(kmer) <=0)
                continue;
            for(int g=0; g<mKmerPositions[kmer].size();g++) {
                long gplong = gp2long(shift(mKmerPositions[kmer][g], i));
                if(abs(gplong - topgp[t]) <= 2){
                    for(int m=i;m<seqlen && m<i+KMER;m++)
                        mask[m]=1;
                }
            }
        }
        vector<int> mismatches = vector<int>();
        for(int i=0;i<seqlen;i++) {
            if(mask[i]==0)
                mismatches.push_back(i);
        }
        if(mismatches.size()<10) {
            MatchResult* mr = new MatchResult();
            mr->sequence = Sequence(sequence);
            mr->startGP = long2gp(topgp[t]);
            mr->mismatches = mismatches;
            delete mask;
            return mr;
        }
        delete mask;
        mask = NULL;
    }

    return NULL;
}

void Matcher::makeMask(unsigned char* mask, unsigned char flag, int seqlen, int start, int kmerSize) {
    for(int i=start;i<seqlen && i<start+kmerSize;i++)
        mask[i]=max(mask[i], flag);
}

int Matcher::base2num(char c) {
    switch(c){
        case 'A':
            return 0;
            break;
        case 'T':
            return 1;
            break;
        case 'C':
            return 2;
            break;
        case 'G':
            return 3;
            break;
        default:
            return -1;
    }
}

unsigned int Matcher::makeKmer(string & seq, int pos, bool& valid) {
    unsigned int kmer = 0;
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
                valid = false;
                return 0;
        }
        // not the tail
        if(i<KMER-1)
            kmer = kmer << 2;
    }
    valid = true;
    return kmer;
}

long Matcher::gp2long(const GenePos& gp){
    long ret = gp.contig;
    return (ret<<32) + gp.position;
}

GenePos Matcher::long2gp(const long val){
    GenePos gp;
    gp.position = (val & 0xFFFFFFFF);
    gp.contig = val >> 32;
    return gp;
}

GenePos Matcher::shift(const GenePos& gp, int i){
    GenePos gpNew;
    gpNew.contig = gp.contig;
    gpNew.position = gp.position - i;
    return gpNew;
}

void Matcher::printStat() {
    cout<<"mUniquePos:"<<mUniquePos<<endl;
    cout<<"mDupePos:"<<mDupePos<<endl;
}

bool Matcher::test() {
    return true;
}