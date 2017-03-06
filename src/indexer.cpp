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

vector<SeqMatch> Indexer::mapRead(Read* r) {
    map<long, int> kmerStat;
    kmerStat[0]=0;
    string seq = r->mSeq.mStr;
    const int step = 1;
    int seqlen = seq.length();
    // first pass, we only want to find if this seq can be partially aligned to the target
    for(int i=0; i< seqlen - KMER + 1; i += step) {
        long kmer = makeKmer(seq, i);
        if(kmer < 0)
            continue;
        // no match
        if(mKmerPos.count(kmer) <=0 ){
            kmerStat[0]++;
            continue;
        }
        GenePos gp = mKmerPos[kmer];
        // is a dupe
        if(gp.contig < 0) {
            for(int g=0; g<mDupeList[gp.position].size();g++) {
                long gplong = gp2long(shift(mDupeList[gp.position][g], i));
                if(kmerStat.count(gplong)==0)
                    kmerStat[gplong] = 1;
                else
                    kmerStat[gplong] += 1;
            }
        } else {
            long gplong = gp2long(shift(gp, i));
            if(kmerStat.count(gplong)==0)
                kmerStat[gplong] = 1;
            else
                kmerStat[gplong] += 1;
        }
    }
    // get 1st and 2nd hit
    long gp1 = 0;
    int count1 = 0;
    long gp2 = 0;
    int count2 = 0;
    map<long, int>::iterator iter;
    for(iter = kmerStat.begin(); iter!=kmerStat.end(); iter++){
        if(iter->first != 0 && iter->second > count1){
            gp2 = gp1;
            count2 = count1;
            gp1 = iter->first;
            count1 = iter->second;
        }  else if(iter->first != 0 && iter->second > count2 ){
            gp2 = iter->first;
            count2 = iter->second;
        }  
    }
    if(count1 < 20){
        // return an null list
        return vector<SeqMatch>();
    }

    unsigned char* mask = new unsigned char[seqlen];
    memset(mask, MATCH_UNKNOWN, sizeof(unsigned char)*seqlen);

    // second pass, make the mask
    for(int i=0; i< seqlen - KMER + 1; i += step) {
        long kmer = makeKmer(seq, i);
        if(kmer < 0 || mKmerPos.count(kmer) <=0)
            continue;
        GenePos gp = mKmerPos[kmer];
        // is a dupe
        if(gp.contig < 0) {
            for(int g=0; g<mDupeList[gp.position].size();g++) {
                long gplong = gp2long(shift(mDupeList[gp.position][g], i));
                if(abs(gplong - gp1) <= 1)
                    makeMask(mask, MATCH_TOP, seqlen, i, KMER);
                else if(abs(gplong - gp2) <= 1)
                    makeMask(mask, MATCH_SECOND, seqlen, i, KMER);
                else if(gplong == 0)
                    makeMask(mask, MATCH_NONE, seqlen, i, KMER);
            }
        } else {
            long gplong = gp2long(shift(gp, i));
            if(abs(gplong - gp1) <= 1)
                    makeMask(mask, MATCH_TOP, seqlen, i, KMER);
            else if(abs(gplong - gp2) <= 1)
                makeMask(mask, MATCH_SECOND, seqlen, i, KMER);
            else if(gplong == 0)
                makeMask(mask, MATCH_NONE, seqlen, i, KMER);
        }
    }


    vector<SeqMatch> result = segmentMask(mask, seqlen, long2gp(gp1), long2gp(gp2));
    delete mask;

    return result;
}

void Indexer::makeMask(unsigned char* mask, unsigned char flag, int seqlen, int start, int kmerSize) {
    for(int i=start;i<seqlen && i<start+kmerSize;i++)
        mask[i]=max(mask[i], flag);
}

vector<SeqMatch> Indexer::segmentMask(unsigned char* mask, int seqlen, GenePos gp1, GenePos gp2) {
    vector<SeqMatch> result;

    const int ALLOWED_GAP = 10;
    const int THRESHOLD_LEN = 30;

    int targets[2] = {MATCH_TOP, MATCH_SECOND};
    GenePos gps[2] = {gp1, gp2};

    cout<<"gp1,"<<gp1.contig<<":"<<gp1.position<<endl;
    cout<<"gp2,"<<gp2.contig<<":"<<gp2.position<<endl;
    for(int i=0;i<seqlen;i++)
        cout<<(int)mask[i];
    cout << endl;

    for(int t=0; t<2; t++){
        int maxStart = -1;
        int maxEnd = -1;

        // get gp1
        int target = targets[t];
        int start = 0;
        int end = 0;
        while(true){
            // get next start
            while(mask[start] != target && start != seqlen-1)
                start++;

            // reach the tail
            if(start >= seqlen-1)
                break;

            if(mask[start] == target){
                end = start+1;
                // get the end
                int g=0;
                while(g<ALLOWED_GAP && end+g<seqlen){
                    if(mask[end+g] > target)
                        break;
                    if(end+g < seqlen && mask[end+g] == target ){
                        end += g+1;
                        g = 0;
                        continue;
                    }
                    g++;
                }
                // left shift to remove the mismatched end
                end--;
                if(end - start > maxEnd - maxStart){
                    maxEnd = end;
                    maxStart = start;
                }
                start++;
            } else {
                // not found
                break;
            }
        }
        if(maxEnd - maxStart >  THRESHOLD_LEN){
            SeqMatch match;
            match.seqStart = maxStart;
            match.seqEnd = maxEnd;
            match.startGP = gps[t];
            result.push_back(match);
        }
    }

    return result;
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
    long ret = gp.contig;
    return (ret<<32) + gp.position;
}

GenePos Indexer::long2gp(const long val){
    GenePos gp;
    gp.position = (val & 0xFFFFFFFF);
    gp.contig = val >> 32;
    return gp;
}

GenePos Indexer::shift(const GenePos& gp, int i){
    GenePos gpNew;
    gpNew.contig = gp.contig;
    gpNew.position = gp.position - i;
    return gpNew;
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