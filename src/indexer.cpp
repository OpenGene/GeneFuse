#include "indexer.h"
#include "util.h"
#include <memory.h>
#include "globalsettings.h"

const int KMER = 16;
// 512M bloom filter
const unsigned int BLOOM_FILTER_SIZE = 1<<29;
const unsigned int BLOOM_FILTER_BITS = BLOOM_FILTER_SIZE - 1;

Indexer::Indexer(string refFile, vector<Fusion>& fusions) {
    mRefFile = refFile;
    mFusions = fusions;
    mReference = new FastaReader(refFile, false);
    mReference->readAll();
    mUniquePos = 0;
    mDupePos = 0;
    mBloomFilter = new unsigned char[BLOOM_FILTER_SIZE];
    memset(mBloomFilter, 0, sizeof(unsigned char) * BLOOM_FILTER_SIZE);
}

Indexer::~Indexer() {
    delete mReference;
    mReference = NULL;
    delete mBloomFilter;
    mBloomFilter = NULL;
}

FastaReader* Indexer::getRef() {
    return mReference;
}

void Indexer::makeIndex() {
    if(mReference == NULL)
        return ;

    map<string, string>& ref = mReference->mAllContigs;
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
        //index forward
        indexContig(ctg, s, 0);
        //index reverse complement
        Sequence seq = ~(Sequence(s));
        indexContig(ctg, seq.mStr, -s.length()+1);
    }
    fillBloomFilter();
    loginfo("mapper indexing done");
}

void Indexer::fillBloomFilter() {
    unordered_map<long, GenePos>::iterator iter;
    for(iter = mKmerPos.begin(); iter!=mKmerPos.end(); iter++) {
        long kmer = iter->first;
        long pos = kmer>>3;
        long bit = kmer & 0x07;
        mBloomFilter[pos] |= (0x1<<bit);
    }
}

void Indexer::indexContig(int ctg, string seq, int start) {
    long kmer = -1;
    for(int i=0; i<seq.length() - KMER; ++i) {
        kmer = makeKmer(seq, i, kmer);
        //cout << kmer << "\t" << seq.substr(i, KMER) << endl;
        if(kmer < 0)
            continue;
        GenePos site;
        site.contig = ctg;
        site.position = i+start;
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
    const int step = 2;
    int seqlen = seq.length();
    // first pass, we only want to find if this seq can be partially aligned to the target
    long kmer = -1;
    for(int i=0; i< seqlen - KMER + 1; i += step) {
        kmer = makeKmer(seq, i, kmer, step);
        if(kmer < 0)
            continue;
        long pos = kmer>>3;
        long bit = kmer & 0x07;
        if( (mBloomFilter[pos] & (0x1<<bit)) == 0) {
            kmerStat[0]++;
            continue;
        }
        GenePos gp = mKmerPos[kmer];
        // is a dupe
        if(gp.contig < 0) {
            // too much keys in this dupe, then skip it
            if(mDupeList[gp.position].size() > GlobalSettings::skipKeyDupThreshold)
                continue;
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
    //TODO: handle small difference caused by INDEL
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
    if(count1 * step < GlobalSettings::majorGeneKeyRequirement || count2 * step < GlobalSettings::minorGeneKeyRequirement){
        // return an null list
        return vector<SeqMatch>();
    }

    unsigned char* mask = new unsigned char[seqlen];
    memset(mask, MATCH_UNKNOWN, sizeof(unsigned char)*seqlen);

    // second pass, make the mask
    kmer = -1;
    for(int i=0; i< seqlen - KMER + 1; i += 1) {
        kmer = makeKmer(seq, i, kmer);
        if(kmer < 0)
            continue;
        long pos = kmer>>3;
        long bit = kmer & 0x07;
        if( (mBloomFilter[pos] & (0x1<<bit)) == 0)
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

    int mismatches = 0;
    for(int i=0; i<seqlen; i++){
        if(mask[i] == MATCH_NONE || mask[i] == MATCH_UNKNOWN )
            mismatches++;
    }

    if(mismatches>GlobalSettings::mismatchThreshold){
        // too many mismatch indicates not a real fusion
        return vector<SeqMatch>();
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
    const int THRESHOLD_LEN = 20;

    int targets[2] = {MATCH_TOP, MATCH_SECOND};
    GenePos gps[2] = {gp1, gp2};

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

    /*if(result.size()>=2 && inRequiredDirection(result)){
        cout<<"gp1,"<<gp1.contig<<":"<<gp1.position<<", ";
        cout<<"gp2,"<<gp2.contig<<":"<<gp2.position<<endl;
        for(int i=0;i<seqlen;i++)
            cout<<(int)mask[i];
        cout << endl;
    }*/

    return result;
}

// this function is to gurantee that all the supporting reads will have same direction
bool Indexer::inRequiredDirection(vector<SeqMatch>& mapping) {
    if(mapping.size()<2)
        return false;

    SeqMatch left = mapping[0];
    SeqMatch right = mapping[1];
    if(left.seqStart > right.seqStart) {
        left = mapping[1];
        right = mapping[0];
    }

    // both are positive, good to go
    if(left.startGP.position >0 && right.startGP.position >0)
        return true;

    // if both are negative, we should use their reverse complement, which will be both positive
    if(left.startGP.position <0 && right.startGP.position <0)
        return false;
    
    // if one is positive, the other is negative, their reverse complement will be the same
    if( mFusions[left.startGP.contig].isReversed() && !mFusions[right.startGP.contig].isReversed()){
        // if left is reversed gene and right is forward gene
        // we should use their reverse complement, which keep the left forward gene
        return false;
    } 
    else if( !mFusions[left.startGP.contig].isReversed() && mFusions[right.startGP.contig].isReversed()){
        // if left is forward gene and right is reversed gene, good to go
        return true;
    }
    else {
        // otherwise, we should keep the left has smaller contig
        if(left.startGP.contig < right.startGP.contig)
            return true;
        // or smaller positive if contig is the same
        if(left.startGP.contig == right.startGP.contig && abs(left.startGP.position) < abs(left.startGP.position))
            return true;
        else
            return false;
    }
    return false;
}

long Indexer::makeKmer(string & seq, int pos, long lastKmer, int step) {
    // else calculate it completely
    long kmer = 0;
    int start = 0;
    // re-use several bits
    if(lastKmer >= 0) {
        kmer = lastKmer;
        start = KMER - step;
        if(step == 1)
            kmer = ((kmer & 0x3FFFFFFF) << 2);
        else if (step == 2)
            kmer = ((kmer & 0x0FFFFFFF) << 2);
        else if (step == 3)
            kmer = ((kmer & 0x03FFFFFF) << 2);
        else if (step == 4)
            kmer = ((kmer & 0x00FFFFFF) << 2);
    }
    for(int i=start;i<KMER;i++){
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
    int temp[2]={gp.position, 0};
    return (ret<<32) | *((long*)temp);
}

GenePos Indexer::long2gp(const long val){
    GenePos gp;
    gp.position = (val & 0x00000000FFFFFFFF);
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
    short contigs[10]={0, 1, 3, 220, -1, 0, 23, 4440, 110, 10};
    int positions[10]={0, 111, 222, -333, 444, 555555, 6, -7777777, 8888, -9999};
    for(int i=0;i<10;i++){
        GenePos gp;
        gp.contig = contigs[i];
        gp.position = positions[i];
        long val = gp2long(gp);
        GenePos gp2 = long2gp(val);
        long val2 = gp2long(gp2);
        cout << gp2.contig << ", " << gp2.position << endl;
        if ((gp.contig == gp2.contig && gp.position == gp2.position && val==val2) == false)
            return false;
    }
    return true;
}