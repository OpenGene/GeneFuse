#include "indexer.h"
#include "util.h"

Indexer::Indexer(string refFile){
    mRefFile = refFile;
    mReference = new FastaReader(refFile);
}

Indexer::~Indexer() {
    delete mReference;
    mReference = NULL;
}

void Indexer::makeIndex(vector<Fusion>& fusions) {
    if(mReference == NULL)
        return ;

    map<string, string> ref = mReference->mAllContigs;
    for(int ctg=0; ctg<fusions.size(); ctg++){
        Gene gene = fusions[ctg].mGene;
        string chr = gene.mChr;
        if(ref.count(chr) == 0) {
            if(ref.count("chr" + chr) >0)
                chr = "chr" + chr;
            else if(ref.count(replace(chr, "chr", "")) >0)
                chr = replace(chr, "chr", "");
            else
                continue;
        }
        string s = ref[chr].substr(gene.mStart, gene.mEnd - gene.mStart);
        str2upper(s);
        indexContig(ctg, s);
    }
}

void Indexer::indexContig(int ctg, string seq) {
    for(int i=0; i<seq.length() - KMER; ++i) {
        unsigned long kmer = makeKmer(seq, i);
        if(kmer < 0)
            continue;
        GenePos site;
        site.contig = ctg;
        site.position = i;
        // this is a dupe
        if(mKmerPos.count(kmer) >0 ){
            GenePos gp = mKmerPos['kmer'];
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
                mKmerPos['kmer'].contig = -1;
                mKmerPos['kmer'].position = mDupeList.size() -1;
            }
        }
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
        kmer = kmer << 2;
    }
    return kmer;
}