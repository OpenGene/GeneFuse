#include "fusion.h"
#include "editdistance.h"
#include <iostream>
#include <fstream>
#include "util.h"
#include <string.h>
#include "builtinfusion.h"
#include <sstream>
#include "globalsettings.h"

Fusion::Fusion(Gene gene){
	mGene = gene;
}

vector<Fusion> Fusion::parseCsv(string filename) {
    ifstream file;
    file.open(filename.c_str(), ifstream::in);
    const int maxLine = 4096;
    char line[maxLine];
    vector<Fusion> fusions;
    Gene workingGene;
    while(file.getline(line, maxLine)){
        // trim \n, \r or \r\n in the tail
        int readed = strlen(line);
        if(readed >=2 ){
            if(line[readed-1] == '\n' || line[readed-1] == '\r'){
                line[readed-1] = '\0';
                if(line[readed-2] == '\r')
                    line[readed-2] = '\0';
            }
        }
        string linestr(line);
        linestr = trim(linestr);
        vector<string> splitted;
        split(linestr, splitted, ",");
        // wrong line
        if(splitted.size()<2)
            continue;
        // comment line
        if(starts_with(splitted[0], "#"))
            continue;
        // gene line
        if(starts_with(splitted[0], ">")){
            if(workingGene.valid()){
                Fusion fusion(workingGene);
                fusions.push_back(fusion);
            }
            workingGene = Gene::parse(linestr);
            continue;
        }
        // position line require id, start, position
        if(splitted.size()<3)
            continue;

        int id = atoi(trim(splitted[0]).c_str());
        int start = atoi(trim(splitted[1]).c_str());
        int end = atoi(trim(splitted[2]).c_str());
        workingGene.addExon(id, start, end);
    }
    // last one
    if(workingGene.valid()){
        Fusion fusion(workingGene);
        fusions.push_back(fusion);
    }
    return fusions;
}

void Fusion::print(){
    mGene.print();
}

void Fusion::printHtml(ofstream& file){
}

string Fusion::pos2str(int pos) {
    return mGene.pos2str(pos);
}

bool Fusion::test() {
    vector<Fusion> fusions = Fusion::parseCsv("testdata/fusions.csv");
    for(int i=0;i<fusions.size();++i) {
        Fusion f = fusions[i];
        if(f.mGene.mName == "ALK") {
            // exon 20
            if(f.pos2str(-30582)!="ALK:exon:20|-chr2:29446222")
                return false;
            // intron 19
            if(f.pos2str(31060)!="ALK:intron:19|+chr2:29446700")
                return false;
        }
        if(f.mGene.mName == "EML4") {
            // exon 6
            if(f.pos2str(95365)!="EML4:exon:6|+chr2:42491855")
                return false;
            // intron 5
            if(f.pos2str(95346)!="EML4:intron:5|+chr2:42491836")
                return false;
        }
    }
    return fusions.size() == 4;
}
