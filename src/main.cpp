// Assume X as sexual chromosome,
// Filtering: male X depth = 1/2 autosome depth, genotype quality = 1/2 autosome GQ 
// Assume diploid


#include <sstream>
#include <fstream>
#include <iostream>
#include <vector>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
//include <boost/algorithm/string.hpp>
#include "strsplit.cpp"
#include "individual.cpp"
#include "parsePed.cpp"

//#define MIN_DP  10
//#define MAX_DP  150
//#define MIN_GQ  70
//#define AB_T    0.15
//#define MIN_ALT 5
//#define AD_MIN_IMPURITY  1 
// integer [1, Inf]
//#define MAX_IMPURITY_SAMPLE  0 
// integer [0,Inf]


int testAlt(std::vector<INDI> & ssvv, int ith_alt);

bool gzvcf(std::string &fname);


std::string vcfile;
std::string pedfile;
std::string outfile;
int 	MIN_DP  = 10;
int 	MAX_DP = 150;
int 	MIN_GQ = 70;
int	MIN_ALT = 5;
double  AB_T = 0.15;
int 	AD_MIN_IMPURITY = 1;
int 	MAX_IMPURITY_SAMPLE = 0;

void show_usage(std::string name)
{
	std::cerr << "Usage: " << name << " <option(s)> \n"
		  << "Options: \n"
		  << "\t--vcf                 Input VCF file or gz VCF file, assume having GT,GQ,DP,AD,ADF,ADR\n"
		  << "\t--ped                 Input PED file, not comment lines allowed\n"
		  << "\t--out		      Output file.\n"
		  << "\t-h,--help             Show help message\n"
		  << "\t--min-dp              The minimum read depth. int [0,Inf); default 10\n"
		  << "\t--max-dp              The maximum read depth. int [0,Inf); default 150\n"
		  << "\t--min-gq              The minimum genotype quality. int [0,Inf); default 70\n"
		  << "\t--min-alt	      The minimum number of ALT reads to be a denovo; int [0, Inf);default 5\n"
		  << "\t--allele-balance      The allele balance for a heterozygote\n"
		  << "\t                        allele balance less than this value or larger than 1-this value\n"
		  << "\t                        will be filtered out. double [0,0.5); default 0.15\n"
		  << "\t--ad-min-impurity     For homoRef, the number of ALT read above which will be considered\n "
		  << "\t                        as impurity sample. int [1,Inf); default 1 (most strigent)\n"
		  << "\t--max-impurity-sample The maximum number of impurity samples allowed, above which the\n" 
		  << "\t                        SNP will not be considered for mutations. int [0, Inf); default 0 (most strigent)\n"
		  << std::endl; 
}


int main(int argc, char * argv[])
{
	if(argc<7)
	{
		show_usage(argv[0]);
		return 1;
	}
	for(int h=1;h<argc;h=h+2)
	{
		std::string arg = argv[h];
		if(arg=="-h" || arg=="--help") {show_usage(argv[0]);return 0;}
		else if(arg=="--vcf") {vcfile=argv[h+1];}
		else if(arg=="--ped") {pedfile=argv[h+1];}
		else if(arg=="--out") {outfile=argv[h+1];}
		else if(arg=="--min-dp"){MIN_DP=std::stoi(argv[h+1]);}
		else if(arg=="--max-dp"){MAX_DP=std::stoi(argv[h+1]);}
		else if(arg=="--min-gq"){MIN_GQ=std::stoi(argv[h+1]);}
		else if(arg=="--min-alt"){MIN_ALT=std::stoi(argv[h+1]);}
		else if(arg=="--allele-balance"){AB_T=std::stod(argv[h+1]);}
		else if(arg=="--ad-min-impurity"){AD_MIN_IMPURITY=std::stoi(argv[h+1]);}
		else if(arg=="--max-impurity-sample"){MAX_IMPURITY_SAMPLE=std::stoi(argv[h+1]);}
		else {std::cerr << "Error:Unrecognised arguments!" << std::endl; return 1;}
	}
	if(vcfile.size()==0){std::cerr << "Error: can't find VCF file!" << std::endl; return 1;}
	if(pedfile.size()==0){std::cerr << "Error: can't find PED file!" << std::endl; return 1;}
	if(outfile.size()==0){std::cerr << "Error: can't fine output file!" << std::endl; return 1;}
	// argument log output
	std::cout << "Running with:\n"
		  << "     --min-dp               " << MIN_DP << "\n"
		  << "     --max-dp               " << MAX_DP << "\n"
		  << "     --min-gq               " << MIN_GQ << "\n"
		  << "     --min-alt	  	  " << MIN_ALT << "\n"
	   	  << "     --allele-balance       " << AB_T << "\n"
		  << "     --ad-min-impurity      " << AD_MIN_IMPURITY << "\n"
		  << "     --max-impurity-sample  " << MAX_IMPURITY_SAMPLE << std::endl;
	//*********************************
	//*********************************
	std::ofstream out(outfile);
	// write output file header
	out << "CHROM\tPOS\tREF\tALT\tCHILD\tCHILD.DP\tCHILD.AD\tCHILD.GQ\tMOTHER\tMOTHER.DP\tMOTHER.AD\tMOTHER.GQ\t"
             << "FATHER\tFATHER.DP\tFATHER.AD\tFATHER.GQ\n";
	//parse ped file
	std::vector<std::string> samples;
	std::vector<std::string> sex;
	std::string father;
	std::string mother;
	pedf2v(pedfile, samples, sex, father, mother);
	// vcf file
	std::stringstream instream;
	if(gzvcf(vcfile))
	{	
		std::ifstream file(vcfile, std::ios_base::in | std::ios_base::binary);
		boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
		inbuf.push(boost::iostreams::gzip_decompressor());
		inbuf.push(file);
		boost::iostreams::copy(inbuf, instream);
	}
	else
	{
		std::ifstream file(vcfile);
		boost::iostreams::copy(file, instream);
	}
	std::string line;
	std::vector<std::string> headerv;
	std::vector<std::string> linev;
	// meta and header lines
	while(std::getline(instream, line))
	{
		if(line[0]=='#' && line[1]=='#') continue;
		if(line[0]=='#' && line[1]=='C')
		{
			mysplit(headerv,line, "\t");
			break;
		}
	}
	// data lines
	// int ll = 0;
	while(std::getline(instream, line))
	{
		// std::cout << ll++ << std::endl;	
		std::vector<INDI> sv;
		mysplit(linev,line, "\t");
		int faid = -1;
		int moid = -1;
		for(int i=9;i<linev.size();i++)
		{
			// std::cout << i << std::endl;
			INDI si(linev[0], linev[1], linev[3], linev[4], linev[8], linev[i], headerv[i]);
			if(headerv[i]==father)
			{
				si.SEX="1";
				si.IDENTITY = "father";
				faid = i-9;
			}
			else if(headerv[i]==mother)
			{
				si.SEX="2";
				si.IDENTITY = "mother";
				moid = i-9;
			}
			else
			{
				for(int j=0; j<samples.size();j++)
				{
					if(headerv[i]==samples[j]){si.SEX=sex[j]; si.IDENTITY="child";}
				}		
			}
			sv.push_back(si); 
		}
		// loop alt alleles
		for(int k=1; k<=sv[0].ALT.size(); k++)
		{
			// std::cout << k << std::endl;
			int re = testAlt(sv, k);
			if(re == -1) continue;
			else //sv[re].echo_me();
			{
			out << sv[re].CHROM << "\t" << sv[re].POS << "\t" << sv[re].REF << "\t" << sv[re].ALT[k-1] << "\t"
			    << sv[re].SNAME << "\t" << sv[re].DP << "\t" << sv[re].AD[0] << "," << sv[re].AD[k] << "\t" << sv[re].GQ << "\t" 
			    << sv[moid].SNAME << "\t" << sv[moid].DP << "\t" << sv[moid].AD[0] << "," << sv[moid].AD[k] << "\t" << sv[moid].GQ << "\t"
			    << sv[faid].SNAME << "\t" << sv[faid].DP << "\t" << sv[faid].AD[0] << "," << sv[faid].AD[k] << "\t" << sv[faid].GQ << "\n";
			}	
		}
	}
	return 0;
}


bool gzvcf(std::string &fname)
{	
	std::vector<std::string> fv;
	mysplit(fv, fname, ".");
	if(fv[fv.size()-1]=="gz") {return true;}
	else if (fv[fv.size()-1]=="vcf") {return false;}
	else {std::cerr << "VCF file?? Check file name, please!" << std::endl; exit(1);}
	
}

int testAlt(std::vector<INDI> & ssvv, int ith_alt) //ith_alt starts from 1
{
	//ssvv[0].echo_me();
	//std::string altHet = ssvv[0].REF + "/" + ssvv[0].ALT[ith_alt-1];
	std::string altHet = "0/" + std::to_string(ith_alt);
	//std::string altHom = ssvv[0].ALT[ith_alt-1] + "/" + ssvv[0].ALT[ith_alt-1];
	std::string altHom = std::to_string(ith_alt) + "/" + std::to_string(ith_alt);
	// std::cout << altHet << "\t" << altHom << std::endl;
	// for autosome, only consider 0/alt
	// for X female, consider 0/alt
	// for X male, consider alt/alt
	std::vector<int> altHetSampleIndex;
	for(int s=0; s<ssvv.size(); s++)
	{
		// ssvv[s].echo_me();
		if(ssvv[s].IDENTITY== "father" || ssvv[s].IDENTITY=="mother")
                { 
			if(!ssvv[s].homoREF()) return -1; // parents must be homo reference
			//std::cout << ssvv[s].SNAME << "homoREF check PASS" << std::endl;
         		if(!(ssvv[s].DPmatch(MIN_DP, MAX_DP) && ssvv[s].GQmatch(MIN_GQ) && ssvv[s].AD[ith_alt]==0)) return -1; // parents DP, GQ, ADalt==0                              
                	// std::cout << ssvv[s].SNAME << ":PASS filtering!" << std::endl;
		}
		if(ssvv[s].CHROM != "X" && ssvv[s].GT==altHet) altHetSampleIndex.push_back(s);
		if(ssvv[s].CHROM == "X" && ssvv[s].SEX=="2" && ssvv[s].GT==altHet) altHetSampleIndex.push_back(s);
		if(ssvv[s].CHROM == "X" && ssvv[s].SEX=="1" && ssvv[s].GT==altHom) altHetSampleIndex.push_back(s);
	}
	// std::cout << altHetSampleIndex[0] << std::endl;
	if(altHetSampleIndex.size()!=1) return -1;
	// std::cout << "Unique mutation filtering PASS" << std::endl;
	int sid = altHetSampleIndex[0];
	if(!(ssvv[sid].DPmatch(MIN_DP, MAX_DP) && ssvv[sid].GQmatch(MIN_GQ) && ssvv[sid].AB_ADF_ADRmatch(ith_alt, AB_T) && ssvv[sid].AD[ith_alt]>=MIN_ALT)) return -1;
	// std::cout << "Child filtering PASS" << std::endl;
	// check impurity samples
	int impurity_sample = 0;
	for(int m=0;m<ssvv.size();m++)
	{
		if(m==sid) continue;
		if(ssvv[m].impurity(ith_alt,AD_MIN_IMPURITY)) impurity_sample++; 
	}
	// std::cout << impurity_sample << std::endl;
	if(impurity_sample > MAX_IMPURITY_SAMPLE) return -1;
	return sid;
}
