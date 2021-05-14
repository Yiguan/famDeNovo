

void str2intv(std::string & ss, std::vector<int> & intv);
std::string unphaseGT(std::string igt);

class INDI
{
public:
	INDI(std::string chrom, std::string pos, std::string ref, std::string alt, std::string format, std::string sdata, std::string sample_name); 
public:
	std::string 		CHROM;
	int	    		POS;
	std::string 		REF;
	std::vector<std::string>	ALT;
	std::string		GT;
	int 			DP;
	int			GQ;
	std::vector<int>	AD;
	std::vector<int>	ADF;
	std::vector<int>	ADR;
	std::string		SNAME;
  	std::string		SEX;
	std::string		IDENTITY;
	void echo_me();
	bool homoREF();
	bool DPmatch(int minDP, int maxDP);
	bool GQmatch(int minGQ);
	bool AB_ADF_ADRmatch(int ith_ALT, double AB_thres);
	bool impurity(int ith_ALT, int ad_min_impurity);
};

INDI::INDI(std::string chrom, std::string pos, std::string ref, std::string alt, std::string format, std::string sdata, std::string sample_name)
{
	CHROM = chrom;
	POS   = std::stoi(pos);
	REF = ref;
	mysplit(ALT, alt, ",");
	std::vector<std::string> FORMAT;
	mysplit(FORMAT, format,":");
	std::vector<std::string> SDATA;
	mysplit(SDATA, sdata, ":");
	for(int i=0;i<FORMAT.size();i++)
	{
		if(FORMAT[i]=="GT") {GT=unphaseGT(SDATA[i]);continue;}
		if(FORMAT[i]=="DP") {DP=std::stoi(SDATA[i]);continue;}
		if(FORMAT[i]=="GQ") {GQ=std::stoi(SDATA[i]);continue;}
		if(FORMAT[i]=="AD") {str2intv(SDATA[i],AD);continue;}
		if(FORMAT[i]=="ADR"){str2intv(SDATA[i],ADR);continue;}
		if(FORMAT[i]=="ADF"){str2intv(SDATA[i],ADF);continue;}
	}
	SNAME = sample_name;
}

void INDI::echo_me()
{
	std::cout << "IDENTITY:" << "\t" << IDENTITY << std::endl;
	std::cout << "SEX:" << "\t" << SEX << std::endl;
	std::cout << "SNAME:" << "\t" << SNAME << std::endl;
	std::cout << "CHROM:" << "\t" << CHROM << std::endl;
	std::cout << "POS:"  << "\t" << POS << std::endl;
	std::cout << "REF:" << "\t" << REF << std::endl;
	std::cout << "GT:" << "\t" << GT << std::endl;
	std::cout << "DP:" << "\t" << DP << std::endl;
	std::cout << "GQ:" << "\t" << GQ << std::endl;
	std::cout << "=======================" << std::endl; 
}

bool INDI::homoREF()
{
	return (GT=="0/0" || GT=="0|0")?true:false;
}

bool INDI::DPmatch(int minDP, int maxDP)
{
	if(SEX=="1" && CHROM=="X")
	{
		return (DP>=minDP/2.0 && DP<=maxDP/2.0)?true:false;
	}else
	{
		return (DP>=minDP && DP<=maxDP)?true:false;
	}	
}



bool INDI::GQmatch(int minGQ)
{
	return (GQ>=minGQ)?true:false;
}

bool INDI::AB_ADF_ADRmatch(int ith_ALT, double AB_thres)
{	// male, X chrom, don't check AB
	if(CHROM=="X" && SEX=="1")
	{
		return (ADF[ith_ALT]>0 && ADR[ith_ALT]>0)?true:false;
	}else
	{
		double ab = 1.0*AD[ith_ALT]/(AD[0]+AD[ith_ALT]);
		return (ab>=AB_thres && ab<=(1-AB_thres) && ADF[ith_ALT]>0 && ADR[ith_ALT]>0)?true:false;	
	}
}


bool INDI::impurity(int ith_ALT, int ad_min_impurity)
{
	return (AD[ith_ALT]>=ad_min_impurity)?true:false;		
}



void str2intv(std::string & ss, std::vector<int> & intv)
{
	std::vector<std::string> strv;
	mysplit(strv,ss, ",");
	for(int i=0;i<strv.size();i++)
	{
		int tt = std::stoi(strv[i]);
		intv.push_back(tt);
	}
}

// unphase GT

std::string unphaseGT(std::string igt)
{
	if(igt[1]=='/' || igt[0]=='.' || igt[2] == '.') return igt;
	else
	{
		std::string tt;
		std::vector <std::string> ttv;
		mysplit(ttv, igt, "|");
		int g1 = std::stoi(ttv[0]);
		int g2 = std::stoi(ttv[1]);
		if(g1>g2) tt = ttv[1] + "/" + ttv[0];
		else tt = ttv[0] + "/" + ttv[1];
		return tt;
	}
}
