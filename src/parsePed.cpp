// read in PED file

bool samev(std::vector<std::string> sv)
{
	for(int i=0; i<sv.size(); i++)
	{
		if(sv[i]!=sv[0]) return false;
	}
	return true;
}




void pedf2v(std::string pedfile, std::vector<std::string> & sampleid, std::vector<std::string> & samplesex, 
		std::string & father, std::string & mother)
{
	std::vector<std::string> fatherv;
	std::vector<std::string> motherv;
	std::ifstream pedf(pedfile);
	std::string line;
	if(pedf.is_open())
	{	
		std::vector<std::string> linev;
		while(getline(pedf, line))
		{
			mysplit(linev, line, "\t");
			sampleid.push_back(linev[1]);
			fatherv.push_back(linev[2]);
			motherv.push_back(linev[3]);
			samplesex.push_back(linev[4]);
		}
		pedf.close();
	}else
	{
		std::cout << "Error: Can't open PED file!" << std::endl;
		exit(1);
	}
	// check unique parents
	if(!samev(fatherv)) {std::cout << "Error: Multiple father in PED?!" << std::endl; exit(1);}
	father = fatherv[0];
	if(!samev(motherv)) {std::cout << "Error: Multiple mother in PED?!" << std::endl; exit(1);}
	mother = motherv[0];
}
