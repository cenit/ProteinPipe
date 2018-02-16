#include <iostream> 
#include <fstream>
#include <vector>
#include <algorithm>
#include <string>

static void show_usage(std::string exe)
{
	std::cerr << "Usage: " << exe << " <option(s)>" << std::endl
		<< "Options:" << std::endl
		<< "\t-h,--help\t\tShow this help message" << std::endl
		<< "\t" << exe << " filename separator atoms (default CA)" << std::endl << std::endl
		<< "\t\tif separator ' ' insert -s" << std::endl
		<< "\t\tif separator '\t' insert -t" << std::endl
		<< "\t\t\tExample : " << exe << " test.pdb ',' N" << std::endl
		<< "\tELSE : list of indices of atoms" << std::endl
		<< "\t\t\tExample : " << exe << " test.txt -s C N CA" << std::endl
		<< std::endl;
}

std::vector<std::string> split(const std::string &txt, const std::string &del)
{
	std::vector<std::string> token;
	std::size_t pos = txt.find_first_of(del), start = 0, end = txt.size();
	while(pos != std::string::npos)
	{
		if(pos)	token.push_back(txt.substr(start, pos));
		start += pos + 1;
		pos = txt.substr(start, end).find_first_of(del);
	}
	if(start != end) token.push_back(txt.substr(start, pos));
	return token;
}

void pdb2xyz(const std::string &filename, const std::string &sep, const std::string &atom, bool check = false)
{
	std::vector<std::string> token;
	std::string row;

	std::ifstream is(filename);
	std::ofstream os(filename + ".xyz");
	int cnt_gap = 0, val_before = 0, val_now;
	while(std::getline(is, row))
		if(row.find("ATOM") == 0 && row.find(atom) != std::string::npos)
		{
			token = split(row, sep);
			val_now = std::stoi(token[5]);
			if(val_now - 1 == val_before)
				++val_before;
			else
				++cnt_gap;
			
			os << token[3] << "\t" << token[6] << "\t" << token[7] << "\t" << token[8] << std::endl;
		}

	is.close();
	os.close();
	if(cnt_gap != 0 && check)
		std::cerr << "WARNING : Found " << cnt_gap << " gaps in atom's sequence! Pay attention!!" << std::endl;
}

bool pair_or(const std::string &row, const std::string &tkn)
{
	return (row.find(tkn) != std::string::npos);
}

bool pair_or(const std::string &row, const std::vector<std::string> &tkn)
{
	bool res = false;
	for(auto &key : tkn)
		res = res || pair_or(row, key);
	return res;
}

void pdb2xyz(const std::string &filename, const std::string &sep, const std::vector<std::string> &atoms, bool check = false)
{
	std::vector<std::string> token;
	std::string row;

	std::ifstream is(filename);
	std::ofstream os(filename + ".xyz");
	int cnt_gap = 0, val_before = 0, val_now;
	while(std::getline(is, row))
		if(row.find("ATOM") == 0 && pair_or(row, atoms))
		{
			token = split(row, sep);
			if(token[2] == atoms[0])
			{
				val_now = std::stoi(token[5]);
				if(val_now - 1 == val_before)
					++val_before;
				else
					++cnt_gap;
			}	
			os << token[3] << "\t" << token[6] << "\t" << token[7] << "\t" << token[8] << std::endl;
		}
	is.close();
	os.close();
	if(cnt_gap != 0 && check)
		std::cerr << "WARNING : Found " << cnt_gap << " gaps in atom's sequence! Pay attention!!" << std::endl;
}

bool fileExists(const std::string &file)
{
	bool ret;
	std::ifstream file_to_check(file.c_str());
	if (file_to_check.is_open())
		ret = true;
	else
		ret = false;
	file_to_check.close();
	return ret;
}

void ErrorFile(const std::string &filename)
{
	std::cerr << "File not found! Given: " << filename << std::endl;
	exit(1);
}

std::string ReadSeparator(const std::string &in)
{
	std::string sep;
	if(in == "-s") sep = " ";
	else if(in == "-t") sep = "\t";
	else sep = in;
	return sep;
}

int main(int argc, char *argv[])
{
	if(argc == 1 || (std::string)argv[1] == "--help" || (std::string)argv[1] == "-h" || argc == 2)	{ show_usage(argv[0]); return 0; }
	
	if(!fileExists(argv[1])) ErrorFile((std::string)argv[1]);
	std::string sep = ReadSeparator((std::string)argv[2]);
	std::string atom = "CA";
	if(argc == 3)
		pdb2xyz(argv[1], sep, atom, true);
	if(argc == 4)
		pdb2xyz(argv[1], sep, argv[3], true);
	if(argc > 4)
	{
		std::vector<std::string> atoms;
		for(int i = 3; i < argc; ++i)
			atoms.push_back(argv[i]);
		pdb2xyz(argv[1], sep, atoms);
	}
	return 0;
}


