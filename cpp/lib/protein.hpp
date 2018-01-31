#if defined(_MSC_VER)
#include <wchar.h>
#define fseeko _fseeki64
#define ftello _ftelli64
#endif

#if defined(__MINGW32__)
#define fseeko fseeko64
#define ftello ftello64
#endif

#if defined(__CYGWIN__)
#define fseeko fseek
#define ftello ftell
#endif

#if defined(__APPLE__) && !defined(__GNUC__)
#define fseeko fseek
#define ftello ftell
#endif


#include <map>
#include <array>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <string>
#include <iterator>
#include <cstring>


std::map<std::string, std::array<float, 3>> AAcolors
{
  {"ALA" , std::array<float, 3>{0.439216f, 0.858824f, 0.576471f}}, // Aquamarine
  {"ARG" , std::array<float, 3>{1.000000f, 0.110000f, 0.680000f}}, // SpicyPink
  {"ASN" , std::array<float, 3>{0.320000f, 0.490000f, 0.460000f}}, // GreenCopper
  {"ASP" , std::array<float, 3>{1.000000f, 0.498039f, 0.000000f}}, // Coral
  {"CYS" , std::array<float, 3>{0.196078f, 0.600000f, 0.800000f}}, // Skyblue
  {"GLY" , std::array<float, 3>{0.850000f, 0.850000f, 0.950000f}}, // Quartz
  {"GLU" , std::array<float, 3>{0.300000f, 0.300000f, 0.100000f}}, // NeonBlue
  {"GLN" , std::array<float, 3>{0.556863f, 0.137255f, 0.137255f}}, // Firebrick
  {"HIS" , std::array<float, 3>{0.858824f, 0.576471f, 0.439216f}}, // Tan
  {"ILE" , std::array<float, 3>{0.800000f, 0.498039f, 0.196078f}}, // Gold
  {"LYS" , std::array<float, 3>{0.576471f, 0.858824f, 0.439216f}}, // GreenYellow
  {"MET" , std::array<float, 3>{0.309804f, 0.184314f, 0.184314f}}, // IndianRed
  {"PHE" , std::array<float, 3>{0.623529f, 0.623529f, 0.372549f}}, // Khaki
  {"PRO" , std::array<float, 3>{0.560784f, 0.560784f, 0.737255f}}, // LightSteelBlue
  {"SER" , std::array<float, 3>{0.196078f, 0.800000f, 0.196078f}}, // LimeGreen
  {"THR" , std::array<float, 3>{0.196078f, 0.800000f, 0.600000f}}, // MediumAquamarine
  {"TRP" , std::array<float, 3>{0.556863f, 0.419608f, 0.137255f}}, // Sienna
  {"TYR" , std::array<float, 3>{1.000000f, 0.500000f, 0.000000f}}, // Orange
  {"VAL" , std::array<float, 3>{0.560784f, 0.737255f, 0.560784f}}, // PaleGreen
};

class Protein
{
  void error(const std::string &, const int &);
public:
  int n;
  float *x, *y, *z;
  std::string *aa;
  Protein();
  Protein(const std::string &);
  ~Protein() = default;
  Protein& operator=(const Protein &);
  Protein(const Protein &);

  void center();
  void normalize();
  void global_normalize();
};
void Protein::error(const std::string &message, const int &n)
{
  if (n)
  {
    std::cerr << "Protein error! " << message << std::endl;
    exit(n);
  }
  else exit(n);
}

Protein::Protein()
{
	this->x = nullptr;
	this->y = nullptr;
	this->z = nullptr;
	this->aa = nullptr;
	this->n = 0;
}


Protein::Protein(const std::string &filename)
{
	// read .xyz file with the first column of AA names
	std::ifstream file(filename);
	if (!file) error("File not found. Given : " + filename, 1);
	file.unsetf(std::ios_base::skipws);
	//count the newlines with an algorithm specialized for counting:
	this->n = std::count(std::istream_iterator<char>(file),
		std::istream_iterator<char>(),
		'\n') + 1;// +1 for the last line

	file.clear();
	file.seekg(0, std::ios::beg);
	file.setf(std::ios_base::skipws);
	this->x = new float[this->n]; this->y = new float[this->n]; this->z = new float[this->n];
	this->aa = new std::string[this->n];
	for (int i = 0; i < this->n; ++i)
	{
		file >> this->aa[i];
		file >> this->x[i];
		file >> this->y[i];
		file >> this->z[i];
	}
	file.close();

  // report output
	std::cout << "Protein name: " << filename << std::endl
			  << "Number of AA read: " << this->n << std::endl;
}

void Protein::center()
{
	float mean_x = std::accumulate(this->x, this->x + this->n, 0.f) / this->n;
	float mean_y = std::accumulate(this->y, this->y + this->n, 0.f) / this->n;
	float mean_z = std::accumulate(this->z, this->z + this->n, 0.f) / this->n;

  // report output
	std::cout << "Mean coordinates:" << std::endl
		      << "x\ty\tz" << std::endl
		      << mean_x << "\t" << mean_y << "\t" << mean_z << std::endl;

	// centering structure
  std::transform(this->x, this->x + this->n, this->x, [&mean_x](float &tmp) {return tmp - mean_x; });
  std::transform(this->y, this->y + this->n, this->y, [&mean_y](float &tmp) {return tmp - mean_y; });
  std::transform(this->z, this->z + this->n, this->z, [&mean_z](float &tmp) {return tmp - mean_z; });
}

void Protein::normalize()
{
	float max_x = std::numeric_limits<float>::min(), max_y = std::numeric_limits<float>::min(), max_z = std::numeric_limits<float>::min();
	for (int i = 0; i < this->n; ++i)
	{
		max_x = (std::fabs(this->x[i]) > max_x) ? std::fabs(this->x[i]) : max_x;
		max_y = (std::fabs(this->y[i]) > max_y) ? std::fabs(this->y[i]) : max_y;
		max_z = (std::fabs(this->z[i]) > max_z) ? std::fabs(this->z[i]) : max_z;
	}

  // report output
  std::cout << "Max coordinates: " << std::endl
            << "x\ty\tz" << std::endl
            << max_x << "\t" << max_y << "\t" << max_z << std::endl;

  // normalize each coordinate
  std::transform(this->x, this->x + this->n, this->x, [&max_x](float &tmp) {return tmp / max_x; });
  std::transform(this->y, this->y + this->n, this->y, [&max_y](float &tmp) {return tmp / max_y; });
  std::transform(this->z, this->z + this->n, this->z, [&max_z](float &tmp) {return tmp / max_z; });
}

void Protein::global_normalize()
{
  float max = std::numeric_limits<float>::min();
  for(int i = 0; i < this->n; ++i)
  {
    max = (std::fabs(this->x[i]) > max) ? std::fabs(this->x[i]) : max;
    max = (std::fabs(this->y[i]) > max) ? std::fabs(this->y[i]) : max;
    max = (std::fabs(this->z[i]) > max) ? std::fabs(this->z[i]) : max;
  }
  // report output
  std::cout << "Global max coordinates: " << std::endl
            << max << std::endl;
  // normalize all coordinates
  std::transform(this->x, this->x + this->n, this->x, [&max](float &tmp) {return tmp / max; });
  std::transform(this->y, this->y + this->n, this->y, [&max](float &tmp) {return tmp / max; });
  std::transform(this->z, this->z + this->n, this->z, [&max](float &tmp) {return tmp / max; });
}

Protein& Protein::operator=(const Protein &p)
{
	if (this->n != 0)
	{
		delete[] this->x;
		delete[] this->y;
		delete[] this->z;
		delete[] this->aa;
	}
	this->n = p.n;
	this->x = new float[p.n];
	this->y = new float[p.n];
	this->z = new float[p.n];
	this->aa = new std::string[p.n];
	std::memcpy(this->x, p.x, sizeof(float)*p.n);
	std::memcpy(this->y, p.y, sizeof(float)*p.n);
	std::memcpy(this->z, p.z, sizeof(float)*p.n);
	std::memcpy(this->aa, p.aa, sizeof(std::string)*p.n);
	return *this;
}
Protein::Protein(const Protein &p)
{
	if (this->n != 0)
	{
		delete[] this->x;
		delete[] this->y;
		delete[] this->z;
		delete[] this->aa;
	}
	this->n = p.n;
	this->x = new float[p.n];
	this->y = new float[p.n];
	this->z = new float[p.n];
	this->aa = new std::string[p.n];
	std::memcpy(this->x, p.x, sizeof(float)*p.n);
	std::memcpy(this->y, p.y, sizeof(float)*p.n);
	std::memcpy(this->z, p.z, sizeof(float)*p.n);
	std::memcpy(this->aa, p.aa, sizeof(std::string)*p.n);
}
Protein protein;
std::map<int, std::vector<int>> cmap;

struct
{
  template<typename T> inline T operator()(const T &x1, const T &x2, const T &y1, const T &y2, T z1 = (T)0., T z2 = (T)0.)
  {
    return std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
  }
} euclidean;

struct
{
  void operator()(const Protein &protein, float thr = 8.f)
  {
    int n_link = 0;
    for (int i = 0; i < protein.n; ++i)
      for (int j = i + 1; j < protein.n; ++j)
      if (euclidean(protein.x[i], protein.x[j], protein.y[i], protein.y[j], protein.z[i], protein.z[j]) < thr)
      {
        cmap[i].push_back(j);
        ++n_link;
      }

    // report output
    std::cout << "Number of links in Cmap with thr = " << thr << " : " << n_link << std::endl;
  }

} contact_map;


