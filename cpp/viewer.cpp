// sizes of window
#define WIDTH 1000 
#define HEIGHT 1000
int LINK;      // 0 - simple line link
			   // 1 - spring link
			   // 2 - cylinder link

#include "object.hpp"
#include "draw3d.hpp"
#include "protein.hpp"


// ----------------------------------------------------------
// Display function
// ----------------------------------------------------------
inline void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear color and depth buffers
	glMatrixMode(GL_MODELVIEW);     // To operate on model-view matrix
									// Clear window and null buffer Z
									// glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
									// Reset transformation
	glLoadIdentity();

	for (int i = 0; i < protein.n; ++i)
	{
		glPushMatrix();
		glRotatef(rotate_x, 1.0, 0.0, 0.0);
		glRotatef(rotate_y, 0.0, 1.0, 0.0);
		glTranslated(transl_x, transl_y, 0.f);
		glScalef(zoom, zoom, zoom);
		sphere.draw(protein.x[i], protein.y[i], protein.z[i], AAcolors[protein.aa[i]][0], AAcolors[protein.aa[i]][1], AAcolors[protein.aa[i]][2]);
		glPopMatrix();
		for (const auto &j : cmap[i])
		{
			glPushMatrix();
			glRotatef(rotate_x, 1.0, 0.0, 0.0);
			glRotatef(rotate_y, 0.0, 1.0, 0.0);
			glTranslated(transl_x, transl_y, 0.f);
			glScalef(zoom, zoom, zoom);
			switch (LINK)
			{
			case 0: link<0>(.005f, protein.x[i], protein.y[i], protein.z[i], protein.x[j], protein.y[j], protein.y[j], 1.f, 1.f, 1.f);
				break;
			case 1: link<1>(.005f, protein.x[i], protein.y[i], protein.z[i], protein.x[j], protein.y[j], protein.y[j], 1.f, 1.f, 1.f);
				break;
			case 2: link<2>(.005f, protein.x[i], protein.y[i], protein.z[i], protein.x[j], protein.y[j], protein.y[j], 1.f, 1.f, 1.f);
				break;
			}
			glPopMatrix();
		}
	}
	glFlush();
	glutSwapBuffers();
}

static void show_usage()
{
	std::cerr << "Usage: Protein Viewer <option(s)>" << std::endl
		<< "Options:" << std::endl
		<< "\t-h,--help\t\tShow this help message" << std::endl
		<< "\t" << "viewer filename[.xyz] [thr] [link]" << std::endl << std::endl
		<< "\t\tfilename must have format AA-name x[coord] y[coord] z[coord]" << std::endl
		<< "\t\tthr is the value of threashold in contact map computation (default : 8)" << std::endl
		<< "\t\tlink is an integer for visualize link (default : 0) , where" << std::endl
		<< "\t\t - 0 : simple line link" << std::endl
		<< "\t\t - 1 : link as spring" << std::endl
		<< "\t\t - 2 : link as cylinder" << std::endl
		<< "\t\t\tExample : viewer test.xyz" << std::endl
		<< "\t\t\tExample : viewer test.txt 12 2" << std::endl
		<< std::endl;
}


int main(int argc, char** argv)
{
	std::string filename;
	float thr;
	int LINK;
	if (argc == 1 || (std::string)argv[1] == "--help" || (std::string)argv[1] == "-h") { show_usage(); return 0; }
	else filename = argv[1];

	switch (argc)
	{
	case 2:
	{
		thr = 8.f;
		LINK = 0;
	}break;
	case 3:
	{
		thr = std::stof(argv[2]);
		LINK = 0;
	}break;
	case 4:
	{
		thr = std::stof(argv[2]);
		LINK = std::stoi(argv[3]);
	}break;
	default: 
	{
		show_usage();
		exit(1);
	}
	}
	

	protein = Protein(filename);
	contact_map(protein, thr);
	protein.center();
	protein.normalize();
	
	draw_window(argc, argv, "Protein Viewer");
	
	return 0;
}
