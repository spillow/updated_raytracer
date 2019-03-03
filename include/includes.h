#include <windows.h>
#include <GL/gl.h>
#include <GL/glut.h>
#include <GL/glui.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "mathobjects.h"
#include <vector>
#include <fstream>
#include <ios>

using namespace std;

using std::endl;
using std::ofstream;

using std::ios_base;

const short int width  = 512;
const short int height = 512;

void MakeGUI();