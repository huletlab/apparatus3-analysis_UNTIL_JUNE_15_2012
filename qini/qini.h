
#include <stdio.h>
#include <stdlib.h>
#include <sys/param.h>
#include <unistd.h>
#include "simpleini/SimpleIni.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <math.h>
#include <string>

#define INI_ERROR -1e11

using namespace std;

double getINI_num( string & inifile, const char *SECTION, const char *KEY);
int
setINI_num (string & inifile, const char *SECTION, const char *KEY, double val);
bool sectionExists( string & inifile, char *SECTION);
