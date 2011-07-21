
#include <stdio.h>
#include <stdlib.h>
#include <sys/param.h>
#include <unistd.h>
#include "/home/pedro/cpp/simpleini/SimpleIni.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <math.h>
#include <string>

using namespace std;

double getINI_num( string & inifile, char *SECTION, char *KEY);
int
setINI_num (string & inifile, char *SECTION, char *KEY, double val);
bool sectionExists( string & inifile, char *SECTION);
