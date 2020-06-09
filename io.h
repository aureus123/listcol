#ifndef _IO_H_
#define _IO_H_

#ifndef VISUALC
#include "mwis_sewell/mwss.h"
#else
extern "C" {
#include "mwis_sewell/mwss.h"
}
#endif

/* for linux users: do not define VISUALC */
#ifndef VISUALC
#include <unistd.h>
#include <sys/times.h>
#else
//#include <windows.h> /* definition of LP conflicts with inner definition :( */
#include <time.h>
#endif

#include <vector>
#include <string>

#define VERBOSE

void set_color(int color);
void show(std::string str);
void warning(std::string str);
void bye(std::string str);
double ECOclock();

#endif
