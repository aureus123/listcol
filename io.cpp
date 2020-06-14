#include "io.h"
#include <iostream>

/* for linux users: do not define VISUALC */
#ifndef VISUALC
#include <unistd.h>
#include <sys/times.h>
#include <conio.h>
#else
#include <time.h>
#include <Windows.h>
#endif

// set_color - change color of text
void set_color(int color)
{
#ifdef VERBOSE
#ifdef VISUALC
	SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE), color);
#else
	textcolor(color);
#endif
#endif
    return;
}

// show - show message
void show(std::string str) {
#ifdef VERBOSE
	std::cout << str << std::endl;
#endif
	return;
}

// warning - show warning
void warning(std::string str) {
#ifdef VERBOSE
	set_color(13);
	std::cout << str << std::endl;
	set_color(7);
#endif
	return;
}

// bye - finish executing and show a message
void bye(std::string str)
{
#ifdef VERBOSE
	set_color(12);
	std::cout << str << std::endl;
	set_color(7);
#endif
	exit(1);
}

double ECOclock() {
#ifndef VISUALC
	// measure user-time: use it on single-threaded system in Linux (more accurate)
	struct tms buf;
	times(&buf);
	return ((double)buf.tms_utime) / (double)sysconf(_SC_CLK_TCK);
#else
	// measure standard wall-clock: use it on Windows 
	return ((double)clock()) / (double)CLOCKS_PER_SEC;
#endif
}