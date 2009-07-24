/*
This class is used for displaying simulation progress and remaining time
Author: Cristea Bogdan
Revision date: 27.08.2007

Usage example:

#include "Progress.cpp"

#define MAX ((int)5e6)

int main(void)
{
  Progress pg;
  pg.set_max(MAX);
  pg(0);
  for(int n=0;n<MAX;n++)
  {
    pg(n+1);
  }
}
 */

#include <iostream>
#include <iomanip>

class Progress //shows simulation progress and remaining time
{
public:
  void operator()(const double percent=0.0)//show progress (percent should vary from 0 to 1)
  {
  	show_etime(percent);
  };
  void set_max(const int max)//set maximum number of iterations
  {
  	max_it=(double)max;
  };
  void operator()(const int iteration=0)//show progress (iteration should vary from 0 to max_it)
  {
  	show_etime(((double)iteration)/max_it);
  };
private:
  void show_etime(const double percent);
  double max_it;//maximum number of iterations
  time_t start,stop,current;
  inline void sec2human(char* rem_time, double sec);//converts seconds to human readable format as string
};

void Progress::show_etime(const double percent)
{
  if(percent==0)
  {
    time(&start);//initialize timer
    current = start;
    std::cout << "0.0%\n";
    fflush(stdout);
    return;
  }
  time(&stop);
  if((difftime(stop, current)<1.0)&&(percent!=1.0)) return;//time diff should be at least 1 sec
  char rem_time[50];//remaining time as string
  double etime = difftime(stop, start);
  sec2human(rem_time, etime/percent-etime);
  std::cout << std::fixed << std::setw(3) << std::setprecision(1) << 100.0*percent << "%\t" << rem_time << " remaining\n";
  fflush(stdout);
  current = stop;
}

void Progress::sec2human(char* rem_time, double sec)
{
  int s = (int)(sec+0.5);//convert to an integer value (round to nearest integer)
  int d = s/86400;
  s -= d*86400;
  int h = s/3600;
  s -= h*3600;
  int m = s/60;
  s -= m*60;
  if(d>0)
    sprintf(rem_time, "%d day, %d hr, %d min, %d sec", d, h, m, s);
  else if(h>0)
    sprintf(rem_time, "%d hr, %d min, %d sec", h, m, s);
  else if(m>0)
    sprintf(rem_time, "%d min, %d sec", m, s);
  else
    sprintf(rem_time, "%d sec", s);
}


