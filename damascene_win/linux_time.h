#ifndef linux_time_h__
#define linux_time_h__


#include <time.h>
#include <windows.h>



int	gettimeofday(struct timeval *tp, void *tzp);
#endif // linux_time_h__
