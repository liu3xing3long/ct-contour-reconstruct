#ifndef common_func_h__
#define common_func_h__


#include <time.h>
#include <windows.h>


namespace CommonFunc
{
	template <typename T> T CLAMP(const T& value, const T& low, const T& high) 
	{
		return value < low ? low : (value > high ? high : value); 
	}


	int gettimeofday(struct timeval *tp, void *tzp);


};



#endif // common_func_h__
