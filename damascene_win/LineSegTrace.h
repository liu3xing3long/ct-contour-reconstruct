#ifndef LineSegTrace_h__
#define LineSegTrace_h__

#include <map>
#include <vector>
using namespace std;



struct SPoint
{
	/**
	*	x, y 坐标, 对应灰度值
	*/
	int nX;
	int nY;
	float fVal;
	bool bTranversed;
	SPoint():nX(0), nY(0), fVal(0.f), bTranversed(false){}
};

class CLineSegTrace
{
public:
	CLineSegTrace(void);
	~CLineSegTrace(void);

public:
	int  initTracePoints(int* fPoints , int width, int height);
	int  initTracePoints(float* fPoints, int width, int height);
	bool traceLineSegs();

	void debugPrintOutput(char* fileSaveName = NULL, int nLengthThreholdLow = 5, int nLengthThreHoldHigh = 1000);

private:
	typedef vector<SPoint>		VPoints;
	/** 点下标-->点*/
	typedef map<int,SPoint>		MPoints;
	/** line index--->lines*/
	typedef map<int, MPoints>	MLines;

	MLines  m_mAllLines;
	MPoints m_vAllPoints;

private:
	int nWidth, nHeight;
};

#endif // LineSegTrace_h__
