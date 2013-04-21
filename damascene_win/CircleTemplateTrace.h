#ifndef CircleTemplateTrace_h__
#define CircleTemplateTrace_h__

#include <queue>
#include <map>


class CCircleTemplateTrace
{
public:
	struct CircleTemplate
	{
		float fX;
		float fY;
		float fR;

		CircleTemplate()
			:fX(0.f),fY(0.f),fR(0.f)
		{

		}

		CircleTemplate(float x, float y, float r)
			:fX(x),fY(y),fR(r)
		{

		}
	};


public:
	CCircleTemplateTrace(void);
	~CCircleTemplateTrace(void);

public:
	void initStartPoint(float circle_x, float circle_y, float radius, int width, int height);
	float findMaxradius(float* fContours, int width, int height);
	void tracePoints(float* fContour);


private:
	float m_fStartX;
	float m_fStartY;
	float m_fRadius;
	int m_iTraceWidth;
	int m_iTraceHeight;


private:
	std::queue<CircleTemplate> m_vCircleVector;
	std::queue<CircleTemplate> m_vRemainVector;

	/*
	 *	已经追踪过的圆心
	 *  int = iWidth+iHeight*width
	 */
	std::map<int, CircleTemplate> m_TracedCircle;
};

#endif // CircleTemplateTrace_h__
