#include "StdAfx.h"
#include "CircleTemplateTrace.h"


CCircleTemplateTrace::CCircleTemplateTrace(void)
	:m_fStartY(0.f), m_fStartX(0.f),m_fRadius(0.f), m_iTraceHeight(0), m_iTraceWidth(0)
{
}


CCircleTemplateTrace::~CCircleTemplateTrace(void)
{
}


template <typename T> T CLAMP(const T& value, const T& low, const T& high) 
{
	return value < low ? low : (value > high ? high : value); 
}

void CCircleTemplateTrace::initStartPoint( float circle_x, float circle_y, float radius, int width, int height )
{
	m_fRadius = radius;
	m_fStartX = circle_x;
	m_fStartY = circle_y;

	m_fStartX = CLAMP<int>(m_fStartX, 0, width);
	m_fStartY = CLAMP<int>(m_fStartY, 0, height);
	m_fRadius = CLAMP<int>(m_fRadius, 0, width>height?width:height);

	m_iTraceWidth = width;
	m_iTraceHeight = height;

	m_vCircleVector.push(CircleTemplate(circle_x, circle_y, radius));
}


//////////////////////////////////////////////////////////////////////////
float CCircleTemplateTrace::findMaxradius(float* fContours, int width, int height)
{
	float* fObvious = (float*)malloc(width* height* sizeof(float));
	memcpy(fObvious, fContours, width* height* sizeof(float));
	int nNonZeroNum = 0;
	float fSum = 0.f;
	for (int iWhole = 0; iWhole <width* height; iWhole++)
	{
		if (fContours[iWhole] > 0.f)
		{
			fSum += fContours[iWhole];
			nNonZeroNum++;
		}
	}

	float fAvg = fSum/nNonZeroNum;

	printf("Avergae Num %f .............\n", fAvg);

	for (int iWhole = 0; iWhole <width* height; iWhole++)
	{
		if (fObvious[iWhole] < fAvg)
		{
			fObvious[iWhole] = 0.f;
		}
	}


	for (int iHei = 0; iHei < height; iHei++)
	{
		for (int iWid = 0; iWid < width; iWid++)
		{
		}
	}


	return 0.f;
}


//////////////////////////////////////////////////////////////////////////
void CCircleTemplateTrace::tracePoints(float* fContour)
{
	while(!m_vCircleVector.empty())
	{
		CircleTemplate temp = m_vCircleVector.front();
		m_vCircleVector.pop();
		int x = (int)temp.fX;
		int y = (int)temp.fY;
		int r = (int)temp.fR;

		for (int i=x-r; i<x+r; i++)
		{
			for (int j=y-r; j<r+r; j++)
			{
				int xIdx = i, yIdx = j;
				xIdx = CLAMP<int>(xIdx, 0, m_iTraceWidth);
				yIdx = CLAMP<int>(yIdx, 0, m_iTraceHeight);

				if (abs(sqrt((xIdx - x)*(xIdx - x) + (yIdx - y) * ((yIdx - y))) - r)< 1.f)
				{

				}
			}
		}
	}

}
