#include "StdAfx.h"
#include "LineSegTrace.h"
#include "common_func.h"
#include "common\LCGRandom.h"
#include "cutil.h"
#include "pgm.cuh"

#define  TRANS_3X3
//#define  TRANS_5X5

CLineSegTrace::CLineSegTrace(void)
	:nWidth(0), nHeight(0)
{
}


CLineSegTrace::~CLineSegTrace(void)
{
}

//////////////////////////////////////////////////////////////////////////
int CLineSegTrace::initTracePoints(int* fPoints , int width, int height)
{
	for (int iHei = 0; iHei < height; iHei++)
	{
		for (int iWid = 0; iWid < width; iWid++)
		{
			int iDx = iHei*width+ iWid;
			if (fPoints[iDx] > 0)
			{/**不是纯黑的点*/
				SPoint pt;
				pt.nX = iWid;
				pt.nY = iHei;
				pt.fVal = (float)fPoints[iDx];
				m_vAllPoints.insert(std::make_pair(iDx, pt));
			}
		}
	}

	nWidth = width;
	nHeight= height;

	return m_vAllPoints.size();
}


//////////////////////////////////////////////////////////////////////////
int CLineSegTrace::initTracePoints( float* fPoints , int width, int height)
{
	for (int iHei = 0; iHei < height; iHei++)
	{
		for (int iWid = 0; iWid < width; iWid++)
		{
			int iDx = iHei*width+ iWid;
			if (fPoints[iDx] > 0.f)
			{/**不是纯黑的点*/
				SPoint pt;
				pt.nX = iWid;
				pt.nY = iHei;
				pt.fVal = fPoints[iDx];
				m_vAllPoints.insert(std::make_pair(iDx, pt));
			}
		}
	}

	nWidth = width;
	nHeight= height;

	return m_vAllPoints.size();
}

//////////////////////////////////////////////////////////////////////////
bool CLineSegTrace::traceLineSegs()
{
	if(m_vAllPoints.empty())
		return false;

	int iLineIdx = 0;

#ifdef TRANS_3X3
	int xTrans[] = {-1, 0, 1};
	int yTrans[] = {-1, 0, 1};

	const int N_SEARCH_DIM = 3;
	const int N_OCT_SEARCH = N_SEARCH_DIM*N_SEARCH_DIM;
#endif

#ifdef TRANS_5X5
	int xTrans[] = {-2,-1,0,1,2};
	int yTrans[] = {-2,-1,0,1,2};
	const int N_SEARCH_DIM = 5;
	const int N_OCT_SEARCH = N_SEARCH_DIM*N_SEARCH_DIM;

#endif

	while (!m_vAllPoints.empty())
	{
		//printf("points left %d................... \n", m_vAllPoints.size());

		/**获取第0个必然存在的元素*/
		MPoints::iterator itPt = m_vAllPoints.begin();
		SPoint pt = itPt->second;

		VPoints vTmpPts;
		MPoints mCurrentLine;
		mCurrentLine.clear();
		vTmpPts.clear();

		vTmpPts.push_back(pt);
		/**递归以此为原点的所有点*/
		while(!vTmpPts.empty())
		{
			SPoint tmp_pt = vTmpPts.back();
			vTmpPts.pop_back();

 			mCurrentLine.insert(make_pair(tmp_pt.nX + tmp_pt.nY*nWidth, tmp_pt));
			/*
			* 八方向查找, 左上、中上、右上； 左中、右中；...
			* 中间点虽然参与了计算，但是是个多余计算，不会被重新添加到点集中
			*/
			//SPoint tl, tm, tr, lm, rm, bl, bm, br;
			/**
			//**top
			*octpt[0].nX = tmp_pt.nX - 1;
			*octpt[0].nY = tmp_pt.nY - 1;
			*octpt[1].nX = tmp_pt.nX;
			*octpt[1].nY = tmp_pt.nY - 1;
			*octpt[2].nX = tmp_pt.nX + 1;
			*octpt[2].nY = tmp_pt.nY - 1;
			//**middle
			*octpt[3].nX = tmp_pt.nX - 1;
			*octpt[3].nY = tmp_pt.nY;
			*octpt[4].nX = tmp_pt.nX;
			*octpt[4].nY = tmp_pt.nY;
			*octpt[5].nX = tmp_pt.nX + 1;
			*octpt[5].nY = tmp_pt.nY;
			//**bottom
			*octpt[6].nX = tmp_pt.nX - 1;
			*octpt[6].nY = tmp_pt.nY + 1;
			*octpt[7].nX = tmp_pt.nX;
			*octpt[7].nY = tmp_pt.nY + 1;
			*octpt[8].nX = tmp_pt.nX + 1;
			*octpt[8].nY = tmp_pt.nY + 1;
			*/
			SPoint octpt[N_OCT_SEARCH];
			int	   idtPt[N_OCT_SEARCH];
			for (int yIdx = 0; yIdx < N_SEARCH_DIM; yIdx++)
			{
				for (int xIdx = 0; xIdx < N_SEARCH_DIM; xIdx++)
				{
					int idx = xIdx + yIdx*N_SEARCH_DIM;
					//ASSERT(idx < N_OCT_SEARCH);
					octpt[idx].nX = tmp_pt.nX + xTrans[xIdx];
					octpt[idx].nY = tmp_pt.nY + yTrans[yIdx];
				}
			}

			for (int i = 0; i< N_OCT_SEARCH; i++)
			{
				octpt[i].nX = CommonFunc::CLAMP<int>(octpt[i].nX,  0, nWidth);
				octpt[i].nY = CommonFunc::CLAMP<int>(octpt[i].nY,  0, nHeight);

				idtPt[i] = octpt[i].nX + octpt[i].nY*nWidth;
				MPoints::iterator itFind = mCurrentLine.find(idtPt[i]);
				MPoints::iterator itFindWhole = m_vAllPoints.find(idtPt[i]);

				/** 不在当前点集，但是是一个valid点*/
				if (itFind == mCurrentLine.end() && itFindWhole != m_vAllPoints.end()){
					vTmpPts.push_back(octpt[i]);
				}
			}
		}

		/**
		* 将当前直线的点从点集中取出
		*/
		MPoints::iterator itBegin=mCurrentLine.begin(), itEnd=mCurrentLine.end();
		for (;itBegin != itEnd; itBegin++)
		{
			MPoints::iterator allPtFind = m_vAllPoints.find(itBegin->first);
			if (allPtFind != m_vAllPoints.end()){
				m_vAllPoints.erase(allPtFind);				
			}
		}

		m_mAllLines.insert(std::make_pair(iLineIdx, mCurrentLine));
		iLineIdx++;
	}

	return true;
}

//////////////////////////////////////////////////////////////////////////
void CLineSegTrace::debugPrintOutput( char* fileSaveName /*= NULL*/, int nLengthThreholdLow /*= 5*/,  int nLengthThreHoldHigh /*= 100*/ )
{
	MLines::iterator itLineBegin = m_mAllLines.begin(), itLineEnd = m_mAllLines.end();
	

	/**
	* 输出图
	*/
#ifdef _COLOR_
	unsigned int* pImageData = new unsigned int[nWidth*nHeight];
	memset(pImageData, 0, nWidth*nHeight*sizeof(unsigned int));
#else
	float* pImageData = new float[nWidth*nHeight];
	memset(pImageData, 0, nWidth*nHeight*sizeof(float));
#endif

	for (; itLineBegin!=itLineEnd; itLineBegin++)
	{
		int iLineIdx = itLineBegin->first;
		MPoints& pts = itLineBegin->second;
		if(pts.size() > nLengthThreholdLow && pts.size() < nLengthThreHoldHigh)
		{
			//printf("Line %d with %d elements\n", iLineIdx, pts.size());
			MPoints::iterator itPtBegin = pts.begin(), itPtEnd = pts.end();
#ifdef _COLOR_
			uint32 uColor = g_RandomGen.Generate();
#else
			float uColor = 1.0;
#endif
			for (;itPtBegin != itPtEnd; itPtBegin++)
			{
				int iDx = itPtBegin->first;
				pImageData[iDx] = uColor;
			}
		}
	}

#ifdef _COLOR_
	cutSavePPM4ub(fileSaveName, (unsigned char*)pImageData, nWidth, nHeight);
#else
	//savePGM(fileSaveName, (float*)pImageData, nWidth, nHeight);
	cutSavePGMf(fileSaveName, pImageData, nWidth, nHeight);
#endif

}

//////////////////////////////////////////////////////////////////////////


