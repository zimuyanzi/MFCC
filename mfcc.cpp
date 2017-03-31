#include "mfcc.h"
#include "math.h"
#include "log4cplus/XtLogger.h"

using namespace log4cplus;
using namespace log4cplus_common;

static const double dPreFactor = 0.9375;  // 预加重因子
static const double dPI        = 3.141592653589; // 圆周率
static const double dFH        = 0.5;
static const double dFL        = 0;

#define CONVERT_AUDIO_TYPE(type)                                           \
	nPointNum = nLen / sizeof(type);                                       \
	type *pNewAudio = (type *)pAudio;                                      \
	for (nIdx = 0; nIdx < nPointNum && nIdx < m_nPointNum; nIdx ++)        \
	{                                                                      \
		m_pAudioData[nIdx] = pNewAudio[nIdx];                              \
	}                                                                      \
	for (nIdx = nPointNum; nIdx < m_nPointNum; nIdx ++)                    \
	{                                                                      \
		m_pAudioData[nIdx] = 0;                                            \
	}                                                                      \

CMFCC::CMFCC(int nPointNum, int nPointByte, int nSampleRate,
	int nDCTX, int nDCTY)
{
	if (is2Power(nPointNum) == false)
	{
		m_nPointNum = m_nFFTPointNumDefalut;
	}
	else
	{
		m_nPointNum = nPointNum;
	}
	m_nPointByte    = nPointByte;
	m_nSampleRate   = nSampleRate;
	m_nDCTX         = nDCTX;
	m_nDCTY         = nDCTY;

	/** 音频数据 */
	m_pAudioData    = new double[m_nPointNum];

	/** hamming窗 */
	m_pHamming      = new double[m_nPointNum];

	/** 梅尔矩阵 */
	m_pMelMatrix    = new double *[m_nDCTY];
	for (int i = 0; i < m_nDCTY; i ++)
	{
		m_pMelMatrix[i] = new double[m_nPointNum];
	}

	/** DCT矩阵 */
	m_pDCT = new double *[m_nDCTX];
	for (int i = 0; i < m_nDCTX; i ++)
	{
		m_pDCT[i]   = new double[m_nDCTY];
	}

	/** 倒谱 */
	m_pStrenWin        = new double[m_nDCTX];

	m_pFFTIn           = NULL;
	m_pFFTOut          = NULL;
}

CMFCC::~CMFCC()
{
	if (m_pAudioData)
	{
		delete [] m_pAudioData;
		m_pAudioData = NULL;
	}
	if (m_pHamming)
	{
		delete [] m_pHamming;
		m_pHamming = NULL;
	}
	for (int i = 0; i < m_nDCTY; i ++)
	{
		if (m_pMelMatrix[i])
		{
			delete [] m_pMelMatrix[i];
			m_pMelMatrix[i] = NULL;
		}
	}
	if (m_pMelMatrix)
	{
		delete [] m_pMelMatrix;
		m_pMelMatrix = NULL;
	}
	for (int i = 0; i < m_nDCTX; i ++)
	{
		if (m_pDCT[i])
		{
			delete [] m_pDCT[i];
			m_pDCT[i] = NULL;
		}
	}
	if (m_pDCT)
	{
		delete [] m_pDCT;
		m_pDCT = NULL;
	}
	if (m_pStrenWin)
	{
		delete [] m_pStrenWin;
		m_pStrenWin = NULL;
	}
}

int CMFCC::init(fft_plan pFFTPlan, fft_complex *pFFTIn, fft_complex *pFFTOut)
{
	if (pFFTPlan == NULL || pFFTIn == NULL || pFFTOut == NULL)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CMFCC::init] param error");
		return -1;
	}

	/** 计算梅尔系数 */
	calMelCoefficient();
	/** 计算hamming窗 */
	calHamming();
	/** 归一化倒谱提升窗口 */
	calStrenWin();
	/** 计算DCT */
	calDCT();

	/** 创建fft相关的 */
	m_fftPlan = pFFTPlan;
	m_pFFTIn = pFFTIn;
	m_pFFTOut = pFFTOut;

	return 0;
}

int CMFCC::calMelCoefficient()
{
	int filterNum = m_nDCTY, frameNum = m_nPointNum, frameSize = m_nPointNum;
	double freMax = m_nSampleRate / 2;//实际最大频率
	double freMin = 0;//实际最小频率
	double melFremax = 1125 * log(1 + freMax / 700);//将实际频率转换成梅尔频率
	double melFremin = 1125 * log(1 + freMin / 700);
	double k = (melFremax - melFremin) / (filterNum + 1);

	double *m = new double[filterNum + 2];
	double *h = new double[filterNum + 2];
	double *f = new double[filterNum + 2];

	for(int i = 0; i < filterNum + 2; i++)
	{
		m[i] = melFremin + k * i;
		h[i] = 700 * (exp(m[i] / 1125) - 1);//将梅尔频率转换成实际频率
		f[i] = floor((frameSize + 1) * h[i] / m_nSampleRate);
	}

	if (m)
	{
		delete [] m;
		m = NULL;
	}
	if (h)
	{
		delete [] h;
		h = NULL;
	}

	/* 初始化 */
	for(int j = 0; j < filterNum; j++)
	{
		for(int k = 0; k < frameNum; k++)
		{
			m_pMelMatrix[j][k] = 0;
		}
	}

	//计算出每个三角滤波器
	for(int j = 0; j < filterNum; j++)
	{
		for(int k = 0; k < frameNum; k++)
		{
			double temp  = 0;
			if ( k < f[j] )
			{
				temp = 0;
			}
			else if (k >= f[j] && k <= f[j + 1])
			{
				if (f[j + 1] - f[j] >= -0.00000000001 && f[j + 1] - f[j] <= 0.00000000001)
				{
					temp = 0;
				}
				else
				{
					temp = (k - f[j]) / (f[j + 1] - f[j]);
				}
			}
			else if (k >= f[j + 1] && k <= f[j + 2])
			{
				if (f[j + 1] - f[j + 2] >= -0.00000000001 && f[j + 1] - f[j + 2] <= 0.00000000001)
				{
					temp = 0;
				}
				else
				{
					temp = (f[j + 2] - k) / (f[j + 2] - f[j + 1]);
				}
			}
			else if (k >= f[j + 2])
			{
				temp = 0;
			}
			m_pMelMatrix[j][k] = temp;
		}
	}

	if (f)
	{
		delete[] f;
		f = NULL;
	}
	return 0;
}

int CMFCC::calHamming()
{
	if (m_nPointNum <= 1)
	{
		return -1;
	}

	for (int i = 0; i < m_nPointNum; i++)
	{
		m_pHamming[i] = (double)(0.54 - 0.46 * cos(2 * dPI * (double)i / ((double)m_nPointNum - 1) ));
	}
	return 0;
}

int CMFCC::calStrenWin()
{
	int i;
	double b = 0.0;
	for(i = 0; i < m_nDCTX; i ++)
	{
		m_pStrenWin[i] = 1 + 6*sin(dPI*(double)(i + 1)/(double)(m_nDCTX));
		if(b < m_pStrenWin[i])
			b = m_pStrenWin[i];
	}
	if (b < 0.000000000001 && b > -0.000000000001)
	{
		return -1;
	}
	for(i = 0; i < m_nDCTX; i ++)
	{
		m_pStrenWin[i] = m_pStrenWin[i]/b;
	}
	return 0;
}

int CMFCC::calDCT()
{
	int nX = 0, nY = 0;
	for(nX = 0; nX < m_nDCTX; nX ++)
	{
		for(nY = 0; nY < m_nDCTY; nY ++)
		{
			m_pDCT[nX][nY] = cos((double)(2*nY + 1)*nX*dPI/(double)(2*m_nDCTY));
		}
	}

	return 0;
}

int CMFCC::preEmphasizing (double *pSample, int nLen, double dFactor)
{
	if (pSample == NULL || nLen < 0)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CMFCC::preEmphasizing] param error");
		return -1;
	}
	for (int i = nLen - 1; i >= 1; i --)
	{
		pSample[i] = pSample[i] - dFactor * pSample[i - 1];
	}
	return 0;
}

int CMFCC::pushAudioData (const unsigned char *pAudio, const int nLen)
{
	if (nLen % m_nPointByte != 0)
	{
		LOG4CPLUS_DEBUG(g_oServerLogger, "[CMFCC::pushAudioData] nLen % m_nPointByte != 0"
			<< ", nLen :" << nLen
			<< ", m_nPointByte :" << m_nPointByte);
		return -1;
	}

	int nPointNum = 0, nIdx = 0;

	if (m_nPointByte == sizeof(short))
	{
		CONVERT_AUDIO_TYPE(short);
	}
	else if(m_nPointByte == sizeof(int))
	{
		CONVERT_AUDIO_TYPE(int);
	}
	else if (m_nPointByte == sizeof(double))
	{
		CONVERT_AUDIO_TYPE(double);
	}
	else
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CMFCC::pushAudioData] data type not support!");
		return -1;
	}
	return 0;
}

int CMFCC::extractMFCCFeature (const unsigned char *pAudio, const int nLen, vector<double> &featureVec)
{
	int nRect = 0;

	/** 转换数据类型 */
	nRect = pushAudioData(pAudio, nLen);
	if (nRect != 0)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CMFCC::extractMFCCFeature] pushAudioData failed!");
		return -1;
	}

	/** 预加重 */
	nRect = preEmphasizing(m_pAudioData, m_nPointNum, dPreFactor);
	if (nRect != 0)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CMFCC::extractMFCCFeature] preEmphasizing failed!");
		return -1;
	}

	/** 加窗 */
	for (int i = 0; i < m_nPointNum; i ++)
	{
		m_pAudioData[i] = m_pAudioData[i] * m_pHamming[i];
	}

	//  矩阵乘积
	//  13x24  24xN   Nx1      13x1
	//  DCT *  Mel *  FFT(x) = MFCC

	/** FFT */
	nRect = FFTW3(m_pAudioData, m_nPointNum);
	if (nRect != 0)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CMFCC::extractMFCCFeature] FFTW3 failed");
		return -1;
	}

	/** mel系数 */
	double *pTempData = new double[m_nDCTY];
	double dSum = 0;
	for (int i = 0; i < m_nDCTY; i ++)
	{
		dSum = 0;
		for (int j = 0; j < (m_nPointNum / 2 + 1); j ++)
		{
			dSum += m_pAudioData[j] * m_pMelMatrix[i][j];
		}
		if (dSum != 0)
		{
			pTempData[i] = log(dSum);
		}
		else
		{
			pTempData[i] = -3000;
		}
	}

	/** DCT */
	for (int i = 1; i < m_nDCTX; i ++)
	{
		dSum = 0;
		for (int j = 0; j < m_nDCTY; j ++)
		{
			dSum += m_pDCT[i][j] * pTempData[j];
		}
		double dResult = dSum * m_pStrenWin[i];
		featureVec.push_back(dResult);
	}

	if (pTempData)
	{
		delete [] pTempData;
		pTempData = NULL;
	}

	if (featureVec.size() == 0)
	{
		return -1;
	}

	return 0;
}

bool CMFCC::is2Power (const int nLen)
{
	return nLen > 0 ? !(nLen & (nLen - 1)) : false;
}

int CMFCC::FFTW3 (double *pSample, int nPointNum)
{
	if (pSample == NULL)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CMFCC::FFTW3] pSample == NULL");
		return -1;
	}

	for (int i = 0; i < m_nPointNum; i ++)
	{
		m_pFFTIn[i][0] = (float)pSample[i];
		m_pFFTIn[i][1] = 0;
	}

	fft_execute(m_fftPlan);

	for (int i = 0; i < m_nPointNum; i ++)
	{
		pSample[i] = m_pFFTOut[i][0]*m_pFFTOut[i][0] + m_pFFTOut[i][1]*m_pFFTOut[i][1];
	}
	return 0;
}

int CMFCC::reset()
{
	for (int i = 0; i < m_nPointNum; i ++)
	{
		m_pFFTIn[i][0] = 0;
		m_pFFTIn[i][1] = 0;
		m_pFFTOut[i][0] = 0;
		m_pFFTOut[i][1] = 0;
	}

	return 0;
}

int CMFCC::extractFrequencyFeature(const unsigned char *pAudio, const int nLen, vector< double > &featureVec)
{
	int nRect = 0;

	/** 转换数据类型 */
	nRect = pushAudioData(pAudio, nLen);
	if (nRect != 0)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CMFCC::extractMFCCFeature] pushAudioData failed!");
		return -1;
	}

	/** 预加重 */
	nRect = preEmphasizing(m_pAudioData, m_nPointNum, dPreFactor);
	if (nRect != 0)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CMFCC::extractMFCCFeature] preEmphasizing failed!");
		return -1;
	}

	/** 加窗 */
	for (int i = 0; i < m_nPointNum; i ++)
	{
		m_pAudioData[i] = m_pAudioData[i] * m_pHamming[i];
	}

	//  矩阵乘积
	//  13x24  24xN   Nx1      13x1
	//  DCT *  Mel *  FFT(x) = MFCC

	/** FFT */
	nRect = FFTW3(m_pAudioData, m_nPointNum);
	if (nRect != 0)
	{
		LOG4CPLUS_WARN(g_oServerLogger, "[CMFCC::extractMFCCFeature] FFTW3 failed");
		return -1;
	}

	/** mel系数 */
	double *pTempData = new double[m_nDCTY];
	for (int i = 0; i < m_nDCTY; i ++)
	{
		double dSum = 0;
		for (int j = 0; j < (m_nPointNum / 2 + 1); j ++)
		{
			dSum += m_pAudioData[j] * m_pMelMatrix[i][j];
		}
		if (dSum != 0)
		{
			pTempData[i] = log(dSum);
		}
		else
		{
			pTempData[i] = 0;
		}
		featureVec.push_back(pTempData[i]);
	}

	if (pTempData)
	{
		delete [] pTempData;
		pTempData = NULL;
	}

	return 0;
}