#ifndef __MFCC_H_FILE__
#define __MFCC_H_FILE__

#include <list>
#include <vector>
#include "fftw/fftw3.h"
#include "vdl/stdafx.h"

using namespace std;

//#define __FFTW_F__

/** 根据链接库的不同调用不同的函数 */
#ifdef __FFTW_F__
	/** 单精度 */
	#define fft_make_planner_thread_safe() fftwf_make_planner_thread_safe()
	#define fft_destroy_plan(a) fftwf_destroy_plan(a)
	#define fft_plan_dft_1d(a, b, c, d, e)  fftwf_plan_dft_1d(a, b, c, d, e)
	#define fft_free(a) fftwf_free(a)
	#define fft_complex fftwf_complex
	#define fft_plan fftwf_plan
	#define fft_malloc fftwf_malloc
	#define fft_execute fftwf_execute
#else
	#define fft_make_planner_thread_safe() fftw_make_planner_thread_safe()
	#define fft_destroy_plan(a) fftw_destroy_plan(a)
	#define fft_plan_dft_1d(a, b, c, d, e)  fftw_plan_dft_1d(a, b, c, d, e)
	#define fft_free(a) fftw_free(a)
	#define fft_complex fftw_complex
	#define fft_plan fftw_plan
	#define fft_malloc fftw_malloc
	#define fft_execute fftw_execute
#endif

class CMFCC
{
public:
	CMFCC(int nPointNum, int nPointByte, int nSampleRate,
		int nDCTX = m_nDCTXDefalut, int nDCTY = m_nDCTYDefalut);
	CMFCC();
	CMFCC(CMFCC &);
	~CMFCC();
public:
	/** 提取mfcc特征 */
	int extractMFCCFeature (const unsigned char *pAudio, const int nLen, vector< double > &featureVec);
	/** 提取特征，只有频率特征 */
	int extractFrequencyFeature(const unsigned char *pAudio, const int nLen, vector< double > &featureVec);
	/** 复位 */
	int reset();
	/** 初始化相关矩阵 */
	int init(fft_plan pFFTPlan, fft_complex *pFFTIn, fft_complex *pFFTOut);

private:
	/** 计算梅尔系数 */
	int calMelCoefficient();
	/** 计算hamming窗 */
	int calHamming();
	/** 归一化倒谱提升窗口 */
	int calStrenWin();
	/** 计算DCT系数 */
	int calDCT();
	/** 预加重 */
	int preEmphasizing (double *pSample, int nLen, double dFactor);
	/** 转换数据类型 */
	int pushAudioData (const unsigned char *pAudio, const int nLen);
	/** 判断是否是2的指数 */
	bool is2Power (const int nLen);
	/** 快速傅里叶变换 */
	int FFTW3 (double *pSample, int nPointNum);

private:
	/** 参数阈值 */
	static const int m_nFFTPointNumDefalut     = 1024;          /** 默认FFT变换点数 */
	static const int m_nSampleRateDefalut      = 48000;         /** 默认采样率 */
	static const int m_nByteDefalut            = 2;             /** 默认每个采样点占用字节 */
	static const int m_nDCTXDefalut            = 13;            /** 默认DCT变换矩阵行数 */
	static const int m_nDCTYDefalut            = 24;            /** 默认DCT变换矩阵列数 */

private:

	int              m_nPointNum;                               /** 采样点数 */
	int              m_nPointByte;                              /** 每个采样点所占字节 */
	int              m_nSampleRate;                             /** 采样率 */
	double *         m_pAudioData;                              /** 音频数据 */

	double **        m_pMelMatrix;                              /** 梅尔矩阵 */
	double *         m_pHamming;                                /** hamming窗 */
	double **        m_pDCT;                                    /** DCT矩阵 */
	double *         m_pStrenWin;                               /** 归一化倒谱提升 */
	int              m_nDCTX;                                   /** DCT变换矩阵行数 */
	int              m_nDCTY;                                   /** DCT变换矩阵列数 */

	fft_plan         m_fftPlan;                                 /** fft plan */
	fft_complex *    m_pFFTIn;                                  /** fft输入数据 */
	fft_complex *    m_pFFTOut;                                 /** fft输出数据 */
};

#endif//#ifndef __MFCC_H_FILE__