#include <opencv2\highgui\highgui.hpp>
#include <math.h>
#include <string.h>
#include <iostream>
using namespace std;

//二维离散变换小波源代码
void DWT(IplImage *pImage, int nLayer)
{
   // 执行条件
   if (pImage)
   {
      if (pImage->nChannels == 1 &&
         pImage->depth == IPL_DEPTH_32F &&
         ((pImage->width >> nLayer) << nLayer) == pImage->width &&
         ((pImage->height >> nLayer) << nLayer) == pImage->height)
         //单通道,32位位深,宽高是n^nLayer
      {
         int     i, x, y, n;
         float   fValue   = 0;
         float   fRadius  = sqrt(2.0f);
         int     nWidth   = pImage->width;
         int     nHeight  = pImage->height;
         int     nHalfW   = nWidth / 2;
         int     nHalfH   = nHeight / 2;
         float **pData    = new float*[pImage->height];
         float  *pCols     = new float[pImage->width];
         float  *pRows  = new float[pImage->height];

         //创建临时矩阵，保存图像处理数据
         for (i = 0; i < pImage->height; i++)
            pData[i] = (float*) (pImage->imageData + pImage->widthStep * i);



         // 多层小波变换
         for (n = 0; n < nLayer; n++, nWidth /= 2, nHeight /= 2, nHalfW /= 2, nHalfH /= 2)
         {
            // 水平变换
            for (y = 0; y < nHeight; y++)
            {
               // 奇偶分离

               //函数原型 void *memcpy(void *dest, const void *src, size_t n);
               memcpy(pCols, pData[y], sizeof(float) * nWidth);

               for (i = 0; i < nHalfW; i++)
               {
                  x = i * 2;
                  pData[y][i] = pCols[x];
                  pData[y][nHalfW + i] = pCols[x + 1];
               }

               /* / 提升小波变换
               //////////////////////////////
               for (i = 0; i < nHalfW - 1; i++)
               {
                  fValue = (pData[y][i] + pData[y][i + 1]) / 2;
                  pData[y][nHalfW + i] -= fValue;
               }

               fValue = (pData[y][nHalfW - 1] + pData[y][nHalfW - 2]) / 2;
               pData[y][nWidth - 1] -= fValue;
               fValue = (pData[y][nHalfW] + pData[y][nHalfW + 1]) / 4;
               pData[y][0] += fValue;

               for (i = 1; i < nHalfW; i++)
               {
                  fValue = (pData[y][nHalfW + i] + pData[y][nHalfW + i - 1]) / 4;
                  pData[y][i] += fValue;
               }
               */ ///////////////////////////////////////////

               // 频带系数
               for (i = 0; i < nHalfW; i++)
               {
                  pData[y][i] *= fRadius;
                  pData[y][nHalfW + i] /= fRadius;
               }
            }

            // 垂直变换
            for (x = 0; x < nWidth; x++)
            {
               // 奇偶分离
               for (i = 0; i < nHalfH; i++)
               {
                  y = i * 2;
                  pRows[i] = pData[y][x];
                  pRows[nHalfH + i] = pData[y + 1][x];
               }
               for (i = 0; i < nHeight; i++)
               {
                  pData[i][x] = pRows[i];
               }
               /*/ 提升小波变换
               ///////////////////////////////////////////
               for (i = 0; i < nHalfH - 1; i++)
               {
                  fValue = (pData[i][x] + pData[i + 1][x]) / 2;
                  pData[nHalfH + i][x] -= fValue;
               }
               fValue = (pData[nHalfH - 1][x] + pData[nHalfH - 2][x]) / 2;
               pData[nHeight - 1][x] -= fValue;
               fValue = (pData[nHalfH][x] + pData[nHalfH + 1][x]) / 4;
               pData[0][x] += fValue;
               for (i = 1; i < nHalfH; i++)
               {
                  fValue = (pData[nHalfH + i][x] + pData[nHalfH + i - 1][x]) / 4;
                  pData[i][x] += fValue;
               }
               */ //////////////////////////////////////////////////////

               // 频带系数
               for (i = 0; i < nHalfH; i++)
               {
                  pData[i][x] *= fRadius;
                  pData[nHalfH + i][x] /= fRadius;
               }
            }
         }
         delete[] pData;
         delete[] pCols;
         delete[] pRows;
      }
   }
}

// 二维离散小波恢复（单通道浮点图像）
void IDWT(IplImage *pImage, int nLayer)
{
   // 执行条件
   if (pImage)
   {
      if (pImage->nChannels == 1 &&
         pImage->depth == IPL_DEPTH_32F &&
         ((pImage->width >> nLayer) << nLayer) == pImage->width &&
         ((pImage->height >> nLayer) << nLayer) == pImage->height)
      {
         int     i, x, y, n;
         float   fValue   = 0;
         float   fRadius  = sqrt(2.0f);
         int     nWidth   = pImage->width >> (nLayer - 1);
         int     nHeight  = pImage->height >> (nLayer - 1);
         int     nHalfW   = nWidth / 2;
         int     nHalfH   = nHeight / 2;
         float **pData    = new float*[pImage->height];
         float  *pCols     = new float[pImage->width];
         float  *pRows  = new float[pImage->height];
         for (i = 0; i < pImage->height; i++)
         {
            pData[i] = (float*) (pImage->imageData + pImage->widthStep * i);
         }
         // 多层小波恢复
         for (n = 0; n < nLayer; n++, nWidth *= 2, nHeight *= 2, nHalfW *= 2, nHalfH *= 2)
         {
            // 垂直恢复
            for (x = 0; x < nWidth; x++)
            {
               // 频带系数
               for (i = 0; i < nHalfH; i++)
               {
                  pData[i][x] /= fRadius;
                  pData[nHalfH + i][x] *= fRadius;
               }
               /* / 提升小波恢复
               /////////////////////////////////////////////////////
               fValue = (pData[nHalfH][x] + pData[nHalfH + 1][x]) / 4;
               pData[0][x] -= fValue;
               for (i = 1; i < nHalfH; i++)
               {
                  fValue = (pData[nHalfH + i][x] + pData[nHalfH + i - 1][x]) / 4;
                  pData[i][x] -= fValue;
               }
               for (i = 0; i < nHalfH - 1; i++)
               {
                  fValue = (pData[i][x] + pData[i + 1][x]) / 2;
                  pData[nHalfH + i][x] += fValue;
               }
               fValue = (pData[nHalfH - 1][x] + pData[nHalfH - 2][x]) / 2;
               pData[nHeight - 1][x] += fValue;
               */ ////////////////////////////////////////////////
               // 奇偶合并
               for (i = 0; i < nHalfH; i++)
               {
                  y = i * 2;
                  pRows[y] = pData[i][x];
                  pRows[y + 1] = pData[nHalfH + i][x];
               }
               for (i = 0; i < nHeight; i++)   pData[i][x] = pRows[i];
            }
            // 水平恢复
            for (y = 0; y < nHeight; y++)
            {
               // 频带系数
               for (i = 0; i < nHalfW; i++)
               {
                  pData[y][i] /= fRadius;
                  pData[y][nHalfW + i] *= fRadius;
               }
               /* / 提升小波恢复
               //////////////////////////////////////////////////
               fValue = (pData[y][nHalfW] + pData[y][nHalfW + 1]) / 4;
               pData[y][0] -= fValue;
               for (i = 1; i < nHalfW; i++)
               {
                  fValue = (pData[y][nHalfW + i] + pData[y][nHalfW + i - 1]) / 4;
                  pData[y][i] -= fValue;
               }
               for (i = 0; i < nHalfW - 1; i++)
               {
                  fValue = (pData[y][i] + pData[y][i + 1]) / 2;
                  pData[y][nHalfW + i] += fValue;
               }
               fValue = (pData[y][nHalfW - 1] + pData[y][nHalfW - 2]) / 2;
               pData[y][nWidth - 1] += fValue;
               */ ///////////////////////////////////////////////////
               // 奇偶合并
               for (i = 0; i < nHalfW; i++)
               {
                  x = i * 2;
                  pCols[x] = pData[y][i];
                  pCols[x + 1] = pData[y][nHalfW + i];
               }
               memcpy(pData[y], pCols, sizeof(float) * nWidth);
            }
         }
         delete[] pData;
         delete[] pCols;
         delete[] pRows;
      }
   }
} 


int main()
{
	// 小波变换层数
   int nLayer = 2;
   const char* path = "E://images//lena.bmp";

	// 输入彩色图像
	IplImage *pSrc = cvLoadImage(path, 0); //读取为黑白图像
	// 计算小波图象大小
	CvSize size = cvGetSize(pSrc);
	if ((pSrc->width >> nLayer) << nLayer != pSrc->width)	
      size.width = ((pSrc->width >> nLayer) + 1) << nLayer;
	if ((pSrc->height >> nLayer) << nLayer != pSrc->height)
     size.height = ((pSrc->height >> nLayer) + 1) << nLayer;
  
	// 创建小波图象
	IplImage *pWavelet = cvCreateImage(size, IPL_DEPTH_32F, pSrc->nChannels);
	if (pWavelet) {
		// 小波图象赋值
		cvSetImageROI(pWavelet, cvRect(0, 0, pSrc->width, pSrc->height));
		cvConvertScale(pSrc, pWavelet, 1, -128);
		cvResetImageROI(pWavelet);
		// 彩色图像小波变换
		IplImage *pImage = cvCreateImage(cvGetSize(pWavelet), IPL_DEPTH_32F, 1);
		if (pImage) {
			for (int i = 1; i <= pWavelet->nChannels; i++)   
         {
				cvSetImageCOI(pWavelet, i);
				cvCopy(pWavelet, pImage, NULL);
				// 二维离散小波变换
				DWT(pImage, nLayer);
				// 二维离散小波恢复
				//IDWT(pImage, nLayer);
				cvCopy(pImage, pWavelet, NULL);

			}
			cvSetImageCOI(pWavelet, 0);
			cvReleaseImage(&pImage);
		}
		// 小波变换图象
		cvSetImageROI(pWavelet, cvRect(0, 0, pSrc->width, pSrc->height));
		cvConvertScale(pWavelet, pSrc, 1, 128);
		cvResetImageROI(pWavelet);
		cvReleaseImage(&pWavelet);
	}
	// 显示图像pSrc
	cv::namedWindow("result");
	cvShowImage("result",pSrc);
	cv::waitKey(0);
	cvReleaseImage(&pSrc); 
}