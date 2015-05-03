// Filename:
//    cv_mouse1.cpp 
// Purpose:
//    Display mouse position while moving on the picture
// by ausk ，20150503
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
using namespace std;
using namespace cv;

IplImage* src=0;  
IplImage* dst=0;  
void on_mouse( int event, int x, int y, int flags, void* ustc)  
{  
    char temp[16];  
    CvPoint pt;  
    CvFont font;  
    cvInitFont(&font, CV_FONT_HERSHEY_SIMPLEX, 0.5, 0.5, 0, 1, CV_AA);  
          
    if( event == CV_EVENT_MOUSEMOVE )  
    {  
        cvCopy(dst,src);          
        sprintf(temp,"(%d,%d)",x,y);  
        pt = cvPoint(x,y);  
        cvPutText(src,temp, pt, &font, cvScalar(255, 255, 255, 0));  
        cvCircle( src, pt, 2,cvScalar(255,0,0,0) ,CV_FILLED, CV_AA, 0 );  
        cvShowImage( "src", src );  
    }   
    else if( event == CV_EVENT_LBUTTONDOWN )  
    {  
        //cvCopy(dst,src);            
        sprintf(temp,"(%d,%d)",x,y);  
        pt = cvPoint(x,y);  
        cvPutText(src,temp, pt, &font, cvScalar(255, 255, 255, 0));  
        cvCircle( src, pt, 2,cvScalar(255,0,0,0) ,CV_FILLED, CV_AA, 0 );  
        cvCopy(src,dst);  
        cvShowImage( "src", src );  
    }   
}  
  
int main(int argc, char** argv)
{ 
    src=cvLoadImage("E:\\images\\lena.bmp");
    if(!src)
    {
    	printf("Cannot read the image!nPlease check it!\n");
    	return 0;
    }  
    dst=cvCloneImage(src);  
  
    cvNamedWindow("src",1);  
    cvSetMouseCallback( "src", on_mouse, 0 );  
      
    cvShowImage("src",src);  
    cvWaitKey(0);   
    cvDestroyAllWindows();  
    cvReleaseImage(&src);  
    cvReleaseImage(&dst);  
  
    return 0;  
}  