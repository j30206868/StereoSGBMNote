#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/highgui/highgui.hpp"

#include "stereoSGBM.h"

#include <iostream>

using namespace cv;
int main(int argc, char ** argv)
{
	Mat L = imread("im2.png", CV_LOAD_IMAGE_GRAYSCALE), R = imread("im6.png", CV_LOAD_IMAGE_GRAYSCALE);
	Mat disp;
	StereoSGBMParams params;
	int cn = 1;

	params.preFilterCap = 63;
	params.SADWindowSize = 5;
	params.P1 = 8 * cn * params.SADWindowSize*params.SADWindowSize;
	params.P2 = 32 * cn * params.SADWindowSize*params.SADWindowSize;
	params.minDisparity = 0;
	params.numDisparities = 64;
	params.uniquenessRatio = 10;
	params.speckleWindowSize = 100;
	params.speckleRange = 128;
	params.disp12MaxDiff = 10;

	int64 t = getTickCount();
	compute(L, R, disp, params);
	t = getTickCount() - t;
	printf("StereoSGBM comput(): number of disparity: %d | time: %fms\n", params.numDisparities, t * 1000 / getTickFrequency());

	Mat visibleMap;
	disp.convertTo(visibleMap, CV_8U, 255 / (params.numDisparities*16.0));

	//double ratio = 255 / params.numDisparities;
	//for (int i = 0; i < disp.cols*disp.rows; i++)
	//	disp.data[i] *= ratio;

	namedWindow("disp", CV_WINDOW_AUTOSIZE);
	imshow("disp", visibleMap);
	waitKey(0);

	system("PAUSE");
	return 0;
}
