#include <opencv2\opencv.hpp>
#include <opencv2\calib3d\calib3d.hpp>
using namespace cv;

//copy stereoSGBM setting
typedef uchar PixType;
typedef short CostType;
typedef short DispType;
enum { NR = 16, NR2 = NR / 2 };
#define DISP_SHIFT 4

#define MODE_SGBM 2
#define MODE_HH 3


struct StereoSGBMParams
{
	StereoSGBMParams()
	{
		minDisparity = numDisparities = 0;
		SADWindowSize = 0;
		P1 = P2 = 0;
		disp12MaxDiff = 0;
		preFilterCap = 0;
		uniquenessRatio = 0;
		speckleWindowSize = 0;
		speckleRange = 0;
		//mode = StereoSGBM::MODE_SGBM;
		mode = MODE_SGBM;
	}

	StereoSGBMParams(int _minDisparity, int _numDisparities, int _SADWindowSize,
		int _P1, int _P2, int _disp12MaxDiff, int _preFilterCap,
		int _uniquenessRatio, int _speckleWindowSize, int _speckleRange,
		int _mode)
	{
		minDisparity = _minDisparity;
		numDisparities = _numDisparities;
		SADWindowSize = _SADWindowSize;
		P1 = _P1;
		P2 = _P2;
		disp12MaxDiff = _disp12MaxDiff;
		preFilterCap = _preFilterCap;
		uniquenessRatio = _uniquenessRatio;
		speckleWindowSize = _speckleWindowSize;
		speckleRange = _speckleRange;
		mode = _mode;
	}

	int minDisparity;
	int numDisparities;
	int SADWindowSize;
	int preFilterCap;
	int uniquenessRatio;
	int P1;
	int P2;
	int speckleWindowSize;
	int speckleRange;
	int disp12MaxDiff;
	int mode;
};

void compute(InputArray leftarr, InputArray rightarr, OutputArray disparr, StereoSGBMParams params);