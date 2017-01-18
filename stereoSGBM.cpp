#include "stereoSGBM.h"

#include <stdio.h>

static void calcPixelCostBT(const Mat& img1, const Mat& img2, int y,
	int minD, int maxD, CostType* cost,
	PixType* buffer, const PixType* tab,
	int tabOfs, int)
{
	int x, c, width = img1.cols, cn = img1.channels();
	int minX1 = std::max(maxD, 0), maxX1 = width + std::min(minD, 0);
	int minX2 = std::max(minX1 - maxD, 0), maxX2 = std::min(maxX1 - minD, width);
	int D = maxD - minD, width1 = maxX1 - minX1, width2 = maxX2 - minX2;
	const PixType *row1 = img1.ptr<PixType>(y), *row2 = img2.ptr<PixType>(y);
	PixType *prow1 = buffer + width2 * 2, *prow2 = prow1 + width*cn * 2;

	/*
	**	查表	  **
	TAB_OFS = 1024
	TAB_SIZE = 2304
	clipTab[0] ~ clipTab[961] = 0
	clipTab[962] = 1
	clipTab[963] = 2
	clipTab[964] = 3
	......
	clipTab[1087] = 126
	clipTab[1088] = 126
	......
	clipTab[2303] = 126
	*/

	tab += tabOfs;

	for (c = 0; c < cn * 2; c++)
	{
		prow1[width*c] = prow1[width*c + width - 1] =
			prow2[width*c] = prow2[width*c + width - 1] = tab[0];
	}

	int n1 = y > 0 ? -(int)img1.step : 0, s1 = y < img1.rows - 1 ? (int)img1.step : 0;
	int n2 = y > 0 ? -(int)img2.step : 0, s2 = y < img2.rows - 1 ? (int)img2.step : 0;

	/*
	n1, n2 第一列  (y軸) = 0, 其他 = (-) width * cn
		n1,n2用來存取上一個row的資料 row1[x + n1) => (x, y-1)
	s1, s2 最後一列(y軸) = 0, 其他 = (+) width * cn
		s1,s2用來存取下一個row的資料 row1[x + s1) => (x, y+1)
	*/
	if (cn == 1)
	{	//灰階影像
		for (x = 1; x < width - 1; x++)
		{
			//左影像 位置 x 的 pixDiff (查表得到)
			prow1[x] = tab[(row1[x + 1] - row1[x - 1]) * 2 + row1[x + n1 + 1] - row1[x + n1 - 1] + row1[x + s1 + 1] - row1[x + s1 - 1]];
			//右影像 位置 x (從右邊往回數, 因為是 width-1-x) 的 pixDiff
			prow2[width - 1 - x] = tab[(row2[x + 1] - row2[x - 1]) * 2 + row2[x + n2 + 1] - row2[x + n2 - 1] + row2[x + s2 + 1] - row2[x + s2 - 1]];

			//prow為每個pixel開了兩格的空間
			//以 x 來說, prow1[x] 放的是 x 的 cost, prow1[x+width] 放的是 x 的 像素值
			
			//存放 左影像 位置 x 的像素到 prow1
			prow1[x + width] = row1[x];
			//存放 右影像 位置 x 的像素到 prow2
			prow2[width - 1 - x + width] = row2[x];
		}
	}
	else
	{	//彩色影像
		for (x = 1; x < width - 1; x++)
		{
			prow1[x] = tab[(row1[x * 3 + 3] - row1[x * 3 - 3]) * 2 + row1[x * 3 + n1 + 3] - row1[x * 3 + n1 - 3] + row1[x * 3 + s1 + 3] - row1[x * 3 + s1 - 3]];
			prow1[x + width] = tab[(row1[x * 3 + 4] - row1[x * 3 - 2]) * 2 + row1[x * 3 + n1 + 4] - row1[x * 3 + n1 - 2] + row1[x * 3 + s1 + 4] - row1[x * 3 + s1 - 2]];
			prow1[x + width * 2] = tab[(row1[x * 3 + 5] - row1[x * 3 - 1]) * 2 + row1[x * 3 + n1 + 5] - row1[x * 3 + n1 - 1] + row1[x * 3 + s1 + 5] - row1[x * 3 + s1 - 1]];

			prow2[width - 1 - x] = tab[(row2[x * 3 + 3] - row2[x * 3 - 3]) * 2 + row2[x * 3 + n2 + 3] - row2[x * 3 + n2 - 3] + row2[x * 3 + s2 + 3] - row2[x * 3 + s2 - 3]];
			prow2[width - 1 - x + width] = tab[(row2[x * 3 + 4] - row2[x * 3 - 2]) * 2 + row2[x * 3 + n2 + 4] - row2[x * 3 + n2 - 2] + row2[x * 3 + s2 + 4] - row2[x * 3 + s2 - 2]];
			prow2[width - 1 - x + width * 2] = tab[(row2[x * 3 + 5] - row2[x * 3 - 1]) * 2 + row2[x * 3 + n2 + 5] - row2[x * 3 + n2 - 1] + row2[x * 3 + s2 + 5] - row2[x * 3 + s2 - 1]];

			prow1[x + width * 3] = row1[x * 3];
			prow1[x + width * 4] = row1[x * 3 + 1];
			prow1[x + width * 5] = row1[x * 3 + 2];

			prow2[width - 1 - x + width * 3] = row2[x * 3];
			prow2[width - 1 - x + width * 4] = row2[x * 3 + 1];
			prow2[width - 1 - x + width * 5] = row2[x * 3 + 2];
		}
	}

	memset(cost, 0, width1*D*sizeof(cost[0]));

	buffer -= minX2;
	cost -= minX1*D + minD; // simplify the cost indices inside the loop

#if 1
	for (c = 0; c < cn * 2; c++, prow1 += width, prow2 += width)
	{
		int diff_scale = c < cn ? 0 : 2;

		// precompute
		//   v0 = min(row2[x-1/2], row2[x], row2[x+1/2]) and
		//   v1 = max(row2[x-1/2], row2[x], row2[x+1/2]) and

		//做右影像
		for (x = minX2; x < maxX2; x++)
		{
			int v = prow2[x];
			int vl = x > 0 ? (v + prow2[x - 1]) / 2 : v;
			int vr = x < width - 1 ? (v + prow2[x + 1]) / 2 : v;
			int v0 = std::min(vl, vr); v0 = std::min(v0, v);
			int v1 = std::max(vl, vr); v1 = std::max(v1, v);
			buffer[x] = (PixType)v0;
			buffer[x + width2] = (PixType)v1;
		}

		//做左影像
		for (x = minX1; x < maxX1; x++)
		{
			int u = prow1[x];
			int ul = x > 0 ? (u + prow1[x - 1]) / 2 : u;
			int ur = x < width - 1 ? (u + prow1[x + 1]) / 2 : u;
			int u0 = std::min(ul, ur); u0 = std::min(u0, u);
			int u1 = std::max(ul, ur); u1 = std::max(u1, u);

#if CV_SIMD128
			v_uint8x16 _u = v_setall_u8((uchar)u), _u0 = v_setall_u8((uchar)u0);
			v_uint8x16 _u1 = v_setall_u8((uchar)u1);

			for (int d = minD; d < maxD; d += 16)
			{
				v_uint8x16 _v = v_load(prow2 + width - x - 1 + d);
				v_uint8x16 _v0 = v_load(buffer + width - x - 1 + d);
				v_uint8x16 _v1 = v_load(buffer + width - x - 1 + d + width2);
				v_uint8x16 c0 = v_max(_u - _v1, _v0 - _u);
				v_uint8x16 c1 = v_max(_v - _u1, _u0 - _v);
				v_uint8x16 diff = v_min(c0, c1);

				v_int16x8 _c0 = v_load_aligned(cost + x*D + d);
				v_int16x8 _c1 = v_load_aligned(cost + x*D + d + 8);

				v_uint16x8 diff1, diff2;
				v_expand(diff, diff1, diff2);
				v_store_aligned(cost + x*D + d, _c0 + v_reinterpret_as_s16(diff1 >> diff_scale));
				v_store_aligned(cost + x*D + d + 8, _c1 + v_reinterpret_as_s16(diff2 >> diff_scale));
			}
#else
			for (int d = minD; d < maxD; d++)
			{
				int v = prow2[width - x - 1 + d];
				int v0 = buffer[width - x - 1 + d];
				int v1 = buffer[width - x - 1 + d + width2];
				int c0 = std::max(0, u - v1); c0 = std::max(c0, v0 - u);
				int c1 = std::max(0, v - u1); c1 = std::max(c1, u0 - v);

				cost[x*D + d] = (CostType)(cost[x*D + d] + (std::min(c0, c1) >> diff_scale));
			}
#endif
		}
	}
#else
	for (c = 0; c < cn * 2; c++, prow1 += width, prow2 += width)
	{
		for (x = minX1; x < maxX1; x++)
		{
			int u = prow1[x];
#if CV_SSE2
			if (useSIMD)
			{
				__m128i _u = _mm_set1_epi8(u), z = _mm_setzero_si128();

				for (int d = minD; d < maxD; d += 16)
				{
					__m128i _v = _mm_loadu_si128((const __m128i*)(prow2 + width - 1 - x + d));
					__m128i diff = _mm_adds_epu8(_mm_subs_epu8(_u, _v), _mm_subs_epu8(_v, _u));
					__m128i c0 = _mm_load_si128((__m128i*)(cost + x*D + d));
					__m128i c1 = _mm_load_si128((__m128i*)(cost + x*D + d + 8));

					_mm_store_si128((__m128i*)(cost + x*D + d), _mm_adds_epi16(c0, _mm_unpacklo_epi8(diff, z)));
					_mm_store_si128((__m128i*)(cost + x*D + d + 8), _mm_adds_epi16(c1, _mm_unpackhi_epi8(diff, z)));
				}
			}
			else
#endif
			{
				for (int d = minD; d < maxD; d++)
				{
					int v = prow2[width - 1 - x + d];
					cost[x*D + d] = (CostType)(cost[x*D + d] + (CostType)std::abs(u - v));
				}
			}
		}
	}
#endif
}

static void computeDisparitySGBM(const Mat& img1, const Mat& img2,
	Mat& disp1, const StereoSGBMParams& params,
	Mat& buffer)
{
#if CV_SSE2
	static const uchar LSBTab[] =
	{
		0, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
		5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
		6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
		5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
		7, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
		5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
		6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0,
		5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0
	};

	volatile bool useSIMD = checkHardwareSupport(CV_CPU_SSE2);
#endif

	const int ALIGN = 16;
	//const int DISP_SHIFT = StereoMatcher::DISP_SHIFT;
	const int DISP_SCALE = (1 << DISP_SHIFT);
	const CostType MAX_COST = SHRT_MAX;

	int minD = params.minDisparity, maxD = minD + params.numDisparities;
	Size SADWindowSize;
	SADWindowSize.width = SADWindowSize.height = params.SADWindowSize > 0 ? params.SADWindowSize : 5;
	//int ftzero = std::max(params.preFilterCap, 15) | 1;
	int ftzero = 63;
	int uniquenessRatio = params.uniquenessRatio >= 0 ? params.uniquenessRatio : 10;
	int disp12MaxDiff = params.disp12MaxDiff > 0 ? params.disp12MaxDiff : 1;
	//int P1 = params.P1 > 0 ? params.P1 : 2, P2 = std::max(params.P2 > 0 ? (int)params.P2 : 5, (int)(P1 + 1));
	int P1 = params.P1 > 0 ? params.P1 : 2, P2 = params.P2;
	P2 = std::max(P2 > 0 ? P2 : 5, P1 + 1);
	int k, width = disp1.cols, height = disp1.rows;
	int minX1 = std::max(maxD, 0), maxX1 = width + std::min(minD, 0);
	int D = maxD - minD, width1 = maxX1 - minX1;
	int INVALID_DISP = minD - 1, INVALID_DISP_SCALED = INVALID_DISP*DISP_SCALE;
	int SW2 = SADWindowSize.width / 2, SH2 = SADWindowSize.height / 2;
	//bool fullDP = params.mode == StereoSGBM::MODE_HH;
	bool fullDP = params.mode == MODE_HH;
	int npasses = fullDP ? 2 : 1;
	const int TAB_OFS = 256 * 4, TAB_SIZE = 256 + TAB_OFS * 2;
	PixType clipTab[TAB_SIZE];

	/************************************************************
		OpenCV 將 gradient 的值限制在 0 ~ ftzero*2 之間
		因此如果 gradient score 的值為 0 => 則 tab[0] = ftzero
		若 gradient 極大,值為255 => tab[255] 最多只到 ftzero * 2
		若 gradient 極小,值為-255 => tab[-255] 最小只到 0
	/************************************************************/
	for (k = 0; k < TAB_SIZE; k++)
		clipTab[k] = (PixType)(std::min(std::max(k - TAB_OFS, -ftzero), ftzero) + ftzero);

	if (minX1 >= maxX1)
	{
		disp1 = Scalar::all(INVALID_DISP_SCALED);
		return;
	}

	CV_Assert(D % 16 == 0);

	// NR - the number of directions. the loop on x below that computes Lr assumes that NR == 8.
	// if you change NR, please, modify the loop as well.
	int D2 = D + 16, NRD2 = NR2*D2;

	// the number of L_r(.,.) and min_k L_r(.,.) lines in the buffer:
	// for 8-way dynamic programming we need the current row and
	// the previous row, i.e. 2 rows in total
	const int NLR = 2;
	const int LrBorder = NLR - 1;

	// for each possible stereo match (img1(x,y) <=> img2(x-d,y))
	// we keep pixel difference cost (C) and the summary cost over NR directions (S).
	// we also keep all the partial costs for the previous line L_r(x,d) and also min_k L_r(x, k)
	size_t costBufSize = width1*D;
	size_t CSBufSize = costBufSize*(fullDP ? height : 1);
	size_t minLrSize = (width1 + LrBorder * 2)*NR2, LrSize = minLrSize*D2;
	int hsumBufNRows = SH2 * 2 + 2;
	size_t totalBufSize = (LrSize + minLrSize)*NLR*sizeof(CostType)+ // minLr[] and Lr[]
		costBufSize*(hsumBufNRows + 1)*sizeof(CostType)+ // hsumBuf, pixdiff
		CSBufSize * 2 * sizeof(CostType)+ // C, S
		width * 16 * img1.channels()*sizeof(PixType)+ // temp buffer for computing per-pixel cost
		width*(sizeof(CostType)+sizeof(DispType)) + 1024; // disp2cost + disp2

	if (buffer.empty() || !buffer.isContinuous() ||
		buffer.cols*buffer.rows*buffer.elemSize() < totalBufSize)
		buffer.create(1, (int)totalBufSize, CV_8U);

	// summary cost over different (nDirs) directions
	CostType* Cbuf = (CostType*)alignPtr(buffer.ptr(), ALIGN);
	CostType* Sbuf = Cbuf + CSBufSize;
	CostType* hsumBuf = Sbuf + CSBufSize;
	CostType* pixDiff = hsumBuf + costBufSize*hsumBufNRows;

	CostType* disp2cost = pixDiff + costBufSize + (LrSize + minLrSize)*NLR;
	DispType* disp2ptr = (DispType*)(disp2cost + width);
	PixType* tempBuf = (PixType*)(disp2ptr + width);

	// add P2 to every C(x,y). it saves a few operations in the inner loops
	for (k = 0; k < width1*D; k++)
		Cbuf[k] = (CostType)P2;

	/************* note *************/
	/*
	設 params.preFilterCap = 63 則
	int ftzero = std::max(params.preFilterCap, 15) | 1
	= 63 (params.preFilterCap = 62或63 ftzero的結果都是63)
	TAB_OFS = 1024
	TAB_SIZE = 2304
	clipTab[0] ~ clipTab[961] = 0
	clipTab[962] = 1
	clipTab[963] = 2
	clipTab[964] = 3
	......
	clipTab[1087] = 126
	clipTab[1088] = 126
	......
	clipTab[2303] = 126

	(目前合理猜測stereoSGBM 5 way aggregation是指: 左 左上 上 右上 右,
	因為每個 row 在掃完就決定了 bestD, 那能夠利用的就只有 "該row" 跟 "之前算過的資訊" 了)
	NR - the number of direction (NR從頭到尾沒被用過)(如果是 5 way NR = 5)
	NR2 似乎才是實際用的值 (預設NR = 16, NR2 = NR/2)

	dx, dy => x,y loop 時累進的方向(pass1時:+1(遞增) or pass2:-1(從後面算回來))

	maxD = minD + params.numDisparities;     (128 in my experiment, minD = 0)

	width = disp1.cols, height = disp1.rows;
	int minX1 = std::max(maxD, 0)            (minX1 = 128)
	maxX1 = width + std::min(minD, 0);   (maxX1 = 612)

	D  = maxD - minD						 (D = 128)
	D2 = D + 16								 (D2 = D + 16)
	width1 = maxX1 - minX1;              (width1 = 484)

	DISP_SCALE = 16(應該是)

	int SW2 = SADWindowSize.width/2
	SH2 = SADWindowSize.height/2;

	// the number of L_r(.,.) and min_k L_r(.,.) lines in the buffer:
	// for 8-way dynamic programming we need the current row and
	// the previous row, i.e. 2 rows in total
	// 意思就是 aggregation 需要 NLR = 2 個 row 去做 aggregation (這兩個row包含前一個row, 與現在的row, 共兩個row)
	const int NLR = 2;
	const int LrBorder = NLR - 1;

	// for each possible stereo match (img1(x,y) <=> img2(x-d,y))
	// we keep pixel difference cost (C) and the summary cost over NR directions (S).
	// we also keep all the partial costs for the previous line L_r(x,d) and also min_k L_r(x, k)
	*/

	for (int pass = 1; pass <= npasses; pass++)
	{
		int x1, y1, x2, y2, dx, dy;

		if (pass == 1)
		{
			y1 = 0; y2 = height; dy = 1;
			x1 = 0; x2 = width1; dx = 1;
		}
		else
		{
			y1 = height - 1; y2 = -1; dy = -1;
			x1 = width1 - 1; x2 = -1; dx = -1;
		}

		CostType *Lr[NLR] = { 0 }, *minLr[NLR] = { 0 };

		for (k = 0; k < NLR; k++)
		{
			// shift Lr[k] and minLr[k] pointers, because we allocated them with the borders,
			// and will occasionally use negative indices with the arrays
			// we need to shift Lr[k] pointers by 1, to give the space for d=-1.
			// however, then the alignment will be imperfect, i.e. bad for SSE,
			// thus we shift the pointers by 8 (8*sizeof(short) == 16 - ideal alignment)
			Lr[k] = pixDiff + costBufSize + LrSize*k + NRD2*LrBorder + 8;
			memset(Lr[k] - LrBorder*NRD2 - 8, 0, LrSize*sizeof(CostType));
			minLr[k] = pixDiff + costBufSize + LrSize*NLR + minLrSize*k + NR2*LrBorder;
			memset(minLr[k] - LrBorder*NR2, 0, minLrSize*sizeof(CostType));
		}

		for (int y = y1; y != y2; y += dy)
		{
			int x, d;
			DispType* disp1ptr = disp1.ptr<DispType>(y);
			CostType* C = Cbuf + (!fullDP ? 0 : y*costBufSize);
			CostType* S = Sbuf + (!fullDP ? 0 : y*costBufSize);

			if (pass == 1) // compute C on the first pass, and reuse it on the second pass, if any.
			{
				int dy1 = y == 0 ? 0 : y + SH2, dy2 = y == 0 ? SH2 : dy1;
				/*
				只有 pass == 1時才會做
				在y == 0時
				dy1 = 0
				dy2 = SH2 = 2
				在y == rest
				dy1 = y + SH2 = y + 2
				dy2 = dy1     = y + 2
				
				演算法會預先做下兩個row的結果
				*/

				for (k = dy1; k <= dy2; k++)
				{
					/***************************************************
						在此迴圈中, 會計算 SAD 的 cost (Block Matching), 步驟如下
						(此處假設 SADWindowSize = 5, SH2 = 2, SW2 = 2)
						1. 計算第一個 row 中每個pixel 的 pixDiff (calcPixelCostBT)
						2. 將水平的 pixDiff 累加到 hsum 中
						(pixDiff * scale => 若上面會左右邊沒有pixel, 則用第一個來補, scale就是控制的權重, 因此對第一個來說scale都是 3 (SW2+1或SH2+1))
						(hsum[x] = (x-2) ~ (x+2) 的 pixDiff總和 )
						3. 將 hsum 的值加到 C 中 (透過C來把Cost傳遞一下面的row)
						C[x] = hsumAdd[x]
						(現在的 C 在經過後面的 cost aggregation 後, 就會成為下次的 Cprev)

						4. 計算第二個 row 中每個pixel 的 pixDiff
						5. 將水平的 pixDiff 累加到 hsum 中
						6. C[x] = Cprev[x] + hsumAdd - hsumSub
						(hsumAdd 是 y+2 那個 row 的水平累加結果)
						(hsumSub 是 y-3 那個 row 的水平累加結果)
						(加入新的, 減掉最後的, 加快運算速度, 重複利用, 不用每次都重算)
						之後就是不斷重複這個過程

						特別需要注意的是 此程式的流程每個 row 都會通過 semi-global cost aggregation
						因此 Cprev 都是經過 cost aggregation的(這其實也是一條由上往下的cost aggregation路徑)
					/***************************************************/
					//除了第一個row以外, 之後做的都是第y+2個row
					//costBufSize  是 一個 row 的cost buffer size大小 = width1 * D bytes
					//hsumBufNRows 是 水平kernel size + 1 = 5 + 1 
					//			   +1 的話 mod 的結果才會控制在 0 ~ 5
					CostType* hsumAdd = hsumBuf + (std::min(k, height - 1) % hsumBufNRows)*costBufSize;

					if (k < height)
					{
						calcPixelCostBT(img1, img2, k, minD, maxD, pixDiff, tempBuf, clipTab, TAB_OFS, ftzero);

						memset(hsumAdd, 0, D*sizeof(CostType));
						/*
								** 累加 kernel中 水平的pixDif **
							邊緣處理 x=0
								hsumAdd = (x=0) + (x=1) + (x=2) 的 pixDiff (各 disparity 獨立) 
								scale 在 x==0 時 為 SW2 + 1, 因為 kernel 左邊沒有像素(少 SW2 個像素), 所以就直接拿中心點的值來充數
							中間     x=1 ~ width1
								hsumAdd[x] = (x-2) + (x-1) + x + (x+1) + (x+2) 的 pixDiff 
								
						*/
						for (x = 0; x <= SW2*D; x += D)
						{
							int scale = x == 0 ? SW2 + 1 : 1;
							for (d = 0; d < D; d++){
								hsumAdd[d] = (CostType)(hsumAdd[d] + pixDiff[x + d] * scale);
								//printf("x:%d d:%d; hsumAdd[d] = %d\n", x, d, hsumAdd[d]);
							}
						}

						if (y > 0)
						{
							const CostType* hsumSub = hsumBuf + (std::max(y - SH2 - 1, 0) % hsumBufNRows)*costBufSize;
							const CostType* Cprev = !fullDP || y == 0 ? C : C - costBufSize;

							for (x = D; x < width1*D; x += D)
							{
								const CostType* pixAdd = pixDiff + std::min(x + SW2*D, (width1 - 1)*D);
								const CostType* pixSub = pixDiff + std::max(x - (SW2 + 1)*D, 0);

#if CV_SSE2
								if (useSIMD)
								{
									for (d = 0; d < D; d += 8)
									{
										__m128i hv = _mm_load_si128((const __m128i*)(hsumAdd + x - D + d));
										__m128i Cx = _mm_load_si128((__m128i*)(Cprev + x + d));
										hv = _mm_adds_epi16(_mm_subs_epi16(hv,
											_mm_load_si128((const __m128i*)(pixSub + d))),
											_mm_load_si128((const __m128i*)(pixAdd + d)));
										Cx = _mm_adds_epi16(_mm_subs_epi16(Cx,
											_mm_load_si128((const __m128i*)(hsumSub + x + d))),
											hv);
										_mm_store_si128((__m128i*)(hsumAdd + x + d), hv);
										_mm_store_si128((__m128i*)(C + x + d), Cx);
									}
								}
								else
#endif
								{
									for (d = 0; d < D; d++)
									{
										//這次的 horizontal kernel sum = 上次的hsum + 下一個pixDiff - 上次的最後一個pixDiff
										int hv = hsumAdd[x + d] = (CostType)(hsumAdd[x - D + d] + pixAdd[d] - pixSub[d]);
										//這次的Cost = 
										C[x + d] = (CostType)(Cprev[x + d] + hv - hsumSub[x + d]);
									}
								}
							}
						}
						else
						{
							for (x = D; x < width1*D; x += D)
							{
								const CostType* pixAdd = pixDiff + std::min(x + SW2*D, (width1 - 1)*D);
								const CostType* pixSub = pixDiff + std::max(x - (SW2 + 1)*D, 0);

								for (d = 0; d < D; d++)
									hsumAdd[x + d] = (CostType)(hsumAdd[x - D + d] + pixAdd[d] - pixSub[d]);
							}
						}
					}

					/*
						第一列的 C[x] = hsumAdd[x] * scale (因為沒有上兩個row, 所以用第一個row的值補)
					*/
					if (y == 0)
					{
						int scale = k == 0 ? SH2 + 1 : 1;
						for (x = 0; x < width1*D; x++)
							C[x] = (CostType)(C[x] + hsumAdd[x] * scale);
					}
				}

				// also, clear the S buffer
				for (k = 0; k < width1*D; k++)
					S[k] = 0;
			}

			// clear the left and the right borders
			memset(Lr[0] - NRD2*LrBorder - 8, 0, NRD2*LrBorder*sizeof(CostType));
			memset(Lr[0] + width1*NRD2 - 8, 0, NRD2*LrBorder*sizeof(CostType));
			memset(minLr[0] - NR2*LrBorder, 0, NR2*LrBorder*sizeof(CostType));
			memset(minLr[0] + width1*NR2, 0, NR2*LrBorder*sizeof(CostType));

			/*
			[formula 13 in the paper]
			compute L_r(p, d) = C(p, d) +
			min(L_r(p-r, d),
			L_r(p-r, d-1) + P1,
			L_r(p-r, d+1) + P1,
			min_k L_r(p-r, k) + P2) - min_k L_r(p-r, k)
			where p = (x,y), r is one of the directions.
			we process all the directions at once:
			0: r=(-dx, 0)
			1: r=(-1, -dy)
			2: r=(0, -dy)
			3: r=(1, -dy)
			4: r=(-2, -dy)
			5: r=(-1, -dy*2)
			6: r=(1, -dy*2)
			7: r=(2, -dy)
			*/
			for (x = x1; x != x2; x += dx)
			{
				int xm = x*NR2, xd = xm*D2;

				int delta0 = minLr[0][xm - dx*NR2] + P2, delta1 = minLr[1][xm - NR2 + 1] + P2;
				int delta2 = minLr[1][xm + 2] + P2, delta3 = minLr[1][xm + NR2 + 3] + P2;

				CostType* Lr_p0 = Lr[0] + xd - dx*NRD2;
				CostType* Lr_p1 = Lr[1] + xd - NRD2 + D2;
				CostType* Lr_p2 = Lr[1] + xd + D2 * 2;
				CostType* Lr_p3 = Lr[1] + xd + NRD2 + D2 * 3;

				Lr_p0[-1] = Lr_p0[D] = Lr_p1[-1] = Lr_p1[D] =
					Lr_p2[-1] = Lr_p2[D] = Lr_p3[-1] = Lr_p3[D] = MAX_COST;

				CostType* Lr_p = Lr[0] + xd;
				const CostType* Cp = C + x*D;
				CostType* Sp = S + x*D;

#if CV_SSE2
				if (useSIMD)
				{
					__m128i _P1 = _mm_set1_epi16((short)P1);

					__m128i _delta0 = _mm_set1_epi16((short)delta0);
					__m128i _delta1 = _mm_set1_epi16((short)delta1);
					__m128i _delta2 = _mm_set1_epi16((short)delta2);
					__m128i _delta3 = _mm_set1_epi16((short)delta3);
					__m128i _minL0 = _mm_set1_epi16((short)MAX_COST);

					for (d = 0; d < D; d += 8)
					{
						__m128i Cpd = _mm_load_si128((const __m128i*)(Cp + d));
						__m128i L0, L1, L2, L3;

						L0 = _mm_load_si128((const __m128i*)(Lr_p0 + d));
						L1 = _mm_load_si128((const __m128i*)(Lr_p1 + d));
						L2 = _mm_load_si128((const __m128i*)(Lr_p2 + d));
						L3 = _mm_load_si128((const __m128i*)(Lr_p3 + d));

						L0 = _mm_min_epi16(L0, _mm_adds_epi16(_mm_loadu_si128((const __m128i*)(Lr_p0 + d - 1)), _P1));
						L0 = _mm_min_epi16(L0, _mm_adds_epi16(_mm_loadu_si128((const __m128i*)(Lr_p0 + d + 1)), _P1));

						L1 = _mm_min_epi16(L1, _mm_adds_epi16(_mm_loadu_si128((const __m128i*)(Lr_p1 + d - 1)), _P1));
						L1 = _mm_min_epi16(L1, _mm_adds_epi16(_mm_loadu_si128((const __m128i*)(Lr_p1 + d + 1)), _P1));

						L2 = _mm_min_epi16(L2, _mm_adds_epi16(_mm_loadu_si128((const __m128i*)(Lr_p2 + d - 1)), _P1));
						L2 = _mm_min_epi16(L2, _mm_adds_epi16(_mm_loadu_si128((const __m128i*)(Lr_p2 + d + 1)), _P1));

						L3 = _mm_min_epi16(L3, _mm_adds_epi16(_mm_loadu_si128((const __m128i*)(Lr_p3 + d - 1)), _P1));
						L3 = _mm_min_epi16(L3, _mm_adds_epi16(_mm_loadu_si128((const __m128i*)(Lr_p3 + d + 1)), _P1));

						L0 = _mm_min_epi16(L0, _delta0);
						L0 = _mm_adds_epi16(_mm_subs_epi16(L0, _delta0), Cpd);

						L1 = _mm_min_epi16(L1, _delta1);
						L1 = _mm_adds_epi16(_mm_subs_epi16(L1, _delta1), Cpd);

						L2 = _mm_min_epi16(L2, _delta2);
						L2 = _mm_adds_epi16(_mm_subs_epi16(L2, _delta2), Cpd);

						L3 = _mm_min_epi16(L3, _delta3);
						L3 = _mm_adds_epi16(_mm_subs_epi16(L3, _delta3), Cpd);

						_mm_store_si128((__m128i*)(Lr_p + d), L0);
						_mm_store_si128((__m128i*)(Lr_p + d + D2), L1);
						_mm_store_si128((__m128i*)(Lr_p + d + D2 * 2), L2);
						_mm_store_si128((__m128i*)(Lr_p + d + D2 * 3), L3);

						__m128i t0 = _mm_min_epi16(_mm_unpacklo_epi16(L0, L2), _mm_unpackhi_epi16(L0, L2));
						__m128i t1 = _mm_min_epi16(_mm_unpacklo_epi16(L1, L3), _mm_unpackhi_epi16(L1, L3));
						t0 = _mm_min_epi16(_mm_unpacklo_epi16(t0, t1), _mm_unpackhi_epi16(t0, t1));
						_minL0 = _mm_min_epi16(_minL0, t0);

						__m128i Sval = _mm_load_si128((const __m128i*)(Sp + d));

						L0 = _mm_adds_epi16(L0, L1);
						L2 = _mm_adds_epi16(L2, L3);
						Sval = _mm_adds_epi16(Sval, L0);
						Sval = _mm_adds_epi16(Sval, L2);

						_mm_store_si128((__m128i*)(Sp + d), Sval);
					}

					_minL0 = _mm_min_epi16(_minL0, _mm_srli_si128(_minL0, 8));
					_mm_storel_epi64((__m128i*)&minLr[0][xm], _minL0);
				}
				else
#endif
				{
					int minL0 = MAX_COST, minL1 = MAX_COST, minL2 = MAX_COST, minL3 = MAX_COST;

					for (d = 0; d < D; d++)
					{
						int Cpd = Cp[d], L0, L1, L2, L3;

						L0 = Cpd + std::min((int)Lr_p0[d], std::min(Lr_p0[d - 1] + P1, std::min(Lr_p0[d + 1] + P1, delta0))) - delta0;
						L1 = Cpd + std::min((int)Lr_p1[d], std::min(Lr_p1[d - 1] + P1, std::min(Lr_p1[d + 1] + P1, delta1))) - delta1;
						L2 = Cpd + std::min((int)Lr_p2[d], std::min(Lr_p2[d - 1] + P1, std::min(Lr_p2[d + 1] + P1, delta2))) - delta2;
						L3 = Cpd + std::min((int)Lr_p3[d], std::min(Lr_p3[d - 1] + P1, std::min(Lr_p3[d + 1] + P1, delta3))) - delta3;

						Lr_p[d] = (CostType)L0;
						minL0 = std::min(minL0, L0);

						Lr_p[d + D2] = (CostType)L1;
						minL1 = std::min(minL1, L1);

						Lr_p[d + D2 * 2] = (CostType)L2;
						minL2 = std::min(minL2, L2);

						Lr_p[d + D2 * 3] = (CostType)L3;
						minL3 = std::min(minL3, L3);

						Sp[d] = saturate_cast<CostType>(Sp[d] + L0 + L1 + L2 + L3);
					}
					minLr[0][xm] = (CostType)minL0;
					minLr[0][xm + 1] = (CostType)minL1;
					minLr[0][xm + 2] = (CostType)minL2;
					minLr[0][xm + 3] = (CostType)minL3;
				}
			}

			if (pass == npasses)
			{
				for (x = 0; x < width; x++)
				{
					disp1ptr[x] = disp2ptr[x] = (DispType)INVALID_DISP_SCALED;
					disp2cost[x] = MAX_COST;
				}

				for (x = width1 - 1; x >= 0; x--)
				{
					CostType* Sp = S + x*D;
					int minS = MAX_COST, bestDisp = -1;

					if (npasses == 1)
					{
						int xm = x*NR2, xd = xm*D2;

						int minL0 = MAX_COST;
						int delta0 = minLr[0][xm + NR2] + P2;
						CostType* Lr_p0 = Lr[0] + xd + NRD2;
						Lr_p0[-1] = Lr_p0[D] = MAX_COST;
						CostType* Lr_p = Lr[0] + xd;

						const CostType* Cp = C + x*D;

#if CV_SSE2
						if (useSIMD)
						{
							__m128i _P1 = _mm_set1_epi16((short)P1);
							__m128i _delta0 = _mm_set1_epi16((short)delta0);

							__m128i _minL0 = _mm_set1_epi16((short)minL0);
							__m128i _minS = _mm_set1_epi16(MAX_COST), _bestDisp = _mm_set1_epi16(-1);
							__m128i _d8 = _mm_setr_epi16(0, 1, 2, 3, 4, 5, 6, 7), _8 = _mm_set1_epi16(8);

							for (d = 0; d < D; d += 8)
							{
								__m128i Cpd = _mm_load_si128((const __m128i*)(Cp + d)), L0;

								L0 = _mm_load_si128((const __m128i*)(Lr_p0 + d));
								L0 = _mm_min_epi16(L0, _mm_adds_epi16(_mm_loadu_si128((const __m128i*)(Lr_p0 + d - 1)), _P1));
								L0 = _mm_min_epi16(L0, _mm_adds_epi16(_mm_loadu_si128((const __m128i*)(Lr_p0 + d + 1)), _P1));
								L0 = _mm_min_epi16(L0, _delta0);
								L0 = _mm_adds_epi16(_mm_subs_epi16(L0, _delta0), Cpd);

								_mm_store_si128((__m128i*)(Lr_p + d), L0);
								_minL0 = _mm_min_epi16(_minL0, L0);
								L0 = _mm_adds_epi16(L0, *(__m128i*)(Sp + d));
								_mm_store_si128((__m128i*)(Sp + d), L0);

								__m128i mask = _mm_cmpgt_epi16(_minS, L0);
								_minS = _mm_min_epi16(_minS, L0);
								_bestDisp = _mm_xor_si128(_bestDisp, _mm_and_si128(_mm_xor_si128(_bestDisp, _d8), mask));
								_d8 = _mm_adds_epi16(_d8, _8);
							}

							short CV_DECL_ALIGNED(16) bestDispBuf[8];
							_mm_store_si128((__m128i*)bestDispBuf, _bestDisp);

							_minL0 = _mm_min_epi16(_minL0, _mm_srli_si128(_minL0, 8));
							_minL0 = _mm_min_epi16(_minL0, _mm_srli_si128(_minL0, 4));
							_minL0 = _mm_min_epi16(_minL0, _mm_srli_si128(_minL0, 2));

							__m128i qS = _mm_min_epi16(_minS, _mm_srli_si128(_minS, 8));
							qS = _mm_min_epi16(qS, _mm_srli_si128(qS, 4));
							qS = _mm_min_epi16(qS, _mm_srli_si128(qS, 2));

							minLr[0][xm] = (CostType)_mm_cvtsi128_si32(_minL0);
							minS = (CostType)_mm_cvtsi128_si32(qS);

							qS = _mm_shuffle_epi32(_mm_unpacklo_epi16(qS, qS), 0);
							qS = _mm_cmpeq_epi16(_minS, qS);
							int idx = _mm_movemask_epi8(_mm_packs_epi16(qS, qS)) & 255;

							bestDisp = bestDispBuf[LSBTab[idx]];
						}
						else
#endif
						{
							for (d = 0; d < D; d++)
							{
								int L0 = Cp[d] + std::min((int)Lr_p0[d], std::min(Lr_p0[d - 1] + P1, std::min(Lr_p0[d + 1] + P1, delta0))) - delta0;

								Lr_p[d] = (CostType)L0;
								minL0 = std::min(minL0, L0);

								int Sval = Sp[d] = saturate_cast<CostType>(Sp[d] + L0);
								if (Sval < minS)
								{
									minS = Sval;
									bestDisp = d;
								}
							}
							minLr[0][xm] = (CostType)minL0;
						}
					}
					else
					{
						for (d = 0; d < D; d++)
						{
							int Sval = Sp[d];
							if (Sval < minS)
							{
								minS = Sval;
								bestDisp = d;
							}
						}
					}

					for (d = 0; d < D; d++)
					{
						if (Sp[d] * (100 - uniquenessRatio) < minS * 100 && std::abs(bestDisp - d) > 1)
							break;
					}
					if (d < D)
						continue;
					d = bestDisp;
					int _x2 = x + minX1 - d - minD;
					if (disp2cost[_x2] > minS)
					{
						disp2cost[_x2] = (CostType)minS;
						disp2ptr[_x2] = (DispType)(d + minD);
					}

					if (0 < d && d < D - 1)
					{
						// do subpixel quadratic interpolation:
						//   fit parabola into (x1=d-1, y1=Sp[d-1]), (x2=d, y2=Sp[d]), (x3=d+1, y3=Sp[d+1])
						//   then find minimum of the parabola.
						/*******************************************************************
							(Subpixel)
							用 Cost 當作 y 軸座標 (aggregation後的)
							用 d    當作 x 軸座標
							用 (d-1) (d) (d+1) 三個點,使用 parabola fitting,
							找出該 拋物線真正的 minimum 位置 作為最後 subpixel 的結果
							可參考 Sub-pixel Estimation of Local Extrema (不過公式有些不同, 目前找不到OpenCV這個版本的算法從哪來的)
						/*******************************************************************/
						int denom2 = std::max(Sp[d - 1] + Sp[d + 1] - 2 * Sp[d], 1);
						d = d*DISP_SCALE + ((Sp[d - 1] - Sp[d + 1])*DISP_SCALE + denom2) / (denom2 * 2);
					}
					else
						d *= DISP_SCALE;
					disp1ptr[x + minX1] = (DispType)(d + minD*DISP_SCALE);
				}

				for (x = minX1; x < maxX1; x++)
				{
					// we round the computed disparity both towards -inf and +inf and check
					// if either of the corresponding disparities in disp2 is consistent.
					// This is to give the computed disparity a chance to look valid if it is.
					int d1 = disp1ptr[x];
					if (d1 == INVALID_DISP_SCALED)
						continue;
					int _d = d1 >> DISP_SHIFT;
					int d_ = (d1 + DISP_SCALE - 1) >> DISP_SHIFT;
					int _x = x - _d, x_ = x - d_;
					if (0 <= _x && _x < width && disp2ptr[_x] >= minD && std::abs(disp2ptr[_x] - _d) > disp12MaxDiff &&
						0 <= x_ && x_ < width && disp2ptr[x_] >= minD && std::abs(disp2ptr[x_] - d_) > disp12MaxDiff)
						disp1ptr[x] = (DispType)INVALID_DISP_SCALED;
				}
			}

			// now shift the cyclic buffers
			std::swap(Lr[0], Lr[1]);
			std::swap(minLr[0], minLr[1]);
		}
	}
}

void compute(InputArray leftarr, InputArray rightarr, OutputArray disparr, StereoSGBMParams params)
{
	Mat buffer;
	const int DISP_SCALE = (1 << DISP_SHIFT);

	Mat left = leftarr.getMat(), right = rightarr.getMat();
	//CV_ASSERT(left.size() == right.size() && left.type() == right.type() &&
	//		  left.depth() == CV_8U);
	disparr.create(left.size(), CV_16S);
	Mat disp = disparr.getMat();

	computeDisparitySGBM(left, right, disp, params, buffer);
	medianBlur(disp, disp, 3);

	if (params.speckleWindowSize > 0)
		filterSpeckles(disp, (params.minDisparity - 1) * DISP_SCALE, params.speckleWindowSize, DISP_SCALE*params.speckleRange, buffer);
}