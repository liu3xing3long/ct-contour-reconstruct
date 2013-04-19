/* GIMP - The GNU Image Manipulation Program
* Copyright (C) 1995 Spencer Kimball and Peter Mattis
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "stdafx.h"
#include "image_enhance.h"




template <typename T> T CLAMP(const T& value, const T& low, const T& high) 
{
	return value < low ? low : (value > high ? high : value); 
}


void image_Enhance( float* fImage, int width, int height )
{

}


#include <string.h>
#include <math.h>

#define MAX_RETINEX_SCALES    8
#define MIN_GAUSSIAN_SCALE   16
#define MAX_GAUSSIAN_SCALE  250
#define SCALE_WIDTH         150
#define ENTRY_WIDTH           4

#define g_try_malloc malloc
#define g_free	free
#define g_warning printf

gfloat RetinexScales[MAX_RETINEX_SCALES];

// GimpPlugInInfo PLUG_IN_INFO =
// {
//   NULL,  /* init_proc  */
//   NULL,  /* quit_proc  */
//   query, /* query_proc */
//   run,   /* run_proc   */
// };

// /*
//  * Applies the algorithm
//  */
// void
// retinex (GimpDrawable *drawable,
//          GimpPreview  *preview)
// {
//   gint          x, y, width, height;
//   gint          size, bytes;
//   guchar       *src  = NULL;
//   guchar       *psrc = NULL;
//   GimpPixelRgn  dst_rgn, src_rgn;
// 
//   bytes = drawable->bpp;
// 
//   /*
//    * Get the size of the current image or its selection.
//    */
//   if (preview)
//     {
//       src = gimp_zoom_preview_get_source (GIMP_ZOOM_PREVIEW (preview),
//                                           &width, &height, &bytes);
//     }
//   else
//     {
//       if (! gimp_drawable_mask_intersect (drawable->drawable_id,
//                                           &x, &y, &width, &height))
//         return;
// 
//       /* Allocate memory */
//       size = width * height * bytes;
//       src = g_try_malloc (sizeof (guchar) * size);
// 
//       if (src == NULL)
//         {
//           g_warning ("Failed to allocate memory");
//           return;
//         }
// 
//       memset (src, 0, sizeof (guchar) * size);
// 
//       /* Fill allocated memory with pixel data */
//       gimp_pixel_rgn_init (&src_rgn, drawable,
//                            x, y, width, height,
//                            FALSE, FALSE);
//       gimp_pixel_rgn_get_rect (&src_rgn, src, x, y, width, height);
//     }
// 
//   /*
//     Algorithm for Multi-scale Retinex with color Restoration (MSRCR).
//    */
//   psrc = src;
//   MSRCR (psrc, width, height, bytes, preview != NULL);
// 
//   if (preview)
//     {
//       gimp_preview_draw_buffer (preview, psrc, width * bytes);
//     }
//   else
//     {
//       gimp_pixel_rgn_init (&dst_rgn, drawable,
//                            x, y, width, height,
//                            TRUE, TRUE);
//       gimp_pixel_rgn_set_rect (&dst_rgn, psrc, x, y, width, height);
// 
//       gimp_progress_update (1.0);
// 
//       gimp_drawable_flush (drawable);
//       gimp_drawable_merge_shadow (drawable->drawable_id, TRUE);
//       gimp_drawable_update (drawable->drawable_id, x, y, width, height);
//     }
// 
//   g_free (src);
// }
// 

/*
* calculate scale values for desired distribution.
*/
void
	retinex_scales_distribution(gfloat* scales, gint nscales, gint mode, gint s)
{
	if (nscales == 1)
	{ /* For one filter we choose the median scale */
		scales[0] = (gint) s / 2;
	}
	else if (nscales == 2)
	{ /* For two filters whe choose the median and maximum scale */
		scales[0] = (gint) s / 2;
		scales[1] = (gint) s;
	}
	else
	{
		gfloat size_step = (gfloat) s / (gfloat) nscales;
		gint   i;

		switch(mode)
		{
		case RETINEX_UNIFORM:
			for(i = 0; i < nscales; ++i)
				scales[i] = 2. + (gfloat) i * size_step;
			break;

		case RETINEX_LOW:
			size_step = (gfloat) log(gfloat(s - 2.0f)) / (gfloat) nscales;
			for (i = 0; i < nscales; ++i)
				scales[i] = 2. + pow (10, (i * size_step) / log (10.f));
			break;

		case RETINEX_HIGH:
			size_step = (gfloat) log(gfloat(s - 2.0f)) / (gfloat) nscales;
			for (i = 0; i < nscales; ++i)
				scales[i] = s - pow (10, (i * size_step) / log (10.f));
			break;

		default:
			break;
		}
	}
}

/*
* Calculate the coefficients for the recursive filter algorithm
* Fast Computation of gaussian blurring.
*/
void
	compute_coefs3 (gauss3_coefs *c, gfloat sigma)
{
	/*
	* Papers:  "Recursive Implementation of the gaussian filter.",
	*          Ian T. Young , Lucas J. Van Vliet, Signal Processing 44, Elsevier 1995.
	* formula: 11b       computation of q
	*          8c        computation of b0..b1
	*          10        alpha is normalization constant B
	*/
	gfloat q, q2, q3;

	q = 0;

	if (sigma >= 2.5)
	{
		q = 0.98711 * sigma - 0.96330;
	}
	else if ((sigma >= 0.5) && (sigma < 2.5))
	{
		q = 3.97156 - 4.14554 * (gfloat) sqrt ((double) 1 - 0.26891 * sigma);
	}
	else
	{
		q = 0.1147705018520355224609375;
	}

	q2 = q * q;
	q3 = q * q2;
	c->b[0] = (1.57825+(2.44413*q)+(1.4281 *q2)+(0.422205*q3));
	c->b[1] = (        (2.44413*q)+(2.85619*q2)+(1.26661 *q3));
	c->b[2] = (                   -((1.4281*q2)+(1.26661 *q3)));
	c->b[3] = (                                 (0.422205*q3));
	c->B = 1.0-((c->b[1]+c->b[2]+c->b[3])/c->b[0]);
	c->sigma = sigma;
	c->N = 3;

	/*
	g_printerr ("q %f\n", q);
	g_printerr ("q2 %f\n", q2);
	g_printerr ("q3 %f\n", q3);
	g_printerr ("c->b[0] %f\n", c->b[0]);
	g_printerr ("c->b[1] %f\n", c->b[1]);
	g_printerr ("c->b[2] %f\n", c->b[2]);
	g_printerr ("c->b[3] %f\n", c->b[3]);
	g_printerr ("c->B %f\n", c->B);
	g_printerr ("c->sigma %f\n", c->sigma);
	g_printerr ("c->N %d\n", c->N);
	*/
}

void
	gausssmooth (gfloat *in, gfloat *out, gint size, gint rowstride, gauss3_coefs *c)
{
	/*
	* Papers:  "Recursive Implementation of the gaussian filter.",
	*          Ian T. Young , Lucas J. Van Vliet, Signal Processing 44, Elsevier 1995.
	* formula: 9a        forward filter
	*          9b        backward filter
	*          fig7      algorithm
	*/
	gint i,n, bufsize;
	gfloat *w1,*w2;

	/* forward pass */
	bufsize = size+3;
	size -= 1;
	w1 = (gfloat *) g_try_malloc (bufsize * sizeof (gfloat));
	w2 = (gfloat *) g_try_malloc (bufsize * sizeof (gfloat));
	w1[0] = in[0];
	w1[1] = in[0];
	w1[2] = in[0];
	for ( i = 0 , n=3; i <= size ; i++, n++)
	{
		w1[n] = (gfloat)(c->B*in[i*rowstride] +
			((c->b[1]*w1[n-1] +
			c->b[2]*w1[n-2] +
			c->b[3]*w1[n-3] ) / c->b[0]));
	}

	/* backward pass */
	w2[size+1]= w1[size+3];
	w2[size+2]= w1[size+3];
	w2[size+3]= w1[size+3];
	for (i = size, n = i; i >= 0; i--, n--)
	{
		w2[n]= out[i * rowstride] = (gfloat)(c->B*w1[n] +
			((c->b[1]*w2[n+1] +
			c->b[2]*w2[n+2] +
			c->b[3]*w2[n+3] ) / c->b[0]));
	}

	g_free (w1);
	g_free (w2);
}

/*
* This function is the heart of the algo.
* (a)  Filterings at several scales and sumarize the results.
* (b)  Calculation of the final values.
*/
void
	MSRCR (guchar *src, gint width, gint height, gint bytes, gboolean preview_mode)
{

	gint          scale,row,col;
	gint          i,j;
	gint          size;
	gint          pos;
	gint          channel;
	guchar       *psrc = NULL;            /* backup pointer for src buffer */
	gfloat       *dst  = NULL;            /* float buffer for algorithm */
	gfloat       *pdst = NULL;            /* backup pointer for float buffer */
	gfloat       *in, *out;
	gint          channelsize;            /* Float memory cache for one channel */
	gfloat        weight;
	gauss3_coefs  coef;
	gfloat        mean, var;
	gfloat        mini, range, maxi;
	gfloat        alpha;
	gfloat        gain;
	gfloat        offset;
	gdouble       max_preview = 0.0;

	if (!preview_mode)
	{
		//      gimp_progress_init (_("Retinex: filtering"));
		max_preview = 3 * rvals.nscales;
	}

	/* Allocate all the memory needed for algorithm*/
	size = width * height * bytes;
	dst = (gfloat*)g_try_malloc (size * sizeof (gfloat));
	if (dst == NULL)
	{
		g_warning ("Failed to allocate memory");
		return;
	}
	memset (dst, 0, size * sizeof (gfloat));

	channelsize  = (width * height);
	in  = (gfloat *) g_try_malloc (channelsize * sizeof (gfloat));
	if (in == NULL)
	{
		g_free (dst);
		g_warning ("Failed to allocate memory");
		return; /* do some clever stuff */
	}

	out  = (gfloat *) g_try_malloc (channelsize * sizeof (gfloat));
	if (out == NULL)
	{
		g_free (in);
		g_free (dst);
		g_warning ("Failed to allocate memory");
		return; /* do some clever stuff */
	}


	/*
	Calculate the scales of filtering according to the
	number of filter and their distribution.
	*/

	retinex_scales_distribution (RetinexScales,
		rvals.nscales, rvals.scales_mode, rvals.scale);

	/*
	Filtering according to the various scales.
	Summerize the results of the various filters according to a
	specific weight(here equivalent for all).
	*/
	weight = 1./ (gfloat) rvals.nscales;

	/*
	The recursive filtering algorithm needs different coefficients according
	to the selected scale (~ = standard deviation of Gaussian).
	*/
	pos = 0;
	for (channel = 0; channel < 3; channel++)
	{
		for (i = 0, pos = channel; i < channelsize ; i++, pos += bytes)
		{
			/* 0-255 => 1-256 */
			in[i] = (gfloat)(src[pos] + 1.0);
		}
		for (scale = 0; scale < rvals.nscales; scale++)
		{
			compute_coefs3 (&coef, RetinexScales[scale]);
			/*
			*  Filtering (smoothing) Gaussian recursive.
			*
			*  Filter rows first
			*/
			for (row=0 ;row < height; row++)
			{
				pos =  row * width;
				gausssmooth (in + pos, out + pos, width, 1, &coef);
			}

			memcpy(in,  out, channelsize * sizeof(gfloat));
			memset(out, 0  , channelsize * sizeof(gfloat));

			/*
			*  Filtering (smoothing) Gaussian recursive.
			*
			*  Second columns
			*/
			for (col=0; col < width; col++)
			{
				pos = col;
				gausssmooth(in + pos, out + pos, height, width, &coef);
			}

			/*
			Summarize the filtered values.
			In fact one calculates a ratio between the original values and the filtered values.
			*/
			for (i = 0, pos = channel; i < channelsize; i++, pos += bytes)
			{
				dst[pos] += weight * (log (src[pos] + 1.) - log (out[i]));
			}

			//            if (!preview_mode)
			//              gimp_progress_update ((channel * rvals.nscales + scale) /
			//                                    max_preview);
		}
	}
	g_free(in);
	g_free(out);

	/*
	Final calculation with original value and cumulated filter values.
	The parameters gain, alpha and offset are constants.
	*/
	/* Ci(x,y)=log[a Ii(x,y)]-log[ Ei=1-s Ii(x,y)] */

	alpha  = 128.;
	gain   = 1.;
	offset = 0.;

	for (i = 0; i < size; i += bytes)
	{
		gfloat logl;

		psrc = src+i;
		pdst = dst+i;

		logl = log((gfloat)psrc[0] + (gfloat)psrc[1] + (gfloat)psrc[2] + 3.);

		pdst[0] = gain * ((log(alpha * (psrc[0]+1.)) - logl) * pdst[0]) + offset;
		pdst[1] = gain * ((log(alpha * (psrc[1]+1.)) - logl) * pdst[1]) + offset;
		pdst[2] = gain * ((log(alpha * (psrc[2]+1.)) - logl) * pdst[2]) + offset;
	}

	/*  if (!preview_mode)
	gimp_progress_update ((2.0 + (rvals.nscales * 3)) /
	((rvals.nscales * 3) + 3));*/

	/*
	Adapt the dynamics of the colors according to the statistics of the first and second order.
	The use of the variance makes it possible to control the degree of saturation of the colors.
	*/
	pdst = dst;

	compute_mean_var (pdst, &mean, &var, size, bytes);
	mini = mean - rvals.cvar*var;
	maxi = mean + rvals.cvar*var;
	range = maxi - mini;

	if (!range)
		range = 1.0;

	for (i = 0; i < size; i+= bytes)
	{
		psrc = src + i;
		pdst = dst + i;

		for (j = 0 ; j < 3 ; j++)
		{
			gfloat c = 255 * ( pdst[j] - mini ) / range;

			psrc[j] = (guchar) CLAMP<int> (c, 0, 255);
		}
	}

	g_free (dst);
}

/*
* Calculate the average and variance in one go.
*/
void
	compute_mean_var (gfloat *src, gfloat *mean, gfloat *var, gint size, gint bytes)
{
	gfloat vsquared;
	gint i,j;
	gfloat *psrc;

	vsquared = 0;
	*mean = 0;
	for (i = 0; i < size; i+= bytes)
	{
		psrc = src+i;
		for (j = 0 ; j < 3 ; j++)
		{
			*mean += psrc[j];
			vsquared += psrc[j] * psrc[j];
		}
	}

	*mean /= (gfloat) size; /* mean */
	vsquared /= (gfloat) size; /* mean (x^2) */
	*var = ( vsquared - (*mean * *mean) );
	*var = sqrt(*var); /* var */
}

