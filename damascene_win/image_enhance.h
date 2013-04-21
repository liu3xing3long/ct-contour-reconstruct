#ifndef image_enhance_h__
#define image_enhance_h__

typedef int gint;
typedef float gfloat;
typedef unsigned char guchar;
typedef int gint;
typedef double gdouble;
typedef bool gboolean;


typedef struct
{
  gint     scale;
  gint     nscales;
  gint     scales_mode;
  gfloat   cvar;
} RetinexParams;

typedef enum
{
  filter_uniform,
  filter_low,
  filter_high
} FilterMode;

#define RETINEX_UNIFORM 0
#define RETINEX_LOW     1
#define RETINEX_HIGH    2


typedef struct
{
  gint    N;
  gfloat  sigma;
  gdouble B;
  gdouble b[4];
} gauss3_coefs;

/*
 * Private variables.
 */
static RetinexParams rvals =
{
  240,             /* Scale */
  3,               /* Scales */
  RETINEX_UNIFORM, /* Echelles reparties uniformement */
  1.2              /* A voir */
};

/*
 * Declare local functions.
 */
void     retinex_scales_distribution (gfloat       *scales,
                                             gint          nscales,
                                             gint          mode,
                                             gint          s);

void     compute_mean_var            (gfloat       *src,
                                             gfloat       *mean,
                                             gfloat       *var,
                                             gint          size,
                                             gint          bytes);
/*
 * Gauss
 */
void     compute_coefs3              (gauss3_coefs *c,
                                             gfloat        sigma);

void     gausssmooth                 (gfloat       *in,
                                             gfloat       *out,
                                             gint          size,
                                             gint          rowtride,
                                             gauss3_coefs *c);

/*
 * MSRCR = MultiScale Retinex with Color Restoration
 */
void     MSRCR                       (guchar       *src,
                                             gint          width,
                                             gint          height,
                                             gint          bytes,
                                             gboolean      preview_mode);



#endif // image_enhance_h__
