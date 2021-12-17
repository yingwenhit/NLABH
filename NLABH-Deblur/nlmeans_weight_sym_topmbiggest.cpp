/*
 * =============================================================
 * nlemans_weight_sym.c 
 *Input: 
             *X: .................... Input image
             *h: .................... Smooth Parameter
             *nwin: ................. half of patch size 
             *bloc .................. half of block size
 *  patch size: [2*nwin+1, 2*nwin+1],
 %  search window: [2*nbloc+1, 2*nbloc+1]

 * This is a MEX-file for MATLAB.
 * =============================================================
 */

/* Revision: 1.0, change the patch and bloc size to half  */
/*Copyright(C) Xiaoqun zhang 16/04/2007*/

#include "mex.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <ctype.h>
#include <vector>
#include <numeric> // std::iota
#include <algorithm>

using namespace std;

// #define SIZE 7

// struct value_index_pair
// {
//     double value;
//     int index;
// };

// int comp(const void *a, const void *b)
// {
//     return (*(struct value_index_pair *)a).value - (*(struct value_index_pair *)b).value;
// }

template <typename T>
vector<size_t> sort_indexes(const vector<T> &v)
{

  // initialize original index locations
  vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i)
    idx[i] = i;

  // sort indexes based on comparing values in v
  sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) { return v[i1] > v[i2]; });

  return idx;
}

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

/*****************************************/
/**************The gateway routine.****  */
/*****************************************/

void mexFunction(int nlhs, mxArray *plhs[],       //output
                 int nrhs, const mxArray *prhs[]) //input
{
  /*input for C function */
  double *in;     /*input image*/
  double *weight; /*output weight*/
  double h;
  mwIndex nxny, nzmax;
  int nx, ny, newnx, newny;
  int nwin, bloc;
  double a, percent_sparse, v, c, sum, e, dist;
  int i, j, k;
  int *dadr, *dd, x, y, xp, yp, adr, adrp, psize, wsize, ds, realadr, realadrp, li;
  double *w, *ww, *lw;
  /*symetry image*/
  double *in_sym;
  mwIndex *irs, *jcs;
  char s;
  int adrp1, adrp2, adrp3, adrp4;
  int ii, ts;
  //
  /*Output for C function */

  int l;
  if (nrhs < 1)
    mexErrMsgTxt("At least one inputs required.");
  if (nrhs > 6)
    mexErrMsgTxt("Too many inputs required.");
  if (nlhs != 1)
    mexErrMsgTxt("One output required.");

  /*Input Image */
  nx = mxGetM(prhs[0]);
  ny = mxGetN(prhs[0]);
  in = mxGetPr(prhs[0]);
  h = mxGetScalar(prhs[1]);
  nxny = nx * ny;

  if (nrhs < 3)
    nwin = 3;
  else
    nwin = mxGetScalar(prhs[2]);

  if (nrhs < 4)
    bloc = 6;
  else
    bloc = mxGetScalar(prhs[3]);


  // mexPrintf("-------- nwin: %d, bloc: %d --------\n", nwin, bloc);

  /*Other parameters*/
  percent_sparse = (2 * bloc + 1) * (2 * bloc + 1) / (double)nxny;
  nzmax = (mwSize)ceil((double)nxny * (double)nxny * percent_sparse);

  plhs[0] = mxCreateSparse(nxny, nxny, nzmax, mxREAL);
  weight = mxGetPr(plhs[0]);

  /* Create a C pointer to a copy of the output matrix,but be careful that the matrix is tranposed */
  a = (2 * nwin) / 4.; //wy, the alpha
  c = 1;

  /*Sparse weight paramters*/
  irs = mxGetIr(plhs[0]);
  jcs = mxGetJc(plhs[0]);

  /* precompute weights */
  ds = nwin;
  psize = (2 * nwin + 1) * (2 * nwin + 1); // patch size;  patch = [-ds,ds]x[-ds,ds]
  wsize = (2 * bloc + 1) * (2 * bloc + 1); // window size;  patch = [-d,d]x[-d,d]

  w = (double *)malloc(psize * sizeof(double));  //patch weight
  dadr = (int *)malloc(psize * sizeof(int));     //patch index
  lw = (double *)malloc(wsize * sizeof(double)); //local window weight

  /*symstry bord of image*/

  newnx = nx + 2 * ds;
  newny = ny + 2 * ds;

  in_sym = (double *)malloc(newny * newnx * sizeof(double));

  // extend the original image
  for (j = 0, k = 0; j < newny; j++)
  {
    if (j < ds)
      y = ds - j;
    else if (j > ny + ds - 1)
      y = 2 * ny + ds - j - 2;
    else
      y = j - ds;

    for (i = 0; i < newnx; i++, k++)
    {
      if (i <= ds)
        x = ds - i;
      else if (i > nx + ds - 1)
        x = 2 * nx + ds - i - 2;
      else
        x = i - ds;
      in_sym[k] = in[y * nx + x];
    }
  }

  vector<double> vec; // = {5, 31, 9, 11, 8, 21, 9, 7, 4};
  vector<int> tmprealadrp;
  vector<size_t> idx;
  //****** wy, end

  //Precompute patch weight and index
  for (sum = 0., i = 0, x = -ds; x <= ds; x++)
    for (y = -ds; y <= ds; y++, i++)
    {
      dadr[i] = y * newnx + x;
      w[i] = exp(-(double)(x * x + y * y) / (2. * a * a));
      sum += w[i];
    }

  for (i = psize; i--;)
    w[i] /= 2 * sum * h * h;




  // Principal loop
  //****** wy

  k = 0;
  for (y = ds; y < newny - ds; y++)
    for (x = ds; x < newnx - ds; x++)
    {
      adr = y * newnx + x;
      realadr = (y - ds) * nx + x - ds;
      sum = 0.;
      jcs[realadr] = k;

      adrp1 = (y - ds) * nx + x - 1 - ds;
      adrp2 = (y - ds) * nx + x + 1 - ds;
      adrp3 = (y - 1 - ds) * nx + x - ds;
      adrp4 = (y + 1 - ds) * nx + x - ds;

      /* loop on patches */
      // Count non zero weight for realadr columns
      int nc = 0;
      for (xp = MAX(x - bloc, ds); xp <= MIN(x + bloc, nx - 1 + ds); xp++)
        for (yp = MAX(y - bloc, ds); yp <= MIN(y + bloc, ny - 1 + ds); yp++)
        {
          
          adrp = yp * newnx + xp;
          realadrp = (yp - ds) * nx + xp - ds;
          for (i = psize, dist = 0., ww = w, dd = dadr; i--; ww++, dd++)
          {
            v = in_sym[adr + *dd] - in_sym[adrp + *dd];
            dist += *ww * (double)(v * v);
          }
          e = (adrp == adr ? c : exp(-dist));

          //********** wy, order top bloc biggest weight
          // save realadrp & e, and then order by e // tmprealadrp=wsize *
          // if (realadr == 4452)
          // {
          //   mexPrintf("e: %f, xp: %d, yp: %d, realadrp: %d \n", e, xp-ds, yp-ds, realadrp);
          //   /* code */
          // }
          
          vec.push_back(e);
          tmprealadrp.push_back(realadrp);
          nc++;

          if (realadrp == adrp1 | realadrp == adrp2 | realadrp == adrp3 | realadrp == adrp4)
          {
            // if (realadr == 5675)
            // {
            //   mexPrintf("%f, %d\n", e, realadrp);
            // }
            weight[k] = e;
            irs[k] = realadrp;
            k++;
          }
        }

      idx = sort_indexes(vec);
      // if (realadr == 5675)
      // {
      //   for (i=0; i<50; i++) //vec.size()
      //   {
      //     mexPrintf("%f, %d\n", vec[idx[i]], tmprealadrp[idx[i]]);
      //   // cout << a << " ";
      //   }   
      //   // cout << '\n';
      //   /* code */
      //   mexPrintf("------------------------\n");
      // }
      

      ts = MIN(bloc, idx.size());
      // ts = idx.size();
      // ts = 30;
      for (ii = 0; ii < ts; ii++)
      {
        int tmpd = tmprealadrp.at(idx.at(ii));
        if (tmpd == adrp1 | tmpd == adrp2 | tmpd == adrp3 | tmpd == adrp4)
          continue;
        weight[k] = vec.at(idx.at(ii));
        // if (realadr == 5675)
        // {
        //   mexPrintf("%f, %d, %d\n", vec.at(idx.at(ii)), ii, tmpd);
        // }
        irs[k] = tmpd;
        //mexPrintf("---- %d, %d, %f, %d----\n", ii, idx.at(ii), vec.at(idx.at(ii)), irs[k]);
        k++;
      }
      vec.clear();
      tmprealadrp.clear();
      idx.clear();
    }
  jcs[nxny] = k;
  free(dadr);
  free(w);
  free(in_sym);
}
