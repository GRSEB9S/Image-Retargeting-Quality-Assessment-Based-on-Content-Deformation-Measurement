#include "mex.h"
#include "math.h"
// undef needed for LCC compiler   undef需要LCC编译器
#undef EXTERN_C
#include <windows.h>
#include <process.h>    

/* Bspline transformation grid function
 * function [Iout,Tx,Ty]=bspline_transform_2d_double(Ox,Oy,Iin,dx,dy)
 * 
 * Ox, Oy are the grid points coordinates
 * Iin is input image, Iout the transformed output image
 * dx and dy are the spacing of the b-spline knots间距
 *
 * Iout: The transformed image
 * Tx: The transformation field in x direction
 * Ty: The transformation field in y direction
 *
 * This function is an implementation of the b-spline registration
 * algorithm in "D. Rueckert et al. : Nonrigid Registration Using Free-Form 
 * Deformations: Application to Breast MR Images".
 * 
 * We used "Fumihiko Ino et al. : a data distrubted parallel algortihm for 我们使用“Fumihiko Ino等人：
 * nonrigid image registration" for the correct formula's, because 用于非刚性图像配准的数据分布式并行算法”
 * (most) other papers contain errors. 用于正确的公式，因为（大多数）其他论文包含错误。
 *
 *   Function is written by D.Kroon University of Twente (July 2008)
 */

// Variables used to detect if the threads are finished用于检测线程是否完成的变量
// (volatile: reload the variable instead of using the value available in a register)
//（volatile：重新加载变量，而不是使用寄存器中的可用值）
static volatile int WaitForThread1;
static volatile int WaitForThread2;

// Convert 2D matrix index to 1D index  x,y是二维坐标的索引，sizx是大小，返回值就是一维的索引，坐标应该是从零开始
int mindex2(int x, int y, int sizx, int sizy) 
{
    if(x<0) { x=0; }
    if(x>(sizx-1)) { x=sizx-1; }
    if(y<0) { y=0; }
    if(y>(sizy-1)) { y=sizy-1; }
    return y*sizx+x;
}


unsigned __stdcall transformvolume(double **Args)//
{
    double *Bu, *Bv, *ZI, *Tx, *Ty, *dxa, *dya, *ThreadID, *Ox, *Oy, *ZO,*Tmap;//jw
    double *ZOsize_d;
    double *Osize_d;
    double *nlhs_d;
    int ZOsize[2]={0,0};
    int Osize[2]={0,0};

    // Multiple threads, one does the odd the other even indexes
    int offset;
        
    // Location of pixel which will be come the current pixel要到达当前像素的像素原位置
    double Tlocalx;
    double Tlocaly;
    
    // Variables to store 1D index
    int indexO;
    int indexZI;
    int indexZO;
    
    // Grid distance;
    int dx,dy; 
    
    // X,Y coordinates of current pixel
    int x,y;
    
    // Linear interpolation variables
    int xBas[4], yBas[4];
    double perc[4];
    double xCom, yCom;
    double color[4]={0,0,0,0};
    
    // B-spline variables
    int u_index=0; 
    int v_index=0;
    int i, j;
    
    // B-Spline loop variabels
    int l,m;
    int nlhs=0;
    
    // Split input into variables将输入拆分成变量
    Bu=Args[0];
    Bv=Args[1];
    ZOsize_d=Args[2];
    Osize_d=Args[3];
    ZI=Args[4];
    Tx=Args[5];
    Ty=Args[6];
    dxa=Args[7];
    dya=Args[8];
    ThreadID=Args[9];
    Ox=Args[10];
    Oy=Args[11];
    ZO=Args[12];
    nlhs_d=Args[13];
    Tmap=Args[14];//jw
   
    nlhs=(int)nlhs_d[0];
    ZOsize[0] = (int)ZOsize_d[0]; 
    ZOsize[1] = (int)ZOsize_d[1]; 
    Osize[0] = (int)Osize_d[0]; 
    Osize[1] = (int)Osize_d[1]; 
    
    /* Get the spacing of the uniform b-spline grid */
    dx=(int)dxa[0]; dy=(int)dya[0];
  
    
    if(ThreadID[0]==1) { offset = 0; }
    if(ThreadID[0]==2) { offset = 1; }
    
    // Loop through all image pixel coordinates奇偶线程分开做
    for (y=offset; y<ZOsize[1]; y=y+2)
    {
        for (x=0; x<ZOsize[0]; x++)
        {
            // Calculate the indexes need to loop up the B-spline values.取余，网格相对像素
            u_index=x%dx; 
            v_index=y%dy;
            
            i=(int)floor(x/(double)dx); // (first row outside image against boundary artefacts)编号
            j=(int)floor(y/(double)dy);
        
            // This part calculates the coordinates of the pixel这部分计算像素的坐标，这些像素将被转换到当前坐标 
            // which will be transformed to the current x,y pixel.
            Tlocalx=0; Tlocaly=0;
            for(l=0; l<4; l++)
            {
                for(m=0; m<4; m++)
                {    
                     if(((i+l)>=0)&&((i+l)<Osize[0])&&((j+m)>=0)&&((j+m)<Osize[1]))
                     {
                          indexO=mindex2(i+l,j+m,Osize[0],Osize[1]);
                          Tlocalx=Tlocalx+Bu[mindex2(l,u_index,4,dx)]*Bv[mindex2(m,v_index,4,dy)]*Ox[indexO];
                          Tlocaly=Tlocaly+Bu[mindex2(l,u_index,4,dx)]*Bv[mindex2(m,v_index,4,dy)]*Oy[indexO];
                     }
                }
            }            

            // Determine the coordinates of the pixel(s) which will be come the current pixel
            // (using linear interpolation)  
            xBas[0]=(int) floor(Tlocalx); yBas[0]=(int) floor(Tlocaly);
            xBas[1]=xBas[0]+0;      yBas[1]=yBas[0]+1;
            xBas[2]=xBas[0]+1;      yBas[2]=yBas[0]+0;
            xBas[3]=xBas[0]+1;      yBas[3]=yBas[0]+1;

            if(xBas[0]>=0&&xBas[0]<ZOsize[0]&&yBas[0]>=0&&yBas[0]<ZOsize[1]) { indexZI=mindex2(xBas[0],yBas[0],ZOsize[0],ZOsize[1]); color[0]=ZO[indexZI]; } else { color[0]=0;}
            if(xBas[1]>=0&&xBas[1]<ZOsize[0]&&yBas[1]>=0&&yBas[1]<ZOsize[1]) { indexZI=mindex2(xBas[1],yBas[1],ZOsize[0],ZOsize[1]); color[1]=ZO[indexZI]; } else { color[1]=0;}
            if(xBas[2]>=0&&xBas[2]<ZOsize[0]&&yBas[2]>=0&&yBas[2]<ZOsize[1]) { indexZI=mindex2(xBas[2],yBas[2],ZOsize[0],ZOsize[1]); color[2]=ZO[indexZI]; } else { color[2]=0;}
            if(xBas[3]>=0&&xBas[3]<ZOsize[0]&&yBas[3]>=0&&yBas[3]<ZOsize[1]) { indexZI=mindex2(xBas[3],yBas[3],ZOsize[0],ZOsize[1]); color[3]=ZO[indexZI]; } else { color[3]=0;}
  
            // Linear interpolation constants (percentages)
            xCom=Tlocalx-floor(Tlocalx); yCom=Tlocaly-floor(Tlocaly);
            perc[0]=(1-xCom) * (1-yCom);
            perc[1]=(1-xCom) * yCom;
            perc[2]=xCom * (1-yCom);
            perc[3]=xCom * yCom;

            // Set the current pixel value
            indexZO=mindex2(x,y,ZOsize[0],ZOsize[1]);
            ZI[indexZO]=color[0]*perc[0]+color[1]*perc[1]+color[2]*perc[2]+color[3]*perc[3];
            
            // Store transformation field Tx定义
            if(nlhs>1) { Tx[indexZO]=Tlocalx-(double)x; }
            if(nlhs>2) { Ty[indexZO]=Tlocaly-(double)y; }
            if(nlhs>3)
            {                                                           //jw   参照变化就和原图一样
              int     indexTmap=0;
              indexTmap=mindex2(xBas[0],yBas[0],ZOsize[0],ZOsize[1]);  
                if(ZO[indexZO]==1)
                {
                    Tmap[indexTmap]=1;
                }                                                                 //jw
            }
        }
    }
    
       
    // Set the thread finished variables
    if(ThreadID[0]==1) { WaitForThread1 = 0; }
    if(ThreadID[0]==2) { WaitForThread2 = 0; }
    
    // explicit end thread, helps to ensure proper recovery of resources allocated for the thread
    _endthreadex( 0 );
    return 0;
}




// The matlab mex function
void mexFunction( int nlhs, mxArray *plhs[],//输出个数，每个参数具体值
                  int nrhs, const mxArray *prhs[] )//输入，常量
{
    // Ox and Oy are the grid points
    // Zo is the input image
    // Zi is the transformed image
    // dx and dy are the spacing of the b-spline knots
    double *Ox,*Oy, *ZO, *dxa, *dya, *ZI,*Tx,*Ty,*Tmap;//jw
    
    // double pointer array to store all needed function variables
    double **ThreadArgs1,**ThreadArgs2;
    
    HANDLE *ThreadList; // Handles to the worker threads
    
    double ThreadID1[1]={1}; // ID of first Thread 
    double ThreadID2[1]={2}; // ID of second Thread
    
    double nlhs_d[1]={0};
    // Size of input image
    mwSize  ZOsizex, ZOsizey;
    double ZOsize_d[2]={0,0};

    // Size of grid
    mwSize  Osizex, Osizey;
    double Osize_d[2]={0,0};
    
    // B-spline variables
    double u,v;
    int u_index=0; 
    int v_index=0;
    double *Bu, *Bv;
         
    // X,Y coordinates of current pixel
    int x,y;
    // Grid distance;
    int dx,dy; 
    
    // Reserve room for 6 function variables(arrays)14---15
    ThreadArgs1 = (double **)malloc( 15* sizeof( double * ) );  
    ThreadArgs2 = (double **)malloc( 15* sizeof( double * ) );  
    
    // Reserve room for handles of threads in ThreadList
    ThreadList = (HANDLE*)malloc(2* sizeof( HANDLE ));
    
  /* Check for proper number of arguments. */
  if(nrhs!=5) {
    mexErrMsgTxt("Five inputs are required.");
  }
   
  // Get the sizes of the grid获得指针
  Osizex = (mwSize)mxGetM(prhs[0]);  
  Osizey = (mwSize)mxGetN(prhs[0]);
  
  // Create image matrix for the return arguments with the size of input image   
  ZOsizex = (mwSize) mxGetM(prhs[2]);  ZOsizey = (mwSize) mxGetN(prhs[2]);
  plhs[0] = mxCreateDoubleMatrix(ZOsizex,ZOsizey, mxREAL); 
  if(nlhs>1) { plhs[1] = mxCreateDoubleMatrix(ZOsizex,ZOsizey, mxREAL); }
  if(nlhs>2) { plhs[2] = mxCreateDoubleMatrix(ZOsizex,ZOsizey, mxREAL); }
  if(nlhs>3) { plhs[3] = mxCreateDoubleMatrix(ZOsizex,ZOsizey, mxREAL); }//jw，实双精度矩阵
  
  /* Assign pointers to each input. */
  Ox=mxGetPr(prhs[0]);
  Oy=mxGetPr(prhs[1]);
  ZO=mxGetPr(prhs[2]);
  dxa=mxGetPr(prhs[3]);
  dya=mxGetPr(prhs[4]);
  
  /* Get the spacing of the uniform b-spline grid */
  dx=(int)dxa[0]; dy=(int)dya[0];
  
    
  /* Assign pointer to output. */
  ZI = mxGetPr(plhs[0]);
  if(nlhs>1) { Tx = mxGetPr(plhs[1]); }
  if(nlhs>2) { Ty = mxGetPr(plhs[2]); }
  if(nlhs>3) { Tmap = mxGetPr(plhs[3]); }//jw
  
  // Make polynomial look up tables 进行多项式查找表
  Bu=malloc(dx*4*sizeof(double));//预定义内存
  Bv=malloc(dy*4*sizeof(double));
  for (x=0; x<dx; x++)
  {
    u=(x/(double)dx)-floor(x/(double)dx);//公式，floor向下取整直接去掉小数部分
    Bu[mindex2(0,x,4,dx)] = pow((1-u),3)/6;
    Bu[mindex2(1,x,4,dx)] = ( 3*pow(u,3) - 6*pow(u,2) + 4)/6;
    Bu[mindex2(2,x,4,dx)] = (-3*pow(u,3) + 3*pow(u,2) + 3*u + 1)/6;
    Bu[mindex2(3,x,4,dx)] = pow(u,3)/6;
  }
  
  for (y=0; y<dy; y++)
  {
    v=(y/(double)dy)-floor(y/(double)dy);
    Bv[mindex2(0,y,4,dy)] = pow((1-v),3)/6;    //(1-v)^3/6
    Bv[mindex2(1,y,4,dy)] = ( 3*pow(v,3) - 6*pow(v,2) + 4)/6;
    Bv[mindex2(2,y,4,dy)] = (-3*pow(v,3) + 3*pow(v,2) + 3*v + 1)/6;
    Bv[mindex2(3,y,4,dy)] = pow(v,3)/6;
  }
  
  ZOsize_d[0]=ZOsizex;  ZOsize_d[1]=ZOsizey;
  Osize_d[0]=Osizex;  Osize_d[1]=Osizey;
  
  nlhs_d[0]=(double)nlhs;
  
  WaitForThread1 = 1;
  ThreadArgs1[0]=Bu;
  ThreadArgs1[1]=Bv;
  ThreadArgs1[2]=ZOsize_d;
  ThreadArgs1[3]=Osize_d;
  ThreadArgs1[4]=ZI;
  ThreadArgs1[5]=Tx;
  ThreadArgs1[6]=Ty;
  ThreadArgs1[7]=dxa;
  ThreadArgs1[8]=dya;
  ThreadArgs1[9]=ThreadID1;
  ThreadArgs1[10]=Ox;
  ThreadArgs1[11]=Oy;
  ThreadArgs1[12]=ZO;
  ThreadArgs1[13]=nlhs_d;
  ThreadArgs1[14]=Tmap;//jw
  ThreadList[0] = (HANDLE)_beginthreadex( NULL, 0, &transformvolume, ThreadArgs1 , 0, NULL );
    
  WaitForThread2 = 1;
  ThreadArgs2[0]=Bu;
  ThreadArgs2[1]=Bv;
  ThreadArgs2[2]=ZOsize_d;
  ThreadArgs2[3]=Osize_d;
  ThreadArgs2[4]=ZI;
  ThreadArgs2[5]=Tx;
  ThreadArgs2[6]=Ty;
  ThreadArgs2[7]=dxa;
  ThreadArgs2[8]=dya;
  ThreadArgs2[9]=ThreadID2;
  ThreadArgs2[10]=Ox;
  ThreadArgs2[11]=Oy;
  ThreadArgs2[12]=ZO;
  ThreadArgs2[13]=nlhs_d;
  ThreadArgs2[14]=Tmap;//jw
  ThreadList[1] = (HANDLE)_beginthreadex( NULL, 0, &transformvolume, ThreadArgs2 , 0, NULL );

  
  WaitForSingleObject(ThreadList[0], INFINITE);
  WaitForSingleObject(ThreadList[1], INFINITE);
  
  CloseHandle( ThreadList[0] );
  CloseHandle( ThreadList[1] );

}
        

