/*
  CSC D18 - RayTracer code.

  Written Dec. 9 2010 - Jan 20, 2011 by F. J. Estrada
  Freely distributable for adacemic purposes only.

  Uses Tom F. El-Maraghi's code for computing inverse
  matrices. You will need to compile together with
  svdDynamic.c

  You need to understand the code provided in
  this file, the corresponding header file, and the
  utils.c and utils.h files. Do not worry about
  svdDynamic.c, we need it only to compute
  inverse matrices.

  You only need to modify or add code in sections
  clearly marked "TO DO" - remember to check what
  functionality is actually needed for the corresponding
  assignment!

  Last updated: Aug. 2017   - F.J.E.
*/

/*****************************************************************************
* COMPLETE THIS TEXT BOX:
*
* 1) Student Name:Yu Qian 
* 2) Student Name:Yue Gan
*
* 1) Student number:1000823460
* 2) Student number:1000620606
* 
* 1) UtorID qianyu3
* 2) UtorID ganyue
* 
* We hereby certify that the work contained here is our own
*
* __________Yu Qian__________             _______Yue Gan______________
* (sign with your name)            (sign with your name)
********************************************************************************/

#include "utils.h"	// <-- This includes RayTracer.h

// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct pointLS *light_list;
struct textureNode *texture_list;
int MAX_DEPTH;

void buildScene(void)
{
#include "buildscene.c"		// <-- Import the scene definition! 
}

void rtShade(struct object3D *obj, struct point3D *p, struct point3D *n, struct ray3D *ray, int depth, double a, double b, struct colourRGB *col)
{
  // This function implements the shading model as described in lecture. It takes
  // - A pointer to the first object intersected by the ray (to get the colour properties)
  // - The coordinates of the intersection point (in world coordinates)
  // - The normal at the point
  // - The ray (needed to determine the reflection direction to use for the global component, as well as for
  //   the Phong specular component)
  // - The current racursion depth
  // - The (a,b) texture coordinates (meaningless unless texture is enabled)
  //
  // Returns:
  // - The colour for this ray (using the col pointer)
  //

  struct colourRGB tmp_col; // Accumulator for colour components
  double R, G, B;           // Colour for the object in R G and B

  // This will hold the colour as we process all the components of
  // the Phong illumination model
  tmp_col.R = 0;
  tmp_col.G = 0;
  tmp_col.B = 0;

  if (obj->texImg == NULL) // Not textured, use object colour
  {
    R = obj->col.R;
    G = obj->col.G;
    B = obj->col.B;
  }
  else
  {
    // Get object colour from the texture given the texture coordinates (a,b), and the texturing function
    // for the object. Note that we will use textures also for Photon Mapping.
    obj->textureMap(obj->texImg, a, b, &R, &G, &B);
    // printf("in here\n\n\n");
  }

  //////////////////////////////////////////////////////////////
  // TO DO: Implement this function. Refer to the notes for
  // details about the shading model.
  //////////////////////////////////////////////////////////////

  // Be sure to update 'col' with the final colour computed here!

  double lambda = -1;

  struct point3D interToCam;
  struct point3D interToLS;
  struct ray3D *itcRay;
  struct ray3D *itlsRay;
  struct pointLS *tmp_LS = light_list;
  struct object3D *tmp_obj;
  struct point3D tmp_normal = *n;
  struct point3D m;
  struct point3D ms;
  struct ray3D *reflect_ray;
  struct colourRGB reflect_col;


  struct point3D tmp_p, tmp_n;
  double tmp_a, tmp_b;

  struct colourRGB refract_col;
  struct ray3D *refract_ray;
  

  interToCam.px = -ray->d.px;
  interToCam.py = -ray->d.py;
  interToCam.pz = -ray->d.pz;
  interToCam.pw = 1;
  normalize(&interToCam);
  itcRay = newRay(p, &interToCam);

  tmp_col.R = 0;
  tmp_col.G = 0;
  tmp_col.B = 0;

  refract_col.R = 0;
  refract_col.G = 0;
  refract_col.B = 0;
  

  // loop through lightsource
  while (tmp_LS != NULL)
  {
    // Local Component
    interToLS = tmp_LS->p0;
    subVectors(p, &interToLS);

    // See if there is anything inbetween this object and the lightsource
    itlsRay = newRay(p, &interToLS);
    findFirstHit(itlsRay, &lambda, obj, &tmp_obj, &tmp_p, &tmp_n, &tmp_a, &tmp_b);

    // Ambient term
    tmp_col.R += obj->alb.ra * R;
    tmp_col.G += obj->alb.ra * G;
    tmp_col.B += obj->alb.ra * B;

    
    if (lambda > 0 && lambda < 1)
    {
      //  tmp_col.R = 0;
      //  tmp_col.G = 0;
      //  tmp_col.B = 0;
      //  printf("lambda: %f\n", lambda);
    }
    else
    {

      normalize(&itlsRay->d);
      double dot_n_s = dot(n, &itlsRay->d);
      //  printf("lambda: %f\n", lambda);

      if (obj->frontAndBack && dot_n_s < 0)
      {
        dot_n_s = -dot_n_s;
      }

      // Diffuse
      tmp_col.R += obj->alb.rd * R * tmp_LS->col.R * max(0, dot_n_s);
      tmp_col.G += obj->alb.rd * G * tmp_LS->col.G * max(0, dot_n_s);
      tmp_col.B += obj->alb.rd * B * tmp_LS->col.B * max(0, dot_n_s);

      // Specular
      m.px = 2 * dot(&itlsRay->d, n) * n->px;
      m.py = 2 * dot(&itlsRay->d, n) * n->py;
      m.pz = 2 * dot(&itlsRay->d, n) * n->pz;
      m.pw = 1;
      subVectors(&itlsRay->d, &m);
      normalize(&m);

      tmp_col.R += obj->alb.rs * R * tmp_LS->col.R * pow(max(0, dot(&m, &itcRay->d)), obj->shinyness);
      tmp_col.G += obj->alb.rs * G * tmp_LS->col.G * pow(max(0, dot(&m, &itcRay->d)), obj->shinyness);
      tmp_col.B += obj->alb.rs * B * tmp_LS->col.B * pow(max(0, dot(&m, &itcRay->d)), obj->shinyness);
    }

    // Global Component
    if (depth < MAX_DEPTH)
    {
      double n_dot_d = 2 * dot(n, &ray->d);
      ms.px = n_dot_d * n->px;
      ms.py = n_dot_d * n->py;
      ms.pz = n_dot_d * n->pz;
      ms.pw = 1;

      subVectors(&ray->d, &ms);
      // We need to reverse te direction 
      ms.px = -ms.px;
      ms.py = -ms.py;
      ms.pz = -ms.pz;
      normalize(&ms);
      
      // specular reflection
      if (obj->alb.rs != 0)
      {

        reflect_ray = newRay(p, &ms);
        rayTrace(reflect_ray, depth + 1, &reflect_col, obj);

        tmp_col.R += obj->alb.rg * reflect_col.R;
        tmp_col.G += obj->alb.rg * reflect_col.G;
        tmp_col.B += obj->alb.rg * reflect_col.B;

      }

      // is refractive
       if (obj->alpha < 1) {

       }
    }
    tmp_LS = tmp_LS->next;
  }

  // printf("R:%f, G:%f, B:%f \n", tmp_col.R, tmp_col.G, tmp_col.B);
  // Add back to col
  if (tmp_col.R >= 1)
  {
    col->R = 1;
  }
  else
  {
    col->R = tmp_col.R;
  }

  if (tmp_col.G >= 1)
  {
    col->G = 1;
  }
  else
  {
    col->G = tmp_col.G;
  }

  if (tmp_col.B >= 1)
  {
    col->B = 1;
  }
  else
  {
    col->B = tmp_col.B;
  }

  return;
}

void findFirstHit(struct ray3D *ray, double *lambda, struct object3D *Os, struct object3D **obj, struct point3D *p, struct point3D *n, double *a, double *b)
{
 // Find the closest intersection between the ray and any objects in the scene.
 // Inputs:
 //   *ray    -  A pointer to the ray being traced
 //   *Os     -  'Object source' is a pointer toward the object from which the ray originates. It is used for reflected or refracted rays
 //              so that you can check for and ignore self-intersections as needed. It is NULL for rays originating at the center of
 //              projection
 // Outputs:
 //   *lambda -  A pointer toward a double variable 'lambda' used to return the lambda at the intersection point
 //   **obj   -  A pointer toward an (object3D *) variable so you can return a pointer to the object that has the closest intersection with
 //              this ray (this is required so you can do the shading)
 //   *p      -  A pointer to a 3D point structure so you can store the coordinates of the intersection point
 //   *n      -  A pointer to a 3D point structure so you can return the normal at the intersection point
 //   *a, *b  -  Pointers toward double variables so you can return the texture coordinates a,b at the intersection point

 /////////////////////////////////////////////////////////////
 // TO DO: Implement this function. See the notes for
 // reference of what to do in here
 /////////////////////////////////////////////////////////////

 *lambda = -1;
 double tmp_lambda = 0;
 struct point3D tmp_p, tmp_n;
 double tmp_a = 0;
 double tmp_b = 0;
//  printf("tmp_a = %f, tmp_b = %f\n", tmp_a, tmp_b);
printf("aaa");
 struct object3D *tmp_obj = object_list;

 while (tmp_obj != NULL)
 {
   if (tmp_obj != Os)
   {
     tmp_obj->intersect(tmp_obj, ray, &tmp_lambda, &tmp_p, &tmp_n, &tmp_a, &tmp_b);
     if ((*lambda < 0 || tmp_lambda < *lambda) && tmp_lambda > 0)
     {
       *lambda = tmp_lambda;
       *obj = tmp_obj;
       *p = tmp_p;
       *n = tmp_n;
       *a = tmp_a;
       *b = tmp_b;
     }
   }
   tmp_obj = tmp_obj->next;
 }
}

void rayTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os)
{
 // Trace one ray through the scene.
 //
 // Parameters:
 //   *ray   -  A pointer to the ray being traced
 //   depth  -  Current recursion depth for recursive raytracing
 //   *col   - Pointer to an RGB colour structure so you can return the object colour
 //            at the intersection point of this ray with the closest scene object.
 //   *Os    - 'Object source' is a pointer to the object from which the ray 
 //            originates so you can discard self-intersections due to numerical
 //            errors. NULL for rays originating from the center of projection. 
 
 double lambda;		// Lambda at intersection
 double a,b;		// Texture coordinates
 struct object3D *obj;	// Pointer to object at intersection
 struct point3D p;	// Intersection point
 struct point3D n;	// Normal at intersection
 struct colourRGB I;	// Colour returned by shading function

 if (depth>MAX_DEPTH)	// Max recursion depth reached. Return invalid colour.
 {
  col->R=-1;
  col->G=-1;
  col->B=-1;
  return;
 }

 ///////////////////////////////////////////////////////
 // TO DO: Complete this function. Refer to the notes
 // if you are unsure what to do here.
 ///////////////////////////////////////////////////////
 
 findFirstHit(ray, &lambda, Os, &obj, &p, &n, &a, &b);

//  printf("lambda: %f\n", lambda);
 if(lambda > 0){

    rtShade(obj, &p, &n, ray, depth, a, b, &I);
    // Color = rtShade();
    col->R = I.R;
    col->G = I.G;
    col->B = I.B;
  }
  // Else color is background

}

int main(int argc, char *argv[])
{
 // Main function for the raytracer. Parses input parameters,
 // sets up the initial blank image, and calls the functions
 // that set up the scene and do the raytracing.
 struct image *im;	// Will hold the raytraced image
 struct view *cam;	// Camera and view for this scene
 int sx;		// Size of the raytraced image
 int antialiasing;	// Flag to determine whether antialiaing is enabled or disabled
 char output_name[1024];	// Name of the output file for the raytraced .ppm image
 struct point3D e;		// Camera view parameters 'e', 'g', and 'up'
 struct point3D g;
 struct point3D up;
 double du, dv;			// Increase along u and v directions for pixel coordinates
 struct point3D pc,d;		// Point structures to keep the coordinates of a pixel and
				// the direction or a ray
 struct ray3D *ray;		// Structure to keep the ray from e to a pixel
 struct colourRGB col;		// Return colour for raytraced pixels
 struct colourRGB background;   // Background colour
 int i,j;			// Counters for pixel coordinates
 unsigned char *rgbIm;
 struct colourRGB antialiasing_col;

 if (argc<5)
 {
  fprintf(stderr,"RayTracer: Can not parse input parameters\n");
  fprintf(stderr,"USAGE: RayTracer size rec_depth antialias output_name\n");
  fprintf(stderr,"   size = Image size (both along x and y)\n");
  fprintf(stderr,"   rec_depth = Recursion depth\n");
  fprintf(stderr,"   antialias = A single digit, 0 disables antialiasing. Anything else enables antialiasing\n");
  fprintf(stderr,"   output_name = Name of the output file, e.g. MyRender.ppm\n");
  exit(0);
 }
 sx=atoi(argv[1]);
 MAX_DEPTH=atoi(argv[2]);
 if (atoi(argv[3])==0) antialiasing=0; else antialiasing=1;
 strcpy(&output_name[0],argv[4]);

 fprintf(stderr,"Rendering image at %d x %d\n",sx,sx);
 fprintf(stderr,"Recursion depth = %d\n",MAX_DEPTH);
 if (!antialiasing) fprintf(stderr,"Antialising is off\n");
 else fprintf(stderr,"Antialising is on\n");
 fprintf(stderr,"Output file name: %s\n",output_name);

 object_list=NULL;
 light_list=NULL;
 texture_list=NULL;

 // Allocate memory for the new image
 im=newImage(sx, sx);
 if (!im)
 {
  fprintf(stderr,"Unable to allocate memory for raytraced image\n");
  exit(0);
 }
 else rgbIm=(unsigned char *)im->rgbdata;

 ///////////////////////////////////////////////////
 // TO DO: You will need to implement several of the
 //        functions below. For Assignment 2, you can use
 //        the simple scene already provided. But
 //        for Assignment 3 you need to create your own
 //        *interesting* scene.
 ///////////////////////////////////////////////////
 buildScene();		// Create a scene. This defines all the
			// objects in the world of the raytracer

 //////////////////////////////////////////
 // TO DO: For Assignment 2 you can use the setup
 //        already provided here. For Assignment 3
 //        you may want to move the camera
 //        and change the view parameters
 //        to suit your scene.
 //////////////////////////////////////////

 // Mind the homogeneous coordinate w of all vectors below. DO NOT
 // forget to set it to 1, or you'll get junk out of the
 // geometric transformations later on.

 // Camera center is at (0,0,-1)
 e.px=0;
 e.py=0;
 e.pz=-1;
 e.pw=1;

 // To define the gaze vector, we choose a point 'pc' in the scene that
 // the camera is looking at, and do the vector subtraction pc-e.
 // Here we set up the camera to be looking at the origin.
 g.px=0-e.px;
 g.py=0-e.py;
 g.pz=0-e.pz;
 g.pw=1;
 // In this case, the camera is looking along the world Z axis, so
 // vector w should end up being [0, 0, -1]

 // Define the 'up' vector to be the Y axis
 up.px=0;
 up.py=1;
 up.pz=0;
 up.pw=1;

 // Set up view with given the above vectors, a 4x4 window,
 // and a focal length of -1 (why? where is the image plane?)
 // Note that the top-left corner of the window is at (-2, 2)
 // in camera coordinates.
 cam=setupView(&e, &g, &up, -1, -2, 2, 4);

 if (cam==NULL)
 {
  fprintf(stderr,"Unable to set up the view and camera parameters. Our of memory!\n");
  cleanup(object_list,light_list, texture_list);
  deleteImage(im);
  exit(0);
 }

 // Set up background colour here
 background.R=0;
 background.G=0;
 background.B=0;

 // Do the raytracing
 //////////////////////////////////////////////////////
 // TO DO: You will need code here to do the raytracing
 //        for each pixel in the image. Refer to the
 //        lecture notes, in particular, to the
 //        raytracing pseudocode, for details on what
 //        to do here. Make sure you undersand the
 //        overall procedure of raytracing for a single
 //        pixel.
 //////////////////////////////////////////////////////
 du=cam->wsize/(sx-1);		// du and dv. In the notes in terms of wl and wr, wt and wb,
 dv=-cam->wsize/(sx-1);		// here we use wl, wt, and wsize. du=dv since the image is
				// and dv is negative since y increases downward in pixel
				// coordinates and upward in camera coordinates.

 fprintf(stderr,"View parameters:\n");
 fprintf(stderr,"Left=%f, Top=%f, Width=%f, f=%f\n",cam->wl,cam->wt,cam->wsize,cam->f);
 fprintf(stderr,"Camera to world conversion matrix (make sure it makes sense!):\n");
 printmatrix(cam->C2W);
 fprintf(stderr,"World to camera conversion matrix:\n");
 printmatrix(cam->W2C);
 fprintf(stderr,"\n");

 fprintf(stderr,"Rendering row: ");
 for (j=0;j<sx;j++)		// For each of the pixels in the image
 {
  fprintf(stderr,"%d/%d, ",j,sx);
  for (i=0;i<sx;i++)
  {
    ///////////////////////////////////////////////////////////////////
    // TO DO - complete the code that should be in this loop to do the
    //         raytracing!
    ///////////////////////////////////////////////////////////////////

    if (antialiasing)
    {
      antialiasing_col.R = 0;
      antialiasing_col.G = 0;
      antialiasing_col.B = 0;

      for (int ai = -1; ai <= 1; ai += 2)
      {
        for (int aj = -1; aj <= 1; aj += 2)
        {
          // Set up d and p0
          d.px = cam->wl + i * du + du / (2 * ai);
          d.py = cam->wt + j * dv + dv / (2 * aj);
          d.pz = -1;
          d.pw = 0;

          pc.px = cam->wl + i * du + du / (2 * ai);
          pc.py = cam->wt + j * dv + dv / (2 * aj);
          pc.pz = -1;
          pc.pw = 1;

          matVecMult(cam->C2W, &d);
          normalize(&d);
          matVecMult(cam->C2W, &pc);

          ray = newRay(&pc, &d);
          col = background;

          rayTrace(ray, 0, &col, NULL);

          antialiasing_col.R += col.R / 4;
          antialiasing_col.G += col.G / 4;
          antialiasing_col.B += col.B / 4;
        }
      }

      // // setPixel(i,j,I)

      
      *(rgbIm + ((i + (j * sx)) * 3) + 0) = (unsigned char)(255 * antialiasing_col.R);
      *(rgbIm + ((i + (j * sx)) * 3) + 1) = (unsigned char)(255 * antialiasing_col.G);
      *(rgbIm + ((i + (j * sx)) * 3) + 2) = (unsigned char)(255 * antialiasing_col.B);
    }
    else
    {
      // Set up d and p0
      d.px = cam->wl + i * du;
      d.py = cam->wt + j * dv;
      d.pz = -1;
      d.pw = 0;

      pc.px = cam->wl + i * du;
      pc.py = cam->wt + j * dv;
      pc.pz = -1;
      pc.pw = 1;

      matVecMult(cam->C2W, &d);
      normalize(&d);
      matVecMult(cam->C2W, &pc);

      ray = newRay(&pc, &d);
      col = background;

      rayTrace(ray, 0, &col, NULL);

      // setPixel(i,j,I)
      *(rgbIm + ((i + (j * sx)) * 3) + 0) = (unsigned char)(255 * col.R);
      *(rgbIm + ((i + (j * sx)) * 3) + 1) = (unsigned char)(255 * col.G);
      *(rgbIm + ((i + (j * sx)) * 3) + 2) = (unsigned char)(255 * col.B);

    }

  } // end for i
 } // end for j

 fprintf(stderr,"\nDone!\n");

 // Output rendered image
 imageOutput(im,output_name);

 // Exit section. Clean up and return.
 cleanup(object_list,light_list,texture_list);		// Object, light, and texture lists
 deleteImage(im);					// Rendered image
 free(cam);						// camera view
 exit(0);
}

