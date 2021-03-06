/*
 CSC D18 - Path Tracer code.
 
 Derived from the ray tracer starter code. Most function
 names are identical, though in practice the implementation
 should be much simpler!
 
 You only need to modify or add code in sections
 clearly marked "TO DO" - remember to check what
 functionality is actually needed for the corresponding
 assignment!
 
 Last updated: Aug. 2017   - F.J.E.
 */

/*****************************************************************************
 * COMPLETE THIS TEXT BOX:
 *
 * 1) Student Name:
 * 2) Student Name:
 *
 * 1) Student number:
 * 2) Student number:
 *
 * 1) UtorID
 * 2) UtorID
 *
 * We hereby certify that the work contained here is our own
 *
 * ____________________             _____________________
 * (sign with your name)            (sign with your name)
 ********************************************************************************/

#include "utils_path.h" // <-- This includes PathTracer.h
//#define __USE_IS            // Use importance sampling for diffuse materials
//#define __USE_ES            // Use explicit light sampling
//#define __DEBUG            // <-- Use this to turn on/off debugging output

// A couple of global structures and data: An object list, a light list, and the
// maximum recursion depth
struct object3D *object_list;
struct textureNode *texture_list;
unsigned long int NUM_RAYS;
int MAX_DEPTH;

#include "buildScene.c" // Import scene definition

void randomLSPoint(double *area, double *x, double *y, double *z)
{
  struct object3D *tmp_obj = object_list;
  // double x,y,z;
  while (tmp_obj != NULL)
  {
    if (tmp_obj->isLightSource == 1)
    {
      planeSample(tmp_obj, x, y, z);
    }
    tmp_obj = tmp_obj->next;
  }
  area = (double *)5;
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

void PathTrace(struct ray3D *ray, int depth, struct colourRGB *col, struct object3D *Os, int CEL)
{
  // Trace one light path through the scene.
  //
  // Parameters:
  //   *ray   -  A pointer to the ray being traced
  //   depth  -  Current recursion depth for recursive raytracing
  //   *col   - Pointer to an RGB colour structure so you can return the object colour
  //            at the intersection point of this ray with the closest scene object.
  //   *Os    - 'Object source' is a pointer to the object from which the ray
  //            originates so you can discard self-intersections due to numerical
  //            errors. NULL for rays originating from the center of projection.

  double lambda;          // Lambda at intersection
  double a, b;            // Texture coordinates
  struct object3D *obj;   // Pointer to object at intersection
  struct point3D p;       // Intersection point
  struct point3D n;       // Normal at intersection
  double R, G, B;         // Handy in case you need to keep track of some RGB colour value
  double dice;            // Handy to keep a random value
  struct ray3D *next_ray; // For the new ray to be used in recursive calls
  double BRDF;
  double w; //weight for explicit sampling
  //w=(2 *Pi * A_ls * dot(N_o,d) * dot (N_ls, d)) / d^2

  if (depth > MAX_DEPTH) // Max recursion depth reached. Return black (no light coming into pixel from this path).
  {
    col->R = ray->Ir; // These are accumulators, initialized at 0. Whenever we find a source of light these
    col->G = ray->Ig; // get incremented accordingly. At the end of the recursion, we return whatever light
    col->B = ray->Ib; // we accumulated into these three values.
    return;
  }

  ///////////////////////////////////////////////////////
  // TO DO: Complete this function. Refer to the notes
  // if you are unsure what to do here.
  ///////////////////////////////////////////////////////
  findFirstHit(ray, &lambda, Os, &obj, &p, &n, &a, &b);

  if (lambda > 0)
  {

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
      // printf("a:%f, b:%f\n", a, b);
      obj->textureMap(obj->texImg, a, b, &R, &G, &B);
    }
    if (obj->isLightSource == 1)
    {

      // If CEL == 0 means previous step hit a LS
      // if (CEL == 0)
      // {
      //   return;
      // }
      // else
      // {
      // printf("\nRGB before R:%f, G:%f, B:%f \n", col->R, col->G, col->B);
      col->R = ray->R * R;
      col->G = ray->G * G;
      col->B = ray->B * B;

      // printf("RGB only R:%f, G:%f, B:%f \n", ray->R * obj->col.R, ray->G * obj->col.G, ray->B * obj->col.B);
      // printf("RGB plus R:%f, G:%f, B:%f \n", ray->R * obj->col.R + col->R, ray->G * obj->col.G + col->G, ray->B * obj->col.B + col->B);

      // col->R += ray->R * obj->col.R;
      // col->G += ray->G * obj->col.G;
      // col->B += ray->B * obj->col.B;
      // printf("RGB after R:%f, G:%f, B:%f \n", col->R, col->G, col->B);

      // }
    }
    else
    {
      dice = drand48();
      double material;
      //choose the surface material
      helper1(obj->diffPct, obj->reflPct, obj->tranPct, &material);

      if (material == 0)
      {

        //diffuse

        // Do Explicit sampling
        // double x,y,z;
        // double w;
        // double area;

        // randomLSPoint(&area, &x, &y, &z);

        // struct ray3D *tmp_ls_ray;
        // struct point3D *tmp_ls_dir;
        // struct object3D *tmp_obj;
        // struct point3D tmp_p, tmp_n, flip_dir;
        // double tmp_a, tmp_b;
        // double tmp_lambda;

        // tmp_ls_dir = newPoint(x,y,z);
        // subVectors(&p, tmp_ls_dir);

        // tmp_ls_ray = newRay(&p, tmp_ls_dir);

        // findFirstHit(tmp_ls_ray, &tmp_lambda, obj, &tmp_obj, &tmp_p, &tmp_n, &tmp_a, &tmp_b);

        // if (tmp_lambda > 0 && tmp_lambda < 1)
        // {
        //   CEL = 1;
        // }
        // else
        // {
        //   tmp_ls_dir->pw = 0;
        //   normalize(tmp_ls_dir);
        //   tmp_ls_dir->pw = 1;
        //   flip_dir.px = -tmp_ls_dir->px;
        //   flip_dir.py = -tmp_ls_dir->py;
        //   flip_dir.pz = -tmp_ls_dir->pz;
        //   flip_dir.pw = 0;

        //   w = (2 * PI * area * dot(&n, tmp_ls_dir) * dot(&tmp_n, &flip_dir)) / pow(dot(tmp_ls_dir, tmp_ls_dir),2);
        //   // printf("weight:%f\n", w);
        //   col->R += ray->R * obj->col.R * 1 * w;
        //   col->G += ray->G * obj->col.G * 1 * w;
        //   col->B += ray->B * obj->col.B * 1 * w;

        //   // printf("col R:%f, G:%f, B:%f\n", ray->R * obj->col.R * 1 * w, ray->G * obj->col.G * 1 * w, ray->B * obj->col.B * 1 * w);
        //   CEL = 0;
        // }

        struct point3D dir;
        //random sampling
        // uniformSample(&n, &dir);
        //importance sampling
        cosWeightedSample(&n, &dir);
        next_ray = newRay(&p, &dir);
        BRDF = dot(&n, &dir);
        if (BRDF < 0)
        {
          BRDF = -BRDF;
        }
        R = ray->R * R * BRDF;
        G = ray->G * G * BRDF;
        B = ray->B * B * BRDF;
        //update next ray to trace
        next_ray->R = R;
        next_ray->G = G;
        next_ray->B = B;

        //find maximum of R,G,B
        double max;
        maximumRGB(R, G, B, &max);
        dice = drand48() * 0.2;
        if (dice < max)
        {
          PathTrace(next_ray, depth + 1, col, obj, CEL);
        }
        else
        {
          //if we killed the light, return background color
          // col->R=ray->Ir + col->R;
          // col->G=ray->Ig + col->G;
          // col->B=ray->Ib + col->B;

          col->R = ray->Ir;
          col->G = ray->Ig;
          col->B = ray->Ib;

          // col->R = col->R + 0;
          // col->G = col->G + 0;
          // col->B = col->B + 0;
          return;
        }
      }
      else if (material == 1)
      {

        //refraction
        double dot_n_d;
        double n1;
        double n2;
        double r_idx;
        double r_0;
        dot_n_d = dot(&n, &ray->d);

        struct point3D flipNormal;
        flipNormal.px = -n.px;
        flipNormal.py = -n.py;
        flipNormal.pz = -n.pz;
        flipNormal.pw = 1;

        if (dot_n_d > 0)
        {
          n1 = obj->r_index;
          n2 = 1;
        }
        else
        {
          n1 = 1;
          n2 = obj->r_index;
        }

        r_idx = n1 / n2;

        if (r_idx > 0)
        {
          r_0 = ((1 - r_idx) / (1 + r_idx)) * ((1 - r_idx) / (1 + r_idx));
        }
        else
        {
          r_0 = ((r_idx - 1) / (r_idx + 1)) * ((r_idx - 1) / (r_idx + 1));
        }

        double angel = acos(dot(&ray->d, &flipNormal) / (dot(&ray->d, &ray->d) * dot(&flipNormal, &flipNormal)));

        double r_theta = r_0 + (pow((1 - r_0) * (1 - cos(angel)), 5));
        double Rt = 1 - r_theta;

        dice = drand48();
        if (dice <= Rt)
        {
          // Refraction
          struct point3D tmp_dir;
          double c = dot(&flipNormal, &ray->d);
          double sqrt_value = 1 - pow(r_idx, 2) * (1 - pow(c, 2));

          R = ray->R * R;
          G = ray->G * G;
          B = ray->B * B;

          if (sqrt_value >= 0)
          {
            tmp_dir.px = r_idx * ray->d.px + (r_idx * c - sqrt(sqrt_value)) * n.px;
            tmp_dir.py = r_idx * ray->d.py + (r_idx * c - sqrt(sqrt_value)) * n.py;
            tmp_dir.pz = r_idx * ray->d.pz + (r_idx * c - sqrt(sqrt_value)) * n.pz;
            tmp_dir.pw = 1;

            next_ray = newRay(&p, &tmp_dir);
            next_ray->R = R;
            next_ray->G = G;
            next_ray->B = B;

            //find the max
            double max;
            maximumRGB(R, G, B, &max);
            if (dice < max)
            {
              // (next_ray->rayPos)(next_ray, .00001, &next_ray->p0);
              PathTrace(next_ray, depth + 1, col, NULL, 1);
              // col->R = col->R * Rt;
              // col->G = col->G * Rt;
              // col->B = col->B * Rt;
            }
            else
            {
              // col->R = ray->Ir + col->R;
              // col->G = ray->Ig + col->G;
              // col->B = ray->Ib + col->B;

              col->R = ray->Ir;
              col->G = ray->Ig;
              col->B = ray->Ib;

              // col->R = col->R + 0;
              // col->G = col->G + 0;
              // col->B = col->B + 0;
              free(next_ray);
              return;
            }
            free(next_ray);
          }
        }
        else
        {
          // Reflection
          struct point3D ms;
          R = ray->R * R;
          G = ray->G * G;
          B = ray->B * B;
          double n_dot_d = 2 * dot(&n, &ray->d);
          ms.px = n_dot_d * n.px;
          ms.py = n_dot_d * n.py;
          ms.pz = n_dot_d * n.pz;
          ms.pw = 1;
          subVectors(&ray->d, &ms);

          normalize(&ms);
          // We need to reverse te direction
          ms.px = -ms.px;
          ms.py = -ms.py;
          ms.pz = -ms.pz;

          //update next ray to trace
          next_ray = newRay(&p, &ms);
          next_ray->R = R;
          next_ray->G = G;
          next_ray->B = B;

          //find the max
          double max;
          maximumRGB(R, G, B, &max);
          dice = drand48() * 0.2;
          if (dice < max)
          {
            PathTrace(next_ray, depth + 1, col, obj, 1);
          }
          else
          {
            col->R = ray->Ir + col->R;
            col->G = ray->Ig + col->G;
            col->B = ray->Ib + col->B;

            // col->R = ray->Ir;
            // col->G = ray->Ig;
            // col->B = ray->Ib;

            // col->R = col->R + 0;
            // col->G = col->G + 0;
            // col->B = col->B + 0;
            return;
          }
        }
      }
      else
      {
        //reflection
        struct point3D ms;
        R = ray->R * R;
        G = ray->G * G;
        B = ray->B * B;
        double n_dot_d = 2 * dot(&n, &ray->d);
        ms.px = n_dot_d * n.px;
        ms.py = n_dot_d * n.py;
        ms.pz = n_dot_d * n.pz;
        ms.pw = 1;
        subVectors(&ray->d, &ms);

        normalize(&ms);
        // We need to reverse te direction
        ms.px = -ms.px + obj->refl_sig * ((drand48() * 4) - 2);
        ms.py = -ms.py + obj->refl_sig * ((drand48() * 4) - 2);
        ms.pz = -ms.pz + obj->refl_sig * ((drand48() * 4) - 2);

        //update next ray to trace
        next_ray = newRay(&p, &ms);
        next_ray->R = R;
        next_ray->G = G;
        next_ray->B = B;

        //find the max
        double max;
        maximumRGB(R, G, B, &max);
        dice = drand48() * 0.2;
        if (dice < max)
        {
          PathTrace(next_ray, depth + 1, col, obj, 1);
        }
        else
        {
          // col->R=ray->Ir + col->R;
          // col->G=ray->Ig + col->G;
          // col->B=ray->Ib + col->B;

          col->R = ray->Ir;
          col->G = ray->Ig;
          col->B = ray->Ib;

          // col->R = col->R + 0;
          // col->G = col->G + 0;
          // col->B = col->B + 0;
          return;
        }
      }
    }
  }
}

int main(int argc, char *argv[])
{
  // Main function for the path tracer. Parses input parameters,
  // sets up the initial blank image, and calls the functions
  // that set up the scene and do the raytracing.
  struct image *im;       // Will hold the final image
  struct view *cam;       // Camera and view for this scene
  int sx;                 // Size of the  image
  int num_samples;        // Number of samples to use per pixel
  char output_name[1024]; // Name of the output file for the .ppm image file
  struct point3D e;       // Camera view parameters 'e', 'g', and 'up'
  struct point3D g;
  struct point3D up;
  double du, dv;        // Increase along u and v directions for pixel coordinates
  struct point3D pc, d; // Point structures to keep the coordinates of a pixel and
  // the direction or a ray
  struct ray3D *ray;    // Structure to keep the ray from e to a pixel
  struct colourRGB col; // Return colour for pixels
  int i, j, k;          // Counters for pixel coordinates and samples
  double *rgbIm;        // Image is now double precision floating point since we
  // will be accumulating brightness differences with a
  // wide dynamic range
  struct object3D *obj; // Will need this to process lightsource weights
  double *wght;         // Holds weights for each pixel - to provide log response
  double pct, wt;

  time_t t1, t2;
  FILE *f;

  if (argc < 5)
  {
    fprintf(stderr, "PathTracer: Can not parse input parameters\n");
    fprintf(stderr, "USAGE: PathTracer size rec_depth num_samples output_name\n");
    fprintf(stderr, "   size = Image size (both along x and y)\n");
    fprintf(stderr, "   rec_depth = Recursion depth\n");
    fprintf(stderr, "   num_samples = Number of samples per pixel\n");
    fprintf(stderr, "   output_name = Name of the output file, e.g. MyRender.ppm\n");
    exit(0);
  }
  sx = atoi(argv[1]);
  MAX_DEPTH = atoi(argv[2]);
  num_samples = atoi(argv[3]);
  strcpy(&output_name[0], argv[4]);

  fprintf(stderr, "Rendering image at %d x %d\n", sx, sx);
  fprintf(stderr, "Recursion depth = %d\n", MAX_DEPTH);
  fprintf(stderr, "NUmber of samples = %d\n", num_samples);
  fprintf(stderr, "Output file name: %s\n", output_name);

  object_list = NULL;
  texture_list = NULL;

  // Allocate memory for the new image
  im = newImage(sx, sx);
  wght = (double *)calloc(sx * sx, sizeof(double));
  if (!im || !wght)
  {
    fprintf(stderr, "Unable to allocate memory for image\n");
    exit(0);
  }
  else
    rgbIm = (double *)im->rgbdata;
  for (i = 0; i < sx * sx; i++)
    *(wght + i) = 1.0;

  buildScene(); // Create a scene.

  // Mind the homogeneous coordinate w of all vectors below. DO NOT
  // forget to set it to 1, or you'll get junk out of the
  // geometric transformations later on.

  // Camera center
  e.px = 0;
  e.py = 0;
  e.pz = -15;
  e.pw = 1;

  // To define the gaze vector, we choose a point 'pc' in the scene that
  // the camera is looking at, and do the vector subtraction pc-e.
  // Here we set up the camera to be looking at the origin.
  g.px = 0 - e.px;
  g.py = 0 - e.py;
  g.pz = 0 - e.pz;
  g.pw = 1;
  // In this case, the camera is looking along the world Z axis, so
  // vector w should end up being [0, 0, -1]

  // Define the 'up' vector to be the Y axis
  up.px = 0;
  up.py = 1;
  up.pz = 0;
  up.pw = 1;

  // Set up view with given the above vectors, a 4x4 window,
  // and a focal length of -1 (why? where is the image plane?)
  // Note that the top-left corner of the window is at (-2, 2)
  // in camera coordinates.
  cam = setupView(&e, &g, &up, -3, -2, 2, 4);

  if (cam == NULL)
  {
    fprintf(stderr, "Unable to set up the view and camera parameters. Our of memory!\n");
    cleanup(object_list, texture_list);
    deleteImage(im);
    exit(0);
  }

  du = cam->wsize / (sx - 1);  // du and dv. In the notes in terms of wl and wr, wt and wb,
  dv = -cam->wsize / (sx - 1); // here we use wl, wt, and wsize. du=dv since the image is
  // and dv is negative since y increases downward in pixel
  // coordinates and upward in camera coordinates.

  fprintf(stderr, "View parameters:\n");
  fprintf(stderr, "Left=%f, Top=%f, Width=%f, f=%f\n", cam->wl, cam->wt, cam->wsize, cam->f);
  fprintf(stderr, "Camera to world conversion matrix (make sure it makes sense!):\n");
  printmatrix(cam->C2W);
  fprintf(stderr, "World to camera conversion matrix:\n");
  printmatrix(cam->W2C);
  fprintf(stderr, "\n");

  // Update light source weights - will give you weights for each light source that add up to 1
  obj = object_list;
  pct = 0;
  while (obj != NULL)
  {
    if (obj->isLightSource)
      pct += obj->LSweight;
    obj = obj->next;
  }
  obj = object_list;
  while (obj != NULL)
  {
    if (obj->isLightSource)
    {
      obj->LSweight /= pct;
    }
    obj = obj->next;
  }
  fprintf(stderr, "\n");

  NUM_RAYS = 0;

  t1 = time(NULL);

  fprintf(stderr, "Rendering pass... ");
  for (k = 0; k < num_samples; k++)
  {
    fprintf(stderr, "%d/%d, ", k, num_samples);
#pragma omp parallel for schedule(dynamic, 1) private(i, j, pc, wt, ray, col, d)
    for (j = 0; j < sx; j++) // For each of the pixels in the image
    {
      for (i = 0; i < sx; i++)
      {
        // Random sample within the pixel's area
        pc.px = (cam->wl + ((i + (drand48() - .5)) * du));
        pc.py = (cam->wt + ((j + (drand48() - .5)) * dv));
        pc.pz = cam->f;
        pc.pw = 1;

        // Convert image plane sample coordinates to world coordinates
        matVecMult(cam->C2W, &pc);

        // Now compute the ray direction
        memcpy(&d, &pc, sizeof(struct point3D));
        subVectors(&cam->e, &d); // Direction is d=pc-e
        normalize(&d);

        // Create a ray and do the raytracing for this pixel.
        ray = newRay(&pc, &d);

        if (ray != NULL)
        {
          wt = *(wght + i + (j * sx));
          // col.R = 0;
          // col.G = 0;
          // col.B = 0;
          PathTrace(ray, 1, &col, NULL, 1);
          (*(rgbIm + ((i + (j * sx)) * 3) + 0)) += col.R * pow(2, -log(wt));
          (*(rgbIm + ((i + (j * sx)) * 3) + 1)) += col.G * pow(2, -log(wt));
          (*(rgbIm + ((i + (j * sx)) * 3) + 2)) += col.B * pow(2, -log(wt));
          wt += col.R;
          wt += col.G;
          wt += col.B;
          *(wght + i + (j * sx)) = wt;
          free(ray);
        }
      } // end for i
    }   // end for j
    if (k % 25 == 0)
      dataOutput(rgbIm, sx, &output_name[0]); // Update output image every 25 passes
  }                                           // End for k
  t2 = time(NULL);

  fprintf(stderr, "\nDone!\n");

  dataOutput(rgbIm, sx, &output_name[0]);

  fprintf(stderr, "Total number of rays created: %ld\n", NUM_RAYS);
  fprintf(stderr, "Rays per second: %f\n", (double)NUM_RAYS / (double)difftime(t2, t1));

  // Exit section. Clean up and return.
  cleanup(object_list, texture_list); // Object and texture lists
  deleteImage(im);                    // Rendered image
  free(cam);                          // camera view
  free(wght);
  exit(0);
}
