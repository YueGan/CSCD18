/*
  CSC D18 - Assignment 1 - 2D light propagation

  This is the place where you will be doing most of your work for solving this
  assignment. Before you start working here, you shold have become familiar
  with the data structures and functions available to you from light2D.h, and
  with the way a scene is built in buildScene.c

  Go over each of the functions below, and implement the different components
  of the solution in the sections marked

  /************************
  / TO DO:
  ************************ /

  Do not add or modify code outside these sections.

  Details about what needs to be implemented are described in the comments, as
  well as in the assignment handout. You must read both carefully.

  Starter by: F.J. Estrada, Aug. 2017
*/

/****************************************************************************
 * Uncomment the #define below to enable debug code, add whatever you need
 * to help you debug your program between #ifdef - #endif blocks
 * ************************************************************************/
#define __DEBUG_MODE

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

struct ray2D makeLightSourceRay(void)
{
 /*
   This function should return a light ray that has its origin at the light
   source, and whose direction depends on the type of light source.

   For point light sources (which emit light in all directions) the direction
    has to be chosen randomly and uniformly over a unit circle (i.e. any
    direction is equally likely)

   For a laser type light source, the direction of the ray is exactly the same
    as the direction of the lightsource.

   Set the colour of the ray to the same colour as the light source, and
    set the inside_outside flag to 0 (we assume light sources are
    outside objects)

   In either case, the direction vector *must be unit length*.
*/

 /************************************************************************
 *  TO DO: Complete this function so that we can sample rays from the
 *         lightsource for the light propagation process.
 ************************************************************************/

 struct ray2D ray;

 // This creates a dummy ray (which won't go anywhere since the direction
 // vector d=[0 0]!. But it's here so you see which data values you have to
 // provide values for given the light source position, and type.
 // ** REPLACE THE CODE BELOW with code that provides a valid ray that
 //    is consistent with the lightsource.

 ray.p.px=lightsource.l.p.px;			// Ray's origin
 ray.p.py=lightsource.l.p.py;
 ray.d.px=lightsource.l.d.px;			// Ray's direction
 ray.d.py=lightsource.l.d.py;
 ray.inside_out=0;		// Initially 0 since the ray starts outside an object
 ray.monochromatic=0;		// Initially 0 since the ray is white (from lightsource)
 ray.R=lightsource.R;			// Ray colour in RGB must be the same as the lightsource
 ray.G=lightsource.G;
 ray.B=lightsource.B;


 return(ray);			// Currently this returns dummy ray
}

void propagateRay(struct ray2D *ray, int depth)
{
 /*
   This function carries out the light propagation process. It is provided with access
   to a ray data structure, and must perform the following steps (in order!):

   - Check if maximum recursion depth has been reached (in which case, it just returns)
   - Find the *closest* intersection between the ray and objects in the scene. This
     means you have to check against the 4 walls, and any circles added in buildScene,
     determine which intersection is closest, and obtain the intersection point, the
     normal at the intersection, and the lambda value at which the intersection happens.
   - Renders the ray onto the image from its starting point all the way up to the
     intersection point.
   - At the intersection, use the material properties to determine how the propagation
     process proceeds:
         * For mirror materials, compute the mirror reflection direction, create a ray
           along that direction whose origin is the intersection point, and propagate
           that ray
         * For scattering materials, choose a random direction within +- 90 degrees of
           the normal, create a ray with that direction with origin at the intersection
           point, and propagate that ray
         * For refracting materials you will need to propagate two rays - one in the
           mirror reflection direction (just like for reflecting materials), and
           another in the refraction direction. Propagate both rays one after the other!

   NOTE: You should only care about intersections for which lambda is POSITIVE (in front
         of the ray), and greater than zero (e.g. if the ray is reflected from some
         object, you do not care about the intersection with the object itself which will
         have a lambda of *very close* to zero)

   In every case, make sure the ray's direction vector has unit length. You will need to
   complete other functions as part of your work here.
*/

 /*********************************************************************************
  * TO DO: Complete this function to implement light propagation!
  ********************************************************************************/

 // Define your local variables here
 if (depth>=max_depth) return;	 	// Leave this be, it makes sure you don't
					// recurse forever
// printf("recursion depth: %d\n", depth);

 // Step 1 - Find *closest* intersection with the 4 walls (the written part of A1
 //          should help you figure out how to do that.

 double lambda = 0;
 struct point2D p;
 struct point2D n;
 int type;
 double r_idx;

 for (int i = 0; i <=3; i++){

   struct point2D normal;

   // normal = (-dy, dx)
   normal.px = -walls[i].w.d.py;
   normal.py = walls[i].w.d.px;

   struct point2D p1_p0;

   p1_p0.px = walls[i].w.p.px - ray->p.px;
   p1_p0.py = walls[i].w.p.py - ray->p.py;


   if(dot(&p1_p0, &normal) / dot(&ray->d, &normal) > 0){

        if(lambda == 0 || dot(&p1_p0, &normal) / dot(&ray->d, &normal) < lambda){

          lambda = dot(&p1_p0, &normal) / dot(&ray->d, &normal);


          p.px = ray->p.px + lambda * ray->d.px;
          p.py = ray->p.py + lambda * ray->d.py;

          n.px = normal.px;
          n.py = normal.py;

          type = walls[i].material_type;

          normalize(&n);
        }
  }

 }






 // How many walls can the ray intersect? how many walls can the ray intersect in the
 // forward direction?

 // Step 2 - Check for intersection against objects in the object array - you must
 //          complete the intersectRay() function, call it, and obtain the closest
 //          intersection (in the forward ray direction) with objects in the scene.
 //          Note that you must provide variables for intersectRay() to return
 //          the point of intersection, normal at intersection, lambda, material type,
 //          and refraction index for the closest object hit by the ray.




 intersectRay(ray, &p, &n, &lambda, &type, &r_idx);

 // printf("after intersction, p: (%f, %f), lambda:%f\n", p.px, p.py, lambda);

 // Step 3 - Check whether the closest intersection with objects is closer than the
 //          closest intersection with a wall. Choose whichever is closer.

 // Step 4 - Render the ray onto the image. Use renderRay(). Provide renderRay() with
 //          the origin of the ray, and the intersection point (it will then draw a
 //          ray from the origin to the intersection). You also need to provide the
 //          ray's colour.

 // printf("this px: %f, py: %f, nx: %f, ny: %f \n", p.px, p.py, n.px, n.py);

 // renderRay(&ray->p, &ray->d, ray->R, ray->G, ray->B); // Original ray
 renderRay(&ray->p, &p, ray->R, ray->G, ray->B); // p lambda origin color

 // printf("normal: (%f, %f)\n", n.px, n.py);

 struct point2D normalPoint;
 normalPoint.px = p.px + n.px;
 normalPoint.py = p.py + n.py;
 // renderRay(&p, &normalPoint, 0, 0, 1); // normal blue


 // Step 5 - Decide how to handle the ray's bounce at the intersection. You will have
 //          to provide code for 3 cases:
 //          If material type = 0, you have a mirror-reflecting object.
 //                                Create a ray in the mirror reflection direction,
 //                                with the same colour as the incoming ray, and
 //                                with origin at the intersection point.
 //                                Then call propagateRay() recursively to trace it.
 //          if material type = 1, you have a scattering surface.
 //                                Choose a random direction within +- 90 degrees
 //                                from the normal at the intersection. Create a
 //                                ray in this direction, with the same colour as
 //                                the incoming ray, and origin at the intersection,
 //                                then call propagateRay() recursively to trace it.
 //          if material type = 2, you have a refracting (transparent) material.
 // 				   Here you need to process two rays:
 //                                * First, determine how much of the incoming light is
 //                                  reflected and how much is transmitted, using
 //				     Schlick's approximation:
 // 					 R0 = ((n1-n2)/(n1+n2))^2
 // 					 R(theta)=R0+((1-R0)*(1-cos(theta))^5)
 //				     If the ray is travelling from air to the inside
 //                                  of an object, n1=1, n2=object's index of refraction.
 //                                  If the ray is travelling from inside an object
 //                                  back onto air, n1=object's index of refraction, n2=1
 //				     And 'theta' is the angle between the normal and the
 // 				     ray direction.
 //				     R(theta) gives the amount Rs of reflected light,
 //				     1.0-R(theta) gives the amount Rt of transmitted light.
 //                                * Now, make a ray in the mirror-reflection direction
 //				     (same as for material type 0), with the same colour
 //				     as the incoming ray, but with intensity modulated
 //				     by Rs. (e.g. if the incoming's colour is R,G,B,
 //                                  the reflected ray's colour will be R*Rs, G*Rs, B*Rs)
 //				     trace this ray.
 //				   * Make a ray in the refracted-ray direction. The
 //				     angle for the transmitted ray is given by Snell's law
 //				     n1*sin(theta1) = n2*sin(theta2). The colour of the
 //				     transmitted ray is the same as the incoming ray but
 //			             modulated by Rt. Trace this ray.
 //	That's it! you're done!

 double angel = 0;
 double tangentAngel = 0;
 struct ray2D reflection;
 struct point2D topRotation;
 struct point2D bottomRotation;
 struct point2D tangentEq;
 tangentEq.px = -n.py;
 tangentEq.py = n.px;

 struct point2D flipNormal;
 flipNormal.px = -n.px;
 flipNormal.py = -n.py;

 // renderRay(&ray->p, &n, 0.5, 1, 0.2); // normal at center yellow
 // renderRay(&ray->p, &tangentEq, 0.5, 1, 0.2); // normal at center yellow

 if(type == 0){

  //  printf("-------normal: (%f, %f)\n", flipNormal.px, flipNormal.py);

  //  printf("-------inside: (%f)\n", dot(&ray->d, &flipNormal) / (dot(&ray->d, &ray->d) * dot(&flipNormal, &flipNormal)));
    // printf("d^2: %f, n^2: %f, dot: %f", dot(&ray->d, &ray->d), dot(&flipNormal, &flipNormal), dot(&ray->d, &flipNormal));
   angel = acos(dot(&ray->d, &flipNormal) / (dot(&ray->d, &ray->d) * dot(&flipNormal, &flipNormal)));
   tangentAngel = acos(dot(&ray->d, &tangentEq) / (dot(&ray->d, &ray->d) * dot(&tangentEq, &tangentEq)));



  //  printf("angel: %f, tangent Angel: %f\n", angel, tangentAngel);

   if(tangentAngel > PI/2){
     angel = -angel;
   }


   topRotation.px = cos(angel);
   topRotation.py = -sin(angel);

   bottomRotation.px = sin(angel);
   bottomRotation.py = cos(angel);

   reflection.d.px = dot(&topRotation, &n);
   reflection.d.py = dot(&bottomRotation, &n);

   reflection.p.px = p.px;
   reflection.p.py = p.py;

   reflection.R = ray->R;
   reflection.G = ray->G;
   reflection.B = ray->B;

   propagateRay(&reflection, depth+1);

  //  printf("r: (%f, %f)\n", reflection.d.px, reflection.d.py);
 }

 if(type == 1){

  //  printf("-------normal: (%f, %f)\n", flipNormal.px, flipNormal.py);

  //  printf("-------inside: (%f)\n", dot(&ray->d, &flipNormal) / (dot(&ray->d, &ray->d) * dot(&flipNormal, &flipNormal)));
    // printf("d^2: %f, n^2: %f, dot: %f", dot(&ray->d, &ray->d), dot(&flipNormal, &flipNormal), dot(&ray->d, &flipNormal));
   angel = rand() / (RAND_MAX + 1.) * PI/2;
   tangentAngel = acos(dot(&ray->d, &tangentEq) / (dot(&ray->d, &ray->d) * dot(&tangentEq, &tangentEq)));



  //  printf("angel: %f, tangent Angel: %f\n", angel, tangentAngel);

   if(tangentAngel > PI/2){
     angel = -angel;
   }


   topRotation.px = cos(angel);
   topRotation.py = -sin(angel);

   bottomRotation.px = sin(angel);
   bottomRotation.py = cos(angel);

   reflection.d.px = dot(&topRotation, &n);
   reflection.d.py = dot(&bottomRotation, &n);

   reflection.p.px = p.px;
   reflection.p.py = p.py;

   reflection.R = ray->R;
   reflection.G = ray->G;
   reflection.B = ray->B;

   propagateRay(&reflection, depth+1);

  //  printf("r: (%f, %f)\n", reflection.d.px, reflection.d.py);
 }

 if(type == 2){

   angel = acos(dot(&ray->d, &flipNormal) / (dot(&ray->d, &ray->d) * dot(&flipNormal, &flipNormal)));
   tangentAngel = acos(dot(&ray->d, &tangentEq) / (dot(&ray->d, &ray->d) * dot(&tangentEq, &tangentEq)));



  //  printf("angel: %f, tangent Angel: %f\n", angel, tangentAngel);

   if(tangentAngel > PI/2){
     angel = -angel;
   }


   topRotation.px = cos(angel);
   topRotation.py = -sin(angel);

   bottomRotation.px = sin(angel);
   bottomRotation.py = cos(angel);

   reflection.d.px = dot(&topRotation, &n);
   reflection.d.py = dot(&bottomRotation, &n);

   reflection.p.px = p.px;
   reflection.p.py = p.py;

    double r_0;
    r_0 = ((1 - r_idx) / (1 + r_idx)) * ((1 - r_idx) / (1 + r_idx));
    double r_theta = r_0 + (pow((1-r_0)*(1-cos(angel)), 5));
    double Rt = 1 - r_theta;


    reflection.R = ray->R * Rt;
    reflection.G = ray->G * Rt;
    reflection.B = ray->B * Rt;

   propagateRay(&reflection, depth+1);

  //  //          if material type = 2, you have a refracting (transparent) material.
  //  // 				   Here you need to process two rays:
  //  //                                * First, determine how much of the incoming light is
  //  //                                  reflected and how much is transmitted, using
  //  //				     Schlick's approximation:
  //  // 					 R0 = ((n1-n2)/(n1+n2))^2
  //  // 					 R(theta)=R0+((1-R0)*(1-cos(theta))^5)
  //
  //  angel = acos(dot(&ray->d, &flipNormal) / (dot(&ray->d, &ray->d) * dot(&flipNormal, &flipNormal)));
  //  tangentAngel = acos(dot(&ray->d, &tangentEq) / (dot(&ray->d, &ray->d) * dot(&tangentEq, &tangentEq)));
  //
  //  double r_0;
  //
  //  //if from air to inside
  //  if(r_idx > 0){
  //
  //    r_0 = ((1 - r_idx) / (1 + r_idx)) * ((1 - r_idx) / (1 + r_idx));
  //
  //  }
  //  else{
  //
  //    r_0 = ((r_idx - 1) / (r_idx + 1)) * ((r_idx - 1) / (r_idx + 1));
  //
  //  }
  //  double r_theta = r_0 + (pow((1-r_0)*(1-cos(angel)), 5));
  //  double Rt = 1 - r_theta;
  //
  //
  //  //				     If the ray is travelling from air to the inside
  //  //                                  of an object, n1=1, n2=object's index of refraction.
  //  //                                  If the ray is travelling from inside an object
  //  //                                  back onto air, n1=object's index of refraction, n2=1
  //  //				     And 'theta' is the angle between the normal and the
  //  // 				     ray direction.
  //  //				     R(theta) gives the amount Rs of reflected light,
  //  //				     1.0-R(theta) gives the amount Rt of transmitted light.
  //  //                                * Now, make a ray in the mirror-reflection direction
  //  //				     (same as for material type 0), with the same colour
  //  //				     as the incoming ray, but with intensity modulated
  //  //				     by Rs. (e.g. if the incoming's colour is R,G,B,
  //  //                                  the reflected ray's colour will be R*Rs, G*Rs, B*Rs)
  //  //				     trace this ray.
  //  //				   * Make a ray in the refracted-ray direction. The
  //  //				     angle for the transmitted ray is given by Snell's law
  //  //				     n1*sin(theta1) = n2*sin(theta2). The colour of the
  //  //				     transmitted ray is the same as the incoming ray but
  //  //			             modulated by Rt. Trace this ray.
  //
  //  if(tangentAngel > PI/2){
  //    angel = -angel;
  //  }
  //
  //
  //  topRotation.px = cos(angel);
  //  topRotation.py = -sin(angel);
  //
  //  bottomRotation.px = sin(angel);
  //  bottomRotation.py = cos(angel);
  //
  //  reflection.d.px = dot(&topRotation, &n);
  //  reflection.d.py = dot(&bottomRotation, &n);
  //
  //  reflection.p.px = p.px;
  //  reflection.p.py = p.py;
  //
  //  reflection.R = ray->R * Rt;
  //  reflection.G = ray->G * Rt;
  //  reflection.B = ray->B * Rt;
  //
  //  propagateRay(&reflection, depth+1);
  // //  propagateRay(&refraction, depth+1);1

 }

 // renderRay(&ray->p, &reflection.d, 0, 1, 0); // normal at center green
 // renderRay(&p, &reflection.d, 0, 1, 0); // reflection green




}

void intersectRay(struct ray2D *ray, struct point2D *p, struct point2D *n, double *lambda, int *type, double *r_idx)
{
 /*
  This function checks for intersection between the ray and any objects in the objects
  array. The objects are circles, so we are in fact solving for the intersection
  between a ray and a circle.

  For a unit circle centered at the origin, we would have the equation

  x^2 + y^2 = 1

  Using vector notation, with C=[x y]', we get

  ||C||^2 = 1

  A point on the ray is given by p + lambda*d

  Substituting in the equation for the circle we have

  (p + lambda*d)(p + lambda*d) - 1 = 0

  If we expand the product above (here the product of two vectors is a DOT product),
  we can form a quadratic equation

  A*lambda^2 + B*lambda + C = 0

  Which as you know, has a very simple solution.

  Your task is to
  * Figure out A, B, and C, considering that your circles don't necessarily have r=1
  * Figure out how to deal with the fact that circles in the scene are likely
    *not* centered at the origin

  Then implement the code that will find the value(s) of lambda at the intersection(s).

  Note that you can have no intersections, 1 intersection, or 2 intersections

  This function *must* find the closest intersection (if any) and update the value
  of lambda, the intersection point p, the normal n at the intersection,
  the corresponding object's material type (needed outside here to figure out how
  to handle the light's bouncing off this object), and the index of refraction for
  the object (needed if this is a transparent object).

  You must make sure that n is a unit-length vector.
 */

 /**********************************************************************************
  * TO DO: Complete this function to find the closest intersection between the
  *        ray and any objects in the scene, as well as the values at the
  *        intersection that will be needed to determine how to bounce/refract the
  *	   ray.
  * *******************************************************************************/

  double inside = 0.0;
  double lambda1 = 0.0;
  double lambda2 = 0.0;

  // For each circle, calculate the lambda(if any)
  for (int i = 0; i <= MAX_OBJECTS; i++){

    if(objects[i].r > 0){

      // Move the points opposite direction of circle center
      struct point2D tmpRayP = ray->p;
      tmpRayP.px = tmpRayP.px - objects[i].c.px;
      tmpRayP.py = tmpRayP.py - objects[i].c.py;

      // Calculate inside first
      inside = dot(&tmpRayP, &ray->d) * dot(&tmpRayP, &ray->d) -
      dot(&ray->d, &ray->d) * (dot(&tmpRayP, &tmpRayP) - objects[i].r*objects[i].r);


      if(round(inside*10000000)/10000000 > 0.0){

        // If positive, find lambda1 and lambda2, replace lambda if found smaller positive one
        lambda1 = (-(dot(&tmpRayP,  &ray->d)) + sqrt(inside)) / dot(&ray->d, &ray->d);
        lambda2 = (-(dot(&tmpRayP,  &ray->d)) - sqrt(inside)) / dot(&ray->d, &ray->d);

        if((round(lambda1*10000000)/10000000 > 0.0) && (lambda1 < *lambda)){

          *lambda = lambda1;
          //p
          p->px = ray->p.px + lambda1 * ray->d.px;
          p->py = ray->p.py + lambda1 * ray->d.py;
          //n
          n->px = (p->px - objects[i].c.px);
          n->py = (p->py - objects[i].c.py);
          normalize(n);
          //material_type
          *type = objects[i].material_type;

        }

        if((round(lambda2*10000000)/10000000 > 0.0) && (lambda2 < *lambda)){

          *lambda = lambda2;

          //p
          p->px = ray->p.px + lambda2 * ray->d.px;
          p->py = ray->p.py + lambda2 * ray->d.py;
          //n
          n->px = (p->px - objects[i].c.px) / objects[i].r;
          n->py = (p->py - objects[i].c.py) / objects[i].r;
          normalize(n);
          //material_type
          *type = objects[i].material_type;

        }


       }
    }
  }
}
