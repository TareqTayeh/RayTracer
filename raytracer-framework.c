/* Name: Tareq Tayeh
 * Student Number: 250725776
 * Course: CS3388A Computer Graphics
 * Assignment Number: Assignment #4
 * Date: 30 November 2017
 * Program purpose: The purpose of this program is to non-recursively ray-trace images of scenes containing a number of simple generic objects
 					and must be able to must be able to:
					- Render spheres and planes as generic objects
					- Implement shading (ambient, diffuse, and specular reflections)
					- Implement shadowing through the use of shadow rays
 * File: raytracer-framework.c
 */

/*            PURPOSE : Simple framework for ray-tracing

        PREREQUISITES : matrix.h
 */

#include <X11/Xlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "matrix.h"

#define INFINITE_PLANE   0
#define PLANE            1
#define SPHERE           2

#define EPSILON          0.00001
#define INFTY            999999.9
#define N_OBJECTS        1024
#define MAX_INTENSITY    255.0

#define Ex               12.0
#define Ey               12.0
#define Ez               3.0

#define Gx               0.0
#define Gy               0.0
#define Gz               0.0

#define UPx              0.0
#define UPy              0.0
#define UPz              1.0

#define Lx               0.0
#define Ly               0.0
#define Lz               10.0

#define Near             1.0
#define Far              25.0

#define THETA            45.0
#define ASPECT           1.5

#define H                300

#define M_PI		     3.14

typedef struct {
    int width, height ;
} window_t ;

typedef struct {
    dmatrix_t UP ;
    dmatrix_t E ;
    dmatrix_t G ;
    dmatrix_t u, v, n ;
} camera_t ;

typedef struct {
    double r, g, b ;
} color_t ;

typedef struct {
    int type ;
    double (*intersection_function)(dmatrix_t *,dmatrix_t *) ;
    dmatrix_t M, Minv ;
    color_t specular_color, diffuse_color, ambient_color ;
    double reflectivity, specular_coeff, diffuse_coeff, ambient_coeff, f ;
} object_t ;

typedef struct {
    dmatrix_t position ;
    color_t color ;
    color_t intensity ;
} light_t ;

object_t object[N_OBJECTS] ;
int nobjects = 0 ;

Display *InitX(Display *d, Window *w, int *s, window_t *Window) {

    d = XOpenDisplay(NULL) ;
    if(d == NULL) {
        printf("Cannot open display\n") ;
        exit(1) ;
    }
    *s = DefaultScreen(d) ;
    *w = XCreateSimpleWindow(d,RootWindow(d,*s),0,0,Window->width,Window->height,1,BlackPixel(d,*s),WhitePixel(d, *s)) ;
    Atom delWindow = XInternAtom(d,"WM_DELETE_WINDOW",0) ;
    XSetWMProtocols(d,*w,&delWindow,1) ;
    XSelectInput(d,*w,ExposureMask | KeyPressMask) ;
    XMapWindow(d,*w) ;
    return(d) ;
}

void SetCurrentColorX(Display *d, GC *gc, unsigned int r, unsigned int g, unsigned int b) {

    XSetForeground(d,*gc,r << 16 | g << 8 | b) ;
}

void SetPixelX(Display *d, Window w, int s, int i, int j) {

    XDrawPoint(d,w,DefaultGC(d,s),i,j) ;
}

void QuitX(Display *d, Window w) {

    XDestroyWindow(d,w) ;
    XCloseDisplay(d) ;
}

light_t *build_light(light_t *light, dmatrix_t *position, color_t color, color_t intensity) {

    dmat_alloc(&light->position,4,1) ;

    light->position = *position ;
    light->color.r = color.r ;
    light->color.g = color.g ;
    light->color.b = color.b ;
    light->intensity.r = intensity.r ;
    light->intensity.g = intensity.g ;
    light->intensity.b = intensity.b ;
    return light ;
}

window_t *build_window(window_t *Window, int height, float aspect) {

    Window->height = height ;
    Window->width =  aspect*height ;

    return(Window) ;
}

camera_t *build_camera(camera_t *Camera, window_t *Window) {

    dmat_alloc(&Camera->E,4,1) ;

    Camera->E.m[1][1] = Ex ;
    Camera->E.m[2][1] = Ey ;
    Camera->E.m[3][1] = Ez ;
    Camera->E.m[4][1] = 1.0 ;

    dmat_alloc(&Camera->G,4,1) ;

    Camera->G.m[1][1] = Gx ;
    Camera->G.m[2][1] = Gy ;
    Camera->G.m[3][1] = Gz ;
    Camera->G.m[4][1] = 1.0 ;

    dmat_alloc(&Camera->n,4,1) ;
    Camera->n = *dmat_normalize(dmat_sub(&Camera->E,&Camera->G)) ;
    Camera->n.l = 3 ;

    dmat_alloc(&Camera->UP,4,1) ;

    Camera->UP.l = 3 ;

    Camera->UP.m[1][1] = UPx ;
    Camera->UP.m[2][1] = UPy ;
    Camera->UP.m[3][1] = UPz ;
    Camera->UP.m[4][1] = 1.0 ;

    dmat_alloc(&Camera->u,4,1) ;

    Camera->u = *dmat_normalize(dcross_product(&Camera->UP,&Camera->n)) ;
    Camera->v = *dmat_normalize(dcross_product(&Camera->n,&Camera->u)) ;

    return(Camera) ;
}

dmatrix_t *intersection_coordinates(dmatrix_t *e, dmatrix_t *direction, double t) {

    dmatrix_t *intersection ;

    intersection = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
    dmat_alloc(intersection,4,1) ;

    intersection->m[1][1] = e->m[1][1] + direction->m[1][1]*t ;
    intersection->m[2][1] = e->m[2][1] + direction->m[2][1]*t ;
    intersection->m[3][1] = e->m[3][1] + direction->m[3][1]*t ;
    intersection->m[4][1] = 1.0 ;

    return intersection ;
}

double infinite_plane_intersection(dmatrix_t *e, dmatrix_t *d) {

    if (d->m[3][1] >= 0.0) return -1.0 ; else return -1.0*e->m[3][1]/d->m[3][1] ;
}

double plane_intersection(dmatrix_t *e, dmatrix_t *d) {

    double t ;
    dmatrix_t *intersection ;

    if (d->m[3][1] >= 0.0) {
        return -1.0 ;
    }
    else {
        t = -1.0*e->m[3][1]/d->m[3][1] ;
        intersection = intersection_coordinates(e,d,t) ;
        if ((fabs(intersection->m[1][1]) < 1.0) && (fabs(intersection->m[2][1]) < 1.0)) {
            delete_dmatrix(intersection) ;
            return t ;
        }
        else {
            delete_dmatrix(intersection) ;
            return -1.0 ;
        }
    }
}

double sphere_intersection(dmatrix_t *e, dmatrix_t *d) {

    double a = ddot_product(from_homogeneous(d),from_homogeneous(d)) ;
    double b = ddot_product(from_homogeneous(e),from_homogeneous(d)) ;
    double c = ddot_product(from_homogeneous(e),from_homogeneous(e)) - 1.0 ;

    double discriminant = b*b - a*c ;

    if (discriminant < 0.0) {
        return -1.0 ;
    }
    else {
        if (discriminant < EPSILON) {
            return -b/a ;
        }
        else {
            double t1 = -b/a - sqrtl(discriminant)/a ;
            double t2 = -b/a + sqrtl(discriminant)/a ;
            if (t1 < t2) {
                if (t1 > EPSILON) {
                    return t1 ;
                }
                else {
                    return -1.0 ;
                }
            }
            else {
                return t2 ;
            }
        }
    }
}

dmatrix_t *normal_to_surface(object_t *object, dmatrix_t *intersection) {

    dmatrix_t *normal ;

    normal = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
    dmat_alloc(normal,4,1) ;

    switch ((*object).type) {

        case PLANE          :   normal->m[1][1] = 0.0 ;
                                normal->m[2][1] = 0.0 ;
                                normal->m[3][1] = 1.0 ;
                                break ;

        case INFINITE_PLANE :   normal->m[1][1] = 0.0 ;
                                normal->m[2][1] = 0.0 ;
                                normal->m[3][1] = 1.0 ;
                                normal->m[4][1] = 0.0 ;
                                break ;

        case SPHERE         :   normal->m[1][1] = intersection->m[1][1] ;
                                normal->m[2][1] = intersection->m[2][1] ;
                                normal->m[3][1] = intersection->m[3][1] ;
                                normal->m[4][1] = 0.0 ;
                                break ;

        default: printf("No such object type\n") ;

    }
    return normal ;
}

int find_min_hit_time(double t0[N_OBJECTS]) {

    double min_t = INFTY ;
    int position = -1 ;

    for(int i = 0 ; i < nobjects ; i++) {
        if (t0[i] != -1.0) {
            if (t0[i] < min_t) {
                min_t = t0[i] ;
                position = i ;
            }
        }
    }
    return position ;
}

dmatrix_t *ray_direction(camera_t *Camera, window_t *Window, double height, double width, double i, double j) {

    dmatrix_t *d ;

    d = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
    dmat_alloc(d,4,1) ;

    d->m[1][1] = -1.0*Near*Camera->n.m[1][1] +
    width*(2.0*i/Window->width - 1.0)*Camera->u.m[1][1] +
    height*(2.0*j/Window->height - 1.0)*Camera->v.m[1][1] ;

    d->m[2][1] = -1.0*Near*Camera->n.m[2][1] +
    width*(2.0*i/Window->width - 1.0)*Camera->u.m[2][1] +
    height*(2.0*j/Window->height - 1.0)*Camera->v.m[2][1] ;

    d->m[3][1] = -1.0*Near*Camera->n.m[3][1] +
    width*(2.0*i/Window->width - 1.0)*Camera->u.m[3][1] +
    height*(2.0*j/Window->height - 1.0)*Camera->v.m[3][1] ;

    d->m[4][1] = 0.0 ;

    return(d) ;
}

dmatrix_t *vector_to_light_source(dmatrix_t *intersection, dmatrix_t *light_position) {

    dmatrix_t *s, *sn ;

    s = dmat_sub(light_position,intersection) ;
    sn = dmat_normalize(s) ;
    delete_dmatrix(s) ;

    return sn ;
}

dmatrix_t *vector_to_center_of_projection(dmatrix_t *intersection, dmatrix_t *e) {

    dmatrix_t *v, *vn ;

    v = dmat_sub(e,intersection) ;
    vn = dmat_normalize(v) ;
    delete_dmatrix(v) ;

    return vn ;
}

dmatrix_t *vector_to_specular_reflection(dmatrix_t *N, dmatrix_t *S) {

    dmatrix_t *r, *rn ;

    r = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
    dmat_alloc(r,4,1) ;

    double sn = 2.0*ddot_product(N,S) ;

    r->m[1][1] = -1.0*S->m[1][1] + sn*N->m[1][1] ;
    r->m[2][1] = -1.0*S->m[2][1] + sn*N->m[2][1] ;
    r->m[3][1] = -1.0*S->m[3][1] + sn*N->m[3][1] ;
    r->m[4][1] = 0.0 ;

    rn = dmat_normalize(r) ;
    delete_dmatrix(r) ;

    return rn ;
}

color_t color_init(double r, double g, double b) {

    color_t s ;

    s.r = r ;
    s.g = g ;
    s.b = b ;

    return s ;
}

color_t color_mult(double a, color_t c) {

    color_t s ;

    s.r = a*c.r ;
    s.g = a*c.g ;
    s.b = a*c.b ;

    return s ;
}

color_t color_add(color_t c1, color_t c2) {

    color_t s ;

    s.r = c1.r + c2.r ;
    s.g = c1.g + c2.g ;
    s.b = c1.b + c2.b ;

    return s ;
}

/* Author: Tareq Tayeh
 * Date of creation: 28 November 2017
 *
 * Function:shade()
 * Parameter:
 *  In: light *light - The light object
 *      object_t *object - The physical objects to be put on the screen
 *		dmatrix_t *e - Eye Position
 *		dmatrix_t *d - Direction of Ray
 *		color_t color - Color object
 *		color_t background - Background color object
 *  Out: None
 * Returns: color_t - The color of the pixel to be shaded
 * Desc: Computes the color of the pixel to be shaded on the window, and eventually would be used to shade entire window
 */

color_t shade(light_t *light, object_t *object, dmatrix_t *e, dmatrix_t *d, color_t color, color_t background) {

    int h,k;
    double t0[N_OBJECTS], Id, Is;

    dmatrix_t *L, *S, *V, *R, *I, *N;
    dmatrix_t *Te, *Td, *Tdn;
    dmatrix_t *Ts, *Tl, *Ti;

    //Color set to background color
    color.r = MAX_INTENSITY*background.r ;
    color.g = MAX_INTENSITY*background.g ;
    color.b = MAX_INTENSITY*background.b ;

    for(k = 0; k < nobjects; k++){
        Te = dmat_mult(&object[k].Minv,e); //Transformed eye position
        Td = dmat_mult(&object[k].Minv,d); //Transformed direction of the ray
        Tdn = dmat_normalize(Td); //Normalized Td
        t0[k] = ((object[k].intersection_function))(Te,Tdn); //Find the intersection between Te, Td, and the object giving a float value of "t"

        delete_dmatrix(Te);
        delete_dmatrix(Td);
        delete_dmatrix(Tdn);
    }

    h = find_min_hit_time(t0); //The time of the first hit

    //If there is a hit
    if(h != -1){
        Ts = dmat_mult(&object[h].Minv, &light->position); //Transformed Light Source
        Te = dmat_mult(&object[h].Minv,e);
        Td = dmat_mult(&object[h].Minv,d);
        Tdn = dmat_normalize(Td);

        I = intersection_coordinates(Te,Tdn,t0[h]); //e + dt

        N = normal_to_surface(&object[h],I);
        V = vector_to_center_of_projection(I,Te);
        S = vector_to_light_source(I,Ts);
        R = vector_to_specular_reflection(I,S);

        //Diffuse Coef: Id = N dot product S
        Id = ddot_product(N,S);//(N->m[1][1] * S->m[1][1]) +  (N->m[2][1] * S->m[2][1]) + (N->m[3][1] * S->m[3][1]);
      	if (Id < 0 ){ //All shadows are at 0
      		Id = 0;
      	}
        //Specular Coef & Rate of Decay: Is = R dot product V to the power of f(decay coef)
        Is = pow(ddot_product(R,V),object[h].f);

        //Red Channel color: Multiplying the color and intensity with the diffuse, specular and diffuse coefficients and colors
        color.r = MAX_INTENSITY * light->color.r * light->intensity.r * (Id * object[h].diffuse_coeff * object[h].diffuse_color.r
                  + object[h].ambient_coeff * object[h].ambient_color.r + Is * object[h].specular_coeff * object[h].specular_color.r );
        //Blue Channel color: Multiplying the color and intensity with the diffuse, specular and diffuse coefficients and colors
        color.b = MAX_INTENSITY * light->color.b * light->intensity.b * (Id * object[h].diffuse_coeff * object[h].diffuse_color.b
                  + object[h].ambient_coeff * object[h].ambient_color.b + Is * object[h].specular_coeff * object[h].specular_color.b );
        //Green Channel color: Multiplying the color and intensity with the diffuse, specular and diffuse coefficients and colors
        color.g = MAX_INTENSITY * light->color.g * light->intensity.g * (Id * object[h].diffuse_coeff * object[h].diffuse_color.g
                  + object[h].ambient_coeff * object[h].ambient_color.g + Is * object[h].specular_coeff * object[h].specular_color.g );

        //L is the vector of reflectivity, Ln is L normalized
        dmatrix_t *Ln;
        L = (dmatrix_t *)malloc(sizeof(dmatrix_t)) ;
        dmat_alloc(L,4,1) ;

        L->m[1][1] = -1.0*V->m[1][1] + N->m[1][1] ;
        L->m[2][1] = -1.0*V->m[2][1] + N->m[2][1] ;
        L->m[3][1] = -1.0*V->m[3][1] + N->m[3][1] ;
        L->m[4][1] = 0.0 ;

        Ln = dmat_normalize(L);

        //Ti: Transformed Back Intersection, Tl: Transformed Light Reflection
        Ti = dmat_mult(&object[h].M, I);
        Tl = dmat_mult(&object[h].M, L);
        //Preventing self intersection
        Ti->m[1][1] = Ti->m[1][1] + EPSILON * Tl->m[1][1];
        Ti->m[2][1] = Ti->m[2][1] + EPSILON * Tl->m[2][1];
        Ti->m[3][1] = Ti->m[3][1] + EPSILON * Tl->m[3][1];

        //Set Color for Reflections
        color = color_add(color_mult(1-object[h].reflectivity,color),
                          color_mult(object[h].reflectivity,shade(light,object,Ti,Tl,color,background)));

      	delete_dmatrix(S);
      	delete_dmatrix(V);
        delete_dmatrix(R);
      	delete_dmatrix(I);
      	delete_dmatrix(N);
      	delete_dmatrix(Ts);
        delete_dmatrix(L);
        delete_dmatrix(Ln);
        delete_dmatrix(Tl);
        delete_dmatrix(Ti);
    }

    return color ;
}

object_t *build_object(int object_type, dmatrix_t *M, color_t ambient_color, color_t diffuse_color, color_t specular_color, double ambient_coeff, double diffuse_coeff, double specular_coeff, double f, double reflectivity) {

    object_t *object ;

    object = malloc(sizeof(*object));

    object->type = object_type ;
    object->M = *M ;
    dmat_alloc(&object->Minv,4,4) ;
    object->Minv = *dmat_inverse(&object->M) ;

    object->reflectivity = reflectivity ;

    object->specular_color.r = specular_color.r ;
    object->specular_color.g = specular_color.g ;
    object->specular_color.b = specular_color.b ;
    object->specular_coeff = specular_coeff ;
    object->f = f ;

    object->diffuse_color.r = diffuse_color.r ;
    object->diffuse_color.g = diffuse_color.g ;
    object->diffuse_color.b = diffuse_color.b ;
    object->diffuse_coeff = diffuse_coeff ;

    object->ambient_color.r = ambient_color.r ;
    object->ambient_color.g = ambient_color.g ;
    object->ambient_color.b = ambient_color.b ;
    object->ambient_coeff = ambient_coeff ;

    switch (object_type) {

        case SPHERE         :   object->intersection_function = &sphere_intersection ;
                                break ;

        case PLANE          :   object->intersection_function = &plane_intersection ;
                                break ;

        case INFINITE_PLANE :   object->intersection_function = &infinite_plane_intersection ;
                                break ;
    }

    nobjects++ ;
    return(object) ;
}

/* Author: Tareq Tayeh
 * Date of creation: 29 November 2017
 *
 * Function: CreateSphere()
 * Parameter:
 *  In: double ambient_color_red - Value of the red ambient component
 *  	double ambient_color_green - Value of the green ambient component
 *		double ambient_color_blue - Value of the blue ambient component
 *		double diffuse_color_red - Value of the red diffuse component
 *		double diffuse_color_green - Value of the green diffuse component
 *		double diffuse_color_blue - Value of the blue diffuse component
 *		double specular_color_red - Value of the red specular component
 *		double specular_color_green - Value of the green specular component
 *		double specular_color_blue - Value of the blue specular component
 *		double ambient_coeff - Value of the ambient coefficiant (Between 0 and 1)
 *		double diffuse_coeff - Value of the diffuse coefficiant (Between 0 and 1)
 *		double specular_coeff - Value of the specular coefficiant (Between 0 and 1)
 *		double f - To get the ray of decay
 *		double reflectivity - Value of the reflectivity component
 *		double xposition - X Position of the Sphere
 *		double yposition - Y Position of the Sphere
 *		double scale - Scale of the Sphere
 *  Out: None
 * Returns: None
 * Desc: Creates a sphere object that would then be built and stored in the physical objects array
 */
void CreateSphere(double ambient_color_red, double ambient_color_green, double ambient_color_blue,
                  double diffuse_color_red, double diffuse_color_green, double diffuse_color_blue,
                  double specular_color_red, double specular_color_green, double specular_color_blue,
                  double ambient_coeff, double diffuse_coeff, double specular_coeff, double f,
                  double reflectivity, double xposition, double yposition, double scale){

  color_t ambient_color, diffuse_color, specular_color;
  dmatrix_t M;
  dmat_alloc(&M,4,4) ;
  M = *dmat_identity(&M) ;

  M.m[1][4] = xposition ;
  M.m[2][4] = yposition ;
  M.m[3][3] = scale ;

  ambient_color.r = ambient_color_red;
  ambient_color.g = ambient_color_green;
  ambient_color.b = ambient_color_blue;
  diffuse_color.r = diffuse_color_red;
  diffuse_color.g = diffuse_color_green;
  diffuse_color.b = diffuse_color_blue;
  specular_color.r = specular_color_red;
  specular_color.g = specular_color_green;
  specular_color.b = specular_color_blue;


  object[nobjects] = *build_object(SPHERE,&M,ambient_color,diffuse_color,specular_color,ambient_coeff,diffuse_coeff,specular_coeff,f,reflectivity) ;
}

int main() {

    Display *d ;
    Window w ;
    XEvent e ;

    int i, j, s ;

    camera_t Camera ;
    window_t Window ;
    light_t light ;
    dmatrix_t M, light_position ;
    color_t pixel, background, light_intensity, light_color, ambient_color, diffuse_color, specular_color ;
    double height, width, aspect, ambient_coeff, diffuse_coeff, specular_coeff, f, reflectivity ;

    /* Set the background color */

    background.r = 0.0 ;
    background.g = 0.0 ;
    background.b = 0.0 ;

    /* Set up light position, intensity, and color */

    dmat_alloc(&light_position,4,1) ;

    light_position.m[1][1] = Lx ;
    light_position.m[2][1] = Ly ;
    light_position.m[3][1] = Lz ;
    light_position.m[4][1] = 1.0 ;

    light_intensity.r = 1.0 ;
    light_intensity.g = 1.0 ;
    light_intensity.b = 1.0 ;

    light_color.r = 1.0 ;
    light_color.g = 1.0 ;
    light_color.b = 1.0 ;

    light = *build_light(&light,&light_position,light_color,light_intensity) ;

    /* Build display window and synthetic camera */

    Window = *build_window(&Window,H,ASPECT) ;
    Camera = *build_camera(&Camera,&Window) ;

    /* Creating spheres */

    CreateSphere(1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.2, 0.4, 0.4, 10, 0.2, 0.0, 0.0, 2.0); //The Red Sphere
    CreateSphere(0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 0.3, 0.5, 0.2, 20, 0.3, 3.0, 6.0, 1.0); //The Blue Sphere
    CreateSphere(0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.5, 0.2, 0.3, 20, 0.3, -3.0, -6.0, 3.0); //The Green Sphere
    CreateSphere(1.0, 1.0, 0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.4, 0.2, 0.4, 20, 0.3, 0.0, -6.0, 4.0); //The Yellow Sphere
    CreateSphere(0.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.1, 0.3, 0.6, 30, 0.3, -3.0, 0.0, 4.0); //The Cyan Sphere

    /* Create tiled floor with planes */

    for (i = -4 ; i < 4 ; i++) {
        for (j = -4 ; j < 4 ; j++) {

            dmat_alloc(&M,4,4) ;
            M = *dmat_identity(&M) ;

            M.m[1][4] = (double)2*i ;
            M.m[2][4] = (double)2*j ;
            M.m[3][4] = -2.0 ;

            reflectivity = 0.5 ;

            specular_color.r = 1.0 ;
            specular_color.g = 1.0 ;
            specular_color.b = 1.0 ;
            specular_coeff = 0.0 ;
            f = 10.0 ;

            diffuse_color.r = 1.0 ;
            diffuse_color.g = 1.0 ;
            diffuse_color.b = 1.0 ;
            diffuse_coeff = 0.0 ;

            ambient_color.r = (double)(abs(i+j)%2) ;
            ambient_color.g = (double)(abs(i+j)%2) ;
            ambient_color.b = (double)(abs(i+j)%2) ;
            ambient_coeff = 1.0 ;

            object[nobjects] = *build_object(PLANE,&M,ambient_color,diffuse_color,specular_color,ambient_coeff,diffuse_coeff,specular_coeff,f,reflectivity) ;
        }
    }

    /* Set near plane dimensions */

    aspect = ASPECT ;
    height = Near*tan(M_PI/180.0*THETA/2.0) ;
    width = height*aspect ;

    dmatrix_t *direction ;

    d = InitX(d,&w,&s,&Window) ;
    XNextEvent(d, &e) ;

    while (1) {
        XNextEvent(d,&e) ;
        if (e.type == Expose) {

            for (i = 0 ; i < Window.width ; i++) {
                for (j = 0  ; j < Window.height ; j++) {
                    direction = ray_direction(&Camera,&Window,height,width,(double)i,(double)j) ;
                    pixel = shade(&light,object,&Camera.E,direction,pixel,background) ;
                    SetCurrentColorX(d,&(DefaultGC(d,s)),(int)pixel.r,(int)pixel.g,(int)pixel.b) ;
                    SetPixelX(d,w,s,i,Window.height - (j + 1)) ;
                    delete_dmatrix(direction) ;
                }
            }
        }
        if (e.type == KeyPress)
            break ;
        if (e.type == ClientMessage)
            break ;
    }
    QuitX(d,w) ;
}
