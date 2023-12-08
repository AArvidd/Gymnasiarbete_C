#define STBI_MSC_SECURE_CRT
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "string.h"
#include "stdio.h"

#define LENGTH 10000

#define PWIDTH 200
#define PHEIGHT 200
#define CHANNELS 3 // RGB channels

#define CX 0
#define CY 0
#define CZ 0

static int INT[LENGTH];
static float FLOAT[LENGTH];
static int objects = 100;


struct ray{
    float x, y, z;
    float a, b, c;
    int depth;
};

struct sphere{
    float x, y, z;
    float R;
    float reflectivity;
    int color[3];
};

struct ground{
    float z;
    float reflectivity;
    float r1, g1, b1;
    float r2, g2, b2;
    int color[3];
};

struct camera{
    float x, y, z;
    int r, g, b;
};

void castRay(int *out, struct ray ray, struct sphere scene[], struct sphere *fScene, struct ground ground, struct camera camera);
void findIntersection(float* arr, struct ray ray, struct sphere scene[], struct sphere* fScene, struct ground ground);
float findSphereIntersection(struct sphere sphere, struct ray ray);
float findGroundIntersection(struct ground ground, struct ray ray);

void sphereReflect(float *out, float X, float Y, float Z, float a, float b, float c, struct sphere sphere);


void groundReflect(float *out, float a, float b, float c);



int main() {



    //struct ray rayTest = {1, 2, 3, 4, 5, 6, 7};
    //struct sphere sphereTest = {0.1f, 0.2f, 0.3f, 0.4f, 0.5f, 6, 7, 8};
    //struct ground groundTest = {1, 2, 3, 4, 5, 6, 7, 8};



    float direction[PWIDTH * PHEIGHT * 3];

    // Create a blank image
    unsigned char* image = (unsigned char*)malloc(PWIDTH * PHEIGHT * CHANNELS);
    memset(image, 255, PWIDTH * PHEIGHT * PHEIGHT);  // Set all pixels to white

    // Set every other pixel in the 50th row to red
    for (int x = 0; x < PWIDTH; x++) {
        if (x % 2 == 0) {
            image[3 * (x + 49 * PWIDTH) + 0] = 255;  // Red channel
            image[3 * (x + 49 * PWIDTH) + 1] = 0;    // Green channel
            image[3 * (x + 49 * PWIDTH) + 2] = 0;    // Blue channel
        }
    }

    // Save the image as a PNG file
    stbi_write_png("output.png", PWIDTH, PHEIGHT, CHANNELS, image, PWIDTH * CHANNELS);

    // Clean up
    free(image);


    return 0;

}

void castRay(int *out, struct ray ray, struct sphere scene[], struct sphere *fScene, struct ground ground, struct camera camera) {

    float intersection[2];
    findIntersection(intersection, ray, scene, fScene, ground);


    if (intersection[1] == -2) {
        *(out+0) = camera.r;
        *(out+1) = camera.g;
        *(out+2) = camera.b;
        return;
    }


    float x = ray.x + ray.a * intersection[0];
    float y = ray.y + ray.b * intersection[0];
    float z = ray.z + ray.c * intersection[0];


    int direct[3];
    float reflectDir[3];
    float reflectivity;

    if (intersection[1] == -1) {
        direct[0] = ground.color[0];
        direct[1] = ground.color[1];
        direct[2] = ground.color[2];
        groundReflect(reflectDir, ray.a, ray.b, ray.c);
        reflectivity = ground.reflectivity;
    } else {
        struct sphere hit = scene[(int)intersection[1]];
        direct[0] = hit.color[0];
        direct[1] = hit.color[1];
        direct[2] = hit.color[2];
        sphereReflect(reflectDir, x, y, z, ray.a, ray.b, ray.c, hit);
        reflectivity = hit.reflectivity;
    }

    if (ray.depth == 0 || reflectivity == 0) {
        *(out + 0) = direct[0];
        *(out + 1) = direct[1];
        *(out + 2) = direct[2];
        return;
    }


    struct ray next = {x, y, z,  reflectDir[0], reflectDir[1], reflectDir[2], ray.depth - 1};
    int reflect[3];

    castRay(reflect, next, scene, fScene, ground, camera);

    *(out+0) = (int) (direct[0] * (1 - reflectivity) + reflect[0] * reflectivity);
    *(out+1) = (int) (direct[1] * (1 - reflectivity) + reflect[1] * reflectivity);
    *(out+2) = (int) (direct[2] * (1 - reflectivity) + reflect[2] * reflectivity);

}

void findIntersection(float* arr, struct ray ray, struct sphere scene[], struct sphere* fScene, struct ground ground) {
    float intersectionLength[objects];


    for (int i = 0; i < scene + objects - fScene; i++) {
        struct sphere current = scene[i];
        intersectionLength[i] = findSphereIntersection(current, ray);
    }

    float shortest = intersectionLength[0];
    int shortestIndex = 0;

    for (int i = 0; i < objects; i++) {
        if (intersectionLength[i] == -1)
            continue;
        shortest = intersectionLength[i];
        shortestIndex = i;
        break;
    }

    for (int i = 0; i < objects; i++) {
        float test = intersectionLength[i];
        if (test == -1) {
            continue;
        }
        if (test < shortest) {
            shortest = test;
            shortestIndex = i;
        }
    }

    if (shortest == -1) {
        shortest = findGroundIntersection(ground, ray);
        shortestIndex = -1;
    }


    if (shortest == -1) {
        shortest = -2;
        shortestIndex = -2;
    }

    *arr = shortest;
    *(arr+1) = (float)shortestIndex;

}

float findSphereIntersection(struct sphere sphere, struct ray ray) {
    float t1;
    float t2;

    float A = ray.a, B = ray.b, C = ray.c;
    float X = ray.x, Y = ray.y, Z = ray.z;
    float h = sphere.x, k = sphere.y, l = sphere.z, r = sphere.R;

    float a = (A * A + B * B + C * C);
    float b = 2 * (A * (X - h) + B * (Y - k) + C * (Z - l));
    float c = -(r * r - (X * X + Y * Y + Z * Z - 2 * (X * h + Y * k + Z * l) + h * h + k * k + l * l));

    float center = -(b / (2 * a));
    float pm = sqrtf((b * b) - (4 * a * c)) / (2 * a);

    t1 = (center + pm);
    t2 = (center - pm);

    if (t1 == NAN && t2 == NAN)
        return -1;

    if (t1 < 0 || t2 < 0)
        return -1;

    float t = t1;

    if (t2 < t1)
        t = t2;

    return t;
}

float findGroundIntersection(struct ground ground, struct ray ray) {
    float t = (ground.z - ray.z) / ray.c;
    if (t < 0)
        t = -1;
    return t;
}


void sphereReflect(float *out, float X, float Y, float Z, float a, float b, float c, struct sphere sphere) {

    float h, k, l;

    h = X - sphere.x;
    k = Y - sphere.y;
    l = Z - sphere.z;

    float abs = sqrtf((h * h) + (k * k) + (l * l));

    h = h / abs;
    k = k / abs;
    l = l / abs;

    float dot = (a * h) + (b * k) + (c * l);
    *(out+0) = a - (2 * dot * h);
    *(out+1) = b - (2 * dot * k);
    *(out+2) = c - (2 * dot * l);

}


void groundReflect(float *out, float a, float b, float c) {

    float h, k, l;
    h = k = 0;
    l = 1;

    float dot = (a * h) + (b * k) + (c * l);
    *(out+0) = a - (2 * dot * h);
    *(out+1) = b - (2 * dot * k);
    *(out+2) = c - (2 * dot * l);


}

