#define STBI_MSC_SECURE_CRT
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include "string.h"
#include "time.h"


#define PWIDTH  5000
#define PHEIGHT 5000
#define CHANNELS 3 // RGB channels

#define WIDTH 2

static int objects = 0;


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
    int r1, g1, b1;
    int r2, g2, b2;
};

struct camera{
    float x, y, z;
    int r, g, b;
};

struct sphere * defineSphere( struct sphere *fScene, float x, float y, float z, float R, float reflectivity, int r, int g, int b);

void pixelPoint(float width, float *pixelGridX, float *pixelGridY, float *pixelGridZ);

void castRay(int *out, struct ray ray, struct sphere scene[], struct sphere *fScene, struct ground ground, struct camera camera);
void findIntersection(float* arr, struct ray ray, struct sphere scene[], struct sphere* fScene, struct ground ground);
float findSphereIntersection(struct sphere sphere, struct ray ray);
float findGroundIntersection(struct ground ground, struct ray ray);

void sphereReflect(float *out, float X, float Y, float Z, float a, float b, float c, struct sphere sphere);
void groundReflect(float *out, float a, float b, float c);

void getGroundColor(int *out, struct ground ground, float x, float y);

static float pixelGridX[PWIDTH * PHEIGHT];
static float pixelGridY[PWIDTH * PHEIGHT];
static float pixelGridZ[PWIDTH * PHEIGHT];



int main() {

    clock_t time;

    struct ground ground = {-20, 0.25f, 0, 0, 0, 255, 255, 255};

    struct camera main = {0, 0, 0, 0, 0, 0};

    struct sphere scene[100];
    struct sphere *fScene = scene;

    fScene = defineSphere(fScene, 30, 0, 0, 10, 0.5f, 255, 0, 255);
    fScene = defineSphere(fScene, 20, -20, 0, 5, 0.25f, 0, 0, 255);
    fScene = defineSphere(fScene, 50, 0, 20, 10, 0.5f, 255, 0, 0);
    fScene = defineSphere(fScene, 100, 100, 10, 40, 1, 255, 255, 255);
    fScene = defineSphere(fScene, 10, 0, 10, 5, 0.5f, 0, 255, 0);
    fScene = defineSphere(fScene, -10, 0, 0, 5, 0, 255, 255, 255);
    fScene = defineSphere(fScene, 0, 20, 0, 5, 0, 0, 255, 255);
    fScene = defineSphere(fScene, 100, 100, 100, 10, 0.1f, 255, 255, 0);
    fScene = defineSphere(fScene, -100, -100, 0, 10, 0, 0, 0, 255);
    fScene = defineSphere(fScene, 100, -100, 100, 10, 0, 128, 128, 128);

    pixelPoint(WIDTH, pixelGridX, pixelGridY, pixelGridZ);
    

    // Create a blank image
    unsigned char* image = (unsigned char*)malloc(PWIDTH * PHEIGHT * CHANNELS);
    memset(image, 255, PWIDTH * PHEIGHT * CHANNELS);



    // Set every other pixel in the 50th row to red
    for (int y = 0; y < PHEIGHT; y++) {
        for (int x = 0; x < PWIDTH; x++) {
            struct ray current = {main.x, main.y, main.z, pixelGridX[x + y * PWIDTH], pixelGridY[x + y * PWIDTH], pixelGridZ[x + y * PWIDTH], 3};
            int color[3];
            castRay(color, current, scene, fScene, ground, main);


            image[3 * (x + y * PWIDTH) + 0] = color[0];  // Red channel
            image[3 * (x + y * PWIDTH) + 1] = color[1];  // Green channel
            image[3 * (x + y * PWIDTH) + 2] = color[2];  // Blue channel
        }
    }



    // Save the image as a PNG file
    stbi_write_png("outputC.png", PWIDTH, PHEIGHT, CHANNELS, image, PWIDTH * CHANNELS);

    // Clean up
    free(image);

    time = clock();

    printf("%f\n", (float)time / CLOCKS_PER_SEC);

    return 0;

}

struct sphere * defineSphere( struct sphere *fScene, float x, float y, float z, float R, float reflectivity, int r, int g, int b){
    fScene->x = x;
    fScene->y = y;
    fScene->z = z;
    fScene->R = R;
    fScene->reflectivity = reflectivity;
    fScene->color[0] = r;
    fScene->color[1] = g;
    fScene->color[2] = b;

    objects++;

    return ++fScene;
}


void pixelPoint(float width, float *pixelGridX, float *pixelGridY, float *pixelGridZ) {

    float height = PHEIGHT * width / PWIDTH;

    float deltaPixelX = width / PWIDTH;
    float deltaPixelY = height / PHEIGHT;

    for (int pX = 0; pX < PWIDTH; pX++) {
        float xCord = (deltaPixelX / 2) - (width / 2) + (pX * deltaPixelX);

        for (int pY = 0; pY < PHEIGHT; pY++) {

            pixelGridY[pX + pY * PWIDTH] = xCord;
            pixelGridX[pX + pY * PWIDTH] = 1;
            pixelGridZ[pX + pY * PWIDTH] = (height / 2) - (deltaPixelY / 2) - (pY * deltaPixelY);

        }
    }

}


void castRay(int *out, struct ray ray, struct sphere scene[], struct sphere *fScene, struct ground ground, struct camera camera) {

    float intersection[2];
    findIntersection(intersection, ray, scene, fScene, ground);

    //printf("intersection: %f : %f \n", intersection[0], intersection[1]);


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
        getGroundColor(direct, ground, x, y);
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

    for (int i = 0; i < 3; i++) {
        if (*(out + i) > 255){
            *(out+0) = 255;
        }
    }

}

void findIntersection(float* arr, struct ray ray, struct sphere scene[], struct sphere* fScene, struct ground ground) {
    float intersectionLength[objects];


    for (int i = 0; i < objects; i++) {
        struct sphere current = scene[i];
        float test = findSphereIntersection(current, ray);
        intersectionLength[i] = test;
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
    double t1;
    double t2;

    double A = ray.a, B = ray.b, C = ray.c;
    double X = ray.x, Y = ray.y, Z = ray.z;
    double h = sphere.x, k = sphere.y, l = sphere.z, r = sphere.R;

    double a = (A * A + B * B + C * C);
    double b = 2 * (A * (X - h) + B * (Y - k) + C * (Z - l));
    double c = -(r * r - (X * X + Y * Y + Z * Z - 2 * (X * h + Y * k + Z * l) + h * h + k * k + l * l));

    double center = -(b / (2 * a));
    double pm = sqrtf((b * b) - (4 * a * c)) / (2 * a);

    t1 = (center + pm);
    t2 = (center - pm);

    double t;

    if (isnan(t1) || isnan(t2))
        return -1;


    if (t1 < 0.01 || t2 < 0.01)
        return -1;


    t = t1;

    if (t2 < t1 && t2 > 0)
        t = t2;

    return (float)t;
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

void getGroundColor(int *out, struct ground ground, float x, float y){
    float size = 5;
    if (x < 0)
        x = x * -1 + size - 1;
    if (y < 0)
        y = y * -1 + size - 1;

    int X = (fmodf(x, size * 2) - size) < 0;
    int Y = (fmodf(y, size * 2) - size) < 0;

    if (X == Y){
        *(out+0) = ground.r1;
        *(out+1) = ground.g1;
        *(out+2) = ground.b1;
    }else{
        *(out+0) = ground.r2;
        *(out+1) = ground.g2;
        *(out+2) = ground.b2;
    }
}

