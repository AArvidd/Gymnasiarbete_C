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


struct ray{
    float x, y, z;
    float a, b, c;
    int depth;
};

struct sphere{
    float x, y, z;
    float R;
    float reflectivity;
    int r, g, b;
};

struct ground{
    float z;
    float reflectivity;
    float r1, g1, b1;
    float r2, g2, b2;
};





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

