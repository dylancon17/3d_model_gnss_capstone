//
// Created by titus on 2025-10-29.
//

#include <stdio.h>
#include <stdint.h>   // for uint16_t
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct fileBIN {
    char fileName[256]; // char to hold the file name
    long fileSize; // file size
    FILE *filePtr; // file pointer data member
};

struct DTMData {
    // X and Y origin
    double xorig;
    double yorig;

    double x_relativeorig;
    double y_relativeorig;

    int step_size; // Step size (i.e. spatial resolution of the DTM)
    int nrows; // number of rows in the elevation (square/rectangular) dataset
    int ncols; // number of columns in the (square/rectangular) dataset
    int ndata; // number of data points in the dataset (should equal nrows * ncols)
    uint16_t *coords; // data member array to hold the parsed heights

    // "First digit" value, where for 5243.6, the "first digit" would be 3.6
    double x_first_digit;
    double y_first_digit;
    double maxHeight;
};

struct XY {
    double x;
    double y;
};

int readBIN(struct fileBIN *file, const char *fileName) {
    strcpy(file->fileName, fileName);
    file->filePtr = fopen(fileName, "rb");
    if (!file->filePtr) {
        perror("Failed to open file");
        return 1;
    }
    printf("File \"%s\" opened successfully.\n", fileName);
    fseek(file->filePtr, 0, SEEK_END);
    const int file_size = ftell(file->filePtr);
    rewind(file->filePtr);
    const int ndata = file_size / sizeof(uint16_t);

    printf("File size: %u bytes\n", file_size);
    printf("Number of points: %u individual coordinates \n", ndata);

    return ndata;
}

double retrieveFirstDigitDecimal(const double num) {
    double num_subtract;
    if (num<0) num_subtract = floor(-1*num/10) * -10;
    else num_subtract = floor(num/10) * 10;

    double num_first_digit = num - num_subtract;
    if (num_first_digit <0) num_first_digit *= -1;

    printf("num_first_digit: %f\n", num_first_digit);
    return num_first_digit;
}

double roundingAlgorithm(const double input, const double first_digit, const double step_size) {
    if (input < 0) return -1*(step_size * round((-1*input - first_digit)/step_size) + first_digit);
    return step_size * round((input - first_digit)/step_size) + first_digit;
}

struct XY getClosestDTMCoord(const struct XY gps, const struct DTMData dtm) {
    struct XY dtm_coord;
    dtm_coord.x = roundingAlgorithm(gps.x, dtm.x_first_digit, dtm.step_size);
    dtm_coord.y = roundingAlgorithm(gps.y, dtm.y_first_digit, dtm.step_size);
    return dtm_coord;
}

void set_relative_origin(struct DTMData *dtm, double xorig_new, double yorig_new) {
    dtm->x_relativeorig = xorig_new;
    dtm->y_relativeorig = yorig_new;
}

double calcMaxHeight(struct DTMData *dtm) {
    int maxVal = dtm->coords[0];
    for (int i = 0; i < dtm->ndata; i++) {
        if (dtm->coords[i] > maxVal) maxVal = dtm->coords[i];
    }
    double maxHeight = (double)maxVal / 100 + 1020;
    return maxHeight;
}

void initializeDTM(
    struct fileBIN *file,
    const char *fileName,
    struct DTMData *dtm,
    double xorig,
    double yorig,
    int stepsize,
    int ncols) {

    dtm->ndata = readBIN(file, fileName);
    dtm->xorig = xorig;
    dtm->yorig = yorig;
    dtm->step_size = stepsize;
    dtm->ncols = ncols;
    dtm->nrows = dtm->ndata/dtm->ncols;
    dtm->coords = malloc(dtm->ndata * sizeof(uint16_t));
    dtm->x_first_digit = retrieveFirstDigitDecimal(dtm->xorig);
    dtm->y_first_digit = retrieveFirstDigitDecimal(dtm->yorig);
    fread(dtm->coords, sizeof(uint16_t), dtm->ndata, file->filePtr);
    dtm->maxHeight = calcMaxHeight(&dtm);
    fclose(file->filePtr);
}

double retrieveHeightFromIndex(const struct DTMData* dtm, const int  i) {
    double val = dtm->coords[i];
    double height = val / 100 + 1020;
    return height;
}

void outOfBoundsCheck(double x_steps, double y_steps) {
    if (x_steps < 0 || y_steps < 0) {
        fprintf(stderr, "Error: Input coordinate exceeds the DTM bounds.\n");
        exit(EXIT_FAILURE);   // portable error code
    }
}

double retrieveHeight(const struct XY xy, const struct DTMData *dtm) {
    double xrounded = roundingAlgorithm(xy.x, dtm->x_first_digit, dtm->step_size);
    double yrounded = roundingAlgorithm(xy.y, dtm->y_first_digit, dtm->step_size);

    int num_steps_x = (int)((xrounded - dtm->xorig)/ dtm->step_size);
    int num_steps_y = (int)(-1*(yrounded - dtm->yorig)/ dtm->step_size);

    outOfBoundsCheck(num_steps_x, num_steps_y);

    int index = num_steps_y * dtm->ncols + num_steps_x;
    return retrieveHeightFromIndex(dtm, index);
}

// Getters
double getMaxHeight(struct DTMData dtm) {
    return dtm.maxHeight;
}

int main() {
    struct fileBIN f;
    struct DTMData dtm;

    initializeDTM(&f, "data.bin", &dtm, -5243.600, 5657585.200, 5, 220);

    struct XY gps;
    gps.x = -5220.0; // 0 - 219
    gps.y = 5657560; // 0 - 269

    double x = gps.x; //printf("GPS x: %f\n", x);

    double y = gps.y; //printf("GPS y: %f\n", y);

    double test_x = roundingAlgorithm(gps.x, dtm.x_first_digit, dtm.step_size);
    printf("x orig: %.3f, ", dtm.xorig);
    printf("Rounding algorithm for x: %.3f, ", test_x);
    printf("dx: %d\n", (int)(test_x - dtm.xorig));

    int num_steps_x = (int)((test_x - dtm.xorig)/dtm.step_size);

    double test_y = roundingAlgorithm(gps.y, dtm.y_first_digit, dtm.step_size);
    printf("y orig: %.3f, ", dtm.yorig);
    printf("Rounding algorithm for y: %.3f, ", test_y);
    printf("dy: %d\n", (int)(test_y - dtm.yorig));

    int num_steps_y = (int)(-1*(test_y - dtm.yorig)/dtm.step_size);

    printf("Number of steps for x: %d\n", num_steps_x);
    printf("Number of steps for y: %d\n", num_steps_y);

    int indexing = num_steps_y*dtm.ncols + num_steps_x;
    printf("Point number: %d\n", indexing + 1);
    printf("dtm value: %hu\n", dtm.coords[indexing]);

    double testfun = retrieveHeight(gps,&dtm);
    printf("Test retrieveHeight: %.2f meters\n",testfun);

    return 0;
}