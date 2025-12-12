#include "rtklib.h"
#include <stdio.h>
#include <stdint.h>   // for uint16_t
#include <stdlib.h>
#include <string.h>
#include <math.h>

/// <summary>
/// Opens a .bin DSM file, determines its size (in bytes), and computes how many uint16_t elevation samples it contains.
/// </summary>
/// <param name="file">FILE struct pointer</param>
/// <param name="fileName">Name of the .bin file</param>
/// <returns>It returns the number of elevation samples in the .bin file.</returns>
///
int read_BIN(file_BIN* file, const char* fileName)
{

    /* ----------------------------------------------------------
    // Copy the file name string into the file_BIN struct so the
    // struct stores which file it represents.
    // ----------------------------------------------------------- */
    strcpy(file->file_name, fileName);

    /* -----------------------------------------------------------
    // Open the .bin file in binary read mode ("rb").
    // "r" = read, "b" = binary (no newline translation).
    // ----------------------------------------------------------- */
    file->file_ptr = fopen(fileName, "rb");

    /* If file doesn't open, send an error message */
    if (!file->file_ptr) {
        fprintf(stdout, "Failed to open file\n");
        return 1;
    }

    /* Move the file pointer to the end so we can measure the total size of the file (in bytes) */
    fseek(file->file_ptr, 0, SEEK_END);

    /* ftell() tells us the current file position in bytes. Since the pointer is at the end it tells us the file size */
    file->file_size = ftell(file->file_ptr);

    /* Compute the number of elevation samples in the file. Each elevation sample is 2 bytes (which is sizeof(uint16_t) */
    const int n_data = (int)(file->file_size / sizeof(uint16_t));

    /* Now that we know the number of elevation samples, we can */
    rewind(file->file_ptr);

    return n_data;
}

/// <summary>
///
/// </summary>
/// <param name="num"></param>
/// <returns></returns>
double retrieve_first_digit_decimal(const double num)
{
    double num_subtract;
    if (num < 0) num_subtract = floor(-1 * num / 10) * -10;
    else num_subtract = floor(num / 10) * 10;

    double num_first_digit = num - num_subtract;
    if (num_first_digit < 0) num_first_digit *= -1;

    return num_first_digit;
}

/// <summary>
/// Test
/// </summary>
/// <param name="input"></param>
/// <param name="first_digit"></param>
/// <param name="step_size"></param>
/// <returns>test</returns>
double rounding_to_first_digit(const double input, const double first_digit, const int step_size)
{
    if (input < 0) return -1 * (step_size * round((-1 * input - first_digit) / step_size) + first_digit);
    return step_size * round((input - first_digit) / step_size) + first_digit;
}

/// <summary>
///
/// </summary>
/// <param name="EN"></param>
/// <param name="DSM"></param>
/// <returns></returns>
east_north get_closest_coordinate(const east_north* EN, const DSMData* DSM)
{
    east_north closest_EN;
    closest_EN.easting = rounding_to_first_digit(EN->easting, DSM->first_digit.easting, DSM->step_size);
    closest_EN.northing = rounding_to_first_digit(EN->northing, DSM->first_digit.northing, DSM->step_size);

    return closest_EN;
}

double calculate_true_height_meters(const DSMData* DSM, const int index)
{
    const double val = (double)(DSM->heights_array[index]);
    const double height = val / 100 + 1020;
    return height;
}

double calc_max_height(const DSMData* DSM) {
    int max_height_index = 0;
    int max_height_grid_code = 0;
    for (int i = 0; i < DSM->n_data_points; i++) {
        //printf("%d\n",i);
        if (DSM->heights_array[i] > max_height_grid_code) {
            max_height_grid_code = DSM->heights_array[i];
            max_height_index = i;
        }
    }
    double max_height_true = calculate_true_height_meters(DSM, max_height_index);
    printf("Max height: %f\n", max_height_true);
    return max_height_true;
}

/// <summary>
/// This function is used to populate the DSM struct and to load values into it from an opened .bin file.
/// </summary>
/// <param name="file">File information structure that contains an open FILE*</param>
/// <param name="file_name">Name of the DSM .bin file</param>
/// <param name="E_origin_DSM">Input easting origin of the DSM (ex. top left location easting of the DSM)</param>
/// <param name="N_origin_DSM">Input northing origin of the DSM (ex. top left location northing of the DSM)</param>
/// <param name="step_size">DSM raster spatial resolution (ex. 5m spatial resolution)</param>
/// <param name="n_columns">Number of columns in the DSM grid</param>
/// <param name="DSM">Output DSM struct to fill in from the raster data</param>
/// 
/*
file_BIN f;
DSM DSM;
create_WGS_84_ellipsoid();
initialize_dsm(&f, "data.bin", &DSM, -5243.600, 5657585.200, 5, 220);
*/
/// 
void initialize_dsm
(
    const char* file_name, /* Name of the DSM .bin file */
    DSMData* DSM, /* Output DSM struct to fill in from the raster data */
    double E_origin_DSM, /* Input easting origin of the DSM (ex. top left location easting of the DSM) */
    double N_origin_DSM, /* Input northing origin of the DSM (ex. top left location northing of the DSM)*/
    int step_size, /* DSM raster spatial resolution (ex. 5m spatial resolution) */
    int n_columns /* Number of columns in the DSM grid */
)
{
    file_BIN file;
    /* 1. Read how many elevation samples are in the DSM raster dataset.
     *    read_BIN() returns the number of 16-bit integer compressed height values */
    DSM->n_data_points = read_BIN(&file, file_name);

    /* 2. Copy spatial metadata from the raster into the DSM struct.
     *    These metadata values are determined beforehand. */
    DSM->origin_dsm.easting = E_origin_DSM;
    DSM->origin_dsm.northing = N_origin_DSM;
    DSM->step_size = step_size;
    DSM->n_columns = n_columns;

    /* 3. Count the number of rows in the square/rectangle DSM grid. */
    DSM->n_rows = DSM->n_data_points / DSM->n_columns;

    /* 4. Allocate memory into the DSM height array to fit the size of the .bin file.
     *    Each height value is stored as compressed uint16_t values (from 0-65535).
     *    And then read the file to store these values from the .bin file into the heights array. */
    DSM->heights_array = malloc(DSM->n_data_points * sizeof(uint16_t));
    fread(DSM->heights_array, sizeof(uint16_t), DSM->n_data_points, file.file_ptr);

    /* 5. Precompute the "first digit" values of the x/y origin.
     *    This is used to properly index through the compressed height values in heights_array */
    DSM->first_digit.easting = retrieve_first_digit_decimal(DSM->origin_dsm.easting);
    DSM->first_digit.northing = retrieve_first_digit_decimal(DSM->origin_dsm.northing);

    /* 6. Calculate the maximum height in the dataset. */
    DSM->max_dsm_height = calc_max_height(DSM);
}

steps_XY calculate_steps_from_origin(const east_north* point, const DSMData* DSM)
{
    const int steps_from_DSM_origin_E = (int)((point->easting - DSM->origin_dsm.easting) / DSM->step_size);
    const int steps_from_DSM_origin_N = (int)(-1 * (point->northing - DSM->origin_dsm.northing) / DSM->step_size);
    steps_XY steps;
    steps.steps_X = steps_from_DSM_origin_E;
    steps.steps_Y = steps_from_DSM_origin_N;
    return steps;
}

///
/// @param x_steps Number of X steps from the origin of the DSM
/// @param y_steps Number of Y steps from the origin of the DSM
/// @return True means that the point is outside the DSM bounds. False means that the point is within the bounds.
int out_of_bounds_check(int x_steps, int y_steps)
{
    if (x_steps < 0 || y_steps < 0) return 1;
    else return 0;
}

/* set_relative_origin ---------------------------------------------------
* set an origin to traverse along the DSM from
* args   : DSM       *DSM      I   DTM Object, see rtklib.h for definition and rnx2rtkp.c for initial setup
*          double    latitude  I   latitude in degrees
*          double    longitude I   longitude in degrees
* return : void (even if out of bounds, that will be handled in get relative height call, but we can change this if you need)*/
void set_relative_origin
(
    struct DSMData* DSM,
    const lat_long* relative_origin_degrees,
    const UTM_projection* proj,
    const ellipsoid* e,
    int* out_of_bounds
)
{
    project_latitude_longitude_to_UTM
    (
        &DSM->relative_origin_traverse,
        relative_origin_degrees,
        proj,
        e
    );
    const steps_XY steps = calculate_steps_from_origin(&DSM->relative_origin_traverse, DSM);

    *out_of_bounds = out_of_bounds_check(steps.steps_X, steps.steps_Y);
    if (*out_of_bounds == 0) {
        DSM->relative_origin_traverse = get_closest_coordinate(&DSM->relative_origin_traverse, DSM);
        printf("\nRelative origin E: %f\n", DSM->relative_origin_traverse.easting);
        printf("Relative origin N: %f\n", DSM->relative_origin_traverse.northing);
    }
    else printf("WARNING: RELATIVE ORIGIN IS OUTSIDE THE DSM BOUNDS");
}

/* get_relative_height ---------------------------------------------------
* Report height for a location with respect to the relative origin
* args   : DSM       *DSM           I   DSM Object, see rtklib.h for definition and rnx2rtkp.c for initial setup
*          int       *E             I   The number of steps East to take (distance = E * step_size) (negative indicates West)
*          int       *N             I   The number of steps North to take (distance = N * step_size) (negative indicates South)
*          double    *h             I   The height of the DSM at that point
*          int       *out_of_bounds I   1 if there is not DSM data at that point (h will not looked at if this is 1)
* return : void (data returned by setting h and out_of_bounds)*/
void get_relative_height
(
    const DSMData* DSM,
    const int* steps_E,
    const int* steps_N,
    double* h,
    int* out_of_bounds
)
{
    const double traverse_E = DSM->relative_origin_traverse.easting + (*steps_E) * (double)DSM->step_size;
    const double traverse_N = DSM->relative_origin_traverse.northing + (*steps_N) * DSM->step_size;

    east_north traverse = { traverse_E,traverse_N };

    const steps_XY steps = calculate_steps_from_origin(&traverse, DSM);

    *out_of_bounds = out_of_bounds_check(steps.steps_X, steps.steps_Y);

    if (*out_of_bounds == 0) {
        const int index = steps.steps_Y * DSM->n_columns + steps.steps_X;
        *h = calculate_true_height_meters(DSM, index);
        printf("Relative height: %f\n", *h);
    }
    else {
        *h = -1.0f;
        //printf("Traversing point is outside the DSM bounds. Cannot compute relative height");
    }

}

void deallocate_dsm(const DSMData* DSM) {
    free(DSM->heights_array);
    DSM = NULL;
}