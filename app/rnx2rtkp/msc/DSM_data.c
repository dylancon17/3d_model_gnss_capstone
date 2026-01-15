#include "rtklib.h"
#include <stdio.h>
#include <stdint.h>   // for uint16_t
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <conio.h>

//NOTE: PLAYGROUND IS IN LINE 195

void initialize_tiles_dataset
(
    TilesDataset* td,
    int num_tiles_x, 
    int num_tiles_y, 
    int tiles_dimension_x, 
    int tiles_dimension_y, 
    double top_left_tile_origin_x, 
    double top_left_tile_origin_y
)
{
    td->num_tiles_x = num_tiles_x;
    td->num_tiles_y = num_tiles_y;
    td->tiles_dimension_x = tiles_dimension_x;
    td->tiles_dimension_y = tiles_dimension_y;
    td->top_left_tile_origin.easting = top_left_tile_origin_x;
    td->top_left_tile_origin.northing = top_left_tile_origin_y;
    td->x_limit = td->top_left_tile_origin.easting + (td->tiles_dimension_x * td->num_tiles_x);
    td->y_limit = td->top_left_tile_origin.northing - (td->tiles_dimension_y * td->num_tiles_y);
}

int open_BIN(file_BIN* file, const char* fileName)
{

    //Set file to be empty
    memset(file, 0, size(*file));
    file->file_ptr = fopen(fileName, "rb");
    if (!file->file_ptr) {
        perror("\nOpening file unsuccessful\n");
        return -1;
    }

    strcpy(file->file_name, fileName);

    fseek(file->file_ptr, 0, SEEK_END);
    file->file_size = ftell(file->file_ptr);
    rewind(file->file_ptr);

    return (int)(file->file_size / sizeof(uint16_t)); // Return the number of data points in the .bin file
}

int read_BIN_data(file_BIN* file, uint16_t* buffer, int n) 
{
    if (!file || !file->file_ptr || !buffer) return -1;

    size_t read = fread(buffer, sizeof(uint16_t), n, file->file_ptr);
    return (read == (size_t)n) ? 0 : -1;
}

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

    printf("\nReading .bin file\n");

    if (file->file_ptr != NULL) { 
        printf("\nFile is not null\n");
        rewind(file->file_ptr); 
        printf("\nFile pointer has been reset.\n");
    }
    file->file_ptr = fopen(fileName, "rb");

    /* If file doesn't open, send an error message */
    if (!file->file_ptr) {
        fprintf(stdout, "Failed to open file\n");
        return 1;
    }
    else { printf("\nFile opened successfully\n"); }

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
double retrieve_anchor_decimal(const double num)
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
double round_to_anchor_step(const double input, const double anchor, const int step_size)
{
    return step_size * round((input - anchor) / step_size) + anchor;
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
    closest_EN.easting = round_to_anchor_step(EN->easting, DSM->first_digit.easting, DSM->step_size);
    closest_EN.northing = round_to_anchor_step(EN->northing, DSM->first_digit.northing, DSM->step_size);

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
    return max_height_true;
}

east_north round_to_tile_origin(const east_north* input, const TilesDataset* tiles_dataset) {
    double tile_origin_easting = round_to_anchor_step(input->easting, -2000, tiles_dataset->tiles_dimension_x);
    double tile_origin_northing = round_to_anchor_step(input->northing, 3000, tiles_dataset->tiles_dimension_y);

    if (tile_origin_easting > tiles_dataset->x_limit - tiles_dataset->tiles_dimension_x) { tile_origin_easting -= tiles_dataset->tiles_dimension_x; }
    if (tile_origin_northing < tiles_dataset->y_limit + tiles_dataset->tiles_dimension_y) { tile_origin_northing += tiles_dataset->tiles_dimension_y; }
    east_north tile_origin = { tile_origin_easting, tile_origin_northing };

    return tile_origin;
}

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
    printf("\nInitializing dsm\n");
    file_BIN file;
    /* Read how many elevation samples are in the DSM raster dataset. read_BIN() returns the number of 16-bit integer compressed height values */
    int n_samples = open_BIN(&file,file_name);
    if (n_samples < 0) {
        fprintf(stderr, "Failed to open DSM file\n");
    }

    DSM->n_data_points = n_samples;

    /* Copy spatial metadata from the raster into the DSM struct. These metadata values are determined beforehand. */
    DSM->origin_dsm.easting = E_origin_DSM;
    DSM->origin_dsm.northing = N_origin_DSM;
    DSM->step_size = step_size;
    DSM->n_columns = n_columns;
    DSM->n_rows = DSM->n_data_points / DSM->n_columns;
    DSM->tile_size_x = DSM->step_size * DSM->n_columns;
    DSM->tile_size_y = DSM->step_size * DSM->n_rows;

    /* Allocate memory into the DSM height array to fit the size of the .bin file. */
    DSM->heights_array = malloc(DSM->n_data_points * sizeof(uint16_t));
    if (!DSM->heights_array) {
        perror("\nmalloc failed\n");
        fclose(file.file_ptr);
        return;
    }

    if (read_BIN_data(&file, DSM->heights_array, (size_t)DSM->n_data_points) != 0) {
        fprintf(stderr, "Failed to read DSM raster data\n");
        free(DSM->heights_array);
        fclose(file.file_ptr);
        return;
    }

    /* Precompute the anchor decimal values of the x/y origin. */
    DSM->first_digit.easting = retrieve_anchor_decimal(DSM->origin_dsm.easting);
    DSM->first_digit.northing = retrieve_anchor_decimal(DSM->origin_dsm.northing);

    /* Calculate the maximum height in the dataset. */
    DSM->max_dsm_height = calc_max_height(DSM);

    fclose(file.file_ptr);
    file.file_ptr = NULL;

    // PLAYGROUND =========================
    printf("\nStarted playground\n");
    double num = -123.856789;
    char buffer[100]; // Enough space for the converted string

    // Convert double to string with 6 decimal places
    if (snprintf(buffer, sizeof(buffer), "%d", (int)round(num)) < 0) {
        fprintf(stderr, "Error converting double to string.\n");
        return 1;
    }

    printf("\nDouble as string: %s\n", buffer);

    char result[64];

    const char* a = "abc";
    const char* b = "123";

    snprintf(result, sizeof(result), "%s_%s", a, b);
    printf("\nAdding strings together: %s\n", result);

    printf("\nPress 'q' to continue...\n");

    int ch;
    do {
        ch = _getch();   // no Enter required, no echo
    } while (ch != 'q');

    printf("\nContinuing...\n");

    printf("\nEnded playground\n");

    /*
    printf("\nTraverse easting: %f\n", traverse.easting);
    printf("\nTraverse northing: %f\n", traverse.northing);
    printf("\nRounded easting: %f\n", rounded.easting);
    printf("\nRounded northing: %f\n", rounded.northing);
    */

    printf("\nInitialization complete\n");
}

steps_XY calculate_steps_from_tile_corner(const east_north* point, const DSMData* DSM)
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
int out_of_bounds_check_tiles_dataset(east_north* traverse, TilesDataset* tiles_dataset)
{
    if (traverse->easting < tiles_dataset->top_left_tile_origin.easting) return 1;
    else if (traverse->easting > tiles_dataset->x_limit) return 1;

    else if (traverse->northing > tiles_dataset->top_left_tile_origin.northing) return 1;
    else if (traverse->northing < tiles_dataset->y_limit) return 1;

    else return 0; 
}

int out_of_bounds_check(int x_steps, int y_steps, DSMData* DSM)
{
    int max_steps_x = DSM->n_columns;
    int max_steps_y = DSM->n_rows;
    if (x_steps < 0 || y_steps < 0 || x_steps > max_steps_x || y_steps > max_steps_y) return 1;
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
    DSMData* DSM,
    const TilesDataset* tiles_dataset,
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
    const steps_XY steps = calculate_steps_from_tile_corner(&DSM->relative_origin_traverse, DSM);

    *out_of_bounds = out_of_bounds_check(steps.steps_X, steps.steps_Y, DSM);
    int out_of_bounds_tiles_dataset = out_of_bounds_check_tiles_dataset(&DSM->relative_origin_traverse, tiles_dataset);

    if (*out_of_bounds == 0 && out_of_bounds_tiles_dataset == 0) {
        DSM->relative_origin_traverse = get_closest_coordinate(&DSM->relative_origin_traverse, DSM);
    }
    else if (*out_of_bounds == 1 && out_of_bounds_tiles_dataset == 0) { 
        printf("\nTODO: OPEN A NEW TILE FOR DSM COMPUTATIONS\n"); 
        double traverse_E = get_closest_coordinate(&DSM->relative_origin_traverse, DSM).easting;
        double traverse_N = get_closest_coordinate(&DSM->relative_origin_traverse, DSM).northing;

        east_north traverse_coords = { traverse_E,traverse_N };
        east_north traverse_to_tile_origin = round_to_tile_origin(&traverse_coords, tiles_dataset);

        char traverse_E_char[100];
        char traverse_N_char[100];
        if (snprintf(traverse_E_char, sizeof(traverse_E_char), "%d", (int)round(traverse_to_tile_origin.easting)) < 0) {
            fprintf(stderr, "Error converting double to string.\n");
            return 1; 
        }
        if (snprintf(traverse_N_char, sizeof(traverse_N_char), "%d", (int)round(traverse_to_tile_origin.northing)) < 0) {
            fprintf(stderr, "Error converting double to string.\n");
            return 1;
        }

        char file_path[100] = "C:\\capstone\\dsm_tiles\\DSM_CGY_5x5km_res1m\\";
        char file_prefix[100] = "DSM_CGY_5x5km_res1m";
        char new_file_name[100];
        char file_extension[100] = ".bin";

        if (snprintf(new_file_name, sizeof(new_file_name), "%s%s_%s_%s%s", file_path, file_prefix, traverse_E_char, traverse_N_char, file_extension) < 0) {
            fprintf(stderr, "Error converting double to string.\n");
            return 1;
        }
        printf("new_file_name: %s",new_file_name);
        free(DSM->heights_array);
        initialize_dsm(new_file_name, &(DSM), traverse_to_tile_origin.easting, traverse_to_tile_origin.northing, 1, 5000);

        printf("\nPress 'c' to continue...\n");

        int ch;
        do {
            ch = _getch();   // no Enter required, no echo
        } while (ch != 'c');

        printf("\nContinuing...\n");
    }
    else if (out_of_bounds_tiles_dataset == 1) { printf("\WARNING: RELATIVE ORIGIN IS OUTSIDE THE DSM TILES DATASET\n"); }

    /*
    double test_input_x = 12000.0;
    double test_input_y = 5635000.0;
    east_north test_input = { test_input_x, test_input_y };
    east_north test_output = round_to_tile_origin(&test_input, tiles_dataset);

    printf("\ntest output x: %f\n", test_output.easting);
    printf("\ntest output y: %f\n", test_output.northing);
    */
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
    const TilesDataset* tiles_dataset,
    const int* steps_E,
    const int* steps_N,
    double* h,
    int* out_of_bounds
)
{
    const double traverse_E = DSM->relative_origin_traverse.easting + (*steps_E) * (double)DSM->step_size;
    const double traverse_N = DSM->relative_origin_traverse.northing + (*steps_N) * (double)DSM->step_size;

    east_north traverse = { traverse_E,traverse_N };

    east_north rounded = round_to_tile_origin(&traverse, tiles_dataset);

    /*
    printf("\nTraverse easting: %f\n", traverse.easting);
    printf("\nTraverse northing: %f\n", traverse.northing);
    printf("\nRounded easting: %f\n", rounded.easting);
    printf("\nRounded northing: %f\n", rounded.northing);
    */

    const steps_XY steps = calculate_steps_from_tile_corner(&traverse, DSM);

    *out_of_bounds = out_of_bounds_check_tiles_dataset(&traverse, tiles_dataset);
    int out_of_bounds_tiles_dataset = out_of_bounds_check_tiles_dataset(&DSM->relative_origin_traverse, tiles_dataset);

    if (*out_of_bounds == 0 && out_of_bounds_tiles_dataset == 0) {
        const int index = steps.steps_Y * DSM->n_columns + steps.steps_X;
        *h = calculate_true_height_meters(DSM, index);
        //printf("\nRelative height calculated: %f\n", *h);
    }
    else if (*out_of_bounds == 1 && out_of_bounds_tiles_dataset == 0) { 

        /*
        char result[64];

        const char* a = "abc";
        const char* b = "123";

        snprintf(result, sizeof(result), "%s_%s", a, b);
        printf("Adding strings together: %s", result);
        */

        //printf("\nTODO: OPEN A NEW TILE FOR DSM COMPUTATIONS\n"); 
        char traverse_E_char[100];
        char traverse_N_char[100];
        if (snprintf(traverse_E_char, sizeof(traverse_E_char), "%.6f", traverse_E) < 0) {
            fprintf(stderr, "Error converting double to string.\n");
            return 1;
        }

        char file_prefix[100] = "DSM_CGY_5x5km_res1m";
        char new_file[100];
        if (snprintf(new_file, sizeof(new_file), "%s_%s_%s", file_prefix, traverse_E_char, traverse_N_char) < 0) {
            fprintf(stderr, "Error converting double to string.\n");
            return 1;
        }


    }
    else if (out_of_bounds_tiles_dataset == 1) { 
        printf("Traversing point is outside the DSM bounds. Cannot compute relative height");
        *h = -1.0f;
    }
}

void deallocate_dsm(const DSMData* DSM) {
    free(DSM->heights_array);
    DSM = NULL;
}