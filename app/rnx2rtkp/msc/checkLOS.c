#include "rtklib.h"

#include <math.h>


struct DTMData;
struct LineState;

typedef struct DTMData{
    int max_distance_check;

    // Function pointers act like methods
    void (*set_relative_origin)(struct DTMData* self, float latitude, float longitude);
    void (*get_relative_height)(struct DTMData* self, int* E, int* N, float* h, boolean* out_of_bounds);
} DTMData;

typedef struct LineState {
    int E, N; //Current point on line
    float d; // Distance travelled along line (not Euclidean)
    int dE, dN; //Absolute step size to each point
    int sE, sN; // Unit direction to step along line 
    int err, e2; // Error tracking
} LineState;

boolean check_los(float sat_az, float sat_elev, float origin_lat, float origin_long, float origin_height, DTMData* DTM);
void step_along_line(LineState* l);


//void los_update() #TODO set imports
//
//for sat in RTKLibSatObject{// TODO properly parse Sat
//    obstructed = check_los(sat_az, sat_elev, origin_lat, origin_long, origin_height, threshold, DTMClass);
//    if not check_los(sat_az, sat_elev, origin_lat, origin_long, origin_height, threshold, DTMClass) {
//        RTKLibObject.removeobservation(sat);
//    }
//}

boolean check_los(float sat_az, float sat_elev, float origin_lat, float origin_long, float origin_height, DTMData* DTM) {
    
    //Set up DTM call
    DTM->set_relative_origin(DTM, origin_lat, origin_long);
    boolean out_of_bounds = 0;
    
    // Set up height checks
    float max_checked_DTM_height = origin_height;
    float current_DTM_height = 0, sat_height = 0;
    float sat_vertical_slope = tanf(sat_elev);

    // Draws a line from (0,0) to (E1, N1) that's the maximum distance required to check
    int E1 = (int)roundf(DTM->max_distance_check * sinf(sat_az));
    int N1 = (int)roundf(DTM->max_distance_check * cosf(sat_az));

    // Set up Bresenham Algorithm
    int dE = abs(E1), sE = 0 < E1 ? 1 : -1;
    int dN = -abs(N1), sN = 0 < N1 ? 1 : -1;; // Defines directions of steps
    LineState line = {0, 0, 0, dE, dN, sE, sN, dE + dN, 0};


    while (1) {
        // Traverse the DTM
        step_along_line(&line);
        DTM->get_relative_height(DTM, &(line.E), &(line.N), &current_DTM_height, &out_of_bounds);

        // If fully traversed DTM, LOS is clear
        if (out_of_bounds) {return 1;}

        // If a higher height has already been checked (or lower than starting height), LOS is clear at that point. Continue traversing
        if (current_DTM_height < max_checked_DTM_height) {continue;}
        
        //Else, onto the new max height checked
        max_checked_DTM_height = current_DTM_height;

        // Check the height of the satellite based on the distance travelled along the satellite line
        sat_height = sat_vertical_slope * line.d + origin_height;

        // Satellite is lower than DTM, LOS is obstructed
        if (sat_height < current_DTM_height) { return 0; }

        // Have traversed the max distance required to trraverse, LOS is clear
        if (line.E == E1 && line.N == N1) { return 1; }

        // Otherwise continue traversing
    }
}

void step_along_line(LineState* l)
{
    l->e2 = 2 * l->err;

    if (l->e2 >= l->dN) {
        l->err += l->dN;
        l->E += l->sE;
    }
    
    if (l->e2 <= l->dE) {
        l->err += l->dE;
        l->N += l->sN;
    }

    if (l->e2 <= l->dE && l->e2 >= l->dN) {
        l->d += 1.41421356237f; //Root 2
    }
    else {
        l->d += 1;
    }
}