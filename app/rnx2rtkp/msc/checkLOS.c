#define _USE_MATH_DEFINES

#include "rtklib.h"

#include <stdlib.h>

#include <math.h>




typedef struct LineState {
    int E, N; //Current point on line
    double d; // Distance travelled along line (not Euclidean)
    int dE, dN; //Absolute step size to each point
    int sE, sN; // Unit direction to step along line 
    int err, e2; // Error tracking
} LineState;

void step_along_line(LineState* l);
void apply_probability(double* y, double* probability);
void determine_DTM_height_var(double* var, struct DSMData* DSM);
void determine_sat_height_var(double* var, double origin_horizontal_variance, double origin_vertical_variance, double sat_vertical_slope, struct DSMData* DTM);
void soltocov_rtk(sol_t* sol, double* P);

int calc_expected_los(rtk_t* rtk, const nav_t* nav, gtime_t tor, double* rr, double* pos, double* Q, double* probability_sum) {
    double rs[6], dts[2], var;
    int sat;
    int svh[2];

    double r;

    double elev[3], azel[2];

    //TODO this whole section could be done better as sat positions are being calculated twice!
    // WARNING - this uses the time of reception, not time of transmission. This is incorrect! Expected error is 1-4", so likely insignificant?
    int num_possible = 0;

    /* ---- GPS/GAL/BDS/QZSS/SBAS broadcast ephemerides ---- */
    for (int i = 0; i < nav->n; i++) {
        eph_t* e = &nav->eph[i];
        sat = e->sat;


        // Get sat position
        if (!satpos(tor, tor, sat, EPHOPT_BRDC, nav, rs, dts, &var, svh)) continue;

        r = geodist(rs + i * 6, rr, elev); //TODO how is rs indexed
        /* geodist failure check */
        if (r <= 0) {
            continue;
        }

        satazel(pos, elev, azel);

        rtk->ssat[i - 1].obstruction_probability = check_los(azel[0], azel[1], pos[0], pos[1], pos[2], Q[0] + Q[4], Q[8], &(rtk->opt.dtm));
    
        num_possible++;
        *probability_sum += rtk->ssat[i - 1].obstruction_probability;

//        if (rtk->ssat[i - 1].obstruction_probability < rtk->opt.dtm.rejection_threshold) {
//            num_expected++;
//        }
    }

    /* ---- GLONASS broadcast ephemerides ---- */
    for (int i = 0; i < nav->ng; i++) {
        geph_t* g = &nav->geph[i];
        sat = g->sat;

        if (!satpos(tor, tor, sat, EPHOPT_BRDC, nav, rs, dts, &var, svh)) continue;

        r = geodist(rs + i * 6, rr, elev); //TODO how is rs indexed
        /* geodist failure check */
        if (r <= 0) {
            continue;
        }

        satazel(pos, elev, azel);

        rtk->ssat[i - 1].obstruction_probability = check_los(azel[0], azel[1], pos[0], pos[1], pos[2], Q[0] + Q[4], Q[8], &(rtk->opt.dtm));

        num_possible++;
        *probability_sum += rtk->ssat[i - 1].obstruction_probability;


//        if (rtk->ssat[i - 1].obstruction_probability < rtk->opt.dtm.rejection_threshold) {
//            num_expected++;
//        }
    }

    return num_possible;
}

/* los update ---------------------------------------------------
* check and update observations based off of line of sight
* args   : rtl_t     *rtk   I   rtk object. Assumed that position is propagated and initialized. Includes link to DEM
*          obsd_t    *level I   rover and base observations. Rover is all first then base. Can also check receiver flag.
*          int       *sat   I   PRN of common satellite[k]
*          int       *iu    I   index of rover sat[k] in obsd_t
*          int       *ir    I   index of base sat[k] in obsd_t
*          int       *ns    I   number of common satellites
*          double    *rs    I   satellite positions in ECEF
* return : int (number of rejected sats) */
extern int los_update(rtk_t* rtk, const obsd_t* obs, const nav_t* nav, gtime_t tor, int* sat, int* iu, int* ir, int* ns, const double* rs) {
    //fprintf(stderr, "%s\n", "LOS Update\n");

    int i, nrej=0;
    int rej_idx[MAXOBS];
    const double *rr = rtk->sol.rr; // TODO add filtering for if the estimated position is poor

    // Converts ECEF to LLH
    double pos[3];
    ecef2pos(rr, pos);
    
    // Converts ECEF covars to LLH
    double P[9]; // 3x3 ENU
    double Q[9];
    
    soltocov_rtk(&(rtk->sol), P);
    covenu(pos, P, Q);

    lat_long ll = { pos[0] * 180 / M_PI, pos[1] * 180 / M_PI};
    printf("\nRelative origin latitude: %f\n", ll.latitude);
    printf("\nRelative origin longitude: %f\n", ll.longitude);

    int out_of_bounds = 0;
    // Set relative origin here as it is constant for the rest of the update
    // TODO, out of bounds filtering may be required...can be used as an optimization. Must reset the scaling if done so
    set_relative_origin(&(rtk->opt.DSM),&(rtk->opt.tiles_dataset), &ll, &(rtk->opt.UTM), &(rtk->opt.ellip), &out_of_bounds);

    double probability_total = 0;

    int num_possible = calc_expected_los(rtk, nav, tor, rr, pos, Q, &probability_total);

    if (num_possible < 1) {
        return 0;
    }

    int num_not_observed = num_possible;

    double observed_probability_sum = 0;

    double r;
    double probability_of_obstruction = 0;

    double current_DTM_height;

    if (out_of_bounds == 1) { //If origin is out of bounds, don't search farther than that
        return 0;
    }

    gtime_t obs_time, point_time;
    int debug = 0;
    if (*ns > 0) {
        obs_time = obs[0].time;
        /* get current truth time */
        point_time = gpst2time(2258, 160665);
        if (abs(timediff(obs_time, point_time)) < 0.1) {
            debug = 1;
        }
    }
    else {
        return 0;
    }

    for (i = 0;i < *ns && i < MAXOBS;i++) {
        
        probability_of_obstruction = rtk->ssat[sat[i]-1].obstruction_probability;

        // Reject the signal if it's likely that it's obstructed (works for DEM processing options 1 and 2)
        if (probability_of_obstruction > rtk->opt.DSM.rejection_threshold) {

            rej_idx[nrej++] = i;
        }
        if (debug) {
            fprintf(stderr, "Setting ssat at index %d with probability: %lf\n", sat[i] - 1, probability_of_obstruction);
        }

        rtk->ssat[sat[i] - 1].obstruction_probability = probability_of_obstruction;

        if (probability_of_obstruction == 1.0) { // Avoid a divide by 0 error
            probability_of_obstruction = 0.999;
        }

        if (probability_of_obstruction < 0.0) { // Handles uninitilized case by setting no scaling
            probability_of_obstruction == 0.0;
        }

        double scaling = 1.0 / (1.0 - probability_of_obstruction); // TODO divide by 0 risk here
        if (scaling > rtk->opt.DSM.max_noise_scaling)
        {
            scaling = rtk->opt.DSM.max_noise_scaling;
        }

        // Save data. Warning! This will last across epochs unless overwritten (actually ssat is reset every epoch). Shouldn't be an issue unless implementation changed
        //fprintf(stderr, "%lf", scaling);
        if (debug) {
            fprintf(stderr, "Setting ssat at index %d with scaling: %lf\n", sat[i] - 1, scaling);
        }
        rtk->ssat[sat[i] - 1].obstruction_scaling = scaling;

        observed_probability_sum += probability_of_obstruction;
        probability_total -= probability_of_obstruction;
        num_not_observed--;

    }

    double average_probability_error = (observed_probability_sum + (num_not_observed - probability_total)) / num_possible;

    // If our prediction was more than a threshold incorrect, our probability calc is likely incorrect due to an incorrect starting position
    // Therefore, don't reject deweight sats to prevent us from continuing along the wrong path
    // TODO - if this happens consider undoing the previous state update as it's guranteed wrong (or was weighted not enough to fix it)
    if (average_probability_error > rtk->opt.dtm.average_prob_error_max && rtk->opt.dtm.processing_type > 5) {
        // Undo any scaling
        for (i = 0;i < *ns && i < MAXOBS;i++) {
            rtk->ssat[sat[i] - 1].obstruction_scaling = 1.0;
        }
        return 0;
    }


    // If deterministic or probabilistic rejection
    if (rtk->opt.DSM.processing_type > 0 && rtk->opt.DSM.processing_type != 4) {
        // Remove the rejected indices
        for (i = nrej - 1; i >= 0; i--) {
            int idx = rej_idx[i];
            memmove(sat + idx, sat + idx + 1, (*ns - idx - 1) * sizeof(int));
            memmove(iu + idx, iu + idx + 1, (*ns - idx - 1) * sizeof(int));
            memmove(ir + idx, ir + idx + 1, (*ns - idx - 1) * sizeof(int));

            (*ns)--;
        }
    }

    //fprintf(stdout, "%d sats rejected\n", nrej);


    return  nrej;
    
}

//Assumes relative origin already set
extern double check_los(double sat_az, double sat_elev, double origin_lat, double origin_long, double origin_height, double origin_horizontal_variance, double origin_vertical_variance, struct DSMData* DSM, TilesDataset* tiles_dataset, int debug) {
    if (debug) {
        fprintf(stderr, "Checking line of sight for: az: %lf elev: %lf at lat: %lf long: %lf height: %lf with hor var: %lf, vert var: %lf\n",
            sat_az * 180 / M_PI,
            sat_elev * 180 / M_PI,
            origin_lat * 180 / M_PI,
            origin_long * 180 / M_PI,
            origin_height,
            origin_horizontal_variance,
            origin_vertical_variance);
    }

    int out_of_bounds = 0;
    double current_DTM_height = 0, sat_height = 0;
    double sat_vertical_slope = tanf(sat_elev);
    int origin_x = 0, origin_y = 0;
    double distance = 0;
    //If the starting height is below the DEM, use the DEM height, or if the option to use_dem_height_only is set
    get_relative_height(DSM, tiles_dataset, &origin_x, &origin_y, &distance, &current_DTM_height, &out_of_bounds);

    double current_building_height = current_DTM_height;
    int current_plane_checked = 1; //Don't check starting plane
    double height_to_check = current_DTM_height;
    double change_in_height_var = pow(DSM->step_size,2);

    double probability_of_obstruction = -1; // -1 = uninitialized - otherwise a range of 0-1

    if (out_of_bounds) {
        return probability_of_obstruction; //-1
    }

    if (debug) {
        fprintf(stderr, "Origin Height: %lf, DEM Height: %lf, Out of Bounds %d",
            origin_height,
            current_DTM_height,
            out_of_bounds);
    }


    if (out_of_bounds != 1 && (current_DTM_height > origin_height || DSM->use_dem_height_only == 1) && DSM->use_dem_height_only != 2) {
        origin_height = current_DTM_height + DSM->antenna_dem_offset;
        origin_vertical_variance = DSM->vertical_point_variance + DSM->antenna_dem_offset_var;
        if (debug) {
            fprintf(stderr, ", Using DTM Height\n");
        }
    }
    else {
        if (debug) {
            fprintf(stderr, ", Using GPS Height\n");
        }
    }
    // Set up height checks
    double max_checked_DTM_height = origin_height;
    // Only search until a distance where you're guranteed to hit a building, or you've searched a ways
    double max_distance_steps = (DSM->max_dsm_height - origin_height) / sat_vertical_slope;
    if (max_distance_steps > DSM->max_distance) {
        max_distance_steps = DSM->max_distance;
    }
    // Scale to units of steps instead of meters
    max_distance_steps = max_distance_steps / DSM->step_size;

    // Draws a line from (0,0) to (E1, N1) that's the maximum distance required to check
    int E1 = (int)roundf(max_distance_steps * sinf(sat_az)); // Find distance to end point and then reduce by step size to reduce the number of required steps
    int N1 = (int)roundf(max_distance_steps * cosf(sat_az));
    // fprintf(stderr, "Determining End Point max distance: %i sinf: %lf cosf: %lf float solution: %lf rounded float solution: %lf",
    //      DTM->max_distance_check, sinf(sat_az), cosf(sat_az), DTM->max_distance_check* sinf(sat_az), roundf(DTM->max_distance_check * sinf(sat_az)));
    // Set up Bresenham Algorithm
    int dE = abs(E1), sE = 0 < E1 ? 1 : -1; // Each step represents the step size
    int dN = -abs(N1), sN = 0 < N1 ? 1 : -1; // Defines directions of steps
    LineState line = {0, 0, 0, dE, dN, sE, sN, dE + dN, 0};

    double DTM_height_var = 0;
    double sat_height_var = 0;
    double y = 0;

    while (1) {
        // Traverse the DTM
        step_along_line(&line);
        
        get_relative_height(DSM, tiles_dataset, &(line.E), &(line.N), &(line.d), &current_DTM_height, &out_of_bounds);

        //fprintf(stderr, "---DTM Height: %lf, dE: %d: dN %d, out_of_bounds: %d---", current_DTM_height, line.E, line.N, out_of_bounds);
        if (debug) {
            fprintf(stderr, "Line State: E: %d N: %d d: %lf dE: %d dN: %d sE: %d sN: %d err: %d e2: %d Height Comparison: %lf Boundary: %d\n",
                line.E, line.N, line.d, line.dE, line.dN, line.sE, line.sN, line.err, line.e2, current_DTM_height, out_of_bounds);
        }
        // If fully traversed DTM, LOS is clear after that point, return current probability
        if (out_of_bounds == 1) {
            if (debug) {
                fprintf(stderr, "out of bounds %lf\n", probability_of_obstruction);
            }
            return probability_of_obstruction;
        }

        // Have traversed the max distance required to trraverse, LOS is clear after that point. Depends on the number of steps travelled
        if (max_distance_steps < line.d) {
            if (probability_of_obstruction < 0.0) { // If uninitialized, initialize
                probability_of_obstruction = 0.0;
            }

            if (debug) {
                fprintf(stderr, "covered max distance %lf\n", probability_of_obstruction);
            }
            return probability_of_obstruction;
        }
        
        
        if (DSM->processing_type >= -3 && DSM->processing_type <= 6) {
            // If a higher height has already been checked, it's not needed to check it again, it's assumed to be low probability of obstruction. NOTE - this is a simplification for efficiency purposes. It also means only one probability per building is estiimated (assuming the building height is ~ constant)
            if (current_DTM_height <= max_checked_DTM_height + 2) {
                if (probability_of_obstruction < 0.0) { // If uninitialized, initilize
                    probability_of_obstruction = 0.0;
                }
                if (debug) {
                    fprintf(stderr, "Skipping Check as Max Checked is Higher %lf\n", max_checked_DTM_height);
                }
                continue;
            }
            max_checked_DTM_height = current_DTM_height;
            height_to_check = current_DTM_height;
        }

        //Height not changing
        if (DSM->processing_type > 6 || DSM->processing_type < -3) {
            if (current_DTM_height > current_building_height - DSM->building_height_margin && current_DTM_height < current_building_height + DSM->building_height_margin) {
                if (current_plane_checked != 0) { // If already checked at this height don't check again
                    continue;
                }
                height_to_check = current_building_height; // Check the front plane of the building
                current_plane_checked == 1; //Will check plane below
                change_in_height_var = 0.0;
            }
            else { // Height did change
                current_plane_checked = 0; //Reset that height to check
                if (current_DTM_height < current_building_height) { // If height decreases, no need to check again
                    current_building_height = current_DTM_height; //Set the next building as the new height
                    continue;
                }

                height_to_check = (current_building_height + current_DTM_height) / 2; //Interpolate between the heights
                change_in_height_var = pow(current_building_height - current_DTM_height, 2) / 2; // Sample variance
                current_building_height = current_DTM_height; //Set the next building as the new height
            }
        }


        // Check the height of the satellite based on the distance travelled along the satellite line
        
        // Calculate precise distance travelled for a height check. Only done here for efficiency
        // Get the distance in meters instead of the units of steps, required for slope uints
        line.d = sqrt(pow(line.E * DSM->step_size, 2) + pow(line.N * DSM->step_size, 2));
        
        sat_height = sat_vertical_slope * line.d + origin_height;

        // Change distance units back to step size instead of meters
        line.d = line.d / DSM->step_size; 


        if (DSM->processing_type == 1) {
            // if processing_type == 1 

            // Boolean rejection (old)
            // Satellite is lower than DTM, LOS is obstructed
            if (sat_height < current_DTM_height) {
                //fprintf(stderr, "sat height is less than DTM height, %lf < %lf\n", sat_height, current_DTM_height);
                return 1;
            }

            //Else, onto the new max height checked
            continue;
        }

        if (DSM->processing_type > 1 || DSM->processing_type < -1) {
            //Calculate probability of obstruction
            determine_sat_height_var(&sat_height_var, origin_horizontal_variance, origin_vertical_variance, sat_vertical_slope, DSM);

            determine_DTM_height_var(&DTM_height_var, DSM);
            
            y = (height_to_check - sat_height) / sqrt(DTM_height_var + change_in_height_var + sat_height_var); // Change in height var, part of the DTM height var

            apply_probability(&y, &probability_of_obstruction);

            if (debug) {
                fprintf(stderr, "Probability Updated To: %lf, Sat Height: %lf, DTM Height: %lf, Sat Height SD: %lf, DTM Height SD: %lf, y: %lf\n", probability_of_obstruction, sat_height, height_to_check, sqrt(sat_height_var + change_in_height_var), sqrt(DTM_height_var), y);
            }

            // If hit 0.99 round up to 1.
            if (probability_of_obstruction > 0.99) {
                probability_of_obstruction = 1.0f;
                if (debug) {
                    fprintf(stderr, "hit max probability %lf\n", probability_of_obstruction);
                }
                return probability_of_obstruction;
            }

            continue;
        }

        // Otherwise continue traversing


    }
}

void step_along_line(LineState* l) // credit to https://zingl.github.io/bresenham.html
{
    l->e2 = 2 * l->err;

    // Traverse North/South
    if (l->e2 >= l->dN) {
        l->err += l->dN;
        l->E += l->sE;
    }
    
    // Traverse East/West
    if (l->e2 <= l->dE) {
        l->err += l->dE;
        l->N += l->sN;
    }

    // If traversed East/West and North/South add root 2 distance
    if (l->e2 <= l->dE && l->e2 >= l->dN) {
        l->d += 1.41421356237f; //Root 2
    }
    // Otherwise only travelled one axis
    else {
        l->d += 1;
    }
}

void apply_probability(double* y, double* probability) {
    double prob = 0.5 * erfc(-*y / sqrt(2.0)); // CDF approximation

    // If probability has not yet been set
    if (*probability <= 0.0) { // Handles uninitialized and 0 case (though 0 technically handled farther down)
        *probability = prob;
        return;
    }

    //If probability has been set multiply the probabilities. Probability should go up given additional probabilities
    *probability = 1 - (1 - prob) * (1 - *probability);
}

void determine_DTM_height_var(double* var, struct DSMData* DSM) {
    //DTM possible error sources:
    // Height accuracy
    // Coordinate Horizontal accuracy - insignificant as all coordinates along the line are called. Can make the assumption that surfaces are generally flat so slight horizontal error doesn't matter
    // Height errros due to calling the nearest height - insignficant as all coordinates along the line are called. Can make the assumption that surfaces are generally flat so no interpolation error would be introduced anyways
    //*var = pow((double)DSM->step_size,2) + DSM->vertical_point_variance;
    *var = DSM->vertical_point_variance;
}

void determine_sat_height_var(double* var, double origin_horizontal_variance, double origin_vertical_variance, double sat_vertical_slope,DSMData* DSM) {
    //Formula for sat height variance: height = tan(elev) * horizontal distance traversed + origin_height
    //Possible Error Sources:
    // Elevation accuracy - function of sat err, pos err and baseline distance. Assumed to be insignificant due to extremely long baseline
    // Distance Traversed - euclidean distance of starting coordinate and traversed coordinate
    //  Starting Coordinate - Use the horizontal variance of the starting coordinate
    //  Traversed Coordinate - This coordinate can be incorrect due to not perfectly following the line. Bresenham is theoretically accurate to up 1/2 a pixel, arguably more accurate than that. 
    //      Standard deviation will be reflected as 0.68 * Step Size / 2. This is probably overly pessimistic.
    // Origin Height - Will have errors, depends on origin height source
    //  DEM Source - use determine_DTM_height_var + DEM offset var
    //  KF Source - use the estimated filter variance
    *var = origin_vertical_variance;
    double distance_var = pow(0.34 * DSM->step_size,2) + origin_horizontal_variance;
    *var = *var + pow(sat_vertical_slope, 2) * distance_var;

    //fprintf(stderr, "Sat Height Var Calculated Using: Origin Vertical Variance: %lf, Origin Horizontal Variance: %lf, Distance Variance: %lf, Calced Variance, %lf\n", origin_vertical_variance, origin_horizontal_variance, distance_var, *var);


    return;
}

/* solution to covariance - copy from solution.c-------------------*/
void soltocov_rtk(sol_t* sol, double* P)
{
    P[0] = sol->qr[0]; /* xx or ee */
    P[4] = sol->qr[1]; /* yy or nn */
    P[8] = sol->qr[2]; /* zz or uu */
    P[1] = P[3] = sol->qr[3]; /* xy or en */
    P[5] = P[7] = sol->qr[4]; /* yz or nu */        
    P[2] = P[6] = sol->qr[5]; /* zx or ue */
}