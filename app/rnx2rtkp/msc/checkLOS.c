#define _USE_MATH_DEFINES

#include "rtklib.h"

#include <stdlib.h>

#include <math.h>



/* ---------- Floating-point DDA (Amanatides & Woo) data structure ---------- */
typedef struct {
    // Origin in grid units (grid coord = meters / step_size)
    double x0, y0;

    // Direction unit vector in grid units (norm not required because we use t in grid units)
    double dx, dy;

    // Current integer cell index
    int ix, iy;

    // Step signs (+1 or -1)
    int stepX, stepY;

    // Parametric distances (in grid units) to next vertical/horizontal grid boundary
    double tMaxX, tMaxY;

    // Parametric distance (in grid units) between successive vertical/horizontal crossings
    double tDeltaX, tDeltaY;

    // last and current t (parametric distance along ray in grid units)
    double tPrev, t;
} DDARay;


/* ---------- Helper prototypes ---------- */
void init_dda_ray(DDARay* r, double azimuth, double origin_x_grid, double origin_y_grid);
void step_dda_ray(DDARay* r);
double compute_cell_weight(const DDARay* r, int cell_ix, int cell_iy, struct DSMData* DSM);
void determine_DTM_height_var(double* var, struct DSMData* DSM);
void determine_sat_height_var(double* var, double origin_horizontal_variance, double origin_vertical_variance, double sat_vertical_slope, struct DSMData* DTM);
void soltocov_rtk(sol_t* sol, double* P);
double phi_from_standardized(double y); // safe phi wrapper
void reset_scaling(rtk_t* rtk, int* ns, int* sat);


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
extern int los_update(rtk_t* rtk, const obsd_t* obs, int* sat, int* iu, int* ir, int* ns, const double* rs) {
    //fprintf(stderr, "%s", "LOS Update\n");

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
        reset_scaling(rtk, ns, sat);
        return 0;
    }




    int i, nrej=0;
    int rej_idx[MAXOBS];
    const double *spp_rr = rtk->sol.rr; 

    // Converts ECEF to LLH
    double spp_pos[3];
    ecef2pos(spp_rr, spp_pos);
    
    // Converts ECEF covars to LLH
    double spp_P[9]; // 3x3 ENU
    double spp_Q[9];
    
    soltocov_rtk(&(rtk->sol), spp_P);
    covenu(spp_pos, spp_P, spp_Q);

    double kf_rr[3] = { rtk->x[0], rtk->x[1], rtk->x[2] };
    double kf_P[9]; // 3x3 ECEF covariance
    double kf_Q[9];

    double kf_pos[3];
    ecef2pos(kf_rr, kf_pos);

    for (int i = 0; i < 3; i++) 
        for (int j = 0; j < 3; j++) 
            kf_P[i + 3 * j] = rtk->P[i + j * rtk->nx];

    covenu(kf_pos, kf_P, kf_Q);

    lat_long kf_ll = { kf_pos[0] * 180 / M_PI, kf_pos[1] * 180 / M_PI};

    int out_of_bounds = 0;
    set_relative_origin(&(rtk->opt.DSM),&(rtk->opt.tiles_dataset), &kf_ll, &(rtk->opt.UTM), &(rtk->opt.ellip), &out_of_bounds);

    if (out_of_bounds == 1) { //If origin is out of bounds, don't search farther than that. Edge case that won't happen
        reset_scaling(rtk, ns, sat);
        return 0;
    }

    double current_DTM_height;
    int origin_x = 0, origin_y = 0;
    double dummy_distance = 0.0;
    get_relative_height(&(rtk->opt.DSM), &(rtk->opt.tiles_dataset), &origin_x, &origin_y, &dummy_distance, &current_DTM_height, &out_of_bounds);

    if (out_of_bounds == 1) { //If origin is out of bounds, don't search farther than that. Edge case that won't happen
        fprintf(stderr, "Theoretically impossible out of bounds hit");
        reset_scaling(rtk, ns, sat);
        return 0;
    }


    if (abs(current_DTM_height + rtk->opt.DSM.antenna_dem_offset - kf_pos[2]) > 5 && (rtk->opt.DSM.processing_type > 7 || rtk->opt.DSM.processing_type == 6)) {
        // Heights don't match. Incorrect starting height and therefore position
        //Try the single point position instead
        lat_long spp_ll = { spp_pos[0] * 180 / M_PI, spp_pos[1] * 180 / M_PI };

        set_relative_origin(&(rtk->opt.DSM), &(rtk->opt.tiles_dataset), &spp_ll, &(rtk->opt.UTM), &(rtk->opt.ellip), &out_of_bounds);

        if (out_of_bounds == 1) { //If origin is out of bounds, don't search farther than that
            reset_scaling(rtk, ns, sat);
            return 0;
        }

        get_relative_height(&(rtk->opt.DSM), &(rtk->opt.tiles_dataset), &origin_x, &origin_y, &dummy_distance, &current_DTM_height, &out_of_bounds);

        if (out_of_bounds == 1) { //If origin is out of bounds, don't search farther than that
            fprintf(stderr, "Theoretically impossible out of bounds hit");
            reset_scaling(rtk, ns, sat);
            return 0;
        }

        if (abs(current_DTM_height + rtk->opt.DSM.antenna_dem_offset - spp_pos[2]) > 5 && (rtk->opt.DSM.processing_type > 7 || rtk->opt.DSM.processing_type == 6)) {
            //Both positions are incorrect. Just raise the variance and deal with it.
            //reset_scaling(rtk, ns, sat);
            //return 0;
            kf_Q[8] = pow(abs(current_DTM_height - spp_pos[2]), 2);
        }

        /* KF height is wrong → fall back to SPP position and covariance */
        memcpy(kf_rr, spp_rr, sizeof(double) * 3);
        memcpy(kf_pos, spp_pos, sizeof(double) * 3);
        memcpy(kf_P, spp_P, sizeof(double) * 9);
        memcpy(kf_Q, spp_Q, sizeof(double) * 9);
    }

    double traverse_origin_relative_m_x = rtk->opt.DSM.relative_origin_traverse_true.easting - rtk->opt.DSM.relative_origin_traverse.easting;
    double traverse_origin_relative_m_y = rtk->opt.DSM.relative_origin_traverse_true.northing - rtk->opt.DSM.relative_origin_traverse.northing;

    double traverse_origin_relative_grid_x = traverse_origin_relative_m_x / rtk->opt.DSM.step_size;
    double traverse_origin_relative_grid_y = traverse_origin_relative_m_y / rtk->opt.DSM.step_size;


    double e[3], azel[2]; //warning gets overwritten each satellite. Should be fine?
    double r;
    double probability_of_obstruction = 0;


    for (i = 0;i < *ns && i < MAXOBS;i++) {


        //fprintf(stderr, "Checking satellite %d\n", i);

        r = geodist(rs + i * 6, kf_rr, e); //TODO how is rs indexed

            /* geodist failure check */
        if (r <= 0) {
            /* Bad geometry → reject satellite */
            rej_idx[nrej++] = i;
            // Reset scaling just in case...
            rtk->ssat[sat[i] - 1].obstruction_scaling = 1;
            //fprintf(stderr, "Geodist Failed %d\n", i);

            fprintf(stderr, "Geodist Check Failed");

            continue;
        }

        satazel(spp_pos, e, azel);
        
        if (debug) {
            fprintf(stderr, "Requesting probability: %d\n", sat[i]);
        }
        probability_of_obstruction = check_los(azel[0], azel[1], kf_pos[0], kf_pos[1], kf_pos[2], kf_Q[0] + kf_Q[4], max(kf_Q[8],25), &(rtk->opt.DSM), &(rtk->opt.tiles_dataset), traverse_origin_relative_grid_x, traverse_origin_relative_grid_y, debug);
        if (debug) {
            fprintf(stderr, "%lf\n", probability_of_obstruction);
        }
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
            probability_of_obstruction = 0.0;
        }

        if (debug) {
            fprintf(stderr, "Probability after editing, %lf", probability_of_obstruction);
        }

        double scaling = 1.0 / (1.0 - probability_of_obstruction);
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
        rtk->ssat[sat[i] - 1].height_offset = current_DTM_height - kf_pos[2];


    }

    //fprintf(stdout, "%d possible sats\n", *ns);

    // If deterministic or probabilistic rejection
    if (rtk->opt.DSM.processing_type > 0 && rtk->opt.DSM.processing_type != 4 && rtk->opt.DSM.processing_type != 9) {
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

void reset_scaling(rtk_t* rtk, int* ns, int* sat) {
    // Undo any scaling
    for (int i = 0;i < *ns && i < MAXOBS;i++) {
        rtk->ssat[sat[i] - 1].obstruction_scaling = 1.0;
    }
    return;
}

//Assumes relative origin already set
extern double check_los(double sat_az, double sat_elev, double origin_lat, double origin_long, double origin_height, double origin_horizontal_variance, double origin_vertical_variance, struct DSMData* DSM, TilesDataset* tiles_dataset, double traverse_origin_x_grid, double traverse_origin_y_grid, int debug) {
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
    double current_DTM_height = 0;
    double sat_vertical_slope = tan(sat_elev);

    // get starting grid coordinates and DTM height
    int origin_x = 0, origin_y = 0;
    double dummy_distance = 0.0;
    get_relative_height(DSM, tiles_dataset, &origin_x, &origin_y, &dummy_distance, &current_DTM_height, &out_of_bounds);

    if (debug) {
        fprintf(stderr, "Origin Height: %lf, DEM Height: %lf, Out of Bounds %d\n",
            origin_height,
            current_DTM_height,
            out_of_bounds);
    }


    if (out_of_bounds) {
        return -1.0; // out of DSM area
    }


    if (out_of_bounds != 1 && (current_DTM_height > origin_height || DSM->use_dem_height_only == 1) && DSM->use_dem_height_only != 2) {
        origin_height = current_DTM_height + DSM->antenna_dem_offset;
        origin_vertical_variance = DSM->vertical_point_variance + DSM->antenna_dem_offset_var;
        if (debug) {
            fprintf(stderr, ", Using DTM Height\n");
        }
        if (debug) fprintf(stderr, "Using DTM Height as origin height (with antenna offset)\n");
    }
    else {
        if (debug) {
            fprintf(stderr, ", Using GPS Height\n");
        }
    }
    // Only search until a distance where you're guranteed to hit a building, or you've searched a ways
    double max_distance_m = (DSM->max_dsm_height - origin_height) / sat_vertical_slope;
    if (max_distance_m > DSM->max_distance) {
        max_distance_m = DSM->max_distance;
    }

    // convert to number of steps (grid units)
    double max_distance_grid_units = max_distance_m / DSM->step_size;

    // end point in integer grid (approx)
    int E1 = (int)round(max_distance_grid_units * sin(sat_az));
    int N1 = (int)round(max_distance_grid_units * cos(sat_az));

    DDARay ray;
    init_dda_ray(&ray, sat_az, traverse_origin_x_grid, traverse_origin_y_grid);

    // accumulation variables for the three options
    double weighted_numer = 0.0; // numerator = sum(w_i * p_i)
    double weighted_denom = 0.0; // denom = sum(w_i)

    double ema_prob = 0.0;       // for correlation EMA
    int ema_initialized = 0;

    double y_max = -INFINITY;    // for max-exceedance

    double probability_of_obstruction = -1.0; // Uninitialized representation

    // Loop over DDA traversal
    while (1) {
        // move to next cell
        step_dda_ray(&ray);

        //fprintf(stderr, "---DTM Height: %lf, dE: %d: dN %d, out_of_bounds: %d---", current_DTM_height, line.E, line.N, out_of_bounds);
        // if we've exceeded max_distance_grid_units stop as original did
        double d_grid = (ray.t + ray.tPrev) / 2;
        
        if (d_grid > max_distance_grid_units) {
            if (debug) fprintf(stderr, "max distance travelled\n");

            return probability_of_obstruction;
        }

        get_relative_height(DSM, tiles_dataset, &ray.ix, &ray.iy, &d_grid, &current_DTM_height, &out_of_bounds);

        if (out_of_bounds) {
            // If we run off the DSM, behavior: LOS clear beyond DSM
            if (debug) fprintf(stderr, "out_of_bounds during traversal\n");
            return probability_of_obstruction;
        }

        // compute sat height at this distance (meters)
        double sat_height = origin_height + sat_vertical_slope * d_grid * DSM->step_size;

        if (debug) {
            fprintf(stderr, "Stepped along ray to: %d, %d, %lf with sat height: %lf, building height: %lf\n", ray.ix, ray.iy, ray.t, sat_height, current_DTM_height);
        }

        // Option-specific handling:
        // For deterministic processing type 1 (boolean), return 1 if DTM is high
        if (DSM->processing_type == 1) {
            // if processing_type == 1 

            // Boolean rejection (old)
            // Satellite is lower than DTM, LOS is obstructed
            if (sat_height < current_DTM_height) {
                if (debug) {
                    fprintf(stderr, "sat height is less than DTM height, %lf < %lf, exiting\n", sat_height, current_DTM_height);
                }
                return 1.0;
            }
            else {
                probability_of_obstruction = 0.0;
                // else continue
                continue;
            }
        }

        // compute sat height variance
        double sat_height_var = 0.0;
        determine_sat_height_var(&sat_height_var, origin_horizontal_variance, origin_vertical_variance, sat_vertical_slope, DSM);

        // compute DTM height variance
        double DTM_height_var = 0.0;
        determine_DTM_height_var(&DTM_height_var, DSM);

        // compute y and local probability p_i
        double denom_var = DTM_height_var + sat_height_var;
        double y = (current_DTM_height - sat_height) / sqrt(denom_var);
        double p_i_0 = phi_from_standardized(y);

        // compute weight for this cell and scale prob(distance_from_cell_center)
        double dw = compute_cell_weight(&ray, ray.ix, ray.iy, DSM);
        double p_i = p_i_0 * dw;

        

        if (debug) {
            fprintf(stderr, "Checking prob for: %lf, with original prob of %lf and deweighting of %lf using vars of %lf, %lf\n", p_i, p_i_0, dw, sat_height_var, DTM_height_var);
        }

        if (DSM->processing_type >= 7 || DSM->processing_type <= -4) {
            if (p_i > probability_of_obstruction) {
                probability_of_obstruction = p_i;
                if (debug) {
                    fprintf(stderr, "Updating prob to new max\n");
                }
            }
        }
        else {
            if (probability_of_obstruction <= 0.0) {
                probability_of_obstruction = p_i * 0.5; // Set max gain per cell when combining cells
            }

            probability_of_obstruction = 1.0 - (1.0 - probability_of_obstruction) * (1.0 - p_i);
            if (debug) {
                fprintf(stderr, "Combined probability into %lf\n", probability_of_obstruction);
            }
        }
        // Save over processing the data
        if (probability_of_obstruction > 0.99) {
            return 1.0;
        }
    }
} // end traversal loop

void init_dda_ray(DDARay* r, double azimuth, double origin_x_grid, double origin_y_grid)
{
    r->x0 = origin_x_grid;
    r->y0 = origin_y_grid;

    r->dx = sin(azimuth);
    r->dy = cos(azimuth);

    // Map fractional origin to the integer cell whose center is nearest
    const double EPS = 1e-14;
    r->ix = (int)floor(r->x0 + 0.5 + EPS);
    r->iy = (int)floor(r->y0 + 0.5 + EPS);

    r->stepX = (r->dx >= 0.0) ? 1 : -1;
    r->stepY = (r->dy >= 0.0) ? 1 : -1;

    // distance from origin to the cell boundaries at ix +/- 0.5
    double nextGridDistX = (r->stepX > 0) ? ((r->ix + 0.5) - r->x0) : (r->x0 - (r->ix - 0.5));
    double nextGridDistY = (r->stepY > 0) ? ((r->iy + 0.5) - r->y0) : (r->y0 - (r->iy - 0.5));

    if (nextGridDistX < 0.0) nextGridDistX = 0.0; // guard numerical noise
    if (nextGridDistY < 0.0) nextGridDistY = 0.0;

    if (fabs(r->dx) < 1e-12) {
        r->tMaxX = INFINITY;
        r->tDeltaX = INFINITY;
    }
    else {
        r->tMaxX = nextGridDistX / fabs(r->dx);
        r->tDeltaX = fabs(1.0 / r->dx);
    }

    if (fabs(r->dy) < 1e-12) {
        r->tMaxY = INFINITY;
        r->tDeltaY = INFINITY;
    }
    else {
        r->tMaxY = nextGridDistY / fabs(r->dy);
        r->tDeltaY = fabs(1.0 / r->dy);
    }

    r->t = 0.0;
    r->tPrev = 0.0;
}



/* ---------- Step the DDA ray to the next crossed grid boundary ---------- */
void step_dda_ray(DDARay* r)
{
    r->tPrev = r->t;

    if (r->tMaxX < r->tMaxY) {
        r->t = r->tMaxX;
        r->tMaxX += r->tDeltaX;
        r->ix += r->stepX;
    }
    else {
        r->t = r->tMaxY;
        r->tMaxY += r->tDeltaY;
        r->iy += r->stepY;
    }
}

/* ---------- Safe phi (CDF) from standardized y using erfc with guards ---------- */
double phi_from_standardized(double y)
{
    // For large magnitudes, clamp to 0 or 1 for numeric safety
    if (y > 8.0) return 1.0;
    if (y < -8.0) return 0.0;
    return 0.5 * erfc(-y / M_SQRT2);
}

/* ---------- compute DTM (DSM) height variance ---------- */
void determine_DTM_height_var(double* var, struct DSMData* DSM) {
    //DTM possible error sources:
    // Height accuracy
    // Coordinate Horizontal accuracy - insignificant as all coordinates along the line are called. Can make the assumption that surfaces are generally flat so slight horizontal error doesn't matter
    // Height errros due to calling the nearest height - insignficant as all coordinates along the line are called. Can make the assumption that surfaces are generally flat so no interpolation error would be introduced anyways
    // Keep a simple model: DSM point vertical variance
    *var = DSM->vertical_point_variance;
}

void determine_sat_height_var(double* var, double origin_horizontal_variance, double origin_vertical_variance, double sat_vertical_slope, DSMData* DSM) {
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
    // double distance_var = pow(0.34 * DSM->step_size, 2) + origin_horizontal_variance;
    // *var = *var + pow(sat_vertical_slope, 2) * distance_var;
    *var = origin_vertical_variance;
    //fprintf(stderr, "Sat Height Var Calculated Using: Origin Vertical Variance: %lf, Origin Horizontal Variance: %lf, Distance Variance: %lf, Calced Variance, %lf\n", origin_vertical_variance, origin_horizontal_variance, distance_var, *var);
    return;
}


/* ---------- compute per-cell weight for traversal over cell centers ---------- */
double compute_cell_weight(const DDARay* r, int cell_ix, int cell_iy, struct DSMData* DSM)
{
    // dx, dy in grid units
    double vx = cell_ix - r->x0;
    double vy = cell_iy - r->y0;

    // perp distance from ray to cell center (grid units)
    double dist_perp_grid = fabs(vx * r->dy - vy * r->dx); // 2D cross product magnitude

    double dist_perp_m = dist_perp_grid * DSM->step_size;

    // footprint sigma: assume uniform cell -> variance step_size^2 / 12
    double sigma_perp = DSM->step_size / sqrt(12.0);

    // Gaussian footprint weight (meters)
    double w_foot = exp(-(dist_perp_m * dist_perp_m) / (2.0 * sigma_perp * sigma_perp));

    return w_foot;
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