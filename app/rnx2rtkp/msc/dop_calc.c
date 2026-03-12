#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include <direct.h>
#include <errno.h>
#include "rtklib.h"

static void ensure_dir_exists(const char* path)
{
    if (_mkdir(path) == 0) return;
    if (errno == EEXIST) return;
}

static int tm_inv_EN_to_latlon_rad(
    double E, double N,
    const UTM_projection* proj,
    const ellipsoid* e,
    double* lat_rad,
    double* lon_rad)
{
    const double a = e->semi_major_axis;
    const double f = e->flattening;

    if (a <= 0.0 || f <= 0.0) return 0;

    const double k0 = proj->scale_factor;
    const double FE = proj->false_easting;
    const double FN = proj->false_northing;

    const double lon0 = proj->central_meridian * D2R;
    const double lat0 = proj->latitude_of_origin * D2R;

    const double b = a * (1.0 - f);
    const double e2 = (a * a - b * b) / (a * a);
    const double ep2 = (a * a - b * b) / (b * b);

    const double x = E - FE;
    const double y = N - FN;

    const double M = y / k0;

    const double mu = M / (a * (1.0 - e2 / 4.0 - 3.0 * e2 * e2 / 64.0 - 5.0 * e2 * e2 * e2 / 256.0));

    const double e1 = (1.0 - sqrt(1.0 - e2)) / (1.0 + sqrt(1.0 - e2));
    const double J1 = (3.0 * e1 / 2.0 - 27.0 * e1 * e1 * e1 / 32.0);
    const double J2 = (21.0 * e1 * e1 / 16.0 - 55.0 * e1 * e1 * e1 * e1 / 32.0);
    const double J3 = (151.0 * e1 * e1 * e1 / 96.0);
    const double J4 = (1097.0 * e1 * e1 * e1 * e1 / 512.0);

    const double phi1 = mu
        + J1 * sin(2.0 * mu)
        + J2 * sin(4.0 * mu)
        + J3 * sin(6.0 * mu)
        + J4 * sin(8.0 * mu);

    const double sinphi1 = sin(phi1);
    const double cosphi1 = cos(phi1);
    const double tanphi1 = tan(phi1);

    const double N1 = a / sqrt(1.0 - e2 * sinphi1 * sinphi1);
    const double R1 = a * (1.0 - e2) / pow(1.0 - e2 * sinphi1 * sinphi1, 1.5);
    const double T1 = tanphi1 * tanphi1;
    const double C1 = ep2 * cosphi1 * cosphi1;
    const double D = x / (N1 * k0);

    double lat = phi1 - (N1 * tanphi1 / R1) * (
        (D * D) / 2.0
        - (5.0 + 3.0 * T1 + 10.0 * C1 - 4.0 * C1 * C1 - 9.0 * ep2) * pow(D, 4) / 24.0
        + (61.0 + 90.0 * T1 + 298.0 * C1 + 45.0 * T1 * T1 - 252.0 * ep2 - 3.0 * C1 * C1) * pow(D, 6) / 720.0
        );

    double lon = lon0 + (
        D
        - (1.0 + 2.0 * T1 + C1) * pow(D, 3) / 6.0
        + (5.0 - 2.0 * C1 + 28.0 * T1 - 3.0 * C1 * C1 + 8.0 * ep2 + 24.0 * T1 * T1) * pow(D, 5) / 120.0
        ) / cosphi1;

    lat += (lat0 - 0.0);

    *lat_rad = lat;
    *lon_rad = lon;
    return 1;
}

static int run_python_fetch_brdc(
    const char* python_exe,
    const char* script_path,
    const char* outdir,
    int year, int month, int day)
{
    char cmd[2048];

    snprintf(cmd, sizeof(cmd),
        "cmd.exe /V:OFF /C \"\"%s\" \"%s\" "
        "--outdir \"%s\" --year %d --month %d --day %d "
        "1>\"%s\\fetch_brdc.log\" 2>&1\"",
        python_exe, script_path, outdir,
        year, month, day, outdir);

    return system(cmd) == 0;
}

static int en_to_ll(const prcopt_t* opt, double E, double N, double pos_llh[3])
{
    double lat, lon;

    if (!tm_inv_EN_to_latlon_rad(E, N, &opt->UTM, &opt->ellip, &lat, &lon)) {
        return 0;
    }

    pos_llh[0] = lat;
    pos_llh[1] = lon;
    pos_llh[2] = 0.0;
    return 1;
}

static int compute_dop_at_pos(
    const nav_t* nav,
    gtime_t t,
    const double rcv_ecef[3],
    lat_long *ll,
    double elmin_rad,
    int navsys,
    double dop_out[4],
    int* ns_out,
    prcopt_t* popt)
{
    double azel[MAXSAT][2];
    int ns = 0;

    if (ns_out) *ns_out = 0;

    double traverse_origin_relative_m_x = popt->DSM.relative_origin_traverse_true.easting - popt->DSM.relative_origin_traverse.easting;
    double traverse_origin_relative_m_y = popt->DSM.relative_origin_traverse_true.northing - popt->DSM.relative_origin_traverse.northing;

    double traverse_origin_relative_grid_x = traverse_origin_relative_m_x / popt->DSM.step_size;
    double traverse_origin_relative_grid_y = traverse_origin_relative_m_y / popt->DSM.step_size;

    for (int sat = 1; sat <= MAXSAT; sat++) {

        int sys = satsys(sat, NULL);
        if (!(sys & navsys)) continue;

        double rs[6] = { 0 };
        double dts[2] = { 0 };
        double var = 0.0;
        int svh = 0;

        if (!satpos(t, t, sat, 0, nav, rs, dts, &var, &svh)) continue;
        if (svh) continue;

        double e[3];
        double r = geodist(rs, rcv_ecef, e);
        if (r <= 0.0) continue;

        double azel_i[2];
        satazel(ll, e, azel_i);

        if (azel_i[1] < elmin_rad) continue;

        double obstruction_prob = check_los(
            azel_i[0],
            azel_i[1],
            ll->latitude, // Technically fine, but bad practice
            ll->longitude, // Technically fine, but bad practice
            0.0, // Ignored as long as DSM->use_dem_height_only is 1
            2.0, // This variance value has been used and tuned against truth positions
            1.0, // This variance value has been used and tuned against truth positions
            &(popt->DSM),
            &(popt->tiles_dataset),
            traverse_origin_relative_grid_x,
            traverse_origin_relative_grid_y,
            0);

        if (obstruction_prob > 0.5 || obstruction_prob < 0) {
            continue;
        }

        azel[ns][0] = azel_i[0];
        azel[ns][1] = azel_i[1];
        ns++;
    }

    if (ns < 4) return 1;

    dops(ns, azel, elmin_rad, dop_out);
    if (ns_out) *ns_out = ns;
    return 1;
}

int dop_csv(prcopt_t* prcopt,
    gtime_t ts,
    int n,
    char** infile,
    const char* dop_outdir,
    double dop_step_sec,
    double dop_grid_m)
{

    fprintf(stderr, "In dop_csv\n");
    obs_t obs = { 0 };
    nav_t nav0 = { 0 };
    sta_t sta0 = { 0 };

    if (n < 1) return 0;
    if (!readrnx(infile[0], 1, "", &obs, &nav0, &sta0) || obs.n <= 0) {
        freeobs(&obs);
        freenav(&nav0, 0);
        return 0;
    }

    fprintf(stderr, "Read RINEX\n");


    gtime_t obs_start = obs.data[0].time;

    gtime_t utc = gpst2utc(obs_start);
    double ep0[6];
    time2epoch(utc, ep0);

    int year = (int)ep0[0];
    int month = (int)ep0[1];
    int day = (int)ep0[2];

    // fprintf(stderr, "Read Obs Time Successfully\n");


    freeobs(&obs);
    freenav(&nav0, 0);

    const char* python_exe = "python.exe";
    const char* script_path = "C:\\capstone\\3d_model_gnss_capstone\\analysis_scripts\\fetch_brdc.py";
    const char* fetch_outdir = "C:\\capstone\\tmp";

    if (!run_python_fetch_brdc(python_exe, script_path, fetch_outdir, year, month, day)) {
        fprintf(stderr, "Python Fetch BRDC Failed\n");

        return 0;
    }

    nav_t nav = { 0 };
    if (!readrnx("C:\\capstone\\tmp\\brdc.rnx", 0, "", NULL, &nav, NULL)) {
        fprintf(stderr, "Read BRDC Failed\n");
        freenav(&nav, 0);
        return 0;
    }

    gtime_t t0, t1;

    if (ts.time != 0) {
        double ep[6];
        time2epoch(ts, ep);
        ep[3] = 0; ep[4] = 0; ep[5] = 0;
        t0 = utc2gpst(epoch2time(ep));
    }
    else {
        double ep[6] = { (double)year,(double)month,(double)day,0,0,0 };
        t0 = utc2gpst(epoch2time(ep));
    }

    t1 = timeadd(t0, 86400.0);

    // fprintf(stderr, "Set End Time\n");


    const double E0 = -6675.47;
    const double N0 = 5656200.28;
    const double length_E = 2000.0;
    const double length_N = 3500.0;


    ensure_dir_exists(dop_outdir);

    // fprintf(stderr, "Past Dir Exists. Starting at Time %lf %lf with Time Step of: %lf\n", t0.sec, t0.time, dop_step_sec);


    for (gtime_t t = t0; timediff(t, t1) <= 0.0; t = timeadd(t, dop_step_sec)) {
        
        double ep[6];
        time2epoch(gpst2utc(t), ep);
        fprintf(stderr, "%04.0f-%02.0f-%02.0f %02.0f:%02.0f\n",
            ep[0], ep[1], ep[2], ep[3], ep[4]);

        gtime_t utc_t = gpst2utc(t);
        double epu[6];
        time2epoch(utc_t, epu);

        int yyyy = (int)epu[0];
        int mm = (int)epu[1];
        int dd = (int)epu[2];
        int HH = (int)epu[3];
        int MN = (int)epu[4];

        char day_folder[64];
        snprintf(day_folder, sizeof(day_folder), "%04d-%02d-%02d", yyyy, mm, dd);

        char out_folder[1024];
        snprintf(out_folder, sizeof(out_folder), "%s\\%s", dop_outdir, day_folder);
         ensure_dir_exists(out_folder);

        char out_file[1024];
        snprintf(out_file, sizeof(out_file),
            "%s\\%04d-%02d-%02d-%02d-%02d.csv",
            out_folder, yyyy, mm, dd, HH, MN);

        FILE* fp = fopen(out_file, "w");
        if (!fp) {
            fprintf(stderr, "Failed on fp\n");
            continue;
        }
        fprintf(fp, "lat_deg,lon_deg,vdop,hdop,pdop,num_sats,E,N\n");

        for (double NN = N0; NN >= N0 - length_N; NN -= dop_grid_m) {
            for (double EE = E0; EE <= E0 + length_E; EE += dop_grid_m) {

                double pos_llh[3];
                if (!en_to_ll(prcopt, EE, NN, pos_llh)) {
                    fprintf(stderr, "East North to LL failed\n");
                    continue;
                }

                lat_long ll = { pos_llh[0] * 180 / M_PI, pos_llh[1] * 180 / M_PI };

                int out_of_bounds = 0;

                set_relative_origin(&(prcopt->DSM), &(prcopt->tiles_dataset), &ll, &(prcopt->UTM), &(prcopt->ellip), &out_of_bounds);

                if (out_of_bounds == 1) { //If origin is out of bounds
                    fprintf(stderr, "Searching in the wrong area. Bad Configuration\n");
                    continue;
                }

                double current_DTM_height;
                int origin_x = 0, origin_y = 0;
                double dummy_distance = 0.0;
                get_relative_height(&(prcopt->DSM), &(prcopt->tiles_dataset), &origin_x, &origin_y, &dummy_distance, &current_DTM_height, &out_of_bounds);

                // Asssume that the antenna is 1m off the ground
                current_DTM_height = current_DTM_height + 1.0;

                if (out_of_bounds == 1) { //If origin is out of bounds, don't search farther than that. Edge case that won't happen
                    fprintf(stderr, "Theoretically impossible out of bounds hit in dop_calc\n");
                    continue;
                }

                pos_llh[2] = current_DTM_height;

                double rcv_ecef[3];
                pos2ecef(pos_llh, rcv_ecef);

                double dop[4];
                dop[0] = -1.0;
                dop[1] = -1.0;
                dop[2] = -1.0;
                dop[3] = -1.0;

                int ns_used = 0;

                if (!compute_dop_at_pos(&nav, t, rcv_ecef, &ll,
                    prcopt->elmin, prcopt->navsys,
                    dop, &ns_used, prcopt)) {
                    fprintf(stderr, "Compute DOP at POS failed\n");

                    continue;
                }

                fprintf(fp, "%.10f,%.10f,%.8f,%.8f,%.8f,%d,%.3f,%.3f\n",
                    ll.latitude,
                    ll.longitude,
                    dop[3], dop[2], dop[1],
                    ns_used,
                    EE,
                    NN);
            }
        }

        fclose(fp);
    }

    freenav(&nav, 0);
    return 1;
}
