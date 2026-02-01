/*------------------------------------------------------------------------------
* rnx2rtkp.c : read rinex obs/nav files and compute receiver positions
*
*          Copyright (C) 2007-2009 by T.TAKASU, All rights reserved.
*
* version : $Revision: 1.1 $ $Date: 2008/07/17 21:55:16 $
* history : 2007/01/16  1.0 new
*           2007/03/15  1.1 add library mode
*           2007/05/08  1.2 separate from postpos.c
*           2009/01/20  1.3 support rtklib 2.2.0 api
*           2009/12/12  1.4 support glonass
*                           add option -h, -a, -l, -x
*           2010/01/28  1.5 add option -k
*           2010/08/12  1.6 add option -y implementation (2.4.0_p1)
*           2014/01/27  1.7 fix bug on default output time format
*-----------------------------------------------------------------------------*/
#include <stdarg.h>
#include "rtklib.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <direct.h>
#include <errno.h>

static const char rcsid[] = "$Id: rnx2rtkp.c,v 1.1 2008/07/17 21:55:16 ttaka Exp $";

#define PROGNAME    "rnx2rtkp"          /* program name */
#define MAXFILE     8                   /* max number of input files */

/* help text -----------------------------------------------------------------*/
static const char* help[] = {
"",
" usage: rnx2rtkp [option]... file file [...]",
"",
" Read RINEX OBS/NAV/GNAV/HNAV/CLK, SP3, SBAS message log files and ccompute ",
" receiver (rover) positions and output position solutions.",
" The first RINEX OBS file shall contain receiver (rover) observations. For the",
" relative mode, the second RINEX OBS file shall contain reference",
" (base station) receiver observations. At least one RINEX NAV/GNAV/HNAV",
" file shall be included in input files. To use SP3 precise ephemeris, specify",
" the path in the files. The extension of the SP3 file shall be .sp3 or .eph.",
" All of the input file paths can include wild-cards (*). To avoid command",
" line deployment of wild-cards, use \"...\" for paths with wild-cards.",
" Command line options are as follows ([]:default). With -k option, the",
" processing options are input from the configuration file. In this case,",
" command line options precede options in the configuration file.",
"",
" -?        print help",
" -k file   input options from configuration file [off]",
" -o file   set output file [stdout]",
" -ts ds ts start day/time (ds=y/m/d ts=h:m:s) [obs start time]",
" -te de te end day/time   (de=y/m/d te=h:m:s) [obs end time]",
" -ti tint  time interval (sec) [all]",
" -p mode   mode (0:single,1:dgps,2:kinematic,3:static,4:moving-base,",
"                 5:fixed,6:ppp-kinematic,7:ppp-static) [2]",
" -m mask   elevation mask angle (deg) [15]",
" -f freq   number of frequencies for relative mode (1:L1,2:L1+L2,3:L1+L2+L5) [2]",
" -v thres  validation threshold for integer ambiguity (0.0:no AR) [3.0]",
" -b        backward solutions [off]",
" -c        forward/backward combined solutions [off]",
" -i        instantaneous integer ambiguity resolution [off]",
" -h        fix and hold for integer ambiguity resolution [off]",
" -e        output x/y/z-ecef position [latitude/longitude/height]",
" -a        output e/n/u-baseline [latitude/longitude/height]",
" -n        output NMEA-0183 GGA sentence [off]",
" -g        output latitude/longitude in the form of ddd mm ss.ss' [ddd.ddd]",
" -t        output time in the form of yyyy/mm/dd hh:mm:ss.ss [sssss.ss]",
" -u        output time in utc [gpst]",
" -d col    number of decimals in time [3]",
" -s sep    field separator [' ']",
" -r x y z  reference (base) receiver ecef pos (m) [average of single pos]",
" -l lat lon hgt reference (base) receiver latitude/longitude/height (deg/m)",
" -y level  output soltion status (0:off,1:states,2:residuals) [0]",
" -x level  debug trace level (0:off) [0]",
" -dem      use a dem to aid in the solution output",
" -dopout        export DOP grid CSVs [off]",
" -dopstep sec   DOP export timestep in seconds [900]",
" -dopgrid m     DOP grid spacing in meters [500]"
};
/* show message --------------------------------------------------------------*/
extern int showmsg(char* format, ...)
{
    va_list arg;
    va_start(arg, format); vfprintf(stderr, format, arg); va_end(arg);
    fprintf(stderr, "\r");
    return 0;
}
extern void settspan(gtime_t ts, gtime_t te) {}
extern void settime(gtime_t time) {}

/* print help ----------------------------------------------------------------*/
static void printhelp(void)
{
    int i;
    for (i = 0;i < (int)(sizeof(help) / sizeof(*help));i++) fprintf(stderr, "%s\n", help[i]);
    exit(0);
}


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

static int en_to_llh(const prcopt_t* opt, double E, double N, double h, double pos_llh[3])
{
    double lat, lon;

    if (!tm_inv_EN_to_latlon_rad(E, N, &opt->UTM, &opt->ellip, &lat, &lon)) {
        return 0;
    }

    pos_llh[0] = lat;
    pos_llh[1] = lon;
    pos_llh[2] = h;
    return 1;
}

static int compute_dop_at_pos(
    const nav_t* nav,
    gtime_t t,
    const double rcv_ecef[3],
    double elmin_rad,
    int navsys,
    double dop_out[4],
    int* ns_out)
{
    double pos_llh[3];
    double azel[MAXSAT][2];
    int ns = 0;

    if (ns_out) *ns_out = 0;

    ecef2pos(rcv_ecef, pos_llh);

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
        satazel(pos_llh, e, azel_i);

        if (azel_i[1] < elmin_rad) continue;

        azel[ns][0] = azel_i[0];
        azel[ns][1] = azel_i[1];
        ns++;
    }

    if (ns < 4) return 0;

    dops(ns, azel, elmin_rad, dop_out);
    if (ns_out) *ns_out = ns;
    return 1;
}


/* rnx2rtkp main -------------------------------------------------------------*/
int main(int argc, char** argv)
{
    prcopt_t prcopt = prcopt_default;
    solopt_t solopt = solopt_default;
    filopt_t filopt = { "" };
    gtime_t ts = { 0 }, te = { 0 };
    double tint = 0.0, es[] = { 2000,1,1,0,0,0 }, ee[] = { 2000,12,31,23,59,59 }, pos[3];
    int i, j, n, ret;
    char* infile[MAXFILE], * outfile = "";

    prcopt.mode = PMODE_KINEMA;
    prcopt.navsys = SYS_GPS | SYS_GLO;
    prcopt.refpos = 1;
    prcopt.glomodear = 1;
    solopt.timef = 0;
    sprintf(solopt.prog, "%s ver.%s", PROGNAME, VER_RTKLIB);
    sprintf(filopt.trace, "%s.trace", PROGNAME);

    int dop_enable = 0;
    const char* dop_outdir = "C:\\capstone\\tmp";
    double dop_step_sec = 60.0 * 15.0;
    double dop_grid_m = 100.0;
    double dop_h_m = 1100.0;

    /* load options from configuration file */
    for (i = 1;i < argc;i++) {
        if (!strcmp(argv[i], "-k") && i + 1 < argc) {
            resetsysopts();
            if (!loadopts(argv[++i], sysopts)) return -1;
            getsysopts(&prcopt, &solopt, &filopt);
        }
    }

    prcopt.DSM.processing_type = 0; //Manually set DEM default
    prcopt.DSM.heights_array = NULL; //First set the heights array to NULL

    for (i = 1, n = 0;i < argc;i++) {
        if (!strcmp(argv[i], "-o") && i + 1 < argc) outfile = argv[++i];
        else if (!strcmp(argv[i], "-ts") && i + 2 < argc) {
            sscanf(argv[++i], "%lf/%lf/%lf", es, es + 1, es + 2);
            sscanf(argv[++i], "%lf:%lf:%lf", es + 3, es + 4, es + 5);
            ts = epoch2time(es);
        }
        else if (!strcmp(argv[i], "-te") && i + 2 < argc) {
            sscanf(argv[++i], "%lf/%lf/%lf", ee, ee + 1, ee + 2);
            sscanf(argv[++i], "%lf:%lf:%lf", ee + 3, ee + 4, ee + 5);
            te = epoch2time(ee);
        }
        else if (!strcmp(argv[i], "-ti") && i + 1 < argc) tint = atof(argv[++i]);
        else if (!strcmp(argv[i], "-k") && i + 1 < argc) { ++i; continue; }
        else if (!strcmp(argv[i], "-p") && i + 1 < argc) prcopt.mode = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-f") && i + 1 < argc) prcopt.nf = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-m") && i + 1 < argc) prcopt.elmin = atof(argv[++i]) * D2R;
        else if (!strcmp(argv[i], "-v") && i + 1 < argc) prcopt.thresar[0] = atof(argv[++i]);
        else if (!strcmp(argv[i], "-s") && i + 1 < argc) strcpy(solopt.sep, argv[++i]);
        else if (!strcmp(argv[i], "-d") && i + 1 < argc) solopt.timeu = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-dem") && i + 1 < argc) prcopt.DSM.processing_type = atoi(argv[++i]);

        /* -------------------- ADDED (for DOP export) -------------------- */
        else if (!strcmp(argv[i], "-dopout")) dop_enable = 1;
        else if (!strcmp(argv[i], "-dopstep") && i + 1 < argc) dop_step_sec = atof(argv[++i]);
        else if (!strcmp(argv[i], "-dopgrid") && i + 1 < argc) dop_grid_m = atof(argv[++i]);
        /* --------------------------------------------------------------- */

        else if (!strcmp(argv[i], "-b")) prcopt.soltype = 1;
        else if (!strcmp(argv[i], "-c")) prcopt.soltype = 2;
        else if (!strcmp(argv[i], "-i")) prcopt.modear = 2;
        else if (!strcmp(argv[i], "-h")) prcopt.modear = 3;
        else if (!strcmp(argv[i], "-t")) solopt.timef = 1;
        else if (!strcmp(argv[i], "-u")) solopt.times = TIMES_UTC;
        else if (!strcmp(argv[i], "-e")) solopt.posf = SOLF_XYZ;
        else if (!strcmp(argv[i], "-a")) solopt.posf = SOLF_ENU;
        else if (!strcmp(argv[i], "-n")) solopt.posf = SOLF_NMEA;
        else if (!strcmp(argv[i], "-g")) solopt.degf = 1;
        else if (!strcmp(argv[i], "-r") && i + 3 < argc) {
            prcopt.refpos = 0;
            for (j = 0;j < 3;j++) prcopt.rb[j] = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "-l") && i + 3 < argc) {
            prcopt.refpos = 0;
            for (j = 0;j < 3;j++) pos[j] = atof(argv[++i]);
            for (j = 0;j < 2;j++) pos[j] *= D2R;
            pos2ecef(pos, prcopt.rb);
        }
        else if (!strcmp(argv[i], "-y") && i + 1 < argc) solopt.sstat = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-x") && i + 1 < argc) solopt.trace = atoi(argv[++i]);
        else if (*argv[i] == '-') printhelp();
        else if (n < MAXFILE) infile[n++] = argv[i];
    }
    if (n <= 0) {
        showmsg("error : no input file");
        return -2;
    }

    /* Initialize DSM and related objects */
    if (prcopt.DSM.processing_type != 0) {
        /* set DEM options */
        prcopt.DSM.max_distance = 100; // Only search for buildings up to 2km away
        prcopt.DSM.antenna_dem_offset = 2; // Assuming antenna is 2m above DEM //TODO-DC, get a better number
        prcopt.DSM.use_dem_height_only = 0; //Use solved GNSS height as height origin for traverses
        prcopt.DSM.rejection_threshold = 0.9; // Reject sats with more than 50% probability of obstruction
        prcopt.DSM.antenna_dem_offset_var = 1; // 1m^2 variance in vehicle height
        prcopt.DSM.vertical_point_variance = pow(0.15, 2); //15cm accuracy
        prcopt.DSM.max_noise_scaling = 20; // Scale errors to a max of 20 times
        prcopt.DSM.building_height_margin = 1;

        //initialize_tiles_dataset(&prcopt.tiles_dataset, 7, 8, 5000, 5000, -22000.00, 5678000.00);
        //initialize_dsm_tile("C:\\capstone\\dsm_tiles\\DSM_CGY_5x5km_res1m\\DSM_CGY_5x5km_res1m_-7000E_5658000N.bin", &(prcopt.DSM), -7000.000, 5658000.000, 1, 5000);

        // initialize_tiles_dataset(&prcopt.tiles_dataset, 1, 1, 20000, 25000, -17000.00, 5678000.00);
        // initialize_dsm_tile("C:\\capstone\\dsm_tiles\\DTMCombinedNorth_567800_-17000_5N_4E.bin", &(prcopt.DSM), -17000.000, 5678000.000, 1, 20000);

        // initialize_tiles_dataset(&prcopt.tiles_dataset, 1, 1, 25000, 25000, -17000.00, 5673000.00);
        // initialize_dsm_tile("C:\\capstone\\dsm_tiles\\DTMCombinedNorth_5673000_-17000_5N_5E.bin", &(prcopt.DSM), -17000.000, 5673000.000, 1, 25000);

        /* Your tile */
        initialize_tiles_dataset(&prcopt.tiles_dataset, 1, 1, 25000, 25000, -15989.47, 5672949.28);
        initialize_dsm_tile("C:\\capstone\\dsm_tiles\\DTMCombinedNorth.bin", &(prcopt.DSM),
            -15989.47, 5672949.28, 1, 25000);


        // initialize_tiles_dataset(&prcopt.tiles_dataset, 1, 1, 35000, 45000, -22000.00, 5678000.00);
        // initialize_dsm_tile("C:\\capstone\\dsm_tiles\\DTMCombined.bin", &(prcopt.DSM), -22000.000, 5678000.000, 1, 35000);

        /* Calgary 114W 3TM */
        prcopt.UTM.central_meridian = -114.0;
        prcopt.UTM.latitude_of_origin = 0.0;
        prcopt.UTM.scale_factor = 0.9999;
        prcopt.UTM.false_easting = 0.0;
        prcopt.UTM.false_northing = 0.0;

        /* WGS 84 Ellipsoid */
        prcopt.ellip.flattening = 1 / 298.257223563;
        prcopt.ellip.semi_major_axis = 6378137.0;
        prcopt.ellip.semi_minor_axis = -1 * (prcopt.ellip.flattening * prcopt.ellip.semi_major_axis - prcopt.ellip.semi_major_axis);
        double a = prcopt.ellip.semi_major_axis;
        double b = prcopt.ellip.semi_minor_axis;
        prcopt.ellip.first_eccentricity = (a * a - b * b) / (a * a);
        prcopt.ellip.second_eccentricity = (a * a - b * b) / (b * b);
    }

    if (dop_enable)
    {
        obs_t obs = { 0 };
        nav_t nav0 = { 0 };
        sta_t sta0 = { 0 };

        if (n >= 1) {
            if (readrnx(infile[0], 1, "", &obs, &nav0, &sta0) && obs.n > 0) {

                gtime_t obs_start = obs.data[0].time;

                gtime_t utc = gpst2utc(obs_start);
                double ep0[6];
                time2epoch(utc, ep0);

                int year = (int)ep0[0];
                int month = (int)ep0[1];
                int day = (int)ep0[2];

                freeobs(&obs);
                freenav(&nav0, 0);

                {
                    const char* python_exe = "C:\\Users\\cronu\\anaconda3\\python.exe";
                    const char* script_path = "C:\\capstone\\3d_model_gnss_capstone\\analysis_scripts\\fetch_brdc.py";
                    const char* fetch_outdir = "C:\\capstone\\tmp";

                    if (run_python_fetch_brdc(python_exe, script_path, fetch_outdir, year, month, day)) {

                        nav_t nav = { 0 };

                        if (readrnx("C:\\capstone\\tmp\\brdc.rnx", 0, "", NULL, &nav, NULL)) {

                            gtime_t t0, t1;

                            if (ts.time != 0) {
                                double ep[6];
                                time2epoch(ts, ep);
                                ep[3] = 0; ep[4] = 0; ep[5] = 0;
                                t0 = utc2gpst(epoch2time(ep));
                            }
                            else {
                                double ep[6] = { (double)year, (double)month, (double)day, 0, 0, 0 };
                                t0 = utc2gpst(epoch2time(ep));
                            }

                            t1 = timeadd(t0, 86400.0);

                            const double E0 = -15989.47;
                            const double N0 = 5672949.28;
                            const double size_m = 25000.0;

                            ensure_dir_exists(dop_outdir);

                            for (gtime_t t = t0; timediff(t, t1) <= 0.0; t = timeadd(t, dop_step_sec)) {

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
                                if (!fp) continue;

                                fprintf(fp, "lat_deg,lon_deg,vdop,hdop,pdop,num_sats\n");

                                int point_count = 0;

                                for (double N = N0; N <= N0 + size_m; N += dop_grid_m) {
                                    for (double E = E0; E <= E0 + size_m; E += dop_grid_m) {

                                        double pos_llh[3];
                                        if (!en_to_llh(&prcopt, E, N, dop_h_m, pos_llh)) continue;

                                        double rcv_ecef[3];
                                        pos2ecef(pos_llh, rcv_ecef);

                                        double dop[4];
                                        int ns_used = 0;

                                        if (!compute_dop_at_pos(&nav, t, rcv_ecef, prcopt.elmin, prcopt.navsys, dop, &ns_used)) {
                                            continue;
                                        }

                                        fprintf(fp, "%.10f,%.10f,%.8f,%.8f,%.8f,%d\n",
                                            pos_llh[0] * R2D,
                                            pos_llh[1] * R2D,
                                            dop[3], dop[2], dop[1],
                                            ns_used);

                                        point_count++;
                                    }
                                }

                                fclose(fp);
                            }

                            freenav(&nav, 0);
                        }
                    }
                }
            }
            else {
                freeobs(&obs);
                freenav(&nav0, 0);
            }
        }
    }
    ret = postpos(ts, te, tint, 0.0, &prcopt, &solopt, &filopt, infile, n, outfile, "", "");

    free(prcopt.DSM.heights_array);


    if (!ret) fprintf(stderr, "%40s\r", "");
    return ret;
}
