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
" -dopout        export DOP grid to area [off]",
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


/* rnx2rtkp main -------------------------------------------------------------*/
int main(int argc, char** argv)
{
    printf("\n Start of rnx2rtkp main \n");
    prcopt_t prcopt = prcopt_default;
    solopt_t solopt = solopt_default;
    filopt_t filopt = { "" };
    gtime_t ts = { 0 }, te = { 0 };
    double tint = 0.0, es[] = { 2000,1,1,0,0,0 }, ee[] = { 2000,12,31,23,59,59 }, pos[3];
    int i, j, n, ret;
    char* infile[MAXFILE], * outfile = "";

    prcopt.mode = PMODE_KINEMA;
    prcopt.navsys = SYS_ALL;
    prcopt.refpos = 1;
    prcopt.glomodear = 1;
    solopt.timef = 0;
    sprintf(solopt.prog, "%s ver.%s", PROGNAME, VER_RTKLIB);
    sprintf(filopt.trace, "%s.trace", PROGNAME);

    int dop_enable = 0;
    int dsmopt = 1;
    const char* dop_outdir = "C:\\capstone\\tmp";
    double dop_step_sec = 60.0 * 15.0;
    double dop_grid_m = 1.0;

    /* load options from configuration file */
    for (i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "-k") && i + 1 < argc) {
            resetsysopts();
            if (!loadopts(argv[++i], sysopts)) return -1;
            getsysopts(&prcopt, &solopt, &filopt);
        }
    }

    prcopt.DSM.processing_type = 0; //Manually set DEM default
    prcopt.DSM.heights_array = NULL; //First set the heights array to NULL

    for (i = 1, n = 0; i < argc; i++) {
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

        else if (!strcmp(argv[i], "-dopout") && i + 1 < argc) dop_enable = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-dsmopt") && i + 1 < argc) dsmopt = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-dopstep") && i + 1 < argc) dop_step_sec = atof(argv[++i]);
        else if (!strcmp(argv[i], "-dopgrid") && i + 1 < argc) dop_grid_m = atof(argv[++i]);

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
            for (j = 0; j < 3; j++) prcopt.rb[j] = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "-l") && i + 3 < argc) {
            prcopt.refpos = 0;
            for (j = 0; j < 3; j++) pos[j] = atof(argv[++i]);
            for (j = 0; j < 2; j++) pos[j] *= D2R;
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

    prcopt.navsys = SYS_ALL;

    /* Initialize DSM and related objects */
    if (prcopt.DSM.processing_type != 0) {
        /* set DEM options */
        prcopt.DSM.max_distance = 500; // Infinite search distance (throttled wtihin code)
        prcopt.DSM.antenna_dem_offset = 1.1; // Assuming antenna is 2m above DEM //TODO-DC, get a better number
        prcopt.DSM.use_dem_height_only = 0; //Use solved GNSS height as height origin for traverses
        prcopt.DSM.rejection_threshold = 0.9; // Reject sats with more than 50% probability of obstruction
        prcopt.DSM.antenna_dem_offset_var = 0.1; // 1m^2 variance in vehicle height
        prcopt.DSM.vertical_point_variance = pow(0.15, 2); //15cm accuracy 
        prcopt.DSM.max_noise_scaling = 100; // Scale errors to a max of 20 times
        prcopt.DSM.building_height_margin = 1; //Unused
      
        //initialize_tiles_dataset(&prcopt.tiles_dataset, 7, 8, 5000, 5000, -22000.00, 5678000.00);
        //initialize_dsm_tile("C:\\capstone\\dsm_tiles\\DSM_CGY_5x5km_res1m\\DSM_CGY_5x5km_res1m_-7000E_5658000N.bin", &(prcopt.DSM), -7000.000, 5658000.000, 1, 5000);

        // initialize_tiles_dataset(&prcopt.tiles_dataset, 1, 1, 20000, 25000, -17000.00, 5678000.00);
        // initialize_dsm_tile("C:\\capstone\\dsm_tiles\\DTMCombinedNorth_567800_-17000_5N_4E.bin", &(prcopt.DSM), -17000.000, 5678000.000, 1, 20000);

        // initialize_tiles_dataset(&prcopt.tiles_dataset, 1, 1, 25000, 25000, -17000.00, 5673000.00);
        // initialize_dsm_tile("C:\\capstone\\dsm_tiles\\DTMCombinedNorth_5673000_-17000_5N_5E.bin", &(prcopt.DSM), -17000.000, 5673000.000, 1, 25000);

        // initialize_tiles_dataset(&prcopt.tiles_dataset, 1, 1, 25000, 25000, -15989.47, 5672949.28);
        // initialize_dsm_tile("C:\\capstone\\dsm_tiles\\DTMCombinedNorth.bin", &(prcopt.DSM),
        //    -15989.47, 5672949.28, 1, 25000);

        initialize_tiles_dataset(&prcopt.tiles_dataset, 1, 1, 25000, 25000, -15989.47, 5672949.28);
        initialize_dsm_tile("C:\\capstone\\dsm_tiles\\DTMCombinedNorth_5672949.28_-15989.47_5N_5E.bin", &(prcopt.DSM), -15989.47, 5672949.28, 1, 25000);

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

    if (dop_enable > 0) {
        // As all other heights are meaningless
        prcopt.DSM.use_dem_height_only = 1;

        dop_traverse traverse;

        if (dop_enable == 1) { // Downtown
            traverse.N0 = 5656200.28;
            traverse.E0 = -6675.47;
            traverse.length_N = 2000.0;
            traverse.length_E = 3500.0;
            traverse.dop_grid_m = 2.0;
            traverse.dop_step_sec = 30.0 * 60.0;
        }
        else {
            if (dop_enable == 2) { // University
                traverse.N0 = 5660700.28;
                traverse.E0 = -10789.47;
                traverse.length_N = 1000.0;
                traverse.length_E = 2500.0;
                traverse.dop_grid_m = 5.0;
                traverse.dop_step_sec = dop_step_sec;
            }
            else {
                if (dop_enable == 3) { // Calgary
                    traverse.N0 = 5670449.28;
                    traverse.E0 = -13489.47;
                    traverse.length_N = 20000.0;
                    traverse.length_E = 20000.0;
                    traverse.dop_grid_m = 100.0;
                    traverse.dop_step_sec = dop_step_sec;
                }
                else {
                    fprintf(stderr, "Unrecognized area. Crashing program");
                    free(prcopt.DSM.heights_array);
                    return 1;
                }
            }
        }

        dop_csv(&prcopt, ts, n, infile, dop_outdir, &traverse);
    }
    else {
        // No point in processing if our goal is just to generate the dop graphs
        ret = postpos(ts, te, tint, 0.0, &prcopt, &solopt, &filopt, infile, n, outfile, "", "");
    }

    free(prcopt.DSM.heights_array);


    if (!ret) fprintf(stderr, "%40s\r", "");
    return ret;
}
