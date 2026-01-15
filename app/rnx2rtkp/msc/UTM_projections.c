#include "rtklib.h"
#define _USE_MATH_DEFINES
#include <math.h>

void project_latitude_longitude_to_UTM
(
	east_north* EN,
	const lat_long* ll,
	const UTM_projection* proj,
	const ellipsoid* e
)
{

	const double a = e->semi_major_axis;
	const double lat = ll->latitude * M_PI / 180;
	const double lon = ll->longitude * M_PI / 180;
	const double W = sqrt(1 - e->first_eccentricity * sin(lat) * sin(lat));
	const double N = a / W;
	const double eta_squared = e->second_eccentricity * cos(lat) * cos(lat);
	const double delta_lon = lon - proj->central_meridian * M_PI / 180;

	/* Series expansion parameters for y-coordinate projection S(lat) --------------------------
	   A, B, C, D, and E are ellipsoid-specific coefficients */

	const double ee = e->first_eccentricity;

	/* Check: A = 1 - 1/4 x e^2 - 3/64 x e^4 - 5/256 x e^6 - 175/16384 x e^8 */
	const double A = 1.0 - 1.0 / 4.0 * ee - 3.0 / 64.0 * ee * ee - 5.0 / 256.0 * ee * ee * ee - 175.0 / 16384.0 * ee * ee * ee * ee;

	/* Check: B = 3/8 x (e^2 + 1/4 x e^4 + 15/128 x e^6 - 455/4096 x e^8) */
	const double B = 3.0 / 8.0 * (ee + 1.0 / 4.0 * ee * ee + 15.0 / 128.0 * ee * ee * ee - 455.0 / 4096.0 * ee * ee * ee * ee);

	/* Check: C = 15/256 x (e^4 + 3/4 x e^6 - 77/128 x e^8) */
	const double C = 15.0 / 256.0 * (ee * ee + 3.0 / 4.0 * ee * ee * ee - 77.0 / 128.0 * ee * ee * ee * ee);

	/* Check: D = 35/3072 x (e^6 - 41/32 x e^8) */
	const double D = 35.0 / 3072.0 * (ee * ee * ee - 41.0 / 32.0 * ee * ee * ee * ee);

	/* Check: E = -315/131072 x e^8 */
	const double E = -315.0 / 131072.0 * ee * ee * ee * ee;

	/* Check: S(lat) = a[A x lat - B x sin(2 x lat) + C x sin(4 x lat) - D x sin(6 x lat) + E x sin(8 x lat)] */
	const double s_lat = a * (A * lat - B * sin(2.0 * lat) + C * sin(4.0 * lat) - D * sin(6.0 * lat) + E * sin(8.0 * lat));

	/* Series expansion parameters for Gauss-Kruger projection -----------------------------------
	   Note that the powers are turned into several multiplications because it is faster for the program than using pow() */

	const double A1 = cos(lat);
	const double A2 = 1.0 / 2.0 * tan(lat) * cos(lat) * cos(lat);
	const double A3 = 1.0 / 6.0 * cos(lat) * cos(lat) * cos(lat) * (1.0 - tan(lat) * tan(lat) + eta_squared);

	/* Check: A4 = 1/24 x tan(lat) x cos(lat)^4 x (5 - tan(lat)^2 + 9 x eta_squared + 4 eta_squared^2 */
	const double A4 = 1.0 / 24.0 * tan(lat) * cos(lat) * cos(lat) * cos(lat) * cos(lat) * (5.0 - tan(lat) * tan(lat) + 9.0 * eta_squared + 4.0 * eta_squared * eta_squared);

	/* Check: A5 = 1/120 x cos(lat)^5 x (5 - 18 x tan(lat)^2 + tan(lat)^4 + 14 x eta_squared - 58 x eta_squared x tan(lat)^2 */
	const double A5 = 1.0 / 120.0 * cos(lat) * cos(lat) * cos(lat) * cos(lat) * cos(lat) * (5.0 - 18.0 * tan(lat) * tan(lat) + tan(lat) * tan(lat) * tan(lat) * tan(lat) + 14.0 * eta_squared - 58.0 * eta_squared * tan(lat) * tan(lat));

	/* Check: A6 = 1/720 x tan(lat) x cos(lat)^6 x (61 - 58 x tan(lat)^2 + tan(lat)^4 + 270 x eta_squared - 330 x eta_squared x tan(lat)^2)*/
	const double A6 = 1.0 / 720.0 * tan(lat) * cos(lat) * cos(lat) * cos(lat) * cos(lat) * cos(lat) * cos(lat) * (61.0 - 58.0 * tan(lat) * tan(lat) + tan(lat) * tan(lat) * tan(lat) * tan(lat) + 270.0 * eta_squared - 330 * eta_squared * tan(lat) * tan(lat));

	/* Check: x = N x (A1 x dlambda^1 + A3 x dlambda^3 + A5 x dlambda^5) */
	const double x = N * (A1 * delta_lon + A3 * delta_lon * delta_lon * delta_lon + A5 * delta_lon * delta_lon * delta_lon * delta_lon * delta_lon);

	/* Check: y = s_lat + N x (A2 x dlambda^2 + A4 x dlambda^4 + A6 x dlambda^6) */
	const double y = s_lat + N * (A2 * delta_lon * delta_lon + A4 * delta_lon * delta_lon * delta_lon * delta_lon + A6 * delta_lon * delta_lon * delta_lon * delta_lon * delta_lon * delta_lon);

	EN->easting = proj->false_easting + proj->scale_factor * x;
	EN->northing = proj->false_northing + proj->scale_factor * y;
}