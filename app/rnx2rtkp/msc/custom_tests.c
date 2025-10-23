#include "rtklib.h"

boolean test_los_case(
	float sat_az,
	float sat_elev,
	float origin_lat,
	float origin_long,
	float origin_height,
	DTMData* DTM,
	boolean expected_result,
	boolean* success);

// TODO-DS as part of setting up project flow, create a trigger to run tests
extern boolean test_los_summary(DTMData* DTM) {
	boolean success = 1;

	// TODO-TC create actual test cases
	test_los_case(0.0f, 0.5f, 10.0f, 10.0f, 100.0f, DTM, 1, success);
	test_los_case(1.57f, 0.1f, 10.0f, 10.0f, 50.0f, DTM, 0, success);

	return success;
}

boolean test_los_case(
	float sat_az,
	float sat_elev,
	float origin_lat,
	float origin_long,
	float origin_height,
	DTMData* DTM,
	boolean expected_result,
	boolean* success) {

	boolean result = check_los(
		sat_az,
		sat_elev,
		origin_lat,
		origin_long,
		origin_height,
		DTM
	);

	if (result != expected_result) {
		printf("Test failed: sat_az=%f, sat_elev=%f, origin_lat=%f, origin_long=%f, origin_height=%f, expected_result=%d, got_result=%d\n",
			sat_az, sat_elev, origin_lat, origin_long, origin_height, expected_result, result);
		success = 0;
		return 0; // Test failed
	}
	else return 1;
}


