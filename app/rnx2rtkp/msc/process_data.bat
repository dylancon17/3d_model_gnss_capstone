@echo off
REM ===================================================================
REM Usage:
REM   process_data.bat [Dataset] [DEM] [PLOT] [PERFORMANCE] [ANALYZE] [PREFIX]
REM Example:
REM   process_data.bat                     -> ALL_DATA, BOTHDEM, NOPLOT, no prefix (default)
REM   process_data.bat 123 BOTHDEM NOPLOT NOPERFORMANCE NOANALYZE   -> datasets 1,2,3, dem and no dem processing, no rtkplotting, no performance analysis, no python script calls
REM   process_data.bat ALL_DATA NODEM PLOT PERFORMANCE ANALYZE RUN4 -> datasets 1-6, NODEM only, open RTKPLOT, include performance monitoring and python analysis, put run4 in all file names
REM ===================================================================

REM ---------------- Configuration (edit paths if necessary) -----------
set "ROOT=C:\capstone\ToShare"
set "DATAYEAR=23"
set "RTKLIB_EXE=C:\capstone\3d_model_gnss_capstone\app\rnx2rtkp\msc\Release\rnx2rtkp.exe"
set "RTKPLOT_EXE=C:\capstone\3d_model_gnss_capstone\app\rtkplot\rtkplot.exe"
set "CONFIG=C:\capstone\3d_model_gnss_capstone\app\rnx2rtkp\msc\config.conf"
set "wpaProfile=C:\capstone\3d_model_gnss_capstone\app\rnx2rtkp\msc\JustCPUrnx2rtkp.wpaProfile"
REM --------------------------------------------------------------------

REM Enable delayed expansion for variables modified inside loops
setlocal enabledelayedexpansion

REM ------------------ Parse arguments with defaults --------------------
REM Arg order: Dataset, DEM, PLOT, PREFIX
set "ARG_DATASET=%~1"
set "ARG_DEM=%~2"
set "ARG_PLOT=%~3"
set "ARG_PERFORMANCE=%~4"
set "ARG_ANALYZE=%~5"
set "ARG_PREFIX=%~6"

if "!ARG_DATASET!"=="" set "ARG_DATASET=ALL_DATA"
if "!ARG_DEM!"=="" set "ARG_DEM=BOTHDEM"
if "!ARG_PLOT!"=="" set "ARG_PLOT=NOPLOT"
if "!ARG_ANALYZE!"=="" set "ARG_ANALYZE=ANALYZE"
REM ARG_PREFIX default is empty
REM ARG_PERFORMANCE default is empty

REM ------------------ Build dataset list preserving order --------------
set "DATASET_LIST="

if /I "!ARG_DATASET!"=="ALL_DATA" (
    set "DATASET_LIST=1 2 3 4 5 6"
) else (
    REM Parse characters in ARG_DATASET and accept only digits 1..6 in the order they appear
    set "s=!ARG_DATASET!"
    set "i=0"
    :__char_loop
    call set "ch=%%s:~%i%,1%%"
    if "!ch!"=="" goto __char_done
    if "!ch!"=="1" (set "DATASET_LIST=!DATASET_LIST! 1")
    if "!ch!"=="2" (set "DATASET_LIST=!DATASET_LIST! 2")
    if "!ch!"=="3" (set "DATASET_LIST=!DATASET_LIST! 3")
    if "!ch!"=="4" (set "DATASET_LIST=!DATASET_LIST! 4")
    if "!ch!"=="5" (set "DATASET_LIST=!DATASET_LIST! 5")
    if "!ch!"=="6" (set "DATASET_LIST=!DATASET_LIST! 6")
    set /a i+=1
    goto __char_loop
    :__char_done
    REM Trim leading spaces
    for /f "tokens=* delims= " %%D in ("!DATASET_LIST!") do set "DATASET_LIST=%%D"
    REM Fallback to ALL_DATA if nothing was accepted
    if "!DATASET_LIST!"=="" set "DATASET_LIST=1 2 3 4 5 6"
)

REM ------------------ Interpret DEM argument ----------------------------
REM Accept DEM, NODEM, BOTHDEM (default BOTHDEM)
set "DO_DEM_1=0"   REM flag: run WITH -dem
set "DO_DEM_0=0"   REM flag: run WITHOUT -dem

if /I "!ARG_DEM!"=="DEM" (
    set "DO_DEM_1=1"
) else if /I "!ARG_DEM!"=="NODEM" (
    set "DO_DEM_0=1"
) else (
    REM BOTHDEM default
    set "DO_DEM_1=1"
    set "DO_DEM_0=1"
)

REM ------------------ Interpret PLOT argument ---------------------------
set "DO_PLOT=0"
if /I "!ARG_PLOT!"=="PLOT" (
    set "DO_PLOT=1"
)

REM ------------------ Generate one timestamp for entire run -----------
for /f %%a in ('wmic os get localdatetime ^| find "."') do set DTS=%%a
set TIMESTAMP=%DTS:~0,8%_%DTS:~8,6%



REM ---------------- Elevate if needed for performance monitoring --------------------
net session >nul 2>&1
if %errorlevel% neq 0 (
    if /I "%ARG_PERFORMANCE%"=="PERFORMANCE" (
        echo Requesting administrator privileges...
        powershell -Command "Start-Process '%~f0' -ArgumentList '%*' -Verb RunAs"
        exit /b
    )
)


REM ------------------ MAIN loop over datasets --------------------------
REM DATASET_LIST is space-separated and preserves user order
for %%D in (!DATASET_LIST!) do (
    set "SPECIFICDATASET=%%D"

    REM Build per-dataset OUTDIR: ROOT\<dataset>\output
    set "OUTDIR=%ROOT%\!SPECIFICDATASET!\output"

    REM Create output dir if it does not exist (silent)
    if not exist "!OUTDIR!" (
        mkdir "!OUTDIR!" >nul 2>&1
    )

    REM Build ROVER and BASE paths for this dataset (delayed expansion)
    set "ROVER=%ROOT%\!SPECIFICDATASET!\RINEXv3_04\!SPECIFICDATASET!.%DATAYEAR%"
    set "BASE=%ROOT%\!SPECIFICDATASET!\EAGLE_HIGHRATE\!SPECIFICDATASET!\RINEXv3_04\EAGLE_HIGHRATE.%DATAYEAR%"

    REM Prepare custom prefix (if provided)
    set "CUSTOM_PREFIX="
    if not "%ARG_PREFIX%"=="" (
        set "CUSTOM_PREFIX=%ARG_PREFIX%_"
    )

    REM Initialize variables to hold produced .pos paths
    set "OUTPATH_DEM="
    set "OUTPATH_NODEM="

    REM ------------------ DEM run (if enabled) ------------------------
    if "!DO_DEM_1!"=="1" (
        set "FINAL_PREFIX=DEM_!CUSTOM_PREFIX!"
        set "OUTFILE=solution_!TIMESTAMP!.pos"
        set "OUTPATH=!OUTDIR!\!FINAL_PREFIX!!OUTFILE!"
	echo Running RTKLIB. Dataset: !SPECIFICDATASET! DEM Flag: Disabled Output Name: !OUTPATH!
	if /I "!ARG_PERFORMANCE!"=="PERFORMANCE" (
            set "OUTPATHETL=!OUTDIR!\!FINAL_PREFIX!solution_!TIMESTAMP!.etl"
	    wpr -start CPU.Light -filemode
	)
        "%RTKLIB_EXE%" -k "%CONFIG%" -dem -o "!OUTPATH!" "!ROVER!O" "!BASE!O" "!ROVER!N" "!ROVER!G" "!ROVER!H" "!ROVER!J" "!ROVER!C" "!ROVER!Q" "!ROVER!P"
	if /I "!ARG_PERFORMANCE!"=="PERFORMANCE" (
            wpr -stop !OUTPATHETL!
            wpaexporter.exe -i !OUTPATHETL! -profile %wpaProfile% -outputfolder !OUTDIR!
            ren "!OUTDIR!\CPU_Usage_(Precise)_Utilization_by_Process,_Thread,_Stack.csv" "!FINAL_PREFIX!solution_!TIMESTAMP!.perf"
            del /Q !OUTPATHETL!
        )
	echo RTKLIB complete. Dataset: !SPECIFICDATASET! DEM Flag: Disabled Output Name: !OUTPATH!
        if "!DO_PLOT!"=="1" (
            start "" /min "%RTKPLOT_EXE%" "!OUTPATH!" >nul 2>&1
	)
    )


    REM ------------------ NODEM run (if enabled) -----------------------
    if "!DO_DEM_0!"=="1" (
        set "FINAL_PREFIX=!CUSTOM_PREFIX!"
        set "OUTFILE=solution_!TIMESTAMP!.pos"
        set "OUTPATH=!OUTDIR!\!FINAL_PREFIX!!OUTFILE!"
	echo Running RTKLIB. Dataset: !SPECIFICDATASET! DEM Flag: Disabled Output Name: !OUTPATH!
	if /I "!ARG_PERFORMANCE!"=="PERFORMANCE" (
            set "OUTPATHETL=!OUTDIR!\!FINAL_PREFIX!solution_!TIMESTAMP!.etl"
	    wpr -start CPU.Light -filemode
	)
        "%RTKLIB_EXE%" -k "%CONFIG%" -o "!OUTPATH!" "!ROVER!O" "!BASE!O" "!ROVER!N" "!ROVER!G" "!ROVER!H" "!ROVER!J" "!ROVER!C" "!ROVER!Q" "!ROVER!P"
	if /I "!ARG_PERFORMANCE!"=="PERFORMANCE" (
            wpr -stop !OUTPATHETL!
            wpaexporter.exe -i !OUTPATHETL! -profile %wpaProfile% -outputfolder !OUTDIR!
            ren "!OUTDIR!\CPU_Usage_(Precise)_Utilization_by_Process,_Thread,_Stack.csv" "!FINAL_PREFIX!solution_!TIMESTAMP!.perf"
            del /Q !OUTPATHETL!
        )
	echo RTKLIB complete. Dataset: !SPECIFICDATASET! DEM Flag: Disabled Output Name: !OUTPATH!
        if "!DO_PLOT!"=="1" (
            start "" /min "%RTKPLOT_EXE%" "!OUTPATH!" >nul 2>&1
	)
    )

    REM End of per-dataset processing
)

REM Clean up and exit
pause
endlocal
exit /b 0
