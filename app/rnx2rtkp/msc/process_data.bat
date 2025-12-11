@echo off

REM ------------------------ Help Menu -------------------------
if /I "%~1"=="-h" goto :show_help
if /I "%~1"=="--help" goto :show_help

REM ---------------- Configuration (edit paths if necessary) -----------
set "ROOT=C:\capstone\ToShare"
set "DATAYEAR=23"
set "RTKLIB_EXE=C:\capstone\3d_model_gnss_capstone\app\rnx2rtkp\msc\Release\rnx2rtkp_vc.exe"
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
if "!ARG_DEM!"=="" set "ARG_DEM=0"
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


    set "FINAL_PREFIX=!CUSTOM_PREFIX!_!ARG_DEM!"
    set "OUTFILE=solution_!TIMESTAMP!.pos"
    set "OUTPATH=!OUTDIR!\!FINAL_PREFIX!!OUTFILE!"
    echo Running RTKLIB. Dataset: !SPECIFICDATASET! DEM Flag: !ARG_DEM! Output Name: !OUTPATH!
    if /I "!ARG_PERFORMANCE!"=="PERFORMANCE" (
        set "OUTPATHETL=!OUTDIR!\!FINAL_PREFIX!solution_!TIMESTAMP!.etl"
        wpr -start CPU.Light -filemode
    )
    "%RTKLIB_EXE%" -k "%CONFIG%" -dem !ARG_DEM! -o "!OUTPATH!" "!ROVER!O" "!BASE!O" "!ROVER!N" "!ROVER!G" "!ROVER!H" "!ROVER!J" "!ROVER!C" "!ROVER!Q" "!ROVER!P"
    if /I "!ARG_PERFORMANCE!"=="PERFORMANCE" (
        wpr -stop !OUTPATHETL!
        wpaexporter.exe -i !OUTPATHETL! -profile %wpaProfile% -outputfolder !OUTDIR!
        ren "!OUTDIR!\CPU_Usage_(Precise)_Utilization_by_Process,_Thread,_Stack.csv" "!FINAL_PREFIX!solution_!TIMESTAMP!.perf"
        del /Q !OUTPATHETL!
    )
    echo RTKLIB complete. Dataset: !SPECIFICDATASET! DEM Flag: !ARG_DEM! Output Name: !OUTPATH!
    if "!DO_PLOT!"=="1" (
        start "" /min "%RTKPLOT_EXE%" "!OUTPATH!" >nul 2>&1
    )

    REM End of per-dataset processing
)

REM Clean up and exit
pause
endlocal
exit /b 0

:show_help
echo.
echo ===================================================================
echo Usage:
echo   process_data.bat [Dataset] [DEM_FLAG] [PLOT] [PERFORMANCE] [ANALYZE] [PREFIX]
echo.
echo DEM Options:
echo    0 = no DEM processing
echo    1 = boolean observation rejection
echo    2 = probability-threshold observation rejection
echo    3 = deweighting + rejection (probability threshold)
echo    4 = deweighting only
echo    5 = deweighting + rejection + reference satellite selection
echo.
echo PLOT Options:
echo    PLOT       = open RTKLIB plots
echo    NOPLOT     = do not open RTKLIB plots
echo.
echo PERFORMANCE Options:
echo    PERFORMANCE   = run CPU monitoring
echo    NOPERFORMANCE = do not run CPU monitoring
echo.
echo ANALYZE Options:
echo    ANALYZE    = call python analysis script
echo    NOANALYZE  = do not call python analysis script
echo.
echo PREFIX:
echo    Added to all output file names
echo.
echo Examples:
echo   process_data.bat 1 2
echo        - Runs datasets 1,2,3 with DEM=2, no plot, no perf, no analysis
echo.
echo   process_data.bat 123 1 NOPLOT NOPERFORMANCE NOANALYZE
echo        - Runs datasets 1,2,3 with DEM=1, no plot/perf/analysis
echo.
echo   process_data.bat ALL_DATA 2 PLOT PERFORMANCE ANALYZE RUN4
echo        - Runs datasets 1 to 6 with full processing and prefix RUN4
echo ===================================================================
echo.
goto :EOF

