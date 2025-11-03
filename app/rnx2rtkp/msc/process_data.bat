@echo off
REM ===================================================================
REM Usage:
REM   process_data.bat [Dataset] [DEM] [PLOT] [PREFIX]
REM Example:
REM   process.bat                     -> ALL_DATA, BOTHDEM, NOPLOT, no prefix (default)
REM   process.bat 3 DEM PLOT RUN3     -> dataset 3, DEM only, plot, prefix RUN3
REM   process.bat 123 BOTHDEM NOPLOT  -> datasets 1,2,3 in order, both DEM & NODEM, no plot
REM   process.bat ALL_DATA NODEM PLOT -> datasets 1-6, NODEM only, open RTKPLOT minimized/non-blocking
REM   TODO add in argument related to auto building the code
REM   TODO add in arugment for measuring resources (check out vsperfcmd)
REM   TODO add in argument for running the python script
REM ===================================================================

REM ---------------- Configuration (edit paths if necessary) -----------
set "ROOT=C:\capstone\ToShare"
set "DATAYEAR=23"
set "RTKLIB_EXE=C:\capstone\3d_model_gnss_capstone\app\rnx2rtkp\msc\Release\rnx2rtkp.exe"
set "RTKPLOT_EXE=C:\capstone\3d_model_gnss_capstone\app\rtkplot\rtkplot.exe"
set "CONFIG=C:\capstone\3d_model_gnss_capstone\app\rnx2rtkp\msc\config.conf"
REM --------------------------------------------------------------------

REM Enable delayed expansion for variables modified inside loops
setlocal enabledelayedexpansion

REM ------------------ Parse arguments with defaults --------------------
REM Arg order: Dataset, DEM, PLOT, PREFIX
set "ARG_DATASET=%~1"
set "ARG_DEM=%~2"
set "ARG_PLOT=%~3"
set "ARG_PREFIX=%~4"

if "%ARG_DATASET%"=="" set "ARG_DATASET=ALL_DATA"
if "%ARG_DEM%"=="" set "ARG_DEM=BOTHDEM"
if "%ARG_PLOT%"=="" set "ARG_PLOT=NOPLOT"
REM ARG_PREFIX default is empty

REM ------------------ Build dataset list preserving order --------------
set "DATASET_LIST="

if /I "%ARG_DATASET%"=="ALL_DATA" (
    set "DATASET_LIST=1 2 3 4 5 6"
) else (
    REM Parse characters in ARG_DATASET and accept only digits 1..6 in the order they appear
    set "s=%ARG_DATASET%"
    set "i=0"
    :__char_loop
    call set "ch=%%s:~%i%,1%%"
    if "%ch%"=="" goto __char_done
    if "%ch%"=="1" (set "DATASET_LIST=!DATASET_LIST! 1")
    if "%ch%"=="2" (set "DATASET_LIST=!DATASET_LIST! 2")
    if "%ch%"=="3" (set "DATASET_LIST=!DATASET_LIST! 3")
    if "%ch%"=="4" (set "DATASET_LIST=!DATASET_LIST! 4")
    if "%ch%"=="5" (set "DATASET_LIST=!DATASET_LIST! 5")
    if "%ch%"=="6" (set "DATASET_LIST=!DATASET_LIST! 6")
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

if /I "%ARG_DEM%"=="DEM" (
    set "DO_DEM_1=1"
) else if /I "%ARG_DEM%"=="NODEM" (
    set "DO_DEM_0=1"
) else (
    REM BOTHDEM default
    set "DO_DEM_1=1"
    set "DO_DEM_0=1"
)

REM ------------------ Interpret PLOT argument ---------------------------
set "DO_PLOT=0"
if /I "%ARG_PLOT%"=="PLOT" (
    set "DO_PLOT=1"
)

REM ------------------ Generate one timestamp for entire run -----------
REM Format: YYYYMMDD_HHMMSS
for /f "tokens=1-6 delims=/:. " %%a in ("%date% %time%") do (
    set "MM=%%a"
    set "DD=%%b"
    set "YYYY=%%c"
    set "HH=%%d"
    set "MIN=%%e"
    set "SEC=%%f"
)

REM Pad single-digit values with leading zero
if "!HH!"==" " set "HH=00"
if "!MIN!"==" " set "MIN=00"
if "!SEC!"==" " set "SEC=00"

REM Final timestamp
set "TIMESTAMP=%YYYY%%MM%%DD%_%HH%%MIN%%SEC%"


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
        set "DEMFLAG=-dem"
        set "DEM_PREFIX=DEM_"
        set "FINAL_PREFIX=!CUSTOM_PREFIX!!DEM_PREFIX!"
        set "OUTFILE=solution_!TIMESTAMP!.pos"
        set "OUTPATH=!OUTDIR!\!FINAL_PREFIX!!OUTFILE!"
        REM Run rnx2rtkp with -dem, silent. Remove >nul 2>&1 if you want console output.
	echo Running RTKLIB. Dataset: !SPECIFICDATASET! DEM Flag: Enabled Output Name: !OUTPATH!
        "%RTKLIB_EXE%" -k "%CONFIG%" !DEMFLAG! -o "!OUTPATH!" "!ROVER!O" "!BASE!O" "!ROVER!N" "!ROVER!G" "!ROVER!H" "!ROVER!J" "!ROVER!C" "!ROVER!Q" "!ROVER!P"
	echo RTKLIB Complete. Dataset: !SPECIFICDATASET! DEM Flag: Enabled Output Name: !OUTPATH!
        set "OUTPATH_DEM=!OUTPATH!"
    )

    REM ------------------ NODEM run (if enabled) -----------------------
    if "!DO_DEM_0!"=="1" (
        set "FINAL_PREFIX=!CUSTOM_PREFIX!"
        set "OUTFILE=solution_!TIMESTAMP!.pos"
        set "OUTPATH=!OUTDIR!\!FINAL_PREFIX!!OUTFILE!"
	echo Running RTKLIB. Dataset: !SPECIFICDATASET! DEM Flag: Disabled Output Name: !OUTPATH!
        "%RTKLIB_EXE%" -k "%CONFIG%" -o "!OUTPATH!" "!ROVER!O" "!BASE!O" "!ROVER!N" "!ROVER!G" "!ROVER!H" "!ROVER!J" "!ROVER!C" "!ROVER!Q" "!ROVER!P"
	echo RTKLIB complete. Dataset: !SPECIFICDATASET! DEM Flag: Disabled Output Name: !OUTPATH!
        set "OUTPATH_NODEM=!OUTPATH!"
    )

    REM ------------------ RTKPLOT launching logic ----------------------
    REM Only call RTKPLOT if DO_PLOT == 1 (user asked for PLOT)
    REM For BOTHDEM case (both OUTPATH_DEM and OUTPATH_NODEM defined), we launch ONE RTKPLOT
    REM with NODEM first, DEM second so DEM overlays.
    if "!DO_PLOT!"=="1" (
        REM BOTH available: launch combined (NODEM then DEM)
        if defined OUTPATH_DEM if defined OUTPATH_NODEM (
            REM Use start "" /min to launch minimized and non-blocking
            start "" /min "%RTKPLOT_EXE%" "!OUTPATH_NODEM!" "!OUTPATH_DEM!" >nul 2>&1
        ) else (
            REM Only DEM exists
            if defined OUTPATH_DEM (
                start "" /min "%RTKPLOT_EXE%" "!OUTPATH_DEM!" >nul 2>&1
            )
            REM Only NODEM exists
            if defined OUTPATH_NODEM (
                start "" /min "%RTKPLOT_EXE%" "!OUTPATH_NODEM!" >nul 2>&1
            )
        )
    )

    REM End of per-dataset processing
)

REM Clean up and exit
pause
endlocal
exit /b 0
