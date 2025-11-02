@echo on
set "ROOT=C:\capstone\ToShare"
set "SPECIFICDATASET=1"
set "DATAYEAR=23"
set "BASECONVERSIONTIME=Converted_On_20251022_T1349"
set "ROVER=%ROOT%\%SPECIFICDATASET%\RINEXv3_04\%SPECIFICDATASET%.%DATAYEAR%"
set "BASE=%ROOT%\%SPECIFICDATASET%\EAGLE_HIGHRATE\%BASECONVERSIONTIME%\RINEXv3_04\EAGLE_HIGHRATE.%DATAYEAR%"
set "RTKLIB_EXE=C:\capstone\3d_model_gnss_capstone\app\rnx2rtkp\msc\Release\rnx2rtkp.exe"
set "RTKPLOT_EXE=C:\capstone\3d_model_gnss_capstone\app\rtkplot\rtkplot.exe"

set "OUTDIR=%ROOT%\output"
if not exist "%OUTDIR%" (
    echo Output directory not found. Creating: %OUTDIR%
    mkdir "%OUTDIR%"
)

set "CONFIG=C:\capstone\3d_model_gnss_capstone\app\rnx2rtkp\msc\config.conf"

REM === Generate timestamp: YYYYMMDD_HHMMSS ===
for /f "tokens=1-5 delims=/: " %%a in ("%date% %time%") do (
    set YYYY=%%c
    set MM=%%a
    set DD=%%b
    set HH=%%d
    set MIN=%%e
)

REM Remove seconds + invalid characters
set "HH=%HH: =0%"
set "TIMESTAMP=%YYYY%%MM%%DD%_%HH%%MIN%"
set "OUT=%OUTDIR%\solution_%TIMESTAMP%.pos"

%RTKLIB_EXE% -k %CONFIG% -o "%OUT%" "%ROVER%O" "%BASE%O" "%ROVER%N" "%ROVER%G" "%ROVER%H" "%ROVER%J" "%ROVER%C" "%ROVER%Q" "%ROVER%P"

%RTKPLOT_EXE% "%OUT%"

pause
