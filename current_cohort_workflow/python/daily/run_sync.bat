@echo off
setlocal
cd /d "%~dp0"
set "VENV=%~dp0.venv"
set "PY=%VENV%\Scripts\python.exe"

if not exist "%PY%" goto setup
"%PY%" -c "import numpy, imageio_ffmpeg, pandas, matplotlib" >nul 2>&1
if errorlevel 1 goto update
goto run

:update
echo [setup] Installing missing packages...
"%PY%" -m pip install --disable-pip-version-check -q -r "%~dp0requirements.txt" >nul 2>&1
"%PY%" -c "import numpy, imageio_ffmpeg, pandas, matplotlib" >nul 2>&1
if errorlevel 1 goto rebuild
goto run

:rebuild
echo [setup] Python environment is broken or from another PC - rebuilding...
rmdir /s /q "%VENV%"

:setup
echo [setup] First run on this PC: creating Python environment, please wait...
py -3 -m venv "%VENV%" >nul 2>&1
if not exist "%PY%" python -m venv "%VENV%"
if not exist "%PY%" goto nopython
"%PY%" -m pip install --disable-pip-version-check -q -r "%~dp0requirements.txt"
if errorlevel 1 goto pipfail

:run
if /i "%~1"=="daily" goto daily
echo ===== Step 1: sync Python + Saleae + video =====
"%PY%" "%~dp0sync_sessions.py" %*
set "RC=%errorlevel%"
echo.
if "%RC%"=="1" echo Some sessions FAILED or no data folders were found - see the messages above and sync_summary.csv.
if "%RC%"=="2" echo Could not write sync_summary.csv - close it in Excel and run again.
echo.
echo ===== Step 2: quick analysis (tables + figures) =====
"%PY%" "%~dp0quick_analysis.py"
if errorlevel 1 echo Quick analysis failed - see the messages above. If a CSV/PNG is open in another program, close it and run again.
echo.
echo Done. Sync: sync folder.  Analysis: analysis\days\dayN\figs\00_day_dashboard.png.  All sessions: session_log.csv.
pause
exit /b %RC%

:daily
shift
echo ===== End of day: arrange today's files, sync, analysis =====
"%PY%" "%~dp0daily_summary.py" %1 %2 %3 %4
set "RC=%errorlevel%"
echo.
if not "%RC%"=="0" echo Something needs attention - read the messages above (CHECK / FAIL lines).
echo Open: analysis\days\dayN\figs\00_day_dashboard.png and analysis\days\dayN\dayN_summary.txt
pause
exit /b %RC%

:nopython
echo ERROR: Python 3 is not installed. Install it from python.org - tick "Add python.exe to PATH" - then run this again.
pause
exit /b 1

:pipfail
echo ERROR: could not install the Python packages (numpy, imageio-ffmpeg, pandas, matplotlib). Check the internet connection and run this again.
pause
exit /b 1
