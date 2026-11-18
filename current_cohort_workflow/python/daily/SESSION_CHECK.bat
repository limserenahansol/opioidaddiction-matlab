@echo off
rem Read-only check right after each session (IC Capture stopped, Saleae capture saved and exported to CSV).
rem   SESSION_CHECK.bat                               newest finished session of today; when the newest log is a
rem                                                   test run (name with 'test', e.g. the morning valve test before
rem                                                   the first mouse) it is checked too: valve time vs the day's valve
rem   SESSION_CHECK.bat --all                         every session of today
rem   SESSION_CHECK.bat --session 9669mblack_prday3   one session (or test run) by its Python log name
rem   SESSION_CHECK.bat --date 2026-11-13 --all       another day (filed copies)
rem Reports: ..\analysis\session_checks\dayN\<session>_check.txt, ..\analysis\session_checks\test_runs\ and
rem          ..\analysis\session_checks\session_checks.csv (close it in Excel before running, or it is not updated)
rem Exit code: 0 all PASS, 1 something to check, 2 nothing to check, 3 a file is still being written,
rem            4 all PASS but session_checks.csv could not be updated (open in Excel? close it and run again)
setlocal
set "PY=%~dp0.venv\Scripts\python.exe"
if not exist "%PY%" (
    echo Python environment not found: %PY%
    echo Run tools\run_sync.bat once to create it, then try again.
    pause
    exit /b 1
)
"%PY%" -B "%~dp0session_check.py" %*
set "RC=%errorlevel%"
echo.
pause
exit /b %RC%
