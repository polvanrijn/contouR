@echo off & setlocal
REM Set to home directory
cd %userprofile%

REM Unzip the file
powershell.exe -NoP -NonI -Command "Expand-Archive '.\f0_ascii.zip'"

REM move the files
move interpolation.exe f0_ascii
move Accent1Hz.fir f0_ascii
cd %userprofile%/f0_ascii

for /r %%i in (*.f0_ascii) do (
	interpolation %%i% 0 4 1e-006 auto 2
)