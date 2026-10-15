@echo off
cd /d "%~dp0"
python portable_launcher.py
if errorlevel 1 pause
