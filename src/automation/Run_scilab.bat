@echo off
:inicio
echo Rodando script Scilab...
"C:\Users\administrator\AppData\Local\scilab-2025.0.0\bin\WScilex.exe" -f "C:\Users\administrator\Scilab_Model_Run_rev4.sci"
echo Scilab fechou ou crashou. Reiniciando...
timeout /t 5 > nul
goto inicio
