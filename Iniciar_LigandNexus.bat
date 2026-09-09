@echo off
setlocal
title LigandNexus
cd /d "%~dp0"

REM ============================================================================
REM  Lanzador de LigandNexus (Windows). La primera vez prepara un entorno local
REM  e instala las dependencias; las siguientes solo abre la app en el navegador.
REM  Requiere Python 3.10+ instalado (https://www.python.org/downloads/,
REM  marcando "Add Python to PATH").
REM
REM  OJO AL EDITAR: este archivo debe quedar en ASCII puro y con saltos CRLF.
REM  cmd.exe lo lee con la pagina de codigos OEM, asi que una tilde o un guion
REM  de caja rompe el parseo y las variables salen vacias.
REM ============================================================================

set "VENV=.venv"
set "PYEXE=%VENV%\Scripts\python.exe"
set "PUERTO=8501"

REM --- Hay ya una instancia abierta? ------------------------------------------
REM  IMPORTANTE: cerrar la PESTANA del navegador NO detiene el servidor; sigue
REM  vivo ocupando el puerto. Y Streamlit, si encuentra el puerto ocupado, NO
REM  da error: se queda esperando en silencio a que se libere. Sin esta com-
REM  probacion, volver a ejecutar el lanzador dejaba la ventana aparentemente
REM  colgada y apilaba un par de procesos zombis en cada intento.
REM  La deteccion usa una conexion TCP real y no lee netstat, porque el texto
REM  de netstat cambia con el idioma de Windows y sus lineas TIME_WAIT sobre-
REM  viven al cierre del servidor (darian falsos positivos).
powershell -NoProfile -Command "try{$c=New-Object Net.Sockets.TcpClient;$c.Connect('127.0.0.1',%PUERTO%);$c.Close();exit 0}catch{exit 1}"
if not errorlevel 1 goto YA_ABIERTA
goto ARRANCAR


:YA_ABIERTA
echo.
echo  ========================================================================
echo    LigandNexus YA ESTA EN MARCHA en el puerto %PUERTO%.
echo  ========================================================================
echo.
echo    Cerrar la pestana del navegador no detiene la aplicacion: el servidor
echo    sigue corriendo en su propia ventana de consola.
echo.
echo      [ENTER]  Abrir la app que ya esta corriendo.
echo      [R]      Reiniciarla ^(util si acaba de actualizar el codigo^).
echo      [S]      Salir sin hacer nada.
echo.
set /p "ACC=Opcion: "
if /i "%ACC%"=="R" goto REINICIAR
if /i "%ACC%"=="S" exit /b 0
start "" "http://localhost:%PUERTO%"
exit /b 0


:REINICIAR
echo.
echo  [LigandNexus] Deteniendo la instancia anterior...
REM  Se detiene EXACTAMENTE el proceso que ocupa el puerto (y el que lo lanzo),
REM  no cualquier python suelto: asi no se toca otra app del usuario.
powershell -NoProfile -Command "Get-NetTCPConnection -LocalPort %PUERTO% -State Listen -ErrorAction SilentlyContinue | ForEach-Object { $q = $_.OwningProcess; $par = (Get-CimInstance Win32_Process | Where-Object { $_.ProcessId -eq $q }).ParentProcessId; Stop-Process -Id $q -Force -ErrorAction SilentlyContinue; if ($par) { Stop-Process -Id $par -Force -ErrorAction SilentlyContinue } }"
timeout /t 3 /nobreak >nul
powershell -NoProfile -Command "try{$c=New-Object Net.Sockets.TcpClient;$c.Connect('127.0.0.1',%PUERTO%);$c.Close();exit 0}catch{exit 1}"
if not errorlevel 1 goto SIGUE_OCUPADO
echo  [LigandNexus] Instancia anterior detenida.
goto ARRANCAR


:SIGUE_OCUPADO
echo.
echo   ERROR: el puerto %PUERTO% sigue ocupado. Cierre a mano la ventana de
echo          consola de LigandNexus que haya abierta y vuelva a intentarlo.
echo.
pause
exit /b 1


:ARRANCAR
REM  OJO: dentro de un bloque entre parentesis, cmd CUENTA los parentesis que
REM  aparezcan en el texto de un echo. Un "(...)" seguido de un punto abortaba
REM  el script entero con "No se esperaba . en este momento", en TODAS las
REM  ejecuciones (aunque el .venv ya existiera, porque el bloque se analiza
REM  antes de evaluar la condicion). Por eso van escapados como ^( y ^).
if not exist "%PYEXE%" (
    echo [LigandNexus] Preparando el entorno por primera vez...
    where py >nul 2>nul && ( py -3 -m venv "%VENV%" ) || ( python -m venv "%VENV%" )
    if not exist "%PYEXE%" (
        echo.
        echo  ERROR: no se encontro Python. Instalelo desde
        echo         https://www.python.org/downloads/  ^(marque "Add Python to PATH"^).
        echo.
        pause
        exit /b 1
    )
    echo [LigandNexus] Instalando dependencias ^(puede tardar unos minutos^)...
    "%PYEXE%" -m pip install --upgrade pip --trusted-host pypi.org --trusted-host files.pythonhosted.org >nul
    "%PYEXE%" -m pip install -r requirements.txt --trusted-host pypi.org --trusted-host files.pythonhosted.org
    if errorlevel 1 (
        echo  ERROR instalando dependencias. Revise su conexion a internet.
        pause
        exit /b 1
    )
)

echo.
echo  [LigandNexus] Iniciando la aplicacion...
echo  [LigandNexus] El navegador se abrira solo cuando el servidor este listo.
echo.
echo  ------------------------------------------------------------------------
echo    PARA DETENER LA APP: CIERRE ESTA VENTANA ^(o pulse Ctrl+C aqui^).
echo    Cerrar la pestana del navegador NO la detiene.
echo  ------------------------------------------------------------------------
echo.

REM Abrir el navegador SOLO cuando el puerto responda (evita la pagina en
REM blanco por abrir antes de tiempo). Corre en paralelo, oculto.
start "" /b powershell -NoProfile -WindowStyle Hidden -Command ^
  "$ok=$false; for($i=0;$i -lt 120 -and -not $ok;$i++){ try{ $c=New-Object Net.Sockets.TcpClient; $c.Connect('127.0.0.1',%PUERTO%); $ok=$c.Connected; $c.Close() }catch{ Start-Sleep -Milliseconds 500 } }; if($ok){ Start-Process 'http://localhost:%PUERTO%' }"

REM Servidor en primer plano: esta ventana es su consola. Cierrela para detener la app.
"%PYEXE%" -m streamlit run app.py --server.headless true --server.address 127.0.0.1 --server.port %PUERTO%
endlocal
