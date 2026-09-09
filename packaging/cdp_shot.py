"""
Captura de pantalla de una app web (Streamlit) vía Chrome DevTools Protocol.
Espera tiempo REAL para que el websocket de Streamlit termine de renderizar.

Uso:  python cdp_shot.py <url> <salida.png> [espera_seg] [ancho] [alto]
"""
import sys, json, time, base64, subprocess, tempfile, urllib.request, os
from websocket import create_connection

URL = sys.argv[1]
OUT = sys.argv[2]
ESPERA = float(sys.argv[3]) if len(sys.argv) > 3 else 22
ANCHO = int(sys.argv[4]) if len(sys.argv) > 4 else 1440
ALTO = int(sys.argv[5]) if len(sys.argv) > 5 else 2200

CHROME = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
PORT = 9333
prof = tempfile.mkdtemp()
proc = subprocess.Popen([
    CHROME, "--headless=new", "--disable-gpu", "--no-sandbox", "--hide-scrollbars",
    f"--remote-debugging-port={PORT}", "--remote-allow-origins=*",
    f"--user-data-dir={prof}", f"--window-size={ANCHO},{ALTO}",
], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

try:
    # Esperar a que el endpoint de depuración responda.
    ws_url = None
    for _ in range(40):
        try:
            data = json.loads(urllib.request.urlopen(f"http://127.0.0.1:{PORT}/json").read())
            pages = [t for t in data if t.get("type") == "page"]
            if pages:
                ws_url = pages[0]["webSocketDebuggerUrl"]
                break
        except Exception:
            pass
        time.sleep(0.5)
    if not ws_url:
        raise RuntimeError("No se obtuvo el webSocketDebuggerUrl de Chrome")

    ws = create_connection(ws_url, max_size=None)
    _id = [0]
    def cmd(method, params=None):
        _id[0] += 1
        ws.send(json.dumps({"id": _id[0], "method": method, "params": params or {}}))
        while True:
            msg = json.loads(ws.recv())
            if msg.get("id") == _id[0]:
                return msg

    cmd("Page.enable")
    cmd("Runtime.enable")
    cmd("Emulation.setDeviceMetricsOverride",
        {"width": ANCHO, "height": ALTO, "deviceScaleFactor": 1, "mobile": False})
    cmd("Page.navigate", {"url": URL})
    time.sleep(ESPERA)   # tiempo REAL para el render del websocket de Streamlit
    res = cmd("Page.captureScreenshot", {"format": "png", "captureBeyondViewport": True})
    png = base64.b64decode(res["result"]["data"])
    with open(OUT, "wb") as f:
        f.write(png)
    print("OK:", OUT, len(png), "bytes")
    ws.close()
finally:
    proc.terminate()
