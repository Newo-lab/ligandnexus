"""
Captura con interacción: navega, escribe en un input y pulsa un botón por texto,
espera y captura. Para mostrar una vista con resultados.

Uso: python cdp_interact.py <url> <salida.png> <texto_a_escribir> <texto_boton>
"""
import sys, json, time, base64, subprocess, tempfile, urllib.request
from websocket import create_connection

URL, OUT, TEXTO, BOTON = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
ANCHO, ALTO = 1400, 1300
CHROME = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
PORT = 9334
prof = tempfile.mkdtemp()
proc = subprocess.Popen([
    CHROME, "--headless=new", "--disable-gpu", "--no-sandbox", "--hide-scrollbars",
    f"--remote-debugging-port={PORT}", "--remote-allow-origins=*",
    f"--user-data-dir={prof}", f"--window-size={ANCHO},{ALTO}",
], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

try:
    ws_url = None
    for _ in range(40):
        try:
            data = json.loads(urllib.request.urlopen(f"http://127.0.0.1:{PORT}/json").read())
            pages = [t for t in data if t.get("type") == "page"]
            if pages:
                ws_url = pages[0]["webSocketDebuggerUrl"]; break
        except Exception:
            pass
        time.sleep(0.5)
    ws = create_connection(ws_url, max_size=None)
    _id = [0]
    def cmd(method, params=None):
        _id[0] += 1
        ws.send(json.dumps({"id": _id[0], "method": method, "params": params or {}}))
        while True:
            m = json.loads(ws.recv())
            if m.get("id") == _id[0]:
                return m

    cmd("Page.enable"); cmd("Runtime.enable")
    cmd("Emulation.setDeviceMetricsOverride",
        {"width": ANCHO, "height": ALTO, "deviceScaleFactor": 1, "mobile": False})
    cmd("Page.navigate", {"url": URL})
    time.sleep(16)

    # 1) Click en el primer text input, escribir.
    js_box = """(() => {const el=document.querySelector('input[type="text"], input[type="number"]');
        if(!el) return null; el.focus(); const r=el.getBoundingClientRect();
        return JSON.stringify({x:r.x+r.width/2, y:r.y+r.height/2});})()"""
    res = cmd("Runtime.evaluate", {"expression": js_box, "returnByValue": True})
    box = res["result"]["result"].get("value")
    if box:
        b = json.loads(box)
        cmd("Input.dispatchMouseEvent", {"type": "mousePressed", "x": b["x"], "y": b["y"], "button": "left", "clickCount": 1})
        cmd("Input.dispatchMouseEvent", {"type": "mouseReleased", "x": b["x"], "y": b["y"], "button": "left", "clickCount": 1})
        cmd("Input.insertText", {"text": TEXTO})
        time.sleep(0.5)

    # 2) Click en el botón cuyo texto contiene BOTON.
    js_btn = """(() => {const bs=[...document.querySelectorAll('button')];
        const b=bs.find(x=>x.innerText && x.innerText.includes('%s'));
        if(!b) return null; const r=b.getBoundingClientRect();
        return JSON.stringify({x:r.x+r.width/2, y:r.y+r.height/2});})()""" % BOTON
    res2 = cmd("Runtime.evaluate", {"expression": js_btn, "returnByValue": True})
    box2 = res2["result"]["result"].get("value")
    if box2:
        b2 = json.loads(box2)
        cmd("Input.dispatchMouseEvent", {"type": "mousePressed", "x": b2["x"], "y": b2["y"], "button": "left", "clickCount": 1})
        cmd("Input.dispatchMouseEvent", {"type": "mouseReleased", "x": b2["x"], "y": b2["y"], "button": "left", "clickCount": 1})

    time.sleep(9)   # esperar la consulta + render
    shot = cmd("Page.captureScreenshot", {"format": "png", "captureBeyondViewport": True})
    open(OUT, "wb").write(base64.b64decode(shot["result"]["data"]))
    print("OK:", OUT)
    ws.close()
finally:
    proc.terminate()
