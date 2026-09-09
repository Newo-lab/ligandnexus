"""
Punto de arranque para el ejecutable congelado (PyInstaller) de LigandNexus.

Usa streamlit.web.bootstrap.run directamente (no la CLI), para evitar la
re-ejecución del proceso que rompe los .exe de Streamlit. Abre el navegador solo.
"""
import os
import sys
import threading
import time
import webbrowser


def _base_dir():
    # En el .exe congelado los datos van a sys._MEIPASS.
    if getattr(sys, "frozen", False):
        return sys._MEIPASS
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _abrir_navegador(port):
    time.sleep(4)
    webbrowser.open(f"http://localhost:{port}")


def main():
    base = _base_dir()
    app = os.path.join(base, "app.py")
    port = 8501

    from streamlit import config
    config.set_option("server.headless", True)
    config.set_option("server.port", port)
    config.set_option("browser.gatherUsageStats", False)
    config.set_option("global.developmentMode", False)

    threading.Thread(target=_abrir_navegador, args=(port,), daemon=True).start()

    from streamlit.web import bootstrap
    bootstrap.run(app, False, [], {})


if __name__ == "__main__":
    main()
