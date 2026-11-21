"""
TFI Teoría de Control - Control automático de presión de descarga en una estación de bombeo
(versión completa, bloque 4)

Lazo (notación de la cátedra, ver "00 - Diagrama de bloques del lazo"):
    θi -> Σ -> e -> Gc (PI en el PLC) -> u -> Ga (variador + bomba) -> Q_b -> Σ(+,+) <- θd
       -> Gp (tramo de oleoducto) -> θo (presión de descarga)
    θo -> H (transmisor de presión) -> f -> vuelve a Σ con signo menos
Supervisión (fuera del lazo, con instrumentos propios): detección de fugas por balance de caudal
y disparo por alta presión.

Señales:
    θi  referencia: 60 bar (13,6 mA en la escala 4-20 mA del transmisor; el PLC la escala a bar)
    f   señal de medición (lectura del transmisor, en bar después del escalado)
    e   error = θi - f
    u   salida del controlador: velocidad de la bomba [%]
    θd  perturbación de caudal: cambio de demanda o fuga [m3/h]
    θo  presión real de descarga [bar]  (valor nominal = 60 bar)

Modelo:
    Ga: Q_b sigue a K_B * u con constante de tiempo TAU_V
    Gp: C_L * dθo/dt = Q_b - Q_out - Q_f,  Q_out = k_out * (θo - P_SIG)
    H : primer orden con constante TAU_T y ruido
    Gc: u = U0 + Kp*e + Ki*∫e dt - Kv*df/dt (derivativa sobre la medición filtrada; Kv = 0 la inhibe),
        con saturación 0-100 % y anti-windup

Uso:
    python simulacion_presion.py          -> imprime resultados y guarda simulacion_presion.png
    python simulacion_presion.py --acciones -> solo P, solo I, solo D, PI y PID; margen de fase y demora crítica
    python simulacion_presion.py --vivo   -> tablero de control en tiempo real: iniciar/detener, perturbación
                                             con tipo, amplitud y duración, referencia, Kp, Ki y Kv,
                                             todas las señales apiladas y QoS (transitorio o falla)
"""
import sys
import numpy as np
import matplotlib
if "--vivo" not in sys.argv:
    matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------- Proceso: tramo de 80 km entre estaciones (oleoducto de 32") ----------------
L, D, e = 80e3, 0.8128, 0.0127
K_CRUDO, E_ACERO = 15000.0, 2.07e6                       # [bar]
C_L = np.pi / 4 * D**2 * L * (1 / K_CRUDO + D / (e * E_ACERO))   # capacitancia [m3/bar]
P_REF, P_SIG = 60.0, 5.0                                  # θi (= valor nominal) y presión en la estación siguiente [bar]
Q_NOM = 1900.0                                            # caudal nominal [m3/h]
K_OUT = Q_NOM / (P_REF - P_SIG)                           # [m3/h por bar]
P_MAX = 70.0                                              # presión máxima admisible del caño [bar]
BANDA = 0.01 * P_REF                                      # desvío admisible: ±1 % (±0,6 bar)

# ---------------- Actuador (Ga) y elemento de medición (H) ----------------
K_B = 25.0                                                # caudal de la bomba por % de velocidad [m3/h/%]
TAU_V = 10.0 / 3600                                       # aceleración del variador: 10 s [h]
TAU_T = 1.0 / 3600                                        # transmisor con amortiguamiento: 1 s [h]
SIG_P = 0.05                                              # ruido del transmisor [bar]
SIG_Q = 10.0                                              # ruido de los caudalímetros [m3/h]
TS = 1.0 / 3600                                           # ciclo (scan) del PLC: 1 s [h]
U0 = Q_NOM / K_B                                          # velocidad de régimen [%]

# ---------------- Controlador (Gc): PID con la derivativa inhibida = PI ----------------
KP = 10.0                                                 # [% por bar]
KI = 400.0                                                # [% por bar·h]  (tiempo integral Kp/Ki = 1,5 min)
KV = 0.0                                                  # derivativa [% por bar/min]: 0 = inhibida
TAU_D = 5.0 / 3600                                        # filtro de la derivativa: 5 s [h]


class Lazo:
    """Estado del lazo cerrado. Cada llamada a paso() es un ciclo (scan) del PLC."""

    def __init__(self, kp=KP, ki=KI, supervision=True, seed=1, kv=KV):
        self.rng = np.random.default_rng(seed)
        self.kp, self.ki, self.kv, self.supervision = kp, ki, kv, supervision
        self.vn = P_REF                    # valor nominal (θi en bar, después del escalado del PLC)
        self.theta_o, self.Q_b, self.P_med, self.integ = P_REF, Q_NOM, P_REF, 0.0
        self.f_filt = P_REF                # medición filtrada para la acción derivativa
        self.k_out = K_OUT                 # baja si cae la demanda aguas abajo
        self.fuga = 0.0                    # caudal fugado [m3/h]
        self.ataque = 0.0                  # bar que se le restan a la lectura del transmisor
        self.parada, self.motivo = None, None
        self.buf, self.p_sup = [], []
        self.k = 0

    def paso(self):
        t = self.k * TS
        Q_f = self.fuga if self.parada is None else 0.0   # la parada aísla el tramo

        # H: elemento de medición (con posible falsificación de la lectura)
        self.P_med += TS / TAU_T * (self.theta_o - self.P_med)
        f = self.P_med + self.rng.normal(0, SIG_P) - self.ataque

        # Σ y Gc: error y PID (derivativa sobre la medición filtrada) con saturación y anti-windup
        err = self.vn - f
        f_prev = self.f_filt
        self.f_filt += TS / TAU_D * (f - self.f_filt)
        uP, uI = self.kp * err, self.ki * self.integ
        uD = -self.kv * (self.f_filt - f_prev) / (TS * 60)   # por minuto
        u = U0 + uP + uI + uD
        if 0.0 < u < 100.0:
            self.integ += err * TS
        u = min(max(u, 0.0), 100.0)
        if self.parada is not None:
            u = 0.0

        # Ga y Gp: actuador y proceso
        self.Q_b += TS / TAU_V * (K_B * u - self.Q_b)
        Q_out = self.k_out * (self.theta_o - P_SIG)
        self.theta_o += TS * (self.Q_b - Q_out - Q_f) / C_L

        # Supervisión (instrumentos independientes del lazo de control):
        # balance de caudal corregido por empaque de línea, promediado 30 s,
        # y disparo por alta presión real
        self.buf.append((self.Q_b + self.rng.normal(0, SIG_Q)) - (Q_out + self.rng.normal(0, SIG_Q)))
        self.p_sup.append(self.theta_o + self.rng.normal(0, SIG_P))
        self.buf, self.p_sup = self.buf[-30:], self.p_sup[-30:]
        balance = np.mean(self.buf) - C_L * (self.p_sup[-1] - self.p_sup[0]) / ((len(self.p_sup) - 1) * TS + 1e-9)
        if self.supervision and self.parada is None:
            if balance > 150.0 and self.k > 60:
                self.parada, self.motivo = t, "fuga"
            elif self.theta_o > P_MAX:
                self.parada, self.motivo = t, "alta presión"

        self.k += 1
        theta_d = (K_OUT - self.k_out) * (self.theta_o - P_SIG) - Q_f   # caudal que agrega la perturbación
        return dict(t=t, theta_o=self.theta_o, f=f, e=err, u=u, Q_b=self.Q_b, balance=balance,
                    Q_out=Q_out, uP=uP, uI=uI, uD=uD, theta_d=theta_d, n=-self.ataque, theta_i=self.vn)


def simular(horas=1.5, kp=KP, ki=KI, cambio_demanda=0.25, fuga=0.0, t_fuga=0.75,
            ataque=0.0, t_ataque=0.75, supervision=True, seed=1):
    """cambio_demanda: hora en que la demanda aguas abajo baja 15 % (None = no hay).
    ataque: bar que se le restan a la lectura del transmisor desde t_ataque."""
    lazo = Lazo(kp, ki, supervision, seed)
    n = int(horas / TS)
    t = np.arange(n) * TS
    reg = {k: np.zeros(n) for k in ("theta_o", "f", "e", "u", "Q_b", "balance")}
    for k in range(n):
        if cambio_demanda is not None and t[k] >= cambio_demanda:
            lazo.k_out = 0.85 * K_OUT                     # θd: baja la demanda
        lazo.fuga = fuga if t[k] >= t_fuga else 0.0       # θd: fuga
        lazo.ataque = ataque if t[k] >= t_ataque else 0.0 # perturbación en la medición
        s = lazo.paso()
        for key in reg:
            reg[key][k] = s[key]
    return t, reg, lazo.parada, lazo.motivo


def routh_kp_critico():
    """Ganancia proporcional crítica del lazo (P solo) por Routh-Hurwitz."""
    tau_p = C_L / K_OUT                                   # constante de tiempo del tramo [h]
    t1, t2, t3 = tau_p, TAU_V, TAU_T
    a3, a2, a1 = t1 * t2 * t3, t1 * t2 + t1 * t3 + t2 * t3, t1 + t2 + t3
    k_crit = a2 * a1 / a3 - 1                             # ganancia de lazo total
    return tau_p, k_crit, k_crit * K_OUT / K_B           # kp crítico en % por bar


def routh_pi(kp=KP, ki=KI):
    """Primera columna de Routh del lazo con PI: s(τp s+1)(τv s+1)(τt s+1) + (K_B/k_out)(Kp s + Ki) = 0."""
    p = np.polymul(np.polymul(np.polymul([C_L / K_OUT, 1], [TAU_V, 1]), [TAU_T, 1]), [1, 0])
    p = np.polyadd(p, [K_B / K_OUT * kp, K_B / K_OUT * ki])
    filas = [list(p[0::2]), list(p[1::2]) + [0.0]]
    ancho = len(filas[0])
    filas = [r + [0.0] * (ancho - len(r)) for r in filas]
    while len(filas) < len(p):
        a, b = filas[-2], filas[-1]
        filas.append([(b[0] * a[i + 1] - a[0] * b[i + 1]) / b[0] if i + 1 < ancho else 0.0
                      for i in range(ancho)])
    return [r[0] for r in filas], np.roots(p)


def vivo(velocidad=30.0, mostrar=True):
    """Tablero de control en tiempo real (pedido del profesor, 24/10).
    - Iniciar / Detener / Reiniciar: arranca en condiciones normales y se puede congelar en cualquier momento.
    - Perturbación con tipo, amplitud y duración (0 = sostenida).
    - Valor de referencia, Kp, Ki y Kv ajustables en vivo.
    - Todas las señales del lazo, una debajo de la otra.
    - QoS: banda ±1 % y 120 s tolerados fuera de banda; más tiempo, o superar la presión máxima, es falla.
    velocidad: segundos de proceso por segundo de reloj (el tramo responde en minutos)."""
    from matplotlib.animation import FuncAnimation
    from matplotlib.widgets import Slider, Button, RadioButtons, TextBox

    VENTANA = 600                                         # s de proceso visibles
    st = {"corriendo": False, "vel": velocidad, "pert": None, "fuera": None, "falla": False,
          "fallas": 0, "log": []}

    def nuevo_lazo():
        lz = Lazo(kp=s_kp.val, ki=s_ki.val * 60, kv=s_kv.val)   # Ki del tablero en % por bar·min
        lz.vn = s_vn.val
        return lz

    fig = plt.figure(figsize=(14, 10))
    fig.canvas.manager.set_window_title("Tablero del lazo de presión") if hasattr(fig.canvas, "manager") else None
    gs = fig.add_gridspec(9, 2, width_ratios=[1, 2.6], left=0.04, right=0.98, top=0.93, bottom=0.05, hspace=0.35, wspace=0.12)
    paneles = [("θi · referencia [bar]", ["theta_i"]), ("e = θi − f · error [bar]", ["e"]),
               ("u · salida del controlador [%]", ["u"]), ("aportes P · I · D [%]", ["uP", "uI", "uD"]),
               ("Q_b bomba y Q_out salida [m³/h]", ["Q_b", "Q_out"]), ("θd · perturbación de caudal [m³/h]", ["theta_d"]),
               ("θo · presión de descarga [bar]", ["theta_o"]), ("n · perturbación en la medición [bar]", ["n"]),
               ("f · lectura del transmisor [bar]", ["f"])]
    colores = {"uI": "tab:blue", "uD": "tab:orange", "Q_out": "gray", "theta_d": "tab:orange", "n": "tab:orange",
               "theta_i": "tab:blue"}
    axs, lineas = [], {}
    for i, (tit, keys) in enumerate(paneles):
        ax = fig.add_subplot(gs[i, 1], sharex=axs[0] if axs else None)
        ax.set_title(tit, fontsize=8, loc="left", pad=2)
        ax.tick_params(labelsize=7)
        for k in keys:
            lineas[k], = ax.plot([], [], lw=1.0, color=colores.get(k, "black"),
                                 ls="--" if k == "Q_out" else "-", label=k)
        if len(keys) > 1:
            ax.legend(fontsize=6, loc="upper left", ncol=3)
        if i < len(paneles) - 1:
            ax.tick_params(labelbottom=False)
        axs.append(ax)
    ax_o = axs[6]
    banda_sp = ax_o.axhspan(59.4, 60.6, color="tab:blue", alpha=0.1)
    linea_max = ax_o.axhline(P_MAX, color="tab:red", ls=":", lw=0.8)
    axs[-1].set_xlabel("tiempo [min]", fontsize=8)

    # --- controles (columna izquierda) ---
    def eje(y, h=0.025):
        return fig.add_axes([0.10, y, 0.18, h])
    s_vn = Slider(eje(0.88), "VN [bar]", 40, 68, valinit=P_REF, valstep=0.5)
    s_kp = Slider(eje(0.84), "Kp [%/bar]", 0, 200, valinit=KP, valstep=1)
    s_ki = Slider(eje(0.80), "Ki [%/bar·min]", 0, 60, valinit=round(KI / 60, 1), valstep=0.1)
    s_kv = Slider(eje(0.76), "Kv [%·min/bar]", 0, 20, valinit=KV, valstep=0.1)
    for s in (s_vn, s_kp, s_ki, s_kv):
        s.label.set_fontsize(8); s.valtext.set_fontsize(8)
    r_tipo = RadioButtons(fig.add_axes([0.04, 0.60, 0.24, 0.12]),
                          ("Demanda aguas abajo [%]", "Fuga [m³/h]", "Transmisor (n) [bar]"))
    for lbl in r_tipo.labels:
        lbl.set_fontsize(8)
    t_amp = TextBox(fig.add_axes([0.12, 0.555, 0.06, 0.03]), "Amplitud ", initial="15")
    t_dur = TextBox(fig.add_axes([0.24, 0.555, 0.04, 0.03]), "Duración [s] ", initial="0")
    t_tol = TextBox(fig.add_axes([0.24, 0.515, 0.04, 0.03]), "Tolerado fuera de banda [s] ", initial="120")
    b_run = Button(fig.add_axes([0.04, 0.45, 0.075, 0.04]), "Iniciar")
    b_reset = Button(fig.add_axes([0.12, 0.45, 0.075, 0.04]), "Reiniciar")
    b_pert = Button(fig.add_axes([0.04, 0.40, 0.155, 0.04]), "Generar perturbación", color="moccasin")
    b_cut = Button(fig.add_axes([0.20, 0.40, 0.08, 0.04]), "Cortar")
    r_vel = RadioButtons(fig.add_axes([0.20, 0.445, 0.08, 0.05]), ("×10", "×30", "×60"), active=1)
    for lbl in r_vel.labels:
        lbl.set_fontsize(8)
    txt_log = fig.text(0.04, 0.37, "", fontsize=7.5, va="top", family="monospace")
    txt_estado = fig.text(0.04, 0.965, "", fontsize=12, weight="bold")
    fig.text(0.04, 0.94, "Estación de bombeo · PID del PLC sobre variador y bomba · tramo de 80 km", fontsize=9, color="gray")

    lazo = nuevo_lazo()
    hist = {k: [] for k in ["t", "theta_i", "e", "u", "uP", "uI", "uD", "Q_b", "Q_out", "theta_d", "theta_o", "n", "f"]}

    def registrar(msg):
        t = lazo.k * TS * 3600
        st["log"].insert(0, f"{int(t // 60):02d}:{int(t % 60):02d}  {msg}")
        st["log"] = st["log"][:12]
        txt_log.set_text("\n".join(st["log"]))
        print(st["log"][0])

    def tipo_actual():
        return {"Demanda aguas abajo [%]": "demanda", "Fuga [m³/h]": "fuga", "Transmisor (n) [bar]": "medicion"}[r_tipo.value_selected]

    def aplicar_pert(p):
        lazo.k_out, lazo.fuga, lazo.ataque = K_OUT, 0.0, 0.0
        if p:
            if p["tipo"] == "demanda":
                lazo.k_out = K_OUT * (1 - p["amp"] / 100)
            elif p["tipo"] == "fuga":
                lazo.fuga = p["amp"]
            else:
                lazo.ataque = -p["amp"]

    def generar(_=None):
        try:
            amp, dur = float(t_amp.text), max(0.0, float(t_dur.text))
        except ValueError:
            registrar("Amplitud o duración inválida"); return
        st["pert"] = {"tipo": tipo_actual(), "amp": amp, "dur": dur, "t0": lazo.k * TS * 3600}
        aplicar_pert(st["pert"])
        registrar(f"Perturbación {st['pert']['tipo']} {amp:g}, {'sostenida' if dur == 0 else f'{dur:g} s'}")

    def cortar(_=None):
        if st["pert"]:
            st["pert"] = None; aplicar_pert(None); registrar("Perturbación cortada")

    def run(_=None):
        st["corriendo"] = not st["corriendo"]
        b_run.label.set_text("Detener" if st["corriendo"] else "Continuar")
        registrar("Simulación en marcha" if st["corriendo"] else "Simulación detenida")

    def reiniciar(_=None):
        nonlocal lazo
        lazo = nuevo_lazo(); st.update(corriendo=False, pert=None, fuera=None, falla=False, fallas=0, log=[])
        b_run.label.set_text("Iniciar")
        for v in hist.values():
            v.clear()
        for _ in range(90):
            avanzar()
        registrar("Sistema en condiciones normales")

    def actualizar_ganancias(_=None):
        lazo.kp, lazo.ki, lazo.kv, lazo.vn = s_kp.val, s_ki.val * 60, s_kv.val, s_vn.val

    for s in (s_vn, s_kp, s_ki, s_kv):
        s.on_changed(actualizar_ganancias)
    b_run.on_clicked(run); b_reset.on_clicked(reiniciar); b_pert.on_clicked(generar); b_cut.on_clicked(cortar)
    r_vel.on_clicked(lambda v: st.update(vel={"×10": 10, "×30": 30, "×60": 60}[v]))
    fig.canvas.mpl_connect("key_press_event", lambda ev: {" ": run, "p": generar, "c": cortar, "r": reiniciar}.get(ev.key, lambda: None)())

    def qos(s):
        try:
            tol = float(t_tol.text)
        except ValueError:
            tol = 120.0
        banda = lazo.vn * 0.01
        t = s["t"] * 3600 + 1
        if abs(s["theta_o"] - lazo.vn) > banda:
            if st["fuera"] is None:
                st["fuera"] = t; registrar(f"θo sale de la banda ±1 % ({s['theta_o']:.2f} bar)")
            if not st["falla"] and (t - st["fuera"] >= tol or s["theta_o"] > P_MAX):
                st["falla"] = True; st["fallas"] += 1
                registrar("FALLA: supera la presión máxima" if s["theta_o"] > P_MAX
                          else f"FALLA: más de {tol:g} s fuera de banda")
        elif st["fuera"] is not None:
            dur = t - st["fuera"]
            registrar(f"θo vuelve a la banda en {dur:.0f} s" + (" (falla registrada)" if st["falla"]
                                                                   else ": transitoria, QoS cumplida"))
            st["fuera"], st["falla"] = None, False

    def avanzar():
        p = st["pert"]
        if p and p["dur"] > 0 and lazo.k * TS * 3600 >= p["t0"] + p["dur"]:
            st["pert"] = None; aplicar_pert(None); registrar(f"Fin de la perturbación ({p['dur']:g} s)")
        parada_antes = lazo.parada
        s = lazo.paso()
        if lazo.parada is not None and parada_antes is None:
            registrar(f"PARADA de la bomba por {lazo.motivo}")
        qos(s)
        hist["t"].append(s["t"] * 60)
        for k in hist:
            if k != "t":
                hist[k].append(s[k])
        if len(hist["t"]) > VENTANA:
            for v in hist.values():
                del v[0]

    def cuadro(_):
        if st["corriendo"]:
            for _ in range(max(1, int(st["vel"] * 0.1))):     # cada cuadro son 100 ms de reloj
                avanzar()
        if not hist["t"]:
            return []
        for k, ln in lineas.items():
            ln.set_data(hist["t"], hist[k])
        t1 = hist["t"][-1]
        axs[0].set_xlim(max(0.0, t1 - VENTANA / 60), max(VENTANA / 60, t1))
        span_min = [2, 1, 10, 4, 200, 100, 3, 2, 3]          # rango mínimo de cada panel (evita ampliar el ruido)
        for ax, (tit, keys), m in zip(axs, paneles, span_min):
            vals = [v for k in keys for v in hist[k]]
            lo, hi = min(vals), max(vals)
            if keys == ["theta_o"]:
                lo, hi = min(lo, lazo.vn - 2 * lazo.vn * 0.01), max(hi, lazo.vn + 2 * lazo.vn * 0.01)
            c, sp_ = (lo + hi) / 2, max(hi - lo, m) * 1.15
            ax.set_ylim(c - sp_ / 2, c + sp_ / 2)
        b = lazo.vn * 0.01
        if hasattr(banda_sp, "set_height"):                 # matplotlib reciente: Rectangle
            banda_sp.set_y(lazo.vn - b); banda_sp.set_height(2 * b)
        else:                                               # versiones anteriores: Polygon
            banda_sp.set_xy([[0, lazo.vn - b], [0, lazo.vn + b], [1, lazo.vn + b], [1, lazo.vn - b], [0, lazo.vn - b]])
        estado = ("PARADA · " + lazo.motivo if lazo.parada else "FALLA · fuera de QoS" if st["falla"]
                  else "Perturbación transitoria" if st["fuera"] is not None
                  else "Perturbación absorbida" if st["pert"] else "Normal")
        t = hist["t"][-1]
        txt_estado.set_text(f"t = {int(t):02d}:{int(t * 60 % 60):02d}   {estado}   ·   θo = {hist['theta_o'][-1]:.2f} bar"
                            f"   ·   fallas registradas: {st['fallas']}" + ("" if st["corriendo"] else "   (detenida)"))
        txt_estado.set_color("firebrick" if (lazo.parada or st["falla"]) else "darkorange" if st["fuera"] else "black")
        return list(lineas.values())

    reiniciar()
    anim = FuncAnimation(fig, cuadro, interval=100, blit=False, cache_frame_data=False)
    if mostrar:
        plt.show()
    return anim, cuadro, dict(run=run, generar=generar, cortar=cortar, reiniciar=reiniciar,
                              estado=st, lazo=lambda: lazo, t_amp=t_amp, t_dur=t_dur, r_tipo=r_tipo, s_ki=s_ki, s_kp=s_kp)



# ---------------- Iteración 3 (devolución del 31/10): acciones por separado y demoras ----------------
ACCIONES = [("Sin control (lazo abierto)", 0.0, 0.0, 0.0), ("Solo P", 10.0, 0.0, 0.0), ("Solo I", 0.0, 400.0, 0.0),
            ("Solo D", 0.0, 0.0, 5.0), ("PI (elegido)", 10.0, 400.0, 0.0), ("PID (Kv = 5)", 10.0, 400.0, 5.0)]


def respuesta(kp, ki, kv, minutos=40, t_cambio_s=60, seed=1):
    """Baja de 15 % en la demanda en t_cambio_s; devuelve θo, u y e cada 1 s (sin supervisión)."""
    lazo = Lazo(kp, ki, False, seed, kv)
    th, u, e = [], [], []
    for k in range(int(minutos * 60)):
        if k == t_cambio_s:
            lazo.k_out = 0.85 * K_OUT
        s = lazo.paso()
        th.append(s["theta_o"]); u.append(s["u"]); e.append(s["e"])
    return np.array(th), np.array(u), np.array(e)


def comparar_acciones(minutos=40):
    """Error en estado estable, pico, tiempo fuera de banda y movimiento de la bomba para cada acción."""
    filas = []
    for nombre, kp, ki, kv in ACCIONES:
        th, u, e = respuesta(kp, ki, kv, minutos)
        d = th[60:] - P_REF
        fuera, corrida, maximo = np.abs(d) > BANDA, 0, 0
        for x in fuera:
            corrida = corrida + 1 if x else 0
            maximo = max(maximo, corrida)
        filas.append(dict(nombre=nombre, pico=d.max(), final=d[-120:].mean(), max_fuera_s=maximo,
                          desvio_u=u[-600:].std()))
    return filas


def margen_y_demora(kp=KP, ki=KI):
    """Margen de fase del lazo nominal y demora pura que lo lleva al límite de estabilidad."""
    w = np.logspace(-1, 4, 200000)                              # rad/h
    s = 1j * w
    lazo_abierto = (kp + ki / s) * K_B / (TAU_V * s + 1) / K_OUT / (C_L / K_OUT * s + 1) / (TAU_T * s + 1)
    i = np.argmin(np.abs(np.abs(lazo_abierto) - 1))
    mf = np.degrees(np.angle(lazo_abierto[i])) + 180
    return mf, np.radians(mf) / w[i] * 3600                     # grados, segundos


def con_demora(demora_s, kp=KP, ki=KI, minutos=30):
    """Lazo PI sin ruido con una demora pura en la medición (por ejemplo, cerrar el lazo a través de la nube)."""
    th, Qb, Pm, integ = P_REF, Q_NOM, P_REF, 0.0
    cola = [P_REF] * demora_s
    k_out, h = K_OUT, []
    for k in range(minutos * 60):
        if k == 60:
            k_out = 0.85 * K_OUT
        Pm += TS / TAU_T * (th - Pm)
        cola.append(Pm); f = cola.pop(0)
        e = P_REF - f
        u = U0 + kp * e + ki * integ
        if 0 < u < 100:
            integ += e * TS
        u = min(max(u, 0.0), 100.0)
        Qb += TS / TAU_V * (K_B * u - Qb)
        th += TS * (Qb - k_out * (th - P_SIG)) / C_L
        h.append(th)
    return np.array(h)


if __name__ == "__main__":
    if "--vivo" in sys.argv:
        vivo()
        sys.exit()
    if "--acciones" in sys.argv:
        print("Baja de 15 % en la demanda, cada acción por separado:")
        for r in comparar_acciones():
            print(f"  {r['nombre']:<27} pico {r['pico']:+5.2f} bar | error final {r['final']:+5.2f} bar | "
                  f"máx. seguido fuera de banda {r['max_fuera_s']:4d} s | desvío de u {r['desvio_u']:4.2f} %")
        for kp, ki in [(KP, KI), (KP, 0.0), (0.0, KI)]:
            mf, dc = margen_y_demora(kp, ki)
            print(f"Kp = {kp:g}, Ki = {ki:g}: margen de fase {mf:.0f}°, demora crítica {dc:.0f} s")
        for d in (0, 30, 50):
            h = con_demora(d)[60:] - P_REF
            print(f"PI con demora de {d:2d} s: pico {h.max():.2f} bar, oscilación al final ±{np.abs(h[-300:]).max():.2f} bar")
        sys.exit()

    tau_p, k_crit, kp_crit = routh_kp_critico()
    print(f"C_L = {C_L:.2f} m3/bar | constante de tiempo del tramo = {tau_p * 60:.1f} min")
    print(f"Routh (solo P): ganancia de lazo crítica = {k_crit:.0f}  ->  Kp crítico ~ {kp_crit:.0f} % por bar")
    col, raices = routh_pi()
    print(f"Routh (PI, Kp = {KP:.0f}, Ki = {KI:.0f}): primera columna "
          f"{'toda positiva -> estable' if all(c > 0 for c in col) else 'con cambio de signo -> inestable'}; "
          f"polos dominantes {raices[-1].real / 60:.2f} ± {abs(raices[-1].imag) / 60:.2f}j por minuto\n")

    for nombre, kw in [("Lazo abierto (SLA)", dict(kp=0.0, ki=0.0)), ("Solo P (Ki = 0)", dict(ki=0.0)), ("PI", {})]:
        t, r, _, _ = simular(**kw)
        final = r["theta_o"][-600:].mean()
        print(f"{nombre:<19}: demanda baja 15 % -> presión final {final:5.2f} bar "
              f"(error en estado estable {P_REF - final:+.2f} bar; banda admisible ±{BANDA:.1f})")

    t, r, _, _ = simular(fuga=200.0, supervision=False)
    i = int(0.74 / TS)
    print(f"\nFuga de 200 m3/h sin supervisión: presión antes {r['theta_o'][i]:.1f} bar, 20 min después "
          f"{r['theta_o'][i + 1200]:.1f} bar; la bomba pasa de {r['Q_b'][i]:.0f} a {r['Q_b'][i + 1200]:.0f} m3/h")
    t, r, parada, motivo = simular(fuga=200.0)
    print(f"Fuga de 200 m3/h con supervisión: parada por {motivo} a los {(parada - 0.75) * 3600:.0f} s")
    t, r, parada, motivo = simular()
    print(f"Cambio de demanda sin fuga, con supervisión: {'falsa parada por ' + motivo if parada else 'sin parada'}")
    for sup in (False, True):
        t, r, parada, motivo = simular(ataque=15.0, supervision=sup)
        print(f"Transmisor falseado -15 bar, {'con' if sup else 'sin'} disparo independiente: "
              f"presión real máxima {r['theta_o'].max():.1f} bar" +
              (f", parada por {motivo}" if parada else ""))

    esc = [("Solo P: cambio de demanda deja error permanente", dict(ki=0.0)),
           ("PI: el mismo cambio, error nulo", {}),
           ("PI con fuga de 200 m3/h: la presión no la muestra, el balance sí", dict(fuga=200.0)),
           ("Transmisor falseado -15 bar: la presión real sube hasta el disparo", dict(ataque=15.0))]
    fig, ax = plt.subplots(4, 1, figsize=(10, 11), sharex=True)
    for a, (tit, kw) in zip(ax, esc):
        t, r, parada, motivo = simular(**kw)
        a.plot(t * 60, r["theta_o"], lw=1.0, label="θo: presión real")
        a.plot(t * 60, r["f"], lw=0.5, alpha=0.6, label="f: lectura del transmisor")
        a.axhline(P_REF, color="gray", ls="--", lw=0.8)
        a.axhspan(P_REF - BANDA, P_REF + BANDA, color="green", alpha=0.08)
        a.axhline(P_MAX, color="tab:red", ls=":", lw=0.8)
        if parada is not None:
            a.axvline(parada * 60, color="black", lw=1)
        a.set_ylabel("bar")
        a.set_title(tit + (f" (parada por {motivo})" if parada else ""), fontsize=10)
        a.legend(loc="lower left", fontsize=8)
    ax[-1].set_xlabel("tiempo [min]  (gris: θi = 60 bar · verde: banda ±1 % · rojo: máximo 70 bar · negro: parada)")
    fig.tight_layout()
    fig.savefig("simulacion_presion.png", dpi=120)
    print("\nFigura guardada en simulacion_presion.png")
