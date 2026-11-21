# Control automático de presión de descarga en una estación de bombeo de oleoducto

Trabajo Final Integrador de **Teoría de Control**: UTN FRBA, comisión K4521, prof. Mgtr. Omar Oscar Civale, noviembre de 2026.
Alumnos: Pedro Mazzieri y Victoria Yedid.

## De qué se trata

Una estación de bombeo tiene que entregarle el crudo a la estación siguiente con una presión de descarga estable. El trabajo diseña, analiza y simula el lazo cerrado que mantiene esa presión en **60 bar ± 1 %**:

- un controlador **PI en el PLC** de la estación lee la presión con un transmisor 4-20 mA;
- según esa lectura, ajusta la velocidad de la bomba a través de un variador;
- el proceso es un tramo de oleoducto de 80 km.

Fuera del lazo hay una capa de **supervisión** con instrumentos propios: detección de fugas por balance de caudal y disparo por alta presión (PSHH-102, a 70 bar). Hace falta porque hay perturbaciones que el lazo no ve o no puede corregir:

- **Una fuga:** el PI la compensa acelerando la bomba, así que la presión no la delata.
- **Un transmisor falseado:** su error entra por la realimentación y se traslada a la presión real.

Este repositorio tiene el código de la simulación. El informe se entrega por separado.

```
                                                                         θd: perturbación de caudal
                                                                           │
                                                                           ▼
θi ──►(Σ)── e ──► Gc: PI en el PLC ── u ──► Ga: variador + bomba ── Qb ──►(Σ)──► Gp: tramo de 80 km ──┬──► θo ──► carga
       ▲ −                                                                                            │
       └── f ◄── H: transmisor 4-20 mA (n entra en la medición) ◄─────────────────────────────────────┘
```

## Contenido

| Archivo | Qué es |
|---|---|
| `simulacion/tablero_PIC-101.html` | Tablero de control que corre en el navegador sin instalar nada. Es la simulación que se muestra en la presentación. |
| `simulacion/simulacion_presion.py` | Programa en Python con el modelo del lazo. Calcula la estabilidad por Routh, corre los escenarios, compara las acciones P, I y D, mide el margen de fase y la demora que tolera el lazo, y tiene un modo en vivo con matplotlib. |
| `simulacion/TFI_control_presion.ipynb` | Cuaderno de Google Colab con el mismo modelo, paso a paso: modelo físico, sintonía, Routh, los 10 escenarios con todas las señales, cada acción por separado, latencia y un tablero con ipywidgets. |
| `requirements.txt` | Dependencias del programa y del cuaderno. |

El tablero, el programa y el cuaderno usan el mismo modelo y las mismas constantes, así que dan los mismos resultados.

## Cómo correrlo

### Tablero

Abrir `simulacion/tablero_PIC-101.html` en el navegador. Incluye:

- **Escenarios:** los 10 de la tabla de resultados, precargados.
- **Perillas en vivo:** baja de demanda y fuga (θd) y error del transmisor (n).
- **Referencia y constantes:** θi, Kp, Ki y Kv, más casillas para habilitar o deshabilitar las acciones P, I y D.
- **Modo:** manual o automático, con transferencia sin golpe.
- **Control de la ejecución:** por tramos, con pausa automática y la opción de guardar un estado y volver a él.
- **QoS:** banda de ±1 % con 120 s tolerados fuera. El registro de eventos avisa si queda un error permanente.

Atajos de teclado:

| Tecla | Acción |
|---|---|
| Espacio | Inicia o detiene |
| T | Avanza un tramo |
| P | Genera una perturbación |
| C | Corta la perturbación |
| R | Reinicia |

### Programa en Python

```bash
pip install -r requirements.txt
cd simulacion
python simulacion_presion.py              # Routh, resultados principales y gráfico simulacion_presion.png
python simulacion_presion.py --acciones   # solo P, solo I, solo D, PI y PID; margen de fase y demora crítica
python simulacion_presion.py --vivo       # tablero en tiempo real con matplotlib
```

### Cuaderno

En Google Colab: *Archivo › Subir cuaderno* y después *Entorno de ejecución › Ejecutar todo*. También corre en Jupyter con las dependencias de `requirements.txt`.

## Modelo

| Parámetro | Valor |
|---|---|
| Tramo | 80 km, 32 pulgadas, caudal nominal 1900 m³/h |
| Capacitancia del tramo C_L | 4,05 m³/bar |
| Descarga k_out | 34,5 m³/h por bar |
| Constante de tiempo del tramo τp | 7 min |
| Bomba y variador | K_B = 25 m³/h por %, τv = 10 s, velocidad de régimen U0 = 76 % |
| Transmisor | 0-100 bar → 4-20 mA (Kt = 0,16 mA/bar), τt = 1 s, ruido de 0,05 bar |
| Ciclo del PLC | 1 s |
| Controlador | Kp = 10 % por bar, Ki = 6,67 % por bar·min, Kv = 0 (PI) |
| QoS | 60 bar ± 1 % (±0,6 bar), 120 s tolerados fuera de banda, máximo de 70 bar |

## Resultados principales

- **Estabilidad (Routh):**
  - Con solo P, el lazo es estable para Kp < 658 % por bar.
  - Con Kp = 10, se vuelve inestable para Ki ≈ 62,5 % por bar·min.
  - Los polos dominantes del PI están en −0,64 ± 0,67j por minuto.
- **Margen de fase:** 55° con cruce en 1,16 rad/min. El lazo tolera unos 50 s de demora pura; por eso el PI queda en el PLC y la nube recibe los datos en forma asíncrona.

| Escenario | Resultado | QoS |
|---|---|---|
| 1. Baja de demanda del 15 %, PI | Pico de 0,7 bar; vuelve a la banda en 83 s, sin error final | Transitoria |
| 2. La misma baja, solo P | Queda 1,0 bar arriba (error permanente) | Falla a los 120 s |
| 3. Baja de 60 s, solo P | Sale 33 s de la banda y vuelve | Transitoria |
| 4. PID sin perturbación | La D transmite el ruido del transmisor a la bomba | Se cumple, con órdenes innecesarias |
| 5. Ki sobre el límite de Routh (Ki = 90) | Oscilación sostenida de +1,2 / −1,4 bar, bomba saturada | Sale de banda cada minuto |
| 6. Fuga de 200 m³/h | El PI la tapa; la supervisión para la bomba a los 14 s | Falla (parada) |
| 7. Transmisor falseado −15 bar | La presión real sube; disparo a 70 bar a los 6 min 11 s | Falla a los 144 s y parada |
| 8. Sin control (manual al 76 %) | La presión sube hacia 69,7 bar | Falla a los 120 s |
| 9. Solo I | Oscila ±1,4 bar con un período de 7,7 min | Falla en cada ciclo |
| 10. Solo D | Frena la subida y termina como sin control | Falla a los 120 s |

## Notación

Sigue a la cátedra:

| Símbolo | Qué es |
|---|---|
| θi | Referencia |
| θo | Variable controlada |
| f | Señal de medición |
| e = θi − f | Error |
| Gc, Ga, Gp, H | Controlador, actuador, proceso y elemento de medición |
| θd | Perturbación dentro del lazo (el lazo la corrige) |
| n | Perturbación en la realimentación (se propaga) |
| Kv | Constante de la acción derivativa |
