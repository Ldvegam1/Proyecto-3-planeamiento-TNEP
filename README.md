# Expansión de Sistemas de Transmisión (TEP) en Gurobi — Caso Garver (1970)

Implementación en **Python + Gurobi (gurobipy)** de los tres modelos descritos en el documento
*Expansión de Sistemas de Transmisión — Formulación Clásica LP, MILP, MINLP* (P. M. De Oliveira-De Jesus,
Uniandes), originalmente resueltos con el Solver de Excel.

| Script | Modelo | Sección del documento | Rúbrica |
|---|---|---|---|
| `tnep_lp.py` | LP interactivo de Villasana-Garver (transporte + redondeo) | 5 (Ec. 3–12) | 30 % |
| `tnep_minlp.py` | MINLP bilineal `n_ij · B_ij · Δθ` | 7 (Ec. 23–31) | 30 % |
| `tnep_milp_bigm.py` | M(0-1)LP disyuntivo con big-M | 6 (Ec. 13–20) | 40 % |
| `validar_fuerza_bruta.py` | Enumeración exhaustiva + flujo DC (validación) | 4.1 | — |
| `run_all.py` | Corre todo y genera la tabla comparativa | — | — |
| `tnep_data.py` | Lectura de la base de datos, `B^a(n)`, flujo DC, reportes | — | — |

## Requisitos

```bash
pip install gurobipy numpy
```
La licencia restringida que trae `gurobipy` (≤ 2000 variables/restricciones) basta para Garver.
Para sistemas grandes use la licencia académica gratuita de Gurobi.

## Ejecución

```bash
python run_all.py --brute-force        # los tres modelos + validación + comparación
python tnep_lp.py                      # solo LP (muestra las 8 iteraciones)
python tnep_milp_bigm.py               # solo M(0-1)LP big-M
python tnep_minlp.py                   # solo MINLP
python validar_fuerza_bruta.py --orden 2,1,3,4   # numeración de filas igual a la del documento
```
Opciones comunes: `--case <carpeta>` (otro sistema), `--verbose` (log de Gurobi),
`--time-limit <s>`; además `--M` en el big-M y `--theta-max` en el MINLP.

## Base de datos única (escalable)

Los tres scripts leen la **misma** carpeta `data/garver/`:

| Archivo | Columnas | Contenido (Garver) |
|---|---|---|
| `buses.csv` | `bus, PG_MW, PL_MW, tipo` | Tabla 1: despacho y demanda del año horizonte |
| `lines.csv` | `from, to, x_pu, cap_MW, length_mi` | Tabla 2: red existente (una fila = un circuito) |
| `corridors.csv` | `id, from, to, cost_pu, x_pu, cap_MW, nmax, length_mi` | Tabla 3 + reactancias de la Sec. 6.2 |
| `params.csv` | `key, value` | `base_MVA=100`, `ref_bus=5`, `bigM=1`, `theta_max=π` |

Ningún script tiene datos de Garver escritos en el código: todos los conjuntos (barras, líneas,
corredores, alternativas `k = 1..nmax`) se construyen a partir de los CSV. Para otro sistema basta
crear `data/<sistema>/` con los mismos cuatro archivos y ejecutar `--case data/<sistema>`.

## Resultados (caso Garver)

| Método | 6-2 | 6-4 | 6-5 | 3-5 | Costo [pu] | Factible |
|---|---|---|---|---|---|---|
| LP Villasana-Garver (8 LPs) | 4 | 2 | 0 | 1 | 1000 | sí |
| M(0-1)LP big-M (14 binarias) | 4 | 2 | 0 | 1 | 1000 | sí |
| MINLP bilineal (4 enteras) | 4 | 2 | 0 | 1 | 1000 | sí |
| Fuerza bruta (400 combinaciones) | 4 | 2 | 0 | 1 | 1000 | sí |

Ángulos del sistema expandido (rad): θ = [0.106 0.311 0.187 0.296 0.000 0.579], idénticos a la Ec. 34.
En el big-M, Δθ₆₂₄ = 0.268, Δθ₆₄₂ = 0.282, Δθ₃₅₁ = 0.187, igual que en la Sec. 6.2.

Secuencia del LP (función objetivo por iteración): 878.75 → 728.75 → 578.75 → 456.625 → 306.625 →
156.625 → 73.75 → 0, agregando líneas en 6-4, 6-4, 6-2, 6-2, 6-2, 3-5, 6-2.

Los archivos de salida quedan en `results/garver/` (`lp_iteraciones.csv`, `lp_plan.csv`,
`milp_plan.csv`, `milp_angulos.csv`, `minlp_plan.csv`, `comparacion_metodos.csv`, `salida_completa.txt`).

## Notas para el informe técnico

1. **LP, iteración 7.** Gurobi obtiene 73.75 (P_D,62 = 0.4916, resto 0), mientras que la hoja
   `iteración 7` del Excel del curso reporta 74.25 (P_D,62 = 0.489 y P_D,65 = 0.0029). Esa hoja está
   configurada con el motor **GRG Nonlinear** (las iteraciones 2–6 usan Simplex LP), que se detiene
   por tolerancia antes del vértice óptimo. Como el LP es convexo, 73.75 es el óptimo exacto. La
   decisión (agregar una línea en 6-2) y el plan final no cambian.
2. **Convención de signos.** Cada corredor se orienta `from → to` (6→2, 6→4, 6→5, 3→5); el flujo
   positivo sale de `from`. El balance se escribe como "salidas − entradas = PG − PL", que equivale
   a `B^a θ = PG − PL + K·P_D` del documento.
3. **Big-M.** Se usa M = 1 rad como Villasana. Cuando un corredor no se construye, la Ec. 19 impone
   |θi − θj| ≤ M; en sistemas mayores M debe ser mayor que la máxima apertura angular posible
   (por ejemplo, la suma de x·Pmax a lo largo del camino más largo entre los extremos del corredor).
   M demasiado grande debilita la relajación LP; M demasiado pequeño recorta soluciones válidas.
4. **MINLP.** El producto entero × continuo `n_ij · θ` es no convexo. Gurobi (`NonConvex = 2`) lo
   resuelve por ramificación espacial y **certifica el óptimo global** (gap = 0), a diferencia del
   GRG No Lineal de Excel, que solo garantiza óptimos locales. Para ello se acotan los ángulos con
   `theta_max`.
5. **Tamaño de los modelos.** Big-M: Σ nmax binarias + igual número de Δθ y P (54 variables,
   66 restricciones). MINLP: una entera por corredor (20 variables, 24 restricciones). El LP necesita
   tantas resoluciones como líneas nuevas (+1), pero cada LP es pequeño.
6. **Fuerza bruta.** El documento numera las combinaciones con los dígitos en el orden
   (6-4, 6-2, 6-5, 3-5): la fila 226 (`2401`) es factible y la 354 (`4201`) viola la capacidad del
   corredor 6-2, como indica la Sec. 4.1. Hay además una tercera combinación de costo 1000
   (`3301`, fila 290) que también es infactible por 6-2.

## Comparación con los Excel

| Elemento | Excel | Python/Gurobi | ¿Coincide? |
|---|---|---|---|
| LP, iteraciones 1–6 (FO, P_D, θ) | 878.75, 728.75, 578.75, 456.625, 306.625, 156.625 | idénticos | sí |
| LP, iteración 7 | 74.2495 (GRG) | 73.7473 (simplex) | misma decisión |
| LP, plan y θ finales (hoja OPTIMUM) | 4-2-0-1, θ = [0.106 0.311 0.187 0.2965 0 0.5787] | idénticos | sí |
| MILP big-M: w activas, Δθ, costo | w624, w642, w351; 0.2677/0.2822/0.1870; 1000 | idénticos | sí |
| MILP: alternativas k | 6-2 y 6-4 con k = 1..5 (16 binarias) | n_max = 4 de la Tabla 3 (14 binarias) | mismo óptimo* |
| MINLP: n y costo | 4-2-0-1, 1000 | idénticos | sí |
| MINLP: referencia angular | θ₁ = 0 | θ₅ = 0 | θ desplazado 0.106 rad; flujos idénticos |
| MINLP: cota n_max | sin cota superior | n_max de la Tabla 3 | mismo óptimo* |

\* Verificado corriendo el Python con k ≤ 5 y con n sin cota: el resultado sigue siendo 4-2-0-1 con costo 1000.

Observaciones sobre los Excel:
- En `C10_E3_Garver_MINLP` el Solver tiene guardada la restricción `O14 >= 1301` (costo ≥ 1301).
  Si se presiona *Resolver* sin quitarla, buscará otro plan y no el de costo 1000 que muestra la hoja.
- En `C10_E2_Garver_MILP`, la celda W26 (flujo F651 en el balance de la barra 6) usa `G39` (Δθ645)
  en lugar de `G40` (Δθ651). No afecta el óptimo porque el corredor 6-5 no se construye, pero muestra
  la ventaja de generar las ecuaciones a partir de la base de datos, como hacen los scripts.
