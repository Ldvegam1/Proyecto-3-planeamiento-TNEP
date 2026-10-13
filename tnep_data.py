"""
tnep_data.py
============
Base de datos comun y utilidades compartidas por los tres modelos de
expansion de transmision (LP de Villasana-Garver, M(0-1)LP big-M y MINLP).

Formato de la base de datos (una carpeta por caso, p.ej. data/garver/):
    buses.csv      bus, PG_MW, PL_MW, tipo          (despacho y demanda del ano horizonte)
    lines.csv      from, to, x_pu, cap_MW, ...      (red existente; una fila = un circuito)
    corridors.csv  id, from, to, cost_pu, x_pu, cap_MW, nmax, ...   (corredores candidatos)
    params.csv     key, value                        (base_MVA, ref_bus, bigM, theta_max)

Para escalar a un sistema mayor basta con crear otra carpeta con los mismos
cuatro archivos

  * Todo se trabaja en por unidad (pu) sobre base_MVA.
  * Un corredor i-j se orienta desde (i) hacia (j); el flujo
    positivo sale de i y entra a j. Delta_theta = theta_i - theta_j.
  * Flujo DC de un circuito: P = (theta_i - theta_j) / x.
"""

import csv
import os
import numpy as np


# ----------------------------------------------------------------------------
# Lectura de la base de datos
# ----------------------------------------------------------------------------
def _read_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return [row for row in csv.DictReader(fh)]


def load_case(folder):
    """Lee los cuatro CSV del caso y devuelve un diccionario con todo en pu."""
    params = {r["key"]: float(r["value"]) for r in _read_csv(os.path.join(folder, "params.csv"))}
    base = params.get("base_MVA", 100.0)

    buses_raw = _read_csv(os.path.join(folder, "buses.csv"))
    buses = [int(r["bus"]) for r in buses_raw]
    PG = {int(r["bus"]): float(r["PG_MW"]) / base for r in buses_raw}
    PL = {int(r["bus"]): float(r["PL_MW"]) / base for r in buses_raw}

    lines = []
    for idx, r in enumerate(_read_csv(os.path.join(folder, "lines.csv")), start=1):
        lines.append(dict(id=idx, f=int(r["from"]), t=int(r["to"]),
                          x=float(r["x_pu"]), cap=float(r["cap_MW"]) / base))

    corridors = []
    for r in _read_csv(os.path.join(folder, "corridors.csv")):
        corridors.append(dict(id=int(r["id"]), f=int(r["from"]), t=int(r["to"]),
                              cost=float(r["cost_pu"]), x=float(r["x_pu"]),
                              B=1.0 / float(r["x_pu"]), cap=float(r["cap_MW"]) / base,
                              nmax=int(r["nmax"])))

    ref = int(params.get("ref_bus", buses[0]))
    case = dict(name=os.path.basename(os.path.normpath(folder)), base=base,
                buses=buses, PG=PG, PL=PL, lines=lines, corridors=corridors,
                ref=ref, bigM=params.get("bigM", 1.0),
                theta_max=params.get("theta_max", np.pi), params=params)
    _validate(case)
    return case


def _validate(case):
    b = set(case["buses"])
    for l in case["lines"]:
        assert l["f"] in b and l["t"] in b, f"Linea {l} conecta barras inexistentes"
    for c in case["corridors"]:
        assert c["f"] in b and c["t"] in b, f"Corredor {c} conecta barras inexistentes"
    assert case["ref"] in b, "La barra de referencia no existe"
    gen, dem = sum(case["PG"].values()), sum(case["PL"].values())
    assert abs(gen - dem) < 1e-6, f"Despacho ({gen*case['base']:.1f} MW) != demanda ({dem*case['base']:.1f} MW)"


def corridor_label(c):
    return f"{c['f']}-{c['t']}"


# ----------------------------------------------------------------------------
# Matriz de susceptancias aumentada B^a(n) y flujo de carga DC
# ----------------------------------------------------------------------------
def susceptance_matrix(case, n_new=None):
    """B^a(n): red existente + n_k lineas en cada corredor k (Ec. 21 y 23)."""
    n_new = n_new or {}
    idx = {b: i for i, b in enumerate(case["buses"])}
    nb = len(case["buses"])
    B = np.zeros((nb, nb))

    def stamp(f, t, y):
        i, j = idx[f], idx[t]
        B[i, i] += y; B[j, j] += y
        B[i, j] -= y; B[j, i] -= y

    for l in case["lines"]:
        stamp(l["f"], l["t"], 1.0 / l["x"])
    for c in case["corridors"]:
        n = n_new.get(c["id"], 0)
        if n:
            stamp(c["f"], c["t"], n * c["B"])
    return B


def dc_power_flow(case, n_new=None):
    """theta = B^-1 (PG - PL) eliminando la barra de referencia (Sec. 6.1).
    Devuelve (theta dict, flujos existentes, flujos corredores) o None si
    la red queda en isla (B reducida singular)."""
    n_new = n_new or {}
    buses = case["buses"]
    B = susceptance_matrix(case, n_new)
    P = np.array([case["PG"][b] - case["PL"][b] for b in buses])
    keep = [i for i, b in enumerate(buses) if b != case["ref"]]
    Br = B[np.ix_(keep, keep)]
    if np.linalg.matrix_rank(Br) < len(keep):
        return None
    th = np.zeros(len(buses))
    th[keep] = np.linalg.solve(Br, P[keep])
    theta = {b: th[i] for i, b in enumerate(buses)}
    f_lines = {l["id"]: (theta[l["f"]] - theta[l["t"]]) / l["x"] for l in case["lines"]}
    f_corr = {c["id"]: n_new.get(c["id"], 0) * c["B"] * (theta[c["f"]] - theta[c["t"]])
              for c in case["corridors"]}
    return theta, f_lines, f_corr


def check_limits(case, n_new, tol=1e-6):
    """Verifica la factibilidad de un plan con el flujo DC. Devuelve (ok, lista de violaciones)."""
    res = dc_power_flow(case, n_new)
    if res is None:
        return False, ["Red en isla (alguna barra queda desconectada)"]
    theta, fl, fc = res
    viol = []
    for l in case["lines"]:
        if abs(fl[l["id"]]) > l["cap"] + tol:
            viol.append(f"Linea {l['f']}-{l['t']}: {abs(fl[l['id']])*case['base']:.2f} MW > {l['cap']*case['base']:.0f} MW")
    for c in case["corridors"]:
        n = n_new.get(c["id"], 0)
        if abs(fc[c["id"]]) > n * c["cap"] + tol:
            viol.append(f"Corredor {corridor_label(c)} ({n} lineas): {abs(fc[c['id']])*case['base']:.2f} MW > {n*c['cap']*case['base']:.0f} MW")
    return len(viol) == 0, viol


def plan_cost(case, n_new):
    return sum(c["cost"] * n_new.get(c["id"], 0) for c in case["corridors"])


# ----------------------------------------------------------------------------
# Reportes
# ----------------------------------------------------------------------------
def print_header(title):
    print("\n" + "=" * 78)
    print(title)
    print("=" * 78)


def print_plan(case, n_new):
    print(f"{'Corredor':<10}{'Lineas nuevas':>15}{'Costo unit.':>13}{'Costo':>10}")
    for c in case["corridors"]:
        n = n_new.get(c["id"], 0)
        print(f"{corridor_label(c):<10}{n:>15d}{c['cost']:>13.0f}{n*c['cost']:>10.0f}")
    print(f"{'TOTAL':<10}{sum(n_new.values()):>15d}{'':>13}{plan_cost(case, n_new):>10.0f}  pu")


def print_power_flow(case, n_new):
    res = dc_power_flow(case, n_new)
    if res is None:
        print("  (red en isla: no hay flujo DC)")
        return
    theta, fl, fc = res
    base = case["base"]
    print("Angulos de barra (rad):  " + "  ".join(f"th{b}={theta[b]:.3f}" for b in case["buses"]))
    print(f"{'Rama':<14}{'Flujo MW':>10}{'Cap. MW':>10}{'Carga %':>10}")
    for l in case["lines"]:
        print(f"{'L '+str(l['f'])+'-'+str(l['t']):<14}{fl[l['id']]*base:>10.2f}{l['cap']*base:>10.0f}"
              f"{100*abs(fl[l['id']])/l['cap']:>10.1f}")
    for c in case["corridors"]:
        n = n_new.get(c["id"], 0)
        if n:
            cap = n * c["cap"]
            print(f"{'C '+corridor_label(c)+f' (x{n})':<14}{fc[c['id']]*base:>10.2f}{cap*base:>10.0f}"
                  f"{100*abs(fc[c['id']])/cap:>10.1f}")


def write_csv(path, header, rows):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def default_case_path():
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "garver")
