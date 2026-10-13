"""
run_all.py  -  Ejecuta los tres modelos sobre la MISMA base de datos y compara resultados.

    python run_all.py                       # caso Garver
    python run_all.py --case data/otro      # cualquier otro caso con el mismo formato
    python run_all.py --brute-force         # agrega la validacion exhaustiva (solo casos pequenos)
"""

import argparse
import os
import time

import tnep_data as td
from tnep_lp import solve_lp_villasana
from tnep_milp_bigm import solve_milp_bigm
from tnep_minlp import solve_minlp


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default=td.default_case_path())
    ap.add_argument("--brute-force", action="store_true")
    args = ap.parse_args()

    case = td.load_case(args.case)
    cors = case["corridors"]
    summary = []

    t0 = time.perf_counter()
    n_lp, hist = solve_lp_villasana(case)
    summary.append(("LP Villasana-Garver (iterativo)", n_lp, time.perf_counter() - t0, f"{len(hist)} LPs resueltos"))

    t0 = time.perf_counter()
    milp = solve_milp_bigm(case)
    summary.append(("M(0-1)LP big-M", milp["n_new"], time.perf_counter() - t0,
                    f"{milp['nbin']} binarias, {milp['nvars']} var, {milp['ncons']} restr."))

    t0 = time.perf_counter()
    minlp = solve_minlp(case)
    summary.append(("MINLP bilineal", minlp["n_new"], time.perf_counter() - t0,
                    f"{minlp['nint']} enteras, {minlp['nvars']} var, {minlp['ncons']} restr."))

    if args.brute_force:
        from validar_fuerza_bruta import brute_force
        t0 = time.perf_counter()
        bf = brute_force(case)
        summary.append(("Fuerza bruta (validacion)", bf, time.perf_counter() - t0, "enumeracion + flujo DC"))

    td.print_header("COMPARACION DE METODOS")
    hdr = f"{'Metodo':<32}" + "".join(f"{td.corridor_label(c):>6}" for c in cors) + f"{'Costo':>8}{'Factible':>10}{'t [s]':>8}"
    print(hdr)
    rows = []
    for name, n_new, t, info in summary:
        ok, _ = td.check_limits(case, n_new)
        cost = td.plan_cost(case, n_new)
        print(f"{name:<32}" + "".join(f"{n_new[c['id']]:>6d}" for c in cors)
              + f"{cost:>8.0f}{'si' if ok else 'NO':>10}{t:>8.3f}   {info}")
        rows.append([name] + [n_new[c["id"]] for c in cors] + [cost, ok, round(t, 4), info])

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", case["name"])
    td.write_csv(os.path.join(out, "comparacion_metodos.csv"),
                 ["metodo"] + [f"n_{td.corridor_label(c)}" for c in cors] + ["costo_pu", "factible", "tiempo_s", "info"],
                 rows)
    print(f"\nResumen guardado en {os.path.join(out, 'comparacion_metodos.csv')}")


if __name__ == "__main__":
    main()
