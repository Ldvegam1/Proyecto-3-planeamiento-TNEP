"""
tnep_lp.py  -  Metodo LP interactivo de Villasana-Garver (Villasana et al., 1985)
==================================================================================
Implementacion en Gurobi (gurobipy)

En cada iteracion se resuelve el LP (Ec. 3-7):

    min   sum_k c_Dk (P+_Dk + P-_Dk)
    s.a.  P_G - B^a theta + K P+_D - K P-_D = P_L          (balance, Ec. 4)
          P_F = X H theta                                   (flujos DC red "existente", Ec. 5)
          -P_F^max <= P_F <= P_F^max                        (limites, Ec. 6)
          P+_D, P-_D >= 0                                   (Ec. 7)
          theta_ref = 0

Los corredores actuan como "transporte" ficticio (sin 2da ley de Kirchhoff)
con costo unitario c_Dk: el LP indica cuanto flujo falta por mover y por donde.
Procedimiento de redondeo (round-up): se agrega UNA linea en el corredor con
mayor flujo P_Dk = P+_Dk + P-_Dk; esa linea pasa a la red "existente" (entra a
B^a y a los limites P_F^max) y se repite hasta que P+_D = P-_D = 0.

Uso:
    python tnep_lp.py                 # caso Garver (data/garver)
    python tnep_lp.py --case data/otro_sistema
"""

import argparse
import os
import gurobipy as gp
from gurobipy import GRB

import tnep_data as td


def build_and_solve_lp(case, n_added, verbose=False):
    """Construye y resuelve el LP de una iteracion. n_added[k] = lineas ya
    agregadas al corredor k (que ahora forman parte de la red existente)."""
    buses, cors = case["buses"], case["corridors"]

    m = gp.Model("TNEP_LP_Villasana")
    m.Params.OutputFlag = 1 if verbose else 0
    m.Params.Method = 0          # Simplex primal

    # --- Variables -----------------------------------------------------------
    theta = m.addVars(buses, lb=-GRB.INFINITY, ub=GRB.INFINITY, name="theta")
    theta[case["ref"]].lb = theta[case["ref"]].ub = 0.0
    kid = [c["id"] for c in cors]
    Pp = m.addVars(kid, lb=0.0, name="Pplus")    # P+_D (sentido from->to)
    Pm = m.addVars(kid, lb=0.0, name="Pminus")   # P-_D (sentido to->from)

    # --- Circuitos de la red "existente" de esta iteracion -------------------
    #     lineas originales + lineas ya instaladas en los corredores
    circuits = [(f"L{l['f']}-{l['t']}", l["f"], l["t"], 1.0 / l["x"], l["cap"]) for l in case["lines"]]
    for c in cors:
        n = n_added.get(c["id"], 0)
        if n:
            circuits.append((f"N{c['f']}-{c['t']}", c["f"], c["t"], n * c["B"], n * c["cap"]))

    # P_F = X H theta  (Ec. 5) y limites (Ec. 6)
    PF = {}
    for name, f, t, b, cap in circuits:
        PF[name] = b * (theta[f] - theta[t])
        m.addConstr(PF[name] <= cap, name=f"cap_up_{name}")
        m.addConstr(PF[name] >= -cap, name=f"cap_dn_{name}")

    # --- Balance nodal (Ec. 4): (B^a theta)_j - [K(P+ - P-)]_j = PG_j - PL_j ---
    #     K_jk = +1 si j es el extremo 'to' del corredor k, -1 si es el 'from'
    balance = {}
    for j in buses:
        Btheta = gp.quicksum(PF[nm] for nm, f, t, b, cap in circuits if f == j) \
               - gp.quicksum(PF[nm] for nm, f, t, b, cap in circuits if t == j)
        KP = gp.quicksum(Pp[c["id"]] - Pm[c["id"]] for c in cors if c["t"] == j) \
           - gp.quicksum(Pp[c["id"]] - Pm[c["id"]] for c in cors if c["f"] == j)
        balance[j] = m.addConstr(Btheta - KP == case["PG"][j] - case["PL"][j], name=f"bal_{j}")

    # --- Funcion objetivo (Ec. 3) --------------------------------------------
    m.setObjective(gp.quicksum(c["cost"] * (Pp[c["id"]] + Pm[c["id"]]) for c in cors), GRB.MINIMIZE)
    m.optimize()

    if m.Status != GRB.OPTIMAL:
        raise RuntimeError(f"LP no optimo (status {m.Status})")

    return dict(obj=m.ObjVal,
                theta={b: theta[b].X for b in buses},
                Pp={k: Pp[k].X for k in kid},
                Pm={k: Pm[k].X for k in kid},
                lmp={j: balance[j].Pi for j in buses})


def solve_lp_villasana(case, verbose=False, tol=1e-6, max_iter=1000):
    cors = case["corridors"]
    base = case["base"]
    n_added = {c["id"]: 0 for c in cors}
    history = []

    td.print_header(f"MODELO LP INTERACTIVO DE VILLASANA-GARVER  -  caso '{case['name']}'")
    for it in range(1, max_iter + 1):
        sol = build_and_solve_lp(case, n_added, verbose)
        PD = {k: sol["Pp"][k] + sol["Pm"][k] for k in sol["Pp"]}

        print(f"\n--- Iteracion {it} ---")
        print("  Lineas instaladas hasta ahora: " +
              ", ".join(f"{td.corridor_label(c)}:{n_added[c['id']]}" for c in cors))
        print("  theta (rad) = [" + " ".join(f"{sol['theta'][b]:.3f}" for b in case["buses"]) + "]")
        print("  P+_D (pu)   = [" + " ".join(f"{sol['Pp'][c['id']]:.4f}" for c in cors) + "]")
        print("  P-_D (pu)   = [" + " ".join(f"{sol['Pm'][c['id']]:.4f}" for c in cors) + "]")
        print(f"  Funcion objetivo = {sol['obj']:.4f}")

        record = dict(iter=it, obj=sol["obj"], PD=PD, theta=sol["theta"], added=None)

        # Criterio de parada: no hay flujo en corredores ficticios
        if max(PD.values()) <= tol:
            print("  Ya no hay congestion en los corredores -> OPTIMO alcanzado.")
            history.append(record)
            break

        # Procedimiento de redondeo: corredor "mas necesario" con capacidad remanente
        candidates = [c for c in cors if n_added[c["id"]] < c["nmax"] and PD[c["id"]] > tol]
        if not candidates:
            raise RuntimeError("Hay flujo en corredores que ya alcanzaron nmax: problema infactible.")
        best = max(candidates, key=lambda c: PD[c["id"]])
        n_added[best["id"]] += 1
        record["added"] = td.corridor_label(best)
        history.append(record)
        print(f"  Corredor mas necesario: {td.corridor_label(best)} "
              f"(P_D = {PD[best['id']]:.4f} pu = {PD[best['id']]*base:.2f} MW) -> se agrega 1 linea "
              f"(total en el corredor: {n_added[best['id']]})")
    else:
        raise RuntimeError("Se alcanzo el maximo de iteraciones")

    td.print_header("RESULTADO FINAL - LP VILLASANA-GARVER")
    td.print_plan(case, n_added)
    print()
    td.print_power_flow(case, n_added)
    ok, viol = td.check_limits(case, n_added)
    print(f"\nVerificacion flujo DC: {'FACTIBLE' if ok else 'INFACTIBLE'}")
    for v in viol:
        print("   - " + v)
    return n_added, history


def main():
    ap = argparse.ArgumentParser(description="TNEP - LP interactivo de Villasana-Garver (Gurobi)")
    ap.add_argument("--case", default=td.default_case_path(), help="carpeta con buses/lines/corridors/params.csv")
    ap.add_argument("--verbose", action="store_true", help="mostrar log de Gurobi")
    args = ap.parse_args()

    case = td.load_case(args.case)
    n_added, history = solve_lp_villasana(case, args.verbose)

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", case["name"])
    cors = case["corridors"]
    td.write_csv(os.path.join(out, "lp_iteraciones.csv"),
                 ["iteracion", "funcion_objetivo"] + [f"PD_{td.corridor_label(c)}_pu" for c in cors]
                 + [f"theta_{b}" for b in case["buses"]] + ["linea_agregada"],
                 [[h["iter"], round(h["obj"], 4)] + [round(h["PD"][c["id"]], 4) for c in cors]
                  + [round(h["theta"][b], 4) for b in case["buses"]] + [h["added"] or "-"] for h in history])
    td.write_csv(os.path.join(out, "lp_plan.csv"), ["corredor", "lineas", "costo"],
                 [[td.corridor_label(c), n_added[c["id"]], n_added[c["id"]] * c["cost"]] for c in cors])
    print(f"\nResultados guardados en {out}")


if __name__ == "__main__":
    main()
