"""
tnep_minlp.py  -  Modelo MINLP bilineal (entero x angulo) de expansion de transmision
=====================================================================================
Implementacion en Gurobi (gurobipy) del modelo de la Seccion 7 del documento TNEP.

Variables de decision:
    n_ij in Z, 0 <= n_ij <= n_ij^max   numero de lineas nuevas en el corredor i-j
    theta_j                             angulos de barra (theta_ref = 0)
    P_ij                                flujo total del corredor i-j

    min  sum_{ij in Omega} c_ij * n_ij                                    (Ec. 24)
    s.a. sum_m B^a_jm(n) theta_m = PG_j - PL_j                            (Ec. 25, B^a depende de n)
         P_jm = B_jm (theta_j - theta_m)        (existentes)              (Ec. 27)
         P_ij = n_ij * B_ij * (theta_i - theta_j)   (BILINEAL)            (Ec. 26 y 28)
         -P_jm^max <= P_jm <= P_jm^max                                    (Ec. 29)
         -n_ij P_ij^max <= P_ij <= n_ij P_ij^max                          (Ec. 30)
         theta_ref = 0                                                    (Ec. 31)

El balance (25) se escribe en forma expandida: (B theta)_j de la red existente
+ flujos de los corredores, que es exactamente sum_m B^a_jm(n) theta_m.

No hay big-M ni restriccion de unicidad. El producto n_ij * theta es no convexo;
Gurobi lo resuelve con ramificacion espacial (NonConvex = 2), por lo que, a
diferencia del GRG No Lineal de Excel, la solucion obtenida es un OPTIMO GLOBAL
certificado (gap = 0). Para ello requiere cotas finitas en theta (theta_max).

Uso:
    python tnep_minlp.py
    python tnep_minlp.py --case data/otro_sistema --theta-max 1.5708
"""

import argparse
import os
import gurobipy as gp
from gurobipy import GRB

import tnep_data as td


def solve_minlp(case, theta_max=None, verbose=False, time_limit=None):
    buses, lines, cors = case["buses"], case["lines"], case["corridors"]
    tmax = case["theta_max"] if theta_max is None else theta_max

    td.print_header(f"MODELO MINLP BILINEAL  -  caso '{case['name']}'  (|theta| <= {tmax:.4f} rad)")

    m = gp.Model("TNEP_MINLP")
    m.Params.OutputFlag = 1 if verbose else 0
    m.Params.NonConvex = 2        # habilita restricciones cuadraticas no convexas (bilineales)
    m.Params.MIPGap = 0.0
    if time_limit:
        m.Params.TimeLimit = time_limit

    kid = [c["id"] for c in cors]
    C = {c["id"]: c for c in cors}

    # --- Variables -----------------------------------------------------------
    theta = m.addVars(buses, lb=-tmax, ub=tmax, name="theta")
    theta[case["ref"]].lb = theta[case["ref"]].ub = 0.0
    n = m.addVars(kid, vtype=GRB.INTEGER, lb=0, ub={k: C[k]["nmax"] for k in kid}, name="n")
    P = m.addVars(kid, lb=-GRB.INFINITY, name="P")
    PF = m.addVars([l["id"] for l in lines], lb=-GRB.INFINITY, name="PF")

    # --- Objetivo (Ec. 24) ---------------------------------------------------
    m.setObjective(gp.quicksum(C[k]["cost"] * n[k] for k in kid), GRB.MINIMIZE)

    # --- Lineas existentes (Ec. 27 y 29) -------------------------------------
    for l in lines:
        m.addConstr(PF[l["id"]] == (theta[l["f"]] - theta[l["t"]]) / l["x"], name=f"flow_ex_{l['id']}")
        PF[l["id"]].lb, PF[l["id"]].ub = -l["cap"], l["cap"]

    # --- Corredores candidatos: producto bilineal (Ec. 28) y limite (Ec. 30) --
    for k in kid:
        c = C[k]
        m.addQConstr(P[k] == c["B"] * n[k] * (theta[c["f"]] - theta[c["t"]]), name=f"bilin_{k}")
        m.addConstr(P[k] <= c["cap"] * n[k], name=f"cap_up_{k}")
        m.addConstr(P[k] >= -c["cap"] * n[k], name=f"cap_dn_{k}")

    # --- Balance nodal con B^a(n) (Ec. 25) -----------------------------------
    for j in buses:
        Bth = gp.quicksum(PF[l["id"]] for l in lines if l["f"] == j) \
            - gp.quicksum(PF[l["id"]] for l in lines if l["t"] == j) \
            + gp.quicksum(P[k] for k in kid if C[k]["f"] == j) \
            - gp.quicksum(P[k] for k in kid if C[k]["t"] == j)
        m.addConstr(Bth == case["PG"][j] - case["PL"][j], name=f"bal_{j}")

    m.optimize()
    if m.Status not in (GRB.OPTIMAL, GRB.TIME_LIMIT) or m.SolCount == 0:
        raise RuntimeError(f"MINLP sin solucion (status {m.Status})")

    n_new = {k: int(round(n[k].X)) for k in kid}
    print("Solucion entera: " + ",  ".join(f"n_{C[k]['f']}{C[k]['t']} = {n_new[k]}" for k in kid))
    print(f"Costo total de la expansion = {m.ObjVal:.2f} pu   "
          f"(gap = {m.MIPGap:.2e}, nodos B&B = {int(m.NodeCount)}, tiempo = {m.Runtime:.3f} s)")
    print("theta (rad) = [" + " ".join(f"{theta[b].X:.3f}" for b in buses) + "]")

    td.print_header("RESULTADO FINAL - MINLP")
    td.print_plan(case, n_new)
    print()
    td.print_power_flow(case, n_new)
    ok, viol = td.check_limits(case, n_new)
    print(f"\nVerificacion flujo DC: {'FACTIBLE' if ok else 'INFACTIBLE'}")
    for v in viol:
        print("   - " + v)

    return dict(n_new=n_new, obj=m.ObjVal, theta={b: theta[b].X for b in buses},
                P={k: P[k].X for k in kid}, runtime=m.Runtime,
                nvars=m.NumVars, ncons=m.NumConstrs + m.NumQConstrs, nint=m.NumIntVars)


def main():
    ap = argparse.ArgumentParser(description="TNEP - MINLP bilineal (Gurobi NonConvex)")
    ap.add_argument("--case", default=td.default_case_path())
    ap.add_argument("--theta-max", type=float, default=None)
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--time-limit", type=float, default=None)
    args = ap.parse_args()

    case = td.load_case(args.case)
    sol = solve_minlp(case, args.theta_max, args.verbose, args.time_limit)

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", case["name"])
    td.write_csv(os.path.join(out, "minlp_plan.csv"), ["corredor", "lineas", "costo", "flujo_MW"],
                 [[td.corridor_label(c), sol["n_new"][c["id"]], sol["n_new"][c["id"]] * c["cost"],
                   round(sol["P"][c["id"]] * case["base"], 3)] for c in case["corridors"]])
    print(f"\nResultados guardados en {out}")


if __name__ == "__main__":
    main()
