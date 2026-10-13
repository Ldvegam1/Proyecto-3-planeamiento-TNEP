"""
tnep_milp_bigm.py  -  Modelo M(0-1)LP disyuntivo con linealizacion big-M (Villasana, 1984)
==========================================================================================
Implementacion en Gurobi (gurobipy) del modelo de la Seccion 6 del documento TNEP.

Variables de decision:
    w_ijk  in {0,1}   = 1 si se instalan EXACTAMENTE k lineas nuevas en el corredor i-j
    dth_ijk           diferencia angular "auxiliar" asociada a la alternativa k
    P_ijk             flujo del corredor i-j si se elige la alternativa k
    theta_j           angulos de barra (theta_ref = 0)

    min  sum_{ij in Omega} sum_k  k * c_ij * w_ijk                                   (Ec. 13)
    s.a. balance nodal con red existente + flujos P_ijk                              (Ec. 14)
         -P_jm^max <= P_jm <= P_jm^max       (lineas existentes)                     (Ec. 15)
         -k P_ij^max w_ijk <= P_ijk <= k P_ij^max w_ijk                              (Ec. 16)
         P_jm  = B_jm (theta_j - theta_m)                                            (Ec. 17)
         P_ijk = k * B_ij * dth_ijk                                                  (Ec. 18)
         M(sum_k w_ijk - 1) <= theta_i - theta_j - sum_k dth_ijk <= M(1 - sum_k w_ijk)  (Ec. 19, big-M)
         sum_k w_ijk <= 1                    (unicidad)                              (Ec. 20)

Con M = 1 rad (valor usado por Villasana); M se lee de params.csv para poder
ajustarlo en sistemas mayores (debe superar la maxima diferencia angular posible
entre los extremos de un corredor NO construido).

Uso:
    python tnep_milp_bigm.py
    python tnep_milp_bigm.py --case data/otro_sistema --M 2.0
"""

import argparse
import os
import gurobipy as gp
from gurobipy import GRB

import tnep_data as td


def solve_milp_bigm(case, M=None, verbose=False, time_limit=None):
    buses, lines, cors = case["buses"], case["lines"], case["corridors"]
    M = case["bigM"] if M is None else M
    base = case["base"]

    td.print_header(f"MODELO M(0-1)LP BIG-M DE VILLASANA  -  caso '{case['name']}'  (M = {M} rad)")

    m = gp.Model("TNEP_MILP_BigM")
    m.Params.OutputFlag = 1 if verbose else 0
    m.Params.MIPGap = 0.0
    if time_limit:
        m.Params.TimeLimit = time_limit

    # --- Conjuntos de indices (i-j, k) ---------------------------------------
    IJK = [(c["id"], k) for c in cors for k in range(1, c["nmax"] + 1)]
    C = {c["id"]: c for c in cors}

    # --- Variables -----------------------------------------------------------
    theta = m.addVars(buses, lb=-GRB.INFINITY, name="theta")
    theta[case["ref"]].lb = theta[case["ref"]].ub = 0.0
    w = m.addVars(IJK, vtype=GRB.BINARY, name="w")
    dth = m.addVars(IJK, lb=-GRB.INFINITY, name="dtheta")
    P = m.addVars(IJK, lb=-GRB.INFINITY, name="P")
    PF = m.addVars([l["id"] for l in lines], lb=-GRB.INFINITY, name="PF")

    # --- Objetivo (Ec. 13) ---------------------------------------------------
    m.setObjective(gp.quicksum(k * C[i]["cost"] * w[i, k] for i, k in IJK), GRB.MINIMIZE)

    # --- Red existente: Ec. 15 y 17 ------------------------------------------
    for l in lines:
        m.addConstr(PF[l["id"]] == (theta[l["f"]] - theta[l["t"]]) / l["x"], name=f"flow_ex_{l['id']}")
        PF[l["id"]].lb, PF[l["id"]].ub = -l["cap"], l["cap"]

    # --- Corredores candidatos: Ec. 16, 18, 19, 20 ---------------------------
    for c in cors:
        i = c["id"]
        K = range(1, c["nmax"] + 1)
        for k in K:
            m.addConstr(P[i, k] == k * c["B"] * dth[i, k], name=f"flow_new_{i}_{k}")              # (18)
            m.addConstr(P[i, k] <= k * c["cap"] * w[i, k], name=f"cap_up_{i}_{k}")                # (16)
            m.addConstr(P[i, k] >= -k * c["cap"] * w[i, k], name=f"cap_dn_{i}_{k}")               # (16)
        sw = gp.quicksum(w[i, k] for k in K)
        gap = theta[c["f"]] - theta[c["t"]] - gp.quicksum(dth[i, k] for k in K)
        m.addConstr(gap <= M * (1 - sw), name=f"bigM_up_{i}")                                     # (19)
        m.addConstr(gap >= -M * (1 - sw), name=f"bigM_dn_{i}")                                    # (19)
        m.addConstr(sw <= 1, name=f"unicidad_{i}")                                                # (20)

    # --- Balance nodal (Ec. 14): salidas - entradas = PG - PL -----------------
    for j in buses:
        out_ex = gp.quicksum(PF[l["id"]] for l in lines if l["f"] == j)
        in_ex = gp.quicksum(PF[l["id"]] for l in lines if l["t"] == j)
        out_new = gp.quicksum(P[i, k] for i, k in IJK if C[i]["f"] == j)
        in_new = gp.quicksum(P[i, k] for i, k in IJK if C[i]["t"] == j)
        m.addConstr(out_ex - in_ex + out_new - in_new == case["PG"][j] - case["PL"][j], name=f"bal_{j}")

    m.optimize()
    if m.Status not in (GRB.OPTIMAL, GRB.TIME_LIMIT) or m.SolCount == 0:
        raise RuntimeError(f"MILP sin solucion (status {m.Status})")

    # --- Lectura de la solucion ---------------------------------------------
    n_new = {c["id"]: 0 for c in cors}
    print("Variables binarias activas:")
    for i, k in IJK:
        if w[i, k].X > 0.5:
            n_new[i] = k
            print(f"  w_{C[i]['f']}{C[i]['t']}{k} = 1   ({k} linea(s) en el corredor {td.corridor_label(C[i])})")
    print(f"\nCosto total de la expansion = {m.ObjVal:.2f} pu   "
          f"(gap = {m.MIPGap:.2e}, nodos B&B = {int(m.NodeCount)}, tiempo = {m.Runtime:.3f} s)")
    print("Angulos theta (rad): " + "  ".join(f"th{b}={theta[b].X:.3f}" for b in buses))
    print("Angulos auxiliares dtheta_ijk (rad):")
    print("  " + "  ".join(f"dth{C[i]['f']}{C[i]['t']}{k}={dth[i, k].X:.3f}" for i, k in IJK))

    td.print_header("RESULTADO FINAL - M(0-1)LP BIG-M")
    td.print_plan(case, n_new)
    print()
    td.print_power_flow(case, n_new)
    ok, viol = td.check_limits(case, n_new)
    print(f"\nVerificacion flujo DC (theta = B^-1 (PG-PL)): {'FACTIBLE' if ok else 'INFACTIBLE'}")
    for v in viol:
        print("   - " + v)

    return dict(n_new=n_new, obj=m.ObjVal, theta={b: theta[b].X for b in buses},
                dth={(i, k): dth[i, k].X for i, k in IJK}, runtime=m.Runtime,
                nvars=m.NumVars, ncons=m.NumConstrs, nbin=m.NumBinVars)


def main():
    ap = argparse.ArgumentParser(description="TNEP - M(0-1)LP big-M de Villasana (Gurobi)")
    ap.add_argument("--case", default=td.default_case_path())
    ap.add_argument("--M", type=float, default=None, help="constante big-M en rad (por defecto params.csv)")
    ap.add_argument("--verbose", action="store_true")
    ap.add_argument("--time-limit", type=float, default=None)
    args = ap.parse_args()

    case = td.load_case(args.case)
    sol = solve_milp_bigm(case, args.M, args.verbose, args.time_limit)

    out = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results", case["name"])
    cors = {c["id"]: c for c in case["corridors"]}
    td.write_csv(os.path.join(out, "milp_plan.csv"), ["corredor", "lineas", "costo"],
                 [[td.corridor_label(c), sol["n_new"][i], sol["n_new"][i] * c["cost"]] for i, c in cors.items()])
    td.write_csv(os.path.join(out, "milp_angulos.csv"), ["variable", "valor_rad"],
                 [[f"theta_{b}", round(v, 6)] for b, v in sol["theta"].items()] +
                 [[f"dtheta_{cors[i]['f']}{cors[i]['t']}{k}", round(v, 6)] for (i, k), v in sol["dth"].items()])
    print(f"\nResultados guardados en {out}")


if __name__ == "__main__":
    main()
