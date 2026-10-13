"""
validar_fuerza_bruta.py  -  Validacion por enumeracion exhaustiva (Seccion 4.1 del documento TNEP)
==================================================================================================
Enumera todas las combinaciones de lineas nuevas por corredor, ejecuta un flujo
de carga DC para cada una y reporta la de minimo costo que no viola limites.
Solo es viable en sistemas pequenos (prod(nmax+1) combinaciones; Garver: 400),
por eso sirve como referencia independiente para validar LP, MILP y MINLP.

La numeracion de "fila" sigue el orden lexicografico del documento
(0000 = fila 1, ..., ultima = fila prod(nmax+1)). Por defecto los digitos siguen
el orden de corridors.csv; el documento TNEP usa el orden 6-4, 6-2, 6-5, 3-5
(ids 2,1,3,4), que se reproduce con:   python validar_fuerza_bruta.py --orden 2,1,3,4
"""

import argparse
import itertools
import tnep_data as td


def brute_force(case, top=5, orden=None):
    cors = case["corridors"]
    if orden:
        byid = {c["id"]: c for c in cors}
        cors = [byid[i] for i in orden]
    ranges = [range(c["nmax"] + 1) for c in cors]
    total = 1
    for r in ranges:
        total *= len(r)

    td.print_header(f"VALIDACION POR FUERZA BRUTA  -  caso '{case['name']}'  ({total} combinaciones)")
    results = []
    for row, combo in enumerate(itertools.product(*ranges), start=1):
        n_new = {c["id"]: v for c, v in zip(cors, combo)}
        ok, viol = td.check_limits(case, n_new)
        results.append((td.plan_cost(case, n_new), row, combo, ok, viol))

    feas = sorted([r for r in results if r[3]], key=lambda r: (r[0], r[1]))
    if not feas:
        print("No existe ninguna combinacion factible.")
        return None
    best_cost = feas[0][0]

    print(f"Combinaciones factibles: {len(feas)} de {total}")
    print(f"\nCombinaciones con el costo minimo factible ({best_cost:.0f} pu):")
    print(f"{'Fila':>5}  {'Combinacion':<12}{'Costo':>8}  Estado")
    for cost, row, combo, ok, viol in sorted(results, key=lambda r: (r[0], r[1])):
        if abs(cost - best_cost) > 1e-9:
            continue
        code = "".join(str(v) for v in combo)
        print(f"{row:>5}  {code:<12}{cost:>8.0f}  {'FACTIBLE' if ok else 'INFACTIBLE: ' + '; '.join(viol)}")

    print(f"\nLas {top} mejores combinaciones factibles:")
    for cost, row, combo, ok, viol in feas[:top]:
        print(f"  fila {row:>4}  " + "  ".join(f"n{c['f']}{c['t']}={v}" for c, v in zip(cors, combo))
              + f"   costo = {cost:.0f}")

    best = {c["id"]: v for c, v in zip(cors, feas[0][2])}
    print("\nOptimo por fuerza bruta:")
    td.print_plan(case, best)
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", default=td.default_case_path())
    ap.add_argument("--orden", default=None, help="ids de corredores en el orden de los digitos, p.ej. 2,1,3,4")
    args = ap.parse_args()
    orden = [int(x) for x in args.orden.split(",")] if args.orden else None
    brute_force(td.load_case(args.case), orden=orden)


if __name__ == "__main__":
    main()
