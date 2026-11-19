#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gera as figuras e a tabela consolidada do artigo a partir dos CSVs exportados
do Jupiter Analytics (nenhum número é inventado: tudo vem dos arquivos).

Entradas
  data/e1e2/*.csv        varredura do E1+E2 (Scenarios), um CSV por metrica,
                         cujo nome contem a chave da metrica (ex.: DP_pct).
  data/manual_runs.csv   rodadas manuais de E3 e E4 (modelo em data/).

Saidas (pasta atual)
  fig_e1e2.pdf, fig_e2.pdf, fig_rho.pdf, fig_e3e4.pdf,
  tabela_consolidada.tex, tabela_validacao.tex,
  resultados_consolidados.csv, validacao_mmnk.csv

Uso:  python3 make_figures.py [--x-mode auto|rate|delay]
"""
import argparse, csv, glob, os, re, sys
from math import factorial
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

METRICS = ["TP_req_s", "U_pct", "DP_pct", "MRT_ms", "MQT_ms", "MPT_ms"]
S_MS = 50.0                       # tempo medio de servico (ms)
MARK = ["o", "s", "^", "D"]
LSTY = ["-", "--", ":", "-."]
COLR = ["#1f4e9c", "#c0392b", "#2e8b57", "#7d3c98"]

plt.rcParams.update({"font.size": 7, "axes.labelsize": 7, "legend.fontsize": 6,
                     "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                     "axes.titlesize": 7, "lines.linewidth": 1.1,
                     "lines.markersize": 3.2, "axes.grid": True,
                     "grid.alpha": 0.3, "pdf.fonttype": 42})


# ----------------------------------------------------------------- leitura
def num(s):
    s = s.strip().replace(",", ".")
    return float(s) if s not in ("", "NaN", "nan") else float("nan")


def read_sweep(path, mode):
    with open(path, encoding="utf-8-sig", newline="") as f:
        lines = [l for l in f.read().splitlines() if l.strip()]
    delim = ";" if lines[0].count(";") >= lines[0].count(",") else ","
    rows = list(csv.reader(lines, delimiter=delim))
    header, body = rows[0], rows[1:]
    xs = [num(r[0]) for r in body]
    m = mode
    if mode == "auto":
        m = "rate" if ("1/" in header[0] or max(xs) <= 1.0) else "delay"
    lam = [x * 1000.0 if m == "rate" else 1000.0 / x for x in xs]   # req/s
    series = {}
    for j, h in enumerate(header[1:], start=1):
        mm = re.search(r"=\s*([0-9]+)", h)
        if not mm:
            continue
        series[int(mm.group(1))] = [num(r[j]) for r in body]
    return lam, series, m


def load_sweep(folder, mode):
    data, lam_ref = {}, None
    for met in METRICS:
        files = sorted(glob.glob(os.path.join(folder, f"*{met}*.csv")),
                       key=os.path.getmtime)
        if not files:
            print(f"  [aviso] sem CSV para {met} em {folder}")
            continue
        if len(files) > 1:
            print(f"  [aviso] varios CSVs para {met}; usando o mais recente: "
                  f"{os.path.basename(files[-1])}")
        lam, series, m = read_sweep(files[-1], mode)
        data[met] = series
        if lam_ref is None:
            lam_ref = lam
            print(f"  eixo X interpretado como '{m}' -> lambda (req/s) = "
                  f"{[round(v, 2) for v in lam]}")
        elif any(abs(a - b) > 0.5 for a, b in zip(lam, lam_ref)):
            print(f"  [ERRO] eixo X de {met} difere dos demais: confira a varredura")
    return lam_ref, data


def at(data, met, N, lam_target, lam):
    if met not in data or N not in data[met]:
        return None
    for l, v in zip(lam, data[met][N]):
        if abs(l - lam_target) < 0.5:
            return v
    return None


# --------------------------------------------------------------- analitico
def mmnk(lam, N, K=20, S=S_MS):
    mu = 1000.0 / S
    a = lam / mu
    w = [a ** n / factorial(n) if n <= N else
         a ** N / factorial(N) * (lam / (N * mu)) ** (n - N) for n in range(K + 1)]
    Z = sum(w)
    p = [x / Z for x in w]
    Lq = sum((n - N) * p[n] for n in range(N, K + 1))
    Ls = sum(min(n, N) * p[n] for n in range(K + 1))
    le = lam * (1 - p[K])
    return {"DP_pct": 100 * p[K], "TP_req_s": le, "U_pct": 100 * Ls / N,
            "MRT_ms": (Lq + Ls) / le * 1000, "MQT_ms": Lq / le * 1000,
            "MPT_ms": S}


def validate(lam, data, K=20):
    out = []
    for met in METRICS:
        for N, vals in data.get(met, {}).items():
            for l, v in zip(lam, vals):
                t = mmnk(l, N, K)[met]
                out.append([met, N, round(l, 3), v, round(t, 4), v - t,
                            100 * (v - t) / t if abs(t) > 1e-9 else float("nan")])
    with open("validacao_mmnk.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["metrica", "N", "lambda_req_s", "jupiter", "analitico",
                    "dif_abs", "dif_rel_pct"])
        w.writerows(out)
    # balanco de fluxo: TP deve ser lambda*(1-DP): detecta simulacao mal configurada
    if "TP_req_s" in data and "DP_pct" in data:
        worst = 0.0
        for N in data["TP_req_s"]:
            for l, tp, dp in zip(lam, data["TP_req_s"][N], data["DP_pct"].get(N, [])):
                worst = max(worst, abs(tp - l * (1 - dp / 100.0)))
        print(f"\nBalanco de fluxo TP = lambda*(1-DP): max |dif| = {worst:.3f} req/s")
    print("\nValidacao Jupiter x M/M/N/K (K=20):")
    for met in METRICS:
        rel = [abs(r[6]) for r in out if r[0] == met and r[6] == r[6] and abs(r[4]) > 1]
        ab = [abs(r[5]) for r in out if r[0] == met]
        if ab:
            print(f"  {met:9s} max |dif abs| = {max(ab):8.3f}"
                  + (f" | max |dif rel| = {max(rel):5.2f}%" if rel else ""))


# ------------------------------------------------------------------ figuras
def fig_e1e2(lam, data):
    panels = [("U_pct", "Utilização (%)", "(a)"),
              ("TP_req_s", "Vazão (req/s)", "(b)"),
              ("DP_pct", "Prob. de descarte (%)", "(c)"),
              ("MRT_ms", "Tempo de resposta (ms)", "(d)"),
              ("MQT_ms", "Tempo de fila (ms)", "(e)"),
              ("MPT_ms", "Tempo de processamento (ms)", "(f)")]
    panels = [p for p in panels if p[0] in data]
    if not panels:
        return
    fig, axs = plt.subplots(2, 3, figsize=(7.16, 3.9))
    axs = list(axs.ravel())
    for ax in axs[len(panels):]:
        ax.axis("off")
    for ax, (met, yl, tag) in zip(axs, panels):
        for i, N in enumerate(sorted(data[met])):
            ax.plot(lam, data[met][N], marker=MARK[i % 4], ls=LSTY[i % 4],
                    color=COLR[i % 4], label=f"N = {N}")
        ax.set_xlabel("Taxa de chegada $\\lambda$ (req/s)")
        ax.set_ylabel(yl)
        ax.set_xticks([round(v) for v in lam])
        ax.set_title(tag, loc="left", pad=2)
        if met == "MPT_ms":
            ax.set_ylim(40, 60)
    axs[0].legend(loc="best", frameon=True)
    fig.tight_layout(pad=0.4)
    fig.savefig("fig_e1e2.pdf")
    plt.close(fig)
    print("  gerada fig_e1e2.pdf")


def bars(ax, labels, vals, ylabel, title, fmt="{:.1f}"):
    xs = range(len(labels))
    ax.bar(xs, vals, color=["#7f8fa6", "#1f4e9c", "#a0b4d6"][:len(vals)],
           edgecolor="black", linewidth=0.5, width=0.62)
    for x, v in zip(xs, vals):
        ax.text(x, v, fmt.format(v), ha="center", va="bottom", fontsize=5.5)
    ax.set_xticks(list(xs)); ax.set_xticklabels(labels)
    ax.set_ylabel(ylabel); ax.set_title(title, loc="left", pad=2)
    ax.set_ylim(0, max(vals) * 1.18 if max(vals) > 0 else 1)


def fig_e3e4(rows):
    e3 = [r for r in rows if r["exp"] == "E3"]
    e4 = [r for r in rows if r["exp"] == "E4"]
    ok3 = len(e3) == 3 and all(r["DP_pct"] is not None and r["MRT_ms"] is not None for r in e3)
    ok4 = len(e4) == 3 and all(None not in (r["U_pct"], r["MRT_ms"], r["DP_pct"]) for r in e4)
    if not (ok3 or ok4):
        print("  [aviso] E3/E4 incompletos: fig_e3e4.pdf nao gerada")
        return
    n = (2 if ok3 else 0) + (3 if ok4 else 0)
    fig, axs = plt.subplots(1, n, figsize=(7.16, 2.2))
    axs = [axs] if n == 1 else list(axs)
    k = 0
    if ok3:
        e3.sort(key=lambda r: r["K"])
        lab = [f"K={r['K']}" for r in e3]
        bars(axs[k], lab, [r["DP_pct"] for r in e3], "DP (%)", "(a) E3: descarte", "{:.2f}"); k += 1
        bars(axs[k], lab, [r["MRT_ms"] for r in e3], "MRT (ms)", "(b) E3: tempo de resposta"); k += 1
    if ok4:
        order = {"Baixo": 0, "Médio": 1, "Medio": 1, "Alto": 2}
        e4.sort(key=lambda r: order.get(r["label"], 9))
        lab = [r["label"] for r in e4]
        bars(axs[k], lab, [r["U_pct"] for r in e4], "U (%)", "(c) E4: utilização"); k += 1
        bars(axs[k], lab, [r["MRT_ms"] for r in e4], "MRT (ms)", "(d) E4: tempo de resposta"); k += 1
        bars(axs[k], lab, [r["DP_pct"] for r in e4], "DP (%)", "(e) E4: descarte", "{:.2f}"); k += 1
    fig.tight_layout(pad=0.4)
    fig.savefig("fig_e3e4.pdf")
    plt.close(fig)
    print("  gerada fig_e3e4.pdf")



def fig_e2(lam, data):
    """E2: lambda=60 req/s, K=20; barras por N."""
    mets = [("U_pct", "U (%)", "(a) E2: utilização", "{:.1f}"),
            ("MRT_ms", "MRT (ms)", "(b) E2: tempo de resposta", "{:.1f}"),
            ("DP_pct", "DP (%)", "(c) E2: descarte", "{:.2f}")]
    Ns = sorted(data.get("MRT_ms", {}))
    vals = {m[0]: [at(data, m[0], N, 60, lam) for N in Ns] for m in mets}
    if not Ns or any(v is None for m in mets for v in vals[m[0]]):
        print("  [aviso] E2 incompleto (precisa de lambda=60 e U, MRT, DP): fig_e2.pdf nao gerada")
        return
    fig, axs = plt.subplots(1, 3, figsize=(7.16, 2.1))
    for ax, (met, yl, tt, fm) in zip(axs, mets):
        bars(ax, [f"N={N}" for N in Ns], vals[met], yl, tt, fm)
    fig.tight_layout(pad=0.4)
    fig.savefig("fig_e2.pdf")
    plt.close(fig)
    print("  gerada fig_e2.pdf")


def fig_rho(lam, data):
    """MRT e DP em funcao de rho = lambda/(N*mu): efeito de escala."""
    if "MRT_ms" not in data or "DP_pct" not in data:
        return
    mu = 1000.0 / S_MS
    fig, axs = plt.subplots(1, 2, figsize=(3.5, 1.9))
    for ax, (met, yl, tt) in zip(axs, [("MRT_ms", "MRT (ms)", "(a)"),
                                       ("DP_pct", "DP (%)", "(b)")]):
        for i, N in enumerate(sorted(data[met])):
            ax.plot([l / (N * mu) for l in lam], data[met][N], marker=MARK[i % 4],
                    ls=LSTY[i % 4], color=COLR[i % 4], label=f"N = {N}")
        ax.set_xlabel("Intensidade de tráfego $\\rho$")
        ax.set_ylabel(yl)
        ax.set_title(tt, loc="left", pad=2)
    axs[0].legend(loc="best")
    fig.tight_layout(pad=0.4)
    fig.savefig("fig_rho.pdf")
    plt.close(fig)
    print("  gerada fig_rho.pdf")


def write_validation_table(lam, data):
    """Tabela Jupiter x analitico no cenario de referencia (lambda=60, N=4, K=20)."""
    if not lam:
        return
    names = [("TP_req_s", "TP (req/s)", 2), ("U_pct", "U (\\%)", 2), ("DP_pct", "DP (\\%)", 3),
             ("MRT_ms", "MRT (ms)", 2), ("MQT_ms", "MQT (ms)", 2), ("MPT_ms", "MPT (ms)", 2)]
    ref = mmnk(60, 4, 20)
    L = ["\\begin{tabular}{@{}lrrr@{}}", "\\toprule",
         "Métrica & Jupiter & Analítico & Dif. (\\%) \\\\", "\\midrule"]
    n = 0
    for met, nm, d in names:
        v = at(data, met, 4, 60, lam)
        if v is None:
            continue
        t = ref[met]
        rel = 100 * (v - t) / t if abs(t) > 1e-9 else float("nan")
        L.append(f"{nm} & {v:.{d}f} & {t:.{d}f} & {rel:+.1f} \\\\")
        n += 1
    L += ["\\bottomrule", "\\end{tabular}"]
    if n:
        open("tabela_validacao.tex", "w", encoding="utf-8").write("\n".join(L) + "\n")
        print("  gerada tabela_validacao.tex")


# ------------------------------------------------- manual + tabela consolidada
def read_manual(path, lam, data):
    rows = []
    if not os.path.exists(path):
        print(f"  [aviso] {path} nao encontrado")
        return rows
    with open(path, encoding="utf-8-sig", newline="") as f:
        for r in csv.DictReader(f):
            d = {"exp": r["exp"].strip(), "label": r["label"].strip(),
                 "lam": num(r["lambda_req_s"]), "N": int(r["N"]), "K": int(r["K"])}
            for met in METRICS:
                v = (r.get(met) or "").strip()
                d[met] = num(v) if v else None
            # fallback: ponto de referencia vem da varredura (lambda=60, N=4, K=20)
            if d["N"] == 4 and d["K"] == 20 and abs(d["lam"] - 60) < 0.5 and lam:
                for met in METRICS:
                    if d[met] is None:
                        d[met] = at(data, met, 4, 60, lam)
            rows.append(d)
    return rows


def consolidate(lam, data, manual):
    rows = []
    if lam:
        for l in lam:                                    # E1: N=4, K=20
            rows.append({"exp": "E1", "label": f"$\\lambda$={l:g}", "lam": l, "N": 4, "K": 20,
                         **{m: at(data, m, 4, l, lam) for m in METRICS}})
        for N in (2, 8):                                 # E2: lambda=60, K=20
            rows.append({"exp": "E2", "label": f"N={N}", "lam": 60, "N": N, "K": 20,
                         **{m: at(data, m, N, 60, lam) for m in METRICS}})
    for r in manual:
        if r["exp"] == "E3" and r["K"] != 20:
            rows.append(r)
        if r["exp"] == "E4" and r["label"] in ("Baixo", "Alto"):
            rows.append(r)
    return [r for r in rows if any(r[m] is not None for m in METRICS)]


def write_table(rows):
    def f(v, d):
        return "--" if v is None or v != v else f"{v:.{d}f}"
    with open("resultados_consolidados.csv", "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["exp", "label", "lambda_req_s", "N", "K"] + METRICS)
        for r in rows:
            w.writerow([r["exp"], r["label"], r["lam"], r["N"], r["K"]] +
                       [r[m] for m in METRICS])
    L = ["\\begin{tabular}{lrrrrrrrrr}", "\\toprule",
         "Exp. & $\\lambda$ & $N$ & $K$ & TP & U & DP & MRT & MQT & MPT \\\\",
         " & (req/s) & & & (req/s) & (\\%) & (\\%) & (ms) & (ms) & (ms) \\\\",
         "\\midrule"]
    prev = None
    for r in rows:
        if prev and r["exp"] != prev:
            L.append("\\midrule")
        prev = r["exp"]
        L.append(f"{r['exp']} & {r['lam']:g} & {r['N']} & {r['K']} & "
                 f"{f(r['TP_req_s'],1)} & {f(r['U_pct'],1)} & {f(r['DP_pct'],2)} & "
                 f"{f(r['MRT_ms'],1)} & {f(r['MQT_ms'],1)} & {f(r['MPT_ms'],1)} \\\\")
    L += ["\\bottomrule", "\\end{tabular}"]
    with open("tabela_consolidada.tex", "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"  gerada tabela_consolidada.tex ({len(rows)} linhas)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--x-mode", default="auto", choices=["auto", "rate", "delay"])
    ap.add_argument("--sweep-dir", default="data/e1e2")
    ap.add_argument("--manual", default="data/manual_runs.csv")
    a = ap.parse_args()
    print("Lendo varredura E1/E2...")
    lam, data = load_sweep(a.sweep_dir, a.x_mode)
    if lam:
        fig_e1e2(lam, data)
        fig_e2(lam, data)
        fig_rho(lam, data)
        write_validation_table(lam, data)
        validate(lam, data)
    else:
        print("  [aviso] nenhum CSV de varredura encontrado")
    print("Lendo rodadas manuais E3/E4...")
    manual = read_manual(a.manual, lam, data)
    fig_e3e4(manual)
    write_table(consolidate(lam, data, manual))


if __name__ == "__main__":
    sys.exit(main())
