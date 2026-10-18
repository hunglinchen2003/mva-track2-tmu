"""Save the six lead Vina poses and draw them against Lys1002, as in the earlier pocket figure."""
from __future__ import annotations

import json
import re
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

VINA = Path(r"F:\BUB1B篩藥\analysis\vina")
OUT = Path(r"F:\2026_MVA_tract2_submission\analysis\fda_top6_k1002")
DRUGS = [
    ("Pralsetinib", "L00296.pdbqt"),
    ("Lapatinib", "L01167.pdbqt"),
    ("Tucatinib", "L00873.pdbqt"),
    ("Entrectinib", "L00431.pdbqt"),
    ("Capmatinib", "L00508.pdbqt"),
    ("Ponatinib", "L00231.pdbqt"),
]
AA = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
STYLES = {
    "hbond": ("#ca8a04", (0, (3, 2)), 1.6),
    "salt": ("#c2410c", (0, (1, 1.5)), 1.5),
    "pi": ("#15803d", (0, (4, 1.6)), 1.3),
    "halogen": ("#0f766e", (0, (2, 2)), 1.3),
    "pi_cation": ("#7c3aed", (0, (2, 2)), 1.3),
    "hydrophobic": ("#78716c", (0, (1, 2)), 1.1),
}
COLORS = ["#1b6ca8", "#c2410c", "#15803d", "#7c3aed", "#b45309", "#be123c"]


def atoms(path: Path, ligand=False):
    rows = []
    started = False
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith("MODEL"):
            started = True
            continue
        if line.startswith("ENDMDL"):
            break
        if ligand and not started:
            continue
        if not line.startswith(("ATOM", "HETATM")):
            continue
        try:
            rows.append({
                "name": line[12:16].strip(),
                "resn": line[17:20].strip() or "UNL",
                "resi": int(line[22:26]) if line[22:26].strip() else 1,
                "xyz": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])]),
                "elem": (line[76:78].strip() or re.sub(r"[^A-Za-z]", "", line[12:14]) or "C")[:1].upper(),
            })
        except ValueError:
            continue
    return rows


def vina_score(path: Path):
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith("REMARK VINA RESULT:"):
            return float(line.split()[3])
    return None


def write_complex(receptor, ligand, dest: Path):
    lines = []
    serial = 1
    for atom in receptor:
        if atom["elem"] == "H":
            continue
        name = f"{atom['name']:<4}"[:4]
        lines.append(
            f"ATOM  {serial:5d} {name} {atom['resn']:>3} A{atom['resi']:4d}    "
            f"{atom['xyz'][0]:8.3f}{atom['xyz'][1]:8.3f}{atom['xyz'][2]:8.3f}  1.00  0.00          {atom['elem']:>2}"
        )
        serial += 1
    for atom in ligand:
        if atom["elem"] == "H":
            continue
        name = f"{atom['name']:<4}"[:4]
        lines.append(
            f"HETATM{serial:5d} {name} UNL L   1    "
            f"{atom['xyz'][0]:8.3f}{atom['xyz'][1]:8.3f}{atom['xyz'][2]:8.3f}  1.00  0.00          {atom['elem']:>2}"
        )
        serial += 1
    dest.write_text("\n".join(lines) + "\nEND\n", encoding="utf-8")


def parse_plip(path: Path):
    text = path.read_text(errors="ignore")
    chunks = re.split(r"\*\*([^*]+)\*\*", text)
    wanted = {
        "Hydrogen Bonds": "hbond",
        "Salt Bridges": "salt",
        "pi-Stacking": "pi",
        "Halogen Bonds": "halogen",
        "pi-Cation Interactions": "pi_cation",
        "Hydrophobic Interactions": "hydrophobic",
    }
    hits = []
    for i in range(1, len(chunks), 2):
        kind = wanted.get(chunks[i].strip())
        if kind is None:
            continue
        for line in chunks[i + 1].splitlines():
            if not line.startswith("|") or "RESNR" in line or set(line.strip("| ")) <= set("-="):
                continue
            coords = re.findall(r"(-?\d+\.\d+),\s*(-?\d+\.\d+),\s*(-?\d+\.\d+)", line)
            head = re.match(r"\|\s+(\d+)\s+\|\s+([A-Z]{3})\s+\|", line)
            dist_m = re.search(r"\|\s+([0-9]+\.[0-9]+)\s+\|", line)
            if head is None or len(coords) < 2:
                continue
            distance = float(dist_m.group(1)) if dist_m else None
            if kind != "hbond" and distance is not None and distance < 2.0:
                continue
            hits.append({
                "kind": kind,
                "resi": int(head.group(1)),
                "resn": head.group(2),
                "dist": distance,
                "lig": np.array([float(v) for v in coords[-2]]),
                "prot": np.array([float(v) for v in coords[-1]]),
            })
    return hits


def bonds(xyz, lo=0.95, hi=1.75):
    pairs = []
    if len(xyz) < 2:
        return pairs
    d = np.linalg.norm(xyz[:, None, :] - xyz[None, :, :], axis=2)
    for i in range(len(xyz)):
        for j in range(i + 1, len(xyz)):
            if lo <= d[i, j] <= hi:
                pairs.append((i, j))
    return pairs


def project(xyz, origin, x_axis, y_axis):
    rel = np.atleast_2d(xyz) - origin
    return np.column_stack([rel @ x_axis, rel @ y_axis])


def main():
    receptor = [a for a in atoms(VINA / "receptor_mut.pdbqt") if a["resn"] != "UNL"]
    lys = [a for a in receptor if a["resi"] == 302 and a["name"] in {"CA", "CB", "CG", "CD", "CE", "NZ"}]
    nz = next(a["xyz"] for a in receptor if a["resi"] == 302 and a["name"] == "NZ")
    (OUT / "poses").mkdir(parents=True, exist_ok=True)
    (OUT / "complexes").mkdir(parents=True, exist_ok=True)
    comps = []
    for name, pdbqt in DRUGS:
        pose = VINA / "out_mut" / pdbqt
        lig = atoms(pose, ligand=True)
        if not lig:
            raise SystemExit(f"missing pose {pose}")
        pose_copy = OUT / "poses" / f"{name}_mut.pdbqt"
        pose_copy.write_text(pose.read_text(errors="ignore"), encoding="utf-8")
        complex_path = OUT / "complexes" / f"{name}_mut_complex.pdb"
        write_complex(receptor, lig, complex_path)
        plip_dir = OUT / "plip" / name
        plip_dir.mkdir(parents=True, exist_ok=True)
        subprocess.run(
            ["plip", "-f", str(complex_path), "-o", str(plip_dir), "-x", "-t"],
            check=False, capture_output=True, text=True, encoding="utf-8", errors="replace",
        )
        reports = list(plip_dir.glob("*_report.txt"))
        hits = parse_plip(reports[0]) if reports else []
        heavy = [a for a in lig if a["elem"] != "H"]
        d_nz = float(np.min(np.linalg.norm(np.stack([a["xyz"] for a in heavy]) - nz, axis=1)))
        k_hits = [h for h in hits if h["resi"] == 302]
        comps.append({
            "name": name, "kcal": vina_score(pose), "lig": heavy, "plip": hits,
            "d_nz": d_nz, "k1002": k_hits,
        })
        print(name, vina_score(pose), "dNZ", round(d_nz, 2), "K1002", len(k_hits), {h["kind"] for h in hits}, flush=True)

    summary = []
    for comp in comps:
        summary.append({
            "name": comp["name"],
            "vina_kcal": comp["kcal"],
            "min_heavy_to_nz_A": round(comp["d_nz"], 2),
            "contacts": [
                {
                    "kind": h["kind"],
                    "uniprot": h["resi"] + 700,
                    "resn": h["resn"],
                    "dist_A": h["dist"],
                    "to_K1002": h["resi"] == 302,
                }
                for h in comp["plip"]
            ],
        })
    (OUT / "interactions.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    lig_all = np.vstack([np.stack([a["xyz"] for a in c["lig"]]) for c in comps])
    kxyz = np.stack([a["xyz"] for a in lys])
    cloud = np.vstack([lig_all, kxyz])
    origin = cloud.mean(0)
    _, _, Vt = np.linalg.svd(cloud - origin, full_matrices=False)
    x_a, y_a = Vt[0], Vt[1]
    ca = np.stack([a["xyz"] for a in receptor if a["name"] == "CA" and np.linalg.norm(a["xyz"] - nz) < 16])

    fig = plt.figure(figsize=(15.6, 12.4), dpi=170)
    fig.patch.set_facecolor("white")
    overview = fig.add_axes((0.04, 0.70, 0.72, 0.26))
    overview.set_title("A    Six lead AutoDock Vina poses around K1002", loc="left", fontsize=13, fontweight="bold")
    overview.set_aspect("equal")
    overview.axis("off")
    trace = project(ca, origin, x_a, y_a)
    overview.plot(trace[:, 0], trace[:, 1], color="#d5cfc6", lw=2.2, solid_capstyle="round", zorder=1)
    kp = project(kxyz, origin, x_a, y_a)
    for i, j in bonds(kxyz, 1.2, 1.7):
        overview.plot([kp[i, 0], kp[j, 0]], [kp[i, 1], kp[j, 1]], color="#9a3412", lw=2.8, zorder=5)
    overview.scatter(kp[:, 0], kp[:, 1], s=18, c="#9a3412", zorder=6)
    nz_p = project(nz.reshape(1, 3), origin, x_a, y_a)[0]
    overview.annotate("K1002", xy=nz_p, xytext=(6, 8), textcoords="offset points", fontsize=11, fontweight="bold", color="#9a3412")
    for comp, color in zip(comps, COLORS):
        xyz = np.stack([a["xyz"] for a in comp["lig"]])
        p = project(xyz, origin, x_a, y_a)
        for i, j in bonds(xyz):
            overview.plot([p[i, 0], p[j, 0]], [p[i, 1], p[j, 1]], color=color, lw=1.45, zorder=3)
    site = project(cloud, origin, x_a, y_a)
    pad = 3.5
    overview.set_xlim(site[:, 0].min() - pad, site[:, 0].max() + pad)
    overview.set_ylim(site[:, 1].min() - pad, site[:, 1].max() + pad)

    fig.legend(
        handles=[Line2D([0], [0], color=c, lw=2.4, label=comp["name"]) for comp, c in zip(comps, COLORS)],
        loc="upper right", frameon=False, fontsize=9, bbox_to_anchor=(0.98, 0.94),
    )
    for index, (comp, color) in enumerate(zip(comps, COLORS)):
        row_i, col = divmod(index, 3)
        ax = fig.add_axes((0.02 + col * 0.33, 0.36 - row_i * 0.32, 0.31, 0.27))
        ax.set_aspect("equal")
        ax.axis("off")
        xyz = np.stack([a["xyz"] for a in comp["lig"]])
        contact_ids = {302} | {hit["resi"] for hit in comp["plip"]}
        chains = []
        for resi in sorted(contact_ids):
            chain = [a for a in receptor if a["resi"] == resi and a["elem"] != "H" and a["name"] not in {"N", "C", "O"}]
            if len(chain) >= 2:
                chains.append(chain)
        scene = [xyz] + [np.stack([a["xyz"] for a in chain]) for chain in chains]
        cloud_b = np.vstack(scene)
        center = cloud_b.mean(0)
        _, _, Vt = np.linalg.svd(cloud_b - center, full_matrices=False)
        x_b, y_b = Vt[0], Vt[1]
        p = project(xyz, center, x_b, y_b)
        for i, j in bonds(xyz):
            ax.plot([p[i, 0], p[j, 0]], [p[i, 1], p[j, 1]], color="#4b5563", lw=1.7, zorder=3)
        elem_c = [
            "#1d4e89" if a["elem"] == "N" else "#b91c1c" if a["elem"] == "O" else "#0f766e" if a["elem"] in {"F", "S", "Cl", "P"} else "#6b7280"
            for a in comp["lig"]
        ]
        ax.scatter(p[:, 0], p[:, 1], s=16, c=elem_c, zorder=4, linewidths=0)
        labeled = {}
        for chain in chains:
            cxyz = np.stack([a["xyz"] for a in chain])
            cp = project(cxyz, center, x_b, y_b)
            tone = "#9a3412" if chain[0]["resi"] == 302 else "#57534e"
            for i, j in bonds(cxyz, 1.15, 1.7):
                ax.plot([cp[i, 0], cp[j, 0]], [cp[i, 1], cp[j, 1]], color=tone, lw=2.1, zorder=4)
            label = f"{AA.get(chain[0]['resn'], chain[0]['resn'])}{chain[0]['resi'] + 700}"
            labeled[label] = cp.mean(0)
        for hit in comp["plip"]:
            a = project(hit["lig"].reshape(1, 3), center, x_b, y_b)[0]
            b = project(hit["prot"].reshape(1, 3), center, x_b, y_b)[0]
            colr, ls, lw = STYLES.get(hit["kind"], STYLES["hbond"])
            ax.plot([a[0], b[0]], [a[1], b[1]], color=colr, lw=lw, ls=ls, zorder=2)
        center2 = p.mean(0)
        for label, point in labeled.items():
            vec = point - center2
            norm = float(np.linalg.norm(vec))
            shift = np.array([0.0, 2.2]) if norm < 1e-6 else 2.6 * vec / norm
            ax.annotate(
                label, xy=point, xytext=(point[0] + shift[0], point[1] + shift[1]),
                fontsize=8, fontweight="bold", color="#9a3412" if label == "K1002" else "#1c1915",
                ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.92),
                zorder=6,
            )
        shown = np.vstack([p] + [project(np.stack([a["xyz"] for a in chain]), center, x_b, y_b) for chain in chains])
        pad = 3.4
        ax.set_xlim(shown[:, 0].min() - pad, shown[:, 0].max() + pad)
        ax.set_ylim(shown[:, 1].min() - pad, shown[:, 1].max() + pad)
        k_note = "contact" if comp["k1002"] else "no direct contact"
        ax.text(
            0.02, 0.98,
            f"{comp['name']}   {comp['kcal']:.2f} kcal/mol\n{comp['d_nz']:.1f} Å to NZ, {k_note}",
            transform=ax.transAxes, fontsize=8.5, color=color, fontweight="bold", va="top", ha="left",
        )

    fig.text(
        0.02, 0.655,
        "B    PLIP contacts on each Vina pose. Yellow = hydrogen bond, gray = hydrophobic, green = pi-stack, purple = pi-cation. Residue numbers are UniProt.",
        fontsize=11, fontweight="bold",
    )
    fig.legend(handles=[
        Line2D([0], [0], color="#ca8a04", lw=1.5, ls=(0, (3, 2)), label="H-bond"),
        Line2D([0], [0], color="#78716c", lw=1.2, ls=(0, (1, 2)), label="hydrophobic"),
        Line2D([0], [0], color="#15803d", lw=1.3, ls=(0, (4, 1.6)), label="pi-stack"),
        Line2D([0], [0], color="#7c3aed", lw=1.3, ls=(0, (2, 2)), label="pi-cation"),
        Line2D([0], [0], color="#c2410c", lw=1.4, ls=(0, (1, 1.5)), label="salt bridge"),
        Line2D([0], [0], color="#9a3412", lw=2.4, label="K1002"),
    ], loc="lower left", ncol=6, frameon=False, fontsize=9, bbox_to_anchor=(0.02, 0.0))
    png = OUT / "k1002_vina_top6.png"
    pdf = OUT / "k1002_vina_top6.pdf"
    fig.savefig(png, dpi=200, bbox_inches="tight", facecolor="white")
    fig.savefig(pdf, bbox_inches="tight", facecolor="white")
    print("wrote", png)


if __name__ == "__main__":
    main()
