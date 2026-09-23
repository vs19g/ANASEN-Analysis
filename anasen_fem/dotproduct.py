#!/usr/bin/env python3
"""Node-by-node dot product of the physical field with a Ramo weighting field.

Reads the two Elmer VTUs from wires2d.sif and wires2d_weight.sif and writes one
carrying DotProduct = E . E_w. E is in V/m and E_w in 1/m (its electrode is
driven to 1 V), so DotProduct is in V/m^2.

By Shockley-Ramo a charge q drifting at v induces i = -q v . E_w, and for
v = mu E that is -q mu (E . E_w). So this map is the induced-current density up
to the constant -q*mu; --scale folds it in. The signs of q and mu differ
between the electron and ion components, so they are left to the caller.

    python3 dotproduct.py
    python3 dotproduct.py --scale -1.0e-4
"""

import argparse
import sys

import numpy as np

try:
    import meshio
except ImportError:
    sys.exit("meshio is required: pip install meshio")


def find_array(point_data, wanted, label):
    """Look up a point-data array, ignoring case and spacing.

    Elmer lowercases variable names into the VTU and the spelling has varied
    between versions, so match loosely and list what is there on failure.
    """
    def norm(s):
        return "".join(s.lower().split())

    target = norm(wanted)
    for key in point_data:
        if norm(key) == target:
            return point_data[key]
    sys.exit(
        "%s: no point-data array matching %r.\nAvailable arrays: %s"
        % (label, wanted, sorted(point_data))
    )


def as_vectors(arr, label):
    """Return an (N, ncomp) float array, tolerating flat or 2-component data."""
    a = np.asarray(arr, dtype=float)
    if a.ndim == 1:
        sys.exit("%s: expected a vector field, got a scalar array" % label)
    if a.ndim != 2:
        sys.exit("%s: expected a 2D array, got shape %s" % (label, a.shape))
    return a


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--phys", default="wires2d/elfield_anasen_t0001.vtu",
                    help="VTU from the physical solve (wires2d.sif)")
    ap.add_argument("--weight", default="wires2d/elfield_weight_t0001.vtu",
                    help="VTU from the weighting solve (wires2d_weight.sif)")
    ap.add_argument("--out", default="wires2d/dotproduct.vtu",
                    help="output VTU")
    ap.add_argument("--field", default="electric field",
                    help="name of the electric-field array in both inputs")
    ap.add_argument("--scale", type=float, default=None,
                    help="constant (e.g. -q*mu) multiplying the dot product; "
                         "writes an extra ScaledDotProduct array")
    ap.add_argument("--rtol", type=float, default=0.0,
                    help="relative tolerance for the node-position check")
    ap.add_argument("--atol", type=float, default=1e-12,
                    help="absolute tolerance (m) for the node-position check")
    args = ap.parse_args()

    phys = meshio.read(args.phys)
    weight = meshio.read(args.weight)

    # Taken node by node, so both solves must have produced the same node list.
    # If they did not, every value below is meaningless and nothing else checks.
    if phys.points.shape != weight.points.shape:
        sys.exit(
            "Node counts differ: %s has %d nodes, %s has %d. Both solves must "
            "read the same mesh (same wires2d/ directory)."
            % (args.phys, len(phys.points), args.weight, len(weight.points))
        )
    if not np.allclose(phys.points, weight.points, rtol=args.rtol, atol=args.atol):
        worst = np.max(np.linalg.norm(phys.points - weight.points, axis=1))
        sys.exit(
            "Node positions do not match (largest discrepancy %.3e m). The two "
            "solves were run on different meshes." % worst
        )

    E = as_vectors(find_array(phys.point_data, args.field, args.phys), args.phys)
    Ew = as_vectors(find_array(weight.point_data, args.field, args.weight), args.weight)

    # A zero field means the solve never ran: Laplace with no Dirichlet
    # condition applied has the trivial solution and Elmer reports no error.
    # Without this the chain completes and writes zeros, and the only symptom
    # is a blank plot.
    for label, arr, path in (
        ("physical", E, args.phys),
        ("weighting", Ew, args.weight),
    ):
        if not np.any(arr):
            sys.exit(
                "the %s field in %s is identically zero.\n"
                "Elmer returned the trivial solution, which means no boundary "
                "condition was applied. Compare the Target Bodies and Target "
                "Boundaries numbers in the .sif against wires2d/mesh.names -- "
                "ElmerGrid's -autoclean renumbers them." % (label, path)
            )

    # Elmer may pad 2D vectors to three components; compare only the
    # components both arrays actually carry.
    ncomp = min(E.shape[1], Ew.shape[1])
    if E.shape[1] != Ew.shape[1]:
        print("note: component counts differ (%d vs %d), using the first %d"
              % (E.shape[1], Ew.shape[1], ncomp))
    dot = np.einsum("ij,ij->i", E[:, :ncomp], Ew[:, :ncomp])

    out = meshio.Mesh(
        points=phys.points,
        cells=phys.cells,
        point_data=dict(phys.point_data),
        cell_data=dict(phys.cell_data),
    )
    out.point_data["DotProduct"] = dot
    out.point_data["WeightingField"] = Ew
    if args.scale is not None:
        out.point_data["ScaledDotProduct"] = args.scale * dot

    meshio.write(args.out, out)

    finite = dot[np.isfinite(dot)]
    print("nodes            : %d" % len(dot))
    print("components used  : %d" % ncomp)
    print("DotProduct min   : %.6e V/m^2" % finite.min())
    print("DotProduct max   : %.6e V/m^2" % finite.max())
    print("DotProduct mean  : %.6e V/m^2" % finite.mean())
    if len(finite) != len(dot):
        print("WARNING: %d non-finite nodes" % (len(dot) - len(finite)))
    print("wrote %s" % args.out)


if __name__ == "__main__":
    main()
