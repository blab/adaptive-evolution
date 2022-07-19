from __future__ import print_function
from pymol import cmd


def saveSurfaceResidues(output_filename):

    myspace = {'surfaceresidues': []}
    cmd.iterate("(exposed_res_01) and n. CA", 'surfaceresidues.append([chain, resi])', space=myspace)


    with open(f"/Users/katekistler/nextstrain/adaptive-evolution/adaptive_loci_results/fixations_per_site/results/{output_filename}.txt", "w") as f:
        for x in myspace['surfaceresidues']:
            f.write(str(x))
            f.write('\n')

    f.close()


cmd.extend("saveSurfaceResidues", saveSurfaceResidues)
