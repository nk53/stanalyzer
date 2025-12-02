import argparse
import typing as t

import stanalyzer.cli.stanalyzer as sta
import MDAnalysis as mda
import math

ANALYSIS_NAME = 'compressibility_modulus'

T = t.TypeVar('T')
Array1D: t.TypeAlias = list[T]
Array2D: t.TypeAlias = list[Array1D[T]]
FL: t.TypeAlias = list[float]
Tup5: t.TypeAlias = tuple[FL, FL, FL, FL, float]


def header(outfile: sta.FileLike | None = None) -> str:
    """Returns a header string and, if optionally writes it to a file"""
    header_str = "compressibility_modulus"

    print(header_str, file=outfile)

    return header_str


def write_compressibility_modulus(psf: sta.FileRef, traj: sta.FileRefList, out: sta.FileRef,
                                  cntmin: int, cntmax: int, bs: int, temp: float,
                                  debug: bool = False) -> None:
    """Writes Compressibility Modulus to 'out' file"""
    # Inputs:
    # cntmin: First dcd # to start the analysis
    # cntmax: Last dcd # to start the analysis
    # cntmax: L
    u = mda.Universe(psf, traj[cntmin:cntmax+1])

    # fout=sys.argv[2]
    # Need to estimate KA from fluctuation of Box ...
    kB = 1.380649e-23  # J/K
    T = 310.15         # Temperature of the system
    kBT = kB*T         # J=N*m = 10^5 dyn * 10^2 cm
    kBT = kBT*1.0e7    # in dyne

    # -----------------------------------------------
    # Unit conversion factors
    # -----------------------------------------------
    atm2bar = 1.01325   # bar/atm
    bar2dyn_cm2 = 1.e6  # dyn/cm^2 / bar
    ang2cm = 1.e-8      # cm/A

    if debug:
        # surface tension (atm*A) => (dyn/cm)
        #       unit is [atm*A] = atm*(atm2bar * bar2dyn_cm2) * A *(ang2cm)
        #                       = (atm2bar*bar2dyn_cm2*ang2cm) [dyn/cm^2 * cm]
        fact_gamma = atm2bar*bar2dyn_cm2*ang2cm

        # free energy derivative w.r.t. curvature
        navog = 6.0221429e23  # /mol
        dyn2new = 1.e-5        # N/dyn
        ang2met = 1.e-10       # m/A
        joule2cal = 1.0/4.184    # cal/joule
        # definition Energy/Area/Curvature -> J/A^2/A^-1 = J/A (SI unit)
        # unit is [atm*A*A] = fact_gamma*ang2cm [dyn]
        #         [atm*A*A*A/A] = fact_gamma*ang2cm [dyn*A/A]
        #                       = fact_gamma*ang2cm*dyn2new*ang2met [N*m/A]
        #                       = fact_gamma*ang2cm*dyn2new*ang2met [J/A]
        # desired output unit is kcal/mol/A
        #       so multiply navog * 10-3 (kcal/cal) * joule2cal to the above unit conversion
        fact_dfdr = fact_gamma*ang2cm*dyn2new*ang2met*navog*joule2cal*1.e-3  # kcal/mol/A
        print('fact_dfdr:', fact_dfdr)

    def read_data(u: mda.Universe, j: int) -> Tup5:
        boxx: FL = []
        boxy: FL = []
        boxz: FL = []
        sa: FL = []
        norm = 0.0
        for i, ts in enumerate(u.trajectory[j*bs:j*bs+bs]):
            tboxx = ts.dimensions[0]  # boxx = boxy
            tboxy = ts.dimensions[1]  # boxy
            tboxz = ts.dimensions[2]  # boxz
            boxx.append(tboxx)
            boxy.append(tboxy)
            boxz.append(tboxz)
            sa.append(tboxx*tboxy)
            norm += 1.0
        return boxx, boxy, boxz, sa, norm

    # setup 2d array
    def init_array2(n1: int, n2: int) -> Array2D[float]:
        array: Array2D[float] = []
        for i in range(0, n1):
            array.append([])
            for j in range(0, n2):
                array[i].append(0.0)  # initialize
        return array

    def init_array1(n1: int) -> Array1D[float]:
        array: Array1D[float] = []
        for i in range(0, n1):
            array.append(0.0)
        return array

    # calculate average & std of data in an array
    def get_ave_var(array: Array1D[float]) -> tuple[float, float]:
        ndata = len(array)
        ave, var = 0.0, 0.0
        for i in range(0, ndata):
            tmp = array[i]
            ave += tmp
            var += tmp**2
        # normlize
        norm = float(ndata)
        ave /= norm
        var = var/norm-ave**2
        return ave, var

    # -----------------
    #
    # MAIN
    #
    # -----------------
    u = mda.Universe(psf, traj)
    nsys = len(u.trajectory) // bs
    print(nsys)
    # ave. box size, sa, ka for each block
    arboxx = init_array1(nsys)
    vrboxx = init_array1(nsys)  # boxx
    arboxy = init_array1(nsys)
    vrboxy = init_array1(nsys)  # boxy
    arboxz = init_array1(nsys)
    vrboxz = init_array1(nsys)  # boxz
    arsa = init_array1(nsys)
    vrsa = init_array1(nsys)    # SA
    norm = init_array1(nsys)    # number of data points in each block
    rka0 = init_array1(nsys)    # Area compressibility modulus (within a block)
    rka = init_array1(nsys)     # area compressibility modulus (from overall SA)

    # ave. sa,ka for the system
    # Ave. SA and its std from simple average over blocks
    asa0, vsa0 = 0.0, 0.0
    # Ave. KA and its std from simple average over blocks
    aka0, vka0 = 0.0, 0.0
    asa, vsa = 0.0, 0.0  # Ave. SA and its standard error
    aka, vka = 0.0, 0.0
    sout = ""
    # loop over systems
    ntot = 0.0  # total sample size

    for j in range(nsys):
        # read dat
        boxx, boxy, boxz, sa, norm[j] = read_data(u, j)

        # get ave. & var.
        arboxx[j], vrboxx[j] = get_ave_var(boxx)
        arboxy[j], vrboxy[j] = get_ave_var(boxy)
        arboxz[j], vrboxz[j] = get_ave_var(boxz)
        arsa[j], vrsa[j] = get_ave_var(sa)

        tmp = arsa[j]
        #   update average SA & standard deviation from simple average over blocks
        asa0 += tmp
        vsa0 += tmp**2

    # update average SA & variance from the whole data
        asa += norm[j]*tmp
        # neeed additional term: norm[i]*(arsa[i]-arsa)**2
        vsa += norm[j]*vrsa[j]
        ntot += norm[j]

    # update asa,vsa
    asa /= ntot
    for i in range(0, nsys):
        vsa += norm[i]*(arsa[i]-asa)**2  # correction for combined variance
    # total variance V = \sum_i n_i [ V_i(x) + (a_i -a)^2 ] / \sum_i n_i
    vsa /= ntot

    # update asa, ssa (ssa for simple avrage)
    norm0 = float(nsys)
    asa0 /= norm0
    vsa0 = math.sqrt(vsa0/norm0-asa0**2)
    # Calculate KA
    for i in range(0, nsys):
        tmp = arsa[i]
        tmp = kBT*arsa[i]/vrsa[i]*1.e16
        rka0[i] = tmp
        aka0 += tmp
        vka0 += tmp**2  # from simple average over block
        # KA and STE from overall average SA
        tmp = kBT*arsa[i]/(vrsa[i]+(arsa[i]-asa)**2)*1.e16
        rka[i] = tmp
        aka += tmp
        vka += tmp**2    # to get mean KA and its STE
        sout += f'# block{i+1:d}: aboxx= {arboxx[i]:g} (A) vboxx= {vrboxx[i]:g} (A^2)\n' \
                f'# block{i+1:d}: aboxy= {arboxy[i]:g} (A) vboxy= {vrboxy[i]:g} (A^2)\n' \
                f'# block{i+1:d}: aboxz= {arboxz[i]:g} (A) vboxz= {vrboxz[i]:g} (A^2)\n' \
                f'# block{i+1:d}: asa= {arsa[i]:g} (A^2) vsa0= {vrsa[i]:g} vsa= ' \
                f'{vrsa[i]+(arsa[i]-asa)**2:G} (A^4)\n' \
                f'# block{i+1:d}: ka0= {rka0[i]:g} ka= {rka[i]:g} (in dyn/cm or mN/m)\n#\n'

    # update asa, ssa, aka, and ska
    norm0 = float(nsys)
    aka /= norm0
    vka = math.sqrt((vka/norm0-aka**2)/norm0)     # simple average
    aka0 /= norm0
    # average considering overall ave. SA
    vka0 = math.sqrt((vka0/norm0-aka0**2)/norm0)

    # SA from the whole data
    sout += "# KA from whole data\n"
    sout += f'{round(kBT * 1e+16 * asa / vsa, 1):g}\n'
    # Average from blocks
    vsa0 = math.sqrt(vsa0/norm0)  # standard error
    sout += "# Ave. SA & std\n" \
            f'{asa:g} {vsa0:g} {vsa:g} (in A)\n' \
            "# Ave. KA0 & ste , KA & ste\n" \
            f'{aka0:g} {vka0:g} {aka:g} {vka:g} (in dyn/cm or mN/m)\n'

    with sta.resolve_file(out, 'w') as outfile:
        header(outfile)
        print(sout, file=outfile)


def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog=f'stanalyzer {ANALYSIS_NAME}')
    sta.add_project_args(parser, 'psf', 'traj', 'out')
    parser.add_argument('-temp', type=float, help='System Temperature')
    parser.add_argument('-cntmin', type=int, help='first DCD number')
    parser.add_argument('-cntmax', type=int, help='Last DCD number')
    parser.add_argument('-bs', type=int, help='block size in frames')

    return parser


def main(settings: dict | None = None) -> None:
    if settings is None:
        settings = dict(sta.get_settings(ANALYSIS_NAME))

    write_compressibility_modulus(**settings)


if __name__ == '__main__':
    main()
