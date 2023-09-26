import csv
from math import acos, asin, copysign, cos, floor, log, pi, sin, sqrt, tan
from random import random

import numpy as np
from tabulate import tabulate
from tqdm import trange


def LaunchPhoton():
    w = 1.0

    x = 0.0
    y = 0.0
    z = 0.0

    dpx = 0.0
    dpy = 0.0
    dpz = 1.0

    return w, x, y, z, dpx, dpy, dpz


def RandomStep(mu_att):
    s = -log(random()) / mu_att

    return s


def BoundaryCheck(
    w, x, y, z, dpx, dpy, dpz, s, dDR, dDT, d, n_t, n_m, alpha_critical, epsilon=1e-9
):
    while True:
        if dpz < 0:
            sbound = (0 - z) / dpz
        elif dpz > 0:
            sbound = (d - z) / dpz
        else:
            break

        if sbound > s:
            break
        else:
            x += sbound * dpx
            y += sbound * dpy
            z += sbound * dpz

            s -= sbound

            alpha_i = acos(abs(dpz))

            if n_t > n_m:
                if alpha_i > alpha_critical:
                    dpz = -dpz
                    break
            else:
                alpha_t = asin(n_t * sin(alpha_i) / n_m)

                R = 0.5 * (
                    (sin(alpha_i - alpha_t)) ** 2
                    / ((sin(alpha_i + alpha_t)) ** 2 + epsilon)
                    + (tan(alpha_i - alpha_t)) ** 2
                    / ((tan(alpha_i + alpha_t)) ** 2 + epsilon)
                )

                if R < random():
                    if dpz < 0:
                        dDR += w
                    else:
                        dDT += w

                    w = 0

                dpz = -dpz

    return w, x, y, z, dpx, dpy, dpz, s, dDR, dDT


def TakeStep(x, y, z, dpx, dpy, dpz, s):
    x += s * dpx
    y += s * dpy
    z += s * dpz

    return x, y, z


def AbsorbPhoton(w, z, pDW, absFact, dz):
    dw = w * absFact

    if (w - dw) < 0:
        dw = w
    w -= dw

    nz = floor(z / dz)
    pDW[nz] += dw

    return w, pDW


def ScatterPhoton(dpx, dpy, dpz, g):
    if g == 0:
        theta = acos(2 * random() - 1)
    else:
        theta = acos(
            (1 + g**2 - ((1 - g**2) / (1 - g + 2 * g * random())) ** 2) / (2 * g)
        )

    psi = 2 * pi * random()

    if abs(dpz) > 0.9999:
        dpx = sin(theta) * cos(psi)
        dpy = sin(theta) * sin(psi)
        dpz = copysign(cos(theta), dpz)
    else:
        dpx = sin(theta) * (dpx * dpz * cos(psi) - dpy * sin(psi)) / sqrt(
            1 - dpz**2
        ) + dpx * cos(theta)
        dpy = sin(theta) * (dpy * dpz * cos(psi) - dpx * sin(psi)) / sqrt(
            1 - dpz**2
        ) + dpy * cos(theta)
        dpz = -sin(theta) * cos(psi) * sqrt(1 - dpz**2) + dpz * cos(theta)

    return dpx, dpy, dpz


def RandomWalk(mu_att, absFact, g, n_t, n_m, alpha_critical, d, dz, dDR, dDT, pDW):
    w, x, y, z, dpx, dpy, dpz = LaunchPhoton()

    while True:
        s = RandomStep(mu_att)
        w, x, y, z, dpx, dpy, dpz, s, dDR, dDT = BoundaryCheck(
            w, x, y, z, dpx, dpy, dpz, s, dDR, dDT, d, n_t, n_m, alpha_critical
        )

        if w <= 0:
            break

        x, y, z = TakeStep(x, y, z, dpx, dpy, dpz, s)
        w, pDW = AbsorbPhoton(w, z, pDW, absFact, dz)

        if w <= 0:
            break

        dpx, dpy, dpz = ScatterPhoton(dpx, dpy, dpz, g)

    return dDR, dDT, pDW


def read_parameters(filename):
    parameters = [[] for i in range(6)]
    with open(filename) as f:
        csv_reader = csv.reader(f, delimiter=",")
        next(csv_reader)
        for row in csv_reader:
            for i in range(6):
                parameters[i].append(float(row[i]))
    parameters = [np.array(p) for p in parameters]
    return tuple(parameters)


if __name__ == "__main__":
    nWavelengths, gs, mu_ss, mu_as, n_ts, n_ms = read_parameters("parameters.csv")
    d = 1.0
    nwl = nWavelengths.size

    mu_atts = mu_as + mu_ss
    absFacts = mu_as / mu_atts

    alpha_criticals = np.zeros(nwl)
    for n, (n_t, n_m) in enumerate(zip(n_ts, n_ms)):
        if n_t > n_m:
            alpha_criticals[n] = asin(n_m / n_t)
        else:
            alpha_criticals[n] = float("nan")

    N_bins = int(1e3)
    dz = d / N_bins

    pDW = np.zeros((nwl, N_bins))
    dDR = np.zeros(nwl)
    dDT = np.zeros(nwl)

    N_photons = int(1e5)


    for n, (mu_att, absFact, g, n_t, n_m, alpha_critical) in enumerate(
        zip(mu_atts, absFacts, gs, n_ts, n_ms, alpha_criticals)
    ):
        for nph in trange(N_photons):
            dDR[n], dDT[n], pDW[n] = RandomWalk(
                mu_att, absFact, g, n_t, n_m, alpha_critical, d, dz, dDR[n], dDT[n], pDW[n]
            )


    dDR /= N_photons
    dDT /= N_photons
    dDW = np.sum(pDW, axis=1) / N_photons

    print("\n")
    print(
        tabulate(
            {
                "Wavelengths": nWavelengths,
                "Diffuse Transmittance": dDT,
                "Diffuse Reflectance": dDR,
                "Absorption": dDW,
            },
            headers="keys",
        )
    )
    print("\n")

    pz = np.arange(0, d, dz)


    np.savetxt("pz", pz)
    np.savetxt("pDW", pDW)
