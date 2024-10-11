import json
import math
import time
from abc import ABC

import numpy as np
from scipy import optimize


class DiesPerWaferCalculator(ABC):
    def __init__(
        self,
        width,
        height,
        xspacing,
        yspacing,
        waferdiameter,
        edgeexclusionwidth,
        ft_grid,
        searchdepth,
        symmetric,
        ft_ShiftRows,
        ft_ShiftCols,
        ft_ShiftRot,
    ) -> None:
        super().__init__()

        self.width = width
        self.height = height
        self.xspacing = xspacing
        self.yspacing = yspacing
        self.waferdiameter = waferdiameter
        self.edgeexclusionwidth = edgeexclusionwidth
        self.ft_grid = ft_grid
        self.searchdepth = searchdepth
        self.symmetric = symmetric
        self.ft_ShiftRows = ft_ShiftRows
        self.ft_ShiftCols = ft_ShiftCols
        self.ft_ShiftRot = ft_ShiftRot

        self.allow_rotation = True

        self.no_local_search_shift = (
            True  # no need for local search when output is integer
        )

        self.maxiter_grid0 = 1000
        self.maxiter_grid1 = 2000
        self.maxiter_grid2 = 4000

        self.maxiter_shift0 = 1000
        self.maxiter_shift1 = 10000
        self.maxiter_shift2 = 60000

        ewr = self.waferdiameter / 2 - self.edgeexclusionwidth
        self.ewr = ewr
        # %% find best solution
        self.Nx = math.ceil(ewr / self.width)  # spacing ignored
        self.Ny = math.ceil(ewr / self.height)
        self.Nmax = max((self.Nx, self.Ny))
        self.Xoff, self.Yoff = np.meshgrid(
            range(-self.Nx, self.Nx + 1), range(-self.Ny, self.Ny + 1), indexing="ij"
        )  # grids of values [-Nx,...,Nx],[-Ny,...,Ny]
        self.Xoff2, self.Yoff2 = np.meshgrid(
            range(-self.Nmax, self.Nmax + 1),
            range(-self.Nmax, self.Nmax + 1),
            indexing="ij",
        )  # grid for shift&rotate

    def CalculatePositions(self, offsets, ft):
        if ft == 0:
            if self.symmetric:
                return (
                    (offsets[0] > 0.5) * (self.width + self.xspacing) / 2
                    + self.Xoff * (self.width + self.xspacing)
                    - self.width / 2,
                    (offsets[1] > 0.5) * (self.height + self.yspacing) / 2
                    + self.Yoff * (self.height + self.yspacing)
                    - self.height / 2,
                    self.width * np.ones((2 * self.Nx + 1, 2 * self.Ny + 1)),
                    self.height * np.ones((2 * self.Nx + 1, 2 * self.Ny + 1)),
                )
            else:
                return (
                    offsets[0]
                    + self.Xoff * (self.width + self.xspacing)
                    - self.width / 2,
                    offsets[1]
                    + self.Yoff * (self.height + self.yspacing)
                    - self.height / 2,
                    self.width * np.ones((2 * self.Nx + 1, 2 * self.Ny + 1)),
                    self.height * np.ones((2 * self.Nx + 1, 2 * self.Ny + 1)),
                )
        elif ft == 1:
            if self.symmetric:
                if offsets[:1] < 0.5:  # center row unique
                    mirror_offset = np.concatenate((np.flip(offsets[2:]), offsets[1:]))
                else:  # center row copied
                    mirror_offset = np.concatenate(
                        (np.flip(offsets[1:-1]), offsets[1:])
                    )
                xoffsets = (mirror_offset > 0.5) * (self.width + self.xspacing) / 2
                yoffset = (offsets[:1] > 0.5) * (self.height + self.yspacing) / 2
            else:
                xoffsets = (offsets[1:] > 0.5) * (self.width + self.xspacing) / 2
                yoffset = offsets[0]
            return (
                np.tile(xoffsets, (2 * self.Nx + 1, 1))
                + self.Xoff * (self.width + self.xspacing)
                - self.width / 2,
                yoffset + self.Yoff * (self.height + self.yspacing) - self.height / 2,
                self.width * np.ones((2 * self.Nx + 1, 2 * self.Ny + 1)),
                self.height * np.ones((2 * self.Nx + 1, 2 * self.Ny + 1)),
            )
        elif ft == 2:
            if self.symmetric:
                if offsets[:1] < 0.5:  # center col unique
                    mirror_offset = np.concatenate((np.flip(offsets[2:]), offsets[1:]))
                else:  # center col copied
                    mirror_offset = np.concatenate(
                        (np.flip(offsets[1:-1]), offsets[1:])
                    )
                yoffsets = (mirror_offset > 0.5) * (self.height + self.yspacing) / 2
                xoffset = (offsets[:1] > 0.5) * (self.width + self.xspacing) / 2
            else:
                yoffsets = (offsets[1:] > 0.5) * (self.height + self.yspacing) / 2
                xoffset = offsets[0]
            return (
                xoffset + self.Xoff * (self.width + self.xspacing) - self.width / 2,
                np.tile(yoffsets, (2 * self.Ny + 1, 1)).transpose()
                + self.Yoff * (self.height + self.yspacing)
                - self.height / 2,
                self.width * np.ones((2 * self.Nx + 1, 2 * self.Ny + 1)),
                self.height * np.ones((2 * self.Nx + 1, 2 * self.Ny + 1)),
            )
        elif ft == 3:
            if self.symmetric:
                if offsets[:1] < 0.5:  # center row unique
                    xshifted_bool = (
                        np.concatenate(
                            (
                                np.flip(offsets[2 : self.Nmax + 2]),
                                offsets[1 : self.Nmax + 2],
                            )
                        )
                        > 0.5
                    )
                    rotated_bool = (
                        np.concatenate(
                            (
                                np.flip(offsets[self.Nmax + 3 :]),
                                offsets[self.Nmax + 2 :],
                            )
                        )
                        > 0.5
                    )
                else:  # center row copied
                    xshifted_bool = (
                        np.concatenate(
                            (
                                np.flip(offsets[1 : self.Nmax + 1]),
                                offsets[1 : self.Nmax + 2],
                            )
                        )
                        > 0.5
                    )
                    rotated_bool = (
                        np.concatenate(
                            (
                                np.flip(offsets[self.Nmax + 2 : -1]),
                                offsets[self.Nmax + 2 :],
                            )
                        )
                        > 0.5
                    )
                yoffset0 = (
                    (offsets[0] > 0.5)
                    * (
                        (
                            self.height * ~rotated_bool[self.Nmax]
                            + self.width * rotated_bool[self.Nmax]
                        )
                        + self.yspacing
                    )
                    / 2
                )
            else:
                xshifted_bool = offsets[1 : 2 * self.Nmax + 2] > 0.5
                rotated_bool = offsets[2 * self.Nmax + 2 :] > 0.5
                yoffset0 = offsets[0]
            xextent = ~rotated_bool * self.width + rotated_bool * self.height  # 1D
            yextent = rotated_bool * self.width + ~rotated_bool * self.height  # 1D
            xoffset0 = (
                xshifted_bool * self.xspacing / 2 - ~xshifted_bool * xextent / 2
            )  # 1D
            Rsum = np.zeros(2 * self.Nmax + 1)  # 1D
            Rsum[self.Nmax + 1 :] = rotated_bool[self.Nmax : -1].cumsum()
            Rsum[: self.Nmax] = -np.flip(np.flip(rotated_bool[: self.Nmax]).cumsum())
            notRsum = np.arange(-self.Nmax, self.Nmax + 1) - Rsum
            yoffsets = (
                yoffset0
                + self.yspacing * np.arange(-self.Nmax, self.Nmax + 1)
                - (
                    self.height * ~rotated_bool[self.Nmax]
                    + self.width * rotated_bool[self.Nmax]
                )
                / 2
                + Rsum * self.width
                + notRsum * self.height
            )
            R = np.tile(rotated_bool, (2 * self.Nmax + 1, 1))
            return (
                np.tile(xoffset0, (2 * self.Nmax + 1, 1))
                + self.Xoff2 * (~R * self.width + R * self.height + self.xspacing),
                np.tile(yoffsets, (2 * self.Nmax + 1, 1)),
                np.tile(xextent, (2 * self.Nmax + 1, 1)),
                np.tile(yextent, (2 * self.Nmax + 1, 1)),
            )

    def constructV2(
        self, offsets, ft
    ):  # return distance^2 of each corner from center of wafer
        Xcorner, Ycorner, Xextent, Yextent = self.CalculatePositions(offsets, ft)
        return np.stack(
            (
                Xcorner**2 + Ycorner**2,
                (Xcorner + Xextent) ** 2 + Ycorner**2,
                Xcorner**2 + (Ycorner + Yextent) ** 2,
                (Xcorner + Xextent) ** 2 + (Ycorner + Yextent) ** 2,
            )
        )

    def gridwithpartialscore(self, offsets, ft):
        V = self.constructV2(offsets, ft)
        valid = np.all(V <= self.ewr**2, axis=0)
        PS = (
            1
            - np.sum(np.clip(V**0.5 - self.ewr, 0, np.inf), axis=0)
            / 4
            / (self.width**2 + self.height**2) ** 0.5
        )
        partial_score = max(PS[~valid])
        return -(valid.sum() + partial_score)

    def countfulldies(self, offsets, ft):
        V = self.constructV2(offsets, ft)
        valid = np.all(V <= self.ewr**2, axis=0)
        return -valid.sum()

    def Rmax(self, offsets, ft, validmask):
        if self.fittype == 0:
            V = self.constructV2(offsets, ft)
        else:
            self.offsets_test[0] = offsets  # just 1 offset
            V = self.constructV2(self.offsets_test, ft)
        return V[validmask].max()

    def fit(self):
        if self.searchdepth == 0:
            maxiter_grid = self.maxiter_grid0
            maxiter_shift = self.maxiter_shift0
        elif self.searchdepth == 1:
            maxiter_grid = self.maxiter_grid1
            maxiter_shift = self.maxiter_shift1
        elif self.searchdepth == 2:
            maxiter_grid = self.maxiter_grid2
            maxiter_shift = self.maxiter_shift2
        Nfits = np.zeros(4)

        start = time.time()
        start2 = start
        if self.ft_grid:
            if self.symmetric:
                # there are literally only 4 solutions in this case, so this is kinda dumb
                fit0 = optimize.dual_annealing(
                    self.countfulldies,
                    [(0, 1), (0, 1)],
                    args=(0,),
                    maxiter=maxiter_grid,
                )
            else:
                fit0 = optimize.dual_annealing(
                    self.gridwithpartialscore,
                    [
                        (0, (self.width + self.xspacing) / 2),
                        (0, (self.height + self.yspacing) / 2),
                    ],
                    args=(0,),
                    maxiter=maxiter_grid,
                )
            print("Uniform Grid fit output:")
            print(fit0)
            print("Uniform Grid time: {:.2f} sec".format(time.time() - start2))
            Nfits[0] = math.floor(-fit0.fun)
            start2 = time.time()
        if self.ft_ShiftRows:
            if self.symmetric:
                fit1 = optimize.dual_annealing(
                    self.countfulldies,
                    [(0, 1)] + [(0, 1) for _ in range(0, self.Ny + 1)],
                    args=(1,),
                    maxiter=maxiter_shift,
                    no_local_search=self.no_local_search_shift,
                )
            else:
                fit1 = optimize.dual_annealing(
                    self.countfulldies,
                    [(0, (self.height + self.yspacing) / 2)]
                    + [(0, 1) for _ in range(-self.Ny, self.Ny + 1)],
                    args=(1,),
                    maxiter=maxiter_shift,
                    no_local_search=self.no_local_search_shift,
                )
            print("Shift Rows fit output:")
            print(fit1)
            print("Shift Rows time: {:.2f} sec".format(time.time() - start2))
            Nfits[1] = math.floor(-fit1.fun)
            start2 = time.time()
        if self.ft_ShiftCols:
            if self.symmetric:
                fit2 = optimize.dual_annealing(
                    self.countfulldies,
                    [(0, 1)] + [(0, 1) for _ in range(0, self.Nx + 1)],
                    args=(2,),
                    maxiter=maxiter_shift,
                    no_local_search=self.no_local_search_shift,
                )
            else:
                fit2 = optimize.dual_annealing(
                    self.countfulldies,
                    [(0, (self.width + self.xspacing) / 2)]
                    + [(0, 1) for _ in range(-self.Nx, self.Nx + 1)],
                    args=(2,),
                    maxiter=maxiter_shift,
                    no_local_search=self.no_local_search_shift,
                )
            print("Shift Columns fit output:")
            print(fit2)
            print("Shift Columns time: {:.2f} sec".format(time.time() - start2))
            Nfits[2] = math.floor(-fit2.fun)
            start2 = time.time()
        if self.ft_ShiftRot:
            if self.symmetric:
                fit3 = optimize.dual_annealing(
                    self.countfulldies,
                    [(0, 1)] + [(0, 1) for _ in range(2 * self.Nmax + 2)],
                    args=(3,),
                    maxiter=maxiter_shift,
                    no_local_search=self.no_local_search_shift,
                )
            else:
                fit3 = optimize.dual_annealing(
                    self.countfulldies,
                    [(0, (max((self.height, self.width)) + self.xspacing) / 2)]
                    + [(0, 1) for _ in range(4 * self.Nmax + 2)],
                    args=(3,),
                    maxiter=maxiter_shift,
                    no_local_search=self.no_local_search_shift,
                )
            print("Shifted & Rotated fit output:")
            print(fit3)
            print("Shifted & Rotated time: {:.2f} sec".format(time.time() - start2))
            Nfits[3] = math.floor(-fit3.fun)
            start2 = time.time()

        fittype = (
            Nfits.argmax()
        )  # "indices corresponding to the first occurrence are returned"
        self.fittype = fittype
        if fittype == 0:
            fit = fit0
        elif fittype == 1:
            fit = fit1
        elif fittype == 2:
            fit = fit2
        elif fittype == 3:
            fit = fit3
        # fit=[fit1,fit2,fit3,fit4][fittype] #doesnt work bc not all fit# defined
        end = time.time()
        Nfit = math.floor(-fit.fun)
        # %% center solution (if not symmetric)
        fit_offsets = fit.x
        fit_V = self.constructV2(fit_offsets, fittype)
        fit_valid = np.tile(
            np.all(fit_V <= self.ewr**2, axis=0), [4, 1, 1]
        )  # mask of valid points in V2
        if self.symmetric:
            final_offsets = fit.x
            min_diameter = 2 * (fit_V[fit_valid].max() ** 0.5 + self.edgeexclusionwidth)
        else:
            self.offsets_test = fit_offsets

            print("\nCentering output:")
            if fittype == 0:
                fit2 = optimize.minimize(
                    self.Rmax,
                    fit_offsets,
                    bounds=[
                        (0, (self.width + self.xspacing) / 2),
                        (0, (self.height + self.yspacing) / 2),
                    ],
                    args=(fittype, fit_valid),
                )
                centered_offsets = fit2.x
            elif fittype == 1:
                fit2 = optimize.minimize_scalar(
                    self.Rmax,
                    bounds=(0, (self.height + self.yspacing) / 2),
                    args=(fittype, fit_valid),
                )
                centered_offsets = fit_offsets
                centered_offsets[0] = fit2.x
            elif fittype == 2:
                fit2 = optimize.minimize_scalar(
                    self.Rmax,
                    bounds=(0, (self.width + self.xspacing) / 2),
                    args=(fittype, fit_valid),
                )
                centered_offsets = fit_offsets
                centered_offsets[0] = fit2.x
            elif fittype == 3:
                fit2 = optimize.minimize_scalar(
                    self.Rmax,
                    bounds=(0, (max((self.width, self.height)) + self.yspacing) / 2),
                    args=(fittype, fit_valid),
                )
                centered_offsets = fit_offsets
                centered_offsets[0] = fit2.x
            print(fit2)
            final_offsets = centered_offsets
            min_diameter = 2 * (fit2.fun**0.5 + self.edgeexclusionwidth)
        print("\nFinal Offsets:")
        print(final_offsets)
        print("Solution diameter:")
        print(min_diameter)
        area_utilization = (
            Nfit * self.width * self.height * 4 * self.waferdiameter**-2 / np.pi
        )

        Xcorner, Ycorner, Xextent, Yextent = self.CalculatePositions(
            final_offsets, fittype
        )
        # Xcen,Ycen=calcXYcen(final_offsets,fittype)
        V2 = self.constructV2(final_offsets, fittype)

        self.valid = np.all(V2**0.5 <= self.ewr, axis=0)
        self.partial = np.logical_and(
            np.any(V2**0.5 <= self.ewr, axis=0), np.logical_not(self.valid)
        )

        self.final_diameter = min_diameter
        self.Xcorner = Xcorner.astype(float)
        self.Ycorner = Ycorner.astype(float)
        self.Xextent = Xextent.astype(float)
        self.Yextent = Yextent.astype(float)

    def format_json_obj(self):
        json_data = {}
        # User inputs
        json_data["user_inputs"] = {
            "width": float(self.width),
            "height": float(self.height),
            "xspacing": float(self.xspacing),
            "yspacing": float(self.yspacing),
            "edge_exclusion_width": float(self.edgeexclusionwidth),
        }
        # Outputs
        json_data["final_wafer_diameter"] = float(self.final_diameter)
        json_data["fit_type"] = int(self.fittype)

        partial_dies = []
        valid_dies = []
        for idx, _ in np.ndenumerate(self.Xcorner):
            corner_pts = [
                [self.Xcorner[idx], self.Ycorner[idx]],
                [self.Xcorner[idx] + self.Xextent[idx], self.Ycorner[idx]],
                [
                    self.Xcorner[idx] + self.Xextent[idx],
                    self.Ycorner[idx] + self.Yextent[idx],
                ],
                [self.Xcorner[idx], self.Ycorner[idx] + self.Yextent[idx]],
            ]
            if self.partial[idx]:
                partial_dies.append(corner_pts)
            if self.valid[idx]:
                valid_dies.append(corner_pts)
        json_data["num_dies"] = len(valid_dies)
        json_data["valid_dies"] = valid_dies
        json_data["partial_dies"] = partial_dies
        print(len(partial_dies))
        print(len(valid_dies))
        return json.dumps(json_data)
