"""Choose initial lipid poses on the leaflet grid before Packmol optimization."""
from itertools import product
from math import acos, atan2, cos, pi, sin

import numpy as np

from .model import coordinates


class LipidSeedGrid:
    def __init__(self, plan):
        self.plan = plan
        self.step = plan.settings.tolerance / 2
        self.occupied = {}
        for placement in plan.placements:
            self.record(placement.xyz)

    def record(self, xyz):
        # Approximate the occupied volume at one-half of the packing tolerance.
        scaled = xyz / self.step
        grid = np.rint(scaled).astype(int)
        for offset in product((-1,0,1), repeat=3):
            nodes = grid + offset
            distances = np.sum((scaled - nodes) ** 2, axis=1)
            nearby = distances < 3
            for key, depth in zip(nodes[nearby].tolist(), (3 - distances[nearby]).tolist()):
                key = tuple(key)
                self.occupied[key] = max(depth,self.occupied.get(key,0.))

    def candidates(self, mol, side, head_bound, rng):
        """Sample poses around the molecular axis with a small optional tilt."""
        sign = 1 if side == 0 else -1
        xyz = np.array([coordinates(line) for line in mol.atoms])
        xyz -= xyz.mean(axis=0)
        ends = xyz[np.array(mol.tails)-1]
        axis = xyz[mol.head-1] - ends.mean(axis=0)
        axis /= np.linalg.norm(axis)
        cross = ends[-1]-ends[0]
        # Tilt toward the longer tail to bring the tail ends closer in Z.
        if cross @ axis > 0:
            cross = -cross
        cross -= (cross @ axis)*axis
        norm = np.linalg.norm(cross)
        lateral = cross / norm if norm > 1.e-6 else cross
        candidates = []
        tilts = (0., .21, .42) if norm > 1.e-6 else (0.,)
        outer = self.plan.interfaces[side] + min(3, self.plan.settings.water)
        tail_band = self.plan.settings.tail_band
        tail_indices = np.array(mol.tails) - 1
        for tilt in tilts:
            normal = sign*(cos(tilt)*axis+sin(tilt)*lateral)
            beta,theta = atan2(normal[0],-normal[1]),acos(np.clip(normal[2],-1.,1.))
            cb,sb,ct,st = cos(beta),sin(beta),cos(theta),sin(theta)
            phase = rng.uniform(0,2*pi)
            for j in range(12):
                gamma = phase + j*pi/6
                cg,sg = cos(gamma),sin(gamma)
                matrix = [[cb*cg-sb*sg*ct,-sb*cg*ct-cb*sg,sb*st],
                          [cb*sg*ct+sb*cg,cb*cg*ct-sb*sg,-cb*st],
                          [sg*st,cg*st,ct]]
                rotated = xyz @ matrix
                heights = sign*rotated[:,2]
                lo = max(head_bound + .02 - heights[mol.head - 1], -tail_band - heights.min())
                hi = min(tail_band - .02 - heights[tail_indices].max(), outer - heights.max())
                z = max(lo,min(hi,.6-heights.min()))
                shift = -(rotated.max(axis=0)+rotated.min(axis=0))/2
                shift[2] = sign*z
                positioned = rotated+shift
                candidates.append(((beta,gamma,theta), shift, positioned,
                                   50*max(0.,lo-hi)**2 + .2*max(0.,-heights.min()-z)**2))
        return candidates

    def choose(self, x, y, candidates, rng):
        # Randomize ties so isolated molecules have no preferred azimuth.
        rng.shuffle(candidates)
        best = float("inf")
        for angles,offset,relative,penalty in candidates:
            xyz = relative + (x,y,0)
            outside = np.maximum(0,np.abs(xyz[:,:2]) - (self.plan.x/2,self.plan.y/2))
            score = penalty + 10*np.sum(outside**2)
            if score >= best:
                continue
            nodes = np.rint(xyz/self.step).astype(int).tolist()
            score += sum(self.occupied.get(tuple(node),0.) for node in nodes)
            if score < best:
                best = score
                chosen = (x+offset[0],y+offset[1],offset[2],*angles)
                chosen_xyz = xyz
        self.record(chosen_xyz)
        return chosen
