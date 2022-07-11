from pyomo import environ as po
from pyomo import dae as pod
import numpy as np


def create_model():
    model = po.ConcreteModel()

    model.t = pod.ContinuousSet(bounds=(0, 1))

    model.rna = po.Var(model.t, bounds=(0, None))
    model.drna_dt = pod.DerivativeVar(model.rna, wrt=model.t)

    model.

    return model
