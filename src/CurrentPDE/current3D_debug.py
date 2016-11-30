#!/usr/bin/env python
# -*- coding: utf-8 -*-
from ufl import *
set_level(DEBUG)

cell = tetrahedron

scalar = FiniteElement("Lagrange", cell, 1)

sigma1 = Coefficient(scalar)
w1 = TrialFunction(scalar)
w2 = TrialFunction(scalar)
f = Coefficient(scalar)


a1 = (-sigma1*inner(grad(w1), grad(q1))-w1*inner(grad(sigma1), grad(q1)))*(dx(1)+dx(2)+dx(3)+dx(4))
L1 = f*q1*(dx(2)+dx(3))

