#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Proveer excepciones personalizadas para las filtraciones."""

# %% InvalidGeometryError
class InvalidGeometryError(Exception):
    """Excepción para geometrías inválidas"""


# %% NotAPolygonError
class NotAPolygonError(Exception):
    """Excepción para elevar si no el objeto no es un polígono."""


# %% EmptyGeometryError
class EmptyGeometryError(Exception):
    """Excepción para geometrías vacías."""


# %% FiltrationError
class FiltrationError(Exception):
    """Excepción para filtraciones."""
