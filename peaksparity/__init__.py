"""Backward-compatible alias: the package was renamed to ``proteoscan``. ``import proteoscan.x``
keeps working by forwarding to ``proteoscan.x``."""
import importlib as _il
import sys as _sys

import proteoscan as _p

__version__ = _p.__version__
__path__ = _p.__path__  # submodule imports resolve inside proteoscan/


class _Finder:
    @staticmethod
    def find_spec(name, path=None, target=None):
        if name.startswith("proteoscan."):
            real = _il.import_module("proteoscan." + name.split(".", 1)[1])
            _sys.modules[name] = real
            import importlib.machinery as m
            return m.ModuleSpec(name, None)
        return None


_sys.meta_path.insert(0, _Finder())
