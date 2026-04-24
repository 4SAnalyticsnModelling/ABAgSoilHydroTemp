import importlib.util
import sys
from pathlib import Path

_PACKAGE_DIR = Path(__file__).resolve().with_suffix("")
_INIT_FILE = _PACKAGE_DIR / "__init__.py"
_SPEC = importlib.util.spec_from_file_location(
    "ab_ag_soil_hydro_temp",
    _INIT_FILE,
    submodule_search_locations=[str(_PACKAGE_DIR)],
)
if _SPEC is None or _SPEC.loader is None:
    raise ImportError("Could not load ab_ag_soil_hydro_temp package.")

_MODULE = importlib.util.module_from_spec(_SPEC)
sys.modules["ab_ag_soil_hydro_temp"] = _MODULE
_SPEC.loader.exec_module(_MODULE)

main = _MODULE.main
run_model_many_site_years = _MODULE.run_model_many_site_years
run_one_site_year_crop = _MODULE.run_one_site_year_crop


if __name__ == "__main__":
    main()
