import gzip
import importlib.util
import importlib.machinery
import types
import sys
from pathlib import Path

# Create dummy package structure to load file_utils without top-level imports
package_name = "organoverse"
utils_pkg_name = f"{package_name}.utils"

sys.modules.setdefault(package_name, types.ModuleType(package_name))
sys.modules.setdefault(utils_pkg_name, types.ModuleType(utils_pkg_name))

base_path = Path(__file__).resolve().parents[1] / "organoverse" / "utils"

# Load exceptions module
exc_loader = importlib.machinery.SourceFileLoader(f"{utils_pkg_name}.exceptions", str(base_path / "exceptions.py"))
exc_spec = importlib.util.spec_from_loader(exc_loader.name, exc_loader)
exc_module = importlib.util.module_from_spec(exc_spec)
exc_loader.exec_module(exc_module)
sys.modules[f"{utils_pkg_name}.exceptions"] = exc_module

# Load file_utils module
fu_loader = importlib.machinery.SourceFileLoader(f"{utils_pkg_name}.file_utils", str(base_path / "file_utils.py"))
fu_spec = importlib.util.spec_from_loader(fu_loader.name, fu_loader)
file_utils = importlib.util.module_from_spec(fu_spec)
fu_loader.exec_module(file_utils)

is_compressed = file_utils.is_compressed


def test_is_compressed(tmp_path):
    plain = tmp_path / "plain.txt"
    plain.write_text("data")
    assert is_compressed(str(plain)) is False

    gz_path = tmp_path / "data.txt.gz"
    with gzip.open(gz_path, "wb") as fh:
        fh.write(b"data")
    assert is_compressed(str(gz_path)) is True
