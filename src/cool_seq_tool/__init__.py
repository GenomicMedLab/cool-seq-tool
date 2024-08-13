"""The cool_seq_tool package"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("cool_seq_tool")
except PackageNotFoundError:
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError


from cool_seq_tool.app import CoolSeqTool

__all__ = ["CoolSeqTool", "__version__"]
