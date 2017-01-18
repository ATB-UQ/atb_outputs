from typing import Any, Dict, Callable, Optional, Tuple, List, Union
from logging import Logger

Atom = Dict[str, Any]

Ring = Dict[str, Any]

Coordinate = List[float]

MolData = Any

Output_File = str

Output_Files = List[Tuple[str, Output_File]]

Output_Module = Callable[[MolData], Union[Output_File, Output_Files]]

