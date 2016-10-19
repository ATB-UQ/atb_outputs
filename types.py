from typing import Callable, Any, Union, Tuple, Dict, List

MolData = Any

Output_File = str

Output_Files = List[Tuple[str, Output_File]]

Output_Module = Callable[[MolData], Union[Output_File, Output_Files]]
