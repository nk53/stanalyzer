from typing import Any, TypeAlias

# the 'type' keyword was added in 3.11; stanalyzer aims to support 3.10
StrDict: TypeAlias = dict[str, Any]
StrDictList: TypeAlias = list[dict[str, Any]]
