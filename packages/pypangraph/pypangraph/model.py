"""Typed pydantic models for the pangraph JSON format.

These models mirror the JSON schema that is generated from the Rust types (see
``pangraph_schema.py``). ``PangraphData.model_validate_json`` parses and validates
a graph document in one call, checking structure, field types, the ``strand`` enum
and the non-negative integer ranges, so no separate validation step is needed.

The ``Uint`` alias carries the schema's ``minimum: 0`` (uint) constraint. Classes
are defined in dependency order so every referenced model exists at definition
time and no forward-reference rebuild is required.
"""

from typing import Annotated, Literal, Optional

from pydantic import BaseModel, Field

# uint: a non-negative integer, matching the schema's {"type": "integer",
# "format": "uint", "minimum": 0}.
Uint = Annotated[int, Field(ge=0)]


class Sub(BaseModel):
    pos: Uint
    alt: Annotated[str, Field(min_length=1, max_length=1)]


class Del(BaseModel):
    pos: Uint
    len: Uint


class Ins(BaseModel):
    pos: Uint
    seq: str


class Edit(BaseModel):
    subs: list[Sub]
    dels: list[Del]
    inss: list[Ins]


class PangraphPath(BaseModel):
    id: Uint
    nodes: list[Uint]
    tot_len: Uint
    circular: bool
    name: Optional[str] = None
    desc: Optional[str] = None


class PangraphBlock(BaseModel):
    id: Uint
    consensus: str
    alignments: dict[str, Edit]


class PangraphNode(BaseModel):
    id: Uint
    block_id: Uint
    path_id: Uint
    strand: Literal["+", "-"]
    position: tuple[Uint, Uint]


class PangraphData(BaseModel):
    """Root of a pangraph document: the paths, blocks and nodes, each keyed by id."""

    paths: dict[str, PangraphPath]
    blocks: dict[str, PangraphBlock]
    nodes: dict[str, PangraphNode]
