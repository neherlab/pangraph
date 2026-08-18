"""Typed msgspec models for the pangraph JSON format.

These structs mirror the JSON schema that is generated from the Rust types (see
``pangraph_schema.py``). ``msgspec`` decodes a graph document directly into a
``PangraphData`` tree, checking structure, field types, the ``strand`` enum and
the non-negative integer ranges in a single pass, so no separate validation step
is needed.

The ``Uint`` alias carries the schema's ``minimum: 0`` (uint) constraint. Field
order follows dependency order for readability; forward references resolve
because annotations are strings (``from __future__ import annotations``).
"""

from __future__ import annotations

from typing import Annotated, Literal, Optional

import msgspec

# uint: a non-negative integer, matching the schema's {"type": "integer",
# "format": "uint", "minimum": 0}.
Uint = Annotated[int, msgspec.Meta(ge=0)]


class PangraphData(msgspec.Struct):
    """Root of a pangraph document: the paths, blocks and nodes, each keyed by id."""

    paths: dict[str, PangraphPath]
    blocks: dict[str, PangraphBlock]
    nodes: dict[str, PangraphNode]


class PangraphPath(msgspec.Struct):
    id: Uint
    nodes: list[Uint]
    tot_len: Uint
    circular: bool
    name: Optional[str] = None
    desc: Optional[str] = None


class PangraphBlock(msgspec.Struct):
    id: Uint
    consensus: str
    alignments: dict[str, Edit]


class PangraphNode(msgspec.Struct):
    id: Uint
    block_id: Uint
    path_id: Uint
    strand: Literal["+", "-"]
    position: tuple[Uint, Uint]


class Edit(msgspec.Struct):
    subs: list[Sub]
    dels: list[Del]
    inss: list[Ins]


class Sub(msgspec.Struct):
    pos: Uint
    alt: Annotated[str, msgspec.Meta(min_length=1, max_length=1)]


class Del(msgspec.Struct):
    pos: Uint
    len: Uint


class Ins(msgspec.Struct):
    pos: Uint
    seq: str
