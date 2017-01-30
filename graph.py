from operator import itemgetter
from io import StringIO, BytesIO
from numpy import linspace, random
from colorsys import hsv_to_rgb
from math import sqrt, ceil
from typing import Optional, Any, Union, Tuple, Dict

USE_NEIGHBOUR_VALENCES = True

USE_DIFFERENT_ELEMENT_SHAPES = False

OPAQUE = 1.0

RAISE_IF_MISSING_GRAPH_TOOL = False

def graph(data, vertex_text=None, decorate_graph=True, use_random_colors: bool = True, shuffle_equivalence_classes: bool = True):
    try:
        from graph_tool.all import Graph
    except ImportError as e:
        if RAISE_IF_MISSING_GRAPH_TOOL:
            raise
        else:
            from sys import stderr
            stderr.write('Package graph_tool could not be imported. Error was: {0}'.format(e))
            return None

    g = Graph(directed=False)

    vertex_types = g.new_vertex_property("string")
    g.vertex_properties['type'] = vertex_types

    if decorate_graph:
        vertex_equivalence_classes = g.new_vertex_property("vector<double>")
        g.vertex_properties['equivalence_class'] = vertex_equivalence_classes

        if USE_DIFFERENT_ELEMENT_SHAPES:
            vertex_shapes = g.new_vertex_property("int")
            g.vertex_properties['shape'] = vertex_shapes

        equivalence_classes = [
            atom['equivalenceGroup']
            for atom in list(data.atoms.values())
        ]

        unique_equivalence_classes = {}
        assigned_equivalence_classes = {}
        next_equivalence_class = 0
        for (atom_index, atom) in sorted(list(data.atoms.items()), key=itemgetter(0)):
            if atom['equivalenceGroup'] == -1:
                unique_equivalence_classes[atom_index] = next_equivalence_class
                next_equivalence_class += 1
            else:
                if atom['equivalenceGroup'] not in assigned_equivalence_classes:
                    assigned_equivalence_classes[atom['equivalenceGroup']] = next_equivalence_class
                    next_equivalence_class += 1
                unique_equivalence_classes[atom_index] = assigned_equivalence_classes[atom['equivalenceGroup']]

        range_unique_equivalence_classes = list(range(len(set(unique_equivalence_classes.values()))))

        equivalence_class_permutation = dict(
            list(zip(
                range_unique_equivalence_classes,
                (random.permutation(range_unique_equivalence_classes) if shuffle_equivalence_classes else range_unique_equivalence_classes), # pylint: disable=no-member
            )),
        )

        unique_elements = set([atom['type'] for atom in list(data.atoms.values())])
        symbol_index_for_element = dict(
            list(map(
                lambda n_element: (n_element[1], n_element[0] % 15), # There are only 15 different symbols
                enumerate(unique_elements),
            )),
        )

        hues = linspace( # pylint: disable=no-member
            0.0,
            1.0,
            len(
                set(
                    list(unique_equivalence_classes.values()),
                ),
            ) // 2 + 1,
        ).tolist()

        colours = list(map(
            lambda c: hsv_to_rgb(*c),
            [(H, S, 0.5) for H in hues for S in (0.3, 0.9)],
        ))

        # Shuffle the colours around the graph
        colours = (
            [
                # Numpy cast all tuples to lists; needs to be reverted
                tuple(colour)
                for colour in random.permutation(colours).tolist()
            ]
            if use_random_colors
            else colours
        )

    if vertex_text == 'element_valence':
        vertex_text_fct = lambda atom_index, atom: '{element}{valence}'.format(
            element=atom['type'],
            valence=len(atom['conn']) if USE_NEIGHBOUR_VALENCES else '',
        )
    elif vertex_text == 'element_equivalence':
        vertex_text_fct = lambda atom_index, atom: '{element}{equivalence}'.format(
            element=atom['type'],
            equivalence=unique_equivalence_classes[atom_index],
        )
    elif vertex_text == 'name_equivalence':
        vertex_text_fct = lambda atom_index, atom: '{symbol}({equivalence})'.format(
            symbol=atom['symbol'],
            equivalence=unique_equivalence_classes[atom_index],
        )
    else:
        raise Exception('Unvalid vertex_text')

    vertices = {}
    for (atom_index, atom) in sorted(list(data.atoms.items()), key=itemgetter(0)):
        v = g.add_vertex()
        vertex_types[v] = vertex_text_fct(atom_index, atom)
        vertices[atom_index] = v

        if decorate_graph:
            vertex_equivalence_classes[v] = colours[unique_equivalence_classes[atom_index]] + (OPAQUE,)
            if USE_DIFFERENT_ELEMENT_SHAPES:
                #vertex_shapes[v] = symbol_index_for_element[atom['type']]
                pass

    for (i, j) in [bond['atoms'] for bond in data.bonds]:
        g.add_edge(vertices[i], vertices[j])

    return g

def graph_img(data, molecule_graph: Optional[Any] = None, return_pos: bool = False, pos: Optional[Any] = None, **kwargs: Dict[str, Any]) -> Union[None, str, Tuple[str, Any, Any]]:
    try:
        from graph_tool.draw import graph_draw
    except ImportError as e:
        if RAISE_IF_MISSING_GRAPH_TOOL:
            raise
        else:
            from sys import stderr
            stderr.write('Package graph_tool could not be imported. Error was: {0}'.format(e))
            return None

    if molecule_graph is None:
        molecule_graph = graph(data, **kwargs)

    if molecule_graph is not None:
        io = StringIO()

        output_pos = graph_draw(
            molecule_graph,
            pos=pos,
            vertex_text=molecule_graph.vertex_properties['type'],
            vertex_font_size=10.0,
            vertex_size=35.0,
            output=io,
            fmt='svg',
            output_size=tuple([150 * ceil(sqrt(len(data.atoms)))]*2),
            vertex_fill_color=molecule_graph.vertex_properties['equivalence_class'],
            **(
                dict(vertex_shape=molecule_graph.vertex_properties['shape'])
                if USE_DIFFERENT_ELEMENT_SHAPES
                else {}
            )
        )

        if not return_pos:
            return io.getvalue()
        else:
            return (io.getvalue(), molecule_graph, output_pos)
    else:
        return None

def graph_gt(data, **kwargs):
    molecule_graph = graph(data, vertex_text='element_valence', decorate_graph=False)

    if molecule_graph is not None:
        io = BytesIO()
        molecule_graph.save(io)
        return io.getvalue()
    else:
        return None
