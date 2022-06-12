#!/usr/bin/env python3


# Małe grafy proste reprezentowane jako para złożona ze zbioru wierzchołków i zbioru krawędzi.
class Graph:

    # Tworzy graf o podanej reprezentacji tekstowej (domyślnie: graf zerowy).
    def __init__( self, text = "?" ):
        self.__vertices, self.__edges, k, i = set( range( ord( text[0] ) - 63 ) ), set(), 0, 1
        for v in range( 1, len( self.__vertices ) ):
            for u in range( v ):
                if k == 0:
                    c, i, k = ord( text[i] ) - 63, i + 1, 6
                k -= 1
                if (c & (1 << k)) != 0:
                    self.__edges.add( (u, v) )

    # OPERACJA PODSTAWOWA
    # Zwraca liczbę wierzchołków grafu.
    def number_of_vertices( self ): 
        return len( self.__vertices )

    # OPERACJA PODSTAWOWA
    # Zwraca zbiór wierzchołków grafu.
    def vertices( self ):
        return self.__vertices

    # OPERACJA PODSTAWOWA
    # Zwraca stopień wierzchołka.
    def vertex_degree( self, vertex ):
        return len( [u + v - vertex for u, v in self.__edges if u == vertex or v == vertex] )

    # OPERACJA PODSTAWOWA
    # Zwraca sąsiedztwo podanego wierzchołka.
    def vertex_neighbors( self, vertex ):
        return set( u + v - vertex for u, v in self.__edges if u == vertex or v == vertex )

    # OPERACJA PODSTAWOWA
    # Dodaje do grafu podany wierzchołek.
    def add_vertex( self, vertex ):
        self.__vertices.add( vertex )

    # OPERACJA PODSTAWOWA
    # Usuwa z grafu podany wierzchołek i wszystkie incydentne z nim krawędzie.
    def delete_vertex( self, vertex ):
        self.__vertices.discard( vertex )
        self.__edges -= set( (u, v) for u, v in self.__edges if u == vertex or v == vertex )

    # OPERACJA PODSTAWOWA
    # Zwraca liczbę krawędzi grafu.
    def number_of_edges( self ):
        return len( self.__edges )

    # OPERACJA PODSTAWOWA
    # Zwraca zbiór krawędzi grafu.
    def edges( self ):
        return self.__edges

    # OPERACJA PODSTAWOWA
    # Zwraca informację o tym, czy podane wierzchołki z sobą sąsiadują.
    def is_edge( self, u, v ):
        return (min( u, v ), max( u, v )) in self.__edges

    # OPERACJA PODSTAWOWA
    # Dodaje do grafu podaną krawędź.
    def add_edge( self, u, v ):
        self.__edges.add( (min( u, v ), max( u, v )) )

    # OPERACJA PODSTAWOWA
    # Usuwa z grafu podaną krawędź.
    def delete_edge( self, u, v ):
        self.__edges.discard( (min( u, v ), max( u, v )) )

    # Sprawdza, czy podane grafy są identyczne.
    def __eq__( self, graph ):
        return self.__vertices == graph.vertices() and self.__edges == graph.edges()

    # Sprawdza, czy podane grafy są różne.
    def __ne__( self, graph ):
        return self.__vertices != graph.vertices() or self.__edges != graph.edges()

    # OPERACJA DODATKOWA
    # Zwraca dopełnienie grafu.
    def complement( self ):
        g = Graph()
        for v in self.__vertices:
            g.add_vertex( v )
            for u in self.__vertices:
                if u < v and (u, v) not in self.__edges:
                    g.add_edge( u, v )
        return g

    # OPERACJA DODATKOWA
    # Zwraca liczbę składowych spójności grafu.
    def connected_components( self ):
        c, cs = 0, {v: 0 for v in self.__vertices}
        def mark_component( vertex, mark ):
            cs[vertex] = mark
            for v in self.__vertices:
                if cs[v] == 0 and self.is_edge( vertex, v ):
                    mark_component( v, mark )
        for v in self.__vertices:
            if cs[v] == 0:
                c += 1
                mark_component( v, c )
        return c

    # OPERACJA DODATKOWA
    # Zwraca graf dwudzielny pełny Kn,m.
    @staticmethod
    def create_complete_bipartite( n, m ):
        g = Graph()
        for i in range( n + m ):
            g.add_vertex( i )
        for i in range( n ):
            for j in range( n, n + m ):
                g.add_edge( i, j )
        return g

    # OPERACJA DODATKOWA
    # Zwraca cykl o podanej liczbie wierzchołków.
    @staticmethod
    def create_cycle( number_of_vertices ):
        g = Graph()
        for i in range( number_of_vertices ):
            g.add_vertex( i )
            g.add_edge( i, (i + 1) % number_of_vertices )
        return g

    # OPERACJA DODATKOWA
    # Zwraca ścieżkę o podanej liczbie wierzchołków.
    @staticmethod
    def create_path( number_of_vertices ):
        g = Graph()
        for i in range( number_of_vertices ):
            g.add_vertex( i )
        for i in range( 1, number_of_vertices ):
            g.add_edge( i - 1, i )
        return g

    # OPERACJA DODATKOWA
    # Zwraca gwiazdę o podanej liczbie wierzchołków.
    @staticmethod
    def create_star( number_of_vertices ):
        g = Graph()
        for i in range( number_of_vertices ):
            g.add_vertex( i )
        for i in range( 1, number_of_vertices ):
            g.add_edge( 0, i )
        return g

    # OPERACJA DODATKOWA
    # Zwraca koło o podanej liczbie wierzchołków.
    @staticmethod
    def create_wheel( number_of_vertices ):
        g = Graph()
        for i in range( number_of_vertices ):
            g.add_vertex( i )
        for i in range( 0, number_of_vertices - 1 ):
            g.add_edge( i, number_of_vertices - 1 )
            g.add_edge( i, (i + 1) % (number_of_vertices - 1) )
        return g

    # OPERACJA DODATKOWA
    # Zwraca ciąg stopni grafu.
    def degree_sequence( self ):
        return sorted( (self.vertex_degree( v ) for v in self.__vertices), reverse = True )

    # OPERACJA DODATKOWA
    # Zwraca graf powstały wskutek ściągnięcia podanej krawędzi.
    def edge_contraction( self, u, v ):
        g = Graph()
        for v0 in self.__vertices:
            if v0 != u:
                g.add_vertex( v0 )
        for u0, v0 in self.__edges:
            if (u0, v0) != sorted( (u, v) ):
                g.add_edge( u0 if u0 != u else v, v0 if v0 != u else v )
        return g

    # OPERACJA DODATKOWA
    # Zwraca podgraf indukowany przez podany zbiór.
    def induced_subgraph( self, vertices ):
        g = Graph()
        for v in self.__vertices & vertices:
            g.add_vertex( v )
        for u, v in self.__edges:
            if u in vertices and v in vertices:
                g.add_edge( u, v )
        return g

    # OPERACJA DODATKOWA
    # Zwraca informację o tym, czy graf jest dwudzielny.
    def is_bipartite( self ):
        cs = {v: 0 for v in self.__vertices}
        def color_component( vertex, color ):
            cs[vertex] = color
            for v in self.__vertices:
                if self.is_edge( v, vertex ):
                    if cs[v] == color:
                        return False
                    if cs[v] == 0:
                        if not color_component( v, 3 - color ):
                            return False
            return True
        for v in self.__vertices:
            if cs[v] == 0:
                if not color_component( v, 1 ):
                    return False
        return True

    # OPERACJA DODATKOWA
    # Zwraca informację o tym, czy graf jest pełnym dwudzielnym.
    def is_complete_bipartite( self ):
        vs = self.vertex_neighbors( next( iter( self.__vertices ) ) )
        for u in self.__vertices:
            for v in self.__vertices:
                if u != v:
                    if self.is_edge( u, v ):
                        if (u in vs and v in vs) or (u not in vs and v not in vs):
                            return False
                    else:
                        if (u in vs and v not in vs) or (u not in vs and v in vs):
                            return False
        return True

    # OPERACJA DODATKOWA
    # Zwraca informację o tym, czy graf jest drzewem.
    def is_tree( self ):
        cs = set()
        def mark_component( vertex ):
            cs.add( vertex )
            for v in self.__vertices:
                if v not in cs and self.is_edge( vertex, v ):
                    mark_component( v )
        mark_component( next( iter( self.__vertices ) ) )
        return len( cs ) == self.number_of_vertices() and self.number_of_edges() == self.number_of_vertices() - 1

    # OPERACJA DODATKOWA
    # Zwraca liczbę cykli C3 zawartych w grafie.
    def number_of_triangles( self ):
        return sum( sum( self.is_edge( v0, u ) and self.is_edge( v1, u ) for u in self.__vertices ) for v0, v1 in self.__edges ) // 3

    # OPERACJA DODATKOWA
    # Zwraca graf, w którym wszystkie wierzchołki stopnia 2 zostały zastąpione krawędziami łączącymi ich sąsiadów.
    def smoothing( self ):
        g = Graph()
        for v in self.__vertices:
            g.add_vertex( v )
        for u, v in self.__edges:
            g.add_edge( u, v )
        vs = [v for v in g.vertices() if g.vertex_degree( v ) == 2]
        while vs:
            u0, u1 = g.vertex_neighbors( vs[0] )
            g.delete_vertex( vs[0] )
            g.add_edge( u0, u1 )
            vs = [v for v in g.vertices() if g.vertex_degree( v ) == 2]
        return g

    # OPERACJA DODATKOWA
    # Zwraca kwadrat grafu.
    def square( self ):
        g = Graph()
        for v in self.__vertices:
            g.add_vertex( v )
            for u in self.__vertices:
                if self.is_edge( u, v ):
                    g.add_edge( u, v )
                elif u != v:
                    for u0 in self.__vertices:
                        if self.is_edge( u0, v ) and self.is_edge( u0, u ):
                            g.add_edge( u, v )
        return g
