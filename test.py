#!/usr/bin/env python3

import graphs
import sys
import time
import tqdm
import zipfile


# Tester struktur grafowych. Sposoby użycia:
#
#   (1) ./test.py [-h|--help]
#
#   Wyświetla krótki opis i kończy działanie.
#
#   (2) ./test.py -t [ST|ST OP]
#
#   Użyte bez parametrów ST i OP wykonuje autotest testera. Użyte z parametrem ST sprawdza, czy klasa ST prawidłowo implementuje wszystkie operacje podstawowe,
#   a z parametrami ST i OP sprawdza, czy klasa ST prawidłowo implementuje operację dodatkową OP. Dopuszczalne wartości parametru ST można odczytać z tablicy
#   GRAPH_STRUCTURES, a parametru OP ze słownika GRAPH_OPERATIONS.      
#
# UWAGA! Tester korzysta z modułu tqdm, który nie jest częścią biblioteki standardowej Pythona i musi być osobno zainstalowany. Opis procesu instalacji można
# znaleźć w dokumentacji tego modułu, dostępnej na stronie https://github.com/tqdm/tqdm.


# Drukuje wskazany komunikat błędu i kończy działanie testera.
def print_error_and_quit( message ):
    print( f"{sys.argv[0]}: {message}" )
    exit( 0 )

# Test operacji wyznaczania dopełnienia grafu (complement).
def test_of_operation_complement( structure, g6_sequence ):
    g6_sequence = tqdm.tqdm( g6_sequence )
    try:
        for g6 in g6_sequence:
            g1, g2 = graphs.Graph( g6 ), structure( g6 )
            if g1.complement() != g2.complement():
                g6_sequence.close()
                print_error_and_quit( f"błędny wynik operacji dla grafu {g6}" )
    except Exception as e:
        print_error_and_quit( f"podczas testu na grafie {g6} wystąpił wyjątek {e}" )

# Test operacji wyznaczania liczby składowych spójności grafu (connected_components).
def test_of_operation_connected_components( structure, g6_sequence ):
    g6_sequence = tqdm.tqdm( g6_sequence )
    try:
        for g6 in g6_sequence:
            g1, g2 = graphs.Graph( g6 ), structure( g6 )
            if g1.connected_components() != g2.connected_components():
                g6_sequence.close()
                print_error_and_quit( f"błędny wynik operacji dla grafu {g6}" )
    except Exception as e:
        print_error_and_quit( f"podczas testu na grafie {g6} wystąpił wyjątek {e}" )

# Test operacji konstrukcji pełnych grafów dwudzielnych (create_complete_bipartite).
def test_of_operation_create_complete_bipartite( structure, g6_sequence ):
    for n in range( 2, 17 ):
        for n1 in range( 1, n ):
            try:
                g1, g2 = graphs.Graph.create_complete_bipartite( n1, n - n1 ) , structure.create_complete_bipartite( n1, n - n1 )
                if g1 != g2:
                    print_error_and_quit( f"błędny wynik operacji dla parametrów {n1}, {n - n1}" )
            except Exception as e:
                print_error_and_quit( f"podczas testu parametrów {n1}, {n - n1} wystąpił wyjątek {e}" )
            try:
                g1, g2 = graphs.Graph.create_complete_bipartite( n - n1, n1 ) , structure.create_complete_bipartite( n - n1, n1 )
                if g1 != g2:
                    print_error_and_quit( f"błędny wynik operacji dla parametrów {n - n1}, {n1}" )
            except Exception as e:
                print_error_and_quit( f"podczas testu parametrów {n - n1}, {n1} wystąpił wyjątek {e}" )

# Test operacji konstrukcji cykli (create_cycle).
def test_of_operation_create_cycle( structure, g6_sequence ):
    for n in range( 3, 17 ):
        try:
            g1, g2 = graphs.Graph.create_cycle( n ) , structure.create_cycle( n )
            if g1 != g2:
                print_error_and_quit( f"błędny wynik operacji dla parametru {n}" )
        except Exception as e:
            print_error_and_quit( f"podczas testu parametru {n} wystąpił wyjątek {e}" )

# Test operacji konstrukcji ścieżek (create_path).
def test_of_operation_create_path( structure, g6_sequence ):
    for n in range( 1, 17 ):
        try:
            g1, g2 = graphs.Graph.create_path( n ) , structure.create_path( n )
            if g1 != g2:
                print_error_and_quit( f"błędny wynik operacji dla parametru {n}" )
        except Exception as e:
            print_error_and_quit( f"podczas testu parametru {n} wystąpił wyjątek {e}" )

# Test operacji konstrukcji ścieżek (create_star).
def test_of_operation_create_star( structure, g6_sequence ):
    for n in range( 1, 17 ):
        try:
            g1, g2 = graphs.Graph.create_star( n ) , structure.create_star( n )
            if g1 != g2:
                print_error_and_quit( f"błędny wynik operacji dla parametru {n}" )
        except Exception as e:
            print_error_and_quit( f"podczas testu parametru {n} wystąpił wyjątek {e}" )

# Test operacji konstrukcji kół (create_wheel).
def test_of_operation_create_wheel( structure, g6_sequence ):
    for n in range( 4, 17 ):
        try:
            g1, g2 = graphs.Graph.create_wheel( n ) , structure.create_wheel( n )
            if g1 != g2:
                print_error_and_quit( f"błędny wynik operacji dla parametru {n}" )
        except Exception as e:
            print_error_and_quit( f"podczas testu parametru {n} wystąpił wyjątek {e}" )

# Test operacji wyznaczania ciągu stopni grafu (degree_sequence).
def test_of_operation_degree_sequence( structure, g6_sequence ):
    g6_sequence = tqdm.tqdm( g6_sequence )
    try:
        for g6 in g6_sequence:
            g1, g2 = graphs.Graph( g6 ), structure( g6 )
            if g1.degree_sequence() != g2.degree_sequence():
                g6_sequence.close()
                print_error_and_quit( f"błędny wynik operacji dla grafu {g6}" )
    except Exception as e:
        print_error_and_quit( f"podczas testu na grafie {g6} wystąpił wyjątek {e}" )

# Test operacji ściągania krawędzi (edge_contraction).
def test_of_operation_edge_contraction( structure, g6_sequence ):
    g6_sequence = tqdm.tqdm( g6_sequence )
    try:
        for g6 in g6_sequence:
            g1, g2 = graphs.Graph( g6 ), structure( g6 )
            for e1 in g1.edges():
                if g1.edge_contraction( *e1 ) != g2.edge_contraction( *e1 ):
                    g6_sequence.close()
                    print_error_and_quit( f"błędny wynik operacji dla grafu {g6} i krawędzi {e1}" )
    except Exception as e:
        print_error_and_quit( f"podczas testu na grafie {g6} i krawędzi {e1} wystąpił wyjątek {e}" )

# Test operacji wyznaczania podgrafu indukowanego (induced_subgraph).
def test_of_operation_induced_subgraph( structure, g6_sequence ):
    g6_sequence = tqdm.tqdm( g6_sequence )
    try:
        for g6 in g6_sequence:
            g1, g2 = graphs.Graph( g6 ), structure( g6 )
            s = g1.vertices().copy()
            for u in g1.vertices():
                s.discard( u )
                if g1.induced_subgraph( s ) != g2.induced_subgraph( s ):
                    g6_sequence.close()
                    print_error_and_quit( f"błędny wynik operacji dla grafu {g6} i zbioru {s}" )
    except Exception as e:
        print_error_and_quit( f"podczas testu na grafie {g6} i zbiorze {s} wystąpił wyjątek {e}" )

# Test operacji sprawdzania dwudzielności grafu (is_bipartite).
def test_of_operation_is_bipartite( structure, g6_sequence ):
    g6_sequence = tqdm.tqdm( g6_sequence )
    try:
        for g6 in g6_sequence:
            g1, g2 = graphs.Graph( g6 ), structure( g6 )
            if g1.is_bipartite() != g2.is_bipartite():
                g6_sequence.close()
                print_error_and_quit( f"błędny wynik operacji dla grafu {g6}" )
    except Exception as e:
        print_error_and_quit( f"podczas testu na grafie {g6} wystąpił wyjątek {e}" )

# Test operacji sprawdzania, czy graf jest pełnym dwudzielnym (is_complete_bipartite).
def test_of_operation_is_complete_bipartite( structure, g6_sequence ):
    g6_sequence = tqdm.tqdm( g6_sequence )
    try:
        for g6 in g6_sequence:
            g1, g2 = graphs.Graph( g6 ), structure( g6 )
            if g1.is_complete_bipartite() != g2.is_complete_bipartite():
                g6_sequence.close()
                print_error_and_quit( f"błędny wynik operacji dla grafu {g6}" )
    except Exception as e:
        print_error_and_quit( f"podczas testu na grafie {g6} wystąpił wyjątek {e}" )

# Test operacji sprawdzania, czy graf jest drzewem (is_tree).
def test_of_operation_is_tree( structure, g6_sequence ):
    g6_sequence = tqdm.tqdm( g6_sequence )
    try:
        for g6 in g6_sequence:
            g1, g2 = graphs.Graph( g6 ), structure( g6 )
            if g1.is_tree() != g2.is_tree():
                g6_sequence.close()
                print_error_and_quit( f"błędny wynik operacji dla grafu {g6}" )
    except Exception as e:
        print_error_and_quit( f"podczas testu na grafie {g6} wystąpił wyjątek {e}" )

# Test operacji wyznaczania liczby trójkątów w grafie (number_of_triangles).
def test_of_operation_number_of_triangles( structure, g6_sequence ):
    g6_sequence = tqdm.tqdm( g6_sequence )
    try:
        for g6 in g6_sequence:
            g1, g2 = graphs.Graph( g6 ), structure( g6 )
            if g1.number_of_triangles() != g2.number_of_triangles():
                g6_sequence.close()
                print_error_and_quit( f"błędny wynik operacji dla grafu {g6}" )
    except Exception as e:
        print_error_and_quit( f"podczas testu na grafie {g6} wystąpił wyjątek {e}" )

# Test operacji wygładzania grafu (smoothing).
def test_of_operation_smoothing( structure, g6_sequence ):
    g6_sequence = tqdm.tqdm( g6_sequence )
    try:
        for g6 in g6_sequence:
            g1, g2 = graphs.Graph( g6 ), structure( g6 )
            # print( f"{g6}" )
            if g1.smoothing() != g2.smoothing():
                g6_sequence.close()
                print_error_and_quit( f"błędny wynik operacji dla grafu {g6}" )
    except Exception as e:
        print_error_and_quit( f"podczas testu na grafie {g6} wystąpił wyjątek {e}" )

# Test operacji wyznaczania kwadratu grafu (square).
def test_of_operation_square( structure, g6_sequence ):
    g6_sequence = tqdm.tqdm( g6_sequence )
    try:
        for g6 in g6_sequence:
            g1, g2 = graphs.Graph( g6 ), structure( g6 )
            if g1.square() != g2.square():
                g6_sequence.close()
                print_error_and_quit( f"błędny wynik operacji dla grafu {g6}" )
    except Exception as e:
        print_error_and_quit( f"podczas testu na grafie {g6} wystąpił wyjątek {e}" )


# Struktury grafów, które tester akceptuje.
GRAPH_STRUCTURES = [ "AdjacencyMatrix", "AdjacencyList", "EdgesList", "IncidenceMatrix" ]

# Dodatkowe operacje, które tester akceptuje, wraz z przypisanymi im funkcjami testującymi.
GRAPH_OPERATIONS = {
    "complement":                test_of_operation_complement,
    "connected_components":      test_of_operation_connected_components,
    "create_complete_bipartite": test_of_operation_create_complete_bipartite,
    "create_cycle":              test_of_operation_create_cycle,
    "create_path":               test_of_operation_create_path,
    "create_star":               test_of_operation_create_star,
    "create_wheel":              test_of_operation_create_wheel,
    "degree_sequence":           test_of_operation_degree_sequence,
    "edge_contraction":          test_of_operation_edge_contraction,
    "induced_subgraph":          test_of_operation_induced_subgraph,
    "is_bipartite":              test_of_operation_is_bipartite,
    "is_complete_bipartite":     test_of_operation_is_complete_bipartite,
    "is_tree":                   test_of_operation_is_tree,
    "number_of_triangles":       test_of_operation_number_of_triangles,
    "smoothing":                 test_of_operation_smoothing,
    "square":                    test_of_operation_square,
}

# Wczytuje do pamięci zestawy grafów testowych.
def read_graphs():
    print( "Przygotowuję dane testowe. Proszę czekać." )
    with zipfile.ZipFile( "graphs.zip" ) as f:
        g6_sequence = "\n".join( f.open( name ).read().decode() for name in f.namelist() ).split()
    print( f"Przygotowano {len( g6_sequence )} grafów." )
    return g6_sequence

# Wykonuje test wskazanej struktury.
def test_structure( structure, g6_sequence ):
    g6_sequence = tqdm.tqdm( g6_sequence )
    for g6 in g6_sequence:

        # Test konwersji z formatu g6.
        try:
            g, h = graphs.Graph( g6 ), structure( g6 )
        except Exception as e:
            g6_sequence.close()
            print_error_and_quit( f"podczas tworzenia grafu {g6} wystąpił wyjątek {e}" )

        # Test funkcji number_of_vertices().
        try:
            if g.number_of_vertices() != h.number_of_vertices():
                g6_sequence.close()
                print_error_and_quit( f"błąd funkcji number_of_vertices() dla grafu {g6}" )
        except Exception as e:
            g6_sequence.close()
            print_error_and_quit( f"podczas testu funkcji number_of_vertices() dla grafu {g6} wystąpił wyjątek {e}" )

        # Test funkcji vertices().
        try:
            if g.vertices() != h.vertices():
                g6_sequence.close()
                print_error_and_quit( f"błąd funkcji vertices() dla grafu {g6}" )
        except Exception as e:
            g6_sequence.close()
            print_error_and_quit( f"podczas testu funkcji vertices() dla grafu {g6} wystąpił wyjątek {e}" )

        # Test funkcji number_of_edges().
        try:
            if g.number_of_edges() != h.number_of_edges():
                g6_sequence.close()
                print_error_and_quit( f"błąd funkcji number_of_edges() dla grafu {g6}" )
        except Exception as e:
            g6_sequence.close()
            print_error_and_quit( f"podczas testu funkcji number_of_edges() dla grafu {g6} wystąpił wyjątek {e}" )

        # Test funkcji edges().
        try:
            if g.edges() != h.edges():
                g6_sequence.close()
                print_error_and_quit( f"błąd funkcji edges() dla grafu {g6}" )
        except Exception as e:
            g6_sequence.close()
            print_error_and_quit( f"podczas testu funkcji edges() dla grafu {g6} wystąpił wyjątek {e}" )

        # Test funkcji is_edge().
        try:
            for u in range( 16 ):
                for v in range( u ):
                    if g.is_edge( u, v ) != h.is_edge( u, v ):
                        g6_sequence.close()
                        print_error_and_quit( f"błąd funkcji is_edge( {u}, {v} ) dla grafu {g6}" )
        except Exception as e:
            g6_sequence.close()
            print_error_and_quit( f"podczas testu funkcji is_edge() dla grafu {g6} wystąpił wyjątek {e}" )

        # Test funkcji vertex_degree().
        try:
            for u in range( g.number_of_vertices() ):
                if g.vertex_degree( u ) != h.vertex_degree( u ):
                    g6_sequence.close()
                    print_error_and_quit( f"błąd funkcji vertex_degree( {u} ) dla grafu {g6}" )
        except Exception as e:
            g6_sequence.close()
            print_error_and_quit( f"podczas testu funkcji vertex_degree() dla grafu {g6} wystąpił wyjątek {e}" )

        # Test funkcji vertex_neighbors().
        try:
            for u in range( g.number_of_vertices() ):
                if g.vertex_neighbors( u ) != h.vertex_neighbors( u ):
                    g6_sequence.close()
                    print_error_and_quit( f"błąd funkcji vertex_neighbors( {u} ) dla grafu {g6}" )
        except Exception as e:
            g6_sequence.close()
            print_error_and_quit( f"podczas testu funkcji vertex_neighbors() dla grafu {g6} wystąpił wyjątek {e}" )

        # Test funkcji add_vertex()/delete_vertex().
        try:
            for u in range( 16 ):
                if u in g.vertices():
                    g.delete_vertex( u )
                    h.delete_vertex( u )
                    if g != h:
                        g6_sequence.close()
                        print_error_and_quit( f"błąd funkcji delete_vertex( {u} ) dla grafu {g6}" )
                else:
                    g.add_vertex( u )
                    h.add_vertex( u )
                    if g != h:
                        g6_sequence.close()
                        print_error_and_quit( f"błąd funkcji add_vertex( {u} ) dla grafu {g6}" )
        except Exception as e:
            g6_sequence.close()
            print_error_and_quit( f"podczas testu funkcji add_vertex() / delete_vertex() dla grafu {g6} wystąpił wyjątek {e}" )

        # Test funkcji add_edge()/delete_edge().
        try:
            g, h = graphs.Graph( g6 ), structure( g6 )
            for u in range( g.number_of_vertices() ):
                for v in range( u ):
                    if g.is_edge( u, v ):
                        g.delete_edge( u, v )
                        h.delete_edge( u, v )
                        if g != h:
                            g6_sequence.close()
                            print_error_and_quit( f"błąd funkcji delete_edge( {u}, {v} ) dla grafu {g6}" )
                        g.add_edge( u, v )
                        h.add_edge( u, v )
                        if g != h:
                            g6_sequence.close()
                            print_error_and_quit( f"błąd funkcji add_edge( {u}, {v} ) dla grafu {g6}" )
                    else:
                        g.add_edge( u, v )
                        h.add_edge( u, v )
                        if g != h:
                            g6_sequence.close()
                            print_error_and_quit( f"błąd funkcji add_edge( {u}, {v} ) dla grafu {g6}" )
                        g.delete_edge( u, v )
                        h.delete_edge( u, v )
                        if g != h:
                            g6_sequence.close()
                            print_error_and_quit( f"błąd funkcji delete_edge( {u}, {v} ) dla grafu {g6}" )
        except Exception as e:
            g6_sequence.close()
            print_error_and_quit( f"podczas testu funkcji add_edge() / delete_edge() dla grafu {g6} wystąpił wyjątek {e}" )

if len( sys.argv ) == 1 or (len( sys.argv ) == 2 and sys.argv[1] in [ "-h", "--help" ]):
    print( f"Sposób użycia: {sys.argv[0]} [-h|--help|-t [STRUKTURA [OPERACJA]]]" )
    print( f"Wiecej informacji znajdziesz w pliku {sys.argv[0]}" )
    exit( 0 )
elif len( sys.argv ) == 2 and sys.argv[1] == "-t":
    g6_sequence = read_graphs()
    for o, f in GRAPH_OPERATIONS.items():
        print( f"Rozpoczynam test operacji {o}." )
        s, _, e = time.time(), f( graphs.Graph, g6_sequence ), time.time()
        print( f"Test operacji {o} zakończony sukcesem. Czas trwania testu to {e - s:.2f}s." )
    print( "Rozpoczynam test struktury graphs.Graph." )
    s, _, e = time.time(), test_structure( graphs.Graph, g6_sequence ), time.time()
    print( f"Test struktury graphs.Graph zakończony sukcesem. Czas trwania testu to {e - s:.2f}s." )
elif len( sys.argv ) == 3 and sys.argv[1] == "-t" and sys.argv[2] in GRAPH_STRUCTURES:
    try:
        import simple_graphs
        print( f"Rozpoczynam test struktury {sys.argv[2]}." )
        s, _, e = time.time(), test_structure( eval( f"simple_graphs.{sys.argv[2]}" ), read_graphs() ), time.time()
        print( f"Test struktury {sys.argv[2]} zakończony sukcesem. Czas trwania testu to {e - s:.2f}s." )
    except ImportError:
        print_error_and_quit( "moduł simple_graphs nie istnieje" )
    except AttributeError:
        print_error_and_quit( f"moduł simple_graphs nie zawiera struktury {sys.argv[2]}" )
elif len( sys.argv ) == 4 and sys.argv[1] == "-t" and sys.argv[2] in GRAPH_STRUCTURES and sys.argv[3] in GRAPH_OPERATIONS:
    try:
        import simple_graphs
        print( f"Rozpoczynam test operacji {sys.argv[3]} dla struktury {sys.argv[2]}." )
        s, _, e = time.time(), GRAPH_OPERATIONS[sys.argv[3]]( eval( f"simple_graphs.{sys.argv[2]}" ), read_graphs() ), time.time()
        print( f"Test operacji {sys.argv[3]} zakończony sukcesem. Czas trwania testu to {e - s:.2f}s." )
    except ImportError:
        print_error_and_quit( "moduł simple_graphs nie istnieje" )
    except AttributeError:
        print_error_and_quit( f"moduł simple_graphs nie zawiera struktury {sys.argv[2]} lub ta struktura nie zawiera funkcji {sys.argv[3]}" )
else:
    print_error_and_quit( f"błędne parametry, wykonaj {sys.argv[0]} -h aby uzyskać pomoc" )
