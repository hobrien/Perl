#!/Users/HeathOBrien/anaconda/bin/python



class Graph(object):
    """ A Python Class
A simple Python graph class, demonstrating the essential 
facts and functionalities of graphs.
    """

    def __init__(self, graph_dict={}):
        """ initializes a graph object """
        self.__graph_dict = graph_dict
        self.__edge_labels = {}
        self.__connected_verticies = set()

    def __getitem__(self, item):
        """returns list of vertices adjacent to item"""
        graph = self.__graph_dict
        labels = self.__edge_labels
        if len(item) == 1:
          return graph[item]
        else:
          return labels[item]

    def __str__(self):
        res = "vertices: "
        for k in self.__graph_dict:
            res += str(k) + " "
        res += "\nedges: "
        for edge in self.__generate_edges():
            res += str(edge) + " "
        return res

    def vertices(self):
        """ returns the vertices of a graph """
        return list(self.__graph_dict.keys())

    def connected_vertices(self):
        """ returns the vertices of a graph """
        return list(self.__connected_verticies)

    def edges(self):
        """ returns the edges of a graph """
        return self.__generate_edges()

    def children(self, vertex):
        return self.__generate_children(vertex)

    def add_vertex(self, vertex):
        """ If the vertex "vertex" is not in 
            self.__graph_dict, a key "vertex" with an empty
            list as a value is added to the dictionary. 
            Otherwise nothing has to be done. 
        """
        if vertex not in self.__graph_dict:
            self.__graph_dict[vertex] = []

    def add_edge(self, edge, label=''):
        """ assumes that edge is of type set, tuple or list; 
            between two vertices can be multiple edges! 
        """
        edge = set(edge)
        (vertex1, vertex2) = tuple(edge)
        if vertex1 in self.__graph_dict:
            self.__graph_dict[vertex1].append(vertex2)
        else:
            self.__graph_dict[vertex1] = [vertex2]
        if not vertex1 == vertex2:
            if vertex2 in self.__graph_dict:
                self.__graph_dict[vertex2].append(vertex1)
            else:
                self.__graph_dict[vertex2] = [vertex1]
        if label:
            self.__edge_labels[(vertex1, vertex2)] = label
            self.__edge_labels[(vertex2, vertex1)] = label
            
    def __generate_edges(self):
        """ A static method generating the edges of the 
            graph "graph". Edges are represented as sets 
            with one (a loop back to the vertex) or two 
            vertices 
        """
        edges = []
        for vertex in self.__graph_dict:
            for neighbour in self.__graph_dict[vertex]:
                if (vertex, neighbour) in self.__edge_labels:
                    label = self.__edge_labels[(vertex, neighbour)]
                else:
                    label =''
                if (label, {vertex, neighbour}) not in edges:
                    edges.append((label, {vertex, neighbour}))
        return edges

    def __generate_children(self, vertex):
        """ A static method generating a list of edge labels for the children
        of a node.  
        """
        edges = {}
        for neighbour in self.__graph_dict[vertex]:
             if neighbour > vertex and (vertex, neighbour) in self.__edge_labels:
                    edges[self.__edge_labels[(vertex, neighbour)]] = neighbour
        return edges

    def find_path(self, start_vertex, end_vertex, path=[]):
        """ find a path from start_vertex to end_vertex 
            in graph """
        graph = self.__graph_dict
        path = path + [start_vertex]
        if start_vertex == end_vertex:
            return path
        if start_vertex not in graph:
            return None
        for vertex in graph[start_vertex]:
            if vertex not in path:
                extended_path = self.find_path(vertex, 
                                               end_vertex, 
                                               path)
                if extended_path: 
                    return extended_path
        return None

    def find_all_paths(self, start_vertex, end_vertex, path=[]):
        """ find all paths from start_vertex to 
            end_vertex in graph """
        graph = self.__graph_dict 
        path = path + [start_vertex]
        if start_vertex == end_vertex:
            return [path]
        if start_vertex not in graph:
            return []
        paths = []
        for vertex in graph[start_vertex]:
            if vertex not in path:
                extended_paths = self.find_all_paths(vertex, 
                                                     end_vertex, 
                                                     path)
                for p in extended_paths: 
                    paths.append(p)
        return paths

    def vertex_degree(self, vertex):
        """ The degree of a vertex is the number of edges connecting
            it, i.e. the number of adjacent vertices. Loops are counted 
            double, i.e. every occurence of vertex in the list 
            of adjacent vertices. """ 
        adj_vertices =  self.__graph_dict[vertex]
        degree = len(adj_vertices) + adj_vertices.count(vertex)
        return degree
 
    def find_isolated_vertices(self):
        """ returns a list of isolated vertices. """
        graph = self.__graph_dict
        isolated = []
        for vertex in graph:
            if not graph[vertex]:
                isolated += [vertex]
        return isolated

    def find_leaf(self):
        """ returns a list of leaves (verticies of degree 1). """
        graph = self.__graph_dict
        leaves = []
        for vertex in graph:
            if self.vertex_degree(vertex) == 1:
                leaves += [vertex]
        return leaves

    def delta(self):
        """ the minimum degree of the vertices """
        min = 100000000
        for vertex in self.__graph_dict:
            vertex_degree = self.vertex_degree(vertex)
            if vertex_degree < min:
                min = vertex_degree
        return min
        
    def Delta(self):
        """ the maximum degree of the vertices """
        max = 0
        for vertex in self.__graph_dict:
            vertex_degree = self.vertex_degree(vertex)
            if vertex_degree > max:
                max = vertex_degree
        return max

    def degree_sequence(self):
        """ calculates the degree sequence """
        seq = []
        for vertex in self.__graph_dict:
            seq.append(self.vertex_degree(vertex))
        seq.sort(reverse=True)
        return tuple(seq)         
                              
    def density(self):
        """ method to calculate the density of a graph """
        g = self.__graph_dict
        V = len(g.keys())
        E = len(self.edges())
        return 2.0 * E / (V *(V - 1))

    def is_connected(self, 
                     start_vertex=None):
        """ determines if the graph is connected """
        gdict = self.__graph_dict        
        vertices = gdict.keys() 
        if not start_vertex:
            # choose a vertex from graph as a starting point
            start_vertex = vertices[0]
        self.__connected_verticies.add(start_vertex)
        #vertices_encountered.add(start_vertex)
        if len(self.__connected_verticies) != len(vertices):
        #if len(vertices_encountered) != len(vertices):
            for vertex in gdict[start_vertex]:
                if vertex not in self.__connected_verticies:
                #if vertex not in vertices_encountered:
                    if self.is_connected(vertex):
                    #if self.is_connected(vertices_encountered, vertex):
                        return True
        else:
            return True
        return False

    @staticmethod
    def erdoes_gallai(dsequence):
        """ Checks if the condition of the Erdoes-Gallai inequality 
            is fullfilled 
            dsequence has to be a valid degree sequence
        """
        if sum(dsequence) % 2:
            # sum of sequence is odd
            return False
        for k in range(1,len(dsequence) + 1):
            left = sum(dsequence[:k])
            right =  k * (k-1) + sum([min(x,k) for x in dsequence[k:]])
            if left > right:
                return False
        return True
        
class DiGraph(object):
    """As above, but for directed graphs. This means that reciprocal relationships aren't
  recorded. IE, if graph_dict['A'] = 'B', this does not mean that 'A' is in graph_dict['B']
    """
    def __init__(self, graph_dict={}):
        """ initializes a graph object """
        self.__graph_dict = graph_dict
        self.__edge_labels = {}
        self.__connected_verticies = set()

    def __getitem__(self, item):
        """returns list of vertices adjacent to item"""
        graph = self.__graph_dict
        labels = self.__edge_labels
        if len(item) == 1:
          return graph[item]
        else:
          return labels[item]

    def __str__(self):
        res = "vertices: "
        for k in self.__graph_dict:
            res += str(k) + " "
        res += "\nedges: "
        for edge in self.__generate_edges():
            res += str(edge) + " "
        return res

    def vertices(self):
        """ returns the vertices of a graph """
        return list(self.__graph_dict.keys())

    def connected_vertices(self):
        """ returns the vertices of a graph """
        return list(self.__connected_verticies)

    def edges(self):
        """ returns the edges of a graph """
        return self.__generate_edges()

    def children(self, vertex):
        return self.__generate_children(vertex)

    def add_vertex(self, vertex):
        """ If the vertex "vertex" is not in 
            self.__graph_dict, a key "vertex" with an empty
            list as a value is added to the dictionary. 
            Otherwise nothing has to be done. 
        """
        if vertex not in self.__graph_dict:
            self.__graph_dict[vertex] = []

    def add_edge(self, edge, label=''):
        """ assumes that edge is of type tuple or list; 
            between two vertices can be multiple edges! 
        """
        (vertex1, vertex2) = tuple(edge)
        if vertex1 in self.__graph_dict:
            self.__graph_dict[vertex1].append(vertex2)
        else:
            self.__graph_dict[vertex1] = [vertex2]
        if not vertex2 in self.__graph_dict:
            self.add_vertex(vertex2)
        if label:
            self.__edge_labels[(vertex1, vertex2)] = label
            
    def __generate_edges(self):
        """ A static method generating the edges of the 
            graph "graph". Edges are represented as sets 
            with one (a loop back to the vertex) or two 
            vertices 
        """
        edges = []
        for vertex in self.__graph_dict:
            for neighbour in self.__graph_dict[vertex]:
                if (vertex, neighbour) in self.__edge_labels:
                    label = self.__edge_labels[(vertex, neighbour)]
                else:
                    label =''
                if (label, {vertex, neighbour}) not in edges:
                    edges.append((label, {vertex, neighbour}))
        return edges

    def __generate_children(self, vertex):
        """ A static method generating a list of edge labels for the children
        of a node.  
        """
        edges = {}
        for neighbour in self.__graph_dict[vertex]:
             if neighbour > vertex and (vertex, neighbour) in self.__edge_labels:
                    edges[self.__edge_labels[(vertex, neighbour)]] = neighbour
        return edges

    def find_path(self, start_vertex, end_vertex, path=[]):
        """ find a path from start_vertex to end_vertex 
            in graph """
        graph = self.__graph_dict
        path = path + [start_vertex]
        if start_vertex == end_vertex:
            return path
        if start_vertex not in graph:
            return None
        for vertex in graph[start_vertex]:
            if vertex not in path:
                extended_path = self.find_path(vertex, 
                                               end_vertex, 
                                               path)
                if extended_path: 
                    return extended_path
        return None

    def find_all_paths(self, start_vertex, end_vertex, path=[]):
        """ find all paths from start_vertex to 
            end_vertex in graph """
        graph = self.__graph_dict 
        path = path + [start_vertex]
        if start_vertex == end_vertex:
            return [path]
        if start_vertex not in graph:
            return []
        paths = []
        for vertex in graph[start_vertex]:
            if vertex not in path:
                extended_paths = self.find_all_paths(vertex, 
                                                     end_vertex, 
                                                     path)
                for p in extended_paths: 
                    paths.append(p)
        return paths

    def vertex_degree(self, vertex):
        """ The degree of a vertex is the number of edges connecting
            it, i.e. the number of adjacent vertices. Loops are counted 
            double, i.e. every occurence of vertex in the list 
            of adjacent vertices. """ 
        adj_vertices =  self.__graph_dict[vertex]
        degree = len(adj_vertices) + adj_vertices.count(vertex)
        return degree
 
    def find_isolated_vertices(self):
        """ returns a list of isolated vertices. """
        graph = self.__graph_dict
        isolated = []
        for vertex in graph:
            if not graph[vertex]:
                isolated += [vertex]
        return isolated

    def find_leaf(self):
        """ returns a list of leaves (verticies of degree 1). """
        graph = self.__graph_dict
        leaves = []
        for vertex in graph:
            if self.vertex_degree(vertex) == 1:
                leaves += [vertex]
        return leaves

    def delta(self):
        """ the minimum degree of the vertices """
        min = 100000000
        for vertex in self.__graph_dict:
            vertex_degree = self.vertex_degree(vertex)
            if vertex_degree < min:
                min = vertex_degree
        return min
        
    def Delta(self):
        """ the maximum degree of the vertices """
        max = 0
        for vertex in self.__graph_dict:
            vertex_degree = self.vertex_degree(vertex)
            if vertex_degree > max:
                max = vertex_degree
        return max

    def degree_sequence(self):
        """ calculates the degree sequence """
        seq = []
        for vertex in self.__graph_dict:
            seq.append(self.vertex_degree(vertex))
        seq.sort(reverse=True)
        return tuple(seq)         
                              
    def density(self):
        """ method to calculate the density of a graph """
        g = self.__graph_dict
        V = len(g.keys())
        E = len(self.edges())
        return 2.0 * E / (V *(V - 1))

    def is_connected(self, 
                     start_vertex=None):
        """ determines if the graph is connected """
        gdict = self.__graph_dict        
        vertices = gdict.keys() 
        if not start_vertex:
            # choose a vertex from graph as a starting point
            start_vertex = vertices[0]
        self.__connected_verticies.add(start_vertex)
        #vertices_encountered.add(start_vertex)
        if len(self.__connected_verticies) != len(vertices):
        #if len(vertices_encountered) != len(vertices):
            for vertex in gdict[start_vertex]:
                if vertex not in self.__connected_verticies:
                #if vertex not in vertices_encountered:
                    if self.is_connected(vertex):
                    #if self.is_connected(vertices_encountered, vertex):
                        return True
        else:
            return True
        return False

    @staticmethod
    def erdoes_gallai(dsequence):
        """ Checks if the condition of the Erdoes-Gallai inequality 
            is fullfilled 
            dsequence has to be a valid degree sequence
        """
        if sum(dsequence) % 2:
            # sum of sequence is odd
            return False
        for k in range(1,len(dsequence) + 1):
            left = sum(dsequence[:k])
            right =  k * (k-1) + sum([min(x,k) for x in dsequence[k:]])
            if left > right:
                return False
        return True


if __name__ == "__main__":

    g = { "a" : ["d"],
          "b" : ["c"],
          "c" : ["b", "c", "d", "e"],
          "d" : ["a", "c"],
          "e" : ["c"],
          "f" : []
        }


    graph = Graph(g)

    print("Vertices of graph:")
    print(graph.vertices())

    print("Edges of graph:")
    print(graph.edges())

    print("Add vertex:")
    graph.add_vertex("z")

    print("Vertices of graph:")
    print(graph.vertices())
 
    print("Add an edge:")
    graph.add_edge({"a","z"})
    
    print("Vertices of graph:")
    print(graph.vertices())

    print("Edges of graph:")
    print(graph.edges())

    print('Adding an edge {"x","y"} with new vertices:')
    graph.add_edge({"x","y"})
    print("Vertices of graph:")
    print(graph.vertices())
    print("Edges of graph:")
    print(graph.edges())
    
    print "Graph: ", graph
    print "Path from a to e: ", graph.find_path('a','e')
    print "All paths from a to e: ", graph.find_all_paths('a','e')
    print "Degree of c: ", graph.vertex_degree('c')
    print "Isolated verticies: ", graph.find_isolated_vertices()
    print "Adjacent nodes of c: ", graph['c']
    print "delta: ", graph.delta()
    print "Delta: ", graph.Delta()
    print "Degree seq: ", graph.degree_sequence()
    print graph.erdoes_gallai(graph.degree_sequence())
    
    complete_graph = { 
        "a" : ["b","c"],
        "b" : ["a","c"],
        "c" : ["a","b"]
    }

    isolated_graph = { 
        "a" : [],
        "b" : [],
        "c" : []
    }


    print(graph)
    print(graph.is_connected())
    print(graph.density())

    print "%%%%%%%%%%%%%%%%%%%%%%%%%"
    graph = Graph(complete_graph)
    print "Complete Graph:\n", graph
    print(graph.is_connected())
    print graph.connected_vertices()
    print(graph.density())

    graph = Graph(isolated_graph)
    print "Isolated Graph:\n", graph
    print(graph.is_connected())
    print graph.connected_vertices()
    print(graph.density())
    extra_edges = 0
    while not graph.is_connected():
        for vertex in graph.vertices():
            if vertex not in graph.connected_vertices():
                graph.add_edge((vertex, graph.vertices()[0]))
                extra_edges += 1
                break
    print extra_edges
    
    graph = Graph()
    for x in range(1,10):
      graph.add_vertex(x)
    print graph
    print "attempt to make graph with labeled edges"
    graph = Graph(g)
    print graph
    graph.add_edge({"h","i"}, '1')
    print graph
    print g
