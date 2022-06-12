//#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "structmember.h"

typedef struct ListItemStruct {
    unsigned short value;
    struct ListItemStruct* nextItem;
} ListItem;

typedef struct {
    PyObject_HEAD
    unsigned short vertices;
    ListItem* neighbors[16];
} AdjacencyList;

static PyObject* AdjacencyList__new__(PyTypeObject* type, PyObject* args, PyObject* kdws){
    AdjacencyList* self;
    self = (AdjacencyList*)type->tp_alloc(type, 0);
    return (PyObject*) self;
};

static void AdjacencyList_dealloc(AdjacencyList* self){
    for(int i = 0; i < 16; i++){
        while(self->neighbors[i]){
            ListItem* toRemove = self->neighbors[i];
            self->neighbors[i] = toRemove->nextItem;
            free(toRemove);
        }
    }
    Py_TYPE(self)->tp_free((PyObject*) self);
}

int insert_edge(AdjacencyList* self, int v, int u){
    //first vertex
    if (self->neighbors[v] == NULL || self->neighbors[v]->value > u){
        ListItem* newEdge = (ListItem*) malloc(sizeof(ListItem*));
        if (!newEdge) return -1;
        newEdge->value = u;
        newEdge->nextItem = self->neighbors[v];
        self->neighbors[v] = newEdge;
    } else {
        ListItem* current = self->neighbors[v];
        while(current->nextItem != NULL && current->nextItem->value <= u) {
            current = current->nextItem;
        }
        // ignore if edge already exists
        if (u != current->value){
            ListItem* newEdge = (ListItem*) malloc(sizeof(ListItem*));
            if (!newEdge) return -1;
            newEdge->value = u;
            newEdge->nextItem = current->nextItem;
            current->nextItem = newEdge;
        }
    }
    //second vertex
    if (self->neighbors[u] == NULL || self->neighbors[u]->value > v){
        ListItem* newEdge = (ListItem*) malloc(sizeof(ListItem*));
        if (!newEdge) return -1;
        newEdge->value = v;
        newEdge->nextItem = self->neighbors[u];
        self->neighbors[u] = newEdge;
    } else {
        ListItem* current = self->neighbors[u];
        while(current->nextItem != NULL && current->nextItem->value <= v) {
            current = current->nextItem;
        }
        // ignore if edge already exists
        if (v != current->value){
            ListItem* newEdge = (ListItem*) malloc(sizeof(ListItem*));
            if (!newEdge) return -1;
            newEdge->value = v;
            newEdge->nextItem = current->nextItem;
            current->nextItem = newEdge;
        }
    }
    return 0;
}

static int AdjacencyList__init__(AdjacencyList* self, PyObject* args, PyObject* kwds){
    char* text = "?";
    if(!PyArg_ParseTuple(args, "|s", &text)){
        return NULL;
    }
    int length = (int) text[0] - 63; // 63 is "?"
    
    self->vertices = 0;
    unsigned short position = 0x8000;
    for (int i = 0; i < length; i++){
        self->vertices = self->vertices | position;
        position = (position >> 1);
    }

    for (int i = 0; i < 16; i++){
        self->neighbors[i] = NULL;
    }
    int k = 0;
    int i = 1;
    int c;
    for (int v = 1; v < length; v++){
        for (int u = 0; u < v; u++){
            if (k == 0){
                c = (int) text[i] - 63;
                i++;
                k = 6;
            }
            k--;
            if ((c & (1 << k)) != 0){
                insert_edge(self, v, u);
            }
        }
    }

    return 0;
}

static PyMemberDef AdjacencyList_Members[] = {
    {"vertices", T_USHORT, offsetof(AdjacencyList, vertices), 0, "Vertices in the graph"},
    {"neighbors", T_UINT, offsetof(AdjacencyList, neighbors), 0, "Neighbours of vertices"},
    {NULL}
};

static PyObject* number_of_vertices(AdjacencyList* self){
    int counter = 0;
    unsigned short position = 0x8000;
    for (int i = 0; i < 16; i++){
        if(self->vertices & position) 
            counter++;
        position = position >> 1;
    }
    return Py_BuildValue("i", counter);
}

static PyObject* vertices(AdjacencyList* self){
    PyObject* py_set = PySet_New(NULL);
    unsigned short position = 0x8000;
    for (int i = 0; i < 16; i++){
        if((self->vertices & position)){
            PyObject* py_int = Py_BuildValue("i", i);
            PySet_Add(py_set, py_int);
            Py_DECREF(py_int);
        }
        position = position >> 1;
    }
    return py_set;
}

static PyObject* vertex_degree(AdjacencyList* self, PyObject* vertex){
    int counter = 0;
    int v;

    if(!PyArg_ParseTuple(vertex, "i", &v)) 
        return NULL;
    
    ListItem* current = self->neighbors[v];
    while(current){
        counter++;
        current = current->nextItem;
    }

    return Py_BuildValue("i", counter);
}

static PyObject* vertex_neighbors(AdjacencyList* self, PyObject* vertex){
    int v;
    if (!PyArg_ParseTuple(vertex, "i", &v)) 
        return NULL;

    PyObject* py_set = PySet_New(NULL);
    ListItem* next = self->neighbors[v];
    while(next){
        PyObject* py_int = Py_BuildValue("i", next->value);
        PySet_Add(py_set, py_int);
        Py_DECREF(py_int);
        next = next->nextItem;
    }

    return py_set;
}

static PyObject* add_vertex(AdjacencyList* self, PyObject* vertex){
    int v;
    if (!PyArg_ParseTuple(vertex, "i", &v)) return NULL;

    unsigned short position;
    position = 0x8000 >> v;
    self->vertices |= position;

    Py_RETURN_NONE;
}

static PyObject* delete_vertex(AdjacencyList* self, PyObject* vertex){
    int v;
    if (!PyArg_ParseTuple(vertex, "i", &v)) 
        return NULL;

    unsigned short number = ~(0x8000 >> v);
    self->vertices &= number; // xor would toggle for already unset bit
    ListItem* next = self->neighbors[v];
    if(next){ // ignore if already null
        self->neighbors[v] = NULL;
        while(next){
            ListItem* nextItem = next->nextItem;
            free(next);
            next = nextItem;
        }
    }
    for(int i = 0; i < 16; i++){
        if(i != v){
            ListItem* current = self->neighbors[i];
            if(current && current->value == v){
                self->neighbors[i] = current->nextItem;
                free(current);
            } else {
                while(current) {
                    if (current->nextItem) {
                        if(current->nextItem->value == v){
                            ListItem* del = current->nextItem;
                            current->nextItem = del->nextItem;
                            free(del);
                        }
                    }
                    current = current->nextItem;
                }
            }
        }
    }

    Py_RETURN_NONE;
}

static PyObject* number_of_edges(AdjacencyList* self){
    int counter = 0;
    for (int i = 0; i < 16; i++){
        ListItem* current = self->neighbors[i];
        while(current){
            counter++;
            current = current->nextItem;
        }
    }

    return Py_BuildValue("i", counter/2);
}

static PyObject* edges(AdjacencyList* self){
    PyObject* py_set = PySet_New(NULL);
    for(int i = 0; i < 16; i++){
        ListItem* current = self->neighbors[i];
        while(current){
            PyObject* py_tuple;
            if (i < current->value)
                py_tuple = Py_BuildValue("(ii)", i, current->value);
            else
                py_tuple = Py_BuildValue("(ii)", current->value, i);
            
            PySet_Add(py_set, py_tuple);
            Py_DECREF(py_tuple);

            current = current->nextItem;
        }
    }

    return py_set;
}

static PyObject* is_edge(AdjacencyList* self, PyObject* uv){
    int u;
    int v;
    if (!PyArg_ParseTuple(uv, "ii", &u, &v)) 
        return NULL;

    ListItem* current = self->neighbors[u];
    while(current){
        if(current->value == v){
            Py_RETURN_TRUE;
        }
        current = current->nextItem;
    }

    Py_RETURN_FALSE;
}

static PyObject* add_edge(AdjacencyList* self, PyObject* uv){
    int u;
    int v;
    if (!PyArg_ParseTuple(uv, "ii", &u, &v))
        return NULL;

    if(u == v)
        Py_RETURN_NONE;

    unsigned short position;
    position = 0x8000;
    // at least one of them isn't a vertice
    if(((position >> u) & self->vertices) == 0 || ((position >> v) & self->vertices) == 0) 
        return NULL;

    int res = insert_edge(self, u, v);
    if (res == -1)
        return NULL;

    Py_RETURN_NONE;
}

static PyObject* delete_edge(AdjacencyList* self, PyObject* uv){
    int u;
    int v;
    if (!PyArg_ParseTuple(uv, "ii", &u, &v))
        return NULL;

    if (u == v)
        Py_RETURN_NONE;

    unsigned short position;
    position = 0x8000;
    // return NULL if one of them isn't a vertice
    if(((position >> u) & self->vertices) == 0 || ((position >> v) & self->vertices) == 0)
        return NULL;

    ListItem* current = self->neighbors[u];
    if(current->value == v){
        current = current->nextItem;
        free(self->neighbors[u]);
        self->neighbors[u] = current;
    } else{
        while(current){
            if(current->nextItem){
                if(current->nextItem->value == v){
                    ListItem* toRemove = current->nextItem;
                    current->nextItem = toRemove->nextItem;
                    free(toRemove);
                }
            }
            current = current->nextItem;
        }
    }

    current = self->neighbors[v];
    if(current->value == u){
        current = current->nextItem;
        free(self->neighbors[v]);
        self->neighbors[v] = current;
    } else{
        while(current){
            if(current->nextItem){
                if(current->nextItem->value == u){
                    ListItem* toRemove = current->nextItem;
                    current->nextItem = toRemove->nextItem;
                    free(toRemove);
                }
            }
            current = current->nextItem;
        }
    }

    Py_RETURN_NONE;
}

int* v2deg(AdjacencyList* g, int* vs){
    // vertex degree calculation
    int counter = 0;
    int degree[16];
    int i;
    for(i=0;i<16;i++){
        // check if the i vertex doesn't have neighbors or doesn't exist
        if(g->neighbors[i] == NULL){
            degree[i]=0;
            continue;
        }
        //calculate the degree of vertices
        ListItem* current = g->neighbors[i];
        degree[i] = 1;
        while(current->nextItem != NULL) {
            current = current->nextItem;
            degree[i]+=1;
        }
        if(degree[i] == 2){
            counter += 1;
        }
    }
    //create dynamic array for 2nd degree vertices
    if (vs == NULL && counter != 0)
        vs = malloc(counter * sizeof(int));
    else if (counter == 0){
        free(vs);
        return NULL;
    }
    else
    {   
        free(vs);
        vs = malloc(counter * sizeof(int));
    }
    counter = 0;
    //find 2nd degree vertices and store them in an array
    for(i=0;i<16;i++){
        if(degree[i] == 2){
            vs[counter]=i;
            counter += 1;
            }
    }
    return vs;
}

static PyObject* smoothing(AdjacencyList* self){
    AdjacencyList* g;
    int i;
    PyTypeObject* type = Py_TYPE(self);
    g = (AdjacencyList*)type->tp_alloc(type, 0);

    //clone vertices to g
    g->vertices = self->vertices;
    for (int i = 0; i < 16; i++){
        g->neighbors[i] = NULL;
    }

    //clone edges to g
    for(i=0;i<16;i++){
        ListItem* current = self->neighbors[i];
        //if vertex doesn't have any neighbors
        if(current == NULL){
            g->neighbors[i] = NULL;
        }
        //if vertex has 1 or more neighbors
        else
        {
            ListItem* prevItem = (ListItem*) malloc(sizeof(ListItem*));
            prevItem->value = current->value;
            g->neighbors[i] = prevItem;
            if(current->nextItem == NULL)
            {
                g->neighbors[i]->nextItem = NULL;
                continue;
            }
            current = current->nextItem;
            while(current != NULL){
                ListItem* newItem = (ListItem*) malloc(sizeof(ListItem*));
                newItem->value = current->value;
                prevItem->nextItem = newItem;
                prevItem = newItem;
                current = current->nextItem;
            }
            prevItem->nextItem = NULL;
        }
    }

    //create dynamic array for 2nd degree vertices
    int* vs = NULL;
    vs = v2deg(g, vs);

    while(vs){
        int vert = vs[0];
        //find 2nd degree vertices
        ListItem* next = g->neighbors[vert];
        short u0 = next->value;
        short u1 = next->nextItem->value;
        
        //delete the vertex from the list
        unsigned short position = ~(0x8000 >> vert);
        g->vertices &= position;
        if(next){ // ignore if already null
            g->neighbors[vert] = NULL;
            while(next){
                ListItem* nextItem = next->nextItem;
                free(next);
                next = nextItem;
            }
        }
        for(int i = 0; i < 16; i++){
            if(i != vert){
                ListItem* current = g->neighbors[i];
                if(current && current->value == vert){
                    g->neighbors[i] = current->nextItem;
                    free(current);
                } else {
                    while(current) {
                        if (current->nextItem) {
                            if(current->nextItem->value == vert){
                                ListItem* toRemove = current->nextItem;
                                current->nextItem = toRemove->nextItem;
                                free(toRemove);
                            }
                        }
                        current = current->nextItem;
                    }
                }
            }
            //insert an edge connecting the deleted vertex's neighbors
            insert_edge(g, u0, u1);

            //update the dynamic array with 2nd degree vertices
            vs = v2deg(g, vs);
            }
        }

    return (PyObject*)g;
}

static PyMethodDef AdjacencyList_Methods[] = {
    {"number_of_vertices", (PyCFunction) number_of_vertices, METH_NOARGS, "Returns number of vertices"},
    {"vertices", (PyCFunction) vertices, METH_NOARGS, "Returns set of graph's vertices"},
    {"vertex_degree", (PyCFunction) vertex_degree, METH_VARARGS, "Returns degree of given vertex"},
    {"vertex_neighbors", (PyCFunction) vertex_neighbors, METH_VARARGS, "Returns neighbors of given vertex"},
    {"add_vertex", (PyCFunction) add_vertex, METH_VARARGS, "Adds vertex"},
    {"delete_vertex", (PyCFunction) delete_vertex, METH_VARARGS, "Removes vertex and all edges to it"},
    {"number_of_edges", (PyCFunction) number_of_edges, METH_NOARGS, "Returns number of edges"},
    {"edges", (PyCFunction) edges, METH_NOARGS, "Returns set of graph's edges"},
    {"is_edge", (PyCFunction) is_edge, METH_VARARGS, "Checks if an edge exists"},
    {"add_edge", (PyCFunction) add_edge, METH_VARARGS, "Adds edge to the graph"},
    {"delete_edge", (PyCFunction) delete_edge, METH_VARARGS, "Removes edge from graph"},
    {"smoothing", (PyCFunction) smoothing, METH_NOARGS, "Returns a smooth graph (2nd degree vertices replaced with edges connecting their neigbors)"},
    {NULL, NULL, 0, NULL}
};

static PyTypeObject AdjacencyListType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "adjacencylist.AdjacencyList",
    .tp_basicsize = sizeof(AdjacencyList),
    .tp_itemsize = 0,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_doc = PyDoc_STR("Simple graph stored as an Adjacency List"),
    .tp_new = AdjacencyList__new__,
    .tp_init = (initproc) AdjacencyList__init__,
    .tp_dealloc = (destructor) AdjacencyList_dealloc,
    .tp_members = AdjacencyList_Members,
    .tp_methods = AdjacencyList_Methods,
};


static PyModuleDef AdjacencyList_Module = {
    PyModuleDef_HEAD_INIT,
    "adjacencylist",
    "Module for simple graph functionality",
    -1,
};

PyMODINIT_FUNC PyInit_simple_graphs(void){
    PyObject* module;
    if (PyType_Ready(&AdjacencyListType) < 0) 
        return NULL;

    module = PyModule_Create(&AdjacencyList_Module);
    if (module == NULL)
        return NULL;

    Py_INCREF(&AdjacencyListType);
    if(PyModule_AddObject(module, "AdjacencyList", (PyObject*) &AdjacencyListType) < 0) {
        Py_DECREF(&AdjacencyListType);
        Py_DECREF(module);
    }
    return module;
}