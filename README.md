# Spruce
A Fast, Space-saving Structure for In-Memory Dynamic Graph Storage

This is a demo of Spruce for ease of use, which contains the core data structure of Spruce and supports the following apis:

### Basic Interface
```C++
bool InsertEdge(SpruceUniverse &spruce, WeightedEdge edge);
```
InsertEdge function inserts a directed, weighted edge into spruce. If the connected verteices are not in it, they will be automatically inserted.
```C++
bool DeleteEdge(SpruceUniverse &spruce, uint64_t from_node_id, uint64_t to_node_id);
```
DeleteEdge function deletes a directed, weighted edge into spruce. If the edge's related vertices are not connected by any edges after deletion, they will be deleted as well.
```C++
bool UpdateEdge(SpruceUniverse &spruce, WeightedEdge edge);
```
UpdateEdge function tries to update an edge in the graph. If the edge already existed, then its weight will be updated. If the edge are not in the graph, then it will be inserted as a new edge. If the edge's weight < 0 (which is viewed as a delete_flag), the edge will be deleted from the graph.
```C++
bool get_neighbours(SpruceUniverse &spruce, uint64_t from_node_id,std::vector<WeightedOutEdge> &neighbours);
```
get_neighbours function gets all the adjacency edges of a vertex and return it using a vector.
```C++
uint64_t GetDegree(SpruceUniverse &spruce, uint64_t from_node_id);
```
Get a vertex's degree.

                               
