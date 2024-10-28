# Spruce

This is a demo of [Spruce](https://dl.acm.org/doi/10.1145/3639282) for ease of use, which contains the core data structure of Spruce and supports the following APIs:

### Required Dependencies

- [CMake](https://cmake.org/)
- [Junction Hash](https://github.com/preshing/junction)

### Compliation

```shell
cmake -DCMAKE_BUILD_TYPE=Release -G "CodeBlocks - Unix Makefiles" -S path/to/Spruce -B path/to/output/directory
```

### Basic Interface
```C++
bool InsertEdge(SpruceUniverse &spruce, WeightedEdge edge);
```
InsertEdge function inserts a directed, weighted edge into Spruce. If the connected vertices are not in it, they will be automatically inserted.
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
get_neighbours function gets all the adjacency edges of a vertex and returns it using a vector.
```C++
uint64_t GetDegree(SpruceUniverse &spruce, uint64_t from_node_id);
```
The getDegree function gets a vertex's degree.

### Notes
For experiments with GFE-Driver, please refer to [gfe_driver_spruce](https://github.com/Stardust-SJF/gfe_driver).
