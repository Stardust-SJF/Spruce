//
// Created by sjf on 3/15/2022.
//

#include "graph_preprocess.h"

// The convert of EdgeList and G500 do not change the graph types. That is to say,
// directed/undirected graph should be recognized by graph systems.
// Generally speaking, graph500 by default consider files as undirected graphs (each edge stored once)
// while spruce can change the option "undirected flag".

#include <iostream>
#include <fstream>
#include <unordered_map>
#include <unordered_set>
#include <string>
#include <filesystem>
#include <iomanip>

bool ConvertToUndirectedGraph(std::string inputFileName) {
    std::ifstream input(inputFileName);
    if (!input.is_open()) {
        std::cerr << "Cannot open input file: " << inputFileName << std::endl;
        return false;
    }
    else {
        std::cerr << "Input file: " << inputFileName << std::endl;
    }

    std::unordered_map<uint32_t, std::unordered_set<uint32_t>> graph;

    int a, b;
    while (input >> a >> b) {
        if (graph[a].find(b) == graph[a].end() && graph[b].find(a) == graph[b].end()) {
            graph[a].insert(b);
        }
    }
    input.close();

    std::string outputFileName = inputFileName + ".el";

    if (std::filesystem::exists(outputFileName)) {
        std::cerr << "Output file: " << outputFileName << " already exists. Aborting conversion." << std::endl;
        return false;
    }

    std::ofstream output(outputFileName);
    if (!output.is_open()) {
        std::cerr << "Cannot create output file: " << outputFileName << std::endl;
        return false;
    }

    for (const auto& [vertex, neighbors] : graph) {
        for (const auto& neighbor : neighbors) {
            output << vertex << " " << neighbor << std::endl;
        }
    }

    output.close();

    return true;
}


int ConvertEdgeListToG500(char* input_file) {
    FILE *input, *output;
    packed_edge edge;
    char *output_filename;
    std::cout << "Converting EdgeList to g500...." << std::endl;
    input = fopen(input_file, "r");
    if (input == NULL) {
        perror("Error opening input file");
        return 1;
    }

    output_filename = strdup(input_file);
    if (output_filename == NULL) {
        perror("Error duplicating filename");
        return 1;
    }

    char *dot = strrchr(output_filename, '.');
    if (dot) {
        *dot = '\0';  // Null-terminate the string at the dot
    }

    output_filename = static_cast<char*>(realloc(output_filename, strlen(output_filename) + strlen(".g500") + 1));
    if (output_filename == NULL) {
        perror("Error allocating memory for output filename");
        return 1;
    }

    strcat(output_filename, ".g500");

    output = fopen(output_filename, "wb");
    if (output == NULL) {
        perror("Error opening output file");
        return 1;
    }

    char line[256];
    uint64_t v0, v1;
    char weight_str[50];
    while (fgets(line, sizeof line, input) != NULL) {
        if (strstr(input_file, "dota") != NULL) {
            // ignore weight
            if (sscanf(line, "%" SCNu64 " %" SCNu64 " %s", &v0, &v1, &weight_str) == 3) {
                edge.v0_low = v0 & 0xFFFFFFFF;
                edge.v1_low = v1 & 0xFFFFFFFF;
                edge.high = ((v0 >> 32) & 0xFFFF) | ((v1 >> 32) & 0xFFFF) << 16;
                fwrite(&edge, sizeof(packed_edge), 1, output);
            }
        }
        else {
            if (sscanf(line, "%" SCNu64 " %" SCNu64, &v0, &v1) == 2) {
                edge.v0_low = v0 & 0xFFFFFFFF;
                edge.v1_low = v1 & 0xFFFFFFFF;
                edge.high = ((v0 >> 32) & 0xFFFF) | ((v1 >> 32) & 0xFFFF) << 16;
                fwrite(&edge, sizeof(packed_edge), 1, output);
            }
        }
    }

    fclose(input);
    fclose(output);
    free(output_filename);

    return 0;
}

int ConvertG500toEdgeList(char* input_file) {
    FILE *input, *output;
    packed_edge edge;
    char *output_filename;

    std::cout << "Input file: " << input_file << std::endl;
    std::cout << "Converting g500 to edgelist ... " << std::endl;


    input = fopen(input_file, "rb");
    if (input == NULL) {
        perror("Error opening input file");
        return 1;
    }

    output_filename = strdup(input_file);
    if (output_filename == NULL) {
        perror("Error duplicating filename");
        return 1;
    }

    char *dot = strrchr(output_filename, '.');
    if (dot) {
        *dot = '\0';  // Null-terminate the string at the dot
    }

    output_filename = static_cast<char*>(realloc(output_filename, strlen(output_filename) + strlen(".el") + 1));
    if (output_filename == NULL) {
        perror("Error allocating memory for output filename");
        return 1;
    }

    strcat(output_filename, ".el");

    output = fopen(output_filename, "w");
    if (output == NULL) {
        perror("Error opening output file");
        return 1;
    }

    while (fread(&edge, sizeof(packed_edge), 1, input) == 1) {
        uint64_t v0 = ((uint64_t)(edge.high & 0xFFFF) << 32) | edge.v0_low;
        uint64_t v1 = ((uint64_t)(edge.high >> 16) << 32) | edge.v1_low;
        fprintf(output, "%" PRId64 " %" PRId64 "\n", v0, v1);
    }

    fclose(input);
    fclose(output);
    free(output_filename);

    return 0;
}

int ConvertToMetisFormat(const std::string &graph_path){
    //Only for undirected graph and node id start from 1; Else need to edit code

    // Init input and output stream
    std::ifstream input_file(graph_path, std::ios::in);
    std::string output_path = graph_path + ".mf";
    std::ofstream output_file;
    output_file.open(output_path, std::ios::out);
    if (!output_file|!input_file)  

        std::cout << "File Error!! " << std::endl;

    // Read Graph Prologues
    std::string temp_line;
    int nodes_num, edges_num, cnt = 0;
    while (cnt++ < 4) {
        getline(input_file, temp_line);
        if(cnt == 3){
            //Read number of nodes
            std::stringstream ss(temp_line);
            ss>>temp_line;
            ss>>temp_line;
            ss>>nodes_num;
            ss>>temp_line;
            ss>>edges_num;
        }
        std::cout << temp_line << std::endl;
    }


    // Read Graph Data(Sorted)
    int from_node_id, to_node_id, former_from_node_id;
    std::vector<std::vector<int>> node_edges(edges_num);    //Use edges_num for count, because some ids may not exist in txt
    int max_node_id = 0; //it decides the size of graph
    //Init the first line.
    node_edges[0].push_back(nodes_num);    //use last to_node_id as vertex num
    node_edges[0].push_back(edges_num);
    former_from_node_id = 0;   // former node id set to invalid for initialization
    while (input_file >> from_node_id) {
        input_file >> to_node_id;
        if(from_node_id < to_node_id){
            node_edges[from_node_id].push_back(to_node_id);
            node_edges[to_node_id].push_back(from_node_id);
        }
        max_node_id = (max_node_id < to_node_id) ? to_node_id : max_node_id;
    }
    node_edges[0][0] = max_node_id; //bug fixed: use last to_node_id as vertex num

    int i, j;
    for(i = 0; i < max_node_id + 1; i++){
        if(node_edges[i].size() != 0) {
            for (j = 0; j + 1 < node_edges[i].size(); j++) {
                output_file << std::to_string(node_edges[i][j]);
                output_file << " ";
            }
            output_file << std::to_string(node_edges[i][j]);
        }
        output_file << "\n";
    }

    input_file.close();
    output_file.close();

    return nodes_num;
}


// For GFE dataset, convert to Ligra symmetrical graphs
int ConvertToLigraFormat(const std::string &graph_path) {
    std::cout << "Converting edge list to Ligra Format..." << std::endl;
    // Init input and output stream
    std::ifstream input_file(graph_path, std::ios::in);
    std::string output_path = graph_path + ".li";
    std::ofstream output_file;
    std::string twitter_2010 = "twitter-2010";
    std::string hyperlink = "Hyperlink";
    std::string graph500 = "graph500";
    output_file.open(output_path, std::ios::out);
    if (!input_file)
        std::cout << "File Error!! " << std::endl;

    // Read Graph Prologues
    std::string temp_line;
    int nodes_num, edges_num, cnt = 0;
    if(graph_path.find(twitter_2010) != std::string::npos || graph_path.find(hyperlink) != std::string::npos ||
            graph_path.find(graph500) != std::string::npos || graph_path.find("GFE") != std::string::npos || graph_path.find("g500")!=std::string::npos){
        nodes_num = 41652230;
        edges_num = 2044000000; //only use edges_num that's enough
    }
    else {
        while (cnt++ < 4) {
            getline(input_file, temp_line);
            if (cnt == 3) {
                //Read number of nodes
                std::stringstream ss(temp_line);
                ss >> temp_line;
                ss >> temp_line;
                ss >> nodes_num;
                ss >> temp_line;
                ss >> edges_num;
            }
            std::cout << temp_line << std::endl;
        }
    }
    long int vertex_num, edge_num;
    if(graph_path.find("amazon")!=std::string::npos){
        vertex_num = 548551;
        edge_num = 925872;
    }
    else if(graph_path.find("LiveJournal")!=std::string::npos) {
        edge_num = 68993773;
        vertex_num = 4847571;
    }
    else if (graph_path.find("twitter")!=std::string::npos) {
        edge_num = 1468365182;
        vertex_num = 41652230;
    }
    else if (graph_path.find("Hyperlink")!=std::string::npos) {
        edge_num = 2043203933;
        vertex_num = 101717775;
    }
    else if (graph_path.find("500-26")!=std::string::npos) {
        vertex_num = 67108863;
        edge_num = 1051922853;
    }
    else if (graph_path.find("uniform-26")!=std::string::npos) {
        vertex_num = 49207468;
        edge_num = 1051922853;
    }
    else if (graph_path.find("uniform-24")!=std::string::npos) {
        vertex_num = 13306413;
        edge_num = 260379520;
    }
    else if (graph_path.find("500-24")!=std::string::npos) {
        vertex_num = 16777215;
        edge_num = 270000000;
    }
    else {
        vertex_num = 16777216;
        edge_num = edges_num;
    }


    // Read Graph Data(Sorted)
    uint32_t from_node_id, to_node_id, former_from_node_id;
    std::vector<std::vector<uint32_t>> node_edges(edge_num);    //Use edges_num for count, because some ids may not exist in txt
    int max_node_id = 0; //it decides the size of graph
    //Init the first line.
//    node_edges[0].push_back(nodes_num);    //use last to_node_id as vertex num
//    node_edges[0].push_back(edges_num);
    former_from_node_id = 0;   // former node id set to invalid for initialization
    long int temp_count = 0;
    double temp_double;
    if (graph_path.find("GFE") != std::string::npos || graph_path.find("g500")!=std::string::npos || graph_path.find("ungraph")!=std::string::npos) {
        while (input_file >> from_node_id) {
            input_file >> to_node_id;
            if (graph_path.find("dota") != std::string::npos) {
                // dota is a weighted graph, ignore graph weight;
                input_file >> temp_double;
            }
            node_edges[from_node_id].push_back(to_node_id);
            node_edges[to_node_id].push_back(from_node_id);
            max_node_id = (max_node_id < to_node_id) ? to_node_id : max_node_id;
            max_node_id = (max_node_id < from_node_id) ? from_node_id : max_node_id;
            temp_count += 2;
            if (temp_count % 1000000000 == 0) {
                std::cout << "reading" << temp_count / 1000000000 << "  billion edges" << std::endl;
            }
        }
    }
    else {
        while (input_file >> from_node_id) {
            input_file >> to_node_id;
            node_edges[from_node_id].push_back(to_node_id);
//            node_edges[to_node_id].push_back(from_node_id);
            max_node_id = (max_node_id < to_node_id) ? to_node_id : max_node_id;
            max_node_id = (max_node_id < from_node_id) ? from_node_id : max_node_id;
            temp_count++;
            if (temp_count % 1000000000 == 0) {
                std::cout << "reading" << temp_count / 1000000000 << "  billion edges" << std::endl;
            }
        }
    }
//    node_edges[0][0] = max_node_id; //bug fixed: use last to_node_id as vertex num

    //Write output file
    output_file << "AdjacencyGraph\n";
    output_file << max_node_id + 1 << "\n"; //begin from 0
    output_file << temp_count << "\n";
//    output_file.close();            // no need to close?
//    std::cout << node_edges[0][1] <<std::endl;

    //like CSR
    int i, j;
    long int offset = 0;
    //write offset
    //node : node_edges[max_node_id] may have edges
//    output_file.open(output_path, std::ios::app | std::ios::out);
    for(i = 0; i < max_node_id + 1; i++){
        output_file << std::to_string(offset) << "\n";
        offset += node_edges[i].size();
    }
    std::cout << offset <<std::endl;
    //write edge
    for(i = 0; i < max_node_id + 1; i++){
        if(node_edges[i].size() != 0) {
            for (j = 0; j < node_edges[i].size(); j++) {
                output_file << std::to_string(node_edges[i][j]);
                output_file << "\n";
            }
        }
    }

    input_file.close();
    output_file.close();

    return nodes_num;
}

int ConvertToGFEFormat(const std::string &graph_path) {
    std::ifstream input_file(graph_path, std::ios::in);
//    std::string output_path = graph_path + ".li";
    std::ofstream output_file;
    std::string twitter_2010 = "twitter-2010";
    std::string hyperlink = "Hyperlink";
    std::string graph500 = "graph500";
//    output_file.open(output_path, std::ios::out);
    if (!output_file|!input_file)  
        std::cout << "File Error!! " << std::endl;

    // Read Graph Prologues
    std::string temp_line;
    int nodes_num, edges_num, cnt = 0;
    if(graph_path.find(twitter_2010) != std::string::npos || graph_path.find(hyperlink) != std::string::npos ||
       graph_path.find(graph500) != std::string::npos || graph_path.find("GFE") != std::string::npos || graph_path.find("Orkut") != std::string::npos){
        nodes_num = 41652230;
        edges_num = 2044000000; //only use edges_num that's enough
    }
    else {
        while (cnt++ < 4) {
            getline(input_file, temp_line);
            if (cnt == 3) {
                //Read number of nodes
                std::stringstream ss(temp_line);
                ss >> temp_line;
                ss >> temp_line;
                ss >> nodes_num;
                ss >> temp_line;
                ss >> edges_num;
            }
            std::cout << temp_line << std::endl;
        }
    }
    long int vertex_num, edge_num;
    if(graph_path.find("amazon")!=std::string::npos){
        vertex_num = 548551;
        edge_num = 925872;
    }
    else if(graph_path.find("LiveJournal")!=std::string::npos) {
        edge_num = 68993773;
        vertex_num = 4847571;
    }
    else if (graph_path.find("twitter")!=std::string::npos) {
        edge_num = 1468365182;
        vertex_num = 41652230;
    }
    else if (graph_path.find("Hyperlink")!=std::string::npos) {
        edge_num = 2043203933;
        vertex_num = 101717775;
    }
    else if (graph_path.find("500-26")!=std::string::npos) {
        vertex_num = 67108863;
        edge_num = 1051922853;
    }
    else if (graph_path.find("uniform-26")!=std::string::npos) {
        vertex_num = 49207468;
        edge_num = 1051922853;
    }
    else if (graph_path.find("uniform-24")!=std::string::npos) {
        vertex_num = 13306413;
        edge_num = 260379520;
    }
    else if (graph_path.find("500-24")!=std::string::npos) {
        vertex_num = 16777215;
        edge_num = 260379520;
    }
    else {
        edge_num = 2044000000;
        vertex_num = 16777216;
    }


    // Read Graph Data(Sorted)
    uint32_t from_node_id, to_node_id, former_from_node_id;
    std::vector<std::vector<uint32_t>> node_edges(edge_num);    //Use edges_num for count, because some ids may not exist in txt
    std::set<uint32_t> nodes;
    int max_node_id = 0; //it decides the size of graph
    //Init the first line.
//    node_edges[0].push_back(nodes_num);    //use last to_node_id as vertex num
//    node_edges[0].push_back(edges_num);
    former_from_node_id = 0;   // former node id set to invalid for initialization
    long int temp_count = 0;
    long int edge_count = 0;
    // for directed graph
    while (input_file >> from_node_id) {
        input_file >> to_node_id;
        node_edges[from_node_id].push_back(to_node_id);
        nodes.insert(from_node_id);
        nodes.insert(to_node_id);
//        node_edges[to_node_id].push_back(from_node_id);
        max_node_id = (max_node_id < to_node_id) ? to_node_id : max_node_id;
        max_node_id = (max_node_id < from_node_id) ? from_node_id : max_node_id;
        temp_count += 2;
        edge_count += 1;
        if (temp_count % 1000000000 == 0) {
            std::cout << "reading" << temp_count / 1000000000 << "  billion edges" << std::endl;
        }
    }
//    node_edges[0][0] = max_node_id; //bug fixed: use last to_node_id as vertex num

    //output 3 files


    std::string::size_type iPos = graph_path.find_last_of('/') + 1;
    std::string filename = graph_path.substr(iPos, graph_path.length() - iPos);
    std::cout << filename << std::endl;

    // Get filename
    std::string name = filename.substr(0, filename.rfind("."));
    std::cout << name << std::endl;
    std::string output_path_prefix = "/mnt/wd_ultrastar/Spruce/Dataset/GFEDataset/";
    std::string output_name_p = name + ".properties";
    std::string output_name_e = name + ".e";
    std::string output_name_v = name + ".v";

    // Write properties
    std::string output_file_p = output_path_prefix + output_name_p;
    std::string output_file_v = output_path_prefix + output_name_v;
    std::string output_file_e = output_path_prefix + output_name_e;
    output_file.open(output_file_p, std::ios::out);
    if (!output_file)  
        std::cout << "File Error!! " << std::endl;

    std::cout << "Writing " << output_file_p << "..." << std::endl;
    output_file << "# Properties file of " << name << "\n";
    output_file << "\n";
    output_file << "# Filename of graph on local filesystem\n";
    output_file << "graph." << name << ".vertex-file = " << output_name_v << "\n";
    output_file << "graph." << name << ".edge-file = " << output_name_e << "\n";
    output_file << "\n";
    output_file << "# Graph metadata for reporting purposes" << "\n";
    output_file << "graph." << name << ".meta.vertices = " << std::to_string(nodes.size()) << "\n";
    output_file << "graph." << name << ".meta.edges = " << std::to_string(edge_count) << "\n";
    output_file << "\n";
    output_file << "# Properties describing the graph format" << "\n";
    output_file << "graph." << name << ".directed = false" << "\n";
    output_file << "\n";
    output_file << "# List of supported algorithms on the graph" << "\n";
    output_file << "graph." << name << ".algorithms = bfs, cdlp, lcc, pr, wcc" << "\n";
    output_file << "\n\n";
    output_file << "#\n";
    output_file << "# Per-algorithm properties describing the input parameters to each algorithm\n";
    output_file << "#\n";
    output_file << "\n# Parameters for BFS\n";
    output_file << "graph." << name << ".bfs.source-vertex = 1" << "\n";
    output_file << "\n# Parameters for CDLP\n";
    output_file << "graph." << name << ".cdlp.max-iterations = 10" << "\n";
    output_file << "\n# No parameters for LCC\n";
    output_file << "\n# Parameters for PR\n";
    output_file << "graph." << name << ".pr.damping-factor = 0.85" << "\n";
    output_file << "graph." << name << ".pr.num-iterations = 10" << "\n";
    output_file << "\n# No parameters for WCC\n";
    output_file.close();

    //Write vertex file:
    std::cout << "Writing " << output_file_v << "..." << std::endl;
    output_file.open(output_file_v, std::ios::out);
    if (!output_file)  
        std::cout << "File Error!! " << std::endl;
    for(auto it=nodes.begin(); it!=nodes.end(); it++){
        output_file << std::to_string(*it) << "\n";
    }
    output_file.close();

    //Write edge file
    std::cout << "Writing " << output_file_e << "..." << std::endl;
    output_file.open(output_file_e, std::ios::out);
    if (!output_file)  
        std::cout << "File Error!! " << std::endl;
    for (int i = 0; i < node_edges.size(); i++) {
        if(node_edges[i].size() != 0) {
            for (int j = 0; j < node_edges[i].size(); j++) {
                output_file << std::to_string(i) << " " << std::to_string(node_edges[i][j]) << "\n";
            }
        }
    }
    output_file.close();
    return nodes_num;
}


int AnalyzeDataset(const std::string &graph_path) {
    std::cout << "Analyzing Dataset..." << std::endl;
    // Init input and output stream
    std::ifstream input_file(graph_path, std::ios::in);
    std::string output_path = graph_path + ".analyze";
    std::ofstream output_file;
    std::string twitter_2010 = "twitter-2010";
    std::string hyperlink = "Hyperlink";
    std::string graph500 = "graph500";

    if (!input_file)  
        std::cout << "File Error!! " << std::endl;

    // Read Graph Prologues
    std::string temp_line;
    int nodes_num, edges_num, cnt = 0;
    if(graph_path.find(twitter_2010) != std::string::npos || graph_path.find(hyperlink) != std::string::npos ||
       graph_path.find(graph500) != std::string::npos || graph_path.find("GFE") != std::string::npos || graph_path.find("g500")!=std::string::npos
       || graph_path.find("com-lj.ungraph.") != std::string::npos){
        nodes_num = 41652230;
        edges_num = 2044000000; //only use edges_num that's enough
    }
    else {
        while (cnt++ < 4) {
            getline(input_file, temp_line);
            if (cnt == 3) {
                //Read number of nodes
                std::stringstream ss(temp_line);
                ss >> temp_line;
                ss >> temp_line;
                ss >> nodes_num;
                ss >> temp_line;
                ss >> edges_num;
            }
            std::cout << temp_line << std::endl;
        }
    }
    long int vertex_num, edge_num;
    if(graph_path.find("amazon")!=std::string::npos){
        vertex_num = 548551;
        edge_num = 925872;
    }
    else if(graph_path.find("LiveJournal")!=std::string::npos) {
        edge_num = 68993773;
        vertex_num = 4847571;
    }
    else if (graph_path.find("twitter")!=std::string::npos) {
        edge_num = 1468365182;
        vertex_num = 41652230;
    }
    else if (graph_path.find("Hyperlink")!=std::string::npos) {
        edge_num = 2043203933;
        vertex_num = 101717775;
    }
    else if (graph_path.find("500-26")!=std::string::npos) {
        vertex_num = 67108863;
        edge_num = 1051922853;
    }
    else if (graph_path.find("uniform-26")!=std::string::npos) {
        vertex_num = 49207468;
        edge_num = 1051922853;
    }
    else if (graph_path.find("uniform-24")!=std::string::npos) {
        vertex_num = 13306413;
        edge_num = 260379520;
    }
    else if (graph_path.find("500-24")!=std::string::npos) {
        vertex_num = 16777215;
        edge_num = 270000000;
    }
    else {
        vertex_num = 16777216;
        edge_num = edges_num;
    }



    // Read Graph Data(Sorted)
    uint32_t from_node_id, to_node_id, former_from_node_id;
    std::vector<std::vector<uint32_t>> node_edges(edge_num);    //Use edges_num for count, because some ids may not exist in txt
    int max_node_id = 0; //it decides the size of graph
    //Init the first line.
//    node_edges[0].push_back(nodes_num);    //use last to_node_id as vertex num
//    node_edges[0].push_back(edges_num);
    former_from_node_id = 0;   // former node id set to invalid for initialization
    long int temp_count = 0;
    double temp_double;
    if (graph_path.find("GFE") != std::string::npos || graph_path.find("g500")!=std::string::npos || graph_path.find("ungraph")!=std::string::npos) {
        while (input_file >> from_node_id) {
            input_file >> to_node_id;
            if (graph_path.find("dota") != std::string::npos) {
                // dota is a weighted graph, ignore graph weight;
                input_file >> temp_double;
            }
            node_edges[from_node_id].push_back(to_node_id);
            node_edges[to_node_id].push_back(from_node_id);
            max_node_id = (max_node_id < to_node_id) ? to_node_id : max_node_id;
            max_node_id = (max_node_id < from_node_id) ? from_node_id : max_node_id;
            temp_count += 2;
            if (temp_count % 1000000000 == 0) {
                std::cout << "reading" << temp_count / 1000000000 << "  billion edges" << std::endl;
            }
        }
    }
    else {
        while (input_file >> from_node_id) {
            input_file >> to_node_id;
            node_edges[from_node_id].push_back(to_node_id);
//            node_edges[to_node_id].push_back(from_node_id);
            max_node_id = (max_node_id < to_node_id) ? to_node_id : max_node_id;
            max_node_id = (max_node_id < from_node_id) ? from_node_id : max_node_id;
            temp_count++;
            if (temp_count % 1000000000 == 0) {
                std::cout << "reading" << temp_count / 1000000000 << "  billion edges" << std::endl;
            }
        }
    }
//    node_edges[0][0] = max_node_id; //bug fixed: use last to_node_id as vertex num

    output_file.open(output_path, std::ios::out | std::ios::trunc);

    //Write output file
    output_file << "AdjacencyGraph\n";
    output_file << max_node_id + 1 << "\n"; //begin from 0
    output_file << temp_count << "\n";

    //initialize analyze data
    int i, j;
    long int offset = 0;
    long long max_degree = 0;
    long long max_gap = 0;
    long long total_degrees = 0;
    long long total_gaps = 0;
    long long temp_gap = 0;
    long long vertex_num_cnt = 0;
    std::vector<long long> degree_numbers;

    //write offset
    //node : node_edges[max_node_id] may have edges
    i = 0;
    while (!node_edges[i].size()) {
        i++;
    }
    for(i; i < max_node_id + 1; i++){
//        output_file << std::to_string(offset) << "\n";
        if (!node_edges[i].size()) {
            temp_gap++;
        }
        else {
            vertex_num_cnt++;
            temp_gap++;
            if (node_edges[i].size() > max_degree) {
                max_degree = node_edges[i].size();
            }
            total_degrees += node_edges[i].size();
            if (temp_gap > max_gap) {
                max_gap = temp_gap;
            }
            degree_numbers.push_back(node_edges[i].size());
            temp_gap = 0;
        }
    }
//    std::cout << offset <<std::endl;
    //write edge
    std::cout << graph_path << "\n";
    output_file << graph_path << "\n";
    std::cout << "Max Gap: " << max_gap << "\n";
    output_file << "Max Gap: " << max_gap << "\n";
    std::cout << "Average Gap: " << std::fixed << std::setprecision(3) << (double)max_node_id/(double)vertex_num_cnt << std::endl;
    output_file << "Average Gap: " << std::fixed << std::setprecision(3) << (double)max_node_id/(double)vertex_num_cnt << "\n";
    std::cout << "Max Degree: " << max_degree << "\n";
    output_file << "Max Degree: " << max_degree << "\n";
    std::cout << "Average Degree: " << std::fixed << std::setprecision(3) << (double)total_degrees/(double)vertex_num_cnt << std::endl;
    output_file << "Average Degree: " << std::fixed << std::setprecision(3) << (double)total_degrees/(double)vertex_num_cnt << "\n";

    std::sort(degree_numbers.begin(), degree_numbers.end());
    long long index_1_percentile = degree_numbers.size() - degree_numbers.size() / 100 - 1;
    long long value_1_percentile = degree_numbers[index_1_percentile];
    long long index_0_1_percentile = degree_numbers.size() - degree_numbers.size() / 1000 - 1;
    long long value_0_1_percentile = degree_numbers[index_0_1_percentile];

    std::cout << "Top 1 percentile: " << value_1_percentile << "\n";
    std::cout << "Top 0.1 percentile: " << value_0_1_percentile << "\n";
    output_file << "Top 1 percentile: " << value_1_percentile << "\n";
    output_file << "Top 0.1 percentile: " << value_0_1_percentile << "\n";

    return 0;
}

int RemovetopthousandDataset(const std::string &graph_path) {
    std::cout << "Removing top 0.1 percent vertices from Dataset..." << std::endl;
    // Init input and output stream
    std::ifstream input_file(graph_path, std::ios::in);
    std::string output_path = graph_path + ".re.el";
    std::ofstream output_file;
    std::string twitter_2010 = "twitter-2010";
    std::string hyperlink = "Hyperlink";
    std::string graph500 = "graph500";

    if (!input_file)  
        std::cout << "File Error!! " << std::endl;

    // Read Graph Prologues
    std::string temp_line;
    int nodes_num, edges_num, cnt = 0;
    if(graph_path.find(twitter_2010) != std::string::npos || graph_path.find(hyperlink) != std::string::npos ||
       graph_path.find(graph500) != std::string::npos || graph_path.find("GFE") != std::string::npos || graph_path.find("g500")!=std::string::npos
       || graph_path.find("com-lj.ungraph.") != std::string::npos){
        nodes_num = 41652230;
        edges_num = 2044000000; //only use edges_num that's enough
    }
    else {
        while (cnt++ < 4) {
            getline(input_file, temp_line);
            if (cnt == 3) {
                //Read number of nodes
                std::stringstream ss(temp_line);
                ss >> temp_line;
                ss >> temp_line;
                ss >> nodes_num;
                ss >> temp_line;
                ss >> edges_num;
            }
            std::cout << temp_line << std::endl;
        }
    }
    long int vertex_num, edge_num;
    if(graph_path.find("amazon")!=std::string::npos){
        vertex_num = 548551;
        edge_num = 925872;
    }
    else if(graph_path.find("LiveJournal")!=std::string::npos) {
        edge_num = 68993773;
        vertex_num = 4847571;
    }
    else if (graph_path.find("twitter")!=std::string::npos) {
        edge_num = 1468365182;
        vertex_num = 41652230;
    }
    else if (graph_path.find("Hyperlink")!=std::string::npos) {
        edge_num = 2043203933;
        vertex_num = 101717775;
    }
    else if (graph_path.find("500-26")!=std::string::npos) {
        vertex_num = 67108863;
        edge_num = 1051922853;
    }
    else if (graph_path.find("uniform-26")!=std::string::npos) {
        vertex_num = 49207468;
        edge_num = 1051922853;
    }
    else if (graph_path.find("uniform-24")!=std::string::npos) {
        vertex_num = 13306413;
        edge_num = 260379520;
    }
    else if (graph_path.find("500-24")!=std::string::npos) {
        vertex_num = 16777215;
        edge_num = 270000000;
    }
    else {
        vertex_num = 16777216;
        edge_num = edges_num;
    }



    // Read Graph Data(Sorted)
    uint32_t from_node_id, to_node_id, former_from_node_id;
    std::vector<std::vector<uint32_t>> node_edges(edge_num);    //Use edges_num for count, because some ids may not exist in txt
    std::vector<std::vector<uint32_t>> node_edges_write(edge_num);
    int max_node_id = 0; //it decides the size of graph
    //Init the first line.
//    node_edges[0].push_back(nodes_num);    //use last to_node_id as vertex num
//    node_edges[0].push_back(edges_num);
    former_from_node_id = 0;   // former node id set to invalid for initialization
    long int temp_count = 0;
    double temp_double;
    if (graph_path.find("GFE") != std::string::npos || graph_path.find("g500")!=std::string::npos || graph_path.find("ungraph")!=std::string::npos) {
        while (input_file >> from_node_id) {
            input_file >> to_node_id;
            if (graph_path.find("dota") != std::string::npos) {
                // dota is a weighted graph, ignore graph weight;
                input_file >> temp_double;
            }
            node_edges[from_node_id].push_back(to_node_id);
            node_edges[to_node_id].push_back(from_node_id);
            node_edges_write[from_node_id].push_back(to_node_id);
            max_node_id = (max_node_id < to_node_id) ? to_node_id : max_node_id;
            max_node_id = (max_node_id < from_node_id) ? from_node_id : max_node_id;
            temp_count += 2;
            if (temp_count % 1000000000 == 0) {
                std::cout << "reading" << temp_count / 1000000000 << "  billion edges" << std::endl;
            }
        }
    }
    else {
        while (input_file >> from_node_id) {
            input_file >> to_node_id;
            node_edges[from_node_id].push_back(to_node_id);
//            node_edges[to_node_id].push_back(from_node_id);
            node_edges_write[from_node_id].push_back(to_node_id);
            max_node_id = (max_node_id < to_node_id) ? to_node_id : max_node_id;
            max_node_id = (max_node_id < from_node_id) ? from_node_id : max_node_id;
            temp_count++;
            if (temp_count % 1000000000 == 0) {
                std::cout << "reading" << temp_count / 1000000000 << "  billion edges" << std::endl;
            }
        }
    }
//    node_edges[0][0] = max_node_id; //bug fixed: use last to_node_id as vertex num

    output_file.open(output_path, std::ios::out | std::ios::trunc);

    //Write output file
//    output_file << "AdjacencyGraph\n";
//    output_file << max_node_id + 1 << "\n"; //begin from 0
//    output_file << temp_count << "\n";

    //initialize analyze data
    int i, j;
    long int offset = 0;
    long long max_degree = 0;
    long long max_gap = 0;
    long long total_degrees = 0;
    long long total_gaps = 0;
    long long temp_gap = 0;
    long long vertex_num_cnt = 0;
    std::vector<long long> degree_numbers;

    //write offset
    //node : node_edges[max_node_id] may have edges
    i = 0;
    while (!node_edges[i].size()) {
        i++;
    }
    for(i; i < max_node_id + 1; i++){
//        output_file << std::to_string(offset) << "\n";
        if (!node_edges[i].size()) {
            temp_gap++;
        }
        else {
            vertex_num_cnt++;
            temp_gap++;
            if (node_edges[i].size() > max_degree) {
                max_degree = node_edges[i].size();
            }
            total_degrees += node_edges[i].size();
            if (temp_gap > max_gap) {
                max_gap = temp_gap;
            }
            degree_numbers.push_back(node_edges[i].size());
            temp_gap = 0;
        }
    }
//    std::cout << offset <<std::endl;
    //write edge
//    std::cout << graph_path << "\n";
//    output_file << graph_path << "\n";
//    std::cout << "Max Gap: " << max_gap << "\n";
//    output_file << "Max Gap: " << max_gap << "\n";
//    std::cout << "Average Gap: " << std::fixed << std::setprecision(3) << (double)max_node_id/(double)vertex_num_cnt << std::endl;
//    output_file << "Average Gap: " << std::fixed << std::setprecision(3) << (double)max_node_id/(double)vertex_num_cnt << "\n";
//    std::cout << "Max Degree: " << max_degree << "\n";
//    output_file << "Max Degree: " << max_degree << "\n";
//    std::cout << "Average Degree: " << std::fixed << std::setprecision(3) << (double)total_degrees/(double)vertex_num_cnt << std::endl;
//    output_file << "Average Degree: " << std::fixed << std::setprecision(3) << (double)total_degrees/(double)vertex_num_cnt << "\n";

    std::sort(degree_numbers.begin(), degree_numbers.end());
    long long index_1_percentile = degree_numbers.size() - degree_numbers.size() / 100 - 1;
    long long value_1_percentile = degree_numbers[index_1_percentile];
    long long index_0_1_percentile = degree_numbers.size() - degree_numbers.size() / 1000 - 1;
    long long value_0_1_percentile = degree_numbers[index_0_1_percentile];

    std::cout << "Top 1 percentile: " << value_1_percentile << "\n";
    std::cout << "Top 0.1 percentile: " << value_0_1_percentile << "\n";
//    output_file << "Top 1 percentile: " << value_1_percentile << "\n";
//    output_file << "Top 0.1 percentile: " << value_0_1_percentile << "\n";
    for(i = 0; i < max_node_id + 1; i++){
//        output_file << std::to_string(offset) << "\n";
        if (node_edges[i].size() <= value_0_1_percentile) {
            continue;
        }
        else {
            for (j = 0; j < node_edges[i].size(); j++) {
                auto k = node_edges[i][j];

                auto it = std::find(node_edges_write[i].begin(), node_edges_write[i].end(), k);
                if (it != node_edges_write[i].end()) {
                    node_edges_write[i].erase(it);  
                }
                it = std::find(node_edges_write[k].begin(), node_edges_write[k].end(), i);
                if (it != node_edges_write[k].end()) {
                    node_edges_write[k].erase(it);  
                }
//
//                output_file << std::to_string(i);
//                output_file << " ";
//                output_file << std::to_string(node_edges[i][j]);
//                output_file << "\n";
            }
        }
    }
    for(i = 0; i < max_node_id + 1; i++){
            for (j = 0; j < node_edges_write[i].size(); j++) {
                output_file << std::to_string(i);
                output_file << " ";
                output_file << std::to_string(node_edges_write[i][j]);
                output_file << "\n";
            }
    }
    output_file.close();
    return 0;
}
