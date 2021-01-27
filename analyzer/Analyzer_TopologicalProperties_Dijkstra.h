/*------------------------------------------------------------------------------
    ooo      L   attice-based  |
  o\.|./o    e   xtensible     | LeMonADE: An Open Source Implementation of the
 o\.\|/./o   Mon te-Carlo      |           Bond-Fluctuation-Model for Polymers
oo---0---oo  A   lgorithm and  |
 o/./|\.\o   D   evelopment    | Copyright (C) 2018, 2021 by
  o/.|.\o    E   nvironment    | LeMonADE Principal Developers (Ron Dockhorn)
    ooo                        |
--------------------------------------------------------------------------------

This file is part of LeMonADE.

LeMonADE is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

LeMonADE is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with LeMonADE.  If not, see <http://www.gnu.org/licenses/>.

------------------------------------------------------------------------------*/

#ifndef ANALYZER_TOPOLOGICAL_PROPERTIES_DIJKSTRA
#define ANALYZER_TOPOLOGICAL_PROPERTIES_DIJKSTRA

#include <iostream>
#include <string>
#include <map>
#include <utility>
#include <functional>
#include <queue>
#include <vector>

#include <LeMonADE/utility/Vector3D.h>
#include <LeMonADE/analyzer/AbstractAnalyzer.h>
#include <LeMonADE/utility/ResultFormattingTools.h>

#include "Histogram1D.h"

/*******************************************************************************
 * CLASS DEFINITION (implementation of methods below)
 * ****************************************************************************/

template<class IngredientsType>
class Analyzer_TopologicalProperties_Dijkstra : public AbstractAnalyzer {
public:


    //constuctor
    Analyzer_TopologicalProperties_Dijkstra(const IngredientsType& ing);

    //initializes the groups. called explicitly or by Taskmanager::init()
    virtual void initialize();

    //does all the calculations esp in every !mcs
    virtual bool execute();

    //write out your results, in this case to standard output
    virtual void cleanup();

    typedef uint32_t NodeIdx;
    typedef std::map<NodeIdx, int> Nodelist;
    typedef std::map<NodeIdx, Nodelist> Graph;
    //typedef std::pair<NodeIdx, Nodelist> Node;
    typedef std::pair<int, NodeIdx> Edge; //! (tentative distance, idxNode)
    typedef std::vector<NodeIdx> NodeVector;

    void dijkstra(Graph &graph, NodeIdx source, Nodelist &distance);

    void findCentersOfTree(Graph &graph, Nodelist degreeNodes, NodeVector &centerNodes);

private:

    //holds a reference of the complete system
    const IngredientsType& ingredients;

    Histogram1D* HG_TopoDistance;
    Histogram1D* HG_FunctionalityTopoDistance;

    int nrOfSamples;

    //only used to make sure you initialize your groups before you do things
    bool initialized;
};

//template<class IngredientsType>
//bool operator>(const typename Analyzer_TopologicalProperties_Dijkstra<IngredientsType>::Edge& r, const typename Analyzer_TopologicalProperties_Dijkstra<IngredientsType>::Edge& k) {
//    return r.first > k.first;
//}


/*****************************************************************************
 * IMPLEMENTATION OF METHODS
 * **************************************************************************/

/* ****************************************************************************
 * constructor. only initializes some variables
 * ***************************************************************************/
template<class IngredientsType>
Analyzer_TopologicalProperties_Dijkstra<IngredientsType>::Analyzer_TopologicalProperties_Dijkstra(const IngredientsType& ing)
: ingredients(ing), initialized(false), nrOfSamples(0) {
}

/* **********************************************************************
 * initialize()
 *
 * this is called in the beginning only once - need for setting up your analyzer
 * **********************************************************************/
template<class IngredientsType>
void Analyzer_TopologicalProperties_Dijkstra<IngredientsType>::initialize() {
    std::cout << "init histogram" << std::endl;
    //HG_TopoDistance = new Histogram1D(0-0.5,ingredients.getMolecules().size()-0.5,ingredients.getMolecules().size());
    HG_TopoDistance = new Histogram1D(0 - 0.5, 20000 - 0.5, 20000);
    HG_FunctionalityTopoDistance = new Histogram1D(0 - 0.5, 20000 - 0.5, 20000);
    //set the initialized tag to true
    initialized = true;
    nrOfSamples = 0;

}

/* ***********************************************************************
 * execute()
 * this is where the calculation happens
 * ***********************************************************************/
template<class IngredientsType>
bool Analyzer_TopologicalProperties_Dijkstra<IngredientsType>::execute() {
    //check if analyzer have been initialized. if not, exit and explain
    if (initialized == false) {
        std::stringstream errormessage;
        errormessage << "Analyzer_TopologicalProperties_Dijkstra::execute()...Analyzer_TopologicalProperties_Dijkstra not initialized\n"
                << "Use Analyzer_TopologicalProperties_Dijkstra::initialize() or Taskmanager::init()\n";

        throw std::runtime_error(errormessage.str());
    }

    Graph treeGraph;


    Nodelist degreeNode; // Degree of node

    for (int k = 0; k < ingredients.getMolecules().size(); k++) {
        for (int l = 0; l < ingredients.getMolecules().getNumLinks(k); l++) {
            // connect from idxAA to idxBB with value ZZ
            // all (symmetric) connections
            treeGraph[k][ingredients.getMolecules().getNeighborIdx(k, l)] = 1;

        }
        // fill with degree of nodes
        degreeNode[k] = ingredients.getMolecules().getNumLinks(k);
    }


    //finding the center of the tree (either one or two)
    NodeVector centerNodes;
    findCentersOfTree(treeGraph, degreeNode, centerNodes);

    //select only one center (if structure has two)
    //calculate the topological distance of center to all nodes and fill histogram
    NodeVector::iterator itv = centerNodes.begin();
    //for (itv = centerNodes.begin(); itv != centerNodes.end(); ++itv)
    {

        nrOfSamples++;

        NodeIdx startNode = (*itv);
        std::cout << std::endl << "startNode: " << startNode << std::endl;

        Nodelist topo_distance;
        dijkstra(treeGraph, startNode, topo_distance);

        std::cout << "startNode -> node => distance" << std::endl;

        Nodelist::iterator it;
        for (it = topo_distance.begin(); it != topo_distance.end(); ++it) {
            std::cout << startNode << " -> " << it->first << "  => " << it->second << std::endl;

            // add value to the histogram
            HG_TopoDistance->AddValue(it->second);

            // add value to the histogram for number links (functionality)
            for (int o = 0; o < degreeNode[it->first]; o++) {
                HG_FunctionalityTopoDistance->AddValue(it->second);
            }
        }
    }
    

    return true;

}

template<class IngredientsType>
void Analyzer_TopologicalProperties_Dijkstra<IngredientsType>::findCentersOfTree(Graph &graph, Nodelist degreeNodes, NodeVector &centerNodes) {

    std::queue <NodeIdx> queueLeaves; // Queue for algorithm holding (reduced) leaves
    std::queue <NodeIdx> queueNextLeaves; // Queue for algorithm holding next (reduced) leaves
    std::queue <NodeIdx> tmpOldLeaves; // Queue for algorithm holding (reduced) leaves

    // Fill and start from leaves

    Nodelist::iterator it;
    for (it = degreeNodes.begin(); it != degreeNodes.end(); ++it) {
        int degreeOfNode = it->second;
        NodeIdx labelNode = it->first;

        if (degreeOfNode == 1) {
            queueLeaves.push(labelNode);
        }
    }

    NodeIdx Seperator = -1;

    queueLeaves.push(Seperator);

    //  while ( (queueLeaves.front()!=Seperator)||(queueLeaves.size() > 2) ) {
    while ((queueLeaves.size() > 0)) {
        NodeIdx node = queueLeaves.front();
        queueLeaves.pop();

        if (node == Seperator) {
            if (queueNextLeaves.empty())
                break;

            std::cout << "next iteration: " << node << std::endl;
            tmpOldLeaves = queueNextLeaves;
            queueLeaves = queueNextLeaves;

            if (queueNextLeaves.size() > 0)
                queueLeaves.push(Seperator);

            queueNextLeaves = std::queue <NodeIdx>();

            continue;
        }


        int dN_old = degreeNodes[node];

        degreeNodes[node]--;

        // new subgraph of all neighbors
        Nodelist tempgraph = graph[node];

        // connected nodes

        std::cout << "connected nodes from: " << node << std::endl;

        Nodelist::iterator it;
        for (it = tempgraph.begin(); it != tempgraph.end(); ++it) {

            std::cout << it->first << " with degree: " << degreeNodes[it->first];

            if (degreeNodes[it->first] > 0) {
                degreeNodes[it->first]--;

                std::cout << " reduced to degree: " << degreeNodes[it->first];
            }

            if (degreeNodes[it->first] == 1) {
                //queueLeaves.push(it->first);
                queueNextLeaves.push(it->first);

                std::cout << " added to queue the idx: " << (it->first);
            }

            std::cout << std::endl;
        }

        /*  std::cout << "myqueue contains: ";
                        while (!queueLeaves.empty())
                        {
                          std::cout << ' ' << queueLeaves.front();
                          queueLeaves.pop();
                        }
                        std::cout << '\n';
         */
    }

    std::cout << "center queue contains: ";
    while (!tmpOldLeaves.empty()) {
        centerNodes.push_back(tmpOldLeaves.front());
        std::cout << ' ' << tmpOldLeaves.front();
        tmpOldLeaves.pop();
    }
    std::cout << '\n';
}

template<class IngredientsType>
void Analyzer_TopologicalProperties_Dijkstra<IngredientsType>::dijkstra(Graph &graph, NodeIdx source, Nodelist &distance) {

    distance.clear(); //! clear all tentative distance information

    //! list of all nodes to be visit and addressed already with least distance on top
    std::priority_queue<Edge, std::vector<Edge>, std::greater<Edge> > queueNode;

    // starting node with tentative distance=0, starting node = source
    queueNode.push(Edge(0, source));

    while (!queueNode.empty()) {

        //get the element with least tentative distance
        Edge tmped = queueNode.top();

        // access the node index
        NodeIdx tmpnl = tmped.second;

        // removes the top element
        queueNode.pop();

        // if the node never visited before
        if (distance.count(tmpnl) == 0) {

            // tentative distance to the recent node
            int dist = tmped.first;

            // set the tentative distance to the node
            distance[tmpnl] = dist;

            // new subgraph of all neighbors
            Nodelist tempgraph = graph[tmpnl];

            Nodelist::iterator it;
            for (it = tempgraph.begin(); it != tempgraph.end(); ++it) {
                int distint = it->second;
                NodeIdx distlabel = it->first;
                queueNode.push(Edge(dist + distint, distlabel));
            }

        }
    }

}

/****************************************************************************
 * cleanup
 * this is where you get the final value and write it to file, std output, or 
 * whatever
 * this is called in the end only once - need for write output to a file or stdout
 * *************************************************************************/
template<class IngredientsType>
void Analyzer_TopologicalProperties_Dijkstra<IngredientsType>::cleanup() {

    // print results into a file
    // get the filename and path
    std::string filenameGeneral = ingredients.getName();
    // delete the .bfm in the name
    filenameGeneral.erase(ingredients.getName().length() - 4, ingredients.getName().length());

    // construct a list
    std::vector < std::vector<double> > tmpResultsHG_ZCOM;

    // we have 3 columns and row
    uint32_t columns = 6;
    uint32_t rows = 1;

    // we have columns
    tmpResultsHG_ZCOM.resize(columns);

    // we have rows
    for (int i = 0; i < columns; i++)
        tmpResultsHG_ZCOM[i].resize(0);

    // fill the histogram
    for (int bin = 0; bin < HG_TopoDistance->GetNrBins(); bin++) {
        if (HG_TopoDistance->GetNrInBin(bin) != 0) {
            tmpResultsHG_ZCOM[0].push_back(HG_TopoDistance->GetRangeInBin(bin));
            tmpResultsHG_ZCOM[1].push_back(HG_TopoDistance->GetNrInBin(bin) / (1.0 * nrOfSamples));
            tmpResultsHG_ZCOM[2].push_back(HG_TopoDistance->GetNrInBinNormiert(bin) / HG_TopoDistance->GetIntervallThickness());
            tmpResultsHG_ZCOM[3].push_back(HG_TopoDistance->GetCumulativeNrInBin(bin) / (1.0 * nrOfSamples));
            tmpResultsHG_ZCOM[4].push_back(HG_TopoDistance->GetCumulativeNrInBinNormiert(bin) / HG_TopoDistance->GetIntervallThickness());
            tmpResultsHG_ZCOM[5].push_back(HG_TopoDistance->GetNrInBin(bin));

        }
    }

    std::stringstream comment;
    comment << " File produced by analyzer Analyzer_TopologicalProperties_Dijkstra" << std::endl
            << " Histogram of topological distance from center with " << ingredients.getMolecules().size() << " monomers" << std::endl
            << " Total counts in all bins: " << HG_TopoDistance->GetNrCounts() << std::endl
            << " Total number of samples: " << nrOfSamples << std::endl
            << " Total number of bins " << HG_TopoDistance->GetNrBins() << std::endl
            << " Interval [" << HG_TopoDistance->GetRangeInBin(0) << " ; " << HG_TopoDistance->GetRangeInBin(HG_TopoDistance->GetNrBins()) << "]" << std::endl
            << " Interval thickness dI = " << HG_TopoDistance->GetIntervallThickness() << std::endl
            << std::endl
            << " d ... topological distance" << std::endl
            << " <c(d)> ... average counts (number monomers) at topological distance d (ensemble average)" << std::endl
            << " <cnorm(d)> ... average normalized counts by all counts at topological distance d (ensemble average)" << std::endl
            << " <s(d)> ... cumulative average counts at topological distance d = sum_d <c(d)> " << std::endl
            << " <snorm(d)> ... cumulative average normalized counts by all counts at topological distance d = sum_d <cnorm(d)> " << std::endl
            << " samples(d) ... total number of samples at topological distance d" << std::endl
            << "d <c(d)> <cnorm(d)> <s(d)> <snorm(d)> samples(d)";

    //new filename
    std::string filenameHG_ZCOM = filenameGeneral + "_TopologicalDistance.dat";

    ResultFormattingTools::writeResultFile(filenameHG_ZCOM, this->ingredients, tmpResultsHG_ZCOM, comment.str());


    // output functionality per topological distance


    // construct a list
    std::vector < std::vector<double> > tmpResultsHG_FuncDist;

    // we have 3 columns and row
    uint32_t columnsFuncDist = 5;
    //uint32_t rowsFuncDist = 1;

    // we have columns
    tmpResultsHG_FuncDist.resize(columnsFuncDist);

    // we have rows
    for (int i = 0; i < columns; i++)
        tmpResultsHG_FuncDist[i].resize(0);

    // fill the histogram
    for (int bin = 0; bin < HG_FunctionalityTopoDistance->GetNrBins(); bin++) {
        if (HG_FunctionalityTopoDistance->GetNrInBin(bin) != 0) {
            tmpResultsHG_FuncDist[0].push_back(HG_FunctionalityTopoDistance->GetRangeInBin(bin));
            tmpResultsHG_FuncDist[1].push_back(HG_TopoDistance->GetNrInBin(bin) / (1.0 * nrOfSamples));
            tmpResultsHG_FuncDist[2].push_back(HG_FunctionalityTopoDistance->GetNrInBin(bin) / HG_TopoDistance->GetNrInBin(bin));
            tmpResultsHG_FuncDist[3].push_back(HG_TopoDistance->GetNrInBin(bin));
            tmpResultsHG_FuncDist[4].push_back(HG_FunctionalityTopoDistance->GetNrInBin(bin));

        }
    }

    std::stringstream commentFuncDist;
    commentFuncDist << " File produced by analyzer Analyzer_TopologicalProperties_Dijkstra" << std::endl
            << " Histogram of functionality at topological distance from center with " << ingredients.getMolecules().size() << " monomers" << std::endl
            << " Total counts in all bins: " << HG_TopoDistance->GetNrCounts() << std::endl
            << " Total number of samples: " << nrOfSamples << std::endl
            << " Total number of bins " << HG_TopoDistance->GetNrBins() << std::endl
            << " Interval [" << HG_TopoDistance->GetRangeInBin(0) << " ; " << HG_TopoDistance->GetRangeInBin(HG_TopoDistance->GetNrBins()) << "]" << std::endl
            << " Interval thickness dI = " << HG_TopoDistance->GetIntervallThickness() << std::endl
            << std::endl
            << " d ... topological distance" << std::endl
            << " <c(d)> ... average counts (number monomers) at topological distance d (ensemble average)" << std::endl
            << " <f(d)> ... average functionality of monomers at topological distance d (ensemble average)" << std::endl
            << " samplesc(d) ... total number of samples of <c(d)> at topological distance d" << std::endl
            << " samplesf(d) ... total number of samples of <f(d)> at topological distance d" << std::endl
            << " d <c(d)> <f(d)> samplesc(d) samplesf(d)";

    //new filename
    std::string filenameHG_FuncDist = filenameGeneral + "_TopologicalDistanceFunctionality.dat";

    ResultFormattingTools::writeResultFile(filenameHG_FuncDist, this->ingredients, tmpResultsHG_FuncDist, commentFuncDist.str());


    delete HG_TopoDistance;
    delete HG_FunctionalityTopoDistance;
}



#endif /*ANALYZER_TOPOLOGICAL_PROPERTIES_DIJKSTRA*/


