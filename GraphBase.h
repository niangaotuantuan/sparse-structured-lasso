/*
 * BBox.cpp
 *
 *  Created on: May. 02nd, 2015
 *  Last modified on: May. 20th, 2015
 *      Author: Yanran Li
 */

#ifndef GRAPHBASE_H_
#define GRAPHBASE_H_
#include <iostream>
class GraphBase {
  protected:
    int numVertPotentials;
    int numEdgePotentials;
    float *vertPotentials;
    float *edgePotentials;
  public:
    GraphBase():vertPotentials(NULL), edgePotentials(NULL), numVertPotentials(0), numEdgePotentials(0){}
    GraphBase(int verticePotentialSize, int edgePotentialSize);
    virtual ~GraphBase();
    virtual void Generate();
    virtual void Save(const std::string fname);
    virtual void LoadByRow(const std::string fname);
};

#endif /* GRAPHBASE_H_ */
