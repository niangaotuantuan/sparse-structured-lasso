/*
 * GraphBase.cpp
 *
 *  Created on: May. 02nd, 2015
 *  Last modified on: May. 21st, 2015
 *      Author: Yanran Li
 */

#include "GraphBase.h"
#include <iostream>

GraphBase::GraphBase(int verticePotentialSize, int edgePotentialSize)
: vertPotentials(NULL), edgePotentials(NULL)
{
  numVertPotentials = verticePotentialSize;
  numEdgePotentials = edgePotentialSize;
}

GraphBase::~GraphBase() {
}

void GraphBase::Generate(){
}

void GraphBase::LoadByRow(const std::string fname){
}

void GraphBase::Save(const std::string fname){
}

