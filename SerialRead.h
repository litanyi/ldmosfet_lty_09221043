#ifndef SERIALREAD_H
#define SERIALREAD_H
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <iomanip>
//jaumin
#include "PatchTopology.h"
#include "PatchGeometry.h"
#include "Patch.h"
#include "NodeVariable.h"
#include "NodeData.h"
#include "VariableDatabase.h"
#include "Utilities.h"

void Get_MaxDop(const std::string filename, hier::Patch<NDIM> & patch, double *max_doping)
{
    double max_dop=1;
    double max_index=0;
    ////////////////////////////////////////////////////////
    std::fstream input(filename.c_str());

    if(!input) {
        TBOX_ERROR("fail to open file: "<< filename << endl);
    }

    string line;

    getline(input,line);

    line.clear();

    while( getline(input, line) )
    {
        stringstream buf(line);
        int node_id;
        double dop_node;

        //read and convert data
        buf >> node_id;
        buf >> dop_node;

        max_index = node_id;
        if(max_dop<dop_node)
            max_dop = dop_node;
    }

    for(int index=0; index<=max_index; index++){
        max_doping[index] = max_dop * 1e-6;
    }

    input.close();
}

void SerialRead(const std::string filename, hier::Patch<NDIM> & patch, double *dop)
{

    // 获取当前网格片关联的网格片几何对象和网格片拓扑对象.
    tbox::Pointer<hier::PatchGeometry<NDIM> > patch_geo =
        patch.getPatchGeometry();
    tbox::Pointer<hier::PatchTopology<NDIM> > patch_top =
        patch.getPatchTopology();

    ////////////////////////////////////////////////////////
    std::fstream input(filename.c_str());

    if(!input) {
        TBOX_ERROR("fail to open file: "<< filename << endl);
    }

    string line;

    getline(input,line);

    line.clear();

    while( getline(input, line) )
    {
        stringstream buf(line);
        int node_id;
        double dop_node;

        //read and convert data
        buf >> node_id;
        buf >> dop_node;

        dop[node_id] = dop_node;
    }
    input.close();
    return;
}
#endif // SERIALREAD_H
