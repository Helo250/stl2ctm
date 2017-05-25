//-----------------------------------------------------------------------------
// Product:     OpenCTM tools
// File:        stl.cpp
// Description: Implementation of the STL file format importer/exporter.
//-----------------------------------------------------------------------------
// Copyright (c) 2009-2010 Marcus Geelnard
//
// This software is provided 'as-is', without any express or implied
// warranty. In no event will the authors be held liable for any damages
// arising from the use of this software.
//
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
//
//     1. The origin of this software must not be misrepresented; you must not
//     claim that you wrote the original software. If you use this software
//     in a product, an acknowledgment in the product documentation would be
//     appreciated but is not required.
//
//     2. Altered source versions must be plainly marked as such, and must not
//     be misrepresented as being the original software.
//
//     3. This notice may not be removed or altered from any source
//     distribution.
//-----------------------------------------------------------------------------

#include <stdexcept>
#include <fstream>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include "stl.h"

#ifdef _MSC_VER
typedef unsigned int uint32;
#else
#include <stdint.h>
typedef uint32_t uint32;
#endif

#define STLA 0
#define STLB 1
#define COR3_MAX 200000
#define FACE_MAX 200000
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
using namespace std;

bool stl = STLB;
float volume;
float surface_area;
float max_x = 0.0;
float max_y = 0.0;
float max_z = 0.0;
float min_x = 0.0;
float min_y = 0.0;
float min_z = 0.0;
float model_x = 0.0;
float model_y = 0.0;
float model_z = 0.0;
// unsigned char


/// Read a 32-bit integer, endian independent.
static uint32 ReadInt32(istream &aStream)
{
  unsigned char buf[4];
  aStream.read((char *) buf, 4);
  return ((uint32) buf[0]) | (((uint32) buf[1]) << 8) |
         (((uint32) buf[2]) << 16) | (((uint32) buf[3]) << 24);
}

/// Write a 32-bit integer, endian independent.
static void WriteInt32(ostream &aStream, uint32 aValue)
{
  unsigned char buf[4];
  buf[0] = aValue & 255;
  buf[1] = (aValue >> 8) & 255;
  buf[2] = (aValue >> 16) & 255;
  buf[3] = (aValue >> 24) & 255;
  aStream.write((char *) buf, 4);
}

/// Read a Vector3, endian independent.
static Vector3 ReadVector3(istream &aStream)
{
  union {
    uint32 i;
    float  f;
  } val;
  Vector3 result;
  val.i = ReadInt32(aStream);
  result.x = val.f;
  val.i = ReadInt32(aStream);
  result.y = val.f;
  val.i = ReadInt32(aStream);
  result.z = val.f;
  return result;
}

/// Write a Vector3, endian independent.
static void WriteVector3(ostream &aStream, Vector3 aValue)
{
  union {
    uint32 i;
    float  f;
  } val;
  val.f = aValue.x;
  WriteInt32(aStream, val.i);
  val.f = aValue.y;
  WriteInt32(aStream, val.i);
  val.f = aValue.z;
  WriteInt32(aStream, val.i);
}

/// Vertex class used when reading and joining the triangle vertices.
class SortVertex {
  public:
    float x, y, z;
    uint32 mOldIndex;

    bool operator<(const SortVertex &v) const
    {
      return (x < v.x) || ((x == v.x) && ((y < v.y) || ((y == v.y) && (z < v.z))));
    }
};

string remove_space(string& str){
    string buff(str);
    char space = ' ';
    str.assign(buff.begin() + buff.find_first_not_of(space),
        buff.begin() + buff.find_last_not_of(space) + 1);
    return str;
}

float SignedVolumeOfTriangle(SortVertex p1, SortVertex p2, SortVertex p3) {
    float v321 = p3.x * p2.y * p1.z;
    float v231 = p2.x * p3.y * p1.z;
    float v312 = p3.x * p1.y * p2.z;
    float v132 = p1.x * p3.y * p2.z;
    float v213 = p2.x * p1.y * p3.z;
    float v123 = p1.x * p2.y * p3.z;
    return (1.0f/6.0f)*(-v321 + v231 + v312 - v132 - v213 + v123);
}

float SingedAreaOfTriangle(SortVertex p1, SortVertex p2, SortVertex p3) {
  float a = (p2.y-p1.y)*(p3.z-p1.z) - (p3.y-p1.y)*(p2.z-p1.z);
  float b = (p2.z-p1.z)*(p3.x-p1.x) - (p2.x-p1.x)*(p3.z-p1.z);
  float c = (p2.x-p1.x)*(p3.y-p1.y) - (p2.y-p1.y)*(p3.x-p1.x);
  return (1.0f/2.0f)*sqrt(a*a + b*b + c*c);
}

/// Import an STL file from a file.
void Import_STL(const char * aFileName, Mesh * aMesh)
{
  // Clear the mesh
  aMesh->Clear();

  // Open the input file
  ifstream f(aFileName, ios_base::in | ios_base::binary);
  if(f.fail())
    throw runtime_error("Could not open input file.");

  // Get the file size
  f.seekg(0, ios_base::end);
  uint32 fileSize = (uint32) f.tellg();
  string line;
  string::size_type position;
  f.seekg(0, ios_base::beg);
  if(fileSize < 84)
    throw runtime_error("Invalid format - not a valid STL file.");
  std::getline(f, line);

  // parse the file is binary or ascii
  position = line.find("solid");
  if (position != line.npos) {
    std::getline(f, line);
    line = remove_space(line);
    if (line.find("facet") == 0) {
      stl = STLA;
    }
  }

  if (stl == STLB){
    f.seekg(0, ios_base::beg);
    // Read header (80 character comment + triangle count)
    char comment[81];
    f.read(comment, 80);
    comment[80] = 0;
    aMesh->mComment = string(comment);
    uint32 triangleCount = ReadInt32(f);
    if(fileSize != (84 + triangleCount * 50))
      throw runtime_error("Invalid format - not a valid STL file.");

    if(triangleCount > 0)
    {
      // Read all the triangle data
      vector<SortVertex> vertices;
      vertices.resize(triangleCount * 3);
      for(uint32 i = 0; i < triangleCount; ++ i)
      {
        // Skip the flat normal
        f.seekg(12, ios_base::cur);

        // Read the three triangle vertices
        for(uint32 j = 0; j < 3; ++ j)
        {
          Vector3 v = ReadVector3(f);
          uint32 index = i * 3 + j;
          vertices[index].x = v.x;
          vertices[index].y = v.y;
          vertices[index].z = v.z;
          vertices[index].mOldIndex = index;
          max_x = MAX(max_x, v.x);
          max_y = MAX(max_y, v.y);
          max_z = MAX(max_z, v.z);
          min_x = MIN(min_x, v.x);
          min_y = MIN(min_y, v.y);
          min_z = MIN(min_z, v.z);
        }
        surface_area += SingedAreaOfTriangle(vertices[i*3], vertices[i*3+1], vertices[i*3+2]);
        volume += SignedVolumeOfTriangle(vertices[i*3], vertices[i*3+1], vertices[i*3+2]);
        // Ignore the two fill bytes
        f.seekg(2, ios_base::cur);
      }
      // caculate the model x y z and volume
      model_x = max_x - min_x;
      model_y = max_y - min_y;
      model_z = max_z - min_z;
      volume = std::abs(volume);
      printf("\nThe model_x: %.0f mm\nThe model_y: %.0f mm\nThe model_z: %.0f mm\n", model_x, model_y, model_z);
      printf("The model surface_area: %.0f mm2\n", surface_area);
      printf("The model volume: %.0f mm3\n", volume);
      // cout << "\nThe model_x: " << setprecision(2) << model_x/10 << " cm" << endl;
      // cout << "The model_y: " << setprecision(2) << model_y/10 << " cm" << endl;
      // cout << "The model_z: " << setprecision(2) << model_z/10 << " cm" << endl;
      // cout << "The volume: " << setprecision(2) << volume << " cm3" << endl;
      // Make sure that no redundant copies of vertices exist (STL files are full
      // of vertex duplicates, so remove the redundancy), and store the data in
      // the mesh object
      sort(vertices.begin(), vertices.end());
      aMesh->mVertices.resize(vertices.size());
      aMesh->mIndices.resize(vertices.size());
      SortVertex * firstEqual = &vertices[0];
      int vertIdx = -1;
      for(uint32 i = 0; i < vertices.size(); ++ i)
      {
        if((i == 0) ||
           (vertices[i].z != firstEqual->z) ||
           (vertices[i].y != firstEqual->y) ||
           (vertices[i].x != firstEqual->x))
        {
          firstEqual = &vertices[i];
          ++ vertIdx;
          aMesh->mVertices[vertIdx] = Vector3(firstEqual->x, firstEqual->y, firstEqual->z);
        }
        // volume += std::abs(SignedVolumeOfTriangle(vertices[i], vertices[i+1], vertices[i+2]));
        aMesh->mIndices[vertices[i].mOldIndex] = vertIdx;
      }
      aMesh->mVertices.resize(vertIdx + 1);
    }
  }
  else{
    vector<SortVertex> vertices;
    uint32 index = 0;
    bool loop = false;

    // string line;
    try{
      while(!f.eof()){
          Vector3 result;
          std::getline(f, line);
          line = remove_space(line);
          // if (line.find("facet")==0){
          //     sscanf(line.c_str(), "%*s %e %e %e", &result.x, &result.y, &result.z);
          // }
          // else
          if (line.find("outer") == 0){
              loop = true;
          }
          else if (line.find("vertex") == 0 && loop){
            sscanf(line.c_str(), "%*s %e %e %e", &result.x, &result.y, &result.z);
            vertices.resize(index + 1);
            vertices[index].x = result.x;
            vertices[index].y = result.y;
            vertices[index].z = result.z;
            vertices[index].mOldIndex = index;
            max_x = MAX(max_x, result.x);
            max_y = MAX(max_y, result.y);
            max_z = MAX(max_z, result.z);
            min_x = MIN(min_x, result.x);
            min_y = MIN(min_y, result.y);
            min_z = MIN(min_z, result.z);
            if (index % 3 == 0 && index > 0){
              volume += SignedVolumeOfTriangle(vertices[index-3], vertices[index-2], vertices[index-1]);
              surface_area += SingedAreaOfTriangle(vertices[index-3], vertices[index-2], vertices[index-1]);
            }
            index ++;
          }
          else if (line.find("endloop") == 0){
            loop = false;
          }
      }
      // caculate the model x y z and volume
      model_x = (max_x - min_x);
      model_y = (max_y - min_y);
      model_z = (max_z - min_z);
      volume = std::abs(volume);
      printf("\nThe model_x: %.0f mm\nThe model_y: %.0f mm\nThe model_z: %.0f mm\n", model_x, model_y, model_z);
      printf("The model surface_area: %.0f mm2\n", surface_area);
      printf("The model volume: %.0f mm3\n", volume);
      // cout << "\nThe model_x: " << setprecision(2) << model_x/10 << " cm" << endl;
      // cout << "The model_y: " << setprecision(2) << model_y/10 << " cm" << endl;
      // cout << "The model_z: " << setprecision(2) << model_z/10 << " cm" << endl;
      // cout << "The volume: " << setprecision(2) << volume << " cm3" << endl;
      // Make sure that no redundant copies of vertices exist (STL files are full
      // of vertex duplicates, so remove the redundancy), and store the data in
      // the mesh object
      sort(vertices.begin(), vertices.end());
      aMesh->mVertices.resize(vertices.size());
      aMesh->mIndices.resize(vertices.size());
      SortVertex * firstEqual = &vertices[0];
      int vertIdx = -1;
      for(uint32 i = 0; i < vertices.size(); ++ i)
      {
        if((i == 0) ||
           (vertices[i].z != firstEqual->z) ||
           (vertices[i].y != firstEqual->y) ||
           (vertices[i].x != firstEqual->x))
        {
          firstEqual = &vertices[i];
          ++ vertIdx;
          aMesh->mVertices[vertIdx] = Vector3(firstEqual->x, firstEqual->y, firstEqual->z);
        }
        aMesh->mIndices[vertices[i].mOldIndex] = vertIdx;
      }
      aMesh->mVertices.resize(vertIdx + 1);
    }
    catch(exception &e){
      // pass
    }
  }
  // Close the input file
  f.close();
}

/// Export an STL file to a file.
void Export_STL(const char * aFileName, Mesh * aMesh, Options &aOptions)
{
  // Open the output file
  ofstream f(aFileName, ios_base::out | ios_base::binary);
  if(f.fail())
    throw runtime_error("Could not open output file.");

  // Write header (80-character comment + triangle count)
  char comment[80];
  for(uint32 i = 0; i < 80; ++ i)
  {
    if(i < aMesh->mComment.size())
      comment[i] = aMesh->mComment[i];
    else
      comment[i] = 0;
  }
  f.write(comment, 80);
  uint32 triangleCount = aMesh->mIndices.size() / 3;
  WriteInt32(f, triangleCount);

  // Write the triangle data
  for(uint32 i = 0; i < triangleCount; ++ i)
  {
    // Get the triangle vertices
    Vector3 v1 = aMesh->mVertices[aMesh->mIndices[i * 3]];
    Vector3 v2 = aMesh->mVertices[aMesh->mIndices[i * 3 + 1]];
    Vector3 v3 = aMesh->mVertices[aMesh->mIndices[i * 3 + 2]];

    // Calculate the triangle normal
    Vector3 n1 = v2 - v1;
    Vector3 n2 = v3 - v1;
    Vector3 n = Normalize(Cross(n1, n2));

    // Write the triangle normal
    WriteVector3(f, n);

    // Coordinates
    WriteVector3(f, v1);
    WriteVector3(f, v2);
    WriteVector3(f, v3);

    // Set the two fill bytes to zero
    f.put(0);
    f.put(0);
  }

  // Close the output file
  f.close();
}
