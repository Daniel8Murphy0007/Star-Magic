#ifndef MODEL_LOADER_H
#define MODEL_LOADER_H

#include <string>
#include <vector>

struct MeshData
{
    std::vector<float> vertices;
    std::vector<float> normals;
    std::vector<float> texcoords;
    std::vector<unsigned int> indices;
};

bool loadOBJ(const std::string &path, MeshData &mesh);
void exportOBJ(const std::string &path, const MeshData &mesh);

#endif // MODEL_LOADER_H
