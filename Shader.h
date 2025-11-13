#ifndef SHADER_H
#define SHADER_H

#include <string>

class Shader
{
public:
    Shader() {}
    ~Shader() {}
    bool loadFromFile(const std::string &vertPath, const std::string &fragPath) { return true; }
    void use() {}
    void setUniform(const std::string &name, float value) {}
};

#endif // SHADER_H
