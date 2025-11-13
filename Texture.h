#ifndef TEXTURE_H
#define TEXTURE_H

#include <string>

class Texture
{
public:
    Texture() {}
    ~Texture() {}
    bool loadFromFile(const std::string &path) { return true; }
    void bind() {}
    void unbind() {}
};

#endif // TEXTURE_H
