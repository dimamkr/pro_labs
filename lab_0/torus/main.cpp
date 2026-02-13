#include <gmsh.h>

int main()
{
    gmsh::initialize();

    gmsh::model::add("thor");

    gmsh::option::setNumber("Mesh.MeshSizeFactor", 0.1);

    gmsh::model::occ::addTorus(0, 0, 0, 20, 5);

    gmsh::model::occ::synchronize();

    gmsh::model::mesh::generate(3);

    gmsh::fltk::run();
}