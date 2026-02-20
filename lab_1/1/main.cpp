#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
#include <array>
#include <filesystem>

#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkTetra.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkSmartPointer.h>

#include <gmsh.h>

// Класс расчётной точки
class CalcNode
{
    friend class CalcMesh;

protected:
    double x, y, z;
    double initX, initY, initZ;
    double smth;
    double vx, vy, vz;

public:
    CalcNode() : x(0.0), y(0.0), z(0.0), initX(0.0), initY(0.0), initZ(0.0), smth(0.0), vx(0.0), vy(0.0), vz(0.0) {}
    CalcNode(double x, double y, double z, double smth, double vx, double vy, double vz)
        : x(x), y(y), z(z), initX(x), initY(y), initZ(z), smth(smth), vx(vx), vy(vy), vz(vz) {}
};

// Класс элемента сетки (тетраэдр)
class Element
{
    friend class CalcMesh;

protected:
    unsigned long nodesIds[4];
};

// Класс расчётной сетки
class CalcMesh
{
public:
    std::vector<CalcNode> nodes;
    std::vector<Element> elements;
    double centerX, centerY, centerZ;

    CalcMesh(const std::vector<double> &nodesCoords, const std::vector<size_t> &tetrsPoints)
        : centerX(0.0), centerY(0.0), centerZ(0.0)
    {
        // Создаём узлы
        nodes.resize(nodesCoords.size() / 3);
        for (size_t i = 0; i < nodes.size(); i++)
        {
            double x = nodesCoords[i * 3];
            double y = nodesCoords[i * 3 + 1];
            double z = nodesCoords[i * 3 + 2];
            // Скалярное поле: координата x
            double smth = x;
            nodes[i] = CalcNode(x, y, z, smth, 0.0, 0.0, 0.0);
        }

        // Создаём тетраэдры
        elements.resize(tetrsPoints.size() / 4);
        for (size_t i = 0; i < elements.size(); i++)
        {
            elements[i].nodesIds[0] = tetrsPoints[i * 4] - 1;
            elements[i].nodesIds[1] = tetrsPoints[i * 4 + 1] - 1;
            elements[i].nodesIds[2] = tetrsPoints[i * 4 + 2] - 1;
            elements[i].nodesIds[3] = tetrsPoints[i * 4 + 3] - 1;
        }
    }

    void step(double t, double dt, int step)
    {
        updateCoords(t, dt);

        snapshot(step);
    }

    void updateCoords(double t, double dt)
    {
        double omega = 2.0 * M_PI * t;

        for (CalcNode &node : nodes)
        {
            double rx0 = node.initX;
            double ry0 = node.initY;
            double r = std::sqrt(rx0 * rx0 + ry0 * ry0);
            double theta0 = std::atan2(ry0, rx0);

            double theta = theta0 + omega * t;

            node.x = r * cos(theta);
            node.y = r * sin(theta);

            node.z += 10 * cos(5 * t) * dt;

            node.vx = -omega * r * sin(theta);
            node.vy = omega * r * cos(theta);
            node.vz = 10 * cos(5 * t);
        }
    }

    // Запись текущего состояния в VTU-файл
    void snapshot(unsigned int snap_number)
    {
        vtkSmartPointer<vtkUnstructuredGrid> grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
        vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();

        auto smth = vtkSmartPointer<vtkDoubleArray>::New();
        smth->SetName("smth");

        auto vel = vtkSmartPointer<vtkDoubleArray>::New();
        vel->SetName("velocity");
        vel->SetNumberOfComponents(3);

        for (const CalcNode &node : nodes)
        {
            pts->InsertNextPoint(node.x, node.y, node.z);
            double v[3] = {node.vx, node.vy, node.vz};
            vel->InsertNextTuple(v);
            smth->InsertNextValue(node.smth);
        }

        grid->SetPoints(pts);
        grid->GetPointData()->AddArray(vel);
        grid->GetPointData()->AddArray(smth);

        for (const Element &elem : elements)
        {
            auto tetra = vtkSmartPointer<vtkTetra>::New();
            tetra->GetPointIds()->SetId(0, elem.nodesIds[0]);
            tetra->GetPointIds()->SetId(1, elem.nodesIds[1]);
            tetra->GetPointIds()->SetId(2, elem.nodesIds[2]);
            tetra->GetPointIds()->SetId(3, elem.nodesIds[3]);
            grid->InsertNextCell(tetra->GetCellType(), tetra->GetPointIds());
        }

        std::string fileName = "output/tetr3d-step-" + std::to_string(snap_number) + ".vtu";
        vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        writer->SetFileName(fileName.c_str());
        writer->SetInputData(grid);
        writer->Write();
    }
};

int main()
{
    // Параметры моделирования
    //------
    const double totalTime = 4.0; // общее время (сек)
    const int numSteps = 500;     // количество кадров
    const double dt = totalTime / numSteps;

    const unsigned int GMSH_TETR_CODE = 4;
    //------

    gmsh::initialize();
    gmsh::model::add("custom_motion_object");

    try
    {
        gmsh::merge("../shell.stl");
    }
    catch (...)
    {
        gmsh::logger::write("Could not load STL mesh: bye!");
        gmsh::finalize();
        return -1;
    }

    // Восстанавливаем геометрию и строим объёмную сетку
    double angle = 40;
    bool forceParametrizablePatches = false;
    bool includeBoundary = true;
    double curveAngle = 180;
    gmsh::model::mesh::classifySurfaces(angle * M_PI / 180., includeBoundary, forceParametrizablePatches, curveAngle * M_PI / 180.);
    gmsh::model::mesh::createGeometry();

    std::vector<std::pair<int, int>> surfaces;
    gmsh::model::getEntities(surfaces, 2);
    std::vector<int> surfaceTags;
    for (auto &s : surfaces)
        surfaceTags.push_back(s.second);
    int loop = gmsh::model::geo::addSurfaceLoop(surfaceTags);
    gmsh::model::geo::addVolume({loop});
    gmsh::model::geo::synchronize();

    // Задаём размер элементов
    int f = gmsh::model::mesh::field::add("MathEval");
    gmsh::model::mesh::field::setString(f, "F", "8");
    gmsh::model::mesh::field::setAsBackgroundMesh(f);

    gmsh::model::mesh::generate(3);

    // Извлекаем узлы
    std::vector<double> nodesCoord;
    std::vector<size_t> nodeTags;
    std::vector<double> parametricCoord;
    gmsh::model::mesh::getNodes(nodeTags, nodesCoord, parametricCoord);

    // Извлекаем тетраэдры
    std::vector<int> elementTypes;
    std::vector<std::vector<size_t>> elementTags;
    std::vector<std::vector<size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags);

    std::vector<size_t> *tetrNodesTags = nullptr;
    for (size_t i = 0; i < elementTypes.size(); i++)
    {
        if (elementTypes[i] == GMSH_TETR_CODE)
        {
            tetrNodesTags = &elementNodeTags[i];
            break;
        }
    }

    if (!tetrNodesTags)
    {
        cerr << "No tetrahedra found. Exiting." << endl;
        gmsh::finalize();
        return -2;
    }

    cout << "Model has " << nodeTags.size() << " nodes and "
         << tetrNodesTags->size() / 4 << " tetrahedra." << endl;

    // Проверка индексов (как в оригинале)
    for (size_t i = 0; i < nodeTags.size(); i++)
        assert(i == nodeTags[i] - 1);
    assert(tetrNodesTags->size() % 4 == 0);

    // Создаём расчётную сетку
    CalcMesh mesh(nodesCoord, *tetrNodesTags);

    // Освобождаем gmsh
    gmsh::finalize();

    if (std::filesystem::exists("output"))
    {
        std::filesystem::remove_all("output");
    }
    std::filesystem::create_directories("output");

    // Записываем начальное состояние
    mesh.snapshot(0);

    double t = 0.0;
    for (int step = 1; step <= numSteps; ++step)
    {
        t += dt;
        mesh.step(t, dt, step);
        std::cout << "step: " << step << std::endl;
    }

    cout << "Done. Generated " << numSteps + 1 << " VTU files." << endl;
    return 0;
}