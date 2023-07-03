#include "Scene.h"

#include "Numeric/PCGSolver.h"

namespace asa
{
namespace
{
////////////////////////////////////////////////
// Add any reusable classes or functions HERE //
////////////////////////////////////////////////
}  // namespace

// advection
void Fluid2::fluidAdvection(const float dt)
{
    {
        // Ink advection HERE

        auto inkRGB_ = this->inkRGB;

        for (size_t i = 0; i < this->grid.getSize().x; i++) 
        {
            for (size_t j = 0; j < this->grid.getSize().y; j++) 
            {
                Vector2 x = grid.getCellPos(Index2(i, j));
        
                auto vx = (velocityX.getValue(i + 1, j) + velocityX.getValue(i, j)) * 0.5f;
                auto vy = (velocityY.getValue(i, j + 1) + velocityY.getValue(i, j)) * 0.5f;
                Vector2 v(vx, vy);
        
                Vector2 x_1 = x - (dt * v);
                Vector2 indices = grid.getCellIndex(x_1);

                indices.x = clamp(indices.x, 0, grid.getSize().x - 1);
                indices.y = clamp(indices.y, 0, grid.getSize().y - 1);
        
                auto aa = inkRGB_.getValue(floor(indices.x), floor(indices.y));
                auto ab = inkRGB_.getValue(floor(indices.x), floor(indices.y) + 1);
                auto ba = inkRGB_.getValue(floor(indices.x) + 1, floor(indices.y));
                auto bb = inkRGB_.getValue(floor(indices.x) + 1, floor(indices.y) + 1);

                auto s = indices.x - floor(indices.x);
                auto t = indices.y - floor(indices.y);
        
                auto newInk = bilerp(aa, ab, ba, bb, s, t);
        
                inkRGB.setValue(i, j, newInk);
            }
        }
    }

    {
        auto velocityX_ = velocityX;
        
        // Velocity X
        for (size_t i = 0; i < velocityX.getSize().x; i++) 
        {
            for (size_t j = 0; j < velocityX.getSize().y; j++) 
            {
                Vector2 x = grid.getFacePosX(Index2(i, j));

                auto vx = velocityX.getValue(i, j);
                auto vy = velocityY.getValue(grid.getFaceIndexY(x).x, grid.getFaceIndexY(x).y);
                Vector2 v(vx, vy);

                Vector2 x_1 = x - (dt * v);
                Vector2 indices = grid.getFaceIndexX(x_1);

                indices.x = clamp(indices.x, 0, velocityX.getSize().x - 1);
                indices.y = clamp(indices.y, 0, velocityX.getSize().y - 1);

                auto aa = velocityX_.getValue(floor(indices.x), floor(indices.y));
                auto ab = velocityX_.getValue(floor(indices.x), floor(indices.y) + 1);
                auto ba = velocityX_.getValue(floor(indices.x) + 1, floor(indices.y));
                auto bb = velocityX_.getValue(floor(indices.x) + 1, floor(indices.y) + 1);

                auto s = indices.x - floor(indices.x);
                auto t = indices.y - floor(indices.y);

                auto newVelocityX = bilerp(aa, ab, ba, bb, s, t);

                velocityX.setValue(i, j, newVelocityX);
            }
        }

        auto velocityY_ = velocityY;

        // Velocity Y
        for (size_t i = 0; i < velocityY.getSize().x; i++) 
        {
            for (size_t j = 0; j < velocityY.getSize().y; j++) 
            {
                Vector2 x = grid.getFacePosX(Index2(i, j));

                auto vx = velocityX.getValue(grid.getFaceIndexX(x).x, grid.getFaceIndexX(x).y);
                auto vy = velocityY.getValue(i, j);
                Vector2 v(vx, vy);

                Vector2 x_1 = x - (dt * v);
                Vector2 indices = grid.getFaceIndexY(x_1);

                indices.x = clamp(indices.x, 0, velocityY.getSize().x - 1);
                indices.y = clamp(indices.y, 0, velocityY.getSize().y - 1);

                auto aa = velocityY_.getValue(floor(indices.x), floor(indices.y));
                auto ab = velocityY_.getValue(floor(indices.x), floor(indices.y) + 1);
                auto ba = velocityY_.getValue(floor(indices.x) + 1, floor(indices.y));
                auto bb = velocityY_.getValue(floor(indices.x) + 1, floor(indices.y) + 1);

                auto s = indices.x - floor(indices.x);
                auto t = indices.y - floor(indices.y);

                auto newVelocityY = bilerp(aa, ab, ba, bb, s, t);

                velocityY.setValue(i, j, newVelocityY);
            }
        }
    }
}

void Fluid2::fluidEmission()
{
    if (Scene::testcase >= Scene::SMOKE) 
    {
        uint lowYIndex = 0.2f * inkRGB.getSize().y;

        for (size_t i = 0; i < inkRGB.getSize().x; i++) 
        {
            Vector3 newColor = i < inkRGB.getSize().x * 0.5f ? Vector3(1, 0, 0) : Vector3(0, 1, 0);

            for (size_t j = 0; j < lowYIndex; j++) 
            {
                inkRGB.setValue(Index2(i, j), newColor);
            }
        }

        for (size_t i = 0; i < velocityY.getSize().x; i++) 
        {
            for (size_t j = 0; j < lowYIndex + 1; j++) 
            {
                velocityY.setValue(i, j, 1);
            }
        }
    }
}

void Fluid2::fluidVolumeForces(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) 
    {
        // Aplicar la gravedad
        for (size_t i = 0; i < this->velocityY.getSize().x; i++) 
        {
            for (size_t j = 0; j < this->velocityY.getSize().y; j++) 
            {                
                velocityY.setValue(i, j, velocityY.getValue(i, j) + dt / Scene::kDensity * Scene::kGravity);
            }
        }
    }
}

void Fluid2::fluidViscosity(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) 
    {
        for (size_t i = 1; i < this->grid.getSize().x; i++) 
        {
            for (size_t j = 1; j < this->grid.getSize().y; j++) 
            {
                float stepVelocityX = dt / Scene::kDensity *
                    ((velocityX.getValue(i + 1, j) - 2 * velocityX.getValue(i, j) + velocityX.getValue(i - 1, j) / pow(grid.getDx().x, 2)) +
                     (velocityX.getValue(i, j + 1) - 2 * velocityX.getValue(i, j) + velocityX.getValue(i, j - 1) / pow(grid.getDx().y, 2)));

                velocityX.setValue(i, j, velocityX[i, j] + stepVelocityX);
                
                float stepVelocityY = dt / Scene::kDensity *
                    ((velocityY.getValue(i + 1, j) - 2 * velocityY.getValue(i, j) + velocityY.getValue(i - 1, j) / pow(grid.getDx().x, 2)) +
                     (velocityY.getValue(i, j + 1) - 2 * velocityY.getValue(i, j) + velocityY.getValue(i, j - 1) / pow(grid.getDx().y, 2)));

                velocityY.setValue(i, j, velocityY[i, j] + stepVelocityY);
            }
        }
    }
}

void Fluid2::fluidPressureProjection(const float dt)
{
    if (Scene::testcase >= Scene::SMOKE) 
    {
        // Incompressibility / Pressure term HERE
        // 
        // Las presiones y tinta se quedán en el centro de celda, lo que se desplazan son las velocidades en x e y
        // 
        // Para las variables de contorno usar solidas, restando -1 a cada elemento que esté pegado a la pared de la matriz, -2 en el caso de las esquinas

        // Make walls normal velocities = 0, todas las las velocidades paredes en las matrices seteadas a 0

        for (size_t i = 0; i < this->grid.getSizeFacesX().x; i++) 
        {
            velocityX[i, 0] = 0;
            velocityX[i, velocityX.getSize()] = 0;
        
            velocityX[0, i] = 0;
            velocityX[velocityX.getSize(), i] = 0;
        
            velocityY[i, 0] = 0;
            velocityY[i, velocityY.getSize()] = 0;
        
            velocityY[0, i] = 0;
            velocityY[velocityY.getSize(), i] = 0;
        }
        
        // Compute RHS
        
        //std::vector<double> b(0);
        //b.resize(grid.getSize().x * grid.getSize().y + 1);
        //
        //
        //for (size_t i = 1; i < this->grid.getSize().x + 1; i++) 
        //{
        //    for (size_t j = 1; j < this->grid.getSize().y + 1; j++) 
        //    {
        //        
        //
        //        b[velocityX.getLinearIndex(i, j)] = velocityX[i, j] - velocityX[i - 1, j - 1];
        //
        //        // En b falta añadir las velocidades en y pero en que indices?
        //
        //        //Usar array2 getlinearindex funcion?
        //
        //    }
        //}
        //
        //std::vector<double> result(0);
        //result.resize(grid.getSize().x * grid.getSize().y);
        //
        //SparseMatrix<double> A(grid.getSize().x * grid.getSize().y, 5);
        //
        //for (size_t i = 0; i < this->grid.getSize().x * this->grid.getSize().y + 1; i++) 
        //{
        //    A.set_element(i, i, 2);
        //
        //    if (i + 1 > grid.getSize().x + 1) 
        //    {
        //        A.add_to_element(i, i, -1);
        //    } 
        //    else 
        //    {
        //        A.set_element(i + 1, i, -1);
        //    }
        //
        //    if (i + 1 > grid.getSize().y + 1) 
        //    {
        //        A.add_to_element(i, i, -1);
        //    } 
        //    else 
        //    {
        //        A.set_element(i, i + 1, -1);
        //    }
        //}
        //
        //PCGSolver<double> solver;
        //solver.set_solver_parameters(0.000001, 200);
        //
        //double residual_out;
        //int iterations_out;
        //
        //// Falta añadir las partes que multiplican a A y b antes de resolver la ecuación
        //
        //solver.solve(A, b, result, residual_out, iterations_out);
        //
        //
        //
        ////Aplicar solucion a cada una de las velocidades
        //// Aplicar en los dos ejes por separado
        //// Aplicar los gradientes de presion a los resultados
        //// No hace falta aplicar velocidades a los bordes?
        //
        //for (size_t i = 1; i < this->grid.getSize().x - 1; i++) 
        //{
        //    for (size_t j = 1; j < this->grid.getSize().y - 1; j++) 
        //    {
        //        // Usar array2 getlinearindex funcion?
        //
        //        //b[i + j] = ;
        //    }
        //}
    }
}
}  // namespace asa
